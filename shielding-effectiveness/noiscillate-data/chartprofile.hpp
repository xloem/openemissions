#pragma once

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TROOT.h>
#include <TSystem.h>

#include "types.hpp"
#include "stats.hpp"

        // TODO: use a 2D histogram? user may understand signal better

template <typename _Profile>
class RootChartProfile
{
public:
  using Profile = _Profile;
  using Scalar = typename Profile::Scalar;

  RootChartProfile(char const * name, TVirtualPad * pad = 0)
  : _hidden(true),
    _pad(pad)
  {
    _graph.SetTitle(name);
    _graph.GetXaxis()->SetTitle("Frequency (Hz)");
    _graph.GetYaxis()->SetTitle("Power (dB)");
    _graph.SetEditable(kFALSE);
    _graph.SetHighlight(kTRUE);
  }

  void init(Int_t count, bool hidden = false, TVirtualPad * pad = 0)
  {
    _hidden = hidden;
    _graph.Set(count);
    if (pad)
    {
      _pad = pad;
    }
    if (!_pad && !hidden)
    {
      if (!gPad)
      {
        gROOT->MakeDefCanvas();
      }
      _pad = gPad;
      _pad->SetTitle(_graph.GetName());
    }
  }

  void resetMetrics()
  {
    _maxDelta = 0;
    _maxError = 0;
    _minValue = Eigen::NumTraits<Scalar>::infinity();
    _maxValue = -_minValue;
    _errorStats_v.clear();
    _errorStats_dB.clear();
  }

  void prepPoint(size_t idx, Scalar xmin, Scalar x, Scalar xmax)
  {
    _graph.SetPoint(idx, x, -50);
    _graph.SetPointError(idx, x - xmin, xmax - x, 50, 50);
  }

  void prepPoints(Int_t count, Scalar start, Scalar hopSize, bool hidden = false)
  {
    init(count, hidden);
    auto delta = hopSize / 2;
    for (Int_t i = 0; i < count; ++ i)
    {
      auto x = start + hopSize * i;
      prepPoint(i, x - delta, x, x + delta);
    }
  }

  void prepPoints(Profile & profile, bool hidden = false)
  {
    init(profile.size(), hidden);
    auto it = profile.begin();
    Int_t idx = 0;
    
    while (it != profile.end())
    {
      Scalar freql = 0, freqh = 0;
      auto freq = *it;
      // - [ ] prepPoint each one individually
      //       set only if not begin, only if not about end
      //       if not set, freql/h will be zero, set to other high/low value
      //       below is mostly to be removed or barely started
      if (it != profile.begin())
      {
        -- it;
        freql = (freq - *it) / 2;
        ++ it;
      }
      ++ it;
      if (it != profile.end())
      {
        freqh = (*it - freq) / 2;
      }
      if (freql == 0)
      {
        freql = freqh;
      }
      if (freqh == 0)
      {
        freqh = freql;
      }
      // freqh, freql, freq set for this point.  *it points to next point
      prepPoint(idx, freql, freq, freqh);
      ++ idx;
    }
  }

  void paint(Profile & profile)
  {
    auto count = _graph.GetN();
    _drawMode = "AL";
    resetMetrics();
    for (size_t i = 0; i < count; ++ i)
    {
      try
      {
        setPointInternal(profile, i);
      }
      catch (std::out_of_range)
      {
        //_drawMode = "AP";
      }
    }
    if (!_hidden)
    {
      _pad->cd();
      _graph.Draw(_drawMode);
      if (_drawMode[1] == 'L')
      {
        _graph.SetMinimum(_minValue);
        _graph.SetMaximum(_maxValue);
      }
      else
      {
        _graph.SetMinimum(-100);
      }
      _graph.Paint(_drawMode);
      _pad->Update();
      gSystem->ProcessEvents();
    }
  }

  void finalize()
  {
    _drawMode = "AL";
    _graph.SetMinimum(_minValue);
    _graph.SetMaximum(_maxValue);
  }

  void paint(Profile & profile, Int_t idx)
  {
    setPointInternal(profile, idx);
    if (!_hidden)
    {
      _pad->cd();
      _graph.Paint(_drawMode);
      _pad->Update();
    }
    gSystem->ProcessEvents();
    if (!_hidden)
    {
      if (_pad->GetCanvasImp() == nullptr)
      {
        throw std::runtime_error("graph closed");
      }
    }
  }

  TVirtualPad * & pad()
  {
    return _pad;
  }

  Scalar maxDelta() const { return _maxDelta; }
  Scalar maxDeltaFreq() const { return _maxDeltaFreq; }
  Scalar maxError() const { return _maxError; }
  Scalar maxErrorFreq() const { return _maxErrorFreq; }
  Scalar maxValue() const { return _maxValue; }
  Scalar maxValueFreq() const { return _maxValueFreq; }
  Scalar minValue() const { return _minValue; }
  Scalar minValueFreq() const { return _minValueFreq; }
  StatsAccumulator<Scalar> const & errorStats_raw() const { return _errorStats_v; }
  StatsAccumulator<Scalar> const & errorStats_dB() const { return _errorStats_dB; }

private:
  void setPointInternal(Profile & profile, Int_t idx)
  {
    auto xs = _graph.GetX();
    auto ys = _graph.GetY();
    auto eyls = _graph.GetEYlow();
    auto eyhs = _graph.GetEYhigh();
    auto freq = xs[idx];

    auto rawY = profile.bestMag(freq);
    auto rawVar = profile.bestMagVariance(freq);
    auto rawErr = sqrt(rawVar) * 3;
    FixedVector<Double_t,3> unproc, proc;
    unproc[0] = 0;
    unproc[1] = -rawErr;
    unproc[2] = rawErr;
    proc = (unproc.array() + rawY).log10() * 20;
    proc.tail(2).array() -= proc[0];
    proc[1] = -proc[1];
    auto delta = abs(proc[0] - ys[idx]);
    if (delta > _maxDelta)
    {
      _maxDelta = delta;
      _maxDeltaFreq = freq;
    }
    if (proc[0] > _maxValue)
    {
      _maxValue = proc[0];
      _maxValueFreq = freq;
    }
    if (proc[0] < _minValue)
    {
      _minValue = proc[0];
      _minValueFreq = freq;
    }
    if (proc[1] > _maxError)
    {
      _maxError = proc[1];
      _maxErrorFreq = freq;
    }
    if (proc[2] > _maxError)
    {
      _maxError = proc[2];
      _maxErrorFreq = freq;
    }
    _errorStats_v.add(rawErr);
    _errorStats_dB.add(proc[1]);
    _errorStats_dB.add(proc[2]);
    ys[idx] = proc[0];
    eyls[idx] = proc[1];
    eyhs[idx] = proc[2];
  }

  bool _hidden;
  TVirtualPad * _pad;
  TGraphAsymmErrors _graph;
  char const * _drawMode;

  Scalar _maxDelta, _maxDeltaFreq;
  Scalar _maxError, _maxErrorFreq;
  Scalar _maxValue, _maxValueFreq;
  Scalar _minValue, _minValueFreq;
  StatsAccumulator<Scalar> _errorStats_v;
  StatsAccumulator<Scalar> _errorStats_dB;
};
