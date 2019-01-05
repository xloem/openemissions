// 1-1-profile-env.cpp
//
// Produces a histogram noise profile for an environment, using a software-defined radio.
//
// The output histograms are of the raw values received from the radio,
// so are expected to have a minimal mean but a variance equal to the noise power.

#include <fstream>
#include <iostream>
#include <random>

#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TSysEvtHandler.h>
#include <TSystem.h>

#include "datafeeder.hpp"
#include "processornoiseprofile.hpp"
#include "soapylive.hpp"

using Scalar = double;

class T11ProfileEnv : public TApplication
{
public:
  T11ProfileEnv(Int_t * argc, char ** argv)
  : TApplication("T11ProfileEnv", argc, argv, 0, -1),
    _sigQuit(*this, kSigQuit),
    _sigInterrupt(*this, kSigInterrupt),
    _sigTermination(*this, kSigTermination),
    _minTuneFreq(200000000),
    _maxTuneFreq(1400000000),
    _sampleRate(2400000),
    _gain(43),
    _hopDistance(0),
    _randGen(_randDev()),
    _terminate(false)
  {
    GetOptions(argc, argv);
    
    _canvas.SetTitle("Recording Profile");
    _graph.SetTitle("Recording Profile");
    _graph.GetXaxis()->SetTitle("Frequency (Hz)");
    _graph.GetYaxis()->SetTitle("Power (dB)");
    _graph.SetEditable(kFALSE);
    //_graph.SetHighlight(kTRUE);

    _sigQuit.Add();
    _sigInterrupt.Add();
    _sigTermination.Add();
  }

  virtual void HandleException(Int_t sig = -1) override
  {
    std::cerr << "-exception-" << std::endl;
    _terminate = true;
  }

  virtual void Run(Bool_t retrn) override
  {
    Int_t status = 0;
    _terminate = false;
    std::cerr << std::fixed;
    std::cerr << "Device construction args: '" << _soapyArgs << "'." << std::endl;
    std::cerr << "Recording at " << _sampleRate << " Hz in hops of " << _hopDistance << " Hz from " << _minTuneFreq << " Hz through " << _maxTuneFreq << " Hz." << std::endl;

    SoapyLive data(_soapyArgs.c_str(), _gain, _sampleRate);
    auto ident = data.ident();
    std::cerr << "Device identification string: '" << ident << "'." << std::endl;
    ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE> processor(_fname, {data.epsilon(), data.midval()}, _env, _emits, {"SoapySDR", ident});

    for (Scalar freq = _minTuneFreq; freq <= _maxTuneFreq; freq += _hopDistance)
    {
      _TuningFreq obj;
      obj.freq = freq;//data.tune(freq);
      obj.idx = _tuneFreqs.size();
      _tuneFreqs.push_back(obj);
    }
    _graph.Set(_tuneFreqs.size());
    for (size_t i = 0; i < _tuneFreqs.size(); ++ i)
    {
      _graph.SetPoint(i, _tuneFreqs[i].freq, 0);
      _graph.SetPointError(i, -_sampleRate / 2, _sampleRate / 2, 0, 0);
    }

    _drawMode = "AL";
    for (auto & freq : _tuneFreqs)
    {
      try
      {
        setPoint(processor, freq.freq, freq.idx);
      }
      catch (std::out_of_range)
      {
        _drawMode = "AP";
      }
    }
    _graph.Draw(_drawMode);
    _graph.SetMinimum(-100);
    _graph.Draw(_drawMode);
    _canvas.Update();

    try
    {

      std::cerr << "Recording environment " << processor.environment() << "." << std::endl;
      auto & emits = processor.emitters();
      std::cerr << emits.size() << " emitter(s) active" << (emits.empty() ? "." : ": ");
      for (auto & emit : emits)
      {
        std::cerr << " " << emit;
      }
      std::cerr << std::endl;
  
      DataFeeder<Scalar, decltype(data), decltype(processor)> feeder(data, processor);
      HeapVector<Complex> buffer(_sampleRate / 8);
      RecBufMeta bufMeta;
  
      std::shuffle(_tuneFreqs.begin(), _tuneFreqs.end(), _randGen);
      auto curFreq_it = _tuneFreqs.begin();
      curFreq_it->freq = data.tune(curFreq_it->freq);
  
      while (true)
      {
        data.readMany(buffer, bufMeta);
        gSystem->ProcessEvents();
        if (buffer.size() == 0) continue;
        feeder.add(buffer.array().real(), bufMeta);
  
        auto freq = curFreq_it->freq;
        auto freqIdx = curFreq_it->idx;
  
        ++ curFreq_it;
        if (curFreq_it == _tuneFreqs.end())
        {
          std::shuffle(_tuneFreqs.begin(), _tuneFreqs.end(), _randGen);
          curFreq_it = _tuneFreqs.begin();
        }
        // tune first to provide more time for radio to settle
        //    TODO might need to change code elsewhere to ensure this helps
        curFreq_it->freq = data.tune(curFreq_it->freq);
  
        // process chart
        // TODO: use a 2D histogram? user may understand signal better
        setPoint(processor, freq, freqIdx);
        if (gPad)
        {
          _graph.Draw(_drawMode);
          _canvas.Update();
        }
        gSystem->ProcessEvents();
        if (_canvas.GetCanvasImp() == nullptr) _terminate = true;
  
        // write file if a full sweep has completed
        if (curFreq_it == _tuneFreqs.begin() || _terminate)
        {
          std::cerr << "-sto";
          processor.write();
          std::cerr << "re-";
          gSystem->ProcessEvents();
          if (_terminate) return;
          _drawMode = "AL";
        }
      }
    }
    catch(...)
    {
      std::cerr << "-error-" << std::endl;
      std::cerr << "-sto";
      processor.write();
      std::cerr << "re-";
      throw;
    }
    gSystem->ProcessEvents();
    processor.write();
    gSystem->ProcessEvents();
    Terminate(0);
  }

  template <typename T>
  void setPoint(T & processor, Scalar freq, size_t idx)
  {
    auto rawY = processor.bestMag(freq);
    auto rawVar = processor.bestMagVariance(freq);
    auto rawErr = sqrt(rawVar) * 3;
    FixedVector<Double_t,3> unproc, proc;
    unproc[0] = 0;
    unproc[1] = -rawErr;
    unproc[2] = rawErr;
    proc = (unproc.array() + rawY).log10() * 20;
    proc.tail(2).array() -= proc[0];
    _graph.GetY()[idx] = proc[0];
    _graph.GetEYlow()[idx] = proc[1];
    _graph.GetEYhigh()[idx] = proc[2];
  }

  virtual void GetOptions(Int_t *argc, char **argv) override
  {
    bool setEnv = false;
    try
    {
      for (int i = 1; i < *argc; ++ i)
      {
        if (argv[i] == nullptr)
        {
          continue;
        }
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help") || !strcmp(argv[i], "-?"))
        {
          argv[i] = nullptr;
          PrintHelp();
          return;
        }
        else if (!strcmp(argv[i], "-o"))
        {
          argv[i] = nullptr;
          if (++ i >= *argc) throw std::invalid_argument("Output filename missing.");
          std::cerr << "Writing output to '" << argv[i] << "'." << std::endl;
          _fname = argv[i];
          argv[i] = nullptr;
        }
        else if (!strcmp(argv[i], "-d"))
        {
          argv[i] = nullptr;
          if (++ i >= *argc) throw std::invalid_argument("SoapySDR args missing.");
          _soapyArgs = argv[i];
          argv[i] = nullptr;
        }
        else if (!strcmp(argv[i], "-f"))
        {
          argv[i] = nullptr;
          if (++ i >= *argc) throw std::invalid_argument("Min freq missing.");
          _minTuneFreq = std::stod(argv[i]);
          argv[i] = nullptr;
          if (++ i >= *argc) throw std::invalid_argument("Max freq missing.");
          _maxTuneFreq = std::stod(argv[i]);
          argv[i] = nullptr;
        }
        else if (!strcmp(argv[i], "-r"))
        {
          argv[i] = nullptr;
          if (++ i >= *argc) throw std::invalid_argument("Samplerate missing.");
          _sampleRate = std::stoll(argv[i]);
          argv[i] = nullptr;
        }
        else if (!strcmp(argv[i], "-g"))
        {
          argv[i] = nullptr;
          if (++ i >= *argc) throw std::invalid_argument("Gain missing.");
          _gain = std::stod(argv[i]);
          argv[i] = nullptr;
        }
        else if (!strcmp(argv[i], "-h"))
        {
          argv[i] = nullptr;
          if (++ i >= *argc) throw std::invalid_argument("Hop distance missing.");
          _hopDistance = std::stod(argv[i]);
          if (_hopDistance < 1) throw std::invalid_argument("Hop distance must be positive.");
          argv[i] = nullptr;
        }
        else if (!setEnv)
        {
          _env = std::stoll(argv[i]);
          argv[i] = nullptr;
          setEnv = true;
        }
        else
        {
          _emits.push_back(std::stoll(argv[i]));
          argv[i] = nullptr;
        }
      }

      if (!_fname.size())
      {
        throw std::invalid_argument("Output filename is required.");
      }

      if (!setEnv)
      {
        throw std::invalid_argument("Environment number is required.");
      }

      if (_hopDistance == 0)
      {
        _hopDistance = _sampleRate / 2;
      }
    }
    catch(std::invalid_argument e)
    {
      std::cerr << "Invalid argument: " << e.what() << std::endl;
      PrintHelp();
    }
    catch(std::out_of_range)
    {
      std::cerr << "Value out of range." << std::endl;
      Terminate(-1);
    }
  }

private:
  void PrintHelp()
  {
    std::cerr                                                                                    << std::endl;
    std::cerr << "Records from an SDR and stores a histogram profile of the values received,"    << std::endl;
    std::cerr << "useful for statistical analysis of the noise."                                 << std::endl;
    std::cerr                                                                                    << std::endl;
    std::cerr << "Usage: 1-1-prof-env [options ...] -o a.noisep env [emit1 ... emitN]"           << std::endl;
    std::cerr                                                                                    << std::endl;
    std::cerr << "            env : An int representing the environment and setup"               << std::endl;
    std::cerr << "          emitX : Ints representing all active emitter setups"                 << std::endl;
    std::cerr                                                                                    << std::endl;
    std::cerr << "Options:"                                                                      << std::endl;
    std::cerr << "    -o a.noisep : create or update the file a.noisep"                          << std::endl;
    std::cerr << "   -d \"args...\" : SoapySDR device construction args (default: \"\")"         << std::endl;
    std::cerr << " -f minHz maxHz : hop between the passed range of tuning frequencies in hertz" << std::endl;
    std::cerr << "                  (default: 2000000000 1400000000)"                            << std::endl;
    std::cerr << "        -r rate : samplerate to record at (default: 2400000)"                  << std::endl;
    std::cerr << "      -g gaindB : dB gain for hardware amplifier (default: 43)"                << std::endl;
    std::cerr << "       -h hopHz : distance between frequencies to tune to"                     << std::endl;
    std::cerr << "                  (default: samplerate / 2)"                                   << std::endl;
    std::cerr                                                                                    << std::endl;
    Terminate(-1);
  }


  class THaltSignalHandler : public TSignalHandler
  {
  public:
    THaltSignalHandler(T11ProfileEnv & target, ESignals sig)
    : TSignalHandler(sig),
      target(target)
    { }

    virtual Bool_t Notify() override
    {
      target.HandleException();
    }

  private:
    T11ProfileEnv & target;
  };

  struct _TuningFreq
  {
    Scalar freq;
    size_t idx;
  };

  THaltSignalHandler _sigQuit, _sigInterrupt, _sigTermination;

  std::string _fname;
  long long _env;
  std::vector<int64_t> _emits;

  char const * _drawMode;

  double _minTuneFreq;
  double _maxTuneFreq;
  long long _sampleRate;
  double _gain;
  double _hopDistance;
  std::string _soapyArgs;

  std::random_device _randDev;
  std::mt19937 _randGen;
  std::vector<_TuningFreq> _tuneFreqs{};

  bool _terminate;

  TCanvas _canvas;
  TGraphAsymmErrors _graph;
};

//ClassImp(T11ProfileEnv)

int main(int argc, char ** argv)
{
  T11ProfileEnv app(&argc, argv);
  app.Run(false);
}
