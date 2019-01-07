// 1-1-profile-env.cpp
//
// Produces a histogram noise profile for an environment, using a software-defined radio.
//
// The output histograms are of the raw values received from the radio,
// so are expected to have a minimal mean but a variance equal to the noise power.

#include <fstream>
#include <iostream>
#include <random>

#include "datafeeder.hpp"
#include "processornoiseprofile.hpp"
#include "rootapplication.hpp"
#include "soapylive.hpp"

using Scalar = double;

class T11ProfileEnv : public RootApplication
{
public:
  T11ProfileEnv(Int_t * argc, char ** argv)
  : RootApplication(
      "1-1-prof-env",
      "Records from an SDR and stores a histogram profile of the values received,\n"
      "useful for statistical analysis of the noise.",
      argc, argv),
    _sweepCount(0),
    _minTuneFreq(200000000),
    _maxTuneFreq(1400000000),
    _sampleRate(2400000),
    _gain(43),
    _hopDistance(0),
    _settleSecs(0.03125),
    _randGen(_randDev()),
    _sweepNum(1)
  { }

  virtual void Run(Bool_t retrn) override
  {
    Int_t status = 0;
    _terminate = false;
    std::cerr << std::fixed;
    std::cerr << "Device construction args: '" << _soapyArgs << "'." << std::endl;
    std::cerr << "Recording at " << _sampleRate << " Hz in hops of " << _hopDistance << " Hz from " << _minTuneFreq << " Hz through " << _maxTuneFreq << " Hz." << std::endl;

    SoapyLive data(_soapyArgs.c_str(), _gain, _sampleRate, _settleSecs);
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
    if (quiet < 1)
    {
      canvas.reset(new TCanvas());
      canvas->SetTitle("Recording Profile");
    }
    graph.Set(_tuneFreqs.size());
    for (size_t i = 0; i < _tuneFreqs.size(); ++ i)
    {
      graph.SetPoint(i, _tuneFreqs[i].freq, 0);
      graph.SetPointError(i, -_sampleRate / 2, _sampleRate / 2, 0, 0);
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
    if (quiet < 1)
    {
      graph.Draw(_drawMode);
      graph.SetMinimum(-100);
      graph.Draw(_drawMode);
      canvas->Update();
    }

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
      HeapVector<Complex> buffer(_sampleRate / 16);//512);////8);
      RecBufMeta bufMeta;
  
      std::shuffle(_tuneFreqs.begin(), _tuneFreqs.end(), _randGen);
      auto curFreq_it = _tuneFreqs.begin();
      curFreq_it->freq = data.tune(curFreq_it->freq);

      _maxDelta = 0;
      _maxError = 0;
      _minValue = Eigen::NumTraits<Scalar>::infinity();
      _maxValue = -_minValue;

      while (true)
      {
        data.readMany(buffer, bufMeta);
        tick();
        if (buffer.size() == 0) continue;
        feeder.add(buffer.array().real(), bufMeta);
  
        auto freq = curFreq_it->freq;
        auto freqIdx = curFreq_it->idx;

        if (quiet < 1)
        {
          std::cerr << "\rSweep #" << _sweepNum << ": " << (curFreq_it - _tuneFreqs.begin()) * 100 / (_tuneFreqs.size() - 1) << "% (last @" << freq << " Hz)   " << std::flush;
        }
  
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
        if (quiet < 1 && gPad)
        {
          graph.Draw(_drawMode);
          canvas->Update();
        }
        tick();
        if (quiet < 1 && canvas->GetCanvasImp() == nullptr) _terminate = true;
  
        // write file if a full sweep has completed
        if (curFreq_it == _tuneFreqs.begin() || _terminate)
        {
          std::cerr << "-sto";
          processor.write();
          std::cerr << "re-";
          tick();
          std::cerr << std::endl;
          std::cerr << "Max Delta = " << _maxDelta << " dB at " << _maxDeltaFreq << " Hz" << std::endl;
          std::cerr << "Max Error = " << _maxError << " dB at " << _maxErrorFreq << " Hz" << std::endl;
          std::cerr << "Max Value = " << _maxValue << " dB at " << _maxValueFreq << " Hz" << std::endl;
          std::cerr << "Min Value = " << _minValue << " dB at " << _minValueFreq << " Hz" << std::endl;

          StatsAccumulator<Scalar> errorStats_v;
          StatsAccumulator<Scalar> errorStats_dB;
          for (auto & f : _tuneFreqs)
          {
            auto v = processor.bestMag(f.freq);
            auto e = sqrt(processor.bestMagVariance(f.freq)) * 3;
            errorStats_v.add(e);
            auto e2 = graph.GetEYhigh()[f.idx] - graph.GetEYlow()[f.idx];
            errorStats_dB.add(e2);
            //assert(log10(v + e) * 20 - log10(v - e) * 20 == e2);
          }
          std::cerr << "Error_dB Mean = " << errorStats_dB.mean() << " dB" << std::endl;
          std::cerr << "Error_raw Mean = " << std::scientific << errorStats_v.mean() << std::fixed << std::endl;
          std::cerr << "Error_dB Variance = " << errorStats_dB.variance() << " dB" << std::endl;
          std::cerr << "Error_raw Variance = " << std::scientific << errorStats_v.variance() << std::fixed << std::endl;


          if (_sweepCount > 0 && _sweepNum >= _sweepCount) _terminate = true;
          if (_terminate) return;
          _drawMode = "AL";
          ++ _sweepNum;
          _maxDelta = 0;
          _maxError = 0;
          _minValue = Eigen::NumTraits<Scalar>::infinity();
          _maxValue = -_minValue;
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
    tick();
    processor.write();
    tick();
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
    auto delta = abs(proc[0] - graph.GetY()[idx]);
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
    if (-proc[1] > _maxError)
    {
      _maxError = -proc[1];
      _maxErrorFreq = freq;
    }
    if (proc[2] > _maxError)
    {
      _maxError = proc[2];
      _maxErrorFreq = freq;
    }
    graph.GetY()[idx] = proc[0];
    graph.GetEYlow()[idx] = proc[1];
    graph.GetEYhigh()[idx] = proc[2];
  }

  void GetOptions(Int_t *argc, char **argv) override
  {
    regArg({{"-o"},{"a.noisep"}}, {"create or update the file a.noisep"},
      [this](std::string, std::string fname)
      {
        std::cerr << "Writing output to '" << fname << "'." << std::endl;
        _fname = fname;
      }, 1 // required
    );
    regArg({{},{"env"}}, {"an int representing the environment and setup"},
      [this](std::string, std::string env)
      {
        _env = std::stoll(env);
      }, 1 // required
    );
    regArg({{},{"emitX ..."}}, {"ints representing all active emitter setups"},
      [this](std::string, std::string emit)
      {
        _emits.push_back(std::stoll(emit));
      }
    );
    regArg({{"-n"},{"sweepCount"}}, {"stop after this many sweeps", "(default: until interrupted)"},
      [this](std::string, std::string sweeps)
      {
        _sweepCount = std::stoull(sweeps);
      }
    );
    regArg({{"-d"},{"\"args...\""}}, {"SoapySDR device construction args", "(default: \"\")"},
      [this](std::string, std::string args)
      {
        _soapyArgs = args;
      }
    );
    regArg({{"-f"},{"minHz"},{"maxHz"}}, {"range of frequencies to sweep in Hz", "(default: 200000000 1400000000)"},
      [this](std::string, std::string min, std::string max)
      {
        _minTuneFreq = std::stod(min);
        _maxTuneFreq = std::stod(max);
      }
    );
    regArg({{"-r"},{"rate"}}, {"sample rate to record at", "(default: 2400000)"},
      [this](std::string, std::string rate)
      {
        _sampleRate = std::stoll(rate);
      }
    );
    regArg({{"-g"},{"gaindB"}}, {"dB gain for hw amplifier", "(default: 43)"},
      [this](std::string, std::string gain)
      {
        _gain = std::stod(gain);
      }
    );
    regArg({{"-p"},{"hopHz"}}, {"distanec between tuning frequencies", "(default: samplerate / 2)"},
      [this](std::string, std::string hop)
      {
        _hopDistance = std::stod(hop);
      }
    );
    regArg({{"-s"},{"settleSecs"}}, {"time to wait after tuning for state to settle", "(default: 0.03125)"},
      [this](std::string, std::string settle)
      {
        _settleSecs = std::stod(settle);
      }
    );
    //regArg({{"-w"},{"dwellSecs"}}, {"time to dwell on a frequency before hopping to antoher", "(default: ?)"});
    //regArg({{"-w"},{"dBdelta"}}, {"maximum dB error to wait for before hopping"});
    RootApplication::GetOptions(argc, argv);
    if (_hopDistance == 0)
    {
      _hopDistance = _sampleRate / 2;
    }
  }

private:
  struct _TuningFreq
  {
    Scalar freq;
    size_t idx;
  };

  std::string _fname;
  long long _env;
  std::vector<int64_t> _emits;

  char const * _drawMode;

  unsigned long long _sweepCount;
  double _minTuneFreq;
  double _maxTuneFreq;
  long long _sampleRate;
  double _gain;
  double _hopDistance;
  double _settleSecs;
  std::string _soapyArgs;

  std::random_device _randDev;
  std::mt19937 _randGen;
  std::vector<_TuningFreq> _tuneFreqs{};

  bool _terminate;
  size_t _sweepNum;

  Scalar _maxDelta, _maxDeltaFreq;
  Scalar _maxError, _maxErrorFreq;
  Scalar _maxValue, _maxValueFreq;
  Scalar _minValue, _minValueFreq;
};

//ClassImp(T11ProfileEnv)

int main(int argc, char ** argv)
{
  T11ProfileEnv app(&argc, argv);
  app.GetOptions(&argc, argv);
  app.Run(false);
}
