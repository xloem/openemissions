// 1-1-profile-env.cpp
//
// Produces a histogram noise profile for an environment, using a software-defined radio.
//
// The output histograms are of the raw values received from the radio,
// so are expected to have a minimal mean but a variance equal to the noise power.

#include <fstream>
#include <iostream>
#include <random>

#include "chartprofile.hpp"
#include "datafeeder.hpp"
#include "processornoiseprofile.hpp"
#include "rootapplication.hpp"
#include "soapylive.hpp"

using Scalar = double;

class T11ProfileEnv : public RootApplication
{
public:
  using Profile = ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE>;

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
    _sweepNum(1),
    _chart("Noise Power over Frequency")
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
  }

  virtual void Run(Bool_t retrn) override
  {
    std::unique_ptr<TCanvas> canvas;
    Int_t status = 0;
    terminate = false;
    std::cerr << std::fixed;
    std::cerr << "Device construction args: '" << _soapyArgs << "'." << std::endl;
    std::cerr << "Recording at " << _sampleRate << " Hz in hops of " << _hopDistance << " Hz from " << _minTuneFreq << " Hz through " << _maxTuneFreq << " Hz." << std::endl;

    SoapyLive data(_soapyArgs.c_str(), _gain, _sampleRate, _settleSecs);
    auto ident = data.ident();
    std::cerr << "Device identification string: '" << ident << "'." << std::endl;
    Profile processor(_fname, {data.epsilon(), data.midval()}, _env, _emits, {"SoapySDR", ident});

    Scalar freq = _minTuneFreq;
    for (Int_t i = 0; freq <= _maxTuneFreq; ++ i, freq = _minTuneFreq + i * _hopDistance)
    {
      _TuningFreq obj;
      obj.freq = freq;//data.tune(freq);
      obj.idx = _tuneFreqs.size();
      _tuneFreqs.push_back(obj);
    }
    if (quiet < 1)
    {
      canvas.reset(new TCanvas("1-1-prof-env"));
      _chart.pad() = &*canvas;
      canvas->SetTitle("Recording Profile");
    }
    _chart.prepPoints(_tuneFreqs.size(), _minTuneFreq, _hopDistance, quiet > 0);

    _chart.paint(processor);

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

      _chart.resetMetrics();

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
        _chart.paint(processor, freqIdx);
  
        // write file if a full sweep has completed
        if (curFreq_it == _tuneFreqs.begin() || terminate)
        {
          std::cerr << "-sto";
          processor.write();
          std::cerr << "re-";
          _chart.finalize();
          tick();
          std::cerr << std::endl;
          std::cerr << "Max Delta = " << _chart.maxDelta() << " dB at " << _chart.maxDeltaFreq() << " Hz" << std::endl;
          std::cerr << "Max Error = " << _chart.maxError() << " dB at " << _chart.maxErrorFreq() << " Hz" << std::endl;
          std::cerr << "Max Value = " << _chart.maxValue() << " dB at " << _chart.maxValueFreq() << " Hz" << std::endl;
          std::cerr << "Min Value = " << _chart.minValue() << " dB at " << _chart.minValueFreq() << " Hz" << std::endl;

          std::cerr << "Error_dB Mean = " << _chart.errorStats_dB().mean() << " dB" << std::endl;
          std::cerr << "Error_raw Mean = " << std::scientific << _chart.errorStats_raw().mean() << std::fixed << std::endl;
          std::cerr << "Error_dB Variance = " << _chart.errorStats_dB().variance() << " dB" << std::endl;
          std::cerr << "Error_raw Variance = " << std::scientific << _chart.errorStats_raw().variance() << std::fixed << std::endl;


          if (_sweepCount > 0 && _sweepNum >= _sweepCount) terminate = true;
          if (terminate) return;

          ++ _sweepNum;

          _chart.resetMetrics();
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

  void GetOptions(Int_t *argc, char **argv) override
  {
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

  size_t _sweepNum;

  RootChartProfile<Profile> _chart;
};

//ClassImp(T11ProfileEnv)

int main(int argc, char ** argv)
{
  T11ProfileEnv app(&argc, argv);
  app.GetOptions(&argc, argv);
  app.Run(false);
}
