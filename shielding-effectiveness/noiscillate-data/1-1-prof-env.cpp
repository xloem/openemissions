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
#include <TGraphAsymmErrors.h>
#include <TSysEvtHandler.h>

#include "datafeeder.hpp"
#include "processornoiseprofile.hpp"
#include "soapylive.hpp"

using Scalar = double;

class T11ProfileEnv : public TApplication, TSignalHandler
{
public:
  T11ProfileEnv(Int_t * argc, char ** argv)
  : TApplication("T11ProfileEnv", argc, argv),
    TSignalHandler(ESignals(kSigInterrupt | kSigQuit | kSigTermination | kSigUser1)),
    _randGen(_randDev()),
    _terminate(false)
  {
    for (Scalar freq = _minTuneFreq; freq <= _maxTuneFreq; freq += _hopDistance)
    {
      _TuningFreq obj;
      obj.freq = freq;
      obj.idx = _tuneFreqs.size();
      _tuneFreqs.push_back(obj);
    }
    _graph.Set(_tuneFreqs.size());
    for (size_t i = 0; i < _tuneFreqs.size(); ++ i)
    {
      _graph.SetPoint(i, _tuneFreqs[i].freq, (Double_t)NAN);
      _graph.SetPointError(i, -_sampleRate / 2, _sampleRate / 2, (Double_t)NAN, (Double_t)NAN);
    }

    _graph.SetTitle("Recording Profile");
    _graph.GetXaxis()->SetTitle("Frequency (Hz)");
    _graph.GetYaxis()->SetTitle("Power (dB)");
  }

  void Run(Bool_t retrn)
  {
    _terminate = false;
    SetSignalHandler(this);
    std::cerr << "Device construction args: '" << _soapyArgs << "'." << std::endl;
    std::cerr << "Recording at " << _sampleRate << " Hz in hops of " << _hopDistance << " Hz from " << _minTuneFreq << " Hz through " << _maxTuneFreq << " Hz." << std::endl;

    SoapyLive data(_soapyArgs.c_str(), _gain, _sampleRate);
    ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE> processor(_fname, {data.epsilon()});

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
    data.tune(curFreq_it->freq);

    for (data.readMany(buffer, bufMeta); buffer.size(); data.readMany(buffer, bufMeta))
    {
      auto preprocessed = buffer.array().real().eval();
      feeder.add(preprocessed, bufMeta);

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
      data.tune(curFreq_it->freq);

      // process chart
      // TODO: use a 2D histogram? user may understand signal better
      auto rawY = processor.bestMag(freq);
      auto rawVar = processor.bestMagVariance(freq);
      auto rawErr = sqrt(rawVar) * 3;
      StackVector<Double_t,3> unproc, proc;
      unproc[0] = 0;
      unproc[1] = -rawErr;
      unproc[2] = rawErr;
      proc = (unproc.array() + rawY).log10() * 20;
      proc.tail(2).array() -= proc[0];
      _graph.GetY()[freqIdx] = proc[0];
      _graph.GetEYlow()[freqIdx] = proc[1];
      _graph.GetEYhigh()[freqIdx] = proc[2];
      _graph.Draw("ALP");

      // write file if a full sweep has completed
      if (curFreq_it == _tuneFreqs.begin() || _terminate)
      {
        std::cerr << "-sto";
        processor.write();
        std::cerr << "re-";
        if (retrn || _terminate) return;
      }
    }
  }

  void GetOptions(Int_t *argc, char **argv)
  {
    bool setEnv = false;
    _gain = 43;
    _sampleRate = 2400000;
    _minTuneFreq =  200000000;
    _maxTuneFreq = 1400000000;
    _hopDistance = 0;
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
        if (!setEnv)
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
      std::cerr << e.what() << std::endl;
      PrintHelp();
    }
    catch(std::out_of_range)
    {
      std::cerr << "Value out of range." << std::endl;
      Terminate(-1);
    }
  }

  // system signal handler
  Bool_t Notify()
  {
    _terminate = true;
  }

private:
  void PrintHelp()
  {
    std::cerr << "\
Records from an SDR and stores a histogram profile of the values received,\
useful for statistical analysis of the noise.\
\
Usage: 1-1-prof-env [options ...] -o a.noisep env [emit1 ... emitN]\
\
            env : An int representing the environment and setup\
          emitX : Ints representing all active emitter setups\
\
Options:\
    -o a.noisep : create or update the file a.noisep\
   -d \"args...\" : SoapySDR device construction args (default: \"\")\
 -f minHz maxHz : hop between the passed range of tuning frequencies in hertz\
                  (default: 2000000000:1400000000)\
        -r rate : samplerate to record at (default: 2400000)\
      -g gaindB : dB gain for hardware amplifier (default: 43)\
       -h hopHz : distance between frequencies to tune to\
                  (default: samplerate / 2)\
\
";
    Terminate(-1);
  }


  struct _TuningFreq
  {
    Scalar freq;
    size_t idx;
  };

  std::string _fname;
  long long _env;
  std::vector<long long> _emits;

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

  TGraphAsymmErrors _graph;
};

//ClassImp(T11ProfileEnv)

int main(int argc, char ** argv)
{
  T11ProfileEnv app(&argc, argv);
  app.Run(false);
}
