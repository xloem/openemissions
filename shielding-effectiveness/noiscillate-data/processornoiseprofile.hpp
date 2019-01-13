#pragma once

#include <fstream>
#include <set>
#include <vector>

#include <time.h>

#include "recbufmeta.hpp"
#include "stats.hpp"

template <typename _Scalar, typename _StatsAccumulator, StatsStatistic STATISTIC = STATS_MEAN>
class ProcessorNoiseProfile
{
  static int64_t constexpr INVALIDINT64 = -0x8000000000000000;
public:
  using StatsAccumulator = _StatsAccumulator;
  using Scalar = _Scalar;

  ProcessorNoiseProfile(std::string fname, StatsAccumulator const & initial, int64_t env = INVALIDINT64, std::vector<int64_t> emits = {}, std::vector<std::string> srcDesc = {})
  : _initial(initial),
    _totalSamps(0),
    _fname(fname),
    _startTime(time(0)),
    _env(env),
    _emits(emits),
    _srcDescs(srcDesc)
  {
    std::ifstream f(fname);
    f.ignore(1);
    if (f.good())
    {
      f.seekg(0, std::ios_base::beg);
      uint8_t ver;
      f.read((char*)&ver, sizeof(ver));
      switch(ver)
      {
      case 0:
        _startTime = 0;
        break;
      case 1:
        {
          f.read((char*)&_env, sizeof(_env));
          if (env != INVALIDINT64 && _env != env) throw std::invalid_argument(std::string("File is for environment #") + std::to_string(_env));
          uint64_t emitterCount;
          f.read((char*)&emitterCount, sizeof(emitterCount));
          _emits.resize(emitterCount);
          if (emits.size() && _emits.size() != emits.size()) throw std::invalid_argument("Emitter mismatch in file.");
          for (size_t i = 0; i < _emits.size(); ++ i)
          {
            f.read((char*)&_emits[i], sizeof(_emits[i]));
            if (emits.size() && _emits[i] != emits[i]) throw std::invalid_argument("Emitter mismatch in file.");
          }
          f.read((char*)&_startTime, sizeof(_startTime));
          f.read((char*)&_stopTime, sizeof(_stopTime));
          uint64_t srcDescCount;
          f.read((char*)&srcDescCount, sizeof(srcDescCount));
          auto it = _srcDescs.begin();
          for (size_t i = 0; i < srcDescCount; ++ i)
          {
            uint64_t len;
            f.read((char*)&len, sizeof(len));
            std::string desc;
            desc.resize(len);
            f.read(&desc[0], len);
            if (it == _srcDescs.end())
            {
              it = _srcDescs.insert(it, std::move(desc));
              ++ it;
            }
            else if (desc != *it)
            {
              throw std::invalid_argument("File source description contains differing '" + desc + "'");
            }
            else
            {
              ++ it;
            }
          }
        }
        break;
      default:
        throw std::invalid_argument("unsupported file version");
      }
      uint64_t count;
      f.read((char*)&count, sizeof(count));
      for (uint64_t i = 0; i < count; ++ i)
      {
        Scalar freq;
        f.read((char*)&freq, sizeof(freq));
        _binsByFrequency.emplace(freq, f);
      }
    }
  }

  void write()
  {
    std::string tmpname = _fname + "_";
    std::ofstream f(tmpname);
    uint8_t ver = 1;
    uint64_t count;
    _stopTime = time(0);
    f.write((char*)&ver, sizeof(ver));
    f.write((char*)&_env, sizeof(_env));
    count = _emits.size();
    f.write((char*)&count, sizeof(count));
    for (auto & emit : _emits)
    {
      f.write((char*)&emit, sizeof(emit));
    }
    f.write((char*)&_startTime, sizeof(_startTime));
    f.write((char*)&_stopTime, sizeof(_stopTime));
    count = _srcDescs.size();
    f.write((char*)&count, sizeof(count));
    for (auto & desc : _srcDescs)
    {
      count = desc.size();
      f.write((char*)&count, sizeof(count));
      f.write(&desc[0], count);
    }
    count = _binsByFrequency.size();
    f.write((char*)&count, sizeof(count));
    for (auto & bin : _binsByFrequency)
    {
      Scalar freq = bin.first;
      f.write((char*)&freq, sizeof(freq));
      bin.second.write(f);
    }
    f.close();
    std::rename(tmpname.c_str(), _fname.c_str());
  }

  template <typename Predicate, typename ... Profiles>
  void fromStatsPredicate(Predicate predicate, Profiles & ... profiles)
  {
    _binsByFrequency.clear();
    std::set<Scalar> sharedFreqs;
    bool initial = true;
    for (auto * profile : {&profiles ...})
    {
      std::set<Scalar> freqs;
      for (auto pair : profile->_binsByFrequency)
      {
        auto freq = pair.first;
        if (initial)
        {
          sharedFreqs.insert(freq);
        }
        else
        {
          if (sharedFreqs.count(freq))
          {
            freqs.insert(freq);
          }
        }
      }
      if (initial)
      {
        _totalSamps = profile->_totalSamps;
        _startTime = profile->_startTime;
        _stopTime = profile->_stopTime;
        initial = false;
      }
      else
      {
        if (_totalSamps > profile->_totalSamps)
        {
          _totalSamps = profile->_totalSamps;
        }
        if (_startTime > profile->_startTime)
        {
          _startTime = profile->_startTime;
        }
        if (_stopTime < profile->_stopTime)
        {
          _stopTime = profile->_stopTime;
        }
        sharedFreqs = std::move(freqs);
      }
    }
    if (!sharedFreqs.size())
    {
      throw std::invalid_argument("profiles share no frequencieds");
    }

    for (auto freq : sharedFreqs)
    {
      _binsByFrequency.emplace(freq, predicate(profiles._binsByFrequency[freq]...));
    }
  }

  size_t bufferSizeMin() const
  {
    return 1;
  }

  size_t bufferSizeMultiple() const
  {
    return 1;
  }

  size_t bufferSizeMax() const
  {
    return std::numeric_limits<size_t>::max();
  }

  template <class DataFeeder, class Derived>
  void process(DataFeeder & feeder, Eigen::DenseBase<Derived> const & chunk, RecBufMeta const & meta)
  {
    _binsByFrequency.emplace(meta.freq, _initial).first->second.add(chunk);
    _totalSamps += chunk.size();
  }

  uint64_t periodsConsidered() const
  {
    return _totalSamps;
  }

  Scalar bestPeriod() const
  {
    return 1;
  }

  Scalar bestFrequency() const
  {
    return 1;
  }

  Scalar bestMag(Scalar & freq)
  {
    auto ret = _binsByFrequency.at(freq).template get<STATISTIC>();
    //std::cerr << freq << " -> " << std::scientific << ret << std::fixed << std::endl;
    return ret;
  }

  Scalar bestMagVariance(Scalar & freq)
  {
    auto & bin = _binsByFrequency.at(freq);
    auto ret = StatsDistributionSampling<Scalar, STATISTIC>(bin.fakeInfinitePopulation(), bin.size()).variance();
    //std::cerr << freq << " x " << bin.size() << " -> " << std::scientific << ret << std::fixed << std::endl;
    return ret;
  }

  decltype(auto) distribution(Scalar & freq)
  {
    auto & bin = _binsByFrequency.at(freq);
    return StatsDistributionSampling<Scalar, STATISTIC>(bin.fakeInfinitePopulation(), bin.size());
  }

  time_t startTime() const { return _startTime; }
  time_t stopTime() const { return _stopTime; }
  int64_t environment() const { return _env; }
  std::vector<int64_t> const & emitters() const { return _emits; }
  std::vector<std::string> const & srcDescs() { return _srcDescs; }

  class FrequencyIterator
  {
    using MapIterator = typename std::map<Scalar, StatsAccumulator>::const_iterator;
  public:
    FrequencyIterator(MapIterator && it)
    : _wrapped(std::move(it))
    { }

    FrequencyIterator & operator--() { -- _wrapped; return *this; }
    FrequencyIterator & operator--(int) { auto ret(_wrapped); -- _wrapped; return ret; }
    FrequencyIterator & operator++() { ++ _wrapped; return *this; }
    FrequencyIterator & operator++(int) { auto ret(_wrapped); ++ _wrapped; return ret; }
    Scalar operator*() { return (*_wrapped).first; }
    bool operator!=(FrequencyIterator const & it) { return _wrapped != it._wrapped; }
    bool operator==(FrequencyIterator const & it) { return _wrapped == it._wrapped; }
  private:
    MapIterator _wrapped;
  };

  FrequencyIterator begin() { return {_binsByFrequency.cbegin()}; }
  FrequencyIterator end() { return {_binsByFrequency.cend()}; }
  size_t size() { return _binsByFrequency.size(); }

private:
  StatsAccumulator _initial;
  std::map<Scalar, StatsAccumulator> _binsByFrequency;
  uint64_t _totalSamps;
  std::string _fname;
  uint64_t _startTime; // UTC unix timestamp
  uint64_t _stopTime;
  int64_t _env;
  std::vector<int64_t> _emits;
  std::vector<std::string> _srcDescs;
};
