#pragma once

#include "types.hpp"
#include "recbufmeta.hpp"

#include <memory>

#include <SoapySDR/Device.hpp>
#include <SoapySDR/Errors.hpp>
#include <SoapySDR/Version.hpp>

#include <TSystem.h>
#include <TUUID.h>

class SoapyLive
{
public:
  SoapyLive(const char * deviceStr, Scalar dB, Scalar sampleRate, Scalar settleSecs = 0.0)
  : _device(SoapySDR::Device::make(deviceStr)),
    _settleSecs(settleSecs),
    _skip(0),
    _sampleTime(0)
  {
    _device->setGainMode(SOAPY_SDR_RX, 0, false);
    
    _device->setSampleRate(SOAPY_SDR_RX, 0, sampleRate);
    _device->setGain(SOAPY_SDR_RX, 0, dB);

    // TODO: hardcoded to CS8, but generic code should go in e.g. emap or freesdr
    std::cerr << "WARNING: hardcoded to 8 bits for RTL-SDR." << std::endl;
    std::cerr << "         precision will be lost from better radios." << std::endl;
    _stream = _device->setupStream(SOAPY_SDR_RX, "CS8", {0}/*, {{"bufflen", "1048576"}}*/);

    if (_device->activateStream(_stream) != 0)
    {
      throw std::runtime_error("failed to activate stream");
    }

    _freq = _device->getFrequency(SOAPY_SDR_RX, 0);
    _rate = _device->getSampleRate(SOAPY_SDR_RX, 0);
    _gain = _device->getGain(SOAPY_SDR_RX, 0);
  }

  std::string ident()
  {
    std::string ret = "driver=" + _device->getDriverKey() + ",hardware=" + _device->getHardwareKey();
    SoapySDR::Kwargs info = SoapySDR::Device::enumerate(_device->getHardwareInfo())[0];
    info["driver"] = _device->getDriverKey();
    info["hardware"] = _device->getHardwareKey();
    SoapySDR::Kwargs extra;
    extra["samplerate"] = std::to_string(_rate);
    extra["gain"] = std::to_string(_gain);
    extra["settle"] = std::to_string(_settleSecs);
    extra["soapysdrapi"] = SoapySDR::getAPIVersion();
    extra["soapysdrlib"] = SoapySDR::getLibVersion();
    SysInfo_t sysinfo;
    if (gSystem->GetSysInfo(&sysinfo) == 0)
    {
      extra["syscpu"] = sysinfo.fCpuType;
      extra["sysmodel"] = sysinfo.fModel;
      extra["sysos"] = sysinfo.fOS;
    }
    return SoapySDR::KwargsToString(info) + "\n" + SoapySDR::KwargsToString(extra);
  }

  Scalar tune(Scalar freq)
  {
    _device->setFrequency(SOAPY_SDR_RX, 0, freq);
    _skip = _settleSecs * _rate;
    _freq = _device->getFrequency(SOAPY_SDR_RX, 0);
    return _freq;
  }

  size_t recommendedMTU() const
  {
    return _device->getStreamMTU(_stream);
  }

  template <typename Derived>
  void readMany(Eigen::MatrixBase<Derived> const & _vec, RecBufMeta & meta)
  {
    Eigen::MatrixBase<Derived> & vec = const_cast<Eigen::MatrixBase<Derived> &>(_vec);
    int flags;
    long long timeNs;

    meta.rate = _rate;
    meta.freq = _freq;

    while (_skip > 0)
    {
      _buf.resize(2 * _skip);
      void * buf = _buf.data();
      auto result = _device->readStream(_stream, &buf, _skip, flags, timeNs, 2000000);
      /*if (result == SOAPY_SDR_TIMEOUT)
      {
        std::cerr << "device timeout" << std::endl;
        vec.derived().conservativeResize(0);
        meta.sampleTime = _sampleTime;
        return;
      }
      else */if (result < 0)
      {
        throw std::runtime_error(std::string("SoapySDR: ") + SoapySDR::errToStr(result));
      }
      _sampleTime += result;
      _skip -= result;
    }

    _buf.resize(2 * _vec.size());
    void * buf = _buf.data();

    meta.sampleTime = _sampleTime;

    // TODO: use global cast function to skip some of this extra stuff
    int count = _device->readStream(_stream, &buf, _vec.size(), flags, timeNs, 10000000);
    if (count < 0) throw std::runtime_error(std::string("stream read error: ") + SoapySDR::errToStr(count));
    Eigen::Map<Eigen::Array<int8_t, 2, Eigen::Dynamic>> eigenRawData(&_buf[0], 2, count);
    auto eigenScalarData = (eigenRawData.cast<Scalar>() + 0.5) / 127.5;
    vec.derived().conservativeResize(eigenScalarData.cols());
    vec.real() = eigenScalarData.row(0);
    vec.imag() = eigenScalarData.row(1);
    _sampleTime += vec.size();
  }

  static constexpr Scalar epsilon() {
    return 1 / 127.5;
  }

  static constexpr Scalar midval() {
    return 0.5 / 127.5;
  }

  ~SoapyLive()
  {
    _device->deactivateStream(_stream);
    _device->closeStream(_stream);
    SoapySDR::Device::unmake(_device); _device = nullptr;
  }

private:
  SoapySDR::Device * _device;
  SoapySDR::Stream * _stream;

  Scalar _settleSecs;
  Scalar _rate;
  Scalar _gain;
  Scalar _freq;
  uint64_t _sampleTime;
  uint64_t _skip;

  std::vector<int8_t> _buf;
};

/*
 * How am I going to manage discarding samples after tuning?
 * I can calculate how many samples to discard, but which piece of code will manage it?
 *
 * hmmm
 *
 * well, I could have the main function do it.  then it would need to let the processor or the
 * data feeder know how many it skipped.  that's not very reusable.
 * I guess the data feeder might do it: high reuse
 * I could also put it in the sample source: then it would need to report the sample time to let
 * the user know how many samples were lost sometimes, to predict timing
 *
 * sounds more conceptually relevent than the data feeder
 *
 * yeah i guess so
 */
