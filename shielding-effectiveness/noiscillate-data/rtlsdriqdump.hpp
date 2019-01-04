#pragma once

#include "types.hpp"
#include "recbufmeta.hpp"

class RtlSdrIQDump
{
public:
  RtlSdrIQDump(std::istream & stream)
  : _stream(stream),
    _sampleTime(0),
    _freq(0)
  { }

  Scalar tune(Scalar freq)
  {
    return _freq = freq;
  }
  
  template <typename Derived>
  void readMany(Eigen::MatrixBase<Derived> const & _vec, RecBufMeta & meta)
  {
    std::vector<uint8_t> val(2 * _vec.size());
    _stream.read(reinterpret_cast<char*>(&val[0]), val.size());
    Eigen::Map<Eigen::Array<uint8_t, 2, Eigen::Dynamic>> eigenRawData(&val[0], 2, _stream.gcount() / 2);
    auto eigenScalarData = eigenRawData.cast<Scalar>() / 127.5 - 1.0;
    Eigen::MatrixBase<Derived> & vec = const_cast<Eigen::MatrixBase<Derived> &>(_vec);
    vec.derived().resize(eigenScalarData.cols());
    vec.real() = eigenScalarData.row(0);
    vec.imag() = eigenScalarData.row(1);
    meta.rate = 0;
    meta.freq = _freq;
    meta.sampleTime = _sampleTime;
    _sampleTime += eigenScalarData.cols();
  }

  static constexpr Scalar epsilon() { return 1 / 127.5; }

private:
  std::istream & _stream;
  uint64_t _sampleTime;
  Scalar _freq;
};

// TODO NEXT: make a generic class for reading and writing raw data files, I guess
