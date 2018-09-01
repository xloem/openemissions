#pragma once

#include <vector>

class RtlSdrIQDump
{
public:
  RtlSdrIQDump(std::istream & stream)
  : stream(stream)
  { }
  
  Complex readOne()
  {
    std::complex<uint8_t> val;
    stream.read(reinterpret_cast<char*>(&val), 2);
    if (stream.eof())
      throw std::runtime_error("EOF");
    return (Complex(val.real(), val.imag()) - Complex(127.5, 127.5)) / 127.5;
  }

  template <typename Derived>
  void readMany(Eigen::MatrixBase<Derived> const & _vec)
  {
    std::vector<uint8_t> val(2 * _vec.size());
    stream.read(reinterpret_cast<char*>(&val[0]), val.size());
    Eigen::Map<Eigen::Array<uint8_t, 2, Eigen::Dynamic>> eigenRawData(&val[0], 2, stream.gcount() / 2);
    auto eigenScalarData = eigenRawData.cast<Scalar>() / 127.5 - 1.0;
    Eigen::MatrixBase<Derived> & vec = const_cast<Eigen::MatrixBase<Derived> &>(_vec);
    vec.derived().resize(eigenScalarData.cols());
    vec.real() = eigenScalarData.row(0);
    vec.imag() = eigenScalarData.row(1);
  }

  static constexpr Scalar epsilon() { return 1 / 127.5; }

private:
  std::istream & stream;
};
