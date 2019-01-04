#pragma once

#include <Eigen/Core>

// TODO: provide for writing to the file
// TODO: merge with file.hpp; the classes can basically be just combined directly I think

// stores raw samples inline
// virtual functions for subclasses to implement specific formats
template <typename Scalar>
class RawLogFile
{
public:
  RawLogFile(std::uniqe_ptr<std::iostream> && stream, bool bigEndian = false)
  : _stream(std::move(stream)),
    _bigEndian(bigEndian)
  { }

  template <typename Derived>
  void readMany(Eigen::MatrixBase<Derived> const &vec)
  {
    _readBuf.conservativeResize(_vec.size());
    _stream.read(reinterpret_cast<char*>(&_readBuf[0]), _readBuf.size() * sizeof(Scalar));
    _readBuf.conservativeResize(_stream.gcount() / sizeof(Scalar));
    if (nativeBigEndian() != _bigEndian) {
      // swap bytes
      // TODO vectorize ? could use a map and reorder columns
      for (size_t i = 0; i < _readBuf.size(); ++ i) {
        Scalar temp = _readBuf[i];
        uint8_t * outPtr = reinterpret_cast<uint8_t*>(&_readBuf[i]);
        uint8_t * inPtr = reinterpret_cast<uint8_t*>(&temp);
        for (size_t j = 0; j < sizeof(Scalar); ++ j) {
          outPtr[j] = inPtr[sizeof(Scalar) - 1 - j];
        }
      }
    }
    vec.derived().conservativeResize(_readBuf.size());
    castFromTo(_readBuf, vec);
  }

  static constexpr ::Scalar epsilon()
  {
    using basic_t = complex_removed_t<Scalar>;
    if (real_limits<Scalar>::is_integer) {
      return ::Scalar(2) / (real_limits<Scalar>::max - real_limits<Scalar>::lowest);
    } else {
      throw std::logic_error("file type is floating point");
    }
  }

private:
  std::unique_ptr<std::iostream> _stream;
  bool _bigEndian;

  HeapVector<Scalar> _readBuf;
};

