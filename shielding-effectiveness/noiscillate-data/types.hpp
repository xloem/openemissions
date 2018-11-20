#pragma once

#include <Eigen/Core>

// TODO: higher precision could be attained by using [u]int64_t and changing code to never leave (-1.0, 1.0)
// this would prevent the error that results when accumulating averages over a long duration
//      perhaps measure how long this duration would be before deciding to do this?
using Scalar = double;
using Complex = std::complex<Scalar>;

template <typename Type, size_t count, int options = 0>
using HeapVectors = Eigen::Matrix<Type, Eigen::Dynamic, count, options>;
template <typename Type, int options = 0>
using HeapVector = HeapVectors<Type, 1, options>;
template <typename Type, size_t maxsize, size_t count, int options = 0>
using StackVectors = Eigen::Matrix<Type, Eigen::Dynamic, count, options, maxsize>;
template <typename Type, size_t maxsize, int options = 0>
using StackVector = StackVectors<Type, maxsize, 1, options>;

// allow casting of all complex types in eigen
// TODO: submit patch to eigen? alternatively make specializations of std::complex?
namespace Eigen { namespace internal {
  template <typename OldType, typename NewType>
  struct cast_impl<std::complex<OldType>, std::complex<NewType>>
  {
    EIGEN_DEVICE_FUNC
    static inline std::complex<NewType> run(const std::complex<OldType>& x)
    {
      return {static_cast<NewType>(x.real()), static_cast<NewType>(x.image())};
    }
  };
} }

template <typename T> struct real_limits : public std::numeric_limits<T> {}; 
template <typename T> struct real_limits<std::complex<T>> : public std::numeric_limits<T> {};

template <typename Old, typename New>
void castFromTo(Eigen::MatrixBase<Old> const & src, Eigen::MatrixBase<New> const & dst)
{
  auto ret = const_cast<Eigen::MatrixBase<New> &>(dst);

  if (! real_limits<New>::is_integer) {
    // new is floating
    if (! real_limits<Old>::is_integer) {
      // new and old are both floating
      ret = src.template cast<New>();
    } else {
      // new is floating, old is integral
      ret = src.template cast<New>() / (static_cast<New>(real_limits<Old>::max - real_limits<Old>::lowest) / 2) -
        (real_limits<Old>::is_signed ? 0 : 1.0);
    }
  } else {
    // new is integral
    if ( real_limits<Old>::is_integer) {
      // new and old are both integral
      // TODO
      throw std::logic_error("unimplemented int-int vector cast");
    } else {
      // new is integral, old is floating
      ret = ((src + (real_limits<New>::is_signed ? 1.0 : 0)) * (real_limits<New>::max - real_limits<New>::lowest) / 2)
          .template cast<New>();
    }
  }
}

inline bool nativeBigEndian()
{
  union test_t {
    uint32_t num;
    uint8_t bytes[4];
  };
  return test_t({0x01234567}).bytes[0] == 0x01;
}
