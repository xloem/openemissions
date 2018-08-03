#pragma once

#include <complex>

#include <Eigen/Core>
#include <vector>

using scalar = double;
using cmplx = std::complex<scalar>;

using cvec = Eigen::Matrix<cmplx,Eigen::Dynamic,1>;
using vec = Eigen::Matrix<scalar,Eigen::Dynamic,1>;
using mat = vec;

template <typename V>
inline constexpr typename V::Scalar sum(V const & v) { return v.sum(); }

inline vec abs(cvec const &v) { return v.cwiseAbs(); }

template <typename V, typename V2>
inline constexpr void concat(V & v, V2 const & xtra)
{
  v.resize(xtra.size() + v.size());
  v.tail(xtra.size()) = xtra;
}

template <typename V, typename S>
inline constexpr V convert(std::vector<S> & s)
{
  return Eigen::Map<Eigen::Matrix<S,Eigen::Dynamic,1>>(s.data(), s.size());
}
