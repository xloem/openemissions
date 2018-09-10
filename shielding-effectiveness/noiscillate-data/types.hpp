#pragma once

#include <Eigen/Core>

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
