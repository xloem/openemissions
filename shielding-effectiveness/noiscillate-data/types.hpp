#pragma once

#include <Eigen/Core>

using Scalar = double;
using Complex = std::complex<Scalar>;


template <typename Type, size_t count>
using HeapVectors = Eigen::Matrix<Type, Eigen::Dynamic, count>;
template <typename Type>
using HeapVector = HeapVectors<Type, 1>;
template <typename Type, size_t maxsize, size_t count>
using StackVectors = Eigen::Matrix<Type, Eigen::Dynamic, count, 0, maxsize>;
template <typename Type, size_t maxsize>
using StackVector = StackVectors<Type, maxsize, 1>;

