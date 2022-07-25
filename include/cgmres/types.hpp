#ifndef TYPES_HPP_
#define TYPES_HPP_

#include "thirdparty/eigen/Eigen/Core"

namespace cgmres {

using Scalar = double;

template <int rows, int cols>
using Matrix = Eigen::Matrix<Scalar, rows, cols>;

template <int size>
using Vector = Eigen::Matrix<Scalar, size, 1>;

template <class MatrixType>
using MatrixBase = Eigen::MatrixBase<MatrixType>;

} // namespace cgmres

#endif // TYPES_HPP_