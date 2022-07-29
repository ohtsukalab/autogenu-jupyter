#ifndef CGMRES__TYPES_HPP_
#define CGMRES__TYPES_HPP_

#include "cgmres/thirdparty/eigen/Eigen/Core"

namespace cgmres {

using Scalar = double;

template <int rows, int cols>
using Matrix = Eigen::Matrix<Scalar, rows, cols>;

template <int size>
using Vector = Eigen::Matrix<Scalar, size, 1>;

template <class MatrixType>
using MatrixBase = Eigen::MatrixBase<MatrixType>;

using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

using VectorXi = Eigen::Matrix<int, Eigen::Dynamic, 1>;

template <class MatrixType>
using Map = Eigen::Map<MatrixType>;

} // namespace cgmres

#endif // CGMRES__TYPES_HPP_