#ifndef CGMRES__TYPES_HPP_
#define CGMRES__TYPES_HPP_

#include "cgmres/thirdparty/eigen/Eigen/Core"

#define EIGEN_STACK_ALLOCATION_LIMIT 0 // this macro allows unlimited stack memory to avoid OBJECT_ALLOCATED_ON_STACK_IS_TOO_BIG static_assertion error.

namespace cgmres {

///
/// @brief Alias of double.  
///
using Scalar = double;

///
/// @brief Alias of Eigen::Matrix.  
///
template <int rows, int cols>
using Matrix = Eigen::Matrix<Scalar, rows, cols>;

///
/// @brief Alias of Eigen::Vector.  
///
template <int size>
using Vector = Eigen::Matrix<Scalar, size, 1>;

///
/// @brief Alias of Eigen::MatrixBase.  
///
template <class MatrixType>
using MatrixBase = Eigen::MatrixBase<MatrixType>;

///
/// @brief Alias of Eigen::MatrixXd (dynamic-size matrix).
///
using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

///
/// @brief Alias of Eigen::VectorXd (dynamic-size vector).  
///
using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

///
/// @brief Alias of Eigen::VectorXi (dynamic-size integer vector).  
///
using VectorXi = Eigen::Matrix<int, Eigen::Dynamic, 1>;

///
/// @brief Alias of Eigen::Map.  
///
template <class MatrixType>
using Map = Eigen::Map<MatrixType>;

} // namespace cgmres

#endif // CGMRES__TYPES_HPP_