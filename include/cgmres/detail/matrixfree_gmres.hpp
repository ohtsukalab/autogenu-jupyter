#ifndef CGMRES__MATRIXFREE_GMRES_HPP_
#define CGMRES__MATRIXFREE_GMRES_HPP_

#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "cgmres/types.hpp"

#include "cgmres/detail/macros.hpp"


namespace cgmres {
namespace detail {

template <typename LinearProblem, int _kmax>
class MatrixFreeGMRES {
public:
  static constexpr int dim = LinearProblem::dim;
  static constexpr int kmax = std::min(dim, _kmax);

  MatrixFreeGMRES()
    : hessenberg_mat_(Matrix<kmax+1, kmax+1>::Zero()), 
      basis_mat_(Matrix<dim, kmax+1>::Zero()), 
      b_vec_(Vector<dim>::Zero()), 
      givens_c_vec_(Vector<kmax+1>::Zero()), 
      givens_s_vec_(Vector<kmax+1>::Zero()), 
      g_vec_(Vector<kmax+1>::Zero()) {
    static_assert(dim > 0);
    static_assert(kmax > 0);
    static_assert(dim >= kmax);
  }

  ~MatrixFreeGMRES() = default;

  template <typename... LinearProblemArgs>
  int solve(LinearProblem& linear_problem, 
            LinearProblemArgs... linear_problem_args, 
            Vector<dim>& linear_problem_solution) {
    // Initializes vectors for QR factrization by Givens rotation.
    givens_c_vec_.setZero();
    givens_s_vec_.setZero();
    g_vec_.setZero();
    // Generates the initial basis of the Krylov subspace.
    linear_problem.eval_b(linear_problem_args..., linear_problem_solution, b_vec_);
    g_vec_.coeffRef(0) = b_vec_.template lpNorm<2>();
    basis_mat_.col(0) = b_vec_ / g_vec_.coeff(0);
    // k : the dimension of the Krylov subspace at the current iteration.
    int k = 0;
    for (; k<kmax; ++k) {
      linear_problem.eval_Ax(linear_problem_args..., basis_mat_.col(k), 
                             basis_mat_.col(k+1));
      for (int j=0; j<=k; ++j) {
        hessenberg_mat_.coeffRef(k, j) = basis_mat_.col(k+1).dot(basis_mat_.col(j));
        basis_mat_.col(k+1).noalias() -= hessenberg_mat_.coeff(k, j) * basis_mat_.col(j);
      }
      hessenberg_mat_.coeffRef(k, k+1) = basis_mat_.col(k+1).template lpNorm<2>();
      if (std::abs(hessenberg_mat_.coeff(k, k+1)) < std::numeric_limits<double>::epsilon()) {
        break;
      }
      else {
        basis_mat_.col(k+1).array() /= hessenberg_mat_.coeff(k, k+1);
      }
      // Givens Rotation for QR factrization of the least squares problem.
      for (int j=0; j<k; ++j) {
        givensRotation(hessenberg_mat_.row(k), j);
      }
      const Scalar nu = std::sqrt(hessenberg_mat_.coeff(k, k)*hessenberg_mat_.coeff(k, k)
                                  +hessenberg_mat_.coeff(k, k+1)*hessenberg_mat_.coeff(k, k+1));
      if (nu) {
        givens_c_vec_.coeffRef(k) = hessenberg_mat_.coeff(k, k) / nu;
        givens_s_vec_.coeffRef(k) = - hessenberg_mat_.coeff(k, k+1) / nu;
        hessenberg_mat_.coeffRef(k, k) = givens_c_vec_.coeff(k) * hessenberg_mat_.coeff(k, k) 
                                          - givens_s_vec_.coeff(k) * hessenberg_mat_.coeff(k, k+1);
        hessenberg_mat_.coeffRef(k, k+1) = 0.0;
        givensRotation(g_vec_, k);
      }
      else {
        throw std::runtime_error("Lose orthogonality of the basis of the Krylov subspace");
      }
    }
    // Computes solution_vec by solving hessenberg_mat_ * y = g_vec.
    for (int i=k-1; i>=0; --i) {
      Scalar tmp = g_vec_.coeff(i);
      for (int j=i+1; j<k; ++j) {
        tmp -= hessenberg_mat_.coeff(j, i) * givens_c_vec_.coeff(j);
      }
      givens_c_vec_.coeffRef(i) = tmp / hessenberg_mat_.coeff(i, i);
    }
    for (int i=0; i<dim; ++i) {
      Scalar tmp = 0.0;
      for (int j=0; j<k; ++j) { 
        tmp += basis_mat_.coeff(i, j) * givens_c_vec_.coeff(j);
      }
      linear_problem_solution.coeffRef(i) += tmp;
    }
    return k;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Matrix<kmax+1, kmax+1> hessenberg_mat_;
  Matrix<dim, kmax+1> basis_mat_;
  Vector<dim> b_vec_;
  Vector<kmax+1> givens_c_vec_, givens_s_vec_, g_vec_;

  template <typename VectorType>
  inline void givensRotation(const MatrixBase<VectorType>& column_vec, 
                             const int i_column) const {
    const Scalar tmp1 = givens_c_vec_.coeff(i_column) * column_vec.coeff(i_column) 
                        - givens_s_vec_.coeff(i_column) * column_vec.coeff(i_column+1);
    const Scalar tmp2 = givens_s_vec_.coeff(i_column) * column_vec.coeff(i_column) 
                        + givens_c_vec_.coeff(i_column) * column_vec.coeff(i_column+1);
    CGMRES_EIGEN_CONST_CAST(VectorType, column_vec).coeffRef(i_column) = tmp1;
    CGMRES_EIGEN_CONST_CAST(VectorType, column_vec).coeffRef(i_column+1) = tmp2;
  }

};

} // namespace detail
} // namespace cgmres

#endif // CGMRES__MATRIXFREE_GMRES_HPP_