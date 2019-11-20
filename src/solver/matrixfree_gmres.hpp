// The matrix-free generalized minimal residual (GMRES) method, which supportes 
// solving the nonlienar problem, used in the C/GMRES method. This program is 
// witten with reference to "C. T. Kelly, Iterative Methods for Linear and 
// Nonlinear Equations, Frontiers in Apllied Mathematics, SIAM (1995)"

#ifndef MATRIXFREE_GMRES_H
#define MATRIXFREE_GMRES_H

#include <iostream>
#include <cmath>
#include <limits>
#include "linear_algebra.hpp"

// Serves the matrix-free GMRES method, which supports for solving the 
// nonlinear problem by using the GMRES method that solves a linear problem 
// Ax = b in a short computational time. This class allocates vectors and 
// matrices used in matrix-free GMRES and computes the solution of GMRES. 
// You have to define MatrixFreeGMRES with template paramter. The first 
// argment of the template paramters is the LinearProblemGenerator which has
// bFunc() and AxFunc(). Second and subsequent arguments are corresponds to
// the argments of bFunc(..., const double* double*) and 
// AxFunc(..., const double*, double*). For example, if you want to solve
// the problem provided by class 'Newton' that has 
// bFunc(const double, const double* const double* const double*)
// and
// AxFunc(const double, const double* const double* const double*),
// you have to define MatrixFreeGMRES<Newton, const double, const double*>.
template <class LinearProblemGenerator, typename... LinearProblemArgs>
class MatrixFreeGMRES {
public:
  // Constructs MatrixFreeGMRES with setting dimension of the solution 
  // and that of the Krylov subspace zero, and sets nullptr for all vectors 
  // and all matrices.
  // MatrixFreeGMRES();
  MatrixFreeGMRES()
    : dim_linear_problem_(0), 
      kmax_(0), 
      hessenberg_mat_(nullptr), 
      basis_mat_(nullptr), 
      b_vec_(nullptr), 
      givens_c_vec_(nullptr), 
      givens_s_vec_(nullptr), 
      g_vec_(nullptr) {
  }

  // Constructs MatrixFreeGMRES with setting dimension of the solution 
  // dim_linear_problem and that of the Krylov subspace as kmax, and allocate 
  // all vectors and all matrices used in the matrix-free GMRES.
  // MatrixFreeGMRES(const int dim_linear_problem, const int kmax);
  MatrixFreeGMRES(const int dim_linear_problem, const int kmax)
    : dim_linear_problem_(dim_linear_problem), 
      kmax_(kmax), 
      hessenberg_mat_(linearalgebra::NewMatrix(kmax+1, kmax+1)), 
      basis_mat_(linearalgebra::NewMatrix(kmax+1, dim_linear_problem)), 
      b_vec_(linearalgebra::NewVector(dim_linear_problem)), 
      givens_c_vec_(linearalgebra::NewVector(kmax+1)), 
      givens_s_vec_(linearalgebra::NewVector(kmax+1)), 
      g_vec_(linearalgebra::NewVector(kmax+1)) {
    if (kmax > dim_linear_problem) {
      kmax_ = dim_linear_problem;
    }
  }

  // Destructs MatrixFreeGMRES with freeing memory of vectors and matrices. 
  // If there are no allocations, this method does not free memory.
  // ~MatrixFreeGMRES();
  ~MatrixFreeGMRES() {
    linearalgebra::DeleteMatrix(hessenberg_mat_);
    linearalgebra::DeleteMatrix(basis_mat_);
    linearalgebra::DeleteVector(b_vec_);
    linearalgebra::DeleteVector(givens_c_vec_);
    linearalgebra::DeleteVector(givens_s_vec_);
    linearalgebra::DeleteVector(g_vec_);
  }

  // Sets dimensions of the solution and that of the Krylov subspace and 
  // reallocates all vectors and alld matrices used in the matrix-free GMRES.
  // void setParameters(const int dim_linear_problem, const int kmax);
  void setParameters(const int dim_linear_problem, const int kmax) {
    linearalgebra::DeleteMatrix(hessenberg_mat_);
    linearalgebra::DeleteMatrix(basis_mat_);
    linearalgebra::DeleteVector(b_vec_);
    linearalgebra::DeleteVector(givens_c_vec_);
    linearalgebra::DeleteVector(givens_s_vec_);
    linearalgebra::DeleteVector(g_vec_);
    dim_linear_problem_ = dim_linear_problem;
    kmax_ = kmax;
    if (kmax > dim_linear_problem) {
      kmax_ = dim_linear_problem;
    }
    hessenberg_mat_ = linearalgebra::NewMatrix(kmax+1, kmax+1);
    basis_mat_ = linearalgebra::NewMatrix(kmax+1, dim_linear_problem);
    b_vec_ = linearalgebra::NewVector(dim_linear_problem);
    givens_c_vec_ = linearalgebra::NewVector(kmax+1);
    givens_s_vec_ = linearalgebra::NewVector(kmax+1);
    g_vec_ = linearalgebra::NewVector(kmax+1);
  }

  // Solves the matrix-free GMRES and generates solution_update_vector, 
  // which is a solution of the matrix-free GMRES.
  void solveLinearProblem(LinearProblemGenerator& linear_problem_generator,
                          LinearProblemArgs... linear_problem_args,
                          double* solution_vec) {
    // Initializes vectors for QR factrization by Givens rotation.
    // Set givens_c_vec_, givens_s_vec_, g_vec_ as zero.
    for (int i=0; i<kmax_+1; ++i) {
      givens_c_vec_[i] = 0;
      givens_s_vec_[i] = 0;
      g_vec_[i] = 0;
    }
    // Generates the initial basis of the Krylov subspace.
    linear_problem_generator.bFunc(linear_problem_args..., solution_vec, 
                                   b_vec_);
    g_vec_[0] = std::sqrt(linearalgebra::SquaredNorm(dim_linear_problem_, 
                                                     b_vec_));
    // basis_mat_[0] = b_vec_ / g_vec[0]
    for (int i=0; i<dim_linear_problem_; ++i) {
      basis_mat_[0][i] = b_vec_[i] / g_vec_[0];
    }
    // k : the dimension of the Krylov subspace at the current iteration.
    int k;
    for (k=0; k<kmax_; ++k) {
      linear_problem_generator.AxFunc(linear_problem_args..., basis_mat_[k], 
                                      basis_mat_[k+1]);
      for (int j=0; j<=k; ++j) {
        hessenberg_mat_[k][j] = linearalgebra::InnerProduct(
            dim_linear_problem_, basis_mat_[k+1], basis_mat_[j]);
        // basis_mat_[k+1] -= hessenberg_mat_[k][j] * basis_mat_[j];
        for (int i=0; i<dim_linear_problem_; ++i) {
          basis_mat_[k+1][i] -= hessenberg_mat_[k][j] * basis_mat_[j][i];
        }
      }
      hessenberg_mat_[k][k+1] = std::sqrt(linearalgebra::SquaredNorm(
          dim_linear_problem_, basis_mat_[k+1]));
      if (std::abs(hessenberg_mat_[k][k+1]) 
          < std::numeric_limits<double>::epsilon()) {
        std::cout << "The modified Gram-Schmidt breakdown at k = " << k 
                  << std::endl;
        break;
      }
      else {
        // basis_mat_[k+1] = basis_mat_[k+1] / hessenberg_mat_[k][k+1];
        for (int i=0; i<dim_linear_problem_; ++i) {
          basis_mat_[k+1][i] = basis_mat_[k+1][i] / hessenberg_mat_[k][k+1];
        }
      }
      // Givens Rotation for QR factrization of the least squares problem.
      for (int j=0; j<k; ++j) {
        givensRotation(hessenberg_mat_[k], j);
      }
      double nu = std::sqrt(hessenberg_mat_[k][k]*hessenberg_mat_[k][k]
                            +hessenberg_mat_[k][k+1]*hessenberg_mat_[k][k+1]);
      if (nu) {
        givens_c_vec_[k] = hessenberg_mat_[k][k] / nu;
        givens_s_vec_[k] = - hessenberg_mat_[k][k+1] / nu;
        hessenberg_mat_[k][k] = givens_c_vec_[k] * hessenberg_mat_[k][k] 
                                - givens_s_vec_[k] * hessenberg_mat_[k][k+1];
        hessenberg_mat_[k][k+1] = 0;
        givensRotation(g_vec_,k);
      }
      else {
        std::cout << "Lose orthogonality of the basis of the Krylov subspace" 
                  << std::endl;
      }
    }
    // Computes solution_vec by solving hessenberg_mat_ * y = g_vec.
    for (int i=k-1; i>=0; --i) {
      double tmp = g_vec_[i];
      for (int j=i+1; j<k; ++j) {
        tmp -= hessenberg_mat_[j][i] * givens_c_vec_[j];
      }
      givens_c_vec_[i] = tmp / hessenberg_mat_[i][i];
    }
    for (int i=0; i<dim_linear_problem_; ++i) {
      double tmp = 0;
      for (int j=0; j<k; ++j) { 
        tmp += basis_mat_[j][i] * givens_c_vec_[j];
      }
      solution_vec[i] += tmp;
    }
  }

  // Prohibits copy constructors.
  MatrixFreeGMRES(const MatrixFreeGMRES&) = delete;
  MatrixFreeGMRES& operator=(const MatrixFreeGMRES&) = delete;

private:
  int dim_linear_problem_, kmax_;
  double **hessenberg_mat_, **basis_mat_;
  double *b_vec_, *givens_c_vec_, *givens_s_vec_, *g_vec_;

  // Applies the Givens rotation for i_column element and i_column+1 
  // element of column_vec, which is a column vector of a matrix.
  // inline void givensRotation(double* column_vec, const int i_column);
  inline void givensRotation(double* column_vec, const int i_column) {
  double tmp1 = givens_c_vec_[i_column] * column_vec[i_column] 
                - givens_s_vec_[i_column] * column_vec[i_column+1];
  double tmp2 = givens_s_vec_[i_column] * column_vec[i_column] 
                + givens_c_vec_[i_column] * column_vec[i_column+1];
  column_vec[i_column] = tmp1;
  column_vec[i_column+1] = tmp2;
}
};

#endif // MATRIXFREE_GMRES_H