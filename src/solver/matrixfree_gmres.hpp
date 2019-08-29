// The matrix-free generalized minimal residual (GMRES) method, which supportes 
// solving the nonlienar problem, used in the C/GMRES method. This program is 
// witten with reference to "C. T. Kelly, Iterative Methods for Linear and 
// Nonlinear Equations, Frontiers in Apllied Mathematics, SIAM (1995)"

#ifndef MATRIXFREE_GMRES_H
#define MATRIXFREE_GMRES_H

#include <iostream>
#include <limits>
#include "nmpc_model.hpp"
#include "linear_algebra.hpp"

// Serves the matrix-free GMRES method, which supports for solving the 
// nonlinear problem by using the GMRES method that solves a linear problem 
// Ax = b in a short computational time. This class allocates vectors and 
// matrices used in matrix-free GMRES and computes the solution of GMRES. 
// You have to inherit this class and override bFunc() and axFunc() which 
// generates vectors corresponding to b and Ax in Ax=b, respectively.
class MatrixFreeGMRES {
public:
  // Constructs MatrixFreeGMRES with setting dimension of the solution 
  // and that of the Krylov subspace zero, and sets nullptr for all vectors 
  // and all matrices.
  MatrixFreeGMRES();

  // Constructs MatrixFreeGMRES with setting dimension of the solution 
  // dim_solution and that of the Krylov subspace as kmax, and allocate 
  // all vectors and all matrices used in the matrix-free GMRES.
  MatrixFreeGMRES(const int dim_solution, const int kmax);

  // Destructs MatrixFreeGMRES with freeing memory of vectors and matrices. 
  // If there are no allocations, this method does not free memory.
  ~MatrixFreeGMRES();

  // Sets dimensions of the solution and that of the Krylov subspace and 
  // reallocates all vectors and alld matrices used in the matrix-free GMRES.
  void setGMRESParameters(const int dim_solution, const int kmax);

  // Solves the matrix-free GMRES and generates solution_update_vector, 
  // which is a solution of the matrix-free GMRES.
  void solveGMRES(const double time, const double* state_vec, 
                  const double* current_solution_vec, double* solution_update_vec);

private:
  bool allocation_flag_;
  int dim_solution_, kmax_;
  double **hessenberg_mat_, **basis_mat_;
  double *b_vec_, *givens_c_vec_, *givens_s_vec_, *g_vec_;

  // Applies the Givens rotation for i_column element and i_column+1 
  // element of column_vec, which is a column vector of a matrix.
  inline void givensRotation(double* column_vec, const int i_column);

  // Generates a vector corresponding to b in Ax=b.
  virtual void bFunc(const double time, const double* state_vec, 
                      const double* current_solution_vec, 
                      double* equation_error_vec) = 0;

  // Generates a vector corresponding to Ax in Ax=b.  
  virtual void axFunc(const double time, const double* state_vec, 
                      const double* current_solution_vec, 
                      const double* direction_vec, double* ax_vec) = 0;
};

#endif // MATRIXFREE_GMRES_H