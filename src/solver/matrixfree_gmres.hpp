//
// The matrix-free generalized minimal residual (GMRES) method, a numerical solver for nonlienar programming problem, used in the C/GMRES method.
// This program is witten with reference to "C. T. Kelly, Iterative Methods for Linear and Nonlinear Equations, Frontiers in Apllied Mathematics, SIAM (1995)"
//

#ifndef MATRIXFREE_GMRES_H
#define MATRIXFREE_GMRES_H

#include <iostream>
#include "nmpc_model.hpp"
#include "linear_funcs.hpp"

// Serves the matrix-free GMRES method, which supports solving the nonlinear problem by using the GMRES method that solves a linear problem Ax = b in a short computational time. 
// This class allocates vectors and matrices used in matrix-free GMRES and computes the solution of GMRES. 
// You have to inherit this class and override bFunc() and axFunc() which generates vectors corresponding to b and Ax in Ax=b, respectively.
class MatrixFreeGMRES {
public:
    // Sets dim_equation_ and kmax_ zero, and all vectors and matrices having nullptr.
    MatrixFreeGMRES();

    // Sets dimensions of the solution vector and the Krylov subspace and allocate vectors and matrices used in the matrix-free GMRES.
    MatrixFreeGMRES(const int dim_equation, const int kmax);

    // Free memory of vectors and matrices.
    ~MatrixFreeGMRES();

    // Sets dimensions of the solution vector and the Krylov subspace and reallocate vectors and matrices used in the matrix-free GMRES.
    void setGMRESParams(const int dim_equation, const int kmax);

    // Solves the GMRES and generates solution_update_vector, which is a solution of the matrix-free GMRES.
    void forwardDifferenceGMRES(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec);


private:
    bool allocation_flag_;
    int dim_equation_, kmax_;
    double **hessenberg_mat_, **basis_mat_;
    double *b_vec_, *givens_c_vec_, *givens_s_vec_, *g_vec_;

    // Applies the Givens rotation for i_column element and i_column+1 element of column_vec, which is a column vector of a matrix.
    inline void givensRotation(double* column_vec, const int i_column);

    // Generates a vector corresponding to b in Ax=b.
    virtual void bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* equation_error_vec) = 0;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    virtual void axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec) = 0;
};

#endif