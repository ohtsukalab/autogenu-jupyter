//
// The matrix-free generalized minimal residual (GMRES) method, a numerical solver for nonlienar programming problem, used in the C/GMRES method.
// This program is witten with reference to "C. T. Kelly, Iterative Methods for Linear and Nonlinear Equations, Frontiers in Apllied Mathematics, SIAM (1995)"
//

#ifndef MATRIXFREE_GMRES_H
#define MATRIXFREE_GMRES_H

#include <pick_model.hpp>
#include <iostream>
#include <eigen3/Eigen/Core>


// Serves the matrix-free GMRES method, which supports solving the nonlinear problem by using the GMRES method that solves a linear problem Ax = b in a short computational time.
// You have to inherit this class and override bFunc() and axFunc() which generates vectors corresponding to b and Ax in Ax=b, respectively.
// This class allocates vectors and matrices used in matrix-free GMRES and computes the solution of GMRES.
class MatrixFreeGMRES {
private:
    int dim_equation_, max_dim_krylov_;
    Eigen::MatrixXd hessenberg_mat_, basis_mat_;
    Eigen::VectorXd b_vec_, givens_c_vec_, givens_s_vec_, g_vec_;

    // Applies the Givens rotation for i_column element and i_column+1 element of  column_vec, which is a column vector of a matrix.
    inline void givensRotation(Eigen::Ref<Eigen::VectorXd> column_vec, const int i_column);

    // Generates a vector corresponding to b in Ax=b.
    virtual void bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec) = 0;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    virtual void axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec) = 0;

public:
    // Sets parameters in the matrix-free GMRES and allocate vectors and matrices used in the matrix-free GMRES.
    MatrixFreeGMRES(const int dim_equation, const int dim_krylov);

    // Resets parameters in the matrix-free GMRES and reallocate vectors and matrices used in the matrix-free GMRES.
    void resetParameters(const int dim_equation, const int dim_krylov);

    // Solves the GMRES and generates solution_update_vector, which is a solution of the matrix-free GMRES.
    void forwardDifferenceGMRES(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> solution_update_vec);
};


#endif