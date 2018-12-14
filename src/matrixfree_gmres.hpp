/*
 *   This supportes the forward-difference GMRES method, a numerical solver for nonlienar programming problem.
 *   This program is witten with reference to "C. T. Kelly, Iterative Methods for Linear and Nonlinear Equations, Frontiers in Apllied Mathematics, SIAM (1995)"
 */

#ifndef MATRIXFREE_GMRES_H
#define MATRIXFREE_GMRES_H

#include <iostream>
#include <eigen3/Eigen/Core>
#include "nmpc_model.hpp"


class MatrixFreeGMRES {
private:
    int dim_equation_, max_dim_krylov_;
    Eigen::MatrixXd hessenberg_mat_, basis_mat_;
    Eigen::VectorXd current_error_vec_, givens_c_vec_, givens_s_vec_, g_vec_;

    inline void givensRotation(Eigen::Ref<Eigen::VectorXd> column_vec, const int i_column);
    virtual void nonlinearEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec) = 0;
    virtual void forwardDifferenceEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec) = 0;
public:
    MatrixFreeGMRES(const int dim_equation, const int dim_krylov);
    void resetParameters(const int dim_equation, const int dim_krylov);
    void forwardDifferenceGMRES(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> solution_update_vec);
};

#endif