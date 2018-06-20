/*
 *   This supportes the newton_gmres method, a numerical solver for nonlienar optimization problem.
 *   This program is witten with reference to "C. T. Kelly, Iterative Methods for Linear and Nonlinear Equations, Frontiers in Apllied Mathematics, SIAM (1995)"
 */

#ifndef NEWTON_GMRES_H
#define NEWTON_GMRES_H

#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "numerical_integrator.hpp"
#include "nmpc_model.hpp"

class newton_gmres : public matrixfree_gmres{
private:
/* --------------------------------------------------
 * n        the dimension of the nonlinear equation
 * x1       a temporary vector
 * hdir     increment of differential derivative
 * Func()   nonlinear function
 * DhFunc() product of partially derivative of the nonlinear function and a vector
 * -------------------------------------------------- */
    Eigen::MatrixXd xtau, lmd;
    Eigen::VectorXd u, u1, du, hu, hu1, tmp;

    nmpc_model model;
    int dimeq, dv, dimx, dimu, dimuc, imax;
    double hdir, tau, htau, tf;

public:
    newton_gmres(const double length_horizon, const int division_num, const int itr_max, const double h_diff, const int k_max);
    void solvenmpc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> s);
    void Func(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu);
    void DhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu);
};

#endif