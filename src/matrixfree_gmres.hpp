/*
 *   This supportes the newton_gmres method, a numerical solver for nonlienar optimization problem.
 *   This program is witten with reference to "C. T. Kelly, Iterative Methods for Linear and Nonlinear Equations, Frontiers in Apllied Mathematics, SIAM (1995)"
 */

#ifndef MATRIXFREE_GMRES_H
#define MATRIXFREE_GMRES_H

#include "nmpc_model.hpp"
#include <iostream>
#include <eigen3/Eigen/Core>


class matrixfree_gmres {
private:
/* --------------------------------------------------
 * n        the dimension of the nonlinear equation
 * kmax     the maximum iteration numebr of the gmres
 * rho      the convergence radius
 * h        Hessian matrix
 * v        Hessian matrix
 * err   vector composed of the erros in gmres
 * Func()   nonlinear function
 * DhFunc() product of partially derivative of the nonlinear function and a vector
 * -------------------------------------------------- */

    int n, kmax;
    nmpc_model model;
    Eigen::MatrixXd h, v;
    Eigen::VectorXd err, r, c, s, g, y;
public:
    matrixfree_gmres(const int division_num, const int k_max);
    void initgmres(const int dim);
    void givappj(const int j, Eigen::Ref<Eigen::VectorXd> v);
    void fdgmres(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> u1);
    virtual void Func(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu) = 0;
    virtual void DhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& du, Eigen::Ref<Eigen::VectorXd> dhu) = 0;
};

#endif