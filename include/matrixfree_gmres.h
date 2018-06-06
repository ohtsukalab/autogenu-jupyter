/*
 *   This suppoerts a solver for matrixfree gmres method
 *   This program is witten with reference to "C. T. Kelly, Iterative methods for linear and nonlinear equations, Frontiers in Apllied Mathematics, SIAM (1995)".
 */

#ifndef MATRIXFREE_GMRES_H
#define MATRIXFREE_GMRES_H

#include <iostream>
#include <eigen3/Eigen/Core>


class matrixfree_gmres {
private:

/* --------------------------------------------------
 * n        the dimension of the nonlinear equation
 * kmax     the maximum iteration numebr of the gmres
 * h        Hessian matrix
 * v        Hessian matrix
 * errvec   vector composed of the erros in gmres
 * Func()   nonlinear function
 * DhFunc() product of partially derivative of the nonlinear function and a vector
 * -------------------------------------------------- */

    int n, kmax;
    Eigen::MatrixXd h, v;
    Eigen::VectorXd err;
public:
    matrixfree_gmres(const int dim, const int k_max);
    ~matrixfree_gmres();
    void fdgmres(const double t, const Eigen::VectorXd& x0, const Eigen::MatrixXd& , Eigen::VectorXd& s);
    void fdgmres(const double t, const Eigen::VectorXd& x0, const Eigen::MatrixXd& x, Eigen::VectorXd& s, const int dim, const int kmax, void (*func()), void void (*dhfunc()));
    virtual void Func(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& s) = 0;
    virtual void DhFunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& s) = 0;
};

#endif