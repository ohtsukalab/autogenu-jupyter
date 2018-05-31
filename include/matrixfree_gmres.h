/*
 *   This suppoerts a solver for matrixfree gmres method
 *   This program is witten with reference to "C. T. Kelly, Iterative methods for linear and nonlinear"
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
    Eigen::VectorXd errvec;
public:
    matrixfree_gmres(const int dimx, const int k_max);
    ~matrixfree_gmres();
    void fdgmres(const double t, const Eigen::VectorXd& x0, const Eigen::MatrixXd& , Eigen::VectorXd& s);
    virtual void Func(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& s) = 0;
    virtual void DhFunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& s) = 0;
};

#endif