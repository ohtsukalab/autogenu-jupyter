#ifndef INIT_CGMRES_H
#define INIT_CGMRES_H


#include "matrixfree_gmres.hpp"
#include <eigen3/Eigen/Core>


class init_cgmres : public matrixfree_gmres{
private:
    nmpc_model model;
    int imax;
    double r, hdir;
    Eigen::VectorXd lmd, u1, hu1, hu2, du;

public:
    init_cgmres(const int itr_max, const double r_tol, const double h_dir, const int k_max);
    void solvenocp(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, Eigen::Ref<Eigen::VectorXd> u);
    void Hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu);
    void Func(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> b);
    void DhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu);
};

#endif