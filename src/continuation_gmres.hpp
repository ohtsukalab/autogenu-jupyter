/*
 *   This suppoerts a solver for continuation/GMRES(C/GMRES) method for nonlinear model predictive control
 *   This program is witten with reference to "T. Ohtsuka, "
 */

#ifndef CONTINUATION_GMRES_H
#define CONTINUATION_GMRES_H


#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "numerical_integrator.hpp"
#include "nmpc_model.hpp"

class continuation_gmres : virtual public matrixfree_gmres, virtual public numerical_integrator{
private:
    nmpc_model model;
    int dimx, dimu, dimuc, dv, dimeq, kmax;
    double hdir, tf_len, a, z, ts;
    Eigen::MatrixXd xtau, lmd;
    Eigen::VectorXd u, u1, hu, hu1, hu2, du, xs, tmp;
public:
    continuation_gmres(const double T_f, const double alpha, const int division_num, const int k_max, const double h_dir, const double zeta);
    void init_cgmres(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, const double r_tol, const int itr_max);
    void fdgmres_init(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> u1, const int n, const int kmax);
    void initFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> s);
    void initDhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu);
    void solvenmpc(const double t, const double ht, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> s);
    void Hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu);
    void Func(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> b);
    void DhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu);
};

#endif