#include "newton_gmres.hpp"


newton_gmres::newton_gmres(const double length_horizon, const int division_num, const int itr_max, const double h_diff, const int k_max) : matrixfree_gmres(division_num, k_max)
{
    tf = length_horizon;
    dv = division_num;
    htau = tf / dv;
    hdir = h_diff;
    dimx = model.dimx;
    dimu = model.dimu;
    dimuc = model.dimu + model.dimc;
    dimeq = dimuc * dv;
    imax = itr_max;

    initgmres(dimeq);

    xtau.resize(dimx, dv+1);
    lmd.resize(dimx, dv+1);
    tmp.resize(dimx);
    u.resize(dimeq);
    u1.resize(dimeq);
    hu.resize(dimeq);
    hu1.resize(dimeq);
    du.resize(dimeq);

    for(int i=0; i<dimeq; i++)
        u(i) = 0.0;

}

void newton_gmres::solvenmpc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> s)
{
    int i;

    for(i=0; i<imax; i++) {
        fdgmres(t, x, u, du);
        u += du;
    }
    s = u.segment(0, dimu);
}


void newton_gmres::Func(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu)
{
    int i;
    xtau.col(0) = x;
    for(i=0, tau=t; i<dv; i++, tau+=htau) {
        model.statefunc(tau, xtau.col(i), u.segment(i*dimuc, dimuc), tmp);
        xtau.col(i+1) = xtau.col(i) + htau * tmp;
    }
    model.phix(t, xtau.col(dv), lmd.col(dv));

    for(i=dv-1, tau=t+tf; i>=0; i--, tau-=htau) {
        model.hxfunc(tau, xtau.col(i), u.segment(i*dimuc, dimuc), lmd.col(i+1), tmp);
        lmd.col(i) = lmd.col(i+1) + htau * tmp;
        model.hufunc(t, xtau.col(i), u.segment(i*dimuc, dimuc), lmd.col(i+1), hu.segment(i*dimuc, dimuc));
    }
}

void newton_gmres::DhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu)
{
    Func(t, x, u, hu);

    u1 = u + hdir * v;
    Func(t, x, u1, hu1);

    dhu = (hu1 - hu) / hdir;
}