#include "continuation_gmres.hpp"



continuation_gmres::continuation_gmres(const double T_f, const double alpha, const int division_num, const int k_max, const double h_dir, const double zeta) : matrixfree_gmres(k_max)
{
    dimx = model.dimx;
    dimu = model.dimu;
    dimuc = model.dimu + model.dimc;
    dv = division_num;
    hdir = h_dir;
    tf_len = T_f;
    a = alpha;
    z = zeta;
    kmax = k_max;
    dimeq = dimuc * dv;

    initgmres(dimeq);
    xtau.resize(dimx, dv+1);
    lmd.resize(dimx, dv+1);
    tmp.resize(dimx);
    u.resize(dimeq);
    u1.resize(dimeq);
    hu.resize(dimeq);
    hu1.resize(dimeq);
    hu2.resize(dimeq);
    du.resize(dimeq);
    xs.resize(dimx);
    uinit.resize(dimuc);

    for(int i=0; i<dimeq; i++){
        du(i) = 0.0;
    }
}


void continuation_gmres::init(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, const double r_tol, const int itr_max)
{
    Eigen::VectorXd uinit(dimuc);

    init_cgmres initializer(itr_max, r_tol, hdir, kmax);
    initializer.solvenocp(t, x0, u0, uinit);

    for(int i=0; i<dv; i++)
        u.segment(i*dimuc, dimuc) = uinit;    
}

void continuation_gmres::solvenmpc(const double t, const double ht, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> s)
{
    fdgmres(t, x, u, du);
    u = u + ht * du;

    s = u.segment(0, dimuc);
}


void continuation_gmres::Hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu)
{
    int i;
    double T, tau, htau;

    T = tf_len * (1.0 - std::exp(- a * t));
    htau = T / dv;

    xtau.col(0) = x;
    for(i=0, tau=t; i<dv; i++, tau+=htau){
        model.statefunc(tau, xtau.col(i), u.segment(i*dimuc, dimu), tmp);
        xtau.col(i+1) = xtau.col(i) + htau * tmp;
    }
    model.phix(tau, xtau.col(dv), lmd.col(dv));

    for(i=dv-1; i>=0; i--, tau-=htau){
        model.hxfunc(tau, xtau.col(i), u.segment(i*dimuc, dimuc), lmd.col(i+1), tmp);
        lmd.col(i) = lmd.col(i+1) + htau * tmp;
        model.hufunc(tau, xtau.col(i), u.segment(i*dimuc, dimuc), lmd.col(i+1), hu.segment(i*dimuc, dimuc));
    }
}


void continuation_gmres::Func(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> b)
{
    ts = t + hdir;
    model.statefunc(t, x, u.segment(0, dimu), tmp);
    xs = x + hdir * tmp;

    Hufunc(t, x, u, hu);
    Hufunc(ts, xs, u, hu1);
    b = - z * hu - (hu1 - hu) / hdir;

    u1 = u + hdir * du;
    Hufunc(ts, xs, u1, hu2);
    b -= (hu2 - hu1) / hdir;
}


void continuation_gmres::DhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu)
{
    u1 = u + hdir * v;
    Hufunc(ts, xs, u1, hu2);   
    dhu = (hu2 - hu1) / hdir;
}

