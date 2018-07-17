#include "init_cgmres.hpp"


init_cgmres::init_cgmres(const int itr_max, const double r_tol, const double h_dir, const int k_max) : matrixfree_gmres(k_max)
{
    lmd.resize(model.dimx);
    u1.resize(model.dimu+model.dimc);
    hu1.resize(model.dimu+model.dimc);
    hu2.resize(model.dimu+model.dimc);
    du.resize(model.dimu+model.dimc);

    hdir = h_dir;
    r = r_tol;
    imax = itr_max;
    initgmres(model.dimu+model.dimc);

    for(int i=0; i<(model.dimu+model.dimc); i++)
        du(i) = 0.0;
}


void init_cgmres::solvenocp(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, Eigen::Ref<Eigen::VectorXd> u)
{
    int i;
    
    u = u0;
    Hufunc(t, x0, u, hu1);
    while(hu1.norm() < r && i > imax){
        fdgmres(t, x0, u, du);
        u += du;
        Hufunc(t, x0, u, hu1);
        i++;
    }
}


void init_cgmres::Hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu)
{
    model.phixfunc(t, x, lmd);
    model.hufunc(t, x, u, lmd, hu);
}


void init_cgmres::Func(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> b)
{
    Hufunc(t, x, u, hu1);
    b = - hu1;
}


void init_cgmres::DhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu)
{
    Hufunc(t, x, u, hu1);
    u1 = u + hdir * v;
    Hufunc(t, x, u1, hu2);

    dhu = (hu2 - hu1) / hdir;
}