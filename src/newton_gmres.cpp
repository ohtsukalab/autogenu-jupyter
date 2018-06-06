#include <eigen3/Eigen/Core>
#include "newton_gmres.h"


newton_gmres::newton_gmres(const nmpc_model m, const int division_num, const double h_dir, const int newtonitr_max, const double rho, const int k_max, const double T)
{
    dimx = m.dimx;
    dimu = m.dimu;
    dv = division_num;
    hdir = h_dir;
    dimeq = dimu * (dv-1);
    dtau = T / division_num;
    model = m;

    u.resize(dimeq);
    u1.resize(dimeq);    
    hu.resize(dimeq);
    hu1.resize(dimeq);
    du.resize(dimeq);
    x.resize(dimx,dv);
    lmd.resize(dimx, dv);
    u = 0;
    u1 = 0;
    hu = 0;
    hu1 = 0;
    du = 0;
    x = 0;
    lmd = 0;
}


void newton_gmers::solver(const double t, const Eigen::VectorXd& x, Eigen::VectorXd& s)
{
    for(i=0; i<i_max; i++){
        fdgmres(t, x0, du, u);
        u += du;
        if(err.norm() < rho)
            break;
    }
    s = u.segment(i*dimu, dimu);
}


void newton_gmres::Hufunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::MatrixXd& hu)
{
    int i;
    double tau;

    x.col(0) = x0;

    for(i=0, tau=t; i<dv; i++){
        model.statefunc(tau, x.col(i), u.segment(i*dimu, dimu), tmp);
        x.col(i+1) = x.col(i) + dtau * x.col(i+1);
        tau += dtau;
    }
    model.phix(tau, x, lmd.col(dv));
    for(; i>0; i--){
        model.hxfunc(tau, x.col(i), u.segment(i*dimu, dimu), lmd.col(i+1), tmp);
        lmd.col(i) = lmd.col(i+1) + dtau * tmp;
        model.hufunc(taut, x.col(i), u.segment(i*dimu, dimu), lmd.col(i+1), hu.segment(i*dimu, dimu));
    }
}


void newton_gmers::Func(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& s)
{
    Hufunc(t, x0, u, s);
}

void newton_gmers::DhFunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& s)
{
    u1 = u + hdir * v;
    Hufunc(t, x0, u1, hu1);
    s = (hu1 - hu) / h;
}
