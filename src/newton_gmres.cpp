#include <eigen3/Eigen/Core>
#include "matrixfree_gmres"
#include "newton_gmres.h"


newton_gmres::newton_gmres(const int dim_x, const int dim_u, const int division_num, const double h_dir, const int newtonitr_max, const double rho, const int k_max, const double T) : matrixfree_gmres(dim_u * (division_num-1), k_max)
{
    dimx = dim_x;
    dimu = dim_u;
    dv = division_num;
    hdir = h_dir;
    dimeq = dimu * (dv-1);
    dtau = T / division_num;

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


void newton_gmers::nsolgm(const double t, const Eigen::VectorXd& x, Eigen::VectorXd& s)
{
    double t;
    for(i=0; i<i_max; i++){
        fdgmres(t, x0, du, u);
        u += du;
        if(errvec.norm() < rho)
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
        statefunc(tau, x.col(i), u.segment(i*dimu, dimu), tmp);
        x.col(i+1) = x.col(i) + dtau * x.col(i+1);
        tau += dtau;
    }
    phix(tau, x, lmd.col(dv));
    for(; i>0; i--){
        hxfunc(tau, x.col(i), u.segment(i*dimu, dimu), lmd.col(i+1), tmp);
        lmd.col(i) = lmd.col(i+1) + dtau * tmp;
        hufunc(taut, x.col(i), u.segment(i*dimu, dimu), lmd.col(i+1), hu.segment(i*dimu, dimu));
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
