#include <eigen3/Eigen/Core>
#include "continuation_gmres.h"

class continuation_gmres : public matrixfree_gmres {
private:
    int dimx, dimu, dv, dimeq;
    double hdir, Tf, a, z, ts;
    Eigen::MatrixXd x, lmd;
    Eigen::VectorXd u, u1, hu, hu1, du, xs, x1s;
public:
    continuation_gmres(const nmpc_model m, const int division_num, const int k_max, const double h_dir, const double T_f, const double alpha, const double zeta) : matrixfree_gmres(m.dimu * (division_num-1), k_max);
    ~continuation_gmres();
    void init(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, const double r_tol, const double newtonitr_max);
    void solver(const double t, const Eigen::VectorXd& x, Eigen::VectorXd& u);
    void Hufunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& hu);
    void Func(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& s);
    void DhFunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& s);
};


continuation_gmres::continuation_gmres(const nmpc_model m, const int division_num, const int k_max, const double h_dir, const double T_f, const double alpha, const double zeta)
{
    dimx = m.dimx;
    dimu = m.dimu;
    dv = division_num;
    hdir = h_dir;
    Tf = T_f;
    a = alpha;
    z = zeta;
    dimeq = dimu * (dv-1);

    u.resize(dimeq);
    u1.resize(dimeq);    
    hu.resize(dimeq);
    hu1.resize(dimeq);
    du.resize(dimeq);
    x.resize(dimx,dv);
    lmd.resize(dimx, dv);
    xs.resize(dimx);
    x1s.resize(dimx);

    u = 0;
    u1 = 0;
    hu = 0;
    hu1 = 0;
    du = 0;
    x = 0;
    lmd = 0;    
}


void continuation_gmres::init(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, const double r_tol, const double newtonitr_max)
{
    Eigen::VectorXd lmd0(dimx), utmp(dimu);
    int i;
    
    utmp = u0;
    for(i=0; i<newtonitr_max; i++){
        model.phix(t, x0, lmd0);
        model.hufunc(t, x0, utmp, lmd0, hu0);
        fdgmres(t, x0, x, )
    }

    for(i=0; i<dv; i++){
        u.segment(i*dimu, dimu) = utmp;
    }
}

void continuation_gmres::initfunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, const double r_tol, const double newtonitr_max)
{

}


void continuation_gmres::initdhfunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, const double r_tol, const double newtonitr_max)
{

}


void continuation_gmres::solver(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& s)
{
    double T, dtau;

    fdgmres(t, x0, u, du);
    u += du * ht;
    s = u.segment(i*dimu, dimu);
}


void newton_gmres::Hufunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::MatrixXd& hu)
{
    int i;
    double T, tau, dtau;

    T = T_f * (1 - std::exp(- alpha * t));
    dtau = T / division_num;

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

void continuation_gmres::Func(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& s)
{
    ts = t + hdir;
    statefcun(t, x0, u, x1s);
    xs = x0 + x1s * hdir;

    Hufunc(t, x0, u, hu);
    Hufunc(ts, xs, u, hu1);
    s = ((1 - hdir * z) * hu - hu1) / hdir;
}

void continuation_gmres::DhFunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& s)
{
    u1 = u + hdir * v;
    Hufunc(ts, xs, u, hu);
    Hufunc(ts, xs, u1, hu1);
    s = (hu1 - hu) / hdir;
}
