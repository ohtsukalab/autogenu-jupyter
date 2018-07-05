#include "newton_gmres.hpp"


newton_gmres::newton_gmres(const double length_horizon, const int division_num, const double conv_radius, const double h_diff, const int k_max) : matrixfree_gmres(k_max)
{
    tf = length_horizon;
    dv = division_num;
    htau = tf / dv;
    hdir = h_diff;
    dimx = model.dimx;
    dimu = model.dimu;
    dimuc = model.dimu + model.dimc;
    dimeq = dimuc * dv;
    rho = conv_radius;

    initgmres(dimeq);
    xtau.resize(dimx, dv+1);
    lmd.resize(dimx, dv+1);
    tmp.resize(dimx);
    u.resize(dimeq);
    u1.resize(dimeq);
    u2.resize(dimeq);
    hu.resize(dimeq);
    hu1.resize(dimeq);
    du.resize(dimeq);

    for(int i=0; i<dimeq; i++)
        u(i) = 0.0;
}


void newton_gmres::solvenmpc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> s)
{
    int i;
    double r;

    i=0;
    Func(t, x, u, hu);
    r = hu.norm();

    while(r > rho){
        fdgmres(t, x, u, du);
        linesearch(t, x, u, du);
        Func(t, x, u, hu);
        r = hu.norm();
    }
    s = u.segment(0, dimu);
}


void newton_gmres::linesearch(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> u, const Eigen::VectorXd& du)
{
    int k;
    double alpha, aupper, alower, alpha1, alpha2, cost, cost0, cost1, dcost0;

// bracketing for linesearch
    k=0;
    cost0 = costsearch(t, x, u, du, 0);
    dcost0 = dcostsearch(t, x, u, du, 0);

    cost = cost0;
    cost1 = costsearch(t, x, u, du, h); 
    alpha = h;
    while( (cost1 < cost) && ( (cost1 > cost0 + mu1 * dcost0 * alpha) && (dcostsearch(t, x, u, du, alpha) < mu2 * dcost0)) ){
        cost = cost1;
        alpha += h;
        cost1 = costsearch(t, x, u, du, alpha);
        k++;
    }

    if(k){
        alower = (k-1) * h;
        aupper = (k+1) * h;
    }
    else {
        alower = 0;
        aupper = h;
    }

// golden section method
    alpha1 = alower+(1-r)*(aupper-alower);
    alpha2 = alower+r*(aupper-alower);
    while ( (std::abs(aupper - alower) >= d) && ( (cost1 > cost0 + mu1 * dcost0 * alpha) && (dcostsearch(t, x, u, du, alpha) < mu2 * dcost0)) ) {
        if(costsearch(t, x, u, du, alpha1) < costsearch(t, x, u, du, alpha2)) {
            aupper = alpha2;
            alpha2 = alpha1;
            alpha1 = alower + (1-r)*(aupper-alower);
        }
        else {
            alower = alpha1;
            alpha1 = alpha2;
            alpha2 = alower + r*(aupper-alower);
        }
    }
    u = u + 0.5 * (aupper + alower) * du;
}


double newton_gmres::costsearch(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& du, const double alpha)
{
    int i;
    double cost;

    u1 = u + alpha * du;

    cost=0;
    xtau.col(0) = x;
    for(i=0, tau=t; i<dv; i++, tau+=htau) {
        model.statefunc(tau, xtau.col(i), u1.segment(i*dimuc, dimuc), tmp);
        xtau.col(i+1) = xtau.col(i) + htau * tmp;
        cost += htau * model.stagecost(tau, xtau.col(i), u1.segment(i*dimuc, dimuc));
    }
    cost += model.terminalcost(tau, xtau.col(dv));

    return cost;
}


double newton_gmres::dcostsearch(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& du, const double alpha)
{
    int i;

    u1 = u + alpha * du;

    for(i=0, tau=t; i<dv; i++, tau+=htau) {
        model.Lufunc(tau, x, u1.segment(i*dimuc, dimuc), hu);
    }

    hu = htau * hu;
    return hu.dot(du);
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