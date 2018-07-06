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

    for(int i=0; i<dimeq; i++){
        du(i) = 0.0;
    }
}


void continuation_gmres::solvenmpc(const double t, const double ht, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> s)
{
    fdgmres(t, x, u, du);
    u = u + ht * du;

    s = u.segment(0, dimu);
}


void continuation_gmres::Hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> hu)
{
    int i;
    double T, tau, htau;

    T = tf_len * (1.0 - std::exp(- a * t));
    htau = T / dv;

    xtau.col(0) = x;
    for(i=0, tau=t; i<dv; i++, tau+=htau){
        model.statefunc(tau, xtau.col(i), u.segment(i*dimuc, dimuc), tmp);
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


void continuation_gmres::init_cgmres(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, const double r_tol, const int itr_max)
{
    int i;
    Eigen::VectorXd utmp(u0.size()), lmd0(x0.size()), dutmp(u0.size()), hu0(u0.size());

    i = 0;
    utmp = u0;
    model.phix(t, x0, lmd0);
    model.hufunc(t, x0, u0, lmd0, hu0);

    while((hu0.norm() > r_tol) && (i < 0)){
        fdgmres_init(t, x0, utmp, dutmp, u0.size(), kmax);
        utmp += dutmp;
        model.phix(t, x0, tmp);
        model.hufunc(t, x0, utmp, tmp, hu0);
        i++;
    }
//    std::cout << "end init, err = " << hu0.norm() << std::endl;
//    std::exit(1);
    for(i=0; i<dv; i++)
        u.segment(i*dimu, dimu) = utmp;
}


void continuation_gmres::initFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> s)
{
    Eigen::VectorXd lmd0(x.size());

    model.phix(t, x, lmd0);
    model.hufunc(t, x, u, lmd0, s);
}


void continuation_gmres::initDhFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::Ref<Eigen::VectorXd> dhu)
{
    Eigen::VectorXd u1(dimuc), hu(dimuc), hu1(dimuc);

    initFunc(t, x, u, hu);
    u1 = u + hdir * v;
    initFunc(t, x, u1, hu1);

    dhu = (hu1 - hu) / hdir;
}


void continuation_gmres::fdgmres_init(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> u1, const int n, const int kmax)
{
    int i, j, k;
    double rho, nu;
    Eigen::MatrixXd hmat(kmax+1, kmax+1), vmat(n, kmax+1);
    Eigen::VectorXd errvec(kmax+1), rvec(n), cvec(kmax), svec(kmax), yvec(kmax), gvec(kmax+1);

    
    for(i=0; i<kmax; i++){
        gvec(i) = 0.0;
        cvec(i) = 0.0;
        svec(i) = 0.0;
    }
    initFunc(t, x, u, rvec);
    
    rho = rvec.norm();
    gvec(0) = rho;
    errvec(0) = rho;
    vmat.col(0) = - rvec / rho;

    for(k=0; k<kmax; k++){
        // Modified Gram-Schmidt
        initDhFunc(t, x, u, vmat.col(k), vmat.col(k+1));
        for(j=0; j<=k; j++){
            hmat(j,k) = vmat.col(k+1).dot(vmat.col(j));
            vmat.col(k+1) -= hmat(j,k) * vmat.col(j);
        }


        hmat(k+1,k) = vmat.col(k+1).norm();
        if(hmat(k+1,k) != 0){
            vmat.col(k+1) = vmat.col(k+1) / hmat(k+1,k);
        }
        else {
            std::cout << "The modified Gram-Schmidt breakdown" << std::endl;
            break;
        }

        // Givens Rotation for the Lieast Squares Problem ||beta * e_1  - H_k * y^k||
        for(j=0; j<k; j++)
            givappj(j, cvec, svec, hmat.col(k));

        nu = std::sqrt(hmat(k,k)*hmat(k,k) + hmat(k+1,k)*hmat(k+1,k));
        if(nu != 0) {
            cvec(k) = hmat(k,k) / nu;
            svec(k) = - hmat(k+1,k) / nu;
            hmat(k,k) = cvec(k) * hmat(k,k) - svec(k) * hmat(k+1,k);
            hmat(k+1,k) = 0;
            givappj(k, cvec, svec, gvec);
        }
        else
            std::cout << "error : h(k,k) = h(k+1,k) = 0" << std::endl;

        rho = std::fabs(gvec(k+1));
        errvec(k+1) = rho;
    }

    // solve H_k * y^k = g
    for(i=kmax-1; i>=0; i--) {
        nu = gvec(i);
        for(j=i+1; j<kmax; j++){
            nu -= hmat(i,j) * yvec(j);
        }
        yvec(i) = nu / hmat(i,i);
    }

    for(i=0; i<n; i++){
        u1(i) = 0;
        for(j=0; j<kmax; j++){
            u1(i) += vmat(i,j) * yvec(j);
        }
    }
}