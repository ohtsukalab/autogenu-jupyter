#include "matrixfree_gmres.hpp"


matrixfree_gmres::matrixfree_gmres(const int division_num, const int k_max)
{
    kmax = k_max;
}

void matrixfree_gmres::initgmres(const int dim)
{
    n = dim;
    h.resize(kmax+1, kmax+1);
    v.resize(n, kmax+1);
    err.resize(kmax+1);
    r.resize(n);
    c.resize(kmax);
    s.resize(kmax);
    y.resize(kmax);
    g.resize(kmax+1);
}


void matrixfree_gmres::givappj(const int j, Eigen::Ref<Eigen::VectorXd> v)
{
    double tmp1, tmp2;

    tmp1 = c(j) * v(j) - s(j) * v(j+1);
    tmp2 = s(j) * v(j) + c(j) * v(j+1);

    v(j) = tmp1;
    v(j+1) = tmp2;
}


void matrixfree_gmres::fdgmres(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> u1)
{
    int i, j, k;
    double rho, nu;
    
    for(i=0; i<kmax; i++){
        g(i) = 0.0;
        c(i) = 0.0;
        s(i) = 0.0;
    }
    Func(t, x, u, r);
    rho = r.norm();
    g(0) = rho;
    err(0) = rho;
    v.col(0) = r / rho;

    for(k=0; k<kmax; k++){
        // Modified Gram-Schmidt
        DhFunc(t, x, u, v.col(k), v.col(k+1));
        for(j=0; j<=k; j++){
            h(j,k) = v.col(k+1).dot(v.col(j));
            v.col(k+1) -= h(j,k) * v.col(j);
        }
        h(k+1,k) = v.col(k+1).norm();
        if(h(k+1,k) != 0)
            v.col(k+1) = v.col(k+1) / h(k+1,k);
        else
            std::cout << "arnordi process breaks down" << std::endl;

        // Givens Rotation for the Lieast Squares Problem ||beta * e_1  - H_k * y^k||
        for(j=0; j<k; j++)
            givappj(j, h.col(k));

        nu = std::sqrt(h(k,k)*h(k,k) + h(k+1,k)*h(k+1,k));
        if(nu != 0) {
            c(k) = h(k,k) / nu;
            s(k) = - h(k+1,k) / nu;
            h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k);
            h(k+1,k) = 0;
            givappj(k,g);
        }
        else
            std::cout << "error : h(k,k) = h(k+1,k) = 0" << std::endl;

        rho = std::fabs(g(k+1));
        err(k+1) = rho;
    }

    // solve H_k * y^k = g
    for(i=kmax-1; i>=0; i--) {
        nu = g(i);
        for(j=i+1; j<kmax; j++){
            nu -= h(i,j) * y(j);
        }
        y(i) = nu / h(i,i);
    }

    for(i=0; i<n; i++){
        u1(i) = 0;
        for(j=0; j<kmax; j++){
            u1(i) += v(i,j) * y(j);
        }
    }
}
