#include "matrixfree_gmres.hpp"


matrixfree_gmres::matrixfree_gmres(const int k_max)
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
    c.resize(kmax+1);
    s.resize(kmax+1);
    y.resize(kmax+1);
    g.resize(kmax+1);

    for(int i=0; i<kmax+1; i++)
        err(i) = 0;
}


void matrixfree_gmres::givappj(const int j, Eigen::Ref<Eigen::VectorXd> v)
{
    double tmp1, tmp2;

    tmp1 = c(j) * v(j) - s(j) * v(j+1);
    tmp2 = s(j) * v(j) + c(j) * v(j+1);

    v(j) = tmp1;
    v(j+1) = tmp2;
}

void matrixfree_gmres::givappj(const int j, Eigen::Ref<Eigen::VectorXd> c, Eigen::Ref<Eigen::VectorXd> s, Eigen::Ref<Eigen::VectorXd> v)
{
    double tmp1, tmp2;

    tmp1 = c(j) * v(j) - s(j) * v(j+1);
    tmp2 = s(j) * v(j) + c(j) * v(j+1);

    v(j) = tmp1;
    v(j+1) = tmp2;
}

double matrixfree_gmres::geterr()
{
    return err.norm();
}

void matrixfree_gmres::fdgmres(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> du)
{
    int i, j, k;
    double rho, nu;
    
    for(i=0; i<kmax+1; i++){
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
        if(h(k+1,k) != 0){
            v.col(k+1) = v.col(k+1) / h(k+1,k);
        }
        else {
            std::cout << "The modified Gram-Schmidt breakdown" << std::endl;
            break;
        }


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

//        std::cout << "\n g(k+1) = " << g(k+1) << std::endl;
        rho = std::fabs(g(k+1));
//        std::cout << "rho = " << rho << std::endl;
        err(k+1) = rho;
    }

    // solve H_k * y^k = g
    for(i=k-1; i>=0; i--){
        nu = g(i);
        for(j=i+1; j<k; j++){
            nu -= h(i,j) * c(j);
        }
        c(i) = nu / h(i,i);
    }
    for(i=0; i<n; i++){
        nu=0;
        for(j=0; j<k; j++){
            nu += v(i,j) * c(j);
        }
        du(i) += nu;
    }
}