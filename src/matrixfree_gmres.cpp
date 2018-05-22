#include <Eigen3/Core>
#include <iostream>


class matrixfree_gmres {
private:
    void (*Func)(), (*DhFunc)();
    int dim, kmax;
    double r_tol;
    Eigen::MatrixXd h, v;
    Eigen::VectorXd errvec;
public:
    matrixfree_gmres(const void (*F)(), const void (*DhF)(), const int dimf, const int k_max, const double rtol){
        (*Func)() = F;
        (*DhFunc)() = DhF;
        dim = dimf;
        kmax = k_max;
        r_tol = rtol;
        h.resize(kmax+1,kmax+1);
        v.resize(kmax+1,dim);
        errvec.resize(kmax+1);
    };
    ~matrixfree_gmres();
    void fdgmres(Eigen::VectorXd &b, Eigen::MatrixXd &x, Eigen::VectorXd &err);
};


// s: solution
// x: current s x
// dhfunc: function for newton
// eta: convergence radius
// rho: error vector
void matrixfree_gmres::fdgmres(Eigen::VectorXd &s, const Eigen::MatrixXd &x)
{
    int i, j, k;
    double beta, rho, nu, w1, w2;
    Eigen::VectorXd r(dim), cvec(kmax+1), svec(kmax+1), gvec(kmax+1);

    s = 0;
    Func(x, r);
    v.col(0) = r / r.norm();
    rho = r.norm();
    beta = rho;
    for(k=0; k<kmax; k++){
        // Modified Gram-Schmidt
        DhFunc(x, v.col(k+1));
        for(j=0; j<k; j++){
            h(j,k) = v.col(k+1) * v.col(j);
            v.col(k+1) = v.col(k+1) - h(j,k) * v.col(j);
        }
        h(k+1,k) = v.col(k+1).norm();
        if(h(k+1,k) != 0)
            v.col(k+1) = v.col(k+1) / h(k+1,k);
        else
            std::cout << "fgmres() : breakdown" << std::endl;

        // Givens Rotation for the Lieast Squares Problem ||beta * e_1  - H_k * y^k||
        for(j=0; j<k; j++){
            w1 = cvec[i] * gvec[i] - svec[i] * gvec[i+1];
            w2 = svec[i] * gvec[i] + cvec[i] * gvec[i+1];
            gvec[i] = w1;
            gvec[i+1] = w2;
        }
        nu = std::sqrt(h(k,k)*h(k,k) + h(k+1,k)*h(k+1,k));
        if(nu != 0) {
            cvec[k] = h(k,k) / nu;
            svec[k] = - h(k+1,k) / nu;
            h(k,k) = cvec[k] * h(k,k) - svec[k] * h(k+1,k);
            h(k+1,k) = 0;
            w1 = cvec[k] * gvec[k] - svec[i] * gvec[k+1];
            w2 = svec[k] * gvec[k] + cvec[i] * gvec[k+1];
        }
        else
            std::cout << "error : h(k,k) = h(k+1,k) = 0\n";
        rho = std::fabs(g[k+1]);
        errvec[k+1] = rho;
    }

    // solve H_k * y^k = gvec
    for(i=kmax; i>=0; i--) {
        for(nu = gvec[i], j=i+1; j<k; j++){
            nu -= h(i,j) * cvec[j];
        }
        cvec[i] = nu / h(i,i);
    }
    s = v * cvec;
}

