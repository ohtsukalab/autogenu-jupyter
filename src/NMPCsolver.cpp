#include <iostream>
#include <eigen3/Eigen/Core>


class newton_gmres : public matrixfree_gmres
{
private:
    int dimx, dimu, dv, dimeq;
    double hdir;
    Eigen::VectorXd u;
public:
    newton_gmres(const int dimx, const int dim_u, const int N, const double h) {
        dimx = dim_x;
        dimu = dim_u;
        dv = N;
        hdir = h;
        dimeq = dimu * (dv-1);
        u.resize(dimeq);
    }

};

class c_gmres : public matrixfree_gmres
{
private:
    int dimx, dimu, dv, dimeq;
    double hdir;
public:
    c_gmres(const int dim_x, const int dim_u, const int N, const double h) {
        dimx = dim_x;
        dimu = dim_u
        dv = N;
        dimeq = dv * dimu;
        Eigen::VectorXd U(dimeq);
        Eigen::MatrixXd x(dimx, dv), x1(dimx, dv), lmd(dimx, dv);
        hdir = h;
        
    };
    ~c_gmres();
    friend bfunc(const double t, const Eigen::VectorXd &x);
    void solver(const double t, const Eigen::VectorXd &x);
};

void bfunc(const double t, const Eigen::VectorXd x)
{
    double T, dtau;
    T = tf * (1 - std::exp(-alpha * t));
    dtau = T / dv;

    x.col(0) = x;
    x.col(i+1) = f(x.col(i), u.)
}

void c_gmres::solver(const double t, const Eigen::VectorXd &x)
{
    int i;
    
    for(i=0; i<dv; i++){
        euler(fx, x[i], u[i],  );
    }
    phix();

}
