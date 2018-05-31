#include <eigen3/Eigen/Core>




class numerical_integration{
private:
    int dimx;
public:
    numerial_integration((dim_x){
        dimx = dimx;
        Eigen::VectorXd a(dim_x), k1(dim_x), k2(dim_x), k3(dim_x), k2(dim_x);
    };
    void euler(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y);
    void euler(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const Eigen::VectorXd& x3, const double tau, Eigen::VectorXd &y);
    void runge_kutta(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y);
    void adams(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y);
};

void numerical_integration::euler(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y);
{
    func(t, x1, x2, a);
    y = x1 + a * tau;
}

void numerical_integration::euler(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const Eigen::VectorXd& x3, const double tau, Eigen::VectorXd &y);
{
    func(t, x1, x2, x3, a);
    y = x3 + a * tau;
}


void numerical_integration::runge_kutta(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y)
{
    Eigen::VectorXd (dim_x);
    func(t, x1, x2, tmp);
    k1 = tau * tmp;
    func(t+0.5*tau, x1+0.5*k0, x2, tmp);
    k2 = tau * tmp;
    func(t+0.5*tau, x1+0.5*k1, x2, tmp);
    k3 = tau * tmp;
    func(t+0.5*tau, x1+k3, x2, tmp);
    k4 = tau * tmp;

    y = (k1 + 2*k2 + 2*k3 + k4)/6;
}

void numerical_integration::adams(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y)
{
    
}

