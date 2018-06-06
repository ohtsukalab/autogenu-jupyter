#include <eigen3/Eigen/Core>

class numerical_integration{
private:
    int dimx;
public:
    numerial_integration((dim_x){
        dimx = dimx;
        Eigen::VectorXd a(dim_x), k1(dim_x), k2(dim_x), k3(dim_x), k4(dim_x), k5(dim_x), y1(dim_x), y2(dim_x), y3(dim_x), y4(dim_x), y5(dim_x);
    };
    void euler(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y);
    void euler(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const Eigen::VectorXd& x3, const double tau, Eigen::VectorXd &y);
    void runge_kutta(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y);
    void adams(void (*func()), const double t, const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double tau, Eigen::VectorXd &y);
};