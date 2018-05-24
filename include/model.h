#include <iostream>
#include <eigen3/Eigen/Dense>


class model {
private:
    //parameters of the model
    double m1, m2, l1, l2, d1, d2, g, J1, J2;
    //parameters in the cost function
    double q[4], r[2], sf[4], xf[4]; 
public:
    model();
    ~model();
    void statefunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, Eigen::VectorXd &x1);
    void phix(double t, const Eigen::VectorXd &x, Eigen::VectorXd &lmd);
    void hxfunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &lmd, Eigen::VectorXd &lmd1);
    void hufunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &lmd, Eigen::VectorXd &hu);
};
