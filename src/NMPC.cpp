#include <iostream>
#include <eigen3/Eigen/Dense>


class model {
public:
    model(dimx, dimy)
    ~model();
    void statefunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, Eigen::VectorXd &x1);
    void hxfunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &lmd, Eigen::VectorXd &lmd1);
    void hufunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &lmd, Eigen::VectorXd &hu);
};

void model::statefunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, Eigen::VectorXd &x1)
{
    x1[0] = ();
    x1[1] = ();
}

void model::hxfunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &lmd, Eigen::VectorXd &lmd1)
{
    lmd1[0] = ();
}

void hufunc(double t, const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &lmd, Eigen::VectorXd &hu)
{
    hu[0] = ();
}




int main()
{

    return 0;
}