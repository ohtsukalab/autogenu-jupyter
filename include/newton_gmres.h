#include <iostream>
#include <eigen3/Eigen/Core>


class newton_gmres : public matrixfree_gmres, public model {
private:
    int dimx, dimu, dv, dimeq, i_max;
    double hdir, dtau, rho;
    Eigen::MatrixXd x, lmd;
    Eigen::VectorXd u, u1, hu, hu1, du;
public:
    newton_gmres(const int dim_x, const int dim_u, const int division_num, const double h_dir, const int newtonitr_max, const double rho, const int k_max) : matrixfree_gmres(dim_u * (division_num-1), k_max);
    ~newton_gmres();
    void nsolgm(const Eigen::VectorXd& x, Eigen::VectorXd& s);
    void Func(const Eigen::VectorXd& u, Eigen::VectorXd& s);
    void DhFunc(const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& s);
    void Hufunc(const double t, const Eigen::MatrixXd& u, Eigen::MatrixXd& hu);
};