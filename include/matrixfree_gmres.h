#include <iostream>
#include <eigen3/Eigen/Core>


class matrixfree_gmres {
private:
    int n, kmax;
    Eigen::MatrixXd h, v;
    Eigen::VectorXd errvec;
public:
    matrixfree_gmres(const int dimx, const int k_max);
    ~matrixfree_gmres();
    void fdgmres(const double t, const Eigen::VectorXd& x0, const Eigen::MatrixXd& u, Eigen::VectorXd& s);
    virtual void Func(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, Eigen::VectorXd& s) = 0;
    virtual void DhFunc(const double t, const Eigen::VectorXd& x0, const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& s) = 0;
};