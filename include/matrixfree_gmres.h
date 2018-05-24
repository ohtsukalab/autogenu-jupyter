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
    void fdgmres(Eigen::VectorXd& s, const Eigen::MatrixXd& x);
    virtual void Func(const Eigen::VectorXd& x, Eigen::VectorXd& s) = 0;
    virtual void DhFunc(const Eigen::VectorXd& x, const Eigen::VectorXd& v, Eigen::VectorXd& s) = 0;
};