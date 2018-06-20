#ifndef NUMERICAL_INTEGRATOR_H
#define NUMERICAL_INTEGRATOR_H

#include "nmpc_model.hpp"
#include <eigen3/Eigen/Core>

class numerical_integrator{
private:
    nmpc_model model;
    Eigen::VectorXd tmp, k1, k2, k3, k4;
public:
    numerical_integrator();
    void euler_state(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double htau, Eigen::Ref<Eigen::VectorXd> x1);
    void euler_lmd(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, const double htau, Eigen::Ref<Eigen::VectorXd> lmd1);
    void runge_kutta_gill(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double htau, Eigen::Ref<Eigen::VectorXd> x1);
};

#endif