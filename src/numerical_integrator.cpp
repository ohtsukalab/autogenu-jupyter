#include "numerical_integrator.hpp"


numerical_integrator::numerical_integrator()
{
    tmp.resize(model.dimx);
    k1.resize(model.dimx);
    k2.resize(model.dimx);
    k3.resize(model.dimx);
    k4.resize(model.dimx);
}

void numerical_integrator::euler_state(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double htau, Eigen::Ref<Eigen::VectorXd> x1)
{
    model.statefunc(t, x, u, tmp);
    x1 = x + tmp * htau;
}

void numerical_integrator::euler_lmd(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, const double htau, Eigen::Ref<Eigen::VectorXd> lmd1)
{
    model.hxfunc(t, x, u, lmd, tmp);
    lmd1 = lmd + tmp * htau;
}

// The four-step Runge-Kutta-Gill method
void numerical_integrator::runge_kutta_gill(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double htau, Eigen::Ref<Eigen::VectorXd> x1)
{
    model.statefunc(t, x, u, k1);
    model.statefunc(t+0.5*htau, x + 0.5*htau*k1, u, k2);
    model.statefunc(t+0.5*htau, x + htau * 0.5 * (std::sqrt(2)-1) * k1 + htau * (1 - (1/std::sqrt(2))) * k2, u, k3);
    model.statefunc(t+htau, x - htau * 0.5 * std::sqrt(2) * k2 + htau * (1 + (1/std::sqrt(2))) * k3, u, k4);
    x1 = x + htau*( k1 + (2-std::sqrt(2))*k2 + (2+std::sqrt(2))*k3 + k4 )/6;
}
