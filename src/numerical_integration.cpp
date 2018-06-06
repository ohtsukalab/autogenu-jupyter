#include <eigen3/Eigen/Core>
#include "numerical_integration.h"


/* Euler method for numerical integration
 * func: function of the dynamics, f(t, x1, x2)
 * t: time
 * x1:
 * x2:
 * tau: integration length
 * y: integrated solution
 */
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
    int i;
    double tstep;

    tstep = tau / 4;

    // Adams-Bashforth Predictor
    func(t, x1, x2, k1);
    y1 = y + tstep * k1;
    func(t+tstep, y1, x2, k2);
    y2 = y1 + tstep * (3/2 * k2 - 1/2 * k1);
    func(t+2*tstep, y2, x2, k3);
    y3 = y2 + tstep * (23/12 * k3 - 4/3 * k2 + 5/12 * k1);
    func(t+3*tstep, y3, x2, k4);
    y4 = y3 + tstep * (55/24 * k4 - 59/24 * k3 + 37/24 * k2 - 3/8 * k1);
    func(t+4*tstep, y4, x2, k5);
    y5 = y4 + tstep * (1901/720 * k5 - 1387/360 * k4 + 109/30 * k3 - 637/360 * k2 + 251/720 * k1);

    // Adams-Moulton Corrector
    func(t+4*tstep, y5, x2, k5);
    y4 = y5 + tstep * k5;
    func(t+3*tstep, y4, x2, k4);


    
}

