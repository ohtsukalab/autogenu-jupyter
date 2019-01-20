//
// Supportes numerical integration methods, the Euler method and the four-step Runge-Kutta-Gill method.
//

#ifndef NUMERICAL_INTEGRATOR_H
#define NUMERICAL_INTEGRATOR_H


#include "nmpc_model.hpp"
#include <eigen3/Eigen/Core>


// Supports numerical integration of the state equation of the system described in nmpc_model.hpp for numerical simnulations.
class NumericalIntegrator{
private:
    NMPCModel model_;
    Eigen::VectorXd dx_vec_, k1_vec_, k2_vec_, k3_vec_, k4_vec_;
public:
    // Allocates vectors.
    NumericalIntegrator();

    // Euler method for the state equation.
    Eigen::VectorXd euler(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length);

    // The four-step Runge-Kutta-Gill method for the state equation.
    Eigen::VectorXd rungeKuttaGill(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length);
};


#endif