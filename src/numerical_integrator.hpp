#ifndef NUMERICAL_INTEGRATOR_H
#define NUMERICAL_INTEGRATOR_H


#include "nmpc_model.hpp"
#include <eigen3/Eigen/Core>

class NumericalIntegrator{
private:
    NMPCModel model_;
    Eigen::VectorXd dx_vec_, k1_vec_, k2_vec_, k3_vec_, k4_vec_;
public:
    NumericalIntegrator(const NMPCModel model);

    // Euler method for the state vector
    void eulerState(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec);

    // Euler method for the lambda vector, the Lagrange multiplier of the state equation
    void eulerLambda(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const Eigen::VectorXd& next_lambda_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> current_lambda_vec);

    // The four-step Runge-Kutta-Gill method
    void rungeKuttaGill(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec);
};

#endif