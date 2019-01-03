#include "numerical_integrator.hpp"


NumericalIntegrator::NumericalIntegrator(const NMPCModel model)
{
    model_ = model;

    // allocate vectors
    dx_vec_.resize(model_.dimState());
    k1_vec_.resize(model_.dimState());
    k2_vec_.resize(model_.dimState());
    k3_vec_.resize(model_.dimState());
    k4_vec_.resize(model_.dimState());
}


void NumericalIntegrator::euler(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec)
{
    model_.stateFunc(current_time, current_state_vec, control_input_vec, dx_vec_);
    next_state_vec = current_state_vec + integration_length * dx_vec_;
}


void NumericalIntegrator::rungeKuttaGill(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec)
{
    model_.stateFunc(current_time, current_state_vec, control_input_vec, k1_vec_);
    model_.stateFunc(current_time+0.5*integration_length, current_state_vec+0.5*integration_length*k1_vec_, control_input_vec, k2_vec_);
    model_.stateFunc(current_time+0.5*integration_length, current_state_vec+integration_length*0.5*(std::sqrt(2)-1)*k1_vec_+integration_length*(1-(1/std::sqrt(2)))*k2_vec_, control_input_vec, k3_vec_);
    model_.stateFunc(current_time+integration_length, current_state_vec-integration_length*0.5*std::sqrt(2)*k2_vec_+integration_length*(1+(1/std::sqrt(2)))*k3_vec_, control_input_vec, k4_vec_);
    next_state_vec = current_state_vec + integration_length*(k1_vec_+(2-std::sqrt(2))*k2_vec_+(2+std::sqrt(2))*k3_vec_+k4_vec_)/6;
}
