#include "zero_horizon_ocp_with_input_saturation.hpp"

ZeroHorizonOCPWithInputSaturation::ZeroHorizonOCPWithInputSaturation(
    const InputSaturationSet& input_saturation_set)
  : OptimalControlProblem(),
    input_saturation_set_(input_saturation_set),
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()
                  +2*input_saturation_set.dimSaturation()),
    dim_saturation_(input_saturation_set.dimSaturation()),
    lambda_vec_(linearalgebra::NewVector(model_.dimState())) {
}

ZeroHorizonOCPWithInputSaturation::~ZeroHorizonOCPWithInputSaturation() {
  linearalgebra::DeleteVector(lambda_vec_);
}

void ZeroHorizonOCPWithInputSaturation::computeOptimalityResidual(
    const double time, const double* state_vec,
    const double* solution_vec, double* optimality_residual) {
  model_.phixFunc(time, state_vec, lambda_vec_);
  model_.huFunc(time, state_vec, solution_vec, lambda_vec_, 
                optimality_residual);
  inputsaturationfunctions::addHamiltonianDerivativeWithSaturatedInput(
      input_saturation_set_, solution_vec, 
      &(solution_vec[dim_control_input_and_constraints_+dim_saturation_]),
      optimality_residual);
  inputsaturationfunctions::computeOptimalityResidualForDummyInput(
      input_saturation_set_, 
      &(solution_vec[dim_control_input_and_constraints_]),
      &(solution_vec[dim_control_input_and_constraints_+dim_saturation_]),
      &(optimality_residual[dim_control_input_and_constraints_]));
  inputsaturationfunctions::computeOptimalityResidualForInputSaturation(
      input_saturation_set_, solution_vec,
      &(solution_vec[dim_control_input_and_constraints_]),
      &(optimality_residual[dim_control_input_and_constraints_+dim_saturation_]));
}

void ZeroHorizonOCPWithInputSaturation::computeTerminalCostDerivative(
    const double time, const double* state_vec,
    double* terminal_cost_derivative_vec) {
  model_.phixFunc(time, state_vec, terminal_cost_derivative_vec);
}

int ZeroHorizonOCPWithInputSaturation::dim_saturation() const {
  return dim_saturation_;
}

int ZeroHorizonOCPWithInputSaturation::dim_solution() const {
  return dim_solution_;
}