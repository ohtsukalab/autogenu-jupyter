#include "zero_horizon_ocp.hpp"

ZeroHorizonOCP::ZeroHorizonOCP() 
  : OptimalControlProblem(),
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()),
    lambda_vec_(linearalgebra::NewVector(dim_state_)) {
}

ZeroHorizonOCP::~ZeroHorizonOCP() {
  linearalgebra::DeleteVector(lambda_vec_);
}

void ZeroHorizonOCP::computeOptimalityResidual(
    const double time, const double* state_vec,
    const double* solution_vec, double* optimality_residual) {
  model_.phixFunc(time, state_vec, lambda_vec_);
  model_.huFunc(time, state_vec, solution_vec, lambda_vec_, optimality_residual);
}

void ZeroHorizonOCP::computeTerminalCostDerivative(
    const double time, const double* state_vec,
    double* terminal_cost_derivative_vec) {
  model_.phixFunc(time, state_vec, terminal_cost_derivative_vec);
}

int ZeroHorizonOCP::dim_solution() const {
  return dim_solution_;
}