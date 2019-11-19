#include "single_shooting_ocp.hpp"

SingleShootingOCP::SingleShootingOCP(const double T_f, const double alpha, 
                                     const int N) 
  : OptimalControlProblem(),
    horizon_(T_f, alpha),
    dim_solution_(N*(model_.dimControlInput()+model_.dimConstraints())),
    N_(N),
    dx_vec_(linearalgebra::NewVector(model_.dimState())),
    state_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())),
    lambda_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())) {
}

SingleShootingOCP::SingleShootingOCP(const double T_f, const double alpha, 
                                     const int N, const double initial_time) 
  : OptimalControlProblem(),
    horizon_(T_f, alpha, initial_time),
    dim_solution_(N*(model_.dimControlInput()+model_.dimConstraints())),
    N_(N),
    dx_vec_(linearalgebra::NewVector(model_.dimState())),
    state_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())),
    lambda_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())) {
}

SingleShootingOCP::~SingleShootingOCP() {
  linearalgebra::DeleteVector(dx_vec_);
  linearalgebra::DeleteMatrix(state_mat_);
  linearalgebra::DeleteMatrix(lambda_mat_);
}

void SingleShootingOCP::computeOptimalityResidual(const double time, 
                                                  const double* state_vec, 
                                                  const double* solution_vec,
                                                  double* optimality_residual) {
  double horizon_length = horizon_.getLength(time);
  double delta_tau = horizon_length / N_;
  // Compute the state trajectory over the horizon on the basis of the 
  // time, solution_vec and the state_vec.
  model_.stateFunc(time, state_vec, solution_vec, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    state_mat_[1][i] = state_vec[i] + delta_tau * dx_vec_[i];
  }
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    model_.stateFunc(
        tau, state_mat_[i], 
        &(solution_vec[i*dim_control_input_and_constraints_]), dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      state_mat_[i+1][j] = state_mat_[i][j] + delta_tau * dx_vec_[j];
    }
  }
  // Compute the Lagrange multiplier over the horizon on the basis of 
  // time, solution_vec and the state_vec.
  model_.phixFunc(tau, state_mat_[N_], lambda_mat_[N_]);
  for (int i=N_-1; i>=1; --i, tau-=delta_tau) {
    model_.hxFunc(
        tau, state_mat_[i], 
        &(solution_vec[i*dim_control_input_and_constraints_]), 
        lambda_mat_[i+1], dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat_[i][j] = lambda_mat_[i+1][j] + delta_tau * dx_vec_[j];
    }
  }
  // Compute the erros in optimality over the horizon on the basis of the 
  // control_input_vec and the state_vec.
  model_.huFunc(time, state_vec, solution_vec, lambda_mat_[1], 
                optimality_residual);
  tau = time;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    model_.huFunc(
        tau, state_mat_[i], 
        &(solution_vec[i*dim_control_input_and_constraints_]), 
        lambda_mat_[i+1], 
        &(optimality_residual[i*dim_control_input_and_constraints_]));
  }
}

void SingleShootingOCP::predictStateFromSolution(const double current_time, 
                                                 const double* current_state,
                                                 const double* solution_vec, 
                                                 const double prediction_length,
                                                 double* predicted_state) {
  model_.stateFunc(current_time, current_state, solution_vec, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    predicted_state[i] =  current_state[i] + prediction_length * dx_vec_[i];
  }
}

void SingleShootingOCP::resetHorizonLength(const double T_f, const double alpha, 
                                           const double initial_time) {
  horizon_.resetLength(T_f, alpha, initial_time);
}

void SingleShootingOCP::resetHorizonLength(const double initial_time) {
  horizon_.resetLength(initial_time);
}

int SingleShootingOCP::dim_solution() const {
  return dim_solution_;
}

int SingleShootingOCP::N() const {
  return N_;
}