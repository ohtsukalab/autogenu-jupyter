#include "ms_ocp_with_input_saturation.hpp"

MSOCPWithInputSaturation::MSOCPWithInputSaturation(
    const InputSaturationSet& input_saturation_set, const double T_f, 
    const double alpha, const int N)
  : OptimalControlProblem(),
    horizon_(T_f, alpha),
    input_saturation_set_(input_saturation_set),
    dim_solution_(N*(model_.dimControlInput()+model_.dimConstraints())),
    dim_saturation_(input_saturation_set_.dimSaturation()),
    N_(N),
    dx_vec_(linearalgebra::NewVector(model_.dimState())) {
}

MSOCPWithInputSaturation::MSOCPWithInputSaturation(
    const InputSaturationSet& input_saturation_set, const double T_f, 
    const double alpha, const int N, const double initial_time)
  : OptimalControlProblem(),
    horizon_(T_f, alpha, initial_time),
    input_saturation_set_(input_saturation_set),
    dim_solution_(N*(model_.dimControlInput()+model_.dimConstraints())),
    dim_saturation_(input_saturation_set_.dimSaturation()),
    N_(N),
    dx_vec_(linearalgebra::NewVector(model_.dimState())) {
}

MSOCPWithInputSaturation::~MSOCPWithInputSaturation() {
  linearalgebra::DeleteVector(dx_vec_);
}

void MSOCPWithInputSaturation::
computeOptimalityResidualForControlInputAndConstraints(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double const* const* input_saturation_multipler_mat,
    double* optimality_residual_for_control_input_and_constraints) {
  // Set the length of the horizon and discretize the horizon.
  double horizon_length = horizon_.getLength(time);
  double delta_tau = horizon_length / N_;
  // Compute optimality error for control input and constraints.
  // Compute optimality error for contol input and constraints.
  model_.huFunc(time, state_vec, control_input_and_constraints_seq, 
                lambda_mat[0], optimality_residual_for_control_input_and_constraints);
  inputsaturationfunctions::addHamiltonianDerivativeWithSaturatedInput(
      input_saturation_set_, control_input_and_constraints_seq, 
      input_saturation_multipler_mat[0], 
      optimality_residual_for_control_input_and_constraints
  );
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.huFunc(tau, state_mat[i-1], 
                  &(control_input_and_constraints_seq[i_total]), 
                  lambda_mat[i], 
                  &(optimality_residual_for_control_input_and_constraints[i_total]));
    inputsaturationfunctions::addHamiltonianDerivativeWithSaturatedInput(
        input_saturation_set_, &(control_input_and_constraints_seq[i_total]), 
        input_saturation_multipler_mat[i], 
        &(optimality_residual_for_control_input_and_constraints[i_total])
    );
  }
}

void MSOCPWithInputSaturation::computeOptimalityResidualForStateAndLambda(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double** optimality_residual_for_state, 
    double** optimality_residual_for_lambda) {
  // Set the length of the horizon and discretize the horizon.
  double horizon_length = horizon_.getLength(time);
  double delta_tau = horizon_length / N_;
  // Compute optimality error for state.
  model_.stateFunc(time, state_vec, control_input_and_constraints_seq, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    optimality_residual_for_state[0][i] = 
        state_mat[0][i] - state_vec[i] - delta_tau * dx_vec_[i];
  }
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.stateFunc(tau, state_mat[i-1], 
                     &(control_input_and_constraints_seq[i_total]), dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      optimality_residual_for_state[i][j] = 
          state_mat[i][j] - state_mat[i-1][j] - delta_tau * dx_vec_[j];
    }
  }
  // Compute optimality error for lambda.
  model_.phixFunc(tau, state_mat[N_-1], dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    optimality_residual_for_lambda[N_-1][i] = lambda_mat[N_-1][i] - dx_vec_[i];
  }
  for (int i=N_-1; i>=1; --i, tau-=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.hxFunc(tau, state_mat[i-1], 
                  &(control_input_and_constraints_seq[i_total]), 
                  lambda_mat[i], dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      optimality_residual_for_lambda[i-1][j] = 
          lambda_mat[i-1][j] - lambda_mat[i][j] - delta_tau * dx_vec_[j];
    }
  }
}

void MSOCPWithInputSaturation::computeStateAndLambdaFromOptimalityResidual(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* optimality_residual_for_state,
    double const* const* optimality_residual_for_lambda,
    double** state_mat, double** lambda_mat) {
  // Set the length of the horizon and discretize the horizon.
  double horizon_length = horizon_.getLength(time);
  double delta_tau = horizon_length / N_;
  // Compute the sequence of state under the error for state.
  model_.stateFunc(time, state_vec, control_input_and_constraints_seq, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    state_mat[0][i] = 
        state_vec[i] 
        + delta_tau * dx_vec_[i] + optimality_residual_for_state[0][i];
  }
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.stateFunc(tau, state_mat[i-1], 
                     &(control_input_and_constraints_seq[i_total]), dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      state_mat[i][j] = 
          state_mat[i-1][j] 
          + delta_tau * dx_vec_[j] + optimality_residual_for_state[i][j];
    }
  }
  // Compute the sequence of lambda under the error for lambda.
  model_.phixFunc(tau, state_mat[N_-1], dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    lambda_mat[N_-1][i] = dx_vec_[i] + optimality_residual_for_lambda[N_-1][i];
  }
  for (int i=N_-1; i>=1; --i, tau-=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.hxFunc(tau, state_mat[i-1], 
                  &(control_input_and_constraints_seq[i_total]), 
                  lambda_mat[i], dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat[i-1][j] = 
          lambda_mat[i][j] 
          + delta_tau * dx_vec_[j] + optimality_residual_for_lambda[i-1][j];
    }
  }
}

void MSOCPWithInputSaturation::computeResidualForDummyInputAndInputSaturation(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    double const* const* input_saturation_multiplier_mat, 
    double** errors_for_dummy_input, 
    double** errors_for_input_saturation) {
  for (int i=0; i<N_; ++i) {
    inputsaturationfunctions::computeOptimalityResidualForDummyInput(
        input_saturation_set_, dummy_input_mat[i], 
        input_saturation_multiplier_mat[i], errors_for_dummy_input[i]);
  }
  for (int i=0; i<N_; ++i) {
    inputsaturationfunctions::computeOptimalityResidualForInputSaturation(
        input_saturation_set_,
        &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]),
        dummy_input_mat[i], errors_for_input_saturation[i]);
  }
}

void MSOCPWithInputSaturation::
multiplyResidualForDummyInputAndInputSaturationInverse(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    double const* const* input_saturation_multiplier_mat, 
    double const* const* multiplied_dummy_input_mat, 
    double const* const* multiplied_Lagrange_multiplier_mat, 
    double** resulted_dummy_input_mat, 
    double** resulted_Lagrange_multiplier_mat) {
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      resulted_dummy_input_mat[i][j] = 
          multiplied_Lagrange_multiplier_mat[i][j] / (2*dummy_input_mat[i][j]);
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      resulted_Lagrange_multiplier_mat[i][j] = 
          multiplied_dummy_input_mat[i][j] / (2*dummy_input_mat[i][j]) 
          - ((input_saturation_multiplier_mat[i][j]
                  +input_saturation_set_.quadratic_weight(j)) 
              *resulted_dummy_input_mat[i][j]) / dummy_input_mat[i][j];
    }
  }
}

void MSOCPWithInputSaturation::computeResidualDifferenceForDummyInput(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    const double* control_input_and_constraints_update_seq, 
    double** dummy_residual_difference_mat) {
  for (int i=0; i<N_; ++i) { 
    int i_total = i * dim_control_input_and_constraints_;
    for (int j=0; j<dim_saturation_; ++j) {
      int index_j = input_saturation_set_.index(j);
      dummy_residual_difference_mat[i][j] = 
          ((2*control_input_and_constraints_seq[i_total+index_j] 
              -input_saturation_set_.min(j)-input_saturation_set_.max(j))
              *control_input_and_constraints_update_seq[i_total+index_j]) 
          / (2*dummy_input_mat[i][j]);
    }
  }
}

void MSOCPWithInputSaturation::computeResidualDifferenceForInputSaturation(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    double const* const* input_saturation_multiplier_mat, 
    const double* control_input_and_constraints_update_seq, 
    double** input_saturation_difference_mat) {
  for (int i=0; i<N_; ++i) {
    int i_total = i * dim_control_input_and_constraints_;
    for (int j=0; j<dim_saturation_; ++j) {
      int index_j = input_saturation_set_.index(j);
      input_saturation_difference_mat[i][j] = 
          - ((input_saturation_multiplier_mat[i][j]
                +input_saturation_set_.quadratic_weight(j))
          *(2*control_input_and_constraints_seq[i_total+index_j]
                -input_saturation_set_.min(j)-input_saturation_set_.max(j))
          *control_input_and_constraints_update_seq[i_total+index_j]) 
          / (2*dummy_input_mat[i][j]*dummy_input_mat[i][j]);
    }
  }
}

void MSOCPWithInputSaturation::predictStateFromSolution(
    const double current_time, const double* current_state,
    const double* solution_vec, const double prediction_length,
    double* predicted_state) {
  model_.stateFunc(current_time, current_state, solution_vec, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    predicted_state[i] =  current_state[i] + prediction_length * dx_vec_[i];
  }
}

void MSOCPWithInputSaturation::resetHorizonLength(const double initial_time) {
  horizon_.resetLength(initial_time);
}

void MSOCPWithInputSaturation::resetHorizonLength(const double T_f, 
                                                  const double alpha, 
                                                  const double initial_time) {
  horizon_.resetLength(T_f, alpha, initial_time);
}

int MSOCPWithInputSaturation::dim_solution() const {
  return dim_solution_;
}

int MSOCPWithInputSaturation::dim_saturation() const {
  return dim_saturation_;
}

int MSOCPWithInputSaturation::N() const {
  return N_;
}