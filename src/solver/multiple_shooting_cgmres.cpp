#include "multiple_shooting_cgmres.hpp"


MultipleShootingCGMRES::MultipleShootingCGMRES(
    const double T_f, const double alpha, const int N, const double zeta, 
    const double finite_difference_step, const int kmax) 
  : MatrixFreeGMRES(), 
    model_(), 
    cgmres_initializer_(),
    dim_state_(model_.dimState()), 
    dim_control_input_(model_.dimControlInput()), 
    dim_constraints_(model_.dimConstraints()), 
    dim_control_input_and_constraints_(
        model_.dimControlInput()+model_.dimConstraints()), 
    dim_state_and_lambda_(2*model_.dimState()), 
    dim_control_input_and_constraints_seq_(
        N*dim_control_input_and_constraints_),
    T_f_(T_f), 
    alpha_(alpha), 
    N_(N), 
    zeta_(zeta), 
    finite_difference_step_(finite_difference_step), 
    kmax_(kmax), 
    initial_time_(0), 
    incremented_time_(0), 
    dx_vec_(linearalgebra::NewVector(dim_state_)), 
    incremented_state_vec_(linearalgebra::NewVector(dim_state_)), 
    control_input_and_constraints_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)), 
    incremented_control_input_and_constraints_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_1_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_2_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_3_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_update_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)), 
    state_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    lambda_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    incremented_state_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    incremented_lambda_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    state_error_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    state_error_mat_1_(linearalgebra::NewMatrix(N, dim_state_)), 
    lambda_error_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    lambda_error_mat_1_(linearalgebra::NewMatrix(N, dim_state_)) {
  // Set dimensions and parameters in GMRES.
  setGMRESParameters(dim_control_input_and_constraints_seq_, kmax);
}

MultipleShootingCGMRES::~MultipleShootingCGMRES() {
  linearalgebra::DeleteVector(dx_vec_);
  linearalgebra::DeleteVector(incremented_state_vec_);
  linearalgebra::DeleteVector(control_input_and_constraints_seq_);
  linearalgebra::DeleteVector(incremented_control_input_and_constraints_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_1_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_2_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_3_);
  linearalgebra::DeleteVector(control_input_and_constraints_update_seq_);
  linearalgebra::DeleteMatrix(state_mat_);
  linearalgebra::DeleteMatrix(lambda_mat_);
  linearalgebra::DeleteMatrix(incremented_state_mat_);
  linearalgebra::DeleteMatrix(incremented_lambda_mat_);
  linearalgebra::DeleteMatrix(state_error_mat_);
  linearalgebra::DeleteMatrix(state_error_mat_1_);
  linearalgebra::DeleteMatrix(lambda_error_mat_);
  linearalgebra::DeleteMatrix(lambda_error_mat_1_);
}

void MultipleShootingCGMRES::setInitParameters(
    const double* initial_guess_solution, const double residual_tolerance, 
    const int max_iteration, const double finite_difference_step, 
    const int kmax) {
  cgmres_initializer_.setInitParameters(initial_guess_solution, 
                                        residual_tolerance, max_iteration, 
                                        finite_difference_step, kmax);
}

void MultipleShootingCGMRES::initSolution(const double initial_time, 
                                          const double* initial_state_vec, 
                                          double* optimal_control_input_vec) {
  double initial_solution_vec[dim_control_input_and_constraints_], 
      initial_errors_in_optimality[dim_control_input_and_constraints_], 
      initial_lambda_vec[dim_state_];
  initial_time_ = initial_time;
  cgmres_initializer_.solveOCPForInit(initial_time, initial_state_vec, 
                                      initial_solution_vec, 
                                      initial_errors_in_optimality);
  model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
  for (int i=0; i<N_; ++i) {
    int i_total = i * dim_control_input_and_constraints_;
    for (int j=0; j<dim_control_input_and_constraints_; ++j) {
      control_input_and_constraints_seq_[i_total+j] = initial_solution_vec[j];
    }
    for (int j=0; j<dim_state_; ++j) {
      state_mat_[i][j] = initial_state_vec[j];
    }
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat_[i][j] = initial_lambda_vec[j];
    }
  }
  for (int i=0; i<N_; ++i) {
    int i_total = i * dim_control_input_and_constraints_;
    for (int j=0; j<dim_control_input_and_constraints_; ++j) {
      control_input_and_constraints_error_seq_[i_total+j] 
          = initial_errors_in_optimality[j];
    }
  }
  for (int i=0; i<dim_control_input_; ++i) {
    optimal_control_input_vec[i] = initial_solution_vec[i];
  }
}

void MultipleShootingCGMRES::controlUpdate(const double time, 
                                           const double sampling_period, 
                                           const double* state_vec, 
                                           double* optimal_control_input_vec) {
  // Predict the incremented state.
  incremented_time_ = time + finite_difference_step_;
  model_.stateFunc(time, state_vec, control_input_and_constraints_seq_, 
                   dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    incremented_state_vec_[i] = state_vec[i] 
        + finite_difference_step_ * dx_vec_[i];
  }
  solveGMRES(time, state_vec, control_input_and_constraints_seq_, 
             control_input_and_constraints_update_seq_);
  // Update state_mat_ and lamdba_mat_ by the difference approximation.
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i] = 
        control_input_and_constraints_seq_[i] 
        + finite_difference_step_ 
        * control_input_and_constraints_update_seq_[i];
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      state_error_mat_1_[i][j] = 
          (1-finite_difference_step_*zeta_) * state_error_mat_[i][j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) { 
      lambda_error_mat_1_[i][j] = 
          (1-finite_difference_step_*zeta_) * lambda_error_mat_[i][j];
    }
  }
  computeStateAndLambdaFromErrors(
        incremented_time_, incremented_state_vec_, 
        incremented_control_input_and_constraints_seq_, state_error_mat_1_, 
        lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
  // state_mat_ += 
  //     (sampling_period/finite_difference_step_) 
  //     * (incremented_state_mat_-state_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      state_mat_[i][j] += 
          (sampling_period/finite_difference_step_) 
          * (incremented_state_mat_[i][j]-state_mat_[i][j]);
    }
  }
  // lambda_mat_ += 
  //     (sampling_period/finite_difference_step_) 
  //     * (incremented_lambda_mat_-lambda_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat_[i][j] += 
          (sampling_period/finite_difference_step_) 
          * (incremented_lambda_mat_[i][j]-lambda_mat_[i][j]);
    }
  }
  // Update control_input_and_constraints_seq_
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    control_input_and_constraints_seq_[i] += sampling_period 
        * control_input_and_constraints_update_seq_[i];
  }
  for (int i=0; i<dim_control_input_; ++i) {
    optimal_control_input_vec[i] = control_input_and_constraints_seq_[i];
  }
}

double MultipleShootingCGMRES::getErrorNorm(const double time, 
                                            const double* state_vec) {
  computeErrorsForControlInputAndConstraints(
      time, state_vec, control_input_and_constraints_seq_, state_mat_, 
      lambda_mat_, control_input_and_constraints_error_seq_);
  computeErrorsForStateAndLambda(
      time, state_vec, control_input_and_constraints_seq_, state_mat_, 
      lambda_mat_, state_error_mat_, lambda_error_mat_);
  double Squared_error = linearalgebra::SquaredNorm(
      dim_control_input_and_constraints_seq_, 
      control_input_and_constraints_error_seq_);
  for (int i=0; i<N_; ++i) {
    Squared_error += 
        (linearalgebra::SquaredNorm(dim_state_, state_error_mat_[i]) 
        + linearalgebra::SquaredNorm(dim_state_, lambda_error_mat_[i]));
  }
  return std::sqrt(Squared_error);
}

inline void MultipleShootingCGMRES::computeErrorsForControlInputAndConstraints(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double* errors_for_control_input_and_constraints) {
  // Set the length of the horizon and discretize the horizon.
  double horizon_length = T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
  double delta_tau = horizon_length / N_;
  // Compute optimality error for control input and constraints.
  model_.huFunc(time, state_vec, control_input_and_constraints_seq, 
                lambda_mat[0], errors_for_control_input_and_constraints);
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.huFunc(tau, state_mat[i-1], 
                  &(control_input_and_constraints_seq[i_total]), 
                  lambda_mat[i], 
                  &(errors_for_control_input_and_constraints[i_total]));
  }
}

inline void MultipleShootingCGMRES::computeErrorsForStateAndLambda(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double** errors_for_state, double** errors_for_lambda) {
  // Set the length of the horizon and discretize the horizon.
  double horizon_length = T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
  double delta_tau = horizon_length / N_;
  // Compute optimality error for state.
  model_.stateFunc(time, state_vec, control_input_and_constraints_seq, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    errors_for_state[0][i] = 
        state_mat[0][i] - state_vec[i] - delta_tau * dx_vec_[i];
  }
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.stateFunc(tau, state_mat[i-1], 
                     &(control_input_and_constraints_seq[i_total]), dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      errors_for_state[i][j] = 
          state_mat[i][j] - state_mat[i-1][j] - delta_tau * dx_vec_[j];
    }
  }
  // Compute optimality error for lambda.
  model_.phixFunc(tau, state_mat[N_-1], dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    errors_for_lambda[N_-1][i] = lambda_mat[N_-1][i] - dx_vec_[i];
  }
  for (int i=N_-1; i>=1; --i, tau-=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.hxFunc(tau, state_mat[i-1], 
                  &(control_input_and_constraints_seq[i_total]), 
                  lambda_mat[i], dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      errors_for_lambda[i-1][j] = 
          lambda_mat[i-1][j] - lambda_mat[i][j] - delta_tau * dx_vec_[j];
    }
  }
}

inline void MultipleShootingCGMRES::computeStateAndLambdaFromErrors(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* errors_for_state, 
    double const* const* errors_for_lambda, double** state_mat, 
    double** lambda_mat) {
  // Set the length of the horizon and discretize the horizon.
  double horizon_length = T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
  double delta_tau = horizon_length / N_;
  // Compute the sequence of state under the error for state.
  model_.stateFunc(time, state_vec, control_input_and_constraints_seq, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    state_mat[0][i] = 
        state_vec[i] + delta_tau*dx_vec_[i] + errors_for_state[0][i];
  }
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.stateFunc(tau, state_mat[i-1], 
                     &(control_input_and_constraints_seq[i_total]), dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      state_mat[i][j] = 
          state_mat[i-1][j] + delta_tau*dx_vec_[j] + errors_for_state[i][j];
    }
  }
  // Compute the sequence of lambda under the error for lambda.
  model_.phixFunc(tau, state_mat[N_-1], dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    lambda_mat[N_-1][i] = dx_vec_[i] + errors_for_lambda[N_-1][i];
  }
  for (int i=N_-1; i>=1; --i, tau-=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.hxFunc(tau, state_mat[i-1], 
                  &(control_input_and_constraints_seq[i_total]), 
                  lambda_mat[i], dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat[i-1][j] = 
          lambda_mat[i][j] + delta_tau*dx_vec_[j] + errors_for_lambda[i-1][j];
    }
  }
}

void MultipleShootingCGMRES::bFunc(const double time, const double* state_vec, 
                                   const double* current_solution_vec, 
                                   double* b_vec) {
  computeErrorsForControlInputAndConstraints(
        time, state_vec, current_solution_vec, state_mat_, lambda_mat_, 
        control_input_and_constraints_error_seq_);
  computeErrorsForControlInputAndConstraints(
        incremented_time_, incremented_state_vec_, current_solution_vec, 
        state_mat_, lambda_mat_, control_input_and_constraints_error_seq_1_);
  computeErrorsForStateAndLambda(
        time, state_vec, current_solution_vec, state_mat_, lambda_mat_, 
        state_error_mat_, lambda_error_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      state_error_mat_1_[i][j] = 
          (1-finite_difference_step_*zeta_) * state_error_mat_[i][j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      lambda_error_mat_1_[i][j] = 
          (1-finite_difference_step_*zeta_) * lambda_error_mat_[i][j];
    }
  }
  computeStateAndLambdaFromErrors(
      incremented_time_, incremented_state_vec_, current_solution_vec, 
      state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, 
      incremented_lambda_mat_);
  computeErrorsForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, current_solution_vec, 
      incremented_state_mat_, incremented_lambda_mat_, 
      control_input_and_constraints_error_seq_3_);
  computeErrorsForStateAndLambda(
      incremented_time_, incremented_state_vec_, current_solution_vec, 
      state_mat_, lambda_mat_, state_error_mat_1_, lambda_error_mat_1_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i] = 
      current_solution_vec[i] 
      + finite_difference_step_ * control_input_and_constraints_update_seq_[i];
  }
  computeStateAndLambdaFromErrors(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, state_error_mat_1_, 
      lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);    
  computeErrorsForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, incremented_state_mat_, 
      incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    b_vec[i] = 
        (1/finite_difference_step_-zeta_) 
        * control_input_and_constraints_error_seq_[i] 
        - control_input_and_constraints_error_seq_3_[i] / finite_difference_step_ 
        - (control_input_and_constraints_error_seq_2_[i]
           -control_input_and_constraints_error_seq_1_[i])
        / finite_difference_step_;
  }
}


void MultipleShootingCGMRES::axFunc(const double time, const double* state_vec, 
                                    const double* current_solution_vec, 
                                    const double* direction_vec, 
                                    double* ax_vec) {
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i] = current_solution_vec[i] 
        + finite_difference_step_*direction_vec[i];
  }
  computeStateAndLambdaFromErrors(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, state_error_mat_1_, 
      lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
  computeErrorsForControlInputAndConstraints(
        incremented_time_, incremented_state_vec_, 
        incremented_control_input_and_constraints_seq_, incremented_state_mat_, 
        incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    ax_vec[i] = (control_input_and_constraints_error_seq_2_[i]
        -control_input_and_constraints_error_seq_1_[i]) / finite_difference_step_;
  }
}