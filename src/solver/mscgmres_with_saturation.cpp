#include "mscgmres_with_saturation.hpp"


MSCGMRESWithSaturation::MSCGMRESWithSaturation(
    const ControlInputSaturationSequence saturation_seq, const double T_f, 
    const double alpha, const int N, const double zeta, 
    const double finite_difference_step, const int kmax)
  : MatrixFreeGMRES(),
    model_(), 
    saturation_seq_(saturation_seq), 
    mscgmres_initializer_(saturation_seq),
    dim_state_(model_.dimState()), 
    dim_control_input_(model_.dimControlInput()), 
    dim_constraints_(model_.dimConstraints()), 
    dim_control_input_and_constraints_(
        model_.dimControlInput()+model_.dimConstraints()),
    dim_control_input_and_constraints_seq_(
        N*(model_.dimControlInput()+model_.dimConstraints())), 
    dim_saturation_(saturation_seq.dimSaturation()), 
    dim_saturation_seq_(N*saturation_seq.dimSaturation()), 
    N_(N), 
    initial_time_(0), 
    T_f_(T_f), 
    alpha_(alpha), 
    zeta_(zeta), 
    finite_difference_step_(finite_difference_step), 
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
    state_mat_1_(linearalgebra::NewMatrix(N, dim_state_)), 
    lambda_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    lambda_mat_1_(linearalgebra::NewMatrix(N, dim_state_)), 
    incremented_state_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    incremented_lambda_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    state_error_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    state_error_mat_1_(linearalgebra::NewMatrix(N, dim_state_)), 
    lambda_error_mat_(linearalgebra::NewMatrix(N, dim_state_)), 
    lambda_error_mat_1_(linearalgebra::NewMatrix(N, dim_state_)), 
    dummy_input_mat_(linearalgebra::NewMatrix(N, dim_saturation_)), 
    dummy_input_mat_1_(linearalgebra::NewMatrix(N, dim_saturation_)), 
    saturation_Lagrange_multiplier_mat_(
        linearalgebra::NewMatrix(N, dim_saturation_)), 
    saturation_Lagrange_multiplier_mat_1_(
        linearalgebra::NewMatrix(N, dim_saturation_)), 
    incremented_saturation_Lagrange_multiplier_mat_(
        linearalgebra::NewMatrix(N, dim_saturation_)), 
    dummy_error_mat_(linearalgebra::NewMatrix(N, dim_saturation_)), 
    dummy_error_mat_1_(linearalgebra::NewMatrix(N, dim_saturation_)), 
    saturation_error_mat_(linearalgebra::NewMatrix(N, dim_saturation_)), 
    saturation_error_mat_1_(linearalgebra::NewMatrix(N, dim_saturation_)), 
    dummy_update_mat_(linearalgebra::NewMatrix(N, dim_saturation_)), 
    saturation_update_mat_(linearalgebra::NewMatrix(N, dim_saturation_)) {
  // Set dimensions and parameters in GMRES.
  setGMRESParameters(dim_control_input_and_constraints_seq_, kmax);
}

MSCGMRESWithSaturation::~MSCGMRESWithSaturation() {
  linearalgebra::DeleteVector(dx_vec_);
  linearalgebra::DeleteVector(incremented_state_vec_);
  linearalgebra::DeleteVector(control_input_and_constraints_seq_);
  linearalgebra::DeleteVector(
      incremented_control_input_and_constraints_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_1_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_2_);
  linearalgebra::DeleteVector(control_input_and_constraints_error_seq_3_);
  linearalgebra::DeleteVector(control_input_and_constraints_update_seq_);
  linearalgebra::DeleteMatrix(state_mat_);
  linearalgebra::DeleteMatrix(state_mat_1_);
  linearalgebra::DeleteMatrix(lambda_mat_);
  linearalgebra::DeleteMatrix(lambda_mat_1_);
  linearalgebra::DeleteMatrix(incremented_state_mat_);
  linearalgebra::DeleteMatrix(incremented_lambda_mat_);
  linearalgebra::DeleteMatrix(state_error_mat_);
  linearalgebra::DeleteMatrix(state_error_mat_1_);
  linearalgebra::DeleteMatrix(lambda_error_mat_);
  linearalgebra::DeleteMatrix(lambda_error_mat_1_);
  linearalgebra::DeleteMatrix(dummy_input_mat_);
  linearalgebra::DeleteMatrix(dummy_input_mat_1_);
  linearalgebra::DeleteMatrix(saturation_Lagrange_multiplier_mat_);
  linearalgebra::DeleteMatrix(saturation_Lagrange_multiplier_mat_1_);
  linearalgebra::DeleteMatrix(
      incremented_saturation_Lagrange_multiplier_mat_);
  linearalgebra::DeleteMatrix(dummy_error_mat_);
  linearalgebra::DeleteMatrix(dummy_error_mat_1_);
  linearalgebra::DeleteMatrix(saturation_error_mat_);
  linearalgebra::DeleteMatrix(saturation_error_mat_1_);
  linearalgebra::DeleteMatrix(dummy_update_mat_);
  linearalgebra::DeleteMatrix(saturation_update_mat_);
}

void MSCGMRESWithSaturation::setInitParameters(
    const double* initial_guess_solution, 
    const double* initial_guess_Lagrange_multiplier, 
    const double residual_tolerance, const int max_iteration, 
    const double finite_difference_step, const int kmax) {
  mscgmres_initializer_.setInitParameters(initial_guess_solution, 
                                          initial_guess_Lagrange_multiplier, 
                                          residual_tolerance, max_iteration, 
                                          finite_difference_step, kmax);
}

void MSCGMRESWithSaturation::initSolution(const double initial_time, 
                                          const double* initial_state_vec, 
                                          double* optimal_control_input_vec) {
  double initial_solution_vec[dim_control_input_and_constraints_+2*dim_saturation_], 
      initial_lambda_vec[dim_state_], 
      initial_error_for_control_input_and_constraint[dim_control_input_and_constraints_], 
      initial_error_for_dummy_input[dim_saturation_], 
      initial_error_for_saturation[dim_saturation_];
  initial_time_ = initial_time;
  mscgmres_initializer_.solveOCPForInit(
      initial_time, initial_state_vec, initial_solution_vec, 
      initial_error_for_control_input_and_constraint, 
      initial_error_for_dummy_input, initial_error_for_saturation);
  model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_control_input_and_constraints_; ++j) {
      control_input_and_constraints_seq_[j] = initial_solution_vec[j];
    }
    for (int j=0; j<dim_saturation_; ++j) {
      dummy_input_mat_[i][j] = 
          initial_solution_vec[dim_control_input_and_constraints_+j];
    }
    for (int j=0; j<dim_saturation_; ++j) {
      saturation_Lagrange_multiplier_mat_[i][j] = 
          initial_solution_vec[dim_control_input_and_constraints_+dim_saturation_+j];
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
      control_input_and_constraints_error_seq_[i_total+j] = 
          initial_error_for_control_input_and_constraint[j];
    }
    for (int j=0; j<dim_saturation_; ++j) {
      dummy_error_mat_[i][j] = initial_error_for_dummy_input[j];
    }
    for (int j=0; j<dim_saturation_; ++j) {
      saturation_error_mat_[i][j] = initial_error_for_saturation[j];
    }
  }
  for (int i=0; i<dim_control_input_; ++i) {
    optimal_control_input_vec[i] = initial_solution_vec[i];
  }
}

void MSCGMRESWithSaturation::controlUpdate(const double time, 
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
  // Update dummy_input_mat_ and saturation_Lagrange_multiplier_mat_.
  computeErrorsDifferenceForDummyInput(
      control_input_and_constraints_seq_, dummy_input_mat_, 
      control_input_and_constraints_update_seq_, dummy_update_mat_);
  computeErrorsDifferenceForSaturation(
      control_input_and_constraints_seq_, dummy_input_mat_, 
      saturation_Lagrange_multiplier_mat_, 
      control_input_and_constraints_update_seq_, saturation_update_mat_);
  for (int i=0; i<N_; ++i) { 
    for (int j=0; j<dim_saturation_; ++j) {
      dummy_input_mat_[i][j] += 
          sampling_period * (dummy_error_mat_1_[i][j]-dummy_update_mat_[i][j]);
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      saturation_Lagrange_multiplier_mat_[i][j] += 
          sampling_period 
          * (saturation_error_mat_1_[i][j]-saturation_update_mat_[i][j]);
    }
  }
  // Update control_input_and_constraints_seq_.
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    control_input_and_constraints_seq_[i] += 
        sampling_period * control_input_and_constraints_update_seq_[i];
  }
  for (int i=0; i<dim_control_input_; ++i) {
    optimal_control_input_vec[i] = control_input_and_constraints_seq_[i];
  }
}

double MSCGMRESWithSaturation::getErrorNorm(const double time, 
                                            const double* state_vec) {
  // Computes errors in optimality for control input and constraints.
  double *errors_for_control_input_and_constraints;
  errors_for_control_input_and_constraints = 
      linearalgebra::NewVector(dim_control_input_and_constraints_seq_);
  computeErrorsForControlInputAndConstraints(
      time, state_vec, control_input_and_constraints_seq_, state_mat_, 
      lambda_mat_, saturation_Lagrange_multiplier_mat_, 
      errors_for_control_input_and_constraints);
  // Computes errors in optimality for state and Lagrange multiplier for 
  // state equation.
  double **errors_for_state_mat, **errors_for_lambda_mat;
  errors_for_state_mat = linearalgebra::NewMatrix(N_, dim_state_);
  errors_for_lambda_mat = linearalgebra::NewMatrix(N_, dim_state_);
  computeErrorsForStateAndLambda(
      time, state_vec, control_input_and_constraints_seq_, state_mat_, 
      lambda_mat_, errors_for_state_mat, errors_for_lambda_mat);
  // Computes errors in optimality for dummy input, constraints about the 
  // control input saturation.
  double **errors_for_dummy_input, **errors_for_input_saturation;
  errors_for_dummy_input = linearalgebra::NewMatrix(N_, dim_saturation_);
  errors_for_input_saturation = linearalgebra::NewMatrix(N_, dim_saturation_);
  computeErrorsForDummyInputAndSaturation(
      control_input_and_constraints_seq_, dummy_input_mat_, 
      saturation_Lagrange_multiplier_mat_, errors_for_dummy_input, 
      errors_for_input_saturation);
  // Computes total error norm.
  double squared_error_norm = linearalgebra::SquaredNorm(
      dim_control_input_and_constraints_seq_, 
      errors_for_control_input_and_constraints);    
  for (int i=0; i<N_; ++i) {
    squared_error_norm += 
        (linearalgebra::SquaredNorm(
            dim_state_, errors_for_state_mat[i]) 
        +linearalgebra::SquaredNorm(
            dim_state_, errors_for_lambda_mat[i])
        +linearalgebra::SquaredNorm(
            dim_saturation_, errors_for_dummy_input[i])
        +linearalgebra::SquaredNorm(
            dim_saturation_, errors_for_input_saturation[i]));
  }
  linearalgebra::DeleteVector(errors_for_control_input_and_constraints);
  linearalgebra::DeleteMatrix(errors_for_state_mat);
  linearalgebra::DeleteMatrix(errors_for_lambda_mat);
  linearalgebra::DeleteMatrix(errors_for_dummy_input);
  linearalgebra::DeleteMatrix(errors_for_input_saturation);
  return std::sqrt(squared_error_norm);
}

inline void MSCGMRESWithSaturation::computeErrorsForControlInputAndConstraints(
      const double time, const double* state_vec, 
      const double* control_input_and_constraints_seq, 
      double const* const* state_mat, double const* const* lambda_mat, 
      double const* const* saturation_Lagrange_multiplier_mat,
      double* errors_for_control_input_and_constraints) {
  // Set and discretize the horizon.
  double horizon_length = T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
  double delta_tau = horizon_length / N_;
  // Compute optimality error for contol input and constraints.
  model_.huFunc(time, state_vec, control_input_and_constraints_seq, 
                lambda_mat[0], errors_for_control_input_and_constraints);
  condensingfunctions::addHamiltonianDerivativeWithConstrainedInput(
      saturation_seq_, control_input_and_constraints_seq, 
      saturation_Lagrange_multiplier_mat[0], 
      errors_for_control_input_and_constraints
  );
  double tau = time + delta_tau;
  for (int i=1; i<N_; ++i, tau+=delta_tau) {
    int i_total = i * dim_control_input_and_constraints_;
    model_.huFunc(tau, state_mat[i-1], 
                  &(control_input_and_constraints_seq[i_total]), 
                  lambda_mat[i], 
                  &(errors_for_control_input_and_constraints[i_total]));
    condensingfunctions::addHamiltonianDerivativeWithConstrainedInput(
        saturation_seq_, &(control_input_and_constraints_seq[i_total]), 
        saturation_Lagrange_multiplier_mat[i], 
        &(errors_for_control_input_and_constraints[i_total])
    );
  }
}

inline void MSCGMRESWithSaturation::computeErrorsForStateAndLambda(
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

inline void MSCGMRESWithSaturation::computeStateAndLambdaFromErrors(
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

inline void MSCGMRESWithSaturation::computeErrorsForDummyInputAndSaturation(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    double const* const* saturation_Lagrange_multiplier_mat, 
    double** errors_for_dummy_input, 
    double** errors_for_input_saturation) {
  for (int i=0; i<N_; ++i) {
    condensingfunctions::computeErrorsForDummyInput(
        saturation_seq_, dummy_input_mat[i], 
        saturation_Lagrange_multiplier_mat[i], errors_for_dummy_input[i]);
  }
  for (int i=0; i<N_; ++i) {
    condensingfunctions::computeErrorsForSaturation(
        saturation_seq_,
        &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]),
        dummy_input_mat[i], errors_for_input_saturation[i]);
  }
}

inline void MSCGMRESWithSaturation::multiplyErrorsForInputAndSaturationInverse(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    double const* const* saturation_Lagrange_multiplier_mat, 
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
          - ((saturation_Lagrange_multiplier_mat[i][j]
                  +saturation_seq_.quadratic_weight(j)) 
              *resulted_dummy_input_mat[i][j]) / dummy_input_mat[i][j];
    }
  }
}

inline void MSCGMRESWithSaturation::computeErrorsDifferenceForDummyInput(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    const double* control_input_and_constraints_update_seq, 
    double** dummy_difference_mat) {
  for (int i=0; i<N_; ++i) { 
    int i_total = i * dim_control_input_and_constraints_;
    for (int j=0; j<dim_saturation_; ++j) {
      int index_j = saturation_seq_.index(j);
      dummy_difference_mat[i][j] = 
          ((2*control_input_and_constraints_seq[i_total+index_j] 
              -saturation_seq_.min(j)-saturation_seq_.max(j))
              *control_input_and_constraints_update_seq[i_total+index_j]) 
          / (2*dummy_input_mat[i][j]);
    }
  }
}

inline void MSCGMRESWithSaturation::computeErrorsDifferenceForSaturation(
    const double* control_input_and_constraints_seq, 
    double const* const* dummy_input_mat, 
    double const* const* saturation_Lagrange_multiplier_mat, 
    const double* control_input_and_constraints_update_seq, 
    double** saturation_difference_mat) {
  for (int i=0; i<N_; ++i) {
    int i_total = i * dim_control_input_and_constraints_;
    for (int j=0; j<dim_saturation_; ++j) {
      int index_j = saturation_seq_.index(j);
      saturation_difference_mat[i][j] = 
          - ((saturation_Lagrange_multiplier_mat[i][j]
                +saturation_seq_.quadratic_weight(j))
          *(2*control_input_and_constraints_seq[i_total+index_j]
                -saturation_seq_.min(j)-saturation_seq_.max(j))
          *control_input_and_constraints_update_seq[i_total+index_j]) 
          / (2*dummy_input_mat[i][j]*dummy_input_mat[i][j]);
    }
  }
}

void MSCGMRESWithSaturation::bFunc(const double time, const double* state_vec, 
                                   const double* current_solution_vec, 
                                   double* b_vec) {
  computeErrorsForControlInputAndConstraints(
      time, state_vec, current_solution_vec, state_mat_, lambda_mat_, 
      saturation_Lagrange_multiplier_mat_, 
      control_input_and_constraints_error_seq_);
  computeErrorsForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, current_solution_vec, 
      state_mat_, lambda_mat_, saturation_Lagrange_multiplier_mat_, 
      control_input_and_constraints_error_seq_1_);
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
  computeErrorsForStateAndLambda(
      incremented_time_, incremented_state_vec_, current_solution_vec, 
      state_mat_, lambda_mat_, state_error_mat_1_, lambda_error_mat_1_);
  computeErrorsForDummyInputAndSaturation(
      current_solution_vec, dummy_input_mat_, 
      saturation_Lagrange_multiplier_mat_, dummy_error_mat_, 
      saturation_error_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      dummy_input_mat_1_[i][j] = - zeta_ * dummy_error_mat_[i][j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      saturation_Lagrange_multiplier_mat_1_[i][j] = 
          - zeta_ * saturation_error_mat_[i][j];
    }
  }
  multiplyErrorsForInputAndSaturationInverse(
      current_solution_vec, dummy_input_mat_, 
      saturation_Lagrange_multiplier_mat_, dummy_input_mat_1_, 
      saturation_Lagrange_multiplier_mat_1_, dummy_error_mat_1_, 
      saturation_error_mat_1_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_saturation_Lagrange_multiplier_mat_[i][j] 
          = saturation_Lagrange_multiplier_mat_[i][j] 
          + finite_difference_step_ * saturation_error_mat_1_[i][j];
    }
  }
  computeErrorsForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, current_solution_vec, 
      incremented_state_mat_, incremented_lambda_mat_, 
      incremented_saturation_Lagrange_multiplier_mat_, 
      control_input_and_constraints_error_seq_3_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i] = 
        current_solution_vec[i] + finite_difference_step_ 
        * control_input_and_constraints_update_seq_[i];
  }
  computeStateAndLambdaFromErrors(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, state_error_mat_1_, 
      lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
  computeErrorsDifferenceForSaturation(
      current_solution_vec, dummy_input_mat_, 
      saturation_Lagrange_multiplier_mat_, 
      control_input_and_constraints_update_seq_, saturation_update_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_saturation_Lagrange_multiplier_mat_[i][j] = 
          saturation_Lagrange_multiplier_mat_[i][j] 
          - finite_difference_step_ * saturation_update_mat_[i][j];
    }
  }
  computeErrorsForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, incremented_state_mat_, 
      incremented_lambda_mat_, 
      incremented_saturation_Lagrange_multiplier_mat_, 
      control_input_and_constraints_error_seq_2_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    b_vec[i] = 
        (1/finite_difference_step_-zeta_) 
        * control_input_and_constraints_error_seq_[i] 
        - control_input_and_constraints_error_seq_3_[i] 
          / finite_difference_step_ 
        - (control_input_and_constraints_error_seq_2_[i]
           -control_input_and_constraints_error_seq_1_[i])
        / finite_difference_step_;
  }
}

void MSCGMRESWithSaturation::axFunc(const double time, const double* state_vec,
                                    const double* current_solution_vec, 
                                    const double* direction_vec, 
                                    double* ax_vec) {
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i] = 
        current_solution_vec[i] + finite_difference_step_ * direction_vec[i];
  }
  computeStateAndLambdaFromErrors(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, state_error_mat_1_, 
      lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
  computeErrorsDifferenceForSaturation(
      current_solution_vec, dummy_input_mat_, 
      saturation_Lagrange_multiplier_mat_, direction_vec, 
      saturation_update_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_saturation_Lagrange_multiplier_mat_[i][j] = 
          saturation_Lagrange_multiplier_mat_[i][j] 
          - finite_difference_step_ * saturation_update_mat_[i][j];
    }
  }
  computeErrorsForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, incremented_state_mat_, 
      incremented_lambda_mat_, 
      incremented_saturation_Lagrange_multiplier_mat_, 
      control_input_and_constraints_error_seq_2_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    ax_vec[i] = 
        (control_input_and_constraints_error_seq_2_[i]
            -control_input_and_constraints_error_seq_1_[i]) 
        / finite_difference_step_;
  }
}