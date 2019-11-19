#include "ms_continuation_with_input_saturation.hpp"

MSContinuationWithInputSaturation::MSContinuationWithInputSaturation(
    const InputSaturationSet& input_saturation_set, 
    const double T_f, const double alpha, const int N,
    const double finite_difference_increment, const double zeta)
  : ocp_(input_saturation_set, T_f, alpha, N),
    dim_state_(ocp_.dim_state()),
    dim_control_input_(ocp_.dim_control_input()),
    dim_constraints_(ocp_.dim_constraints()),
    dim_control_input_and_constraints_(
        ocp_.dim_control_input()+ocp_.dim_constraints()), 
    dim_saturation_(ocp_.dim_saturation()),
    dim_control_input_and_constraints_seq_(
        N*(ocp_.dim_control_input()+ocp_.dim_constraints())), 
    N_(N),
    finite_difference_increment_(finite_difference_increment),
    zeta_(zeta),
    incremented_time_(0),
    incremented_state_vec_(linearalgebra::NewVector(ocp_.dim_state())),
    incremented_control_input_and_constraints_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_1_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_2_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_3_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    incremented_state_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    incremented_lambda_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    state_residual_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    state_residual_mat_1_(linearalgebra::NewMatrix(N_, dim_state_)),
    lambda_residual_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    lambda_residual_mat_1_(linearalgebra::NewMatrix(N_, dim_state_)),
    incremented_dummy_input_mat_(linearalgebra::NewMatrix(N_, dim_saturation_)),
    incremented_input_sautration_multiplier_mat_(
        linearalgebra::NewMatrix(N_, dim_saturation_)),
    dummy_input_residual_mat_(linearalgebra::NewMatrix(N_, dim_saturation_)), 
    dummy_input_residual_mat_1_(linearalgebra::NewMatrix(N_, dim_saturation_)), 
    input_saturation_residual_mat_(
        linearalgebra::NewMatrix(N_, dim_saturation_)), 
    input_saturation_residual_mat_1_(
        linearalgebra::NewMatrix(N_, dim_saturation_)), 
    dummy_input_difference_mat_(linearalgebra::NewMatrix(N_, dim_saturation_)), 
    input_saturation_multiplier_difference_mat_(
        linearalgebra::NewMatrix(N_, dim_saturation_)) {
}

MSContinuationWithInputSaturation::MSContinuationWithInputSaturation(
    const InputSaturationSet& input_saturation_set, 
    const double T_f, const double alpha, const int N,
    const double initial_time,
    const double finite_difference_increment, const double zeta)
  : ocp_(input_saturation_set, T_f, alpha, N, initial_time),
    dim_state_(ocp_.dim_state()),
    dim_control_input_(ocp_.dim_control_input()),
    dim_constraints_(ocp_.dim_constraints()),
    dim_control_input_and_constraints_(
        ocp_.dim_control_input()+ocp_.dim_constraints()), 
    dim_saturation_(ocp_.dim_saturation()),
    dim_control_input_and_constraints_seq_(
        N*(ocp_.dim_control_input()+ocp_.dim_constraints())), 
    N_(N),
    finite_difference_increment_(finite_difference_increment),
    zeta_(zeta),
    incremented_time_(0),
    incremented_state_vec_(linearalgebra::NewVector(ocp_.dim_state())),
    incremented_control_input_and_constraints_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_1_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_2_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    control_input_and_constraints_residual_seq_3_(
        linearalgebra::NewVector(dim_control_input_and_constraints_seq_)),
    incremented_state_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    incremented_lambda_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    state_residual_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    state_residual_mat_1_(linearalgebra::NewMatrix(N_, dim_state_)),
    lambda_residual_mat_(linearalgebra::NewMatrix(N_, dim_state_)),
    lambda_residual_mat_1_(linearalgebra::NewMatrix(N_, dim_state_)),
    incremented_dummy_input_mat_(linearalgebra::NewMatrix(N_, dim_saturation_)),
    incremented_input_sautration_multiplier_mat_(
        linearalgebra::NewMatrix(N_, dim_saturation_)),
    dummy_input_residual_mat_(linearalgebra::NewMatrix(N_, dim_saturation_)), 
    dummy_input_residual_mat_1_(linearalgebra::NewMatrix(N_, dim_saturation_)), 
    input_saturation_residual_mat_(
        linearalgebra::NewMatrix(N_, dim_saturation_)), 
    input_saturation_residual_mat_1_(
        linearalgebra::NewMatrix(N_, dim_saturation_)), 
    dummy_input_difference_mat_(linearalgebra::NewMatrix(N_, dim_saturation_)), 
    input_saturation_multiplier_difference_mat_(
        linearalgebra::NewMatrix(N_, dim_saturation_)) {
}

MSContinuationWithInputSaturation::~MSContinuationWithInputSaturation() {
  linearalgebra::DeleteVector(incremented_state_vec_);
  linearalgebra::DeleteVector(incremented_control_input_and_constraints_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_residual_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_residual_seq_1_);
  linearalgebra::DeleteVector(control_input_and_constraints_residual_seq_2_);
  linearalgebra::DeleteVector(control_input_and_constraints_residual_seq_3_);
  linearalgebra::DeleteMatrix(incremented_state_mat_);
  linearalgebra::DeleteMatrix(incremented_lambda_mat_);
  linearalgebra::DeleteMatrix(state_residual_mat_);
  linearalgebra::DeleteMatrix(state_residual_mat_1_);
  linearalgebra::DeleteMatrix(lambda_residual_mat_);
  linearalgebra::DeleteMatrix(lambda_residual_mat_1_);
  linearalgebra::DeleteMatrix(incremented_dummy_input_mat_);
  linearalgebra::DeleteMatrix(incremented_input_sautration_multiplier_mat_);
  linearalgebra::DeleteMatrix(dummy_input_residual_mat_);
  linearalgebra::DeleteMatrix(dummy_input_residual_mat_1_);
  linearalgebra::DeleteMatrix(input_saturation_residual_mat_);
  linearalgebra::DeleteMatrix(input_saturation_residual_mat_1_);
  linearalgebra::DeleteMatrix(dummy_input_difference_mat_);
  linearalgebra::DeleteMatrix(input_saturation_multiplier_difference_mat_);
}

void MSContinuationWithInputSaturation::integrateSolution(
    double* control_input_and_constraints_seq, 
    double** state_mat, double** lambda_mat,
    double** dummy_input_mat, 
    double** input_saturation_multiplier_mat,
    const double* control_input_and_constraints_update_seq, 
    const double integration_length) {
  // Update state_mat_ and lamdba_mat_ by the difference approximation.
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i] 
        = control_input_and_constraints_seq[i] 
            + finite_difference_increment_
                * control_input_and_constraints_update_seq[i];
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      state_residual_mat_1_[i][j] = 
          (1-finite_difference_increment_*zeta_) * state_residual_mat_[i][j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) { 
      lambda_residual_mat_1_[i][j] = 
          (1-finite_difference_increment_*zeta_) * lambda_residual_mat_[i][j];
    }
  }
  ocp_.computeStateAndLambdaFromOptimalityResidual(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, state_residual_mat_1_, 
      lambda_residual_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
  // state_mat_ += 
  //     (sampling_period/finite_difference_step_) 
  //     * (incremented_state_mat_-state_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      state_mat[i][j] += (integration_length/finite_difference_increment_) 
                          * (incremented_state_mat_[i][j]-state_mat[i][j]);
    }
  }
  // lambda_mat_ += 
  //     (sampling_period/finite_difference_step_) 
  //     * (incremented_lambda_mat_-lambda_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat[i][j] += (integration_length/finite_difference_increment_) 
                           * (incremented_lambda_mat_[i][j]-lambda_mat[i][j]);
    }
  }
  ocp_.computeResidualDifferenceForDummyInput(
      control_input_and_constraints_seq, dummy_input_mat, 
      control_input_and_constraints_update_seq, dummy_input_difference_mat_);
  ocp_.computeResidualDifferenceForInputSaturation(
      control_input_and_constraints_seq, dummy_input_mat, 
      input_saturation_multiplier_mat, control_input_and_constraints_update_seq,
      input_saturation_multiplier_difference_mat_);
  for (int i=0; i<N_; ++i) { 
    for (int j=0; j<dim_saturation_; ++j) {
      dummy_input_mat[i][j] 
          += integration_length 
              * (dummy_input_residual_mat_1_[i][j]
                    -dummy_input_difference_mat_[i][j]);
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      input_saturation_multiplier_mat[i][j] 
          += integration_length 
              * (input_saturation_residual_mat_1_[i][j] 
                    -input_saturation_multiplier_difference_mat_[i][j]);
    }
  }
  // Update control_input_and_constraints_seq_
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    control_input_and_constraints_seq[i] 
        += integration_length 
            * control_input_and_constraints_update_seq[i];
  }
}

double MSContinuationWithInputSaturation::computeErrorNorm(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq,
    double const* const* state_mat, 
    double const* const* lambda_mat,
    double const* const* dummy_input_mat, 
    double const* const* input_saturation_multiplier_mat) {
  ocp_.computeOptimalityResidualForControlInputAndConstraints(
      time, state_vec, control_input_and_constraints_seq, state_mat, lambda_mat, 
      input_saturation_multiplier_mat, 
      control_input_and_constraints_residual_seq_);
  ocp_.computeOptimalityResidualForStateAndLambda(
      time, state_vec, control_input_and_constraints_seq, state_mat, lambda_mat,
      state_residual_mat_, lambda_residual_mat_);
  ocp_.computeResidualForDummyInputAndInputSaturation(
      control_input_and_constraints_seq, dummy_input_mat, 
      input_saturation_multiplier_mat, dummy_input_residual_mat_, 
      input_saturation_residual_mat_);
  double squared_error_norm 
      = linearalgebra::SquaredNorm(dim_control_input_and_constraints_seq_, 
                                   control_input_and_constraints_residual_seq_);
  for (int i=0; i<N_; ++i) {
    squared_error_norm += linearalgebra::SquaredNorm(dim_state_, 
                                                     state_residual_mat_[i]);
  }
  for (int i=0; i<N_; ++i) {
    squared_error_norm += linearalgebra::SquaredNorm(dim_state_, 
                                                     lambda_residual_mat_[i]);
  }
  for (int i=0; i<N_; ++i) {
    squared_error_norm 
        += linearalgebra::SquaredNorm(dim_saturation_, 
                                      dummy_input_residual_mat_[i]);
  }
  for (int i=0; i<N_; ++i) {
    squared_error_norm 
        += linearalgebra::SquaredNorm(dim_saturation_, 
                                      input_saturation_residual_mat_[i]);
  }
  return std::sqrt(squared_error_norm);
}

void MSContinuationWithInputSaturation::resetHorizonLength(
    const double T_f, const double alpha, const double initial_time) {
  ocp_.resetHorizonLength(T_f, alpha, initial_time);
}

void MSContinuationWithInputSaturation::resetHorizonLength(
    const double initial_time) {
  ocp_.resetHorizonLength(initial_time);
}

void MSContinuationWithInputSaturation::bFunc(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat,
    double const* const* dummy_input_mat, 
    double const* const* input_saturation_multiplier_mat,
    const double* current_control_input_and_constraints_update_seq, 
    double* b_vec) {
  incremented_time_ = time + finite_difference_increment_;
  ocp_.predictStateFromSolution(time, state_vec, 
                                control_input_and_constraints_seq,
                                finite_difference_increment_, 
                                incremented_state_vec_);
  ocp_.computeOptimalityResidualForControlInputAndConstraints(
      time, state_vec, control_input_and_constraints_seq, state_mat, lambda_mat, 
      input_saturation_multiplier_mat,
      control_input_and_constraints_residual_seq_);
  ocp_.computeOptimalityResidualForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, 
      control_input_and_constraints_seq, state_mat, lambda_mat, 
      input_saturation_multiplier_mat,
      control_input_and_constraints_residual_seq_1_);
  ocp_.computeOptimalityResidualForStateAndLambda(
      time, state_vec, control_input_and_constraints_seq, state_mat, lambda_mat, 
      state_residual_mat_, lambda_residual_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      state_residual_mat_1_[i][j] = 
          (1-finite_difference_increment_*zeta_) * state_residual_mat_[i][j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      lambda_residual_mat_1_[i][j] = 
          (1-finite_difference_increment_*zeta_) * lambda_residual_mat_[i][j];
    }
  }
  ocp_.computeStateAndLambdaFromOptimalityResidual(
      incremented_time_, incremented_state_vec_, 
      control_input_and_constraints_seq, 
      state_residual_mat_1_, lambda_residual_mat_1_, 
      incremented_state_mat_, incremented_lambda_mat_);
  ocp_.computeOptimalityResidualForStateAndLambda(
      incremented_time_, incremented_state_vec_, 
      control_input_and_constraints_seq, 
      state_mat, lambda_mat, state_residual_mat_1_, lambda_residual_mat_1_);
  ocp_.computeResidualForDummyInputAndInputSaturation(
      control_input_and_constraints_seq, dummy_input_mat, 
      input_saturation_multiplier_mat, dummy_input_residual_mat_, 
      input_saturation_residual_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_dummy_input_mat_[i][j] 
          = - zeta_ * dummy_input_residual_mat_[i][j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_input_sautration_multiplier_mat_[i][j] 
          = - zeta_ * input_saturation_residual_mat_[i][j];
    }
  }
  ocp_.multiplyResidualForDummyInputAndInputSaturationInverse(
      control_input_and_constraints_seq, dummy_input_mat, 
      input_saturation_multiplier_mat, incremented_dummy_input_mat_,
      incremented_input_sautration_multiplier_mat_, dummy_input_residual_mat_1_,
      input_saturation_residual_mat_1_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_input_sautration_multiplier_mat_[i][j] 
          = input_saturation_multiplier_mat[i][j] 
          + finite_difference_increment_ 
            * input_saturation_residual_mat_1_[i][j];
    }
  }
  ocp_.computeOptimalityResidualForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, 
      control_input_and_constraints_seq, incremented_state_mat_,
      incremented_lambda_mat_, incremented_input_sautration_multiplier_mat_,
      control_input_and_constraints_residual_seq_3_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i]
        = control_input_and_constraints_seq[i] 
            + finite_difference_increment_ 
               * current_control_input_and_constraints_update_seq[i];
  }
  ocp_.computeStateAndLambdaFromOptimalityResidual(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, state_residual_mat_1_,
      lambda_residual_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
  ocp_.computeResidualDifferenceForInputSaturation(
      control_input_and_constraints_seq, dummy_input_mat,
      input_saturation_multiplier_mat, 
      current_control_input_and_constraints_update_seq,
      input_saturation_multiplier_difference_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_input_sautration_multiplier_mat_[i][j] 
          = input_saturation_multiplier_mat[i][j] 
              - finite_difference_increment_ 
                  * input_saturation_multiplier_difference_mat_[i][j];
    }
  }
  ocp_.computeOptimalityResidualForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, incremented_state_mat_,
      incremented_lambda_mat_, incremented_input_sautration_multiplier_mat_,
      control_input_and_constraints_residual_seq_2_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    b_vec[i] = 
        (1/finite_difference_increment_-zeta_) 
        * control_input_and_constraints_residual_seq_[i] 
        - control_input_and_constraints_residual_seq_3_[i] 
          / finite_difference_increment_ 
        - (control_input_and_constraints_residual_seq_2_[i]
           -control_input_and_constraints_residual_seq_1_[i])
        / finite_difference_increment_;
  }
}

void MSContinuationWithInputSaturation::AxFunc(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat,
    double const* const* dummy_input_mat, 
    double const* const* input_saturation_multiplier_mat,
    const double* direction_vec, double* ax_vec) {
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    incremented_control_input_and_constraints_seq_[i]
        = control_input_and_constraints_seq[i] 
            + finite_difference_increment_ * direction_vec[i];
  }
  ocp_.computeStateAndLambdaFromOptimalityResidual(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, state_residual_mat_1_,
      lambda_residual_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
  ocp_.computeResidualDifferenceForInputSaturation(
      control_input_and_constraints_seq, dummy_input_mat,
      input_saturation_multiplier_mat, direction_vec,
      input_saturation_multiplier_difference_mat_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      incremented_input_sautration_multiplier_mat_[i][j] 
          = input_saturation_multiplier_mat[i][j] 
              - finite_difference_increment_ 
                  * input_saturation_multiplier_difference_mat_[i][j];
    }
  }
  ocp_.computeOptimalityResidualForControlInputAndConstraints(
      incremented_time_, incremented_state_vec_, 
      incremented_control_input_and_constraints_seq_, incremented_state_mat_,
      incremented_lambda_mat_, incremented_input_sautration_multiplier_mat_,
      control_input_and_constraints_residual_seq_2_);
  for (int i=0; i<dim_control_input_and_constraints_seq_; ++i) {
    ax_vec[i] 
        = (control_input_and_constraints_residual_seq_2_[i]
            -control_input_and_constraints_residual_seq_1_[i]) 
          / finite_difference_increment_;
  }
}

int MSContinuationWithInputSaturation::dim_state() const {
  return dim_state_;
}

int MSContinuationWithInputSaturation::dim_control_input() const {
  return dim_control_input_;
}

int MSContinuationWithInputSaturation::dim_constraints() const {
  return dim_constraints_;
}

int MSContinuationWithInputSaturation::dim_saturation() const {
  return ocp_.dim_saturation();
}

int MSContinuationWithInputSaturation::dim_condensed_problem() const {
  return dim_control_input_and_constraints_seq_;
}

int MSContinuationWithInputSaturation::N() const {
  return ocp_.N();
}