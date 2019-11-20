#include "ms_cgmres_with_input_saturation.hpp"

MSCGMRESWithInputSaturation::MSCGMRESWithInputSaturation(
    const InputSaturationSet& input_saturation_set, const double T_f, 
    const double alpha, const int N, const double finite_difference_increment,
    const double zeta, const int kmax)
  : continuation_problem_(input_saturation_set, T_f, alpha, N, 
                        finite_difference_increment, zeta),
    mfgmres_(continuation_problem_.dim_condensed_problem(), kmax),
    solution_initializer_(input_saturation_set, finite_difference_increment, 
                          kmax),
    dim_state_(continuation_problem_.dim_state()),
    dim_control_input_(continuation_problem_.dim_control_input()),
    dim_constraints_(continuation_problem_.dim_constraints()),
    dim_saturation_(continuation_problem_.dim_saturation()),
    N_(N),
    control_input_and_constraints_seq_(
        linearalgebra::NewVector(N*(dim_control_input_+dim_constraints_))),
    control_input_and_constraints_update_seq_(
        linearalgebra::NewVector(N*(dim_control_input_+dim_constraints_))),
    initial_control_input_and_constraints_vec_(
        linearalgebra::NewVector(dim_control_input_+dim_constraints_)),
    initial_lambda_vec_(linearalgebra::NewVector(dim_state_)),
    initial_dummy_input_vec_(linearalgebra::NewVector(dim_saturation_)),
    initial_input_saturation_vec_(linearalgebra::NewVector(dim_saturation_)),
    state_mat_(linearalgebra::NewMatrix(N, dim_state_)),
    lambda_mat_(linearalgebra::NewMatrix(N, dim_state_)),
    dummy_input_mat_(linearalgebra::NewMatrix(N, dim_saturation_)),
    input_saturation_multiplier_mat_(
        linearalgebra::NewMatrix(N, dim_saturation_)) {
}

MSCGMRESWithInputSaturation::~MSCGMRESWithInputSaturation() {
  linearalgebra::DeleteVector(control_input_and_constraints_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_update_seq_);
  linearalgebra::DeleteVector(initial_control_input_and_constraints_vec_);
  linearalgebra::DeleteVector(initial_lambda_vec_);
  linearalgebra::DeleteVector(initial_dummy_input_vec_);
  linearalgebra::DeleteVector(initial_input_saturation_vec_);
  linearalgebra::DeleteMatrix(state_mat_);
  linearalgebra::DeleteMatrix(lambda_mat_);
  linearalgebra::DeleteMatrix(dummy_input_mat_);
  linearalgebra::DeleteMatrix(input_saturation_multiplier_mat_);
}

void MSCGMRESWithInputSaturation::controlUpdate(const double time, 
                                           const double* state_vec, 
                                           const double sampling_period, 
                                           double* control_input_vec) {
  mfgmres_.solveLinearProblem(continuation_problem_, time, state_vec, 
                              control_input_and_constraints_seq_, 
                              state_mat_, lambda_mat_, dummy_input_mat_, 
                              input_saturation_multiplier_mat_,
                              control_input_and_constraints_update_seq_);
  continuation_problem_.integrateSolution(
      control_input_and_constraints_seq_, state_mat_, lambda_mat_, 
      dummy_input_mat_, input_saturation_multiplier_mat_,
      control_input_and_constraints_update_seq_, sampling_period);
  getControlInput(control_input_vec);
}

void MSCGMRESWithInputSaturation::getControlInput(
    double* control_input_vec) const {
  for (int i=0; i<dim_control_input_; ++i) {
    control_input_vec[i] = control_input_and_constraints_seq_[i];
  }
}

void MSCGMRESWithInputSaturation::setParametersForInitialization(
    const double* initial_guess_solution, 
    const double newton_residual_tolerance, const int max_newton_iteration) {
  solution_initializer_.setInitialGuessSolution(initial_guess_solution);
  solution_initializer_.setCriterionsOfNewtonTermination(
    newton_residual_tolerance, max_newton_iteration);
}

void MSCGMRESWithInputSaturation::setInitialInputSaturationMultiplier(
    const double initial_input_saturation_multiplier) {
  solution_initializer_.setInitialInputSaturationMultiplier(
      initial_input_saturation_multiplier);
}

void MSCGMRESWithInputSaturation::setInitialInputSaturationMultiplier(
    const double* initial_input_saturation_multiplier) {
  solution_initializer_.setInitialInputSaturationMultiplier(
      initial_input_saturation_multiplier);
}

void MSCGMRESWithInputSaturation::initializeSolution(
    const double initial_time, const double* initial_state_vec) {
  solution_initializer_.computeInitialSolution(
      initial_time, initial_state_vec, 
      initial_control_input_and_constraints_vec_, initial_dummy_input_vec_, 
      initial_input_saturation_vec_);
  solution_initializer_.getInitialLambda(initial_time, initial_state_vec, 
                                         initial_lambda_vec_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_control_input_+dim_constraints_; ++j) {
      control_input_and_constraints_seq_[i*(dim_control_input_+dim_constraints_)+j] 
          = initial_control_input_and_constraints_vec_[j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      state_mat_[i][j] = initial_state_vec[j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat_[i][j] = initial_lambda_vec_[j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      dummy_input_mat_[i][j] = initial_dummy_input_vec_[j];
    }
  }
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_saturation_; ++j) {
      input_saturation_multiplier_mat_[i][j] = initial_input_saturation_vec_[j];
    }
  }
  continuation_problem_.resetHorizonLength(initial_time);
}

double MSCGMRESWithInputSaturation::getErrorNorm(const double time, 
                                            const double* state_vec) {
  return continuation_problem_.computeErrorNorm(
      time, state_vec, control_input_and_constraints_seq_, state_mat_,
      lambda_mat_, dummy_input_mat_, input_saturation_multiplier_mat_);
}