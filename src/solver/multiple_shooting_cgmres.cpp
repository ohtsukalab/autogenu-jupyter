#include "multiple_shooting_cgmres.hpp"

MultipleShootingCGMRES::MultipleShootingCGMRES(
    const double T_f, const double alpha, const int N, 
    const double finite_difference_increment, const double zeta, const int kmax)
  : continuation_problem_(T_f, alpha, N, finite_difference_increment, zeta),
    mfgmres_(continuation_problem_.dim_condensed_problem(), kmax),
    solution_initializer_(finite_difference_increment, kmax),
    dim_state_(continuation_problem_.dim_state()),
    dim_control_input_(continuation_problem_.dim_control_input()),
    dim_constraints_(continuation_problem_.dim_constraints()),
    N_(N),
    control_input_and_constraints_seq_(
        linearalgebra::NewVector(N*(dim_control_input_+dim_constraints_))),
    control_input_and_constraints_update_seq_(
        linearalgebra::NewVector(N*(dim_control_input_+dim_constraints_))),
    initial_control_input_and_constraints_vec_(
        linearalgebra::NewVector(dim_control_input_+dim_constraints_)),
    initial_lambda_vec_(linearalgebra::NewVector(dim_state_)),
    state_mat_(linearalgebra::NewMatrix(N, dim_state_)),
    lambda_mat_(linearalgebra::NewMatrix(N, dim_state_)) {
}

MultipleShootingCGMRES::~MultipleShootingCGMRES() {
  linearalgebra::DeleteVector(control_input_and_constraints_seq_);
  linearalgebra::DeleteVector(control_input_and_constraints_update_seq_);
  linearalgebra::DeleteVector(initial_control_input_and_constraints_vec_);
  linearalgebra::DeleteVector(initial_lambda_vec_);
  linearalgebra::DeleteMatrix(state_mat_);
  linearalgebra::DeleteMatrix(lambda_mat_);
}

void MultipleShootingCGMRES::controlUpdate(const double time, 
                                           const double* state_vec,
                                           const double sampling_period, 
                                           double* control_input_vec) {
  mfgmres_.solveLinearProblem(continuation_problem_, time, state_vec,
                              control_input_and_constraints_seq_,
                              state_mat_, lambda_mat_, 
                              control_input_and_constraints_update_seq_);
  continuation_problem_.integrateSolution(
      control_input_and_constraints_seq_, state_mat_, lambda_mat_, 
      control_input_and_constraints_update_seq_, sampling_period);
  getControlInput(control_input_vec);
}

void MultipleShootingCGMRES::getControlInput(double* control_input_vec) const {
  for (int i=0; i<dim_control_input_; ++i) {
    control_input_vec[i] = control_input_and_constraints_seq_[i];
  }
}

void MultipleShootingCGMRES::setParametersForInitialization(
    const double* initial_guess_solution, 
    const double newton_residual_tolerance,
    const int max_newton_iteration) {
  solution_initializer_.setInitialGuessSolution(initial_guess_solution);
  solution_initializer_.setCriterionsOfNewtonTermination(
      newton_residual_tolerance, max_newton_iteration);
}

void MultipleShootingCGMRES::initializeSolution(
    const double initial_time, const double* initial_state_vec) {
  solution_initializer_.computeInitialSolution(
      initial_time, initial_state_vec, 
      initial_control_input_and_constraints_vec_);
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
  solution_initializer_.getInitialLambda(initial_time, initial_state_vec, 
                                       initial_lambda_vec_);
  for (int i=0; i<N_; ++i) {
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat_[i][j] = initial_lambda_vec_[j];
    }
  }
  continuation_problem_.resetHorizonLength(initial_time);
}

double MultipleShootingCGMRES::getErrorNorm(const double time, 
                                            const double* state_vec) {
  return continuation_problem_.computeErrorNorm(
      time, state_vec, control_input_and_constraints_seq_,state_mat_, 
      lambda_mat_);
}