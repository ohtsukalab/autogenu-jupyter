#include "continuation_gmres.hpp"

ContinuationGMRES::ContinuationGMRES(const double T_f, const double alpha,
                                     const int N, 
                                     const double finite_difference_increment,
                                     const double zeta, const int kmax)
  : continuation_problem_(T_f, alpha, N, finite_difference_increment, zeta),
    mfgmres_(continuation_problem_.dim_solution(), kmax),
    solution_initializer_(finite_difference_increment, kmax),
    dim_control_input_(continuation_problem_.dim_control_input()),
    dim_constraints_(continuation_problem_.dim_constraints()),
    solution_vec_(
        linearalgebra::NewVector(continuation_problem_.dim_solution())),
    solution_update_vec_(
        linearalgebra::NewVector(continuation_problem_.dim_solution())), 
    initial_solution_vec_(
        linearalgebra::NewVector(solution_initializer_.dim_solution())) {
}

ContinuationGMRES::~ContinuationGMRES() {
  linearalgebra::DeleteVector(solution_vec_);
  linearalgebra::DeleteVector(solution_update_vec_);
  linearalgebra::DeleteVector(initial_solution_vec_);
}

void ContinuationGMRES::controlUpdate(const double time, const double* state_vec, 
                                      const double sampling_period, 
                                      double* control_input_vec) {
  mfgmres_.solveLinearProblem(continuation_problem_, time, state_vec, 
                              solution_vec_, solution_update_vec_);
  continuation_problem_.integrateSolution(solution_vec_, solution_update_vec_, 
                                          sampling_period);
  for (int i=0; i<dim_control_input_; ++i) {
    control_input_vec[i] = solution_vec_[i];
  }
}

void ContinuationGMRES::setParametersForInitialization(
    const double* initial_guess_solution, 
    const double newton_residual_tolerance,
    const int max_newton_iteration) {
  solution_initializer_.setInitialGuessSolution(initial_guess_solution);
  solution_initializer_.setCriterionsOfNewtonTermination(
      newton_residual_tolerance, max_newton_iteration);
}

void ContinuationGMRES::initializeSolution(const double initial_time,  
                                           const double* initial_state_vec) {
  solution_initializer_.computeInitialSolution(initial_time, initial_state_vec, 
                                             initial_solution_vec_);
  for (int i=0; i<continuation_problem_.N(); ++i) {
    for (int j=0; j<solution_initializer_.dim_solution(); ++j) {
      solution_vec_[i*(dim_control_input_+dim_constraints_)+j] 
          = initial_solution_vec_[j];
    }
  }
  continuation_problem_.resetHorizonLength(initial_time);
}

void ContinuationGMRES::getControlInput(double* control_input_vec) const {
  for (int i=0; i<dim_control_input_; ++i) {
    control_input_vec[i] = solution_vec_[i];
  }
}

double ContinuationGMRES::getErrorNorm(const double time, 
                                       const double* state_vec) {
  return continuation_problem_.computeErrorNorm(time, state_vec, 
                                                  solution_vec_);
}