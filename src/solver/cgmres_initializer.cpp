#include "cgmres_initializer.hpp"

CGMRESInitializer::CGMRESInitializer(const double finite_difference_increment,
                                     const int kmax,
                                     const double residual_tolerance, 
                                     const int max_newton_iteration)
  : newton_(finite_difference_increment),
    mfgmres_(newton_.dim_solution(), kmax),
    dim_control_input_(newton_.dim_control_input()),
    dim_constraints_(newton_.dim_constraints()),
    dim_solution_(newton_.dim_solution()),
    max_newton_iteration_(max_newton_iteration),
    newton_residual_tolerance_(residual_tolerance),
    initial_guess_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    solution_update_vec_(linearalgebra::NewVector(dim_solution_)) {
}

CGMRESInitializer::CGMRESInitializer(const double finite_difference_increment,
                                     const int kmax) 
  : newton_(finite_difference_increment),
    mfgmres_(newton_.dim_solution(), kmax),
    dim_control_input_(newton_.dim_control_input()),
    dim_constraints_(newton_.dim_constraints()),
    dim_solution_(newton_.dim_solution()),
    max_newton_iteration_(50),
    newton_residual_tolerance_(1e-08),
    initial_guess_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    solution_update_vec_(linearalgebra::NewVector(dim_solution_)) {
}

CGMRESInitializer::~CGMRESInitializer() {
  linearalgebra::DeleteVector(initial_guess_solution_vec_);
  linearalgebra::DeleteVector(solution_update_vec_);
}

void CGMRESInitializer::setCriterionsOfNewtonTermination(
    const double newton_residual_tolerance, 
    const int max_newton_iteration) {
  newton_residual_tolerance_ = newton_residual_tolerance;
  max_newton_iteration_ = max_newton_iteration;
}

void CGMRESInitializer::setInitialGuessSolution(
    const double* initial_guess_solution) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_guess_solution_vec_[i] = initial_guess_solution[i];
  }
}

void CGMRESInitializer::computeInitialSolution(const double initial_time, 
                                               const double* initial_state_vec, 
                                               double* initial_solution_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_solution_vec[i] = initial_guess_solution_vec_[i];
  }
  int num_itr = 0;
  double optimality_error = newton_.errorNorm(initial_time, initial_state_vec, 
                                              initial_solution_vec);
  while (optimality_error > newton_residual_tolerance_ 
         && num_itr < max_newton_iteration_) {
    mfgmres_.solveLinearProblem(newton_, initial_time, initial_state_vec, 
                                initial_solution_vec, solution_update_vec_);
    for (int i=0; i<dim_solution_; ++i) {
      initial_solution_vec[i] += solution_update_vec_[i];
    }
    optimality_error = newton_.errorNorm(initial_time, initial_state_vec, 
                                         initial_solution_vec);
    ++num_itr;
  }
}

void CGMRESInitializer::getInitialLambda(const double initial_time, 
                                         const double* initial_state_vec, 
                                         double* initial_lambda_vec) {
  newton_.getTerminalCostDerivatives(initial_time, initial_state_vec, 
                                     initial_lambda_vec);
}

int CGMRESInitializer::dim_solution() const {
  return dim_solution_;
}