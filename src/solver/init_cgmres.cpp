#include "init_cgmres.hpp"


InitCGMRES::InitCGMRES() 
  : MatrixFreeGMRES(), 
    model_(), 
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()), 
    max_iteration_(0),
    finite_difference_step_(0), 
    residual_tolerance_(0),
    initial_guess_solution_(linearalgebra::NewVector(dim_solution_)),
    solution_update_vec_(linearalgebra::NewVector(dim_solution_)), 
    incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)), 
    lambda_vec_(linearalgebra::NewVector(model_.dimState())), 
    errors_in_optimality_(linearalgebra::NewVector(dim_solution_)), 
    errors_in_optimality_1_(linearalgebra::NewVector(dim_solution_)), 
    errors_in_optimality_2_(linearalgebra::NewVector(dim_solution_)) {
  // Set dimensions and parameters used in GMRES.
  setGMRESParameters(dim_solution_, dim_solution_);
}

InitCGMRES::~InitCGMRES() {
  linearalgebra::DeleteVector(initial_guess_solution_);
  linearalgebra::DeleteVector(solution_update_vec_);
  linearalgebra::DeleteVector(incremented_solution_vec_);
  linearalgebra::DeleteVector(lambda_vec_);
  linearalgebra::DeleteVector(errors_in_optimality_);
  linearalgebra::DeleteVector(errors_in_optimality_1_);
  linearalgebra::DeleteVector(errors_in_optimality_2_);
}

void InitCGMRES::setInitParameters(const double* initial_guess_solution, 
                                   const double residual_tolerance, 
                                   const int max_iteration, 
                                   const double finite_difference_step, 
                                   const int kmax) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_guess_solution_[i] = initial_guess_solution[i];
  }
  residual_tolerance_ = residual_tolerance;
  max_iteration_ = max_iteration;
  finite_difference_step_ = finite_difference_step;
  if(kmax < dim_solution_) {
    setGMRESParameters(dim_solution_, kmax);
  }
  else {
    setGMRESParameters(dim_solution_, dim_solution_);
  }
}

void InitCGMRES::solveOCPForInit(const double initial_time, 
                                 const double* initial_state_vec, 
                                 double* initial_solution_vec, 
                                 double* initial_errors_in_optimality) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_solution_vec[i] = initial_guess_solution_[i];
  }
  computeErrorsInOptimality(initial_time, initial_state_vec, 
                            initial_solution_vec, 
                            initial_errors_in_optimality);
  int iteration_counter = 0;
  while (linearalgebra::SquaredNorm(dim_solution_, 
                                    initial_errors_in_optimality)
             > residual_tolerance_*residual_tolerance_ 
         && iteration_counter < max_iteration_) {
    solveGMRES(initial_time, initial_state_vec, initial_solution_vec, 
               solution_update_vec_);
    for (int i=0; i<dim_solution_; ++i) {
        initial_solution_vec[i] += solution_update_vec_[i];
    }
    computeErrorsInOptimality(initial_time, initial_state_vec, 
                              initial_solution_vec, 
                              initial_errors_in_optimality);
    ++iteration_counter;
  }
}

inline void InitCGMRES::computeErrorsInOptimality(
    const double time, const double* state_vec, 
    const double* solution_vec, double* errors_in_optimality) {
  model_.phixFunc(time, state_vec, lambda_vec_);
  model_.huFunc(time, state_vec, solution_vec, lambda_vec_, 
                errors_in_optimality);
}

void InitCGMRES::bFunc(const double time, const double* state_vec, 
                       const double* solution_vec, double* b_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = solution_vec[i] 
        + finite_difference_step_ * solution_update_vec_[i];
  }
  computeErrorsInOptimality(time, state_vec, solution_vec, 
                            errors_in_optimality_);
  computeErrorsInOptimality(time, state_vec, incremented_solution_vec_, 
                            errors_in_optimality_1_);
  for (int i=0; i<dim_solution_; ++i) {
    b_vec[i] = - errors_in_optimality_[i] 
        - (errors_in_optimality_1_[i]-errors_in_optimality_[i]) 
        / finite_difference_step_;
  }
}

void InitCGMRES::axFunc(const double time, const double* state_vec, 
                        const double* solution_vec, 
                        const double* direction_vec, double* ax_vec) {
  for (int i=0; i<dim_solution_; ++i) { 
    incremented_solution_vec_[i] = solution_vec[i] 
        + finite_difference_step_ * direction_vec[i];
  }
  computeErrorsInOptimality(time, state_vec, incremented_solution_vec_, 
                            errors_in_optimality_1_);
  for (int i=0; i<dim_solution_; ++i) {
    ax_vec[i] = (errors_in_optimality_1_[i]-errors_in_optimality_[i]) 
        / finite_difference_step_;
  }
}