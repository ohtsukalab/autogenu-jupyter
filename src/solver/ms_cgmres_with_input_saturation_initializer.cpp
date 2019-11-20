#include "ms_cgmres_with_input_saturation_initializer.hpp"

MSCGMRESWithInputSaturationInitializer::
MSCGMRESWithInputSaturationInitializer(
    const InputSaturationSet& input_saturation_set,
    const double finite_difference_increment, const int kmax, 
    const double newton_residual_tolerance, const int max_newton_iteration)
  : newton_(finite_difference_increment, input_saturation_set),
    mfgmres_(newton_.dim_solution(), kmax),
    input_saturation_set_(input_saturation_set),
    dim_control_input_(newton_.dim_control_input()),
    dim_constraints_(newton_.dim_constraints()),
    dim_input_saturation_(input_saturation_set.dimSaturation()),
    dim_solution_(newton_.dim_solution()),
    max_newton_iteration_(max_newton_iteration),
    newton_residual_tolerance_(newton_residual_tolerance),
    initial_guess_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    initial_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    solution_update_vec_(linearalgebra::NewVector(dim_solution_)) {
}

MSCGMRESWithInputSaturationInitializer::
MSCGMRESWithInputSaturationInitializer(
    const InputSaturationSet& input_saturation_set,
    const double finite_difference_increment, const int kmax)
  : newton_(finite_difference_increment, input_saturation_set),
    mfgmres_(newton_.dim_solution(), kmax),
    input_saturation_set_(input_saturation_set),
    dim_control_input_(newton_.dim_control_input()),
    dim_constraints_(newton_.dim_constraints()),
    dim_input_saturation_(input_saturation_set.dimSaturation()),
    dim_solution_(newton_.dim_solution()),
    max_newton_iteration_(50),
    newton_residual_tolerance_(1e-08),
    initial_guess_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    initial_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    solution_update_vec_(linearalgebra::NewVector(dim_solution_)) {
}

MSCGMRESWithInputSaturationInitializer::
~MSCGMRESWithInputSaturationInitializer() {
  linearalgebra::DeleteVector(initial_guess_solution_vec_);
  linearalgebra::DeleteVector(initial_solution_vec_);
  linearalgebra::DeleteVector(solution_update_vec_);
}

void MSCGMRESWithInputSaturationInitializer::setCriterionsOfNewtonTermination(
    const double newton_residual_tolerance, const int max_newton_iteration) {
  newton_residual_tolerance_ = newton_residual_tolerance;
  max_newton_iteration_ = max_newton_iteration;
}

void MSCGMRESWithInputSaturationInitializer::setInitialGuessSolution(
    const double* initial_guess_control_input_and_constraints) {
  for (int i=0; i<dim_control_input_+dim_constraints_; ++i) {
    initial_guess_solution_vec_[i] 
        = initial_guess_control_input_and_constraints[i];
  }
  for (int i=0; i<dim_input_saturation_; ++i) {
    initial_guess_solution_vec_[dim_control_input_+dim_constraints_+i]
        = computeDummyInput(
            initial_guess_control_input_and_constraints[
                  input_saturation_set_.index(i)], 
            input_saturation_set_.min(i), input_saturation_set_.max(i));
  }
  for (int i=0; i<dim_input_saturation_; ++i) {
    initial_guess_solution_vec_[dim_control_input_+dim_constraints_
                                +dim_input_saturation_+i]
        = 0.001;
  }
}

void MSCGMRESWithInputSaturationInitializer::
setInitialInputSaturationMultiplier(
    const double initial_input_saturation_multiplier) {
  for (int i=0; i<input_saturation_set_.dimSaturation(); ++i) {
    initial_guess_solution_vec_[dim_control_input_+dim_constraints_
                                +dim_input_saturation_+i]
        = initial_input_saturation_multiplier;
  }
} 

void MSCGMRESWithInputSaturationInitializer::
setInitialInputSaturationMultiplier(
    const double* initial_input_saturation_multiplier) {
  for (int i=0; i<input_saturation_set_.dimSaturation(); ++i) {
    initial_guess_solution_vec_[dim_control_input_+dim_constraints_
                                +dim_input_saturation_+i]
        = initial_input_saturation_multiplier[i];
  }
} 

void MSCGMRESWithInputSaturationInitializer::computeInitialSolution(
    const double initial_time, const double* initial_state_vec, 
    double* initial_control_input_and_constraints_vec, 
    double* initial_dummy_input_vec, double* initial_input_saturation_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_solution_vec_[i] = initial_guess_solution_vec_[i];
  }
  int num_itr = 0;
  double optimality_error = newton_.errorNorm(initial_time, initial_state_vec, 
                                              initial_solution_vec_);
  while (optimality_error > newton_residual_tolerance_ 
         && num_itr < max_newton_iteration_) {
    mfgmres_.solveLinearProblem(newton_, initial_time, initial_state_vec, 
                                initial_solution_vec_, solution_update_vec_);
    for (int i=0; i<dim_solution_; ++i) {
      initial_solution_vec_[i] += solution_update_vec_[i];
    }
    optimality_error = newton_.errorNorm(initial_time, initial_state_vec, 
                                         initial_solution_vec_);
    ++num_itr;
  }
  for (int i=0; i<dim_control_input_+dim_constraints_; ++i) {
    initial_control_input_and_constraints_vec[i] = initial_solution_vec_[i];
  }
  for (int i=0; i<dim_input_saturation_; ++i) {
    initial_dummy_input_vec[i] 
        = initial_solution_vec_[dim_control_input_+dim_constraints_+i];
  }
  for (int i=0; i<dim_input_saturation_; ++i) {
    initial_input_saturation_vec[i] 
        = initial_solution_vec_[dim_control_input_+dim_constraints_+dim_input_saturation_+i];
  }
}

void MSCGMRESWithInputSaturationInitializer::getInitialLambda(
    const double initial_time, const double* initial_state_vec, 
    double* initial_lambda_vec) {
  newton_.getTerminalCostDerivatives(initial_time, initial_state_vec, 
                                     initial_lambda_vec);
}

int MSCGMRESWithInputSaturationInitializer::dim_solution() const {
  return dim_solution_;
}

double MSCGMRESWithInputSaturationInitializer::computeDummyInput(
    const double input, const double min_input, const double max_input) const {
  if (min_input < input && input < max_input) {
    double max_plus_min = max_input + min_input;
    double max_minus_min = max_input - min_input;
    return std::sqrt(
        (max_minus_min*max_minus_min)/4
            -(input-max_plus_min/2)*(input-max_plus_min/2));
  }
  else {
    return (min_input+max_input) / 2;
  }
}