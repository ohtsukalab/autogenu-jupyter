#include "init_mscgmres_with_saturation.hpp"


InitMSCGMRESWithSaturation::InitMSCGMRESWithSaturation(
    const ControlInputSaturationSequence saturation_seq) 
  : MatrixFreeGMRES(), 
    model_(), 
    saturation_seq_(saturation_seq), 
    max_iteration_(0),
    finite_difference_step_(0),
    dim_control_input_and_constraints_(
        model_.dimControlInput()+model_.dimConstraints()), 
    dim_saturation_(saturation_seq.dimSaturation()), 
    dim_solution_(
        model_.dimControlInput()+model_.dimConstraints()
        +2*saturation_seq.dimSaturation()), 
    initial_guess_solution_(linearalgebra::NewVector(dim_solution_)),
    initial_guess_Lagrange_multiplier_(
        linearalgebra::NewVector(dim_saturation_)),
    incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)), 
    solution_update_vec_(linearalgebra::NewVector(dim_solution_)), 
    lambda_vec_(linearalgebra::NewVector(model_.dimState())), 
    errors_in_optimality_(linearalgebra::NewVector(dim_solution_)), 
    errors_in_optimality_1_(linearalgebra::NewVector(dim_solution_)), 
    errors_in_optimality_2_(linearalgebra::NewVector(dim_solution_)) {
  // Set parameters in GMRES.
  setGMRESParameters(dim_solution_, dim_solution_);
}

InitMSCGMRESWithSaturation::~InitMSCGMRESWithSaturation() {
  linearalgebra::DeleteVector(initial_guess_solution_);
  linearalgebra::DeleteVector(initial_guess_Lagrange_multiplier_);
  linearalgebra::DeleteVector(incremented_solution_vec_);
  linearalgebra::DeleteVector(solution_update_vec_);
  linearalgebra::DeleteVector(lambda_vec_);
  linearalgebra::DeleteVector(errors_in_optimality_);
  linearalgebra::DeleteVector(errors_in_optimality_1_);
  linearalgebra::DeleteVector(errors_in_optimality_2_);
}

void InitMSCGMRESWithSaturation::setInitParameters(
    const double* initial_guess_solution, 
    const double* initial_guess_Lagrange_multiplier,
    const double residual_tolerance, const int max_iteration, 
    const double finite_difference_step, const int kmax) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_guess_solution_[i] = initial_guess_solution[i];
  }
  for (int i=0; i<dim_saturation_; ++i) {
    initial_guess_Lagrange_multiplier_[i] = 
        initial_guess_Lagrange_multiplier[i];
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

void InitMSCGMRESWithSaturation::solveOCPForInit(
    const double initial_time, const double* initial_state_vec, 
    double* initial_solution_vec, 
    double* errors_for_control_input_and_constraints, 
    double* errors_for_dummy_input, 
    double* errors_for_control_input_sauturations) {
  // Substitute initial guess solution to solution_vec.
  for (int i=0; i<dim_control_input_and_constraints_; ++i) {
      initial_solution_vec[i] = initial_guess_solution_[i];
  }
  for (int i=0; i<dim_saturation_; ++i) {
    if (initial_solution_vec[saturation_seq_.index(i)] < saturation_seq_.min(i)) {
      initial_solution_vec[dim_control_input_and_constraints_+i] =
          (saturation_seq_.min(i)+saturation_seq_.max(i)) / 2;
    }
    else if (initial_solution_vec[saturation_seq_.index(i)] > saturation_seq_.max(i)) {
      initial_solution_vec[dim_control_input_and_constraints_+i] = 
          (saturation_seq_.min(i)+saturation_seq_.max(i)) / 2;
    }
    else {
      double max_plus_min = saturation_seq_.max(i) + saturation_seq_.min(i);
      double max_minus_min = saturation_seq_.max(i) - saturation_seq_.min(i);
      initial_solution_vec[dim_control_input_and_constraints_+i] = 
          std::sqrt((max_minus_min*max_minus_min)/4-
              (initial_guess_solution_[saturation_seq_.index(i)]-max_plus_min/2)
              *(initial_guess_solution_[saturation_seq_.index(i)]-max_plus_min/2));
    }
  }
  for (int i=0; i<dim_saturation_; ++i) {
    initial_solution_vec[dim_control_input_and_constraints_+dim_saturation_+i] = 
        initial_guess_Lagrange_multiplier_[i];
  }
  computeErrorsInOptimality(initial_time, initial_state_vec, 
                            initial_solution_vec, errors_in_optimality_);
  int itr_counter = 0;
  while (
      linearalgebra::SquaredNorm(dim_solution_, errors_in_optimality_) 
      > residual_tolerance_*residual_tolerance_ && itr_counter < max_iteration_) {
    solveGMRES(initial_time, initial_state_vec, initial_solution_vec, 
                 solution_update_vec_);
    for (int j=0; j<dim_solution_; ++j) {
      initial_solution_vec[j] += solution_update_vec_[j];
    }
    computeErrorsInOptimality(initial_time, initial_state_vec, 
                              initial_solution_vec, errors_in_optimality_);
    ++itr_counter;
  }
  for (int i=0; i<dim_control_input_and_constraints_; ++i) {
    errors_for_control_input_and_constraints[i] = errors_in_optimality_[i];
  }
  for (int i=0; i<dim_saturation_; ++i) {
    errors_for_dummy_input[i] = 
        errors_in_optimality_[dim_control_input_and_constraints_+i];
  }
  for (int i=0; i<dim_saturation_; ++i) {
    errors_for_control_input_sauturations[i] =
        errors_in_optimality_[dim_control_input_and_constraints_+dim_saturation_+i];
  }
}

inline void InitMSCGMRESWithSaturation::computeErrorsInOptimality(
    const double time, const double* state_vec, const double* solution_vec, 
    double* errors_in_optimality) {
  model_.phixFunc(time, state_vec, lambda_vec_);
  model_.huFunc(time, state_vec, solution_vec, lambda_vec_, 
                errors_in_optimality);
  condensingfunctions::addHamiltonianDerivativeWithConstrainedInput(
    saturation_seq_, solution_vec,
    &(solution_vec[dim_control_input_and_constraints_+dim_saturation_]), 
    errors_in_optimality);
  condensingfunctions::computeErrorsForDummyInput(
    saturation_seq_,
    &(solution_vec[dim_control_input_and_constraints_]),
    &(solution_vec[dim_control_input_and_constraints_+dim_saturation_]),
    &(errors_in_optimality[dim_control_input_and_constraints_])
  );
  condensingfunctions::computeErrorsForSaturation(
    saturation_seq_, solution_vec,
    &(solution_vec[dim_control_input_and_constraints_]), 
    &(errors_in_optimality[dim_control_input_and_constraints_+dim_saturation_])
  );
}

void InitMSCGMRESWithSaturation::bFunc(const double time, 
                                       const double* state_vec, 
                                       const double* current_solution_vec, 
                                       double* b_vec) {
  computeErrorsInOptimality(time, state_vec, current_solution_vec, 
                            errors_in_optimality_);
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_step_ * solution_update_vec_[i];
  }
  computeErrorsInOptimality(time, state_vec, incremented_solution_vec_, 
                            errors_in_optimality_1_);
  for (int i=0; i<dim_solution_; ++i) {
    b_vec[i] = - errors_in_optimality_[i] 
        - (errors_in_optimality_1_[i]-errors_in_optimality_[i]) 
        / finite_difference_step_;
  }
}

void InitMSCGMRESWithSaturation::axFunc(const double time, 
                                        const double* state_vec, 
                                        const double* current_solution_vec, 
                                        const double* direction_vec, 
                                        double* ax_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_step_ * direction_vec[i];
  }
  computeErrorsInOptimality(time, state_vec, incremented_solution_vec_, 
                            errors_in_optimality_1_);
  for (int i=0; i<dim_solution_; ++i) {
    ax_vec[i] = (errors_in_optimality_1_[i]-errors_in_optimality_[i]) 
        / finite_difference_step_;
  }
}