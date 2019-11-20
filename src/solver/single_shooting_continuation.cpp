#include "single_shooting_continuation.hpp"

SingleShootingContinuation::SingleShootingContinuation(
    const double T_f, const double alpha, const int N,
    const double finite_difference_increment, const double zeta) 
  : ocp_(T_f, alpha, N),
    dim_state_(ocp_.dim_state()),
    dim_control_input_(ocp_.dim_control_input()),
    dim_constraints_(ocp_.dim_constraints()),
    dim_solution_(N*(dim_control_input_+dim_constraints_)),
    finite_difference_increment_(finite_difference_increment),
    zeta_(zeta),
    incremented_time_(0),
    incremented_state_vec_(linearalgebra::NewVector(ocp_.dim_state())),
    incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_1_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_2_(linearalgebra::NewVector(dim_solution_)) {
}

SingleShootingContinuation::SingleShootingContinuation(
    const double T_f, const double alpha, const int N, 
    const double initial_time,
    const double finite_difference_increment, const double zeta) 
  : ocp_(T_f, alpha, N, initial_time),
    dim_state_(ocp_.dim_state()),
    dim_control_input_(ocp_.dim_control_input()),
    dim_constraints_(ocp_.dim_constraints()),
    dim_solution_(N*(dim_control_input_+dim_constraints_)),
    finite_difference_increment_(finite_difference_increment),
    zeta_(zeta),
    incremented_time_(0),
    incremented_state_vec_(linearalgebra::NewVector(ocp_.dim_state())),
    incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_1_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_2_(linearalgebra::NewVector(dim_solution_)) {
}

SingleShootingContinuation::~SingleShootingContinuation() {
  linearalgebra::DeleteVector(incremented_state_vec_);
  linearalgebra::DeleteVector(incremented_solution_vec_);
  linearalgebra::DeleteVector(optimality_residual_);
  linearalgebra::DeleteVector(optimality_residual_1_);
  linearalgebra::DeleteVector(optimality_residual_2_);
}

void SingleShootingContinuation::integrateSolution(
    double* solution_vec, const double* solution_update_vec, 
    const double integration_length) {
  for (int i=0; i<dim_solution_; ++i) {
    solution_vec[i] += integration_length * solution_update_vec[i];
  }
}

double SingleShootingContinuation::computeErrorNorm(const double time, 
                                                    const double* state_vec,
                                                    const double* solution_vec) {
  ocp_.computeOptimalityResidual(time, state_vec, solution_vec,
                                 optimality_residual_);
  return std::sqrt(
        linearalgebra::SquaredNorm(dim_solution_, optimality_residual_));
}

void SingleShootingContinuation::resetHorizonLength(const double T_f, 
                                                    const double alpha, 
                                                    const double initial_time) {
  ocp_.resetHorizonLength(T_f, alpha, initial_time);
}

void SingleShootingContinuation::resetHorizonLength(const double initial_time) {
  ocp_.resetHorizonLength(initial_time);
}

void SingleShootingContinuation::bFunc(const double time, 
                                       const double* state_vec, 
                                       const double* current_solution_vec, 
                                       const double* current_solution_update_vec, 
                                       double* b_vec) {
  incremented_time_ = time + finite_difference_increment_;
  ocp_.predictStateFromSolution(time, state_vec, current_solution_vec,
                                finite_difference_increment_, 
                                incremented_state_vec_);
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_increment_ * current_solution_update_vec[i];
  }
  ocp_.computeOptimalityResidual(time, state_vec, current_solution_vec, 
                                 optimality_residual_);
  ocp_.computeOptimalityResidual(incremented_time_, incremented_state_vec_, 
                                 current_solution_vec, optimality_residual_1_);
  ocp_.computeOptimalityResidual(incremented_time_, incremented_state_vec_, 
                                 incremented_solution_vec_, 
                                 optimality_residual_2_);
  for (int i=0; i<dim_solution_; ++i) {
    b_vec[i] = (1/finite_difference_increment_-zeta_) * optimality_residual_[i] 
        - optimality_residual_2_[i] / finite_difference_increment_;
  }
}

void SingleShootingContinuation::AxFunc(const double time, 
                                        const double* state_vec, 
                                        const double* current_solution_vec, 
                                        const double* direction_vec, 
                                        double* ax_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_increment_ * direction_vec[i];
  }
  ocp_.computeOptimalityResidual(incremented_time_, incremented_state_vec_, 
                                 incremented_solution_vec_, 
                                 optimality_residual_2_);
  for (int i=0; i<dim_solution_; ++i) {
    ax_vec[i] = 
        (optimality_residual_2_[i]-optimality_residual_1_[i]) 
        / finite_difference_increment_;
  }
}

int SingleShootingContinuation::dim_state() const {
  return dim_state_;
}

int SingleShootingContinuation::dim_control_input() const {
  return dim_control_input_;
}

int SingleShootingContinuation::dim_constraints() const {
  return dim_constraints_;
}

int SingleShootingContinuation::dim_solution() const {
  return dim_solution_;
}

int SingleShootingContinuation::N() const {
  return ocp_.N();
}