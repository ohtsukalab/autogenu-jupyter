#include "continuation_gmres.hpp"


ContinuationGMRES::ContinuationGMRES(const double T_f, const double alpha, 
                                     const int N, const double zeta, 
                                     const double finite_difference_step, 
                                     const int kmax) 
  : MatrixFreeGMRES(), 
    model_(), 
    cgmres_initializer_(),
    dim_state_(model_.dimState()), 
    dim_control_input_(model_.dimControlInput()), 
    dim_constraints_(model_.dimConstraints()), 
    dim_control_input_and_constraints_(
        model_.dimControlInput()+model_.dimConstraints()), 
    dim_solution_(N*(model_.dimControlInput()+model_.dimConstraints())), 
    T_f_(T_f), 
    alpha_(alpha), 
    N_(N), 
    zeta_(zeta), 
    finite_difference_step_(finite_difference_step), 
    kmax_(kmax), 
    initial_time_(0), 
    incremented_time_(0), 
    dx_vec_(linearalgebra::NewVector(dim_state_)), 
    incremented_state_vec_(linearalgebra::NewVector(dim_state_)),
    solution_vec_(linearalgebra::NewVector(dim_solution_)), 
    incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)), 
    errors_in_optimality_(linearalgebra::NewVector(dim_solution_)), 
    errors_in_optimality_1_(linearalgebra::NewVector(dim_solution_)), 
    errors_in_optimality_2_(linearalgebra::NewVector(dim_solution_)), 
    solution_update_vec_(linearalgebra::NewVector(dim_solution_)), 
    state_mat_(linearalgebra::NewMatrix(N+1, dim_state_)), 
    lambda_mat_(linearalgebra::NewMatrix(N+1, dim_state_)) {
  // Set dimensions and parameters used in GMRES.
  setGMRESParameters(dim_solution_, kmax); 
}

ContinuationGMRES::~ContinuationGMRES() {
  linearalgebra::DeleteVector(dx_vec_);
  linearalgebra::DeleteVector(incremented_state_vec_);
  linearalgebra::DeleteVector(solution_vec_);
  linearalgebra::DeleteVector(incremented_solution_vec_);
  linearalgebra::DeleteVector(errors_in_optimality_);
  linearalgebra::DeleteVector(errors_in_optimality_1_);
  linearalgebra::DeleteVector(errors_in_optimality_2_);
  linearalgebra::DeleteVector(solution_update_vec_);
  linearalgebra::DeleteMatrix(state_mat_);
  linearalgebra::DeleteMatrix(lambda_mat_);
}

void ContinuationGMRES::setInitParameters(const double* initial_guess_solution,
                                          const double residual_tolerance, 
                                          const int max_iteration, 
                                          const double finite_difference_step, 
                                          const int kmax) {
  cgmres_initializer_.setInitParameters(initial_guess_solution, 
                                        residual_tolerance, max_iteration, 
                                        finite_difference_step, kmax);
}

void ContinuationGMRES::initSolution(const double initial_time, 
                                     const double* initial_state_vec, 
                                     double* optimal_control_input_vec) {
  double initial_solution_vec[dim_control_input_and_constraints_], 
      initial_errors_in_optimality[dim_control_input_and_constraints_];
  initial_time_ = initial_time;
  cgmres_initializer_.solveOCPForInit(initial_time, initial_state_vec, 
                                      initial_solution_vec, 
                                      initial_errors_in_optimality);
  for (int i=0; i<N_; ++i) {
    int i_total = i * dim_control_input_and_constraints_;
    for (int j=0; j<dim_control_input_and_constraints_; ++j) {
      solution_vec_[i_total+j] = initial_solution_vec[j];
    }
    // Intialize the errors_in_optimality_.
    for (int j=0; j<dim_control_input_and_constraints_; ++j) {
      errors_in_optimality_[i_total+j] = initial_errors_in_optimality[j];
    }
  }
  for (int i=0; i<dim_control_input_; ++i) {
      optimal_control_input_vec[i] = initial_solution_vec[i];
  }
}

void ContinuationGMRES::controlUpdate(const double time, 
                                      const double sampling_period, 
                                      const double* state_vec, 
                                      double* optimal_control_input_vec) {
  // Predict the incremented state.
  incremented_time_ = time + finite_difference_step_;
  model_.stateFunc(time, state_vec, solution_vec_, dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    incremented_state_vec_[i] = state_vec[i] 
        + finite_difference_step_ * dx_vec_[i];
  }
  solveGMRES(time, state_vec, solution_vec_, solution_update_vec_);
  for (int i=0; i<dim_solution_; ++i) {
    solution_vec_[i] += sampling_period * solution_update_vec_[i];
  }
  for (int i=0; i<dim_control_input_; ++i) {
    optimal_control_input_vec[i] = solution_vec_[i];
  }
}

double ContinuationGMRES::getErrorNorm(const double time, 
                                       const double* state_vec) {
  double error_vec[dim_solution_];
  computeErrorsInOptimality(time, state_vec, solution_vec_, error_vec);
  return std::sqrt(linearalgebra::SquaredNorm(dim_solution_, error_vec));
}

void ContinuationGMRES::computeErrorsInOptimality(
    const double time, const double* state_vec, 
    const double* solution_vec, double* errors_in_optimality) {
    // Set the length of the horizon and discretize the horizon.
  double horizon_length = T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
  double delta_tau = horizon_length / N_;

  // Compute the state trajectory over the horizon on the basis of the 
  // time, solution_vec and the state_vec.
  for (int i=0; i<dim_state_; ++i) {
    state_mat_[0][i] = state_vec[i];
  }
  double tau = time;
  for (int i=0; i<N_; ++i, tau+=delta_tau) {
    model_.stateFunc(
        tau, state_mat_[i], 
        &(solution_vec[i*dim_control_input_and_constraints_]), dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      state_mat_[i+1][j] = state_mat_[i][j] + delta_tau * dx_vec_[j];
    }
  }
  // Compute the Lagrange multiplier over the horizon on the basis of 
  // time, solution_vec and the state_vec.
  model_.phixFunc(tau, state_mat_[N_], lambda_mat_[N_]);
  for (int i=N_-1; i>=0; --i, tau-=delta_tau) {
    model_.hxFunc(
        tau, state_mat_[i], 
        &(solution_vec[i*dim_control_input_and_constraints_]), 
        lambda_mat_[i+1], dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat_[i][j] = lambda_mat_[i+1][j] + delta_tau * dx_vec_[j];
    }
  }
  // Compute the erros in optimality over the horizon on the basis of the 
  // control_input_vec and the state_vec.
  tau = time;
  for (int i=0; i<N_; ++i, tau+=delta_tau) {
    model_.huFunc(
        tau, state_mat_[i], 
        &(solution_vec[i*dim_control_input_and_constraints_]), 
        lambda_mat_[i+1], 
        &(errors_in_optimality[i*dim_control_input_and_constraints_]));
  }
}

void ContinuationGMRES::bFunc(const double time, const double* state_vec,
                              const double* current_solution_vec, 
                              double* b_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_step_ * solution_update_vec_[i];
  }
  computeErrorsInOptimality(time, state_vec, current_solution_vec, 
                            errors_in_optimality_);
  computeErrorsInOptimality(incremented_time_, incremented_state_vec_, 
                            current_solution_vec, errors_in_optimality_1_);
  computeErrorsInOptimality(incremented_time_, incremented_state_vec_, 
                            incremented_solution_vec_, errors_in_optimality_2_);
  for (int i=0; i<dim_solution_; ++i) {
    b_vec[i] = (1/finite_difference_step_-zeta_) * errors_in_optimality_[i] 
        - errors_in_optimality_2_[i] / finite_difference_step_;
  }
}

inline void ContinuationGMRES::axFunc(const double time, 
                                      const double* state_vec, 
                                      const double* current_solution_vec, 
                                      const double* direction_vec, 
                                      double* ax_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_step_ * direction_vec[i];
  }
  computeErrorsInOptimality(incremented_time_, incremented_state_vec_, 
                            incremented_solution_vec_, errors_in_optimality_2_);
  for (int i=0; i<dim_solution_; ++i) {
    ax_vec[i] = 
        (errors_in_optimality_2_[i]-errors_in_optimality_1_[i]) 
        / finite_difference_step_;
  }
}