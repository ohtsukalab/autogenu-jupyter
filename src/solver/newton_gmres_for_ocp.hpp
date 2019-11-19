#ifndef NEWTON_GMRES_FOR_OCP_H
#define NEWTON_GMRES_FOR_OCP_H

#include <cmath>
#include "linear_algebra.hpp"

template <class OCPType, typename... OCPConstructorArgs>
class NewtonGMRESForOCP {
public:
  NewtonGMRESForOCP(const double finite_difference_increment, 
                    OCPConstructorArgs... ocp_constructor_args) 
    : ocp_(ocp_constructor_args...), 
      dim_solution_(ocp_.dim_solution()),
      finite_difference_increment_(finite_difference_increment),
      incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)),
      optimality_residual_(linearalgebra::NewVector(dim_solution_)),
      optimality_residual_1_(linearalgebra::NewVector(dim_solution_)) {
  }

  ~NewtonGMRESForOCP() {
    linearalgebra::DeleteVector(incremented_solution_vec_);
    linearalgebra::DeleteVector(optimality_residual_);
    linearalgebra::DeleteVector(optimality_residual_1_);
  }

  double errorNorm(const double time, const double* state_vec, 
                   const double* solution_vec) {
    ocp_.computeOptimalityResidual(time, state_vec, solution_vec, 
                                  optimality_residual_);
    return std::sqrt(
          linearalgebra::SquaredNorm(dim_solution_, optimality_residual_));
  }

  void bFunc(const double time, const double* state_vec, 
             const double* current_solution_vec, 
             const double* current_solution_update_vec, double* b_vec) {
    for (int i=0; i<dim_solution_; ++i) {
      incremented_solution_vec_[i] = current_solution_vec[i] 
          + finite_difference_increment_ * current_solution_update_vec[i];
    }
    ocp_.computeOptimalityResidual(time, state_vec, current_solution_vec, 
                                  optimality_residual_);
    ocp_.computeOptimalityResidual(time, state_vec, incremented_solution_vec_, 
                                  optimality_residual_1_);
    for (int i=0; i<dim_solution_; ++i) {
      b_vec[i] = - optimality_residual_[i] 
          - (optimality_residual_1_[i]-optimality_residual_[i]) 
          / finite_difference_increment_;
    }
  }

  void AxFunc(const double time, const double* state_vec, 
              const double* current_solution_vec, const double* direction_vec,
              double* ax_vec) {
    for (int i=0; i<dim_solution_; ++i) { 
      incremented_solution_vec_[i] = current_solution_vec[i] 
          + finite_difference_increment_ * direction_vec[i];
    }
    ocp_.computeOptimalityResidual(time, state_vec, incremented_solution_vec_, 
                                  optimality_residual_1_);
    for (int i=0; i<dim_solution_; ++i) {
      ax_vec[i] = (optimality_residual_1_[i]-optimality_residual_[i]) 
          / finite_difference_increment_;
    }
  }

  void getTerminalCostDerivatives(const double time, 
                                  const double* state_vec, 
                                  double* terminal_cost_derivative_vec) {
    ocp_.computeTerminalCostDerivative(time, state_vec, 
                                      terminal_cost_derivative_vec);
  }

  int dim_state() const {
    return ocp_.dim_state();
  }

  int dim_control_input() const {
    return ocp_.dim_control_input();
  }

  int dim_constraints() const {
    return ocp_.dim_constraints();
  }

  int dim_solution() const {
    return ocp_.dim_solution();
  }

private:
  OCPType ocp_;
  const int dim_solution_;
  double finite_difference_increment_; 
  double *incremented_solution_vec_, *optimality_residual_, 
      *optimality_residual_1_;

};

#endif // NEWTON_GMRES_FOR_OCP_H