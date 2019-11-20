// Privides the optimal control problem (OCP) with horizon whose length is zero.
// The OCP also consider the constrains provided by InputSaturationSet.

#ifndef ZERO_HORIZON_OCP_WITH_INPUT_SATURATION_H
#define ZERO_HORIZON_OCP_WITH_INPUT_SATURATION_H

#include "input_saturation_set.hpp"
#include "input_saturation_functions.hpp"
#include "optimal_control_problem.hpp"
#include "linear_algebra.hpp"

// Privides the optimal control problem (OCP) with horizon whose length is zero.
// The OCP also consider the constrains provided by InputSaturationSet.
class ZeroHorizonOCPWithInputSaturation final : public OptimalControlProblem {
public:
  // Allocates a vector.
  ZeroHorizonOCPWithInputSaturation(
      const InputSaturationSet& input_saturation_set);

  // Free a vector.
  ~ZeroHorizonOCPWithInputSaturation();

  // Computes the optimaliy residual under time, state_vec, and solution_vec 
  // that represents the control input and Lgrange multiplier with respect to
  // equality constraints. The result is set in optimality_residual.
  void computeOptimalityResidual(const double time, const double* state_vec, 
                                 const double* solution_vec,
                                 double* optimality_residual);

  // Computes the partial derivative of the terminal cost with respect to
  // the state.
  void computeTerminalCostDerivative(const double time, const double* state_vec,
                                     double* terminal_cost_derivative_vec);

  // Returns the number of the control inputs that are constrained by the 
  // INputSaturationSet.
  int dim_saturation() const;

  // Return the dimension of the solution, 
  // i.e., dim_control_input+dim_constraints.
  int dim_solution() const override;

private:
  InputSaturationSet input_saturation_set_;
  int dim_solution_, dim_saturation_;
  double *lambda_vec_;
};

#endif // ZERO_HORIZON_OCP_WITH_INPUT_SATURATION_H