#ifndef ZERO_HORIZON_OCP_WITH_INPUT_SATURATION_H
#define ZERO_HORIZON_OCP_WITH_INPUT_SATURATION_H

#include "input_saturation_set.hpp"
#include "input_saturation_functions.hpp"
#include "optimal_control_problem.hpp"
#include "linear_algebra.hpp"

class ZeroHorizonOCPWithInputSaturation final : public OptimalControlProblem {
public:
  ZeroHorizonOCPWithInputSaturation(
      const InputSaturationSet& input_saturation_set);
  ~ZeroHorizonOCPWithInputSaturation();

  void computeOptimalityResidual(const double time, const double* state_vec, 
                                 const double* solution_vec,
                                 double* optimality_residual);

  void computeTerminalCostDerivative(const double time, const double* state_vec,
                                     double* terminal_cost_derivative_vec);

  int dim_saturation() const;
  int dim_solution() const override;

private:
  InputSaturationSet input_saturation_set_;
  int dim_solution_, dim_saturation_;
  double *lambda_vec_;
};

#endif // ZERO_HORIZON_OCP_WITH_INPUT_SATURATION_H