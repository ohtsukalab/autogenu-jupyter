#ifndef ZERO_HORIZON_OCP_H
#define ZERO_HORIZON_OCP_H

#include "optimal_control_problem.hpp"
#include "linear_algebra.hpp"

class ZeroHorizonOCP final : public OptimalControlProblem {
public:
  ZeroHorizonOCP();
  ~ZeroHorizonOCP();

  void computeOptimalityResidual(const double time, const double* state_vec, 
                                 const double* solution_vec,
                                 double* optimality_residual);

  void computeTerminalCostDerivative(const double time, const double* state_vec,
                                     double* terminal_cost_derivative_vec);

  int dim_solution() const override;

private:
  int dim_solution_;
  double *lambda_vec_;

};

#endif // ZERO_HORIZON_OCP_H