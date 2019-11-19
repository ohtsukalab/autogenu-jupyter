#ifndef SINGLE_SHOOTING_OCP_H
#define SINGLE_SHOOTING_OCP_H

#include "optimal_control_problem.hpp"
#include "time_varying_smooth_horizon.hpp"
#include "linear_algebra.hpp"

class SingleShootingOCP final : public OptimalControlProblem {
public:
  SingleShootingOCP(const double T_f, const double alpha, const int N);
  SingleShootingOCP(const double T_f, const double alpha, const int N,
                    const double initial_time);
  ~SingleShootingOCP();

  void computeOptimalityResidual(const double time, const double* state_vec, 
                                 const double* solution_vec,
                                 double* optimality_residual);

  void predictStateFromSolution(const double current_time, 
                                const double* current_state,
                                const double* solution_vec, 
                                const double prediction_length,
                                double* predicted_state);

  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);
  void resetHorizonLength(const double initial_time);

  int dim_solution() const override;
  int N() const;

private:
  TimeVaryingSmoothHorizon horizon_;
  int dim_solution_, N_;
  double *dx_vec_, **state_mat_, **lambda_mat_;

};

#endif // SINGLE_SHOOTING_OCP_H