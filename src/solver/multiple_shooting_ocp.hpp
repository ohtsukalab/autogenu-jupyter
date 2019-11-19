#ifndef MULTIPLE_SHOOTING_OCP_H
#define MULTIPLE_SHOOTING_OCP_H

#include "optimal_control_problem.hpp"
#include "time_varying_smooth_horizon.hpp"
#include "linear_algebra.hpp"

class MultipleShootingOCP final : public OptimalControlProblem {
public:
  MultipleShootingOCP(const double T_f, const double alpha, const int N);
  MultipleShootingOCP(const double T_f, const double alpha, const int N,
                      const double initial_time);
  ~MultipleShootingOCP();

  void computeOptimalityResidualForControlInputAndConstraints(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double* optimality_redisual_for_control_input_and_constraints);

  void computeOptimalityResidualForStateAndLambda(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double** optimality_residual_for_state, 
    double** optimality_residual_for_lambda);

  void computeStateAndLambdaFromOptimalityResidual(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* optimality_residual_for_state,
    double const* const* optimality_residual_for_lambda,
    double** state_mat, double** lambda_mat);

  void predictStateFromSolution(const double current_time, 
                                const double* current_state,
                                const double* solution_vec, 
                                const double prediction_length,
                                double* predicted_state);

  void resetHorizonLength(const double initial_time);

  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  int dim_solution() const override;
  int N() const;

private:
  TimeVaryingSmoothHorizon horizon_;
  int dim_solution_, N_;
  double *dx_vec_;

};

#endif // MULTIPLE_SHOOTING_OCP_H