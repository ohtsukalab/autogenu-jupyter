// The multiple-shooting two-point boundary-value problem of the
// finite-horizon optimal control problem. Functions for condensing of the 
// solution are also provided.

#ifndef MULTIPLE_SHOOTING_OCP_H
#define MULTIPLE_SHOOTING_OCP_H

#include "optimal_control_problem.hpp"
#include "time_varying_smooth_horizon.hpp"
#include "linear_algebra.hpp"

// This class provides multiple-shooting two-point boundary-value problem of 
// the finite-horizon optimal control problem. Functions for condensing of the 
// solution are also provided.
class MultipleShootingOCP final : public OptimalControlProblem {
public:
  // Constructs MultipleShootingOCP with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  MultipleShootingOCP(const double T_f, const double alpha, const int N);

  // Constructs MultipleShootingOCP with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  MultipleShootingOCP(const double T_f, const double alpha, const int N,
                      const double initial_time);

  // Free vectors and matrices.
  ~MultipleShootingOCP();

  // Computes the optimaliy residual with respect to the control input and the 
  // equality constraints under time, state_vec, and solution_vec 
  // that represents the control input sequence. The result is stored in 
  // optimality_redisual_for_control_input_and_constraints.
  void computeOptimalityResidualForControlInputAndConstraints(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double* optimality_redisual_for_control_input_and_constraints);

  // Computes the optimaliy residual with respect to the state and lambda,
  // the Lagrange multiplier with respect to the state equation
  // under time, state_vec, and solution_vec that represents the control 
  // input sequence. The former is stored in optimality_residual_for_state
  // and the latter optimality_residual_for_lambda.
  void computeOptimalityResidualForStateAndLambda(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* state_mat, double const* const* lambda_mat, 
    double** optimality_residual_for_state, 
    double** optimality_residual_for_lambda);

  // Computes the state and lambda, the Lagrange multiplier with respect to 
  // the state equation from the optimality residual with respect to the state 
  // and Lambda. This function is needed for condensing of the solution of the 
  // multiple-shooting based methods.
  void computeStateAndLambdaFromOptimalityResidual(
    const double time, const double* state_vec, 
    const double* control_input_and_constraints_seq, 
    double const* const* optimality_residual_for_state,
    double const* const* optimality_residual_for_lambda,
    double** state_mat, double** lambda_mat);

  // Predicts the state in the finite future. Under time,state_vec, and 
  // solution_vec that represents the control input sequence., and 
  // prediction_length The result is set in predicted_state.
  void predictStateFromSolution(const double current_time, 
                                const double* current_state,
                                const double* solution_vec, 
                                const double prediction_length,
                                double* predicted_state);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double initial_time);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  // Returns the dimension of the solution, which is equivalent to 
  // N*(dim_control_input+dim_constraints).
  int dim_solution() const override;

  // Returns the grid number of the horizon.
  int N() const;

private:
  TimeVaryingSmoothHorizon horizon_;
  int dim_solution_, N_;
  double *dx_vec_;
};

#endif // MULTIPLE_SHOOTING_OCP_H