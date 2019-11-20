// The multiple-shooting two-point boundary-value problem of the
// finite-horizon optimal control problem. Functions for condensing of the 
// solution are also provided.

#ifndef MSOCP_WITH_INPUT_SATURATION_H
#define MSOCP_WITH_INPUT_SATURATION_H

#include "input_saturation_functions.hpp"
#include "input_saturation_set.hpp"
#include "optimal_control_problem.hpp"
#include "time_varying_smooth_horizon.hpp"
#include "linear_algebra.hpp"

// This class provides multiple-shooting two-point boundary-value problem of 
// the finite-horizon optimal control problem. Functions for condensing of the 
// solution are also provided.
class MSOCPWithInputSaturation final : public OptimalControlProblem {
public:
  // Constructs MSOCPWithInputSaturation with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  InputSaturationSet: The set composed of the input saturation constraints.
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  MSOCPWithInputSaturation(const InputSaturationSet& input_saturation_set, 
                           const double T_f, const double alpha, const int N);

  // Constructs MSOCPWithInputSaturation with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  InputSaturationSet: The set composed of the input saturation constraints.
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  MSOCPWithInputSaturation(const InputSaturationSet& input_saturation_set, 
                           const double T_f, const double alpha, const int N,
                           const double initial_time);

  // Free vectors and matrices.
  ~MSOCPWithInputSaturation();

  // Computes the optimaliy residual with respect to the control input and the 
  // equality constraints under time, state_vec, and solution_vec 
  // that represents the control input sequence. The result is stored in 
  // optimality_redisual_for_control_input_and_constraints.
  void computeOptimalityResidualForControlInputAndConstraints(
      const double time, const double* state_vec, 
      const double* control_input_and_constraints_seq, 
      double const* const* state_mat, double const* const* lambda_mat, 
      double const* const* input_saturation_multiplier_mat,
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

  // Computes optimality residual for dummy input and constraints
  // on the saturation functions for the contorl input. 
  // The resulted errors are assigned in optimality_residual_for_dummy_input and
  // optimality_residual_for_input_saturation
  void computeResidualForDummyInputAndInputSaturation(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* input_saturation_multiplier_mat, 
      double** optimality_residual_for_dummy_input, 
      double** optimality_residual_for_input_saturation);

  // Computes the invers of the matrix of optimality residual for dummy input
  // and constraints on the saturation functions for the control input,
  // multiply that inverse matrix to dummy input matrix and Lagrange multiplier
  // mat. The multiplied matrices have to be assigned in
  // multiplied_dummy_input_mat and multiplied_Lagrange_multiplier_mat and 
  // the resulted matrices are assigned in resulted_dummy_input_mat and in 
  // resulted_Lagrange_multiplier_mat.
  void multiplyResidualForDummyInputAndInputSaturationInverse(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* saturation_Lagrange_multiplier_mat, 
      double const* const* multiplied_dummy_input_mat, 
      double const* const* multiplied_Lagrange_multiplier_mat, 
      double** resulted_dummy_input_mat, 
      double** resulted_Lagrange_multiplier_mat);

  // Computes the difference value of the dummy input corresponding to the
  // difference in the control input and Lagrange multiplier with respect to
  // the constraints that are not condensed. That update quantities have to 
  // be assigned in control_input_and_constraints_update_mat and the result is
  // assigned in dummy_residual_difference_mat.
  void computeResidualDifferenceForDummyInput(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      const double* control_input_and_constraints_update_seq, 
      double** dummy_residual_difference_mat);

  // Computes the difference value of the Lagrange multipliers with respect to
  // the constraints on the condensed saturation functions of the control
  // input corresponding to the difference in the control input and Lagrange 
  // multiplier with respect to the constraints that are not condensed. That 
  // update quantities have to be assigned in 
  // control_input_and_constraints_update_mat and the result is assigned in
  // input_saturation_residual_difference_mat.
  void computeResidualDifferenceForInputSaturation(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* input_saturation_multiplier_mat, 
      const double* control_input_and_constraints_update_seq, 
      double** input_saturation_residual_difference_mat);

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

  // Returns the number of the constrained control input.
  int dim_saturation() const;

  // Returns the grid number of the horizon.
  int N() const;

private:
  TimeVaryingSmoothHorizon horizon_;
  InputSaturationSet input_saturation_set_;
  int dim_solution_, dim_saturation_, N_;
  double *dx_vec_;
};

#endif // MSOCP_WITH_INPUT_SATURATION_H