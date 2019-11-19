#ifndef MSOCP_WITH_INPUT_SATURATION_H
#define MSOCP_WITH_INPUT_SATURATION_H

#include "input_saturation_functions.hpp"
#include "input_saturation_set.hpp"
#include "optimal_control_problem.hpp"
#include "time_varying_smooth_horizon.hpp"
#include "linear_algebra.hpp"

class MSOCPWithInputSaturation final : public OptimalControlProblem {
public:
  MSOCPWithInputSaturation(const InputSaturationSet& input_saturation_set, 
                           const double T_f, const double alpha, const int N);
  MSOCPWithInputSaturation(const InputSaturationSet& input_saturation_set, 
                           const double T_f, const double alpha, const int N,
                           const double initial_time);
  ~MSOCPWithInputSaturation();

  void computeOptimalityResidualForControlInputAndConstraints(
      const double time, const double* state_vec, 
      const double* control_input_and_constraints_seq, 
      double const* const* state_mat, double const* const* lambda_mat, 
      double const* const* input_saturation_multiplier_mat,
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

  void computeResidualForDummyInputAndInputSaturation(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* input_saturation_multiplier_mat, 
      double** optimality_residual_for_dummy_input, 
      double** optimality_residual_for_input_saturation);

  void multiplyResidualForDummyInputAndInputSaturationInverse(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* saturation_Lagrange_multiplier_mat, 
      double const* const* multiplied_dummy_input_mat, 
      double const* const* multiplied_Lagrange_multiplier_mat, 
      double** resulted_dummy_input_mat, 
      double** resulted_Lagrange_multiplier_mat);

  void computeResidualDifferenceForDummyInput(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      const double* control_input_and_constraints_update_seq, 
      double** dummy_residual_difference_mat);

  void computeResidualDifferenceForInputSaturation(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* input_saturation_multiplier_mat, 
      const double* control_input_and_constraints_update_seq, 
      double** input_saturation_residual_difference_mat);

  void predictStateFromSolution(const double current_time, 
                                const double* current_state,
                                const double* solution_vec, 
                                const double prediction_length,
                                double* predicted_state);

  void resetHorizonLength(const double initial_time);

  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  int dim_solution() const override;
  int dim_saturation() const;
  int N() const;

private:
  TimeVaryingSmoothHorizon horizon_;
  InputSaturationSet input_saturation_set_;
  int dim_solution_, dim_saturation_, N_;
  double *dx_vec_;
};

#endif // MSOCP_WITH_INPUT_SATURATION_H