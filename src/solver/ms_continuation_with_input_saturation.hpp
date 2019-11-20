// This class provides the linear problem of the continuation transformation 
// for the multiple-shooting optimal control problem, which is solved in 
// Matrix-free GMRES. 

#ifndef MS_CONTINUATION_WITH_INPUT_SATURATION_H
#define MS_CONTINUATION_WITH_INPUT_SATURATION_H

#include <cmath>
#include "linear_algebra.hpp"
#include "input_saturation_set.hpp"
#include "ms_ocp_with_input_saturation.hpp"

// Linear problem of the continuation transformation for the multiple-shooting 
// optimal control problem, which is solved in Matrix-free GMRES. This class 
// is intended for use with MatrixfreeGMRES class. 
class MSContinuationWithInputSaturation {
public:
  // Constructs MSContinuationWithInputSaturation with setting parameters and 
  // allocates vectors and matrices.
  // Arguments:
  //  InputSaturationSet: The set composed of the input saturation constraints.
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  MSContinuationWithInputSaturation(
      const InputSaturationSet& input_saturation_set, const double T_f, 
      const double alpha, const int N, const double finite_difference_increment,
      const double zeta);

  // Constructs MSContinuationWithInputSaturation with setting parameters and 
  // allocates vectors and matrices.
  // Arguments:
  //  InputSaturationSet: The set composed of the input saturation constraints.
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  MSContinuationWithInputSaturation(
      const InputSaturationSet& input_saturation_set, const double T_f, 
      const double alpha, const int N, const double initial_time, 
      const double finite_difference_increment, const double zeta);

  // Free vectors and matrices.
  ~MSContinuationWithInputSaturation();

  // Integrates the solution for given optimal update vector of the solution 
  // and the integration length.
  void integrateSolution(double* control_input_and_constraints_seq, 
                         double** state_mat, double** lambda_mat,
                         double** dummy_input_mat, 
                         double** input_saturation_multiplier_mat,
                         const double* control_input_and_constraints_update_seq, 
                         const double integration_length);

  // Computes and returns the squared norm of the errors in optimality under 
  // the state_vec and the current solution.
  double computeErrorNorm(const double time, const double* state_vec, 
                          const double* control_input_and_constraints_seq,
                          double const* const* state_mat, 
                          double const* const* lambda_mat,
                          double const* const* dummy_input_mat, 
                          double const* const* input_saturation_multiplier_mat);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double initial_time);

  // Computes a vector correspongin to b in Ax=b. This function is called in
  // MatrixfreeGMRES.
  void bFunc(const double time, const double* state_vec, 
             const double* control_input_and_constraints_seq, 
             double const* const* state_mat, double const* const* lambda_mat,
             double const* const* dummy_input_mat, 
             double const* const* input_saturation_multiplier_mat,
             const double* current_control_input_and_constraints_update_seq, 
             double* b_vec);

  // Computes a vector correspongin to Ax in Ax=b. This function is called in
  // MatrixfreeGMRES.
  void AxFunc(const double time, const double* state_vec, 
              const double* control_input_and_constraints_seq, 
              double const* const* state_mat, double const* const* lambda_mat,
              double const* const* dummy_input_mat, 
              double const* const* input_saturation_multiplier_mat,
              const double* direction_vec, double* ax_vec);

  // Returns the dimension of the state.
  int dim_state() const;

  // Returns the dimension of the control input.
  int dim_control_input() const;

  // Returns the dimension of the equality constraints.
  int dim_constraints() const;

  int dim_saturation() const;

  // Returns the dimension of the solution, which is equivalent to 
  // N*(dim_control_input+dim_constraints).
  int dim_condensed_problem() const;

  int N() const;

  MSContinuationWithInputSaturation(const MSContinuationWithInputSaturation&) 
      = delete;
  MSContinuationWithInputSaturation& operator=(
        const MSContinuationWithInputSaturation&) = delete;

private:
  MSOCPWithInputSaturation ocp_;
  const int dim_state_, dim_control_input_, dim_constraints_, 
      dim_control_input_and_constraints_, dim_saturation_,
      dim_control_input_and_constraints_seq_, N_;
  double finite_difference_increment_, zeta_, incremented_time_; 
  double *incremented_state_vec_, 
      *incremented_control_input_and_constraints_seq_, 
      *control_input_and_constraints_residual_seq_, 
      *control_input_and_constraints_residual_seq_1_, 
      *control_input_and_constraints_residual_seq_2_, 
      *control_input_and_constraints_residual_seq_3_;
  double **incremented_state_mat_, **incremented_lambda_mat_, 
      **state_residual_mat_, **state_residual_mat_1_, 
      **lambda_residual_mat_, **lambda_residual_mat_1_,
      **incremented_dummy_input_mat_, 
      **incremented_input_sautration_multiplier_mat_,
      **dummy_input_residual_mat_, **dummy_input_residual_mat_1_, 
      **input_saturation_residual_mat_, **input_saturation_residual_mat_1_,
      **dummy_input_difference_mat_,
      **input_saturation_multiplier_difference_mat_;
};

#endif // MS_CONTINUATION_WITH_INPUT_SATURATION_H