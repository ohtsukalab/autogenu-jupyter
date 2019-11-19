#ifndef MS_CONTINUATION_WITH_INPUT_SATURATION_H
#define MS_CONTINUATION_WITH_INPUT_SATURATION_H

#include <cmath>
#include "linear_algebra.hpp"
#include "input_saturation_set.hpp"
#include "ms_ocp_with_input_saturation.hpp"

class MSContinuationWithInputSaturation {
public:
  MSContinuationWithInputSaturation(
      const InputSaturationSet& input_saturation_set, const double T_f, 
      const double alpha, const int N, const double finite_difference_increment,
      const double zeta);
  MSContinuationWithInputSaturation(
      const InputSaturationSet& input_saturation_set, const double T_f, 
      const double alpha, const int N, const double initial_time, 
      const double finite_difference_increment, const double zeta);
  ~MSContinuationWithInputSaturation();

  void integrateSolution(double* control_input_and_constraints_seq, 
                         double** state_mat, double** lambda_mat,
                         double** dummy_input_mat, 
                         double** input_saturation_multiplier_mat,
                         const double* control_input_and_constraints_update_seq, 
                         const double integration_length);

  double computeErrorNorm(const double time, const double* state_vec, 
                          const double* control_input_and_constraints_seq,
                          double const* const* state_mat, 
                          double const* const* lambda_mat,
                          double const* const* dummy_input_mat, 
                          double const* const* input_saturation_multiplier_mat);

  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  void resetHorizonLength(const double initial_time);

  void bFunc(const double time, const double* state_vec, 
             const double* control_input_and_constraints_seq, 
             double const* const* state_mat, double const* const* lambda_mat,
             double const* const* dummy_input_mat, 
             double const* const* input_saturation_multiplier_mat,
             const double* current_control_input_and_constraints_update_seq, 
             double* b_vec);

  void AxFunc(const double time, const double* state_vec, 
              const double* control_input_and_constraints_seq, 
              double const* const* state_mat, double const* const* lambda_mat,
              double const* const* dummy_input_mat, 
              double const* const* input_saturation_multiplier_mat,
              const double* direction_vec, double* ax_vec);

  int dim_state() const;
  int dim_control_input() const;
  int dim_constraints() const;
  int dim_saturation() const;
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