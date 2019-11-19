#ifndef MULTIPLE_SHOOTING_CONTINUATION_H
#define MULTIPLE_SHOOTING_CONTINUATION_H

#include <cmath>
#include "linear_algebra.hpp"
#include "multiple_shooting_ocp.hpp"

class MultipleShootingContinuation {
public:
  MultipleShootingContinuation(const double T_f, const double alpha, 
                               const int N,
                               const double finite_difference_increment,
                               const double zeta);
  MultipleShootingContinuation(const double T_f, const double alpha, 
                               const int N,
                               const double initial_time, 
                               const double finite_difference_increment,
                               const double zeta);
  ~MultipleShootingContinuation();

  void integrateSolution(double* control_input_and_constraints_seq, 
                         double** state_mat, double** lambda_mat,
                         const double* control_input_and_constraints_update_seq, 
                         const double integration_length);

  double computeErrorNorm(const double time, const double* state_vec, 
                          const double* control_input_and_constraints_seq,
                          double const* const* state_mat, 
                          double const* const* lambda_mat);

  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  void resetHorizonLength(const double initial_time);

  void bFunc(const double time, const double* state_vec, 
             const double* control_input_and_constraints_seq, 
             double const* const* state_mat, double const* const* lambda_mat,
             const double* current_control_input_and_constraints_update_seq, 
             double* b_vec);

  void AxFunc(const double time, const double* state_vec, 
              const double* control_input_and_constraints_seq, 
              double const* const* state_mat, double const* const* lambda_mat,
              const double* direction_vec, double* ax_vec);

  int dim_state() const;
  int dim_control_input() const;
  int dim_constraints() const;
  int dim_condensed_problem() const;
  int N() const;

  MultipleShootingContinuation(const MultipleShootingContinuation&) = delete;
  MultipleShootingContinuation& operator=(const MultipleShootingContinuation&) 
      = delete;

private:
  MultipleShootingOCP ocp_;
  const int dim_state_, dim_control_input_, dim_constraints_, 
      dim_control_input_and_constraints_, dim_control_input_and_constraints_seq_, 
      N_;
  double finite_difference_increment_, zeta_, incremented_time_; 
  double *incremented_state_vec_, 
      *incremented_control_input_and_constraints_seq_, 
      *control_input_and_constraints_residual_seq_, 
      *control_input_and_constraints_residual_seq_1_, 
      *control_input_and_constraints_residual_seq_2_, 
      *control_input_and_constraints_residual_seq_3_;
  double **incremented_state_mat_, **incremented_lambda_mat_, 
      **state_residual_mat_, **state_residual_mat_1_, 
      **lambda_residual_mat_, **lambda_residual_mat_1_;
};

#endif // MULTIPLE_SHOOTING_CONTINUATION_H