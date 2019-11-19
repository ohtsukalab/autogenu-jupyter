#ifndef SINGLE_SHOOTING_CONTINUATION_H
#define SINGLE_SHOOTING_CONTINUATION_H

#include <cmath>
#include "linear_algebra.hpp"
#include "single_shooting_ocp.hpp"

class SingleShootingContinuation {
public:
  SingleShootingContinuation(const double T_f, const double alpha, const int N,
                             const double finite_difference_increment,
                             const double zeta);
  SingleShootingContinuation(const double T_f, const double alpha, const int N,
                             const double initial_time, 
                             const double finite_difference_increment,
                             const double zeta);
  ~SingleShootingContinuation();

  void integrateSolution(double* solution_vec, 
                         const double* solution_update_vec, 
                         const double integration_length);

  double computeErrorNorm(const double time, const double* state_vec, 
                          const double* solution_vec);

  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  void resetHorizonLength(const double initial_time);

  void bFunc(const double time, const double* state_vec, 
             const double* current_solution_vec, 
             const double* current_solution_update_vec, double* b_vec);

  void AxFunc(const double time, const double* state_vec, 
              const double* current_solution_vec, const double* direction_vec,
              double* ax_vec);

  int dim_state() const;
  int dim_control_input() const;
  int dim_constraints() const;
  int dim_solution() const;
  int N() const;

  // Prohibits copy constructor.
  SingleShootingContinuation(const SingleShootingContinuation&) = delete;
  SingleShootingContinuation& operator=(const SingleShootingContinuation&) 
      = delete;

private:
  SingleShootingOCP ocp_;
  const int dim_state_, dim_control_input_, dim_constraints_, dim_solution_;
  double finite_difference_increment_, zeta_, incremented_time_; 
  double *incremented_state_vec_, *incremented_solution_vec_, 
      *optimality_residual_, *optimality_residual_1_, *optimality_residual_2_;
};

#endif // SINGLE_SHOOTING_CONTINUATION_H 