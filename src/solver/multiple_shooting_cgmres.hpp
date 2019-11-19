#ifndef MULTIPLE_SHOOTING_CGMRES_H
#define MULTIPLE_SHOOTING_CGMRES_H

#include "matrixfree_gmres.hpp"
#include "multiple_shooting_continuation.hpp"
#include "cgmres_initializer.hpp"
#include "linear_algebra.hpp"

class MultipleShootingCGMRES {
public:
  MultipleShootingCGMRES(const double T_f, const double alpha, const int N,
           const double finite_difference_increment,
           const double zeta, const int kmax);
  ~MultipleShootingCGMRES();

  void controlUpdate(const double time, const double* state_vec, 
                     const double sampling_period, double* control_input_vec);

  void getControlInput(double* control_input_vec) const;

  void setParametersForInitialization(const double* initial_guess_solution, 
                                      const double newton_residual_tolerance,
                                      const int max_newton_iteration);

  void initializeSolution(const double initial_time,  
                          const double* initial_state_vec);

  double getErrorNorm(const double time, const double* state_vec);

  MultipleShootingCGMRES(const MultipleShootingCGMRES&) = delete;
  MultipleShootingCGMRES& operator=(const MultipleShootingCGMRES&) = delete;

private:
  MultipleShootingContinuation continuation_problem_;
  MatrixFreeGMRES<MultipleShootingContinuation, const double, const double*, 
                  const double*, double const* const*, double const* const*> 
                  mfgmres_;
  CGMRESInitializer solution_initializer_;
  const int dim_state_, dim_control_input_, dim_constraints_, N_;
  double *control_input_and_constraints_seq_, 
    *control_input_and_constraints_update_seq_, 
    *initial_control_input_and_constraints_vec_, *initial_lambda_vec_;
  double **state_mat_, **lambda_mat_;
};

#endif // MULTIPLE_SHOOTING_CGMRES_H