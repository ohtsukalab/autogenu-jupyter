#ifndef CONTINUATION_GMRES_H
#define CONTINUATION_GMRES_H

#include "nmpc_model.hpp"
#include "matrixfree_gmres.hpp"
#include "cgmres_initializer.hpp"
#include "single_shooting_continuation.hpp"
#include "linear_algebra.hpp"

class ContinuationGMRES {
public:
  ContinuationGMRES(const double T_f, const double alpha, const int N,
                    const double finite_difference_increment,
                    const double zeta, const int kmax);
  ~ContinuationGMRES();

  void controlUpdate(const double time, const double* state_vec, 
                     const double sampling_period, double* control_input_vec);

  void getControlInput(double* control_input_vec) const;

  void setParametersForInitialization(const double* initial_guess_solution, 
                                      const double newton_residual_tolerance,
                                      const int max_newton_iteration);

  void initializeSolution(const double initial_time,  
                          const double* initial_state_vec);

  double getErrorNorm(const double time, const double* state_vec);

  ContinuationGMRES(const ContinuationGMRES&) = delete;
  ContinuationGMRES& operator=(const ContinuationGMRES&) = delete;

private:
  SingleShootingContinuation continuation_problem_;
  MatrixFreeGMRES<SingleShootingContinuation, const double, const double*, 
                  const double*> mfgmres_;
  CGMRESInitializer solution_initializer_;
  const int dim_control_input_, dim_constraints_;
  double *solution_vec_, *solution_update_vec_, *initial_solution_vec_;
};

#endif // CONTINUATION_GMRES_H