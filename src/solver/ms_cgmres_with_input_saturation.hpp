#ifndef MS_CGMRES_WITH_INPUT_SATURATION_H
#define MS_CGMRES_WITH_INPUT_SATURATION_H

#include "matrixfree_gmres.hpp"
#include "ms_continuation_with_input_saturation.hpp"
#include "input_saturation_set.hpp"
#include "ms_cgmres_with_input_saturation_initializer.hpp"
#include "linear_algebra.hpp"

class MSCGMRESWithInputSaturation {
public:
  MSCGMRESWithInputSaturation(const InputSaturationSet& input_saturation_set,
                              const double T_f, const double alpha, const int N,
                              const double finite_difference_increment,
                              const double zeta, const int kmax);
  ~MSCGMRESWithInputSaturation();

  void controlUpdate(const double time, const double* state_vec, 
                     const double sampling_period, double* control_input_vec);

  void getControlInput(double* control_input_vec) const;

  void setParametersForInitialization(const double* initial_guess_solution, 
                                      const double newton_residual_tolerance,
                                      const int max_newton_iteration);

  void initializeSolution(const double initial_time,  
                          const double* initial_state_vec);

  double getErrorNorm(const double time, const double* state_vec);

  MSCGMRESWithInputSaturation(const MSCGMRESWithInputSaturation&) = delete;
  MSCGMRESWithInputSaturation& operator=(const MSCGMRESWithInputSaturation&) 
      = delete;

private:
  MSContinuationWithInputSaturation continuation_problem_;
  MatrixFreeGMRES<MSContinuationWithInputSaturation, const double, const double*, 
                  const double*, double const* const*, double const* const*, 
                  double const* const*, double const* const*> mfgmres_;
  MSCGMRESWithInputSaturationInitializer solution_initializer_;
  const int dim_state_, dim_control_input_, dim_constraints_, dim_saturation_,
            N_;
  double *control_input_and_constraints_seq_, 
         *control_input_and_constraints_update_seq_, 
         *initial_control_input_and_constraints_vec_, *initial_lambda_vec_,
         *initial_dummy_input_vec_, *initial_input_saturation_vec_;
  double **state_mat_, **lambda_mat_, **dummy_input_mat_,
         **input_saturation_multiplier_mat_;
};

#endif // MSCGMRES_WITH_INPUT_SATURATION_H