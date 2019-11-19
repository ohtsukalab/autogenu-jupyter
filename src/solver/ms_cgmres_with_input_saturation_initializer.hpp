#ifndef MS_CGMRES_WITH_INPUT_SATURATION_INITIALIZER_H
#define MS_CGMRES_WITH_INPUT_SATURATION_INITIALIZER_H

#include <cmath>
#include "newton_gmres_for_ocp.hpp"
#include "zero_horizon_ocp_with_input_saturation.hpp"
#include "input_saturation_set.hpp"
#include "matrixfree_gmres.hpp"
#include "linear_algebra.hpp"

class MSCGMRESWithInputSaturationInitializer {
public:
  MSCGMRESWithInputSaturationInitializer(
      const InputSaturationSet& input_saturation_set,
      const double finite_difference_increment, const int kmax, 
      const double newton_residual_tolerance, const int max_newton_iteration);

  MSCGMRESWithInputSaturationInitializer(
      const InputSaturationSet& input_saturation_set,
      const double finite_difference_increment, const int kmax);

  ~MSCGMRESWithInputSaturationInitializer();

  void setCriterionsOfNewtonTermination(const double newton_residual_tolerance, 
                                        const int max_newton_iteration);

  void setInitialGuessSolution(
      const double* initial_guess_control_input_and_constraints);

  void setInitialInputSaturationMultiplier(
      const double initial_input_saturation_multiplier);

  void computeInitialSolution(const double initial_time, 
                              const double* initial_state_vec, 
                              double* initial_control_input_and_constraints_vec,
                              double* initial_dummy_input_vec,
                              double* initial_input_saturation_vec);
  
  void getInitialLambda(const double initial_time, 
                        const double* initial_state_vec, 
                        double* initial_lambda_vec);

  int dim_solution() const;

private:
  NewtonGMRESForOCP<ZeroHorizonOCPWithInputSaturation, 
                    InputSaturationSet> newton_;
  MatrixFreeGMRES<NewtonGMRESForOCP<ZeroHorizonOCPWithInputSaturation, 
                                    InputSaturationSet>, 
                  const double, const double*, const double*> mfgmres_;
  InputSaturationSet input_saturation_set_;
  const int dim_control_input_, dim_constraints_, dim_input_saturation_,
            dim_solution_;
  int max_newton_iteration_;
  double newton_residual_tolerance_;
  double *initial_guess_solution_vec_, *initial_solution_vec_, 
      *solution_update_vec_;

  double computeDummyInput(const double input, 
                           const double min_input, 
                           const double max_input) const;
};

#endif // MS_CGMRES_WITH_INPUT_SATURATION_INITIALIZER_H