#ifndef CGMRES_INITIALIZER_H
#define CGMRES_INITIALIZER_H

#include "newton_gmres_for_ocp.hpp"
#include "zero_horizon_ocp.hpp"
#include "matrixfree_gmres.hpp"
#include "linear_algebra.hpp"

class CGMRESInitializer {
public:
  CGMRESInitializer(const double finite_difference_increment,
                    const int kmax, const double newton_residual_tolerance, 
                    const int max_newton_iteration);
  CGMRESInitializer(const double finite_difference_increment, const int kmax);
  ~CGMRESInitializer();

  void setCriterionsOfNewtonTermination(const double newton_residual_tolerance, 
                                        const int max_newton_iteration);

  void setInitialGuessSolution(const double* initial_guess_solution);

  void computeInitialSolution(const double initial_time, 
                              const double* initial_state_vec, 
                              double* initial_solution_vec);
  
  void getInitialLambda(const double initial_time, 
                        const double* initial_state_vec, 
                        double* initial_lambda_vec);

  int dim_solution() const;

private:
  NewtonGMRESForOCP<ZeroHorizonOCP> newton_;
  MatrixFreeGMRES<NewtonGMRESForOCP<ZeroHorizonOCP>, 
                  const double, const double*, const double*> mfgmres_;
  const int dim_control_input_, dim_constraints_, dim_solution_;
  int max_newton_iteration_;
  double newton_residual_tolerance_;
  double *initial_guess_solution_vec_, *solution_update_vec_;
};

#endif // CGMRES_INITIALIZER_H