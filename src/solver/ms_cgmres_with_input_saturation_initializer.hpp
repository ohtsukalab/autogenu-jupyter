// Computes the initial solution of the C/GMRES and multiple shooting-based
// C/GMRES method with condensing of dummy input for the control input 
// saturation and the Lagrange multiplier for the constraints using the Newton 
// GMRES method.

#ifndef MS_CGMRES_WITH_INPUT_SATURATION_INITIALIZER_H
#define MS_CGMRES_WITH_INPUT_SATURATION_INITIALIZER_H

#include <cmath>
#include "newton_gmres_for_ocp.hpp"
#include "zero_horizon_ocp_with_input_saturation.hpp"
#include "input_saturation_set.hpp"
#include "matrixfree_gmres.hpp"
#include "linear_algebra.hpp"

// Solver of the nonlinear optimal control problem for the initialization 
// of the solution in the C/GMERS method. The main method is 
// computeInitialSolution() which computes the solution of the OCP with horizon 
// whose length is zero by Newton GMRES method. Before using 
// computeInitialSolution() method, you have to set initial guess solution by 
// setInitialGuessSolution().
class MSCGMRESWithInputSaturationInitializer {
public:
  // Sets parameters and allocates vectors. 
  // Arguments:
  //  InputSaturationSet: The set composed of the input saturation constraints.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  kmax: A parameter for the GMRES method. This parameter represents the
  //     dimension of the Krylov subspace and maximum iteration number of the
  //     GMRES method.
  //   newton_residual_tolerance: A convergence criteria for the Newton iteration. 
  //     Newton iteration terminates when the error is less than this value.
  //   max_newton_iteration: Maximum number of the Newton iteration. Newton 
  //     iteration for the initialization terminates when the number of the 
  //     iteration is equal to this value.
  MSCGMRESWithInputSaturationInitializer(
      const InputSaturationSet& input_saturation_set,
      const double finite_difference_increment, const int kmax, 
      const double newton_residual_tolerance, const int max_newton_iteration);

  // Sets parameters and allocates vectors. 
  // Arguments:
  //  InputSaturationSet: The set composed of the input saturation constraints.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  kmax: A parameter for the GMRES method. This parameter represents the
  //     dimension of the Krylov subspace and maximum iteration number of the
  //     GMRES method.
  MSCGMRESWithInputSaturationInitializer(
      const InputSaturationSet& input_saturation_set,
      const double finite_difference_increment, const int kmax);

  // Free vectors and matrices.
  ~MSCGMRESWithInputSaturationInitializer();

  // Sets parameters for Newton iteration. 
  //   newton_residual_tolerance: A convergence criteria for the Newton iteration. 
  //     Newton iteration terminates when the error is less than this value.
  //   max_newton_iteration: Maximum number of the Newton iteration. Newton 
  //     iteration for the initialization terminates when the number of the 
  //     iteration is equal to this value.
  void setCriterionsOfNewtonTermination(const double newton_residual_tolerance, 
                                        const int max_newton_iteration);

  // Sets the initial guess soluion, which is composed by the control input
  // vector and the Lgrange multiplier with respect the equality constraints.
  void setInitialGuessSolution(
      const double* initial_guess_control_input_and_constraints);

  // Sets the initial guess of the Lagrange multiplier with respect to the 
  // constraints on the control input saturation function. The all elements of 
  // the multiplier are filled by initial_input_saturation_multiplier.
  void setInitialInputSaturationMultiplier(
      const double initial_input_saturation_multiplier);

  // Sets the initial guess of the Lagrange multiplier with respect to the 
  // constraints on the control input saturation function. The multiplier 
  // is set by initial_input_saturation_multiplier.
  void setInitialInputSaturationMultiplier(
      const double* initial_input_saturation_multiplier);

  // Solves the optimal control problem with horizon whose length is zero
  // under initial_time and initial_state_vec.
  void computeInitialSolution(const double initial_time, 
                              const double* initial_state_vec, 
                              double* initial_control_input_and_constraints_vec,
                              double* initial_dummy_input_vec,
                              double* initial_input_saturation_vec);

  // Computes the initial lambda, which is the Lagrange multiplier of the state 
  // equation. This corresponds to the partial derivative of the terminal cost
  // with respect to the state.
  void getInitialLambda(const double initial_time, 
                        const double* initial_state_vec, 
                        double* initial_lambda_vec);

  // Returns the dimenstion of the solution, which is equivalent to the 
  // dim_control_input and the dim_constraints.
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

  // Computes and returns the dummy input from input, min_input, and max_input.
  double computeDummyInput(const double input, 
                           const double min_input, 
                           const double max_input) const;
};

#endif // MS_CGMRES_WITH_INPUT_SATURATION_INITIALIZER_H