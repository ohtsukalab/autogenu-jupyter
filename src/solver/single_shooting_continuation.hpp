// This class provides the linear problem of the continuation transformation 
// for the single-shooting optimal control problem, which is solved in 
// Matrix-free GMRES. This program is witten with reference to "T. Ohtsuka A 
// continuation/GMRES method for fast computation of nonlinear receding horizon 
// control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)".

#ifndef SINGLE_SHOOTING_CONTINUATION_H
#define SINGLE_SHOOTING_CONTINUATION_H

#include <cmath>
#include "linear_algebra.hpp"
#include "single_shooting_ocp.hpp"

// Linear problem of the continuation transformation for the single-shooting 
// optimal control problem, which is solved in Matrix-free GMRES. This class 
// is intended for use with MatrixfreeGMRES class. 
class SingleShootingContinuation {
public:
  // Constructs SingleShootingContinuation with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  SingleShootingContinuation(const double T_f, const double alpha, const int N,
                             const double finite_difference_increment,
                             const double zeta);

  // Constructs SingleShootingContinuation with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  SingleShootingContinuation(const double T_f, const double alpha, const int N,
                             const double initial_time, 
                             const double finite_difference_increment,
                             const double zeta);

  // Free vectors and matrices.
  ~SingleShootingContinuation();

  // Integrates the solution for given optimal update vector of the solution 
  // and the integration length.
  void integrateSolution(double* solution_vec, 
                         const double* solution_update_vec, 
                         const double integration_length);

  // Computes and returns the squared norm of the errors in optimality under 
  // the state_vec and the current solution.
  double computeErrorNorm(const double time, const double* state_vec, 
                          const double* solution_vec);

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
             const double* current_solution_vec, 
             const double* current_solution_update_vec, double* b_vec);

  // Computes a vector correspongin to Ax in Ax=b. This function is called in
  // MatrixfreeGMRES.
  void AxFunc(const double time, const double* state_vec, 
              const double* current_solution_vec, const double* direction_vec,
              double* ax_vec);

  // Returns the dimension of the state.
  int dim_state() const;

  // Returns the dimension of the control input.
  int dim_control_input() const;

  // Returns the dimension of the equality constraints.
  int dim_constraints() const;

  // Returns the dimension of the solution, which is equivalent to 
  // N*(dim_control_input+dim_constraints).
  int dim_solution() const;

  // Returns the grid number of the horizon.
  int N() const;

private:
  SingleShootingOCP ocp_;
  const int dim_state_, dim_control_input_, dim_constraints_, dim_solution_;
  double finite_difference_increment_, zeta_, incremented_time_; 
  double *incremented_state_vec_, *incremented_solution_vec_, 
      *optimality_residual_, *optimality_residual_1_, *optimality_residual_2_;
};

#endif // SINGLE_SHOOTING_CONTINUATION_H 