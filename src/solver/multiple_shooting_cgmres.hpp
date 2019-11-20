// The multiple shooting based continuation GMRES (C/GMRES) method, a fast 
// algorithm of nonlinear model predictive control (NMPC). This program is 
// witten with reference to "T. Ohtsuka A continuation/GMRES method for fast 
// computation of nonlinear receding horizon control, Automatica, Vol. 40, 
// No. 4, pp. 563-574 (2004)" and "Y. Shimizu, T. Ohtsuka, M. Diehl, A 
// real‐time algorithm for nonlinear receding horizon control using multiple 
// shooting and continuation/Krylov method, International Journal of Robust 
// and Nonlinear Control, Vol. 19, No. 8, pp. 919-936 (2008)".

#ifndef MULTIPLE_SHOOTING_CGMRES_H
#define MULTIPLE_SHOOTING_CGMRES_H

#include "matrixfree_gmres.hpp"
#include "multiple_shooting_continuation.hpp"
#include "cgmres_initializer.hpp"
#include "linear_algebra.hpp"

// Solver of the nonlinear optimal control problem for NMPC using the 
// multiple shooting-based C/GMRES method, a fast numerical algorithm of NMPC. 
// This solver also supports condensing of the state and Lagrange multipliers 
// for state equation in the linear problem. The main method is
// controlUpdate() that updates the solution of NMPC using the multiple 
// shooting-based C/GMRES method. 
// Before using controlUpdate() method, you have to initialize the solution.
// For this initialization, you are required to set parameters by
// setParametersForInitialization() method and initializeSolution() method. 
// Without these initialization, all components of the solution is zero.
class MultipleShootingCGMRES {
public:
  // Constructs MultipleShootingCGMRES with setting parameters and allocates 
  // vectors and matrices used in the C/GMRES method. 
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  //  kmax: A parameter for the GMRES method. This parameter represents the
  //     dimension of the Krylov subspace and maximum iteration number of the
  //     GMRES method.
  MultipleShootingCGMRES(const double T_f, const double alpha, const int N,
           const double finite_difference_increment,
           const double zeta, const int kmax);

  // Free vectors and matrices.
  ~MultipleShootingCGMRES();

  // Updates the solution by solving the matrix-free GMRES. The optimal control
  // to be applied to the actual system is assigned in control_input_vec.
  void controlUpdate(const double time, const double* state_vec, 
                     const double sampling_period, double* control_input_vec);

  // Initial value of the current optimal control input is assigned 
  // in control_input_vec.
  void getControlInput(double* control_input_vec) const;

  // Sets parameters for the initialization of the solution of the C/GMRES 
  // method. Call before initializes the solutino by initializeSolution().
  // This initializaiton is done by solving an optimal control problem (OCP) 
  // with horizon whose length is zero using the Newton-GMRES method.
  // Argments:
  //   initial_guess_solution: An initial guess solution of the OCP vectors
  //     are composed of a contorl input vector and a Lagrange multiplier for
  //     equality constraints.
  //   newton_residual_tolerance: A convergence criteria for the Newton iteration. 
  //     Newton iteration terminates when the error is less than this value.
  //   max_newton_iteration: Maximum number of the Newton iteration. Newton 
  //     iteration for the initialization terminates when the number of the 
  //     iteration is equal to this value.
  void setParametersForInitialization(const double* initial_guess_solution, 
                                      const double newton_residual_tolerance,
                                      const int max_newton_iteration);

  // Initializes the solution of the C/GMRES method by solving the optimal
  // control problem with the horizon whose length is zero. soltuion_vec_ and 
  // errors_in_optimality_ is fullfilled with the solution of this OCP. The 
  // control input to be applied to the actual system is assigned in 
  // optimal_control_input_vec.
  void initializeSolution(const double initial_time,  
                          const double* initial_state_vec);

  // Returns the squared norm of the optimality residual under time, state_vec, 
  // and the current solution.
  double getErrorNorm(const double time, const double* state_vec);

  // Prohibits copy due to memory allocation.
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