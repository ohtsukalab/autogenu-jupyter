// This class provides the continuation GMRES (C/GMRES) method, a fast 
// algorithm of nonlinear model predictive control (NMPC). This program is 
// witten with reference to "T. Ohtsuka A continuation/GMRES method for fast 
// computation of nonlinear receding horizon control, Automatica, Vol. 40, 
// No. 4, pp. 563-574 (2004)".

#ifndef CONTINUATION_GMRES_H
#define CONTINUATION_GMRES_H

#include "linear_algebra.hpp"
#include "matrixfree_gmres.hpp"
#include "nmpc_model.hpp"
#include "init_cgmres.hpp"

// Solver of the nonlinear optimal control problem for NMPC using the 
// C/GMRES method, a fast numerical algorithm of NMPC. The main method is
// controlUpdate() that updates the solution of NMPC using the C/GMRES method.
// Before using controlUpdate() method, you have to initialize the solution.
// For this initialization, you are required to set parameters by
// setInitParameters() method and initSolution method. Without these 
// initialization, all components of the solution is zero.
class ContinuationGMRES final : virtual public MatrixFreeGMRES {
public:
  // Constructs ContinuationGMRES with setting parameters and allocates vectors 
  // and matrices used in the C/GMRES method. This constructor also calls 
  // the constructor of GMRES and allocations of GMRES.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  //  finite_difference_step: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  kmax: A parameter for the GMRES method. This parameter represents the
  //     dimension of the Krylov subspace and maximum iteration number of the
  //     GMRES method.
  ContinuationGMRES(const double T_f, const double alpha, const int N, 
                    const double zeta, const double finite_difference_step, 
                    const int kmax);

  // Free vectors and matrices.
  ~ContinuationGMRES();

  // Sets parameters for the initialization of the solution of the C/GMRES 
  // method. Call before initializes the solutino by initSolution() method.
  // This initializaiton is done by solving an optimal control problem (OCP) 
  // with horizon whose length is zero using the Newton-GMRES method.
  // Argments:
  //   initial_guess_solution: An initial guess solution of the OCP vectors
  //     are composed of a contorl input vector and a Lagrange multiplier for
  //     equality constraints.
  //   residual_tolerance: A convergence criteria for the Newton iteration. 
  //     Newton iteration terminates when the error is less than this value.
  //   max_iteration: Maximum number of the Newton iteration. Newton iteration
  //     for the initialization terminates when the number of the iteration 
  //     is equal to this value.
  //   finite_difference_step: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //   kmax: A parameter for the GMRES method. This parameter represents the
  //     dimension of the Krylov subspace and maximum iteration number of the
  //     GMRES method.
  void setInitParameters(const double* initial_guess_solution, 
                         const double residual_tolerance, 
                         const int max_iteration, 
                         const double finite_difference_step, 
                         const int kmax);

  // Initializes the solution of the C/GMRES method by solving the optimal
  // control problem with the horizon whose length is zero. soltuion_vec_ and 
  // errors_in_optimality_ is fullfilled with the solution of this OCP. The 
  // control input to be applied to the actual system is assigned in 
  // optimal_control_input_vec.
  void initSolution(const double initial_time, const double* initial_state_vec, 
                    double* optimal_control_input_vec);

  // Updates the solution by solving the matrix-free GMRES. The optimal control
  // to be applied to the actual system is assigned in optima_control_input_vec.
  void controlUpdate(const double time, const double sampling_period, 
                     const double* state_vec, 
                     double* optimal_control_input_vec);

  // Returns the squared norm of the errors in optimality under the state_vec 
  // and the current solution.
  double getErrorNorm(const double time, const double* state_vec);

private:
  NMPCModel model_;
  InitCGMRES cgmres_initializer_;
  int dim_state_, dim_control_input_, dim_constraints_, 
      dim_control_input_and_constraints_, dim_solution_, N_, kmax_;
  double T_f_, alpha_, zeta_, finite_difference_step_, incremented_time_, 
      initial_time_;
  double *dx_vec_, *incremented_state_vec_, *solution_vec_, 
      *incremented_solution_vec_, *errors_in_optimality_, 
      *errors_in_optimality_1_, *errors_in_optimality_2_, 
      *solution_update_vec_;
  double **state_mat_, **lambda_mat_;

  // Computes the erros in optimaliy under time, state_vec, and 
  // solution_vec. The result is set in errors_in_optimality.
  inline void computeErrorsInOptimality(const double time, 
                                        const double* state_vec, 
                                        const double* solution_vec, 
                                        double* errors_in_optimality);

  // Computes a vector correspongin to b in Ax=b based on the firmulation of
  // the C/GMRES method. 
  void bFunc(const double time, const double* state_vec, 
             const double* current_solution_vec, double* b_vec) override;

  // Generates a vector corresponding to Ax in Ax=b with using the forward 
  // difference approximation.
  void axFunc(const double time, const double* state_vec, 
              const double* current_solution_vec, const double* direction_vec, 
              double* ax_vec) override;
};

#endif // CONTINUATION_GMRES_H