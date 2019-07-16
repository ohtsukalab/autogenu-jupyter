// The multiple shooting based continuation GMRES (C/GMRES) method, a fast 
// algorithm of nonlinear model predictive control (NMPC). This program is 
// witten with reference to "T. Ohtsuka A continuation/GMRES method for fast 
// computation of nonlinear receding horizon control, Automatica, Vol. 40, 
// No. 4, pp. 563-574 (2004)" and "Y. Shimizu, T. Ohtsuka, M. Diehl, A 
// real‐time algorithm for nonlinear receding horizon control using multiple 
// shooting and continuation/Krylov method, International Journal of Robust 
// and Nonlinear Control, Vol. 19, No. 8, pp. 919-936 (2008)".

#ifndef MSCGMRES_WITH_SATURATION_H
#define MSCGMRES_WITH_SATURATION_H

#include "linear_algebra.hpp"
#include "matrixfree_gmres.hpp"
#include "nmpc_model.hpp"
#include "condensing_functions.hpp"
#include "control_input_saturation_sequence.hpp"
#include "init_mscgmres_with_saturation.hpp"

// Solver of the nonlinear optimal control problem for NMPC using the 
// multiple shooting-based C/GMRES method, a fast numerical algorithm of NMPC. 
// This solver also supports condensing of the state, Lagrange multipliers 
// for state equation, and constraints on the saturation function of the
// control input in the linear problem. The main method is
// controlUpdate() that updates the solution of NMPC using the multiple 
// shooting-based C/GMRES method. Before using controlUpdate() method, you 
// have to initialize the solution. For this initialization, you are required 
// to set parameters by setInitParameters() method and initSolution method. 
// Without these initialization, all components of the solution is zero.
class MSCGMRESWithSaturation final : virtual public MatrixFreeGMRES {
public:
  // Constructs MSCGMRES with setting parameters and allocates vectors and 
  // matrices used in the multiple shooting-based C/GMRES method. This 
  // constructor also calls the constructor of GMRES and allocations of GMRES.
  // Arguments:
  //  saturatio_seq: Parameters about the control input saturation.
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
  MSCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, 
                         const double T_f, const double alpha, const int N, 
                         const double zeta, const double finite_difference_step,
                         const int kmax);

  // Free vectors and matrices.
  ~MSCGMRESWithSaturation();

  // Sets parameters for the initialization of the solution of the multiple
  // shooting-based C/GMRES method.  Call before initializes the solutino by 
  // initSolution() method. This initializaiton is done by solving an optimal 
  // control problem (OCP) with horizon whose length is zero using the 
  // Newton-GMRES method.
  // Argments:
  //   initial_guess_solution: An initial guess solution of the OCP vectors
  //     are composed of a contorl input vector and a Lagrange multiplier for
  //     equality constraints.
  //   initial_guess_Lagrange_multipliers: An initial guess of Lagrange 
  //     multipliers for condensed constraints on the saturation function
  //     of the control input.
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
                         const double* initial_guess_Lagrange_multipliers,
                         const double residual_tolerance, 
                         const int max_iteration, 
                         const double finite_difference_step, const int kmax);

  // Initializes the solution of the multiple shooting-based C/GMRES method by 
  // solving the optimal control problem with the horizon whose length is 
  // zero. soltuion_vec_ and errors_in_optimality_ is fullfilled with the 
  // solution of this OCP. The control input to be applied to the actual 
  // system is assigned in optimal_control_input_vec.
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
  ControlInputSaturationSequence saturation_seq_;
  InitMSCGMRESWithSaturation mscgmres_initializer_;
  int dim_state_, dim_control_input_, dim_constraints_, 
      dim_control_input_and_constraints_, 
      dim_control_input_and_constraints_seq_, dim_saturation_, 
      dim_saturation_seq_, N_, dim_krylov_;
  int *saturation_index;
  double T_f_, alpha_, zeta_, finite_difference_step_, incremented_time_, 
      initial_time_;
  double *dx_vec_, *incremented_state_vec_, 
      *control_input_and_constraints_seq_, 
      *incremented_control_input_and_constraints_seq_, 
      *control_input_and_constraints_error_seq_, 
      *control_input_and_constraints_error_seq_1_, 
      *control_input_and_constraints_error_seq_2_, 
      *control_input_and_constraints_error_seq_3_, 
      *control_input_and_constraints_update_seq_;
  double **state_mat_, **state_mat_1_, **lambda_mat_1_, **lambda_mat_, 
      **incremented_state_mat_, **incremented_lambda_mat_, **state_error_mat_, 
      **state_error_mat_1_, **lambda_error_mat_, **lambda_error_mat_1_, 
      **dummy_input_mat_, **dummy_input_mat_1_, 
      **saturation_Lagrange_multiplier_mat_, 
      **saturation_Lagrange_multiplier_mat_1_, 
      **incremented_saturation_Lagrange_multiplier_mat_, **dummy_error_mat_, 
      **dummy_error_mat_1_,**saturation_error_mat_, **saturation_error_mat_1_, 
      **dummy_update_mat_, **saturation_update_mat_;

  // Computes errors in optimality for control input and constraints under 
  // time, state_vec, control_input_and_constraints_seq, state_mat, lambda_mat,
  // and saturation_Lagrange_multiplier_mat. The resulted errors are assigned 
  // in errors_for_control_input_and_constraints.
  inline void computeErrorsForControlInputAndConstraints(
      const double time, const double* state_vec, 
      const double* control_input_and_constraints_seq, 
      double const* const* state_mat, double const* const* lambda_mat, 
      double const* const* saturation_Lagrange_multiplier_mat,
      double* errors_for_control_input_and_constraints);

  // Computes errors in optimality for state and the Lagrange multiplier
  // for the state equation under time, state_vec, 
  // control_input_and_constraints_seq, state_mat, and lambda_mat. 
  // The resulted errors are assigned in errors_for_state and errors_for_lambda.
  inline void computeErrorsForStateAndLambda(
      const double time, const double* state_vec, 
      const double* control_input_and_constraints_seq, 
      double const* const* state_mat, double const* const* lambda_mat, 
      double** errors_for_state, double** errors_for_lambda);

  // Computes the sequence of state and the Lagrange multiplier for the state 
  // equation under time, state_vec, control_input_and_constraints_seq, 
  // errors_for_state, and errors_for_lambda. The results are assigned in 
  // state_mat and lambda_mat.
  inline void computeStateAndLambdaFromErrors(
      const double time, const double* state_vec, 
      const double* control_input_and_constraints_seq, 
      double const* const* errors_for_state, 
      double const* const* errors_for_lambda, double** state_mat, 
      double** lambda_mat);

  // Computes errors in optimality for dummy input and constraints
  // on the saturation functions for the contorl input. 
  // The resulted errors are assigned in errors_for_dummy_input and
  // errors_for_input_saturation.
  inline void computeErrorsForDummyInputAndSaturation(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* saturation_Lagrange_multiplier_mat, 
      double** errors_for_dummy_input, 
      double** errors_for_input_saturation);

  // Computes the invers of the matrix of errors in optimality for dummy input
  // and constraints on the saturation functions for the control input,
  // multiply that inverse matrix to dummy input matrix and Lagrange multiplier
  // mat, and returns that result. The multiplied matrices have to be assigned
  // in multiplied_dummy_input_mat and multiplied_Lagrange_multiplier_mat and 
  // the resulted matrices are assigned in resulted_dummy_input_mat and in 
  // resulted_Lagrange_multiplier_mat.
  inline void multiplyErrorsForInputAndSaturationInverse(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* saturation_Lagrange_multiplier_mat, 
      double const* const* multiplied_dummy_input_mat, 
      double const* const* multiplied_Lagrange_multiplier_mat, 
      double** resulted_dummy_input_mat, 
      double** resulted_Lagrange_multiplier_mat);

  // Computes the difference value of the dummy input corresponding to the
  // difference in the control input and Lagrange multiplier with respect to
  // the constraints that are not condensed. That update quantities have to 
  // be assigned in control_input_and_constraints_update_mat and the result is
  // assigned in dummy_difference_mat.
  inline void computeErrorsDifferenceForDummyInput(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      const double* control_input_and_constraints_update_seq, 
      double** dummy_difference_mat);

  // Computes the difference value of the Lagrange multipliers with respect to
  // the constraints on the condensed saturation functions of the control
  // input corresponding to the difference in the control input and Lagrange 
  // multiplier with respect to the constraints that are not condensed. That 
  // update quantities have to be assigned in 
  // control_input_and_constraints_update_mat and the result is assigned in
  // saturation_difference_mat.
  inline void computeErrorsDifferenceForSaturation(
      const double* control_input_and_constraints_seq, 
      double const* const* dummy_input_mat, 
      double const* const* saturation_Lagrange_multiplier_mat, 
      const double* control_input_and_constraints_update_seq, 
      double** saturation_difference_mat);

  // Computes a vector correspongin to b in Ax=b based on the firmulation of
  // the multiple shooting-based C/GMRES method with condensing of state and 
  // Lagrange multiplier for the state equation.
  void bFunc(const double time, const double* state_vec, 
             const double* current_solution_vec, double* b_vec) override;

  // Generates a vector corresponding to Ax in Ax=b with using the forward 
  // difference approximation.
  void axFunc(const double time, const double* state_vec, 
              const double* current_solution_vec, const double* direction_vec,
              double* ax_vec) override;
};

#endif // MSCGMRES_WITH_SATURATION_H