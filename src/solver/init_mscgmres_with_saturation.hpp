// Computes the initial solution of the multiple shooting-based C/GMRES method 
// with condensing of state, lambda, and control input saturation using the 
// Newton GMRES method.

#ifndef INIT_MSCGMRES_WITH_SATURATION_H
#define INIT_MSCGMRES_WITH_SATURATION_H

#include "linear_algebra.hpp"
#include "matrixfree_gmres.hpp"
#include "control_input_saturation_sequence.hpp"
#include "condensing_functions.hpp"

// Solver of the nonlinear optimal control problem for the initialization 
// of the solution in the multiple shooting-based C/GMERS method. The main 
// method is solveOCPForInit() which computes the solution of the OCP with 
// horizon whose length is zero by Newton GMRES method. Before using 
// solveOCPForInit() method, you have to set parameters by setInitParameters().
class InitMSCGMRESWithSaturation final : public MatrixFreeGMRES {
public:
  // Allocates vectors. Also sets parameters of GMRES with setting
  // kmax by the dimension of the solution.
  InitMSCGMRESWithSaturation(
      const ControlInputSaturationSequence saturation_seq);

  // Free vectors/
  ~InitMSCGMRESWithSaturation();

  // Sets parameters and allocates vectors. Also sets parameters of GMRES.
  void setInitParameters(const double* initial_guess_solution, 
                         const double* initial_guess_Lagrange_multiplier,
                         const double residual_tolerance, 
                         const int max_iteration, 
                         const double finite_difference_step, const int kmax);

  // Iterates Newton GMRES by calling solveGMRES() to the OCP whose length 
  // of the horizon is zero, and obtains the solution of the initialization 
  // for the C/GMRES method. The solution is assigned initial_solution_vec.
  // The errors in optimality for control input and constraints that is not 
  // condensed is assigned in errors_for_control_input_and_constraints, errors 
  // in optimality for dummy input is assigned in errors_for_dummy_input, and 
  // errors in optimality for constraints of input saturation is assigned in
  // errors_for_control_input_saturations.
  void solveOCPForInit(const double initial_time, 
                       const double* initial_state_vec, 
                       double* initial_solution_vec, 
                       double* errors_for_control_input_and_constraints, 
                       double* errors_for_dummy_input, 
                       double* errors_for_control_input_sauturations);

private:
  NMPCModel model_;
  ControlInputSaturationSequence saturation_seq_;

  int dim_control_input_and_constraints_, dim_saturation_, dim_solution_, 
      max_iteration_;
  double finite_difference_step_, residual_tolerance_;
  double *initial_guess_solution_, *initial_guess_Lagrange_multiplier_, 
      *incremented_solution_vec_, *solution_update_vec_, *lambda_vec_, 
      *errors_in_optimality_, *errors_in_optimality_1_, 
      *errors_in_optimality_2_;

  // Computes the optimality error vector under current_solution_vec.
  inline void computeErrorsInOptimality(const double time, 
                                        const double* state_vec, 
                                        const double* solution_vec, 
                                        double* errors_in_optimality);

  // Computes a vector correspongin to b in Ax=b for Newton GMRES method.
  void bFunc(const double time, const double* state_vec, 
             const double* current_solution_vec, double* b_vec) override;

  // Generates a vector corresponding to Ax in Ax=b with using the forward 
  // difference approximation for Newton GMRES method.
  void axFunc(const double time, const double* state_vec, 
              const double* current_solution_vec, const double* direction_vec,
              double* ax_vec) override;
};

#endif // INIT_MSCGMRES_WITH_SATURATION_H