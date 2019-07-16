// Computes the initial solution of the C/GMRES and multiple shooting-based
// C/GMRES method using the Newton GMRES method.

#ifndef INIT_CGMRES_H
#define INIT_CGMRES_H

#include "linear_algebra.hpp"
#include "matrixfree_gmres.hpp"

// Solver of the nonlinear optimal control problem for the initialization 
// of the solution in the C/GMERS method. The main method is solveOCPForInit()
// which computes the solution of the OCP with horizon whose length is zero
// by Newton GMRES method. Before using solveOCPForInit() method, you have to 
// set parameters by setInitParameters(). 
class InitCGMRES final : public MatrixFreeGMRES {
public:
  // Allocates vectors. Also sets parameters of GMRES with setting
  // kmax by the dimension of the solution.
  InitCGMRES();

  // Free vectors.
  ~InitCGMRES();

  // Sets parameters and allocates vectors. Also sets parameters of GMRES.
  void setInitParameters(const double* initial_guess_solution, 
                         const double residual_tolerance, 
                         const int max_iteration, 
                         const double finite_difference_step, const int kmax);

  // Iterates Newton GMRES by calling solveGMRES() to the OCP whose length 
  // of the horizon is zero, and obtains the solution of the initialization 
  // for the C/GMRES method. The solution is assigned initial_solution_vec
  // and the errors in optimality is assigned initial_error_in_optimality.
  void solveOCPForInit(const double initial_time, 
                       const double* initial_state_vec, 
                       double* initial_solution_vec, 
                       double* initial_errors_in_optimality);

private:
  NMPCModel model_;
  int dim_solution_, max_iteration_;
  double finite_difference_step_, residual_tolerance_;
  double *initial_guess_solution_, *solution_update_vec_, 
      *incremented_solution_vec_, *lambda_vec_, *errors_in_optimality_, 
      *errors_in_optimality_1_, *errors_in_optimality_2_;

  // Computes the error in optimality under time, state_vec, and solution_vec. 
  // The errors in optimality is assigned in errors_in_optimality.
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

#endif // INIT_CGMRES_H