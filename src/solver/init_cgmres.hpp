//
// Computes the initial solution of the C/GMRES method using the Newton GMRES method.
//

#ifndef INIT_CGMRES_H
#define INIT_CGMRES_H

#include "linear_funcs.hpp"
#include "matrixfree_gmres.hpp"


// Computes the initial solution of the C/GMRES method using the Newton GMRES method.
class InitCGMRES final : public MatrixFreeGMRES{
public:
    // Allocates vectors.
    InitCGMRES();

    // Free vectors.
    ~InitCGMRES();

    // Sets parameters and allocates vectors.
    void setInitParams(const double* initial_guess_solution, const double residual_tolerance, const int max_iteration, const double finite_diff_step, const int kmax);

    // Calls the forwardDifferenceGMRES, solves the GMRES, and obtains the solution of the initialization for the C/GMRES method.
    void solveOCPForInit(const double initial_time, const double* initial_state_vec, double* initial_solution_vec, double* optimality_error_vec);


private:
    NMPCModel model_;
    int dim_solution_, max_iteration_;
    double finite_diff_step_, residual_tolerance_;
    double *initial_guess_solution_, *solution_update_vec_, *incremented_solution_vec_, *lambda_vec_, *error_vec_, *error_vec_1_, *error_vec_2_;

    // Computes the optimality error vector under current_solution_vec.
    inline void computeOptimalityErrors(const double time_param, const double* state_vec, const double* current_solution_vec, double* optimality_vec);

    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec) override;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec) override;
};


#endif