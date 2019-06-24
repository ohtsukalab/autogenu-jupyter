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
    // Sets parameters and allocates vectors.
    InitCGMRES(const double finite_diff_step, const int kmax);

    // Free vectors.
    ~InitCGMRES();

    // Calls the forwardDifferenceGMRES, solves the GMRES, and obtains the solution of the initialization for the C/GMRES method.
    void solve0stepNOCP(const double initial_time, const double* initial_state_vec, const double* initial_guess_vec, const double convergence_radius, const int max_iteration, double* solution_vec);

    // Returns the optimality error vector under current_solution_vec.
    void getOptimalityErrorVec(const double initial_time, const double* initial_state_vec, const double* current_solution_vec, double* error_vec);


private:
    NMPCModel model_;
    int dim_solution_;
    double finite_diff_step_;
    double *solution_update_vec_, *incremented_solution_vec_, *lambda_vec_, *error_vec_, *error_vec_1_, *error_vec_2_;

    // Computes the optimality error vector under current_solution_vec.
    inline void computeOptimalityErrors(const double time_param, const double* state_vec, const double* current_solution_vec, double* optimality_vec);

    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec) override;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec) override;
};


#endif