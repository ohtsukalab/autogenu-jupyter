//
// The continuation GMRES (C/GMRES) method, a fast algorithm of nonlinear model predictive control (NMPC).
// This program is witten with reference to "T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)".
//

#ifndef CONTINUATION_GMRES_H
#define CONTINUATION_GMRES_H

#include "linear_funcs.hpp"
#include "matrixfree_gmres.hpp"
#include "nmpc_model.hpp"
#include "init_cgmres.hpp"


// Solves the nonlinear optimal control problem for NMPC using the C/GMRES method.
class ContinuationGMRES final : virtual public MatrixFreeGMRES{
public:
    // Sets parameters and allocates vectors and matrices.
    ContinuationGMRES(const double horizon_max_length, const double alpha, const int horizon_division_num, const double finite_diff_step, const double zeta, const int kmax);

    // Free vectors and matrices.
    ~ContinuationGMRES();

    // Initializes the solution of the C/GMRES method.
    void initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double convergence_radius, const int max_iteration);

    // Updates the solution by solving the matrix-free GMRES.
    void controlUpdate(const double current_time, const double sampling_period, const double* current_state_vec, double* optimal_control_input_vec);

    // Returns the intial vector of the control input sequence
    void getControlInput(double* control_input_vec) const;

    // Returns the optimality error norm under the current_state_vec and the current solution.
    double getError(const double current_time, const double* current_state_vec);


private:
    NMPCModel model_;
    int dim_state_, dim_control_input_, dim_constraints_, dim_control_input_and_constraints_, dim_solution_, horizon_division_num_, kmax_;

    // initial_time_, T_f_, alpha_ : parameters of the length of the horizon
    // The horizon length at time t is given by T_f_*(1.0-std::exp(-alpha_*(t - initial_time_))).
    double initial_time_, T_f_, alpha_, zeta_, finite_diff_step_, incremented_time_;
    double *dx_vec_, *incremented_state_vec_, *solution_vec_, *incremented_solution_vec_, *optimality_vec_, *optimality_vec_1_, *optimality_vec_2_, *solution_update_vec_;
    double **state_mat_, **lambda_mat_;

    // Computes the optimaliy error vector under current_solution_vec.
    inline void computeOptimalityError(const double time_param, const double* state_vec, const double* current_solution_vec, double* optimality_error_vec);

    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec) override;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec) override;
};

#endif