//
// The multiple shooting based continuation GMRES (C/GMRES) method, a fast algorithm of nonlinear model predictive control (NMPC).
// This program is witten with reference to "T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)" and "Y. Shimizu, T. Ohtsuka, M. Diehl, A real‐time algorithm for nonlinear receding horizon control using multiple shooting and continuation/Krylov method, International Journal of Robust and Nonlinear Control, Vol. 19, No. 8, pp. 919-936 (2008)".
//

#ifndef MULTIPLE_SHOOTING_CGMRES_WITH_SATURATION_H
#define MULTIPLE_SHOOTING_CGMRES_WITH_SATURATION_H


#include "linear_funcs.hpp"
#include "matrixfree_gmres.hpp"
#include "nmpc_model.hpp"
#include "init_cgmres_with_saturation.hpp"
#include "control_input_saturation_sequence.hpp"


// Solves the nonlinear optimal control problem using the mutiple shooting based C/GMRES method with condensing for saturations of the control input.
class MultipleShootingCGMRESWithSaturation final : virtual public MatrixFreeGMRES{
private:
    NMPCModel model_;
    ControlInputSaturationSequence saturation_seq_;

    int dim_state_, dim_control_input_, dim_constraints_, dim_control_input_and_constraints_, dim_control_input_and_constraints_seq_, dim_saturation_, dim_saturation_seq_, horizon_division_num_, dim_krylov_;
    int *saturation_index;
    // initial_time_, horizon_max_length_, alpha_ : parameters of the length of the horizon
    // The horizon length at time t is given by horizon_max_length_*(1.0-std::exp(-alpha_*(t - initial_time_))).
    double initial_time_, horizon_max_length_, alpha_, zeta_, difference_increment_, incremented_time_;

    double *dx_vec_, *incremented_state_vec_, *control_input_and_constraints_seq_, *incremented_control_input_and_constraints_seq_, *control_input_and_constraints_error_seq_, *control_input_and_constraints_error_seq_1_, *control_input_and_constraints_error_seq_2_, *control_input_and_constraints_error_seq_3_, *control_input_and_constraints_update_seq_;

    double **state_mat_, **state_mat_1_, **lambda_mat_1_, **lambda_mat_, **incremented_state_mat_, **incremented_lambda_mat_, **state_error_mat_, **state_error_mat_1_, **lambda_error_mat_, **lambda_error_mat_1_,**dummy_input_mat_, **dummy_input_mat_1_, **saturation_lagrange_multiplier_mat_, 
    **saturation_lagrange_multiplier_mat_1_, **incremented_saturation_lagrange_multiplier_mat_, **dummy_error_mat_, **dummy_error_mat_1_,**saturation_error_mat_, **saturation_error_mat_1_, **dummy_update_mat_, **saturation_update_mat_;


    // Adds partial derivative of the saturation with respect to the optimality error for the control input and other constraints
    inline void addHamiltonianDerivativeWithControlInput(const double* control_input_and_constraints_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_control_input_and_constraints_vec);

    // Computes the optimality error for the dummy input
    inline void computeDummyOptimality(const double* dummy_input_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_dummy);

    // Computes the optimality error of the saturation
    inline void computeSaturationOptimality(const double* control_input_and_constraint_vec, const double* dummy_input_vec, double* optimality_for_saturation);


    // Computes the optimaliy error sequence for the control input and constraints sequence under the current solution.
    inline void computeOptimalityErrorforControlInputAndConstraints(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* state_mat, double const* const* lambda_mat, double const* const* saturation_lagrange_multiplier_mat, double* optimality_for_control_input_and_constraints);

    // Computes the optimaliy error sequence for the state and lambda under the current solution.
    inline void computeOptimalityErrorforStateAndLambda(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* state_mat, double const* const* lambda_mat, double** optimality_for_state, double** optimality_for_lambda);

    // Computes the sequence of state and lambda under the optimality error sequence for state and lambda for condencing.
    inline void computeStateAndLambda(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* optimality_for_state, double const* const* optimality_for_lambda, double** state_mat, double** lambda_mat);

    // Computes the optimality error sequence for the saturation on the function of the control input
    inline void computeOptimalityErrorforSaturation(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, double const* const* saturation_lagrange_multiplier_seq, double** optimality_for_dummy, double** optimality_for_saturation);


    // A function for condensing the variables with respect to the saturation function on the control input
    inline void multiplySaturationErrorInverse(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, double const* const* saturation_lagrange_multiplier_seq, double const* const* multiplied_dummy_input_seq, double const* const* multiplied_lagrange_multiplier_seq, double** resulted_dummy_input_seq, double** resulted_lagrange_multiplier_seq);

    // A functions for condensing the variables with respect to the saturation function on the control input
    inline void computeDummyOptimalityDifference(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, const double* control_input_and_constraints_update_seq, double** dummy_difference_seq);
    
    // A functions for condensing the variables with respect to the saturation function on the control input
    inline void computeSaturationOptimalityDifference(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, double const* const* saturation_lagrange_multiplier_seq, const double* control_input_and_constraints_update_seq, double** saturation_difference_seq);


    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec) override;

    // Computes a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec) override;


public:
    // Sets parameters and allocates vectors and matrices.
    MultipleShootingCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov);

    // Free vectors and matrices.
    ~MultipleShootingCGMRESWithSaturation();


    // Initializes the solution of the multiple shooting based C/GMRES method with condensing of the saturations on the control input
    // with setting all the initial guess Lagrange multipliers of the condensed saturation to 0.
    void initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double convergence_radius, const int max_iteration);

    // with setting all the initial guess Lagrange multipliers of the condensed saturation to initial_guess_lagrange_multiplier.
    void initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration);

    // with setting the initial guess Lagrange multiplier vector of the condensed saturations to initial_guess_lagrange_multiplier.
    void initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double* initial_guess_lagrange_multiplier_vec, const double convergence_radius, const int max_iteration);


    // Updates the solution by solving the matrix-free GMRES.
    void controlUpdate(const double current_time, const double sampling_period, const double* current_state_vec, double* optimal_control_input_vec);

    // Set the control input vector as the intial vector of the current optimal control input sequence
    void getControlInput(double* control_input_vec) const;

    // Returns the optimality error norm under the current_state_vec and the current solution
    double getError(const double current_time, const double* current_state_vec);
};


#endif