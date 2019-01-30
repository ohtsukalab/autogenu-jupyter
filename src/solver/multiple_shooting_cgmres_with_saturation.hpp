//
// The multiple shooting based continuation GMRES (C/GMRES) method, a fast algorithm of nonlinear model predictive control (NMPC).
// This program is witten with reference to "T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)" and "Y. Shimizu, T. Ohtsuka, M. Diehl, A real‐time algorithm for nonlinear receding horizon control using multiple shooting and continuation/Krylov method, International Journal of Robust and Nonlinear Control, Vol. 19, No. 8, pp. 919-936 (2008)".
//

#ifndef MULTIPLE_SHOOTING_CGMRES_WITH_SATURATION_H
#define MULTIPLE_SHOOTING_CGMRES_WITH_SATURATION_H


#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "nmpc_model.hpp"
#include "init_cgmres_with_saturation.hpp"
#include "control_input_saturation_sequence.hpp"


// Solves the nonlinear optimal control problem using the mutiple shooting based C/GMRES method with condensing for saturations of the control input.
// Describe the model of a system to be controlled in nmpc_model.hpp and nmpc_model.cpp, and set saturations on the control input in main.cpp. 
class MultipleShootingCGMRESWithSaturation final : virtual public MatrixFreeGMRES{
private:
    NMPCModel model_;
    ControlInputSaturationSequence saturation_seq_;

    int dim_state_, dim_control_input_, dim_constraints_, dim_control_input_and_constraints_, dim_control_input_and_constraints_seq_, dim_saturation_, dim_saturation_seq_, horizon_division_num_, dim_krylov_;

    // initial_time_, horizon_max_length_, alpha_ : parameters of the length of the horizon
    // The horizon length at time t is given by horizon_max_length_*(1.0-std::exp(-alpha_*(time_param-initial_time_))).
    double initial_time_, horizon_max_length_, alpha_, zeta_, difference_increment_, incremented_time_;

    Eigen::VectorXd dx_vec_, incremented_state_vec_, control_input_and_constraints_seq_, control_input_and_constraints_error_seq_, control_input_and_constraints_error_seq_1_, control_input_and_constraints_error_seq_2_, control_input_and_constraints_error_seq_3_, control_input_and_constraints_update_seq_;

    Eigen::MatrixXd state_mat_, lambda_mat_, incremented_state_mat_, incremented_lambda_mat_, state_error_mat_, state_error_mat_1_, lambda_error_mat_, lambda_error_mat_1_, state_update_mat_, lambda_update_mat_, dummy_input_mat_, saturation_lagrange_multiplier_mat_, dummy_error_mat_, dummy_error_mat_1_, saturation_error_mat_, saturation_error_mat_1_, dummy_update_mat_, saturation_update_mat_;


    // Compute 1step optimality error about the saturation
    // Adds partial derivative of the saturation with respect to the control input
    inline void addHamiltonianDerivativeWithControlInput(const Eigen::VectorXd& control_input_and_constraints_vec, const Eigen::VectorXd& saturation_lagrange_multiplier_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_control_input_and_constraints_vec);

    // Computes the optimality for the dummy input
    inline void computeDummyOptimality(const Eigen::VectorXd& dummy_input_vec, const Eigen::VectorXd& saturation_lagrange_multiplier_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_dummy);

    // Computes the optimality of the saturation
    inline void computeSaturationOptimality(const Eigen::VectorXd& control_input_and_constraint_vec, const Eigen::VectorXd& dummy_input_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_saturation);


    // Compute the optimality error sequence
    // Computes the optimaliy error for control input and constraints under current solution.
    inline void computeOptimalityErrorforControlInputAndConstraints(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& state_mat, const Eigen::MatrixXd& lambda_mat, const Eigen::MatrixXd& saturation_lagrange_multiplier_mat, Eigen::Ref<Eigen::VectorXd> optimality_for_control_input_and_constraints);

    // Computes the optimaliy error for state and lambda under current solution.
    inline void computeOptimalityErrorforStateAndLambda(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& state_mat, const Eigen::MatrixXd& lambda_mat, Eigen::Ref<Eigen::MatrixXd> optimality_for_state, Eigen::Ref<Eigen::MatrixXd> optimality_for_lambda);

    // Computes the sequence of state and lambda under the error for state and lambda for condencing.
    inline void computeStateAndLambda(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& optimality_for_state, const Eigen::MatrixXd& optimality_for_lambda, Eigen::Ref<Eigen::MatrixXd> state_mat, Eigen::Ref<Eigen::MatrixXd> lambda_mat);

    // Compute optimality error for saturation on the control input
    inline void computeOptimalityErrorforSaturation(const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& dummy_input_seq, const Eigen::MatrixXd& saturation_lagrange_multiplier_seq, Eigen::Ref<Eigen::MatrixXd> optimality_for_dummy, Eigen::Ref<Eigen::MatrixXd> optimality_for_saturation);


    // Functions for saturations
    inline void multiplySaturationErrorInverse(const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& dummy_input_seq, const Eigen::MatrixXd& saturation_lagrange_multiplier_seq, const Eigen::MatrixXd& multiplied_dummy_input_seq, const Eigen::MatrixXd& multiplied_lagrange_multiplier_seq, Eigen::Ref<Eigen::MatrixXd> resulted_dummy_input_seq, Eigen::Ref<Eigen::MatrixXd> resulted_lagrange_multiplier_seq);

    inline void computeDummyOptimalityDifference(const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& dummy_input_seq, const Eigen::VectorXd& control_input_and_constraints_update_seq, Eigen::Ref<Eigen::MatrixXd> dummy_difference_seq);
    
    inline void computeSaturationOptimalityDifference(const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& dummy_input_seq, const Eigen::MatrixXd& saturation_lagrange_multiplier_seq, const Eigen::VectorXd& control_input_and_constraints_update_seq, Eigen::Ref<Eigen::MatrixXd> saturation_difference_seq);


    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> b_vec) override;

    // Computes a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> ax_vec) override;


public:
    // Sets parameters and allocates vectors and matrices.
    MultipleShootingCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov);


    // Initializes the solution of the multiple shooting based C/GMRES method with condensing of the saturations on the control input
    // 1: with setting all the initial guess Lagrange multipliers of the condensed saturation to 0.
    // 2: with setting the initial guess Lagrange multiplier vector of the condensed saturations to initial_guess_lagrange_multiplier.
    // 3: with setting all the initial guess Lagrange multipliers of the condensed saturation to initial_guess_lagrange_multiplier.
    void initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const double convergence_radius, const int max_iteration);
    void initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const Eigen::VectorXd& initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration);
    void initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const double initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration);


    // Updates the solution by solving the matrix-free GMRES.
    void controlUpdate(const double current_time, const double sampling_period, const Eigen::VectorXd& current_state_vec, Eigen::Ref<Eigen::VectorXd> optimal_control_input_vec);

    // Returns the intial vector of the control input sequence
    Eigen::VectorXd getControlInput() const;

    // Returns the optimality error norm under the current_state_vec and the current solution.
    double getError(const double current_time, const Eigen::VectorXd& current_state_vec);
};


#endif