//
// The continuation GMRES (C/GMRES) method, a fast algorithm of nonlinear model predictive control (NMPC).
// This program is witten with reference to "T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)".
//

#ifndef CONTINUATION_GMRES_H
#define CONTINUATION_GMRES_H

#include <pick_model.hpp>
#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "init_cgmres.hpp"


// Solves the nonlinear optimal control problem using the C/GMRES method.
// You have to describe the model of a system you want to control in nmpc_model.hpp and nmpc_model.cpp. 
class ContinuationGMRES //final : virtual public MatrixFreeGMRES{
{
private:
    NMPCModel model_;
    int dim_state_, dim_control_input_, dim_constraints_, dim_control_input_and_constraints_, dim_solution_, horizon_division_num_, dim_krylov_;

    // initial_time_, horizon_max_length_, alpha_ : parameters of the length of the horizon
    // The horizon length at time t is given by horizon_max_length_*(1.0-std::exp(-alpha_*(time_param-initial_time_))).
    double initial_time_, horizon_max_length_, alpha_, zeta_, difference_increment_, incremented_time_;
    Eigen::MatrixXd state_mat_, lambda_mat_;
    Eigen::VectorXd dx_vec_, incremented_state_vec_, solution_vec_, optimality_vec_, optimality_vec_1_, optimality_vec_2_, solution_update_vec_;

    // Computes the optimaliy error vector under current_solution_vec.
    inline void computeOptimalityError(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> optimality_vec);

    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec);

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec);

    /* Copied from matrix-free GMRES */
        int dim_equation_, max_dim_krylov_;
        Eigen::MatrixXd hessenberg_mat_, basis_mat_;
        Eigen::VectorXd b_vec_, givens_c_vec_, givens_s_vec_, g_vec_;

        // Applies the Givens rotation for i_column element and i_column+1 element of  column_vec, which is a column vector of a matrix.
        inline void givensRotation(Eigen::Ref<Eigen::VectorXd> column_vec, const int i_column);

        // Solves the GMRES and generates solution_update_vector, which is a solution of the matrix-free GMRES.
        void forwardDifferenceGMRES(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> solution_update_vec);

public:
    // Sets parameters and allocates vectors and matrices.
    // ContinuationGMRES(const NMPCModel model, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov);
    void setSolver(const NMPCModel model, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov);

    // Initializes the solution of the C/GMRES method.
    void initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const double convergence_radius, const int max_iteration);

    // Updates the solution by solving the matrix-free GMRES.
    void controlUpdate(const double current_time, const double sampling_period, const Eigen::VectorXd& current_state_vec, Eigen::Ref<Eigen::VectorXd> optimal_control_input_vec);

    // Returns the optimality error norm under the current_state_vec and the current solution.
    double getError(const double current_time, const Eigen::VectorXd& current_state_vec);
};


#endif