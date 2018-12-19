#ifndef MULTIPLE_SHOOTING_CGMRES_H
#define MULTIPLE_SHOOTING_CGMRES_H


#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "init_cgmres.hpp"

class MultipleShootingCGMRES : virtual public MatrixFreeGMRES{
private:
    NMPCModel model_;
    int dim_state_, dim_control_input_, dim_constraints_, dim_control_input_and_constraints_, dim_state_and_lambda_, dim_control_input_and_constraints_seq_, dim_state_and_lambda_seq_, horizon_division_num_, dim_krylov_;

    double initial_time_, horizon_max_length_, alpha_, zeta_, difference_increment_, incremented_time_;
    Eigen::VectorXd dx_vec_, incremented_state_vec_, control_input_and_constraints_seq_, incremented_control_input_and_constraints_seq_, control_input_and_constraints_error_seq_, control_input_and_constraints_error_seq_1_, control_input_and_constraints_error_seq_2_, control_input_and_constraints_error_seq_3_,control_input_and_constraints_update_seq_;
    Eigen::MatrixXd state_mat_, lambda_mat_, incremented_state_mat_, incremented_lambda_mat_, state_error_mat_, state_error_mat_1_, lambda_error_mat_, lambda_error_mat_1_, state_update_mat_, lambda_update_mat_;

    inline void computeOptimalityErrorforControlInputAndConstraints(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& state_mat, const Eigen::MatrixXd& lambda_mat, Eigen::Ref<Eigen::VectorXd> optimality_for_control_input_and_constraints);
    inline void computeOptimalityErrorforStateAndLambda(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& state_mat, const Eigen::MatrixXd& lambda_mat, Eigen::Ref<Eigen::MatrixXd> optimality_for_state, Eigen::Ref<Eigen::MatrixXd> optimality_for_lambda);
    inline void computeStateAndLambda(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& optimality_for_state, const Eigen::MatrixXd& optimality_for_lambda, Eigen::Ref<Eigen::MatrixXd> state_mat, Eigen::Ref<Eigen::MatrixXd> lambda_mat);
    void nonlinearEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec) override;
    void forwardDifferenceEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec) override;

public:
    MultipleShootingCGMRES(const NMPCModel model, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov);
    void initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const double convergence_radius, const int max_iteration);
    void controlUpdate(const double current_time, const double sampling_period, const Eigen::VectorXd& current_state_vec, Eigen::Ref<Eigen::VectorXd> optimal_control_input_vec);
    double getError(const double current_time, const Eigen::VectorXd& current_state_vec);
};



#endif