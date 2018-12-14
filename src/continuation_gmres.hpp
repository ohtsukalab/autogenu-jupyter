#ifndef CONTINUATION_GMRES_H
#define CONTINUATION_GMRES_H


#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "init_cgmres.hpp"


class ContinuationGMRES final : virtual public MatrixFreeGMRES{
private:
    NMPCModel model_;
    int dim_state_, dim_control_input_, dim_constraints_, dim_1step_solution_, dim_solution_, horizon_division_num_, dim_krylov_;
    double initial_time_, horizon_max_length_, alpha_, zeta_, difference_increment_, incremented_time_;
    Eigen::MatrixXd state_mat_, lambda_mat_;
    Eigen::VectorXd dx_vec_, incremented_state_vec_, solution_vec_, incremented_solution_vec_, optimality_vec_, optimality_vec_1_, optimality_vec_2_, solution_update_vec_;

    inline void computeOptimalityError(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> optimality_vec);
    void nonlinearEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec) override;
    void forwardDifferenceEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec) override;

public:
    ContinuationGMRES(const NMPCModel model, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov);
    void initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const double convergence_radius, const int max_iteration);
    void controlUpdate(const double current_time, const double sampling_period, const Eigen::VectorXd& current_state_vec, Eigen::Ref<Eigen::VectorXd> optimal_control_input_vec);
    double getError(const double current_time, const Eigen::VectorXd& current_state_vec);
};

#endif