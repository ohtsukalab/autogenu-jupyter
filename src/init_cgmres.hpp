//
// Computes the initial solution of the C/GMRES method.
//

#ifndef INIT_CGMRES_H
#define INIT_CGMRES_H


#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"


// Computes the initial solution of the C/GMRES method using the Newton method.
class InitCGMRES final : public MatrixFreeGMRES{
private:
    NMPCModel model_;
    int dim_solution_;
    double difference_increment_;
    Eigen::VectorXd solution_update_vec_, incremented_solution_vec_, lambda_vec_, error_vec_, error_vec_1_, error_vec_2_;

    // Computes the optimality error vector under current_solution_vec.
    inline void computeOptimalityErrors(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> optimality_vec);

    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec) override;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec) override;

public:
    // Sets parameters and allocates vectors.
    InitCGMRES(const NMPCModel model, const double difference_increment, const int dim_krylov);

    // Calls the forwardDifferenceGMRES, solves the GMRES, and obtains the solution of the initialization for the C/GMRES method.
    void solve0stepNOCP(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_vec, const double convergence_radius, const int max_iteration, Eigen::Ref<Eigen::VectorXd> solution_vec);

    // Returns the optimality error vector under current_solution_vec.
    Eigen::VectorXd getOptimalityErrorVec(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec);
};


#endif