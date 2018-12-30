#include "init_cgmres.hpp"


InitCGMRES::InitCGMRES(const NMPCModel model, const double difference_increment, const int dim_krylov) : MatrixFreeGMRES(model.dimControlInput()+model.dimConstraints(), dim_krylov)
{
    // Set parameters.
    model_ = model;
    difference_increment_ = difference_increment;
    dim_solution_ = model_.dimControlInput()+model_.dimConstraints();

    // Allocate vectors.
    solution_update_vec_.resize(dim_solution_);
    incremented_solution_vec_.resize(dim_solution_);
    lambda_vec_.resize(model_.dimState());
    error_vec_.resize(dim_solution_);
    error_vec_1_.resize(dim_solution_);
    error_vec_2_.resize(dim_solution_);

    // Initialize solution of the forward-difference GMRES.
    for(int i=0; i<dim_solution_; i++){
        solution_update_vec_(i) = 0.0;
        error_vec_(i) = 0.0;
    }
}


void InitCGMRES::solve0stepNOCP(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_vec, const double convergence_radius, const int max_iteration, Eigen::Ref<Eigen::VectorXd> solution_vec)
{
    int i=0;
    solution_vec = initial_guess_vec;
    computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);
    while(error_vec_.squaredNorm() > convergence_radius && i < max_iteration){
        forwardDifferenceGMRES(initial_time, initial_state_vec, solution_vec, solution_update_vec_);
        solution_vec += solution_update_vec_;
        computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);
        i++;
    }
}


Eigen::VectorXd InitCGMRES::getOptimalityErrorVec(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec)
{
    Eigen::VectorXd error_vec(dim_solution_);

    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);

    return error_vec;
}


inline void InitCGMRES::computeOptimalityErrors(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> optimality_vec)
{
    model_.phixFunc(time_param, state_vec, lambda_vec_);
    model_.huFunc(time_param, state_vec, current_solution_vec, lambda_vec_, optimality_vec);
}


void InitCGMRES::bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec) 
{
    computeOptimalityErrors(time_param, state_vec, current_solution_vec, error_vec_);
    equation_error_vec = - error_vec_;
}


void InitCGMRES::axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec)
{
    incremented_solution_vec_ = current_solution_vec + difference_increment_ * direction_vec;
    computeOptimalityErrors(time_param, state_vec, incremented_solution_vec_, error_vec_1_);

    forward_difference_error_vec = (error_vec_1_ - error_vec_) / difference_increment_;
}