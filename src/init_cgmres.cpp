#include "init_cgmres.hpp"


inline void InitCGMRES::computeOptimalityErrors(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> optimality_vec)
{
    model_.phixFunc(time_param, state_vec, lambda_vec_);
    model_.huFunc(time_param, state_vec, current_solution_vec, lambda_vec_, optimality_vec);
}


void InitCGMRES::bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> b_vec) 
{
    computeOptimalityErrors(time_param, state_vec, current_solution_vec, error_vec_);
    computeOptimalityErrors(time_param, state_vec, current_solution_vec+difference_increment_*solution_update_vec_, error_vec_1_);

    b_vec = - error_vec_ - (error_vec_1_ - error_vec_)/difference_increment_;
}


void InitCGMRES::axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> ax_vec)
{
    computeOptimalityErrors(time_param, state_vec, current_solution_vec+difference_increment_*direction_vec, error_vec_1_);

    ax_vec = (error_vec_1_ - error_vec_)/difference_increment_;
}


InitCGMRES::InitCGMRES(const double difference_increment, const int max_dim_krylov) : MatrixFreeGMRES(), 
    model_(), 
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()), 
    difference_increment_(difference_increment), 
    solution_update_vec_(Eigen::VectorXd::Zero(dim_solution_)), 
    lambda_vec_(Eigen::VectorXd::Zero(model_.dimState())), 
    error_vec_(Eigen::VectorXd::Zero(dim_solution_)), 
    error_vec_1_(Eigen::VectorXd::Zero(dim_solution_)), 
    error_vec_2_(Eigen::VectorXd::Zero(dim_solution_))
{
    // Set dimensions and parameters in GMRES.
    setGMRESParams(dim_solution_, max_dim_krylov);
}


void InitCGMRES::solve0stepNOCP(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_vec, const double convergence_radius, const int max_iteration, Eigen::Ref<Eigen::VectorXd> solution_vec)
{
    int i=0;
    solution_vec = initial_guess_vec;
    computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);
    while(error_vec_.squaredNorm() > convergence_radius*convergence_radius && i < max_iteration){
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