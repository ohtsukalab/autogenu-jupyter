#include "init_cgmres.hpp"


InitCGMRES::InitCGMRES() : MatrixFreeGMRES(), 
    model_(), 
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()), 
    max_iteration_(0),
    finite_diff_step_(0), 
    residual_tolerance_(0),
    initial_guess_solution_(linearfunc::newVector(dim_solution_)),
    solution_update_vec_(linearfunc::newVector(dim_solution_)), 
    incremented_solution_vec_(linearfunc::newVector(dim_solution_)), 
    lambda_vec_(linearfunc::newVector(model_.dimState())), 
    error_vec_(linearfunc::newVector(dim_solution_)), 
    error_vec_1_(linearfunc::newVector(dim_solution_)), 
    error_vec_2_(linearfunc::newVector(dim_solution_))
{
    setGMRESParams(dim_solution_, dim_solution_);
}


InitCGMRES::~InitCGMRES()
{
    linearfunc::deleteVector(initial_guess_solution_);
    linearfunc::deleteVector(solution_update_vec_);
    linearfunc::deleteVector(incremented_solution_vec_);
    linearfunc::deleteVector(lambda_vec_);
    linearfunc::deleteVector(error_vec_);
    linearfunc::deleteVector(error_vec_1_);
    linearfunc::deleteVector(error_vec_2_);
}


void InitCGMRES::setInitParams(const double* initial_guess_solution, const double residual_tolerance, const int max_iteration, const double finite_diff_step, const int kmax)
{
    for (int i=0; i<dim_solution_; i++) {
        initial_guess_solution_[i] = initial_guess_solution[i];
    }
    residual_tolerance_ = residual_tolerance;
    max_iteration_ = max_iteration;
    finite_diff_step_ = finite_diff_step;
    setGMRESParams(dim_solution_, kmax);
}


void InitCGMRES::solveOCPForInit(const double initial_time, const double* initial_state_vec, double* initial_solution_vec, double* optimality_error_vec)
{
    for (int i=0; i<dim_solution_; i++) {
        initial_solution_vec[i] = initial_guess_solution_[i];
    }
    computeOptimalityErrors(initial_time, initial_state_vec, initial_solution_vec, optimality_error_vec);
    int j = 0;
    while (linearfunc::squaredNorm(dim_solution_, optimality_error_vec)> residual_tolerance_*residual_tolerance_ && j < max_iteration_) {
        forwardDifferenceGMRES(initial_time, initial_state_vec, initial_solution_vec, solution_update_vec_);
        for (int i=0; i<dim_solution_; i++) {
            initial_solution_vec[i] += solution_update_vec_[i];
        }
        computeOptimalityErrors(initial_time, initial_state_vec, initial_solution_vec, optimality_error_vec);
        j++;
    }
}


inline void InitCGMRES::computeOptimalityErrors(const double time_param, const double* state_vec, const double* current_solution_vec, double* optimality_error_vec)
{
    model_.phixFunc(time_param, state_vec, lambda_vec_);
    model_.huFunc(time_param, state_vec, current_solution_vec, lambda_vec_, optimality_error_vec);
}


void InitCGMRES::bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec) 
{
    for (int i=0; i<dim_solution_; i++) {
        incremented_solution_vec_[i] = current_solution_vec[i] + finite_diff_step_*solution_update_vec_[i];
    }
    computeOptimalityErrors(time_param, state_vec, current_solution_vec, error_vec_);
    computeOptimalityErrors(time_param, state_vec, incremented_solution_vec_, error_vec_1_);

    for (int i=0; i<dim_solution_; i++) {
        b_vec[i] = - error_vec_[i] - (error_vec_1_[i] - error_vec_[i])/finite_diff_step_;
    }
}


void InitCGMRES::axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec)
{
    for (int i=0; i<dim_solution_; i++) { 
        incremented_solution_vec_[i] = current_solution_vec[i] + finite_diff_step_*direction_vec[i];
    }
    computeOptimalityErrors(time_param, state_vec, incremented_solution_vec_, error_vec_1_);

    for (int i=0; i<dim_solution_; i++) {
        ax_vec[i] = (error_vec_1_[i] - error_vec_[i])/finite_diff_step_;
    }
}