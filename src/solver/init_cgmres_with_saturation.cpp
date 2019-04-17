#include "init_cgmres_with_saturation.hpp"


inline void InitCGMRESWithSaturation::addHamiltonianDerivativeWithControlInput(const double* control_input_and_constraints_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_control_input_and_constraints_vec)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_control_input_and_constraints_vec[saturation_seq_.index(i)] += [2*control_input_and_constraints_vec(saturation_seq_.index(i)] - saturation_seq_.min(i) - saturation_seq_.max(i)) * saturation_lagrange_multiplier_vec[i];
    }
}


inline void InitCGMRESWithSaturation::computeDummyOptimality(const double* dummy_input_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_dummy)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_dummy[i] = 2 * (saturation_seq_.quadratic_weight(i) + saturation_lagrange_multiplier_vec[i]) * dummy_input_vec[i] - saturation_seq_.dummy_weight(i);
    }
}


inline void InitCGMRESWithSaturation::computeSaturationOptimality(const double* control_input_and_constraint_vec, const double* dummy_input_vec, double* optimality_for_saturation)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_saturation[i] = control_input_and_constraint_vec[saturation_seq_.index(i)] * (control_input_and_constraint_vec[saturation_seq_.index(i)] - saturation_seq_.min(i) - saturation_seq_.max(i)) + saturation_seq_.min(i) * saturation_seq_.max(i) + dummy_input_vec[i] * dummy_input_vec[i];
    }
}


inline void InitCGMRESWithSaturation::computeOptimalityErrors(const double time_param, const double* state_vec, const double* solution_vec, double* optimality_vec)
{
    model_.phixFunc(time_param, state_vec, lambda_vec_);
    model_.huFunc(time_param, state_vec, solution_vec, lambda_vec_, optimality_vec);
    addHamiltonianDerivativeWithControlInput(solution_vec, &(solution_vec[dim_control_input_and_constraints_+dim_saturation_]), optimality_vec);

    computeDummyOptimality(&(solution_vec[dim_control_input_and_constraints_]), &(solution_vec[dim_control_input_and_constraints_+dim_saturation_]), &(optimality_vec[dim_control_input_and_constraints_]));
    computeSaturationOptimality(solution_vec, &(solution_vec[dim_control_input_and_constraints_]), &(optimality_vec[dim_control_input_and_constraints_+dim_saturation_]));
}


void InitCGMRESWithSaturation::bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec) 
{
    computeOptimalityErrors(time_param, state_vec, current_solution_vec, error_vec_);
    for(int i=0; i<dim_solution_; i++){
        incremented_solution_vec_[i] = current_solution_vec[i] + difference_increment_*solution_update_vec_[i];
    }
    computeOptimalityErrors(time_param, state_vec, incremented_solution_vec_, error_vec_1_);

    for(int i=0; i<dim_solution_; i++){
        b_vec[i] = - error_vec_[i] - (error_vec_1_[i] - error_vec_[i])/difference_increment_;
    }
}


void InitCGMRESWithSaturation::axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec)
{
    for(int i=0; i<dim_solution_; i++){
        incremented_solution_vec_[i] = current_solution_vec[i] + difference_increment_*direction_vec[i];
    }
    computeOptimalityErrors(time_param, state_vec, incremented_solution_vec_, error_vec_1_);

    for(int i=0; i<dim_solution_; i++){
        ax_vec[i] = (error_vec_1_[i] - error_vec_[i])/difference_increment_;
    }
}



InitCGMRESWithSaturation::InitCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, const double difference_increment, const int max_dim_krylov) : MatrixFreeGMRES(), 
    model_(), 
    saturation_seq_(saturation_seq), 
    dim_control_input_and_constraints_(model_.dimControlInput()+model_.dimConstraints()), 
    dim_saturation_(saturation_seq.dimSaturation()), 
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()+2*saturation_seq.dimSaturation()), 
    difference_increment_(difference_increment), 
    incremented_solution_vec_(linearfunc::newVector(dim_solution_)), 
    solution_update_vec_(linearfunc::newVector(dim_solution_)), 
    lambda_vec_(linearfunc::newVector(model_.dimState())), 
    error_vec_(linearfunc::newVector(dim_solution_)), 
    error_vec_1_(linearfunc::newVector(dim_solution_)), 
    error_vec_2_(linearfunc::newVector(dim_solution_))
{
    // Set parameters in GMRES.
    setGMRESParams(dim_solution_, max_dim_krylov);
}


InitCGMRESWithSaturation::InitCGMRESWithSaturation()
{
    linearfunc::deleteVector(incremented_solution_vec_);
    linearfunc::deleteVector(solution_update_vec_);
    linearfunc::deleteVector(lambda_vec_);
    linearfunc::deleteVector(error_vec_);
    linearfunc::deleteVector(error_vec_1_);
    linearfunc::deleteVector(error_vec_2_);
}



void InitCGMRESWithSaturation::solve0stepNOCP(const double initial_time, const double* initial_state_vec, const double* initial_guess_control_input_vec, const double* initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration, double* solution_vec)
{
    // Substitute initial guess solution to solution_vec.
    for(int i=0; i<dim_control_input_and_constraints_; i++){
        solution_vec_[i] = initial_guess_control_input_vec[i]
    }
    for(int i=0; i<dim_saturation_; i++){
        solution_vec[dim_control_input_and_constraints_+i] = std::sqrt((saturation_seq_.max(i)-saturation_seq_.min(i))*(saturation_seq_.max(i)-saturation_seq_.min(i))/4 - (initial_guess_control_input_vec[saturation_seq_.index(i)] - (saturation_seq_.max(i)+saturation_seq_.min(i))/2)*(initial_guess_control_input_vec[saturation_seq_.index(i)] - (saturation_seq_.max(i)+saturation_seq_.min(i))/2));
    }
    for(int i=0; i<dim_saturation_; i++){
        solution_vec_[dim_control_input_and_constraints_+dim_saturation_+i] = initial_guess_lagrange_multiplier[i];
    }

    // Solve the 0step nonlinear optimal control problem using Newton GMRES method.
    computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);

    int i=0;
    while(linearfunc::squaredNorm(dim_solution_, error_vec_) > convergence_radius*convergence_radius && i < max_iteration){
        forwardDifferenceGMRES(initial_time, initial_state_vec, solution_vec, solution_update_vec_);
        for(int j=0; j<dim_solution_; j++){
            solution_vec[j] += solution_update_vec_[j];
        }
        computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);
        i++;
    }
}


void InitCGMRESWithSaturation::getControlInputAndConstraintsError(const double initial_time, const double* initial_state_vec, const double* current_solution_vec, double* control_input_and_constraints_error_vec)
{
    double error_vec[dim_solution_];

    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);

    for(int i=0; i<dim_control_input_and_constraints_; i++){
        control_input_and_constraints_error_vec[i] = error_vec[i];
    }
}


void InitCGMRESWithSaturation::getDummyInputError(const double initial_time, const double* initial_state_vec, const double* current_solution_vec, double* dummy_input_error_vec)
{
    double error_vec[dim_solution_];

    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);

    for(int i=0; i<dim_saturation_; i++){
        dummy_input_error_vec[i] = error_vec[dim_control_input_and_constraints_+i];
    }
}


void InitCGMRESWithSaturation::getControlInputSaturationError(const double initial_time, const double* initial_state_vec, const double* current_solution_vec,  double* control_input_saturation_error_vec)
{
    double error_vec[dim_solution_];

    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);

    for(int i=0; i<dim_saturation_; i++){
        control_input_saturation_error_vec[i] = error_vec[dim_control_input_and_constraints_+dim_saturation_+i];
    }
}