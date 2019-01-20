#include "init_cgmres_with_saturation.hpp"


inline void InitCGMRESWithSaturation::addHamiltonianDerivativeWithControlInput(const Eigen::VectorXd& control_input_and_constraints_vec, const Eigen::VectorXd& saturation_lagrange_multiplier_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_control_input_and_constraints_vec)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_control_input_and_constraints_vec(saturation_seq_.index(i)) += (2*control_input_and_constraints_vec(saturation_seq_.index(i)) - saturation_seq_.min(i) - saturation_seq_.max(i)) * saturation_lagrange_multiplier_vec(i);
    }
}


inline void InitCGMRESWithSaturation::computeDummyOptimality(const Eigen::VectorXd& dummy_input_vec, const Eigen::VectorXd& saturation_lagrange_multiplier_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_dummy)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_dummy(i) = 2 * (saturation_seq_.quadratic_weight(i) + saturation_lagrange_multiplier_vec(i)) * dummy_input_vec(i) - saturation_seq_.dummy_weight(i);
    }
}


inline void InitCGMRESWithSaturation::computeSaturationOptimality(const Eigen::VectorXd& control_input_and_constraint_vec, const Eigen::VectorXd& dummy_input_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_saturation)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_saturation(i) = control_input_and_constraint_vec(saturation_seq_.index(i)) * (control_input_and_constraint_vec(saturation_seq_.index(i)) - saturation_seq_.min(i) - saturation_seq_.max(i)) + saturation_seq_.min(i) * saturation_seq_.max(i) + dummy_input_vec(i) * dummy_input_vec(i);
    }
}


inline void InitCGMRESWithSaturation::computeOptimalityErrors(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& solution_vec, Eigen::Ref<Eigen::VectorXd> optimality_vec)
{
    model_.phixFunc(time_param, state_vec, lambda_vec_);
    model_.huFunc(time_param, state_vec, solution_vec.segment(0, dim_control_input_and_constraints_), lambda_vec_, optimality_vec.segment(0, dim_control_input_and_constraints_));
    addHamiltonianDerivativeWithControlInput(solution_vec.segment(0, dim_control_input_and_constraints_), solution_vec.segment(dim_control_input_and_constraints_+dim_saturation_, dim_saturation_), optimality_vec.segment(0, dim_control_input_and_constraints_));

    computeDummyOptimality(solution_vec.segment(dim_control_input_and_constraints_,dim_saturation_), solution_vec.segment(dim_control_input_and_constraints_+dim_saturation_,dim_saturation_), optimality_vec.segment(dim_control_input_and_constraints_, dim_saturation_));
    computeSaturationOptimality(solution_vec.segment(0, dim_control_input_and_constraints_), solution_vec.segment(dim_control_input_and_constraints_, dim_saturation_), optimality_vec.segment(dim_control_input_and_constraints_+dim_saturation_, dim_saturation_));
}


void InitCGMRESWithSaturation::bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> b_vec) 
{
    computeOptimalityErrors(time_param, state_vec, current_solution_vec, error_vec_);
    computeOptimalityErrors(time_param, state_vec, current_solution_vec+difference_increment_*solution_update_vec_, error_vec_1_);

    b_vec = - error_vec_ - (error_vec_1_ - error_vec_)/difference_increment_;
}


void InitCGMRESWithSaturation::axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> ax_vec)
{
    computeOptimalityErrors(time_param, state_vec, current_solution_vec+difference_increment_*direction_vec, error_vec_1_);

    ax_vec = (error_vec_1_ - error_vec_)/difference_increment_;
}



InitCGMRESWithSaturation::InitCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, const double difference_increment, const int max_dim_krylov) : MatrixFreeGMRES(), 
    model_(), 
    saturation_seq_(saturation_seq), 
    dim_control_input_and_constraints_(model_.dimControlInput()+model_.dimConstraints()), 
    dim_saturation_(saturation_seq.dimSaturation()), 
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()+2*saturation_seq.dimSaturation()), 
    difference_increment_(difference_increment), 
    solution_update_vec_(Eigen::VectorXd::Zero(dim_solution_)), 
    lambda_vec_(Eigen::VectorXd::Zero(model_.dimState())), 
    error_vec_(Eigen::VectorXd::Zero(dim_solution_)), 
    error_vec_1_(Eigen::VectorXd::Zero(dim_solution_)), 
    error_vec_2_(Eigen::VectorXd::Zero(dim_solution_))
{
    // Set parameters in GMRES.
    setGMRESParams(dim_solution_, max_dim_krylov);
}


void InitCGMRESWithSaturation::solve0stepNOCP(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_control_input_vec, const Eigen::VectorXd& initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration, Eigen::Ref<Eigen::VectorXd> solution_vec)
{
    // Substitute initial guess solution to solution_vec.
    solution_vec.segment(0, dim_control_input_and_constraints_) = initial_guess_control_input_vec;
    for(int i=0; i<dim_saturation_; i++){
        solution_vec(dim_control_input_and_constraints_+i) = std::sqrt((saturation_seq_.max(i)-saturation_seq_.min(i))*(saturation_seq_.max(i)-saturation_seq_.min(i))/4 - (initial_guess_control_input_vec(saturation_seq_.index(i)) 
            - (saturation_seq_.max(i)+saturation_seq_.min(i))/2)*(initial_guess_control_input_vec(saturation_seq_.index(i)) - (saturation_seq_.max(i)+saturation_seq_.min(i))/2));
    }
    solution_vec.segment(dim_control_input_and_constraints_+dim_saturation_, dim_saturation_) = initial_guess_lagrange_multiplier;

    // Solve the 0step nonlinear optimal control problem using Newton GMRES method.
    computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);

    int i=0;
    while(error_vec_.squaredNorm() > convergence_radius*convergence_radius && i < max_iteration){
        forwardDifferenceGMRES(initial_time, initial_state_vec, solution_vec, solution_update_vec_);
        solution_vec += solution_update_vec_;
        computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);
        i++;
    }
}


Eigen::VectorXd InitCGMRESWithSaturation::getControlInputAndConstraintsError(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec)
{
    Eigen::VectorXd error_vec(dim_solution_);

    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);

    return error_vec.segment(0, dim_control_input_and_constraints_);
}


Eigen::VectorXd InitCGMRESWithSaturation::getDummyInputError(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec)
{
    Eigen::VectorXd error_vec(dim_solution_);

    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);

    return error_vec.segment(dim_control_input_and_constraints_, dim_saturation_);
}


Eigen::VectorXd InitCGMRESWithSaturation::getControlInputSaturationError(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec)
{
    Eigen::VectorXd error_vec(dim_solution_);

    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);

    return error_vec.segment(dim_control_input_and_constraints_+dim_saturation_, dim_saturation_);
}