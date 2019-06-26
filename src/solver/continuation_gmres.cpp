#include "continuation_gmres.hpp"


ContinuationGMRES::ContinuationGMRES(const double T_f, const double alpha, const int horizon_division_num, const double finite_diff_step, const double zeta, const int kmax) : MatrixFreeGMRES(), 
    model_(), 
    cgmres_initializer_(),
    dim_state_(model_.dimState()), 
    dim_control_input_(model_.dimControlInput()), 
    dim_constraints_(model_.dimConstraints()), 
    dim_control_input_and_constraints_(model_.dimControlInput()+model_.dimConstraints()), 
    dim_solution_(horizon_division_num*(model_.dimControlInput()+model_.dimConstraints())), 
    horizon_division_num_(horizon_division_num), 
    kmax_(kmax), 
    initial_time_(0), 
    T_f_(T_f), 
    alpha_(alpha), 
    zeta_(zeta), 
    finite_diff_step_(finite_diff_step), 
    incremented_time_(0), 
    dx_vec_(linearfunc::newVector(dim_state_)), 
    incremented_state_vec_(linearfunc::newVector(dim_state_)),
    solution_vec_(linearfunc::newVector(dim_solution_)), 
    incremented_solution_vec_(linearfunc::newVector(dim_solution_)), 
    optimality_vec_(linearfunc::newVector(dim_solution_)), 
    optimality_vec_1_(linearfunc::newVector(dim_solution_)), 
    optimality_vec_2_(linearfunc::newVector(dim_solution_)), 
    solution_update_vec_(linearfunc::newVector(dim_solution_)), 
    state_mat_(linearfunc::newMatrix(horizon_division_num+1, dim_state_)), 
    lambda_mat_(linearfunc::newMatrix(horizon_division_num+1, dim_state_))
{
    // Set dimensions and parameters in GMRES.
    setGMRESParams(dim_solution_, kmax);
}


ContinuationGMRES::~ContinuationGMRES()
{
    linearfunc::deleteVector(dx_vec_);
    linearfunc::deleteVector(incremented_state_vec_);
    linearfunc::deleteVector(solution_vec_);
    linearfunc::deleteVector(incremented_solution_vec_);
    linearfunc::deleteVector(optimality_vec_);
    linearfunc::deleteVector(optimality_vec_1_);
    linearfunc::deleteVector(optimality_vec_2_);
    linearfunc::deleteVector(solution_update_vec_);
    linearfunc::deleteMatrix(state_mat_);
    linearfunc::deleteMatrix(lambda_mat_);
}


void ContinuationGMRES::setInitParams(const double* initial_guess_solution, const double residual_tolerance, const int max_iteration, const double finite_diff_step, const int kmax)
{
    cgmres_initializer_.setInitParams(initial_guess_solution, residual_tolerance, max_iteration, finite_diff_step, kmax);
}


void ContinuationGMRES::initSolution(const double initial_time, const double* initial_state_vec, double* optimal_control_input_vec)
{
    double initial_solution_vec[dim_control_input_and_constraints_], initial_error_vec[dim_control_input_and_constraints_];
    initial_time_ = initial_time;
    cgmres_initializer_.solveOCPForInit(initial_time, initial_state_vec, initial_solution_vec, initial_error_vec);
    for (int i=0; i<horizon_division_num_; i++) {
        // Intialize the solution_vec_.
        for (int j=0; j<dim_control_input_and_constraints_; j++) {
            solution_vec_[i*dim_control_input_and_constraints_+j] = initial_solution_vec[j];
        }
        // Intialize the optimality_vec_.
        for (int j=0; j<dim_control_input_and_constraints_; j++) {
            optimality_vec_[i*dim_control_input_and_constraints_+j] = initial_error_vec[j];
        }
    }
    for (int i=0; i<dim_control_input_; i++) {
        optimal_control_input_vec[i] = initial_solution_vec[i];
    }
}


void ContinuationGMRES::controlUpdate(const double current_time, const double sampling_period, const double* current_state_vec, double* optimal_control_input_vec)
{
    // Predict the incremented state.
    incremented_time_ = current_time + finite_diff_step_;
    model_.stateFunc(current_time, current_state_vec, solution_vec_, dx_vec_);
    for (int i=0; i<dim_state_; i++) {
        incremented_state_vec_[i] = current_state_vec[i] + finite_diff_step_ * dx_vec_[i];
    }
    forwardDifferenceGMRES(current_time, current_state_vec, solution_vec_, solution_update_vec_);
    for (int i=0; i<dim_solution_; i++) {
        solution_vec_[i] += sampling_period * solution_update_vec_[i];
    }
    for (int i=0; i<dim_control_input_; i++) {
        optimal_control_input_vec[i] = solution_vec_[i];
    }
}


double ContinuationGMRES::getError(const double current_time, const double* current_state_vec)
{
    double error_vec[dim_solution_];
    computeOptimalityError(current_time, current_state_vec, solution_vec_, error_vec);
    return std::sqrt(linearfunc::squaredNorm(dim_solution_, error_vec));
}


void ContinuationGMRES::computeOptimalityError(const double time_param, const double* state_vec, const double* current_solution_vec, double* optimality_error_vec)
{
    // Set and discretize the horizon.
    double horizon_length = T_f_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute the state trajectory over the horizon on the basis of the control_input_vec and the current_state_vec.
    for (int i=0; i<dim_state_; i++) {
        state_mat_[0][i] = state_vec[i];
    }
    double tau=time_param;
    for (int i=0; i<horizon_division_num_; i++, tau+=delta_tau) {
        model_.stateFunc(tau, state_mat_[i], &(current_solution_vec[i*dim_control_input_and_constraints_]), dx_vec_);
        for (int j=0; j<dim_state_; j++) {
            state_mat_[i+1][j] = state_mat_[i][j] + delta_tau * dx_vec_[j];
        }
    }

    // Compute the Lagrange multiplier over the horizon on the basis of the control_input_vec and the current_state_vec.
    model_.phixFunc(tau, state_mat_[horizon_division_num_], lambda_mat_[horizon_division_num_]);
    for (int i=horizon_division_num_-1; i>=0; i--, tau-=delta_tau) {
        model_.hxFunc(tau, state_mat_[i], &(current_solution_vec[i*dim_control_input_and_constraints_]), lambda_mat_[i+1], dx_vec_);
        for (int j=0; j<dim_state_; j++) {
            lambda_mat_[i][j] = lambda_mat_[i+1][j] + delta_tau * dx_vec_[j];
        }
    }

    // Compute the optimality condition over the horizon on the basis of the control_input_vec and the current_state_vec.
    tau = time_param;
    for (int i=0; i<horizon_division_num_; i++, tau+=delta_tau) {
        model_.huFunc(tau, state_mat_[i], &(current_solution_vec[i*dim_control_input_and_constraints_]), lambda_mat_[i+1], &(optimality_error_vec[i*dim_control_input_and_constraints_]));
    }
}


void ContinuationGMRES::bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec)
{
    for (int i=0; i<dim_solution_; i++) {
        incremented_solution_vec_[i] = current_solution_vec[i] + finite_diff_step_*solution_update_vec_[i];
    }
    computeOptimalityError(time_param, state_vec, current_solution_vec, optimality_vec_);
    computeOptimalityError(incremented_time_, incremented_state_vec_, current_solution_vec, optimality_vec_1_);
    computeOptimalityError(incremented_time_, incremented_state_vec_, incremented_solution_vec_, optimality_vec_2_);

    for (int i=0; i<dim_solution_; i++) {
        b_vec[i] = (1/finite_diff_step_ - zeta_) * optimality_vec_[i] - optimality_vec_2_[i]/finite_diff_step_;
    }
}


inline void ContinuationGMRES::axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec)
{
    for (int i=0; i<dim_solution_; i++) {
        incremented_solution_vec_[i] = current_solution_vec[i] + finite_diff_step_ * direction_vec[i];
    }
    computeOptimalityError(incremented_time_, incremented_state_vec_, incremented_solution_vec_, optimality_vec_2_);

    for (int i=0; i<dim_solution_; i++) {
        ax_vec[i] = (optimality_vec_2_[i] - optimality_vec_1_[i])/finite_diff_step_;
    }
}