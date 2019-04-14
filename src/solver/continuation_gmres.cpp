#include "continuation_gmres.hpp"


void ContinuationGMRES::computeOptimalityError(const double time_param, const double* state_vec, const double* current_solution_vec, double* optimality_error_vec)
{
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute the state trajectory over the horizon on the basis of the control_input_vec and the current_state_vec.
    // state_mat_ = state_vec
    for(int i=0; i<dim_state_; i++){
        state_mat_[0][i] = state_vec[i];
    }
    double tau=time_param;
    for(int i=0; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.stateFunc(tau, state_mat_[i], &(current_solution_vec[i*dim_control_input_and_constraints_]), dx_vec_);
        // state_mat[i+1] = state_mat_[i] + delta_tau*dx_vec_
        for(int j=0; j<dim_state_; j++){
            state_mat_[i+1][j] = state_mat_[i][j] + delta_tau*dx_vec_[j];
        }
    }

    // Compute the Lagrange multiplier over the horizon on the basis of the control_input_vec and the current_state_vec.
    model_.phixFunc(tau, state_mat_[horizon_division_num_], lambda_mat_[horizon_division_num_]);
    for(int i=horizon_division_num_-1; i>=0; i--, tau-=delta_tau){
        model_.hxFunc(tau, state_mat_[i], &(current_solution_vec[i*dim_control_input_and_constraints_]), lambda_mat_[i+1], dx_vec_);
        // lambda_mat_[i] = lamdba_mat_[i+1] + delta_tau*dx_vec_
        for(int j=0; j<dim_state_; j++){
            lambda_mat_[i][j] = lambda_mat_[i+1][j] + delta_tau*dx_vec_[j];
        }
    }

    // Compute the optimality condition over the horizon on the basis of the control_input_vec and the current_state_vec.
    tau = time_param;
    for(int i=0; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.huFunc(tau, state_mat_[i], &(current_solution_vec[i*dim_control_input_and_constraints_]), lambda_mat_[i+1], &(optimality_error_vec[i*dim_control_input_and_constraints_]));
    }
}


void ContinuationGMRES::bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec)
{
    // incremented_solution_vec_ = current_solution_vec + difference_incremente_*solution_update_vec_
    for(int i=0; i<dim_solution_; i++){
        incremented_solution_vec_[i] = current_solution_vec[i] + difference_increment_*solution_update_vec_[i];
    }
    computeOptimalityError(time_param, state_vec, current_solution_vec, optimality_vec_);
    computeOptimalityError(incremented_time_, incremented_state_vec_, current_solution_vec, optimality_vec_1_);
    computeOptimalityError(incremented_time_, incremented_state_vec_, incremented_solution_vec_, optimality_vec_2_);

    // b_vec = (1/differenece_incremente_-zeta_)*optimality_vec_ - optimality_vec_2_/difference_increment_
    for(int i=0; i<dim_solution_; i++){
        b_vec[i] = (1/difference_increment_-zeta_)*optimality_vec_[i] - optimality_vec_2_[i]/difference_increment_;
    }
}


inline void ContinuationGMRES::axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec)
{
    // incremented_solution_vec_ = current_solution_vec + difference_increment_*direction_vec
    for(int i=0; i<dim_solution_; i++){
        incremented_solution_vec_[i] = current_solution_vec[i] + difference_increment_*direction_vec[i];
    }
    computeOptimalityError(incremented_time_, incremented_state_vec_, incremented_solution_vec_, optimality_vec_2_);

    // ax_vec = (optimality_vec_2_ - optimality_vec_1_)/difference_increment_
    for(int i=0; i<dim_solution_; i++){
        ax_vec[i] = (optimality_vec_2_[i] - optimality_vec_1_[i])/difference_increment_;
    }
}


ContinuationGMRES::ContinuationGMRES(const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int max_dim_krylov) : MatrixFreeGMRES(), 
    model_(), 
    dim_state_(model_.dimState()), 
    dim_control_input_(model_.dimControlInput()), 
    dim_constraints_(model_.dimConstraints()), 
    dim_control_input_and_constraints_(model_.dimControlInput()+model_.dimConstraints()), 
    dim_solution_(horizon_division_num*(model_.dimControlInput()+model_.dimConstraints())), 
    horizon_division_num_(horizon_division_num), 
    max_dim_krylov_(max_dim_krylov), 
    initial_time_(0), 
    horizon_max_length_(horizon_max_length), 
    alpha_(alpha), 
    zeta_(zeta), 
    difference_increment_(difference_increment), 
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
    setGMRESParams(dim_solution_, max_dim_krylov);
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


void ContinuationGMRES::initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double convergence_radius, const int max_iteration)
{
    double initial_solution_vec[dim_control_input_and_constraints_], initial_error_vec[dim_control_input_and_constraints_];
    InitCGMRES initializer(difference_increment_, max_dim_krylov_);

    initial_time_ = initial_time;
    initializer.solve0stepNOCP(initial_time, initial_state_vec, initial_guess_input_vec, convergence_radius, max_iteration, initial_solution_vec);
    initializer.getOptimalityErrorVec(initial_time, initial_state_vec, initial_solution_vec, initial_error_vec);

    for(int i=0; i<horizon_division_num_; i++){
        // Intialize the solution_vec_.
        // solution_vec_[i] = initial_solution_vec
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            solution_vec_[i*dim_control_input_and_constraints_+j] = initial_solution_vec[j];
        }
        // Intialize the optimality_vec_.
        // optimality_vec_[i] = initial_error_vec
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            optimality_vec_[i*dim_control_input_and_constraints_+j] = initial_error_vec[j];
        }
    }
}


void ContinuationGMRES::controlUpdate(const double current_time, const double sampling_period, const double* current_state_vec, double* optimal_control_input_vec)
{
    // Predict the incremented state.
    incremented_time_ = current_time + difference_increment_;
    model_.stateFunc(current_time, current_state_vec, solution_vec_, dx_vec_);

    // incremented_state_vec_ = current_state_vec + difference_increment_*dx_vec_
    for(int i=0; i<dim_state_; i++){
        incremented_state_vec_[i] = current_state_vec[i] + difference_increment_*dx_vec_[i];
    }

    // Solves the matrix-free GMRES and updates the solution.
    forwardDifferenceGMRES(current_time, current_state_vec, solution_vec_, solution_update_vec_);

    // solution_vec_ += sampling_period*solution_update_vec_
    for(int i=0; i<dim_solution_; i++){
        solution_vec_[i] += sampling_period * solution_update_vec_[i];
    }
    // optimal_control_input_vec = solution_vec
    for(int i=0; i<dim_control_input_; i++){
        optimal_control_input_vec[i] = solution_vec_[i];
    }
}


void ContinuationGMRES::getControlInput(double* control_input_vec) const
{
    // control_input_vec = solution_vec
    for(int i=0; i<dim_control_input_; i++){
        control_input_vec[i] = solution_vec_[i];
    }
}


double ContinuationGMRES::getError(const double current_time, const double* current_state_vec)
{
    double error_vec[dim_solution_];

    computeOptimalityError(current_time, current_state_vec, solution_vec_, error_vec);

    return std::sqrt(linearfunc::squaredNorm(dim_solution_, error_vec));
}