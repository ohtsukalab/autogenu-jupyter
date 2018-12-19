#include "multiple_shooting_cgmres.hpp"



MultipleShootingCGMRES::MultipleShootingCGMRES(const NMPCModel model, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov) : MatrixFreeGMRES((model.dimControlInput()+model.dimConstraints())*horizon_division_num, dim_krylov)
{
    // set dimensions and parameters
    model_ = model;
    dim_state_ = model_.dimState();
    dim_control_input_ = model_.dimControlInput();
    dim_constraints_ = model_.dimConstraints();
    dim_control_input_and_constraints_ = dim_control_input_ + dim_constraints_;
    dim_state_and_lambda_ = 2*dim_state_;
    dim_control_input_and_constraints_seq_ = horizon_division_num * dim_control_input_and_constraints_;
    dim_state_and_lambda_seq_ = horizon_division_num * dim_state_and_lambda_;

    // set parameters for horizon and the C/GMRES
    horizon_max_length_ = horizon_max_length;
    alpha_ = alpha;
    horizon_division_num_ = horizon_division_num;
    difference_increment_ = difference_increment;
    zeta_ = zeta;
    dim_krylov_ = dim_krylov;

    // allocate vectors and matrices
    dx_vec_.resize(dim_state_);
    incremented_state_vec_.resize(dim_state_);
    control_input_and_constraints_seq_.resize(dim_control_input_and_constraints_seq_);
    incremented_control_input_and_constraints_seq_.resize(dim_control_input_and_constraints_seq_);
    control_input_and_constraints_error_seq_.resize(dim_control_input_and_constraints_seq_);
    control_input_and_constraints_error_seq_1_.resize(dim_control_input_and_constraints_seq_);
    control_input_and_constraints_error_seq_2_.resize(dim_control_input_and_constraints_seq_);
    control_input_and_constraints_error_seq_3_.resize(dim_control_input_and_constraints_seq_);
    control_input_and_constraints_update_seq_.resize(dim_control_input_and_constraints_seq_);

    state_mat_.resize(dim_state_, horizon_division_num);
    lambda_mat_.resize(dim_state_, horizon_division_num);
    incremented_state_mat_.resize(dim_state_, horizon_division_num);
    incremented_lambda_mat_.resize(dim_state_, horizon_division_num);
    state_error_mat_.resize(dim_state_, horizon_division_num);
    state_error_mat_1_.resize(dim_state_, horizon_division_num);
    lambda_error_mat_.resize(dim_state_, horizon_division_num);
    lambda_error_mat_1_.resize(dim_state_, horizon_division_num);
    state_update_mat_.resize(dim_state_, horizon_division_num);
    lambda_update_mat_.resize(dim_state_, horizon_division_num);

    // initialize solution of the forward-difference GMRES
    control_input_and_constraints_update_seq_ = Eigen::VectorXd::Zero(dim_control_input_and_constraints_seq_);
    state_update_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num);
    lambda_update_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num);
}


void MultipleShootingCGMRES::initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const double convergence_radius, const int max_iteration)
{
    Eigen::VectorXd initial_control_input_and_constraints_vec(dim_control_input_and_constraints_), initial_control_input_and_constraints_error(dim_control_input_and_constraints_), initial_lambda_vec(dim_state_);
    InitCGMRES initializer(model_, difference_increment_, dim_krylov_);

    initial_time_ = initial_time;

    // intialize the solution
    initializer.solve0stepNOCP(initial_time, initial_state_vec, initial_guess_input_vec, convergence_radius, max_iteration, initial_control_input_and_constraints_vec);
    model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
    for(int i=0; i<horizon_division_num_; i++){
        control_input_and_constraints_seq_.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_) = initial_control_input_and_constraints_vec;
        state_mat_.col(i) = initial_state_vec;
        lambda_mat_.col(i) = initial_lambda_vec;
    }

    // intialize the error
    initial_control_input_and_constraints_error = initializer.getOptimalityErrorVec(initial_time, initial_state_vec, initial_control_input_and_constraints_vec);
    for(int i=0; i<horizon_division_num_; i++){
        control_input_and_constraints_error_seq_.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_) = initial_control_input_and_constraints_error;
    }
    state_error_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num_);
    lambda_error_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num_);
}


void MultipleShootingCGMRES::controlUpdate(const double current_time, const double sampling_period, const Eigen::VectorXd& current_state_vec, Eigen::Ref<Eigen::VectorXd> optimal_control_input_vec)
{
    // predict the incremented state
    incremented_time_ = current_time + difference_increment_;
    model_.stateFunc(current_time, current_state_vec, control_input_and_constraints_seq_.segment(0, dim_control_input_), dx_vec_);
    incremented_state_vec_ = current_state_vec + difference_increment_*dx_vec_;

    forwardDifferenceGMRES(current_time, current_state_vec, control_input_and_constraints_seq_, control_input_and_constraints_update_seq_);
    control_input_and_constraints_seq_ += sampling_period * control_input_and_constraints_update_seq_;

    // update state_mat_ and lamdba_mat
    computeStateAndLambda(incremented_time_, incremented_state_vec_, control_input_and_constraints_seq_, (1-difference_increment_*zeta_)*state_error_mat_, (1-difference_increment_*zeta_)*lambda_error_mat_, incremented_state_mat_, incremented_lambda_mat_);
    state_update_mat_ = (incremented_state_mat_ - state_mat_)/difference_increment_;
    lambda_update_mat_ = (incremented_lambda_mat_ - lambda_mat_)/difference_increment_;
    state_mat_ += sampling_period * state_update_mat_;
    lambda_mat_ += sampling_period * lambda_update_mat_;

    optimal_control_input_vec = control_input_and_constraints_seq_.segment(0, dim_control_input_);
}


double MultipleShootingCGMRES::getError(const double current_time, const Eigen::VectorXd& current_state_vec)
{
    Eigen::VectorXd control_input_and_constraints_error_seq(horizon_division_num_*dim_control_input_and_constraints_); 
    Eigen::MatrixXd state_error_mat(dim_state_, horizon_division_num_), lambda_error_mat(dim_state_, horizon_division_num_);

    computeOptimalityErrorforControlInputAndConstraints(current_time, current_state_vec, control_input_and_constraints_seq_, state_mat_, lambda_mat_, control_input_and_constraints_error_seq);
    computeOptimalityErrorforStateAndLambda(current_time, current_state_vec, control_input_and_constraints_seq_, state_mat_, lambda_mat_, state_error_mat, lambda_error_mat);

    double squared_error = control_input_and_constraints_error_seq.squaredNorm();
    for(int i=0; i<horizon_division_num_; i++){
        squared_error += (state_error_mat.col(i).squaredNorm() + lambda_error_mat.col(i).squaredNorm());
    }

    std::cout << "control_error = " << control_input_and_constraints_error_seq.norm() << ", state_error = " << std::sqrt(squared_error - control_input_and_constraints_error_seq.squaredNorm()) << std::endl;

    if(std::isnan(control_input_and_constraints_error_seq.norm()) || std::isnan(squared_error)){
        std::exit(1);
    }

    return std::sqrt(squared_error);
}


inline void MultipleShootingCGMRES::computeOptimalityErrorforControlInputAndConstraints(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& state_mat, const Eigen::MatrixXd& lambda_mat, Eigen::Ref<Eigen::VectorXd> optimality_for_control_input_and_constraints)
{
    // set the horizon
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    model_.huFunc(time_param, state_vec, control_input_and_constraints_seq.segment(0, dim_control_input_and_constraints_), lambda_mat.col(0), optimality_for_control_input_and_constraints.segment(0, dim_control_input_and_constraints_));
    double tau = time_param+delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.huFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_), lambda_mat.col(i), optimality_for_control_input_and_constraints.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_));
    }
}


inline void MultipleShootingCGMRES::computeOptimalityErrorforStateAndLambda(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& state_mat, const Eigen::MatrixXd& lambda_mat, Eigen::Ref<Eigen::MatrixXd> optimality_for_state, Eigen::Ref<Eigen::MatrixXd> optimality_for_lambda)
{
    // set the horizon
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // compute optimality error for state
    model_.stateFunc(time_param, state_vec, control_input_and_constraints_seq.segment(0, dim_control_input_), dx_vec_);
    optimality_for_state.col(0) = state_mat.col(0) - state_vec - delta_tau * dx_vec_;
    double tau = time_param + delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.stateFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_), dx_vec_);
        optimality_for_state.col(i) = state_mat.col(i) - state_mat.col(i-1) - delta_tau * dx_vec_;
    }

    // compute optimality error for lambda
    model_.phixFunc(tau, state_mat.col(horizon_division_num_-1), dx_vec_);
    optimality_for_lambda.col(horizon_division_num_-1) = lambda_mat.col(horizon_division_num_-1) - dx_vec_;
    for(int i=horizon_division_num_-1; i>=1; i--, tau-=delta_tau){
        model_.hxFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_), lambda_mat.col(i), dx_vec_);
        optimality_for_lambda.col(i-1) = lambda_mat.col(i-1) - lambda_mat.col(i) - delta_tau * dx_vec_;
    }
}


inline void MultipleShootingCGMRES::computeStateAndLambda(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& optimality_for_state, const Eigen::MatrixXd& optimality_for_lambda, Eigen::Ref<Eigen::MatrixXd> state_mat, Eigen::Ref<Eigen::MatrixXd> lambda_mat)
{
    // set the horizon
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // compute state
    model_.stateFunc(time_param, state_vec, control_input_and_constraints_seq.segment(0, dim_control_input_), dx_vec_);
    state_mat.col(0) = state_vec + delta_tau * dx_vec_ + optimality_for_state.col(0);
    double tau = time_param + delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.stateFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_), dx_vec_);
        state_mat.col(i) = state_mat.col(i-1) + delta_tau * dx_vec_ + optimality_for_state.col(i);
    }

    // compute lambda
    model_.phixFunc(tau, state_mat.col(horizon_division_num_-1), dx_vec_);
    lambda_mat.col(horizon_division_num_-1) = dx_vec_ + optimality_for_lambda.col(horizon_division_num_-1);
    for(int i=horizon_division_num_-1; i>=1; i--, tau-=delta_tau){
        model_.hxFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_), lambda_mat.col(i), dx_vec_);
        lambda_mat.col(i-1) = lambda_mat.col(i) + delta_tau * dx_vec_ + optimality_for_lambda.col(i-1);
    }
}


void MultipleShootingCGMRES::nonlinearEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec)
{
    // control_input_and_constraints_error_seq_ : F(U, X, x, t)
    computeOptimalityErrorforControlInputAndConstraints(time_param, state_vec, current_solution_vec, state_mat_, lambda_mat_, control_input_and_constraints_error_seq_);
    // control_input_and_constraints_error_seq_1_: F(U, X, x+hdx, t+ht)
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, control_input_and_constraints_error_seq_1_);

    // state_error_mat_, lambda_error_mat_ : G(U, X, x, t)
    computeOptimalityErrorforStateAndLambda(time_param, state_vec, current_solution_vec, state_mat_, lambda_mat_, state_error_mat_, lambda_error_mat_);
    // state_error_mat_1_, lambda_error_mat_1_ : G(U, X, x+hdx, t+ht)
    computeOptimalityErrorforStateAndLambda(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, state_error_mat_1_, lambda_error_mat_1_);

    // incremented_state_mat_, incremented_lambda_mat_ : K(U, (1-hz)G(U, X, x, t), x+hdx, t+ht)
    computeStateAndLambda(incremented_time_, incremented_state_vec_, current_solution_vec, (1-difference_increment_*zeta_)*state_error_mat_, (1-difference_increment_*zeta_)*lambda_error_mat_, incremented_state_mat_, incremented_lambda_mat_);
    // control_input_and_constraints_error_seq_3_ : F(U, K(U, (1-hz)G(U, X, x, t), x+hdx, t+ht), x+hdx, t+ht)
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_3_);

    incremented_control_input_and_constraints_seq_ = current_solution_vec + difference_increment_ * control_input_and_constraints_update_seq_;
    // incremented_state_mat_, incremented_lambda_mat_ : K(U+h\hat{U}, G(U, X, x+hdx, t+ht), x+hdx, t+ht)
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    // F(U+h\hat{U}, K(U+h\hat{U}, G(U, X, x+hdx, t+ht), x+hdx, t+ht), x+hdx, t+ht)
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);

    equation_error_vec = - (zeta_ - 1/difference_increment_) * control_input_and_constraints_error_seq_ - control_input_and_constraints_error_seq_3_/difference_increment_ - (control_input_and_constraints_error_seq_2_ - control_input_and_constraints_error_seq_1_)/difference_increment_;
}


void MultipleShootingCGMRES::forwardDifferenceEquation(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec)
{
    incremented_control_input_and_constraints_seq_ = current_solution_vec + difference_increment_ * direction_vec;
    // incremented_state_mat_, incremented_lambda_mat_ : K(U+huk, G(U, X, x+hdx, t+ht), x+hdx, t+ht)
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    // F(U+huk, K(U+huk, G(U, X, x+hdx, t+ht), x+hdx, t+ht), x+hdx, t+ht)
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);

    // (F(U+huk, K(U+huk, G(U, X, x+hdx, t+ht), x+hdx, t+ht), x+hdx, t+ht) - F(U, x, x+hdx, t+ht))/h
    forward_difference_error_vec =(control_input_and_constraints_error_seq_2_ - control_input_and_constraints_error_seq_1_)/difference_increment_;
}