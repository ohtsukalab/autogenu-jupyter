#include "multiple_shooting_cgmres.hpp"

void MultipleShootingCGMRES::setDummy(double dummy)
{
    double new_dummy = dummy;
}

void MultipleShootingCGMRES::setSolver(NMPCModel model, double horizon_max_length, double alpha, int horizon_division_num, double difference_increment, double zeta, int dim_krylov)
{
    /* Copied from matrix-free GMRES */
        dim_equation_ = (model.dimControlInput()+model.dimConstraints())*horizon_division_num;
        max_dim_krylov_ = dim_krylov;

        // Allocate matrices and vectors for matrix-free GMRES.
        hessenberg_mat_.resize(max_dim_krylov_+1, max_dim_krylov_+1);
        basis_mat_.resize(dim_equation_, max_dim_krylov_+1);
        b_vec_.resize(dim_equation_);
        givens_c_vec_.resize(max_dim_krylov_+1);
        givens_s_vec_.resize(max_dim_krylov_+1);
        g_vec_.resize(max_dim_krylov_+1);


    // Set dimensions and parameters.
    model_ = model;
    dim_state_ = model_.dimState();
    dim_control_input_ = model_.dimControlInput();
    dim_constraints_ = model_.dimConstraints();
    dim_control_input_and_constraints_ = dim_control_input_ + dim_constraints_;
    dim_state_and_lambda_ = 2*dim_state_;
    dim_control_input_and_constraints_seq_ = horizon_division_num * dim_control_input_and_constraints_;
    dim_state_and_lambda_seq_ = horizon_division_num * dim_state_and_lambda_;

    // Set parameters for horizon and the C/GMRES.
    horizon_max_length_ = horizon_max_length;
    alpha_ = alpha;
    horizon_division_num_ = horizon_division_num;
    difference_increment_ = difference_increment;
    zeta_ = zeta;
    dim_krylov_ = dim_krylov;

    // Allocate vectors and matrices.
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

    // Initialize solution of the forward-difference GMRES.
    control_input_and_constraints_update_seq_ = Eigen::VectorXd::Zero(dim_control_input_and_constraints_seq_);
    state_update_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num);
    lambda_update_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num);
}

// MultipleShootingCGMRES::MultipleShootingCGMRES(
    //     const NMPCModel model, 
    //     const double horizon_max_length, 
    //     const double alpha, 
    //     const int horizon_division_num, 
    //     const double difference_increment, 
    //     const double zeta, const int dim_krylov) : 
    //     MatrixFreeGMRES((model.dimControlInput()+model.dimConstraints())*horizon_division_num, dim_krylov)
    // {
    //     // Set dimensions and parameters.
    //     model_ = model;
    //     dim_state_ = model_.dimState();
    //     dim_control_input_ = model_.dimControlInput();
    //     dim_constraints_ = model_.dimConstraints();
    //     dim_control_input_and_constraints_ = dim_control_input_ + dim_constraints_;
    //     dim_state_and_lambda_ = 2*dim_state_;
    //     dim_control_input_and_constraints_seq_ = horizon_division_num * dim_control_input_and_constraints_;
    //     dim_state_and_lambda_seq_ = horizon_division_num * dim_state_and_lambda_;

    //     // Set parameters for horizon and the C/GMRES.
    //     horizon_max_length_ = horizon_max_length;
    //     alpha_ = alpha;
    //     horizon_division_num_ = horizon_division_num;
    //     difference_increment_ = difference_increment;
    //     zeta_ = zeta;
    //     dim_krylov_ = dim_krylov;

    //     // Allocate vectors and matrices.
    //     dx_vec_.resize(dim_state_);
    //     incremented_state_vec_.resize(dim_state_);
    //     control_input_and_constraints_seq_.resize(dim_control_input_and_constraints_seq_);
    //     incremented_control_input_and_constraints_seq_.resize(dim_control_input_and_constraints_seq_);
    //     control_input_and_constraints_error_seq_.resize(dim_control_input_and_constraints_seq_);
    //     control_input_and_constraints_error_seq_1_.resize(dim_control_input_and_constraints_seq_);
    //     control_input_and_constraints_error_seq_2_.resize(dim_control_input_and_constraints_seq_);
    //     control_input_and_constraints_error_seq_3_.resize(dim_control_input_and_constraints_seq_);
    //     control_input_and_constraints_update_seq_.resize(dim_control_input_and_constraints_seq_);

    //     state_mat_.resize(dim_state_, horizon_division_num);
    //     lambda_mat_.resize(dim_state_, horizon_division_num);
    //     incremented_state_mat_.resize(dim_state_, horizon_division_num);
    //     incremented_lambda_mat_.resize(dim_state_, horizon_division_num);
    //     state_error_mat_.resize(dim_state_, horizon_division_num);
    //     state_error_mat_1_.resize(dim_state_, horizon_division_num);
    //     lambda_error_mat_.resize(dim_state_, horizon_division_num);
    //     lambda_error_mat_1_.resize(dim_state_, horizon_division_num);
    //     state_update_mat_.resize(dim_state_, horizon_division_num);
    //     lambda_update_mat_.resize(dim_state_, horizon_division_num);

    //     // Initialize solution of the forward-difference GMRES.
    //     control_input_and_constraints_update_seq_ = Eigen::VectorXd::Zero(dim_control_input_and_constraints_seq_);
    //     state_update_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num);
    //     lambda_update_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num);
// }

void MultipleShootingCGMRES::initSolution(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_input_vec, const double convergence_radius, const int max_iteration)
{
    Eigen::VectorXd initial_control_input_and_constraints_vec(dim_control_input_and_constraints_), initial_control_input_and_constraints_error(dim_control_input_and_constraints_), initial_lambda_vec(dim_state_);
    InitCGMRES initializer(model_, difference_increment_, dim_krylov_);
    initial_time_ = initial_time;

    // Intialize the solution
    initializer.solve0stepNOCP(initial_time, initial_state_vec, initial_guess_input_vec, convergence_radius, max_iteration, initial_control_input_and_constraints_vec);
    model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
    for(int i=0; i<horizon_division_num_; i++){
        control_input_and_constraints_seq_.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_) = initial_control_input_and_constraints_vec;
        state_mat_.col(i) = initial_state_vec;
        lambda_mat_.col(i) = initial_lambda_vec;
    }

    // Intialize the optimality error.
    initial_control_input_and_constraints_error = initializer.getOptimalityErrorVec(initial_time, initial_state_vec, initial_control_input_and_constraints_vec);
    for(int i=0; i<horizon_division_num_; i++){
        control_input_and_constraints_error_seq_.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_) = initial_control_input_and_constraints_error;
    }
    state_error_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num_);
    lambda_error_mat_ = Eigen::MatrixXd::Zero(dim_state_, horizon_division_num_);
}


void MultipleShootingCGMRES::controlUpdate(const double current_time, const double sampling_period, const Eigen::VectorXd& current_state_vec, Eigen::Ref<Eigen::VectorXd> optimal_control_input_vec)
{
    // Predict the incremented state.
    incremented_time_ = current_time + difference_increment_;
    model_.stateFunc(current_time, current_state_vec, control_input_and_constraints_seq_.segment(0, dim_control_input_), dx_vec_);
    incremented_state_vec_ = current_state_vec + difference_increment_*dx_vec_;

    forwardDifferenceGMRES(current_time, current_state_vec, control_input_and_constraints_seq_, control_input_and_constraints_update_seq_);

    // Update state_mat_ and lamdba_mat_ by the difference approximation.
    computeStateAndLambda(incremented_time_, incremented_state_vec_, control_input_and_constraints_seq_+difference_increment_*control_input_and_constraints_update_seq_, (1-difference_increment_*zeta_)*state_error_mat_, (1-difference_increment_*zeta_)*lambda_error_mat_, incremented_state_mat_, incremented_lambda_mat_);
    state_update_mat_ = (incremented_state_mat_ - state_mat_)/difference_increment_;
    lambda_update_mat_ = (incremented_lambda_mat_ - lambda_mat_)/difference_increment_;
    state_mat_ += sampling_period * state_update_mat_;
    lambda_mat_ += sampling_period * lambda_update_mat_;

    // Update control_input_and_constraints_seq_.
    control_input_and_constraints_seq_ += sampling_period * control_input_and_constraints_update_seq_;

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

    return std::sqrt(squared_error);
}


inline void MultipleShootingCGMRES::computeOptimalityErrorforControlInputAndConstraints(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& state_mat, const Eigen::MatrixXd& lambda_mat, Eigen::Ref<Eigen::VectorXd> optimality_for_control_input_and_constraints)
{
    // Set and discretize the horizon.
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
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute optimality error for state.
    model_.stateFunc(time_param, state_vec, control_input_and_constraints_seq.segment(0, dim_control_input_), dx_vec_);
    optimality_for_state.col(0) = state_mat.col(0) - state_vec - delta_tau * dx_vec_;
    double tau = time_param + delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.stateFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_), dx_vec_);
        optimality_for_state.col(i) = state_mat.col(i) - state_mat.col(i-1) - delta_tau * dx_vec_;
    }

    // Compute optimality error for lambda.
    model_.phixFunc(tau, state_mat.col(horizon_division_num_-1), dx_vec_);
    optimality_for_lambda.col(horizon_division_num_-1) = lambda_mat.col(horizon_division_num_-1) - dx_vec_;
    for(int i=horizon_division_num_-1; i>=1; i--, tau-=delta_tau){
        model_.hxFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_), lambda_mat.col(i), dx_vec_);
        optimality_for_lambda.col(i-1) = lambda_mat.col(i-1) - lambda_mat.col(i) - delta_tau * dx_vec_;
    }
}


inline void MultipleShootingCGMRES::computeStateAndLambda(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_and_constraints_seq, const Eigen::MatrixXd& optimality_for_state, const Eigen::MatrixXd& optimality_for_lambda, Eigen::Ref<Eigen::MatrixXd> state_mat, Eigen::Ref<Eigen::MatrixXd> lambda_mat)
{
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute the sequence of state under the error for state.
    model_.stateFunc(time_param, state_vec, control_input_and_constraints_seq.segment(0, dim_control_input_), dx_vec_);
    state_mat.col(0) = state_vec + delta_tau * dx_vec_ + optimality_for_state.col(0);
    double tau = time_param + delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.stateFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_), dx_vec_);
        state_mat.col(i) = state_mat.col(i-1) + delta_tau * dx_vec_ + optimality_for_state.col(i);
    }

    // Compute the sequence of lambda under the error for lambda.
    model_.phixFunc(tau, state_mat.col(horizon_division_num_-1), dx_vec_);
    lambda_mat.col(horizon_division_num_-1) = dx_vec_ + optimality_for_lambda.col(horizon_division_num_-1);
    for(int i=horizon_division_num_-1; i>=1; i--, tau-=delta_tau){
        model_.hxFunc(tau, state_mat.col(i-1), control_input_and_constraints_seq.segment(i*dim_control_input_and_constraints_, dim_control_input_and_constraints_), lambda_mat.col(i), dx_vec_);
        lambda_mat.col(i-1) = lambda_mat.col(i) + delta_tau * dx_vec_ + optimality_for_lambda.col(i-1);
    }
}


void MultipleShootingCGMRES::bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> equation_error_vec)
{
    computeOptimalityErrorforControlInputAndConstraints(time_param, state_vec, current_solution_vec, state_mat_, lambda_mat_, control_input_and_constraints_error_seq_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, control_input_and_constraints_error_seq_1_);

    computeOptimalityErrorforStateAndLambda(time_param, state_vec, current_solution_vec, state_mat_, lambda_mat_, state_error_mat_, lambda_error_mat_);
    computeOptimalityErrorforStateAndLambda(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, state_error_mat_1_, lambda_error_mat_1_);

    computeStateAndLambda(incremented_time_, incremented_state_vec_, current_solution_vec, (1-difference_increment_*zeta_)*state_error_mat_, (1-difference_increment_*zeta_)*lambda_error_mat_, incremented_state_mat_, incremented_lambda_mat_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_3_);

    incremented_control_input_and_constraints_seq_ = current_solution_vec + difference_increment_ * control_input_and_constraints_update_seq_;
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);

    equation_error_vec = - (zeta_ - 1/difference_increment_) * control_input_and_constraints_error_seq_ - control_input_and_constraints_error_seq_3_/difference_increment_ - (control_input_and_constraints_error_seq_2_ - control_input_and_constraints_error_seq_1_)/difference_increment_;
}


void MultipleShootingCGMRES::axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> forward_difference_error_vec)
{
    incremented_control_input_and_constraints_seq_ = current_solution_vec + difference_increment_ * direction_vec;
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);

    forward_difference_error_vec =(control_input_and_constraints_error_seq_2_ - control_input_and_constraints_error_seq_1_)/difference_increment_;
}

inline void MultipleShootingCGMRES::givensRotation(Eigen::Ref<Eigen::VectorXd> column_vec, const int i_column)
{
    double tmp1, tmp2;

    tmp1 = givens_c_vec_(i_column) * column_vec(i_column) - givens_s_vec_(i_column) * column_vec(i_column+1);
    tmp2 = givens_s_vec_(i_column) * column_vec(i_column) + givens_c_vec_(i_column) * column_vec(i_column+1);

    column_vec(i_column) = tmp1;
    column_vec(i_column+1) = tmp2;
}

void MultipleShootingCGMRES::forwardDifferenceGMRES(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> solution_update_vec)
{
    // Initialize vectors for QR factrization by Givens rotation.
    for(int i=0; i<=max_dim_krylov_; i++){
        givens_c_vec_(i) = 0.0;
        givens_s_vec_(i) = 0.0;
        g_vec_(i) = 0.0;
    }

    // Generate the initial basis of the Krylov subspace.
    bFunc(time_param, state_vec, current_solution_vec, b_vec_);
    g_vec_(0) = b_vec_.norm();
    basis_mat_.col(0) = b_vec_ / g_vec_(0);


    // k : the dimension of the Krylov subspace at the current iteration.
    int k;
    for(k=0; k<max_dim_krylov_; k++){
        axFunc(time_param, state_vec, current_solution_vec, basis_mat_.col(k), basis_mat_.col(k+1));
        for(int j=0; j<=k; j++){
            hessenberg_mat_(j,k) = basis_mat_.col(k+1).dot(basis_mat_.col(j));
            basis_mat_.col(k+1) -= hessenberg_mat_(j,k) * basis_mat_.col(j);
        }
        hessenberg_mat_(k+1,k) = basis_mat_.col(k+1).norm();

        if(hessenberg_mat_(k+1,k) != 0){
            basis_mat_.col(k+1) = basis_mat_.col(k+1) / hessenberg_mat_(k+1,k);
        }
        else {
            std::cout << "The modified Gram-Schmidt breakdown at k=" << k << std::endl;
            break;
        }

        // Givens Rotation for QR factrization of the least squares problem.
        for(int j=0; j<k; j++){
            givensRotation(hessenberg_mat_.col(k), j);
        }
        double nu = std::sqrt(hessenberg_mat_(k,k)*hessenberg_mat_(k,k) + hessenberg_mat_(k+1,k)*hessenberg_mat_(k+1,k));
        if(nu != 0) {
            givens_c_vec_(k) = hessenberg_mat_(k,k) / nu;
            givens_s_vec_(k) = - hessenberg_mat_(k+1,k) / nu;
            hessenberg_mat_(k,k) = givens_c_vec_(k) * hessenberg_mat_(k,k) - givens_s_vec_(k) * hessenberg_mat_(k+1,k);
            hessenberg_mat_(k+1,k) = 0;
            givensRotation(g_vec_,k);
        }
        else{
            std::cout << "Lose orthogonality of the basis of the Krylov subspace!!" << std::endl;
        }
    }

    // Solve hessenberg_mat * y = g_vec and obtain y.
    for(int i=k-1; i>=0; i--){
        double tmp=g_vec_(i);
        for(int j=i+1; j<k; j++){
            tmp -= hessenberg_mat_(i,j) * givens_c_vec_(j);
        }
        givens_c_vec_(i) = tmp / hessenberg_mat_(i,i);
    }
    for(int i=0; i<dim_equation_; i++){
        double tmp=0;
        for(int j=0; j<k; j++){
            tmp += basis_mat_(i,j) * givens_c_vec_(j);
        }
        solution_update_vec(i) += tmp;
    }
}