#include "matrixfree_gmres.hpp"


inline void MatrixFreeGMRES::givensRotation(double* column_vec, const int i_column)
{
    double tmp1, tmp2;

    tmp1 = givens_c_vec_[i_column] * column_vec[i_column] - givens_s_vec_[i_column] * column_vec[i_column+1];
    tmp2 = givens_s_vec_[i_column] * column_vec[i_column] + givens_c_vec_[i_column] * column_vec[i_column+1];

    column_vec[i_column] = tmp1;
    column_vec[i_column+1] = tmp2;
}


MatrixFreeGMRES::MatrixFreeGMRES() : 
    dim_equation_(0), 
    max_dim_krylov_(0), 
    hessenberg_mat_(nullptr), 
    basis_mat_(nullptr), 
    b_vec_(nullptr), 
    givens_c_vec_(nullptr), 
    givens_s_vec_(nullptr), 
    g_vec_(nullptr)
{}


MatrixFreeGMRES::MatrixFreeGMRES(const int dim_equation, const int max_dim_krylov) : 
    dim_equation_(dim_equation), 
    max_dim_krylov_(max_dim_krylov), 
    hessenberg_mat_(linearfunc::newMatrix(max_dim_krylov_+1, max_dim_krylov_+1)), 
    basis_mat_(linearfunc::newMatrix(max_dim_krylov_+1, dim_equation)), 
    b_vec_(linearfunc::newVector(dim_equation)), 
    givens_c_vec_(linearfunc::newVector(max_dim_krylov+1)), 
    givens_s_vec_(linearfunc::newVector(max_dim_krylov+1)), 
    g_vec_(linearfunc::newVector(max_dim_krylov+1))
{}


MatrixFreeGMRES::~MatrixFreeGMRES()
{
    linearfunc::deleteMatrix(hessenberg_mat_);
    linearfunc::deleteMatrix(basis_mat_);
    linearfunc::deleteVector(b_vec_);
    linearfunc::deleteVector(givens_c_vec_);
    linearfunc::deleteVector(givens_s_vec_);
    linearfunc::deleteVector(g_vec_);
}


void MatrixFreeGMRES::setGMRESParams(const int dim_equation, const int max_dim_krylov)
{
    dim_equation_ = dim_equation;
    max_dim_krylov_ = max_dim_krylov;

    hessenberg_mat_ = linearfunc::newMatrix(max_dim_krylov+1, max_dim_krylov+1);
    basis_mat_ = linearfunc::newMatrix(max_dim_krylov+1, dim_equation);
    b_vec_ = linearfunc::newVector(dim_equation);
    givens_c_vec_ = linearfunc::newVector(max_dim_krylov+1);
    givens_s_vec_ = linearfunc::newVector(max_dim_krylov+1);
    g_vec_ = linearfunc::newVector(max_dim_krylov+1);
}


void MatrixFreeGMRES::forwardDifferenceGMRES(const double time_param, const double* state_vec, const double* current_solution_vec, double* solution_update_vec)
{
    // Initialize vectors for QR factrization by Givens rotation.
    // Set givens_c_vec_, givens_s_vec_, g_vec_ as zero.
    for(int i=0; i<max_dim_krylov_+1; i++){
        givens_c_vec_[i] = 0;
        givens_s_vec_[i] = 0;
        g_vec_[i] = 0;
    }

    // Generate the initial basis of the Krylov subspace.
    bFunc(time_param, state_vec, current_solution_vec, b_vec_);
    g_vec_[0] = std::sqrt(linearfunc::squaredNorm(dim_equation_, b_vec_));
    for(int i=0; i<dim_equation_; i++){
        basis_mat_[0][i] = b_vec_[i] / g_vec_[0];
    }

    // : the dimension of the Krylov subspace at the current iteration.
    int k;
    for(k=0; k<max_dim_krylov_; k++){
        axFunc(time_param, state_vec, current_solution_vec, basis_mat_[k], basis_mat_[k+1]);
        for(int j=0; j<=k; j++){
            hessenberg_mat_[k,j] = linearfunc::innerProduct(dim_equation_, basis_mat_[k+1], basis_mat_[j]);
            // basis_mat_.col(k+1) -= hessenberg_mat_(j,k) * basis_mat_.col(j);
            for(int i=0; i<dim_equation_; i++){
                basis_mat_[k+1][i] -= hessenberg_mat_[k,j] * basis_mat_[j][i];
            }
        }
        hessenberg_mat_[k,k+1] = std::sqrt(linearfunc::squaredNorm(dim_equation_, basis_mat_[k+1]));

        if(hessenberg_mat_[k,k+1] != 0){
            // basis_mat_.col(k+1) = basis_mat_.col(k+1) / hessenberg_mat_(k+1,k);
            for(int i=0; i<dim_equation_; i++){
                basis_mat_[k+1][i] = basis_mat_[k+1][i] / hessenberg_mat_[k,k+1];
            }
        }
        else {
            std::cout << "The modified Gram-Schmidt breakdown at k=" << k << std::endl;
            break;
        }

        // Givens Rotation for QR factrization of the least squares problem.
        for(int j=0; j<k; j++){
            givensRotation(hessenberg_mat_[k], j);
        }
        double nu = std::sqrt(hessenberg_mat_[k,k]*hessenberg_mat_[k,k] + hessenberg_mat_[k,k+1]*hessenberg_mat_[k,k+1]);
        if(nu != 0) {
            givens_c_vec_[k] = hessenberg_mat_[k,k] / nu;
            givens_s_vec_[k] = - hessenberg_mat_[k,k+1] / nu;
            hessenberg_mat_[k,k] = givens_c_vec_[k] * hessenberg_mat_[k,k] - givens_s_vec_[k] * hessenberg_mat_[k,k+1];
            hessenberg_mat_[k,k+1] = 0;
            givensRotation(g_vec_,k);
        }
        else{
            std::cout << "Lose orthogonality of the basis of the Krylov subspace!!" << std::endl;
        }
    }

    // Solve hessenberg_mat * y = g_vec and obtain y.
    for(int i=k-1; i>=0; i--){
        double tmp=g_vec_[i];
        for(int j=i+1; j<k; j++){
            tmp -= hessenberg_mat_[j,i] * givens_c_vec_[j];
        }
        givens_c_vec_[i] = tmp / hessenberg_mat_[i,i];
    }
    for(int i=0; i<dim_equation_; i++){
        double tmp=0;
        for(int j=0; j<k; j++){
            tmp += basis_mat_[j,i] * givens_c_vec_[j];
        }
        solution_update_vec[i] += tmp;
    }
}