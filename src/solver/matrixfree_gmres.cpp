#include "matrixfree_gmres.hpp"


MatrixFreeGMRES::MatrixFreeGMRES()
  : allocation_flag_(false), 
    dim_solution_(0), 
    kmax_(0), 
    hessenberg_mat_(nullptr), 
    basis_mat_(nullptr), 
    b_vec_(nullptr), 
    givens_c_vec_(nullptr), 
    givens_s_vec_(nullptr), 
    g_vec_(nullptr) {
}

MatrixFreeGMRES::MatrixFreeGMRES(const int dim_solution, const int kmax)
  : allocation_flag_(true), 
    dim_solution_(dim_solution), 
    kmax_(kmax), 
    hessenberg_mat_(linearalgebra::NewMatrix(kmax+1, kmax+1)), 
    basis_mat_(linearalgebra::NewMatrix(kmax+1, dim_solution)), 
    b_vec_(linearalgebra::NewVector(dim_solution)), 
    givens_c_vec_(linearalgebra::NewVector(kmax+1)), 
    givens_s_vec_(linearalgebra::NewVector(kmax+1)), 
    g_vec_(linearalgebra::NewVector(kmax+1)) {
}

MatrixFreeGMRES::~MatrixFreeGMRES() {
  if (allocation_flag_) {
    // Frees vectors and matrices only if the vectors and matrices are 
    // allocated.
    linearalgebra::DeleteMatrix(hessenberg_mat_);
    linearalgebra::DeleteMatrix(basis_mat_);
    linearalgebra::DeleteVector(b_vec_);
    linearalgebra::DeleteVector(givens_c_vec_);
    linearalgebra::DeleteVector(givens_s_vec_);
    linearalgebra::DeleteVector(g_vec_);
  }
}

void MatrixFreeGMRES::setGMRESParameters(const int dim_solution, 
                                         const int kmax) {
  dim_solution_ = dim_solution;
  kmax_ = kmax;
  if (allocation_flag_) {
    // Frees vectors and matrices for reallocation only if the vectors and 
    // matrices are allocated.
    linearalgebra::DeleteMatrix(hessenberg_mat_);
    linearalgebra::DeleteMatrix(basis_mat_);
    linearalgebra::DeleteVector(b_vec_);
    linearalgebra::DeleteVector(givens_c_vec_);
    linearalgebra::DeleteVector(givens_s_vec_);
    linearalgebra::DeleteVector(g_vec_);
  }
  allocation_flag_ = true;
  hessenberg_mat_ = linearalgebra::NewMatrix(kmax+1, kmax+1);
  basis_mat_ = linearalgebra::NewMatrix(kmax+1, dim_solution);
  b_vec_ = linearalgebra::NewVector(dim_solution);
  givens_c_vec_ = linearalgebra::NewVector(kmax+1);
  givens_s_vec_ = linearalgebra::NewVector(kmax+1);
  g_vec_ = linearalgebra::NewVector(kmax+1);
}

void MatrixFreeGMRES::solveGMRES(const double time, 
                                 const double* state_vec, 
                                 const double* current_solution_vec, 
                                 double* solution_update_vec) {
  // Initializes vectors for QR factrization by Givens rotation.
  // Set givens_c_vec_, givens_s_vec_, g_vec_ as zero.
  for (int i=0; i<kmax_+1; ++i) {
    givens_c_vec_[i] = 0;
    givens_s_vec_[i] = 0;
    g_vec_[i] = 0;
  }
  // Generates the initial basis of the Krylov subspace.
  bFunc(time, state_vec, current_solution_vec, b_vec_);
  g_vec_[0] = std::sqrt(linearalgebra::SquaredNorm(dim_solution_, b_vec_));
  // basis_mat_[0] = b_vec_ / g_vec[0]
  for (int i=0; i<dim_solution_; ++i) {
    basis_mat_[0][i] = b_vec_[i] / g_vec_[0];
  }
  // k : the dimension of the Krylov subspace at the current iteration.
  int k;
  for (k=0; k<kmax_; ++k) {
    axFunc(time, state_vec, current_solution_vec, basis_mat_[k], 
           basis_mat_[k+1]);
    for (int j=0; j<=k; ++j) {
      hessenberg_mat_[k][j] = linearalgebra::InnerProduct(
          dim_solution_, basis_mat_[k+1], basis_mat_[j]);
      // basis_mat_[k+1] -= hessenberg_mat_[k][j] * basis_mat_[j];
      for (int i=0; i<dim_solution_; ++i) {
        basis_mat_[k+1][i] -= hessenberg_mat_[k][j] * basis_mat_[j][i];
      }
    }
    hessenberg_mat_[k][k+1] = std::sqrt(linearalgebra::SquaredNorm(
        dim_solution_, basis_mat_[k+1]));
    if (hessenberg_mat_[k][k+1]) {
      // basis_mat_[k+1] = basis_mat_[k+1] / hessenberg_mat_[k][k+1];
      for (int i=0; i<dim_solution_; ++i) {
        basis_mat_[k+1][i] = basis_mat_[k+1][i] / hessenberg_mat_[k][k+1];
      }
    }
    else {
      std::cout << "The modified Gram-Schmidt breakdown at k = " << k 
                << std::endl;
      break;
    }
    // Givens Rotation for QR factrization of the least squares problem.
    for (int j=0; j<k; ++j) {
      givensRotation(hessenberg_mat_[k], j);
    }
    double nu = std::sqrt(hessenberg_mat_[k][k]*hessenberg_mat_[k][k]
                          +hessenberg_mat_[k][k+1]*hessenberg_mat_[k][k+1]);
    if (nu) {
      givens_c_vec_[k] = hessenberg_mat_[k][k] / nu;
      givens_s_vec_[k] = - hessenberg_mat_[k][k+1] / nu;
      hessenberg_mat_[k][k] = givens_c_vec_[k] * hessenberg_mat_[k][k] 
                              - givens_s_vec_[k] * hessenberg_mat_[k][k+1];
      hessenberg_mat_[k][k+1] = 0;
      givensRotation(g_vec_,k);
    }
    else {
      std::cout << "Lose orthogonality of the basis of the Krylov subspace" 
                << std::endl;
    }
  }
  // Computes solution_update_vec by solving hessenberg_mat_ * y = g_vec.
  for (int i=k-1; i>=0; --i) {
    double tmp = g_vec_[i];
    for (int j=i+1; j<k; ++j) {
      tmp -= hessenberg_mat_[j][i] * givens_c_vec_[j];
    }
    givens_c_vec_[i] = tmp / hessenberg_mat_[i][i];
  }
  for (int i=0; i<dim_solution_; ++i) {
    double tmp = 0;
    for (int j=0; j<k; ++j) { 
      tmp += basis_mat_[j][i] * givens_c_vec_[j];
    }
    solution_update_vec[i] += tmp;
  }
}

inline void MatrixFreeGMRES::givensRotation(double* column_vec, 
                                            const int i_column) {
  double tmp1 = givens_c_vec_[i_column] * column_vec[i_column] 
                - givens_s_vec_[i_column] * column_vec[i_column+1];
  double tmp2 = givens_s_vec_[i_column] * column_vec[i_column] 
                + givens_c_vec_[i_column] * column_vec[i_column+1];
  column_vec[i_column] = tmp1;
  column_vec[i_column+1] = tmp2;
}