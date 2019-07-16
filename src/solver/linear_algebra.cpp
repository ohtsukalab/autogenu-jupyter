#include "linear_algebra.hpp"

double* linearalgebra::NewVector(const int dim) {
  double* vec = new double[dim];
  for (int i=0; i<dim; ++i) {
    vec[i] = 0;
  }
  return vec;
}

void linearalgebra::DeleteVector(double* vec) {
  delete[] vec;
}

double** linearalgebra::NewMatrix(const int dim_row, const int dim_column) {
  double** mat = new double*[dim_row];
  mat[0] = new double[dim_row*dim_column];
  for (int i=1; i<dim_row; ++i) {
    mat[i] = mat[i-1] + dim_column;
  }
  for (int i=0; i<dim_row*dim_column; ++i) {
    mat[0][i] = 0;
  }
  return mat;
}

void linearalgebra::DeleteMatrix(double** mat) {
  delete[] mat[0];
  delete[] mat;
}

double linearalgebra::InnerProduct(const int dim, const double *vec1, 
                                   const double *vec2) {
  double ans = 0;
  for (int i=0; i<dim; ++i) {
    ans += vec1[i] * vec2[i];
  }
  return ans;
}


double linearalgebra::SquaredNorm(const int dim, const double *vec) {
  double ans = 0;
  for (int i=0; i<dim; ++i) {
    ans += vec[i] * vec[i];
  }
  return ans;
}