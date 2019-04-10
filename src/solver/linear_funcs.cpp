#include "linear_funcs.hpp"



double* linearfunc::newVector(const int dim)
{
    double* vec = new double[dim];
    for(int i=0; i<dim; i++){
        vec[i] = 0;
    }
    return vec;
}


void linearfunc::deleteVector(double* vec)
{
    delete[] vec;
}


double** linearfunc::newMatrix(const int dim_column, const int dim_row)
{
    double** mat = new double*[dim_row];
    mat[0] = new double[dim_row*dim_column];
    for(int i=1; i<dim_row; i++){
        mat[i] = mat[i-1] + dim_column;
    }
    for(int i=0; i<dim_row*dim_column; i++){
        mat[0][i] = 0;
    }
    return mat;
}


void linearfunc::deleteMatrix(double** mat)
{
    delete[] mat[0];
    delete[] mat;
}



void linearfunc::substituteVec(const int dim, const double* origin_vec, double* vec)
{
    for(int i=0; i<dim; i++){
        vec[i] = origin_vec[i];
    }
}


void linearfunc::increaseVec(const int dim, const double increase_scalar, const double* increase_vec, double* vec)
{
    for(int i=0; i<dim; i++){
        vec[i] += increase_scalar * increase_vec[i];
    }
}


void linearfunc::sumVec(const int dim, const double* increse_vec, double* vec)
{
    for(int i=0; i<dim; i++){
        vec[i] += increse_vec[i];
    }
}


void linearfunc::sumVec(const int dim, const double* vec1, const double* vec2, double* result_vec)
{
    for(int i=0; i<dim; i++){
        result_vec[i] = vec1[i] + vec2[i];
    }
}


void linearfunc::sumMat(const int dim_low, const int dim_column, const double** mat1, const double** mat2, double** result_mat)
{
    for(int i=0; i<dim_column; i++){
        for(int j=0; j<dim_low; j++){
            result_mat[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
}


void linearfunc::decreaseVec(const int dim, const double decrease_scalar, const double* decrease_vec, double* vec)
{
    for(int i=0; i<dim; i++){
        vec[i] -= decrease_scalar * decrease_vec[i];
    }
}

void linearfunc::subtractVec(const int dim, const double* vec1, const double* vec2, double* result_vec)
{
    for(int i=0; i<dim; i++){
        result_vec[i] = vec1[i] - vec2[i];
    }
}



void linearfunc::multipleVec(const int dim, const double scalar, double* vec)
{
    for(int i=0; i<dim; i++){
        vec[i] = scalar * vec[i];
    }
}


void linearfunc::multipleVec(const int dim, const double scalar, const double* vec, double* result_vec)
{
    for(int i=0; i<dim; i++){
        result_vec[i] = scalar * vec[i];
    }
}


void linearfunc::multipleMat(const int dim_low, const int dim_column, const double scalar, const double** mat, double** result_mat)
{
    for(int i=0; i<dim_column; i++){
        for(int j=0; j < dim_low; j++){
            result_mat[i][j] = scalar * mat[i][j];
        }
    }
}


void linearfunc::divisionVec(const int dim, const double scalar, double* vec)
{
    for(int i=0; i<dim; i++){
        vec[i] = vec[i]/scalar;
    }
}


void linearfunc::divisionVec(const int dim, const double scalar, const double* vec, double* result_vec)
{
    for(int i=0; i<dim; i++){
        result_vec[i] = vec[i]/scalar;
    }
}


void linearfunc::divisionMat(const int dim_low, const int dim_column, const double scalar, const double** mat, double** result_mat)
{
    for(int i=0; i<dim_column; i++){
        for(int j=0; j < dim_low; j++){
            result_mat[i][j] = mat[i][j]/scalar;
        }
    }
}




double linearfunc::innerProduct(const int dim, const double *vec1, const double *vec2)
{
    double ans = 0;
    for(int i=0; i<dim; i++){
        ans += vec1[i] * vec2[i];
    }
    return ans;
}


double linearfunc::squaredNorm(const int dim, const double *vec)
{
    double ans = 0;
    for(int i = 0; i < dim; i++){
        ans += vec[i] * vec[i];
    }
}


void linearfunc::dotProduct(const int dim_vec, const int dim_column, const int dim_low, const double *vec, const double **mat, double *result_vec)
{
    for(int i=0; i<dim_column; i++){
        result_vec[i] = 0;
        for(int j=0; j<dim_low; j++){
            result_vec[i] += mat[i][j] * vec[j];
        }
    }
}


void linearfunc::integrator::(const int dim, const double* initial_vec, const double integration_step, const double* derivative_vector, double* result_vec)
{
    for(int i=0; i<dim, i++){
        result_vec[i] = initial_vec[i] + integration_step * derivative_vector[i];
    }
}