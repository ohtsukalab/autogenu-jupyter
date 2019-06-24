// 
// Functions for linear arithmetic.
// 

#ifndef LINEAR_FUNCS_H
#define LINEAR_FUNCS_H


namespace linearfunc{
// Allocates memory for a vector and set all components zero.
double* newVector(const int dim);

// Free memory of a vector.
void deleteVector(double* vec);

// Allocates memory for a matrix and set all components zero.    
double** newMatrix(const int dim_row, const int dim_column);

//  Free memory of a matrix.
void deleteMatrix(double** mat);

// Returns inner product of two vectors.
double innerProduct(const int dim, const double *vec1, const double *vec2);

// Returns squared norm of a vector.
double squaredNorm(const int dim, const double *vec);
}

#endif