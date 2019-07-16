#ifndef LINEAR_ALGEBRA_H 
#define LINEAR_ALGEBRA_H 

// Functions supporting linear algebra. 
namespace linearalgebra {
// Allocates memory for a vector whose dimension is dim and set all components 
// zero. Then returns the pointer to the vector.
double* NewVector(const int dim);

// Free memory of a vector. 
void DeleteVector(double* vec);

// Allocates memory for a matrix whose dimensions are given by dim_row and 
// dim_column and set all components zero. Then returns the pointer to the 
// matrix.
double** NewMatrix(const int dim_row, const int dim_column);

//  Free memory of a matrix.
void DeleteMatrix(double** mat);

// Returns inner product of vec_1 and vec_2.
double InnerProduct(const int dim, const double *vec1, const double *vec2);

// Returns squared norm of vec.
double SquaredNorm(const int dim, const double *vec);
}

#endif // LINEAR_ALGEBRA