namespace linearfunc{
    double* newVector(const int dim);
    void deleteVector(double* vec);
    double** newMatrix(const int dim_row, const int dim_column);
    void deleteMatrix(double** mat);
    double innerProduct(const int dim, const double *vec1, const double *vec2);
    double squaredNorm(const int dim, const double *vec);
}

