namespace linearfunc{
    double* newVector(const int dim);
    void deleteVector(double* vec);
    double** newMatrix(const int dim_row, const int dim_column);
    void deleteMatrix(double** mat);
    void substituteVec(const int dim, const double* origin_vec, double* vec);
    void increaseVec(const int dim, const double increase_scalar, const double* increase_vec, double* vec);
    void sumVec(const int dim, const double* increse_vec, double* vec);
    void sumVec(const int dim, const double* vec1, const double* vec2, double* result_vec);
    void sumMat(const int dim_low, const int dim_column, const double** mat1, const double** mat2, double** result_mat);
    void decreaseVec(const int dim, const double decrease_scalar, const double* decrease_vec, double* vec);
    void subtractVec(const int dim, const double* vec1, const double* vec2, double* result_vec);
    void multipleVec(const int dim, const double scalar, double* vec);
    void multipleVec(const int dim, const double scalar, const double* vec, double* result_vec);
    void multipleMat(const int dim_low, const int dim_column, const double scalar, const double** mat, double** result_mat);
    void divisionVec(const int dim, const double scalar, double* vec);
    void divisionVec(const int dim, const double scalar, const double* vec, double* result_vec);
    void divisionMat(const int dim_low, const int dim_column, const double scalar, const double** mat, double** result_mat);
    double innerProduct(const int dim, const double *vec1, const double *vec2);
    double squaredNorm(const int dim, const double *vec);
    void dotProduct(const int dim_vec, const int dim_column, const int dim_low, const double* vec, const double **mat, double *result_vec);
    namespace integrator{
        eulerMethod(const int dim, const double* initial_vec, const double integration_step, const double* derivative_vector, double* result_vec);
    }
}

