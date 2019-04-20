//
// Computes the initial solution of the multiple shooting based C/GMRES method with condensing of the saturations on the control input using the Newton GMRES method.
//

#ifndef INIT_CGMRES_WITH_SATURATION_H
#define INIT_CGMRES_WITH_SATURATION_H


#include "linear_funcs.hpp"
#include "matrixfree_gmres.hpp"
#include "control_input_saturation_sequence.hpp"


// Computes the initial solution of the multiple shooting based C/GMRES method with condensing of the saturations on the control input using the Newton GMRES method.
class InitCGMRESWithSaturation final : public MatrixFreeGMRES{
private:
    NMPCModel model_;
    ControlInputSaturationSequence saturation_seq_;

    int dim_control_input_and_constraints_, dim_saturation_, dim_solution_;
    double difference_increment_;
    double *incremented_solution_vec_, *solution_update_vec_, *lambda_vec_, *error_vec_, *error_vec_1_, *error_vec_2_;


    // Adds partial derivative of the saturation with respect to the control input
    inline void addHamiltonianDerivativeWithControlInput(const double* control_input_and_constraints_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_control_input_and_constraints_vec);

    // Computes the optimality for the dummy input
    inline void computeDummyOptimality(const double* dummy_input_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_dummy);

    // Computes the optimality of the saturation
    inline void computeSaturationOptimality(const double* control_input_and_constraint_vec, const double* dummy_input_vec, double* optimality_for_saturation);


    // Computes the optimality error vector under current_solution_vec.
    inline void computeOptimalityErrors(const double time_param, const double* state_vec, const double* solution_vec, double* optimality_vec);

    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec) override;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec) override;


public:
    // Sets parameters and allocates vectors.
    InitCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, const double difference_increment, const int max_dim_krylov);

    // Free vectors/
    ~InitCGMRESWithSaturation();


    // Calls the forwardDifferenceGMRES, solves the GMRES, and obtains the solution 
    void solve0stepNOCP(const double initial_time, const double* initial_state_vec, const double* initial_guess_control_input_vec, const double* initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration, double* solution_vec);

    // Returns the optimality error for the control input and constraints under current_solution_vec.
    void getControlInputAndConstraintsError(const double initial_time, const double* initial_state_vec, const double* current_solution_vec, double* control_input_and_constraints_error_vec);

    // Returns the optimality error for dummy input under current_solution_vec.
    void getDummyInputError(const double initial_time, const double* initial_state_vec, const double* current_solution_vec, double* dummy_input_error_vec);

    // Returns the optimality error for saturation under current_solution_vec.
    void getControlInputSaturationError(const double initial_time, const double* initial_state_vec, const double* current_solution_vec, double* control_input_saturation_error_vec);
};


#endif