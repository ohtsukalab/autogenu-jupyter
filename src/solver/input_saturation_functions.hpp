#ifndef INPUT_SATURATION_FUNCTIONS_H 
#define INPUT_SATURATION_FUNCTIONS_H 

#include "input_saturation_set.hpp"

namespace inputsaturationfunctions {
// Computes the partial derivative of the product of constraints on the control
// input saturation that is condensed and the corresponding Lagrange multiplier
// with respect to the control input vector and adds it to a given errors in 
// optimality. That derivative is added to 
// errors_for_control_input_and_constraints_vec.
void addHamiltonianDerivativeWithSaturatedInput(
    InputSaturationSet& input_saturation_set,
    const double* control_input_and_constraints_vec, 
    const double* input_saturation_multiplier_vec, 
    double* optimality_residual_for_control_input_and_constraints_vec);

// Computes the partial derivative of the Hamiltonian with respect to the 
// dummy input.
void computeOptimalityResidualForDummyInput(
    InputSaturationSet& input_saturation_set,
    const double* dummy_input_vec, 
    const double* input_saturation_multiplier_vec, 
    double* optimality_residual_for_dummy_input);

// Computes the errors in optimality of the condensed constraints on
// the control input saturation function.
void computeOptimalityResidualForInputSaturation(
    InputSaturationSet& input_saturation_set,
    const double* control_input_and_constraint_vec, 
    const double* dummy_input_vec, double* optimality_residual_for_saturation);
}

#endif // INPUT_SATURATION_FUNCTIONS_H