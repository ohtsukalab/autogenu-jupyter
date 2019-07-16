#ifndef CONDENSING_FUNCTIONS_H 
#define CONDENSING_FUNCTIONS_H 

#include "control_input_saturation_sequence.hpp"

namespace condensingfunctions {
// Computes the partial derivative of the product of constraints on the control
// input saturation that is condensed and the corresponding Lagrange multiplier
// with respect to the control input vector and adds it to a given errors in 
// optimality. That derivative is added to 
// errors_for_control_input_and_constraints_vec.
void addHamiltonianDerivativeWithConstrainedInput(
    ControlInputSaturationSequence& saturation_seq,
    const double* control_input_and_constraints_vec, 
    const double* saturation_lagrange_multiplier_vec, 
    double* errors_for_control_input_and_constraints_vec);

// Computes the partial derivative of the Hamiltonian with respect to the 
// dummy input.
void computeErrorsForDummyInput(
    ControlInputSaturationSequence& saturation_seq,
    const double* dummy_input_vec, 
    const double* saturation_lagrange_multiplier_vec, 
    double* errors_for_dummy_input);

// Computes the errors in optimality of the condensed constraints on
// the control input saturation function.
void computeErrorsForSaturation(
    ControlInputSaturationSequence& saturation_seq,
    const double* control_input_and_constraint_vec, 
    const double* dummy_input_vec, double* errors_for_saturation);
}

#endif // CONDENSING_FUNCTIONS_H 