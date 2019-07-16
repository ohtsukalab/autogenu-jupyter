#include "condensing_functions.hpp"


void condensingfunctions::addHamiltonianDerivativeWithConstrainedInput(
    ControlInputSaturationSequence& saturation_seq,
    const double* control_input_and_constraints_vec, 
    const double* saturation_lagrange_multiplier_vec, 
    double* errors_for_control_input_and_constraints_vec) {
  for (int i=0; i<saturation_seq.dimSaturation(); ++i) {
    int index_i = saturation_seq.index(i);
    errors_for_control_input_and_constraints_vec[index_i] +=
        (2*control_input_and_constraints_vec[index_i] 
            -saturation_seq.min(i)-saturation_seq.max(i)) 
        * saturation_lagrange_multiplier_vec[i];
  }
}

void condensingfunctions::computeErrorsForDummyInput(
    ControlInputSaturationSequence& saturation_seq,
    const double* dummy_input_vec, 
    const double* saturation_lagrange_multiplier_vec, 
    double* errors_for_dummy_input) {
 for (int i=0; i<saturation_seq.dimSaturation(); ++i) {
    errors_for_dummy_input[i] = 
        2 * (saturation_seq.quadratic_weight(i)
            +saturation_lagrange_multiplier_vec[i]) 
        * dummy_input_vec[i] - saturation_seq.dummy_weight(i);
  }
}

void condensingfunctions::computeErrorsForSaturation(
    ControlInputSaturationSequence& saturation_seq,
    const double* control_input_and_constraint_vec, 
    const double* dummy_input_vec, double* errors_for_saturation) {
  for (int i=0; i<saturation_seq.dimSaturation(); ++i) {
    int index_i = saturation_seq.index(i);
    double min_i = saturation_seq.min(i);
    double max_i = saturation_seq.max(i);
    errors_for_saturation[i] = 
        control_input_and_constraint_vec[index_i] * 
        (control_input_and_constraint_vec[index_i]-min_i-max_i)
        + min_i * max_i + dummy_input_vec[i] * dummy_input_vec[i];
  }
}