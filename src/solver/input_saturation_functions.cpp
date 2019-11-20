#include "input_saturation_functions.hpp"

void inputsaturationfunctions::addHamiltonianDerivativeWithSaturatedInput(
    InputSaturationSet& input_saturation_set,
    const double* control_input_and_constraints_vec, 
    const double* input_saturation_multiplier_vec, 
    double* optimality_residual_for_control_input_and_constraints_vec) {
  for (int i=0; i<input_saturation_set.dimSaturation(); ++i) {
    int index_i = input_saturation_set.index(i);
    optimality_residual_for_control_input_and_constraints_vec[index_i] +=
        (2*control_input_and_constraints_vec[index_i] 
            -input_saturation_set.min(i)-input_saturation_set.max(i)) 
        * input_saturation_multiplier_vec[i];
  }
}

void inputsaturationfunctions::computeOptimalityResidualForDummyInput(
    InputSaturationSet& input_saturation_set,
    const double* dummy_input_vec, 
    const double* input_saturation_multiplier_vec, 
    double* optimality_residual_for_dummy_input) {
 for (int i=0; i<input_saturation_set.dimSaturation(); ++i) {
    optimality_residual_for_dummy_input[i] 
        = 2 * (input_saturation_set.quadratic_weight(i)
                  +input_saturation_multiplier_vec[i]) 
            * dummy_input_vec[i] - input_saturation_set.dummy_weight(i);
  }
}

void inputsaturationfunctions::computeOptimalityResidualForInputSaturation(
    InputSaturationSet& input_saturation_set,
    const double* control_input_and_constraint_vec, 
    const double* dummy_input_vec, double* optimality_residual_for_saturation) {
  for (int i=0; i<input_saturation_set.dimSaturation(); ++i) {
    int index_i = input_saturation_set.index(i);
    double min_i = input_saturation_set.min(i);
    double max_i = input_saturation_set.max(i);
    optimality_residual_for_saturation[i] = 
        control_input_and_constraint_vec[index_i] * 
        (control_input_and_constraint_vec[index_i]-min_i-max_i)
        + min_i * max_i + dummy_input_vec[i] * dummy_input_vec[i];
  }
}