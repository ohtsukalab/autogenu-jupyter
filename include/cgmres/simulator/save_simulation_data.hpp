#ifndef SAVE_SIMULATION_DATA_H
#define SAVE_SIMULATION_DATA_H

#include <iostream>
#include <fstream>

namespace cgmres {

// Saves state_vec, contorl_input_vec, and error_norm to file streams.
void saveData(const int dim_state, const int dim_control_input, 
              std::ofstream& state_data, std::ofstream& control_input_data, 
              std::ofstream& error_data, const double time_param, 
              const double* state_vec, const double* control_input_vec, 
              const double error_norm);

} // namespace cgmres

#endif // SAVE_SIMULATION_DATA_H