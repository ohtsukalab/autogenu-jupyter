//
// Save data of the numerical simulation.
//

#ifndef SAVE_SIMULATION_DATA_H
#define SAVE_SIMULATION_DATA_H


#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>


namespace nmpcsim {
void makeSaveDir(const std::string dir_name);
void saveData(const int dim_state, const int dim_control_input, std::ofstream& state_data, std::ofstream& control_input_data, std::ofstream& error_data, const double time_param, const double* state_vec, const double* control_input_vec, const double optimality_error);
}


#endif