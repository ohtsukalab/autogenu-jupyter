#include "save_simulation_data_for_windows.hpp"

void nmpcsim::makeSaveDir(const std::string dir_name)
{
    int error = mkdir(dir_name.c_str());
}

void nmpcsim::saveData(const int dim_state, const int dim_control_input, std::ofstream& state_data, std::ofstream& control_input_data, std::ofstream& error_data, const double time_param, const double* state_vec, const double* control_input_vec, const double optimality_error)
{
    for(int i=0; i<dim_state; i++){
        state_data << state_vec[i] << " ";
    }
    state_data << "\n";

    for(int i=0; i<dim_control_input; i++){
        control_input_data << control_input_vec[i] << " ";
    }
    control_input_data << "\n";

    error_data << optimality_error << "\n";
}