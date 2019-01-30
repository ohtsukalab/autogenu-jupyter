#include "save_simulation_data.hpp"


void nmpcsim::makeSaveDir(const std::string dir_name)
{
    int error = mkdir(dir_name.c_str(), 0755);
}


void nmpcsim::saveData(std::ofstream& state_data, std::ofstream& control_input_data, std::ofstream& error_data, const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_vec, const double optimality_error)
{
    for(int i=0; i<state_vec.size(); i++){
        state_data << state_vec(i) << " ";
    }
    state_data << "\n";

    for(int i=0; i<control_input_vec.size(); i++){
        control_input_data << control_input_vec(i) << " ";
    }
    control_input_data << "\n";

    error_data << optimality_error << "\n";
}
