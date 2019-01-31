//
// Save data of the numerical simulation.
//

#ifndef SAVE_SIMULATION_DATA_H
#define SAVE_SIMULATION_DATA_H


#include <fstream>
#include <string>
#include <sys/stat.h>
// #include <direct.h> include this header instead of <sys/stat.h> when you use Windows OS. 
#include <Eigen/Core>


namespace nmpcsim{
    void makeSaveDir(const std::string dir_name);
    void saveData(std::ofstream& state_data, std::ofstream& control_input_data, std::ofstream& error_data, const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_vec, const double optimality_error);
}


#endif