//
// Simulator of the C/GMRES method and multiple shooting based C/GMRES method.
//

#ifndef SIMULATOR_H
#define SIMULATOR_H


#include <iostream>
#include <fstream>
#include <chrono>
#include <eigen3/Eigen/Core>
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "continuation_gmres.hpp"
#include "multiple_shooting_cgmres.hpp"


// Simulates the C/GMRES method and multiple shooting based C/GMRES method.
// Opens file streams and saves simulation data.
// Updates the state using the method in NumericalIntegrator.
class Simulator final : virtual public NumericalIntegrator{
private:
    NMPCModel model_;
    
    void saveData(std::ofstream& state_data, std::ofstream& control_input_data, std::ofstream& error_data, const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_vec, const double optimality_error);
public:
    Simulator(const NMPCModel model);
    // Serves a numerical simulation for the C/GMRES and the multiple shooting based C/GMRES.
    void simulation(ContinuationGMRES cgmres_solver, const Eigen::VectorXd& initial_state_vec, const double start_time, const double end_time, const double sampling_period, const std::string savefile_name);
    void simulation(MultipleShootingCGMRES cgmres_solver, const Eigen::VectorXd& initial_state_vec, const double start_time, const double end_time, const double sampling_period, const std::string savefile_name);
};


#endif