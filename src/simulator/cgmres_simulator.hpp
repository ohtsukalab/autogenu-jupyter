//
// Simulator of the C/GMRES method.
//

#ifndef CGMRES_SIMULATOR_H
#define CGMRES_SIMULATOR_H


#include <iostream>
#include <fstream>
#include <chrono>
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "save_simulation_data.hpp"
#include "continuation_gmres.hpp"


// Simulates NMPC using the C/GMRES method and multiple shooting based C/GMRES method.
// Opens file streams and saves simulation data.
// Updates the state using the method in NumericalIntegrator.
namespace nmpcsim{
    // Perform a numerical simulation for the C/GMRES method.
    void simulation(ContinuationGMRES nmpc_solver, const double* initial_state_vec, const double start_time, const double end_time, const double sampling_period, const std::string savefile_name);
};


#endif