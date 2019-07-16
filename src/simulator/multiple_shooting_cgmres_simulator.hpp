#ifndef MULTIPLE_SHOOTING_CGMRES_SIMULATOR_H
#define MULTIPLE_SHOOTING_CGMRES_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "save_simulation_data.hpp"
#include "multiple_shooting_cgmres.hpp"

namespace nmpcsim {
// Simulates NMPC using the multiple shooting-based C/GMRES method. Opens file 
// streams and saves simulation data into them. 
void simulation(MultipleShootingCGMRES& nmpc_solver, 
                const double* initial_state_vec, const double start_time, 
                const double end_time, const double sampling_period, 
                const std::string save_dir, const std::string savefile_name);
}

#endif // MULTIPLE_SHOOTING_CGMRES_SIMULATOR_H