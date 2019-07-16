#ifndef MSCGMRES_WITH_SATURATION_SIMULATOR_H
#define MSCGMRES_WITH_SATURATION_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "save_simulation_data.hpp"
#include "mscgmres_with_saturation.hpp"

namespace nmpcsim {
// Simulates NMPC using the multiple shooting-based C/GMRES method. Opens file 
// streams and saves simulation data into them. 
void simulation(MSCGMRESWithSaturation& nmpc_solver, 
                const double* initial_state_vec, const double start_time, 
                const double end_time, const double sampling_period, 
                const std::string save_dir, const std::string savefile_name);
}

#endif // MSCGMRES_WITH_SATURATION_SIMULATOR_H