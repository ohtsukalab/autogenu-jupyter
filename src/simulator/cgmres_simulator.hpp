#ifndef CGMRES_SIMULATOR_H
#define CGMRES_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "save_simulation_data.hpp"
#include "continuation_gmres.hpp"
#include "multiple_shooting_cgmres.hpp"
#include "ms_cgmres_with_input_saturation.hpp"

namespace nmpcsim {
// Simulates NMPC using the C/GMRES-based methods. Opens file streams and saves 
// simulation data into them.
// NMPCSolver: the solver class. Select from ContinuationGMRES, 
//             MultipleShootingCGMRES, and MSCGMRESWithSaturation.
template <class NMPCSolver>
void simulation(NMPCSolver& nmpc, 
                const double* initial_state_vec, const double start_time, 
                const double end_time, const double sampling_period, 
                const std::string save_dir, const std::string savefile_name);
}

#endif // CGMRES_SIMULATOR_H