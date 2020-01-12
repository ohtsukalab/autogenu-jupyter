#ifndef CGMRES_SIMULATOR_H
#define CGMRES_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "save_simulation_data.hpp"


namespace cgmres {

// Simulates NMPC using the C/GMRES-based methods. Opens file streams and saves 
// simulation data into them.
// NMPCSolver: the solver class. Select from ContinuationGMRES, 
//             MultipleShootingCGMRES, and MSCGMRESWithInputSaturation.
template <class NMPCSolver>
void simulation(NMPCSolver& nmpc, const double* initial_state_vec, 
                const double start_time, const double end_time, 
                const double sampling_period, const std::string save_dir, 
                const std::string savefile_name) {
  NMPCModel model;
  NumericalIntegrator integrator;
  double current_state_vec[model.dim_state()], next_state_vec[model.dim_state()],
      control_input_vec[model.dim_control_input()];
  std::chrono::system_clock::time_point start_clock, end_clock;

  std::string savefile_header = save_dir + "/" + savefile_name;
  std::ofstream state_data(savefile_header + "_state.dat"), 
    control_input_data(savefile_header + "_control_input.dat"), 
    error_data(savefile_header + "_error.dat"),
    conditions_data(savefile_header + "_conditions.dat");

  double total_time = 0;
  for (int i=0; i<model.dim_state(); i++) {
    current_state_vec[i] = initial_state_vec[i];
  }
  nmpc.initializeSolution(start_time, current_state_vec);
  nmpc.getControlInput(control_input_vec);

  std::cout << "Start simulation" << std::endl;
  for (double current_time=start_time; current_time<end_time; 
       current_time+=sampling_period) {
    // Saves the current datas.
    saveData(model.dim_state(), model.dim_control_input(), state_data, 
             control_input_data, error_data, current_time, current_state_vec, 
             control_input_vec, 
             nmpc.getErrorNorm(current_time, current_state_vec));

    // Computes the next state vector using the 4th Runge-Kutta-Gill method.
    integrator.rungeKuttaGill(current_time, current_state_vec, 
                              control_input_vec, sampling_period, 
                              next_state_vec);

    // Updates the solution and measure the computational time of the update.
    start_clock = std::chrono::system_clock::now();
    nmpc.controlUpdate(current_time, current_state_vec, 
                       sampling_period, control_input_vec);
    end_clock = std::chrono::system_clock::now();

    // Converts the computational time to seconds.
    double step_time = 
        std::chrono::duration_cast<std::chrono::microseconds>(
            end_clock-start_clock).count();
    step_time *= 1e-06;
    total_time += step_time;

    // Updates the state.
    for (int i=0; i<model.dim_state(); i++) {
      current_state_vec[i] = next_state_vec[i];
    }
  }

  // cout the simulation conditions.
  std::cout << "End simulation\n" 
      << "Total CPU time for control update: " << total_time << " [sec]\n" 
      << "sampling time: " << sampling_period << " [sec]" << "\n" 
      << "CPU time for per control update: " 
      << total_time/((int)( (end_time-start_time)/(sampling_period))) 
      << " [sec]" << std::endl;

  // Save simulation conditions.
  conditions_data << "simulation name: " << savefile_name << "\n"
      << "simulation time: " << end_time-start_time << " [sec]\n"
      << "Total CPU time for control update: " << total_time << " [sec]\n"
      << "sampling time: " << sampling_period << " [sec]\n"
      << "CPU time for per control update: " 
      << total_time/((int)((end_time-start_time)/(sampling_period))) 
      << " [sec]\n";

  state_data.close();
  control_input_data.close();
  error_data.close();
  conditions_data.close();
}

} // namespace cgmres

#endif // CGMRES_SIMULATOR_H