#ifndef CGMRES__SIMULATOR_HPP_
#define CGMRES__SIMULATOR_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "cgmres/types.hpp"
#include "cgmres/simulator/integrator.hpp"


namespace cgmres {

template <typename OCP, typename NMPC, typename VectorType>
void simulation(const OCP& ocp, NMPC& nmpc, const VectorType& x0, 
                const Scalar t0, const Scalar tf, const double dt, 
                const std::string save_dir, 
                const std::string savefile_name) {
  VectorX x = x0;
  std::chrono::system_clock::time_point start_clock, end_clock;

  std::string savefile_header = save_dir + "/" + savefile_name;
  std::ofstream x_log(savefile_header + "_x.log"), 
                u_log(savefile_header + "_u.log"), 
                opterr_log(savefile_header + "_opterr.log"),
                conditions_log(savefile_header + "_conditions.log");

  Scalar total_CPU_time = 0.0;
  std::cout << "Start simulation" << std::endl;
  for (Scalar t=t0; t<tf; t+=dt) {
    // take logs
    x_log << x.transpose() << '\n';
    u_log << nmpc.uopt()[0].transpose() << '\n';
    opterr_log << nmpc.optError(t, x) << '\n';

    // Computes the next state vector using the 4th Runge-Kutta method.
    const VectorX x1 = RK4(ocp, t, dt, x, nmpc.uopt()[0]);

    // Updates the solution and measure the computational time of the update.
    start_clock = std::chrono::system_clock::now();
    nmpc.update(t, x);
    end_clock = std::chrono::system_clock::now();

    // Converts the computational time to seconds.
    Scalar CPU_time = 
        std::chrono::duration_cast<std::chrono::microseconds>(
            end_clock-start_clock).count();
    CPU_time *= 1e-03;
    total_CPU_time += CPU_time;

    // Updates the state.
    x = x1;
  }

  // cout the simulation conditions.
  std::cout << "End simulation\n" 
      << "Total CPU time for control update: " << total_CPU_time << " [ms]\n" 
      << "sampling time: " << 1000.0 * dt << " [ms]" << "\n" 
      << "CPU time for per control update: " 
      << total_CPU_time/((int)( (tf-t0)/dt)) << " [ms]" << std::endl;

  // Save simulation conditions.
  conditions_log << "simulation: " << savefile_name << "\n"
      << "simulation time: " << tf-t0 << " [sec]\n"
      << "Total CPU time for control update: " << 0.001 * total_CPU_time << " [s]\n"
      << "sampling time: " << 1000.0 * dt << " [ms]\n"
      << "CPU time for per control update: " 
      << total_CPU_time/((int)( (tf-t0)/dt)) << " [ms]" << std::endl;

  x_log.close();
  u_log.close();
  opterr_log.close();
  conditions_log.close();
}

} // namespace cgmres

#endif // CGMRES__SIMULATOR_HPP_