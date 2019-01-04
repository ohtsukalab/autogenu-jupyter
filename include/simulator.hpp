//
//
// Simulator of the C/GMRES method and multiple shooting based C/GMRES method.
//

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <pick_model.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <eigen3/Eigen/Core>
#include "continuation_gmres.hpp"
#include "multiple_shooting_cgmres.hpp"


// Simulates the C/GMRES method and multiple shooting based C/GMRES method.
// Opens file streams and saves simulation data.

class Simulator 
{
private:
    NMPCModel model_;
    void saveData(std::ofstream& state_data, std::ofstream& control_input_data, std::ofstream& error_data, const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_vec, const double optimality_error);

    /* Copied from numerical integrator */ 
	    Eigen::VectorXd dx_vec_, k1_vec_, k2_vec_, k3_vec_, k4_vec_;
public:
    // Simulator(const NMPCModel model);
    void initModel(NMPCModel model);

    // Serves a numerical simulation for the C/GMRES and the multiple shooting based C/GMRES.
    void simulation(ContinuationGMRES cgmres_solver, const Eigen::VectorXd& initial_state_vec, const double start_time, const double end_time, const double sampling_period, const std::string savefile_name);
    void simulation(MultipleShootingCGMRES cgmres_solver, const Eigen::VectorXd& initial_state_vec, const double start_time, const double end_time, const double sampling_period, const std::string savefile_name);

	/* Copied from numerical integrator */ 
		// Euler method for the state equation.
		void euler(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec);

		// The four-step Runge-Kutta-Gill method for the state equation.
		void rungeKuttaGill(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec);

};


#endif