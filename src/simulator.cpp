#include "simulator.hpp"
#include <ros/ros.h>

void Simulator::saveData(std::ofstream& state_data, std::ofstream& control_input_data, std::ofstream& error_data, const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& control_input_vec, const double optimality_error)
{
    for(int i=0; i<model_.dimState(); i++){
        state_data << state_vec(i) << " ";
    }
    state_data << "\n";

    for(int i=0; i<model_.dimControlInput(); i++){
        control_input_data << control_input_vec(i) << " ";
    }
    control_input_data << "\n";

    error_data << optimality_error << "\n";
}


void Simulator::initModel(NMPCModel model)
{
    model_ = model;

    k1_vec_.resize(model_.dimState());
    k2_vec_.resize(model_.dimState());
    k3_vec_.resize(model_.dimState());
    k4_vec_.resize(model_.dimState());

}

void Simulator::euler(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec)
{
    model_.stateFunc(current_time, current_state_vec, control_input_vec, dx_vec_);
    next_state_vec = current_state_vec + integration_length * dx_vec_;
}


void Simulator::rungeKuttaGill(const double current_time, const Eigen::VectorXd& current_state_vec, const Eigen::VectorXd& control_input_vec, const double integration_length, Eigen::Ref<Eigen::VectorXd> next_state_vec)
{
    model_.stateFunc(current_time, current_state_vec, control_input_vec, k1_vec_);
    model_.stateFunc(current_time+0.5*integration_length, current_state_vec+0.5*integration_length*k1_vec_, control_input_vec, k2_vec_);
    model_.stateFunc(current_time+0.5*integration_length, current_state_vec+integration_length*0.5*(std::sqrt(2)-1)*k1_vec_+integration_length*(1-(1/std::sqrt(2)))*k2_vec_, control_input_vec, k3_vec_);
    model_.stateFunc(current_time+integration_length, current_state_vec-integration_length*0.5*std::sqrt(2)*k2_vec_+integration_length*(1+(1/std::sqrt(2)))*k3_vec_, control_input_vec, k4_vec_);
    next_state_vec = current_state_vec + integration_length*(k1_vec_+(2-std::sqrt(2))*k2_vec_+(2+std::sqrt(2))*k3_vec_+k4_vec_)/6;
}



void Simulator::simulation(ContinuationGMRES cgmres_solver, const Eigen::VectorXd& initial_state_vec, const double start_time, const double end_time, const double sampling_period, const std::string savefile_name)
{
    Eigen::VectorXd current_state_vec(model_.dimState()), next_state_vec(model_.dimState()), control_input_vec(model_.dimControlInput());
    std::chrono::system_clock::time_point start_clock, end_clock;

    std::ofstream state_data(savefile_name + "_state.dat");
    std::ofstream control_input_data(savefile_name + "_control_input.dat");
    std::ofstream error_data(savefile_name + "_error.dat");
    std::ofstream conditions_data(savefile_name + "_conditions.dat");

    double total_time = 0;
    current_state_vec = initial_state_vec;
    std::cout << "Start simulation" << std::endl;
    for(double current_time=start_time; current_time<end_time; current_time+= sampling_period){
        saveData(state_data, control_input_data, error_data, current_time, current_state_vec, control_input_vec, cgmres_solver.getError(current_time, current_state_vec));

        // Compute the next state vector using the 4th Runge-Kutta-Gill method.
        rungeKuttaGill(current_time, current_state_vec, control_input_vec, sampling_period, next_state_vec);

        // Update the solution and measure computational time.
        start_clock = std::chrono::system_clock::now();
        cgmres_solver.controlUpdate(current_time, sampling_period, current_state_vec, control_input_vec);
        end_clock = std::chrono::system_clock::now();

        // Convert computational time to seconds.
        double step_time = std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count();
        step_time *= 1e-06;
        total_time += step_time;

        current_state_vec = next_state_vec;
    }
    std::cout << "End simulation" << std::endl;
    std::cout << "Total CPU time for control update: " << total_time << " [sec]" << std::endl;
    
    // Save simulation conditions.
    conditions_data << "simulation name: " << savefile_name << "\n";
    conditions_data << "simulation time: " << end_time-start_time << " [sec]\n";
    conditions_data << "CPU time for control update (total): " << total_time << " [sec]\n";
    conditions_data << "sampling time: " << sampling_period << " [sec]\n";
    conditions_data << "CPU time for control update (1step): " << total_time/((int)( (end_time-start_time)/(sampling_period))) << " [sec]\n";

    state_data.close();
    control_input_data.close();
    error_data.close();
    conditions_data.close();
}



void Simulator::simulation(MultipleShootingCGMRES cgmres_solver, const Eigen::VectorXd& initial_state_vec, const double start_time, const double end_time, const double sampling_period, const std::string savefile_name)
{

    Eigen::VectorXd current_state_vec(model_.dimState()), next_state_vec(model_.dimState()), control_input_vec(model_.dimControlInput());
    std::chrono::system_clock::time_point start_clock, end_clock;

    std::ofstream state_data(savefile_name + "_state.dat");
    std::ofstream control_input_data(savefile_name + "_control_input.dat");
    std::ofstream error_data(savefile_name + "_error.dat");
    std::ofstream conditions_data(savefile_name + "_conditions.dat");

    double total_time = 0;
    current_state_vec = initial_state_vec;
    std::cout << "Start simulation" << std::endl;
    for(double current_time=start_time; current_time<end_time; current_time+= sampling_period){
        // saveData(state_data, control_input_data, error_data, current_time, current_state_vec, control_input_vec, cgmres_solver.getError(current_time, current_state_vec));

        // Compute the next state vector using the 4th Runge-Kutta-Gill method.
        rungeKuttaGill(current_time, current_state_vec, control_input_vec, sampling_period, next_state_vec);

        // Update the solution and measure computational time.
        start_clock = std::chrono::system_clock::now();
        cgmres_solver.controlUpdate(current_time, sampling_period, current_state_vec, control_input_vec);
        end_clock = std::chrono::system_clock::now();

        // Convert computational time to seconds.
        double step_time = std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count();
        step_time *= 1e-06;
        total_time += step_time;

        current_state_vec = next_state_vec;
    }
    std::cout << "End simulation" << std::endl;
    std::cout << "Total CPU time for control update: " << total_time << " [sec]" << std::endl;

    // Save simulation conditions.
    conditions_data << "simulation name: " << savefile_name << "\n";
    conditions_data << "simulation time: " << end_time-start_time << " [sec]\n";
    conditions_data << "CPU time for control update (total): " << total_time << " [sec]\n";
    conditions_data << "sampling time: " << sampling_period << " [sec]\n";
    conditions_data << "CPU time for control update (1step): " << total_time/((int)( (end_time-start_time)/(sampling_period))) << " [sec]\n";

    state_data.close();
    control_input_data.close();
    error_data.close();
    conditions_data.close();
}