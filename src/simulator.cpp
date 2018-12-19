#include "simulator.hpp"



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


Simulator::Simulator(const NMPCModel model) : NumericalIntegrator(model)
{
    model_ = model;
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
    current_state_vec=initial_state_vec;
    std::cout << "Start simulation" << std::endl;
    for(double current_time=start_time; current_time<end_time; current_time+= sampling_period){
        saveData(state_data, control_input_data, error_data, current_time, current_state_vec, control_input_vec, cgmres_solver.getError(current_time, current_state_vec));

        // compute the next state vector by the 4th Runge-Kutta-Gill method
        rungeKuttaGill(current_time, current_state_vec, control_input_vec, sampling_period, next_state_vec);

        start_clock = std::chrono::system_clock::now();
        cgmres_solver.controlUpdate(current_time, sampling_period, current_state_vec, control_input_vec);
        end_clock = std::chrono::system_clock::now();

        // compute the computational time of the control update by the C/GMRES method
        double step_time = std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count();
        step_time *= 1e-06;
        total_time += step_time;

        current_state_vec = next_state_vec;
    }
    std::cout << "End simulation" << std::endl;
    std::cout << "CPU time: " << total_time << " [sec]" << std::endl;
    
    // save simulation conditions
    conditions_data << savefile_name << "\n";
    conditions_data << "simulation time: " << end_time-start_time << " [sec]\n";
    conditions_data << "CPU time (total): " << total_time << " [sec]\n";
    conditions_data << "sampling time: " << sampling_period << " [sec]\n";
    conditions_data << "CPU time (1step): " << total_time/((int)( (end_time-start_time)/(sampling_period))) << " [sec]\n";

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
    current_state_vec=initial_state_vec;
    std::cout << "Start simulation" << std::endl;
    for(double current_time=start_time; current_time<end_time; current_time+= sampling_period){
        saveData(state_data, control_input_data, error_data, current_time, current_state_vec, control_input_vec, cgmres_solver.getError(current_time, current_state_vec));

        // compute the next state vector by the 4th Runge-Kutta-Gill method
        rungeKuttaGill(current_time, current_state_vec, control_input_vec, sampling_period, next_state_vec);

        start_clock = std::chrono::system_clock::now();
        cgmres_solver.controlUpdate(current_time, sampling_period, current_state_vec, control_input_vec);
        end_clock = std::chrono::system_clock::now();

        // compute the computational time of the control update by the C/GMRES method
        double step_time = std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count();
        step_time *= 1e-06;
        total_time += step_time;

        current_state_vec = next_state_vec;
    }
    std::cout << "End simulation" << std::endl;
    std::cout << "CPU time: " << total_time << " [sec]" << std::endl;
    
    // save simulation conditions
    conditions_data << savefile_name << "\n";
    conditions_data << "simulation time: " << end_time-start_time << " [sec]\n";
    conditions_data << "CPU time (total): " << total_time << " [sec]\n";
    conditions_data << "sampling time: " << sampling_period << " [sec]\n";
    conditions_data << "CPU time (1step): " << total_time/((int)( (end_time-start_time)/(sampling_period))) << " [sec]\n";

    state_data.close();
    control_input_data.close();
    error_data.close();
    conditions_data.close();
}