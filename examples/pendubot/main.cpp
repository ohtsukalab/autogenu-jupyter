#include "nmpc_model.hpp"
// #include "continuation_gmres.hpp"
// #include "cgmres_simulator.hpp"
#include "multiple_shooting_cgmres.hpp"
#include "multiple_shooting_cgmres_simulator.hpp"


int main()
{
    // Define the model in NMPC.
    NMPCModel nmpc_model;

    // Define the solver of C/GMRES.
    // ContinuationGMRES nmpc_solver(1.0, 1.0, 50, 1.0e-06, 1000, 3);
    MultipleShootingCGMRES nmpc_solver(1.0, 1.0, 50, 1.0e-06, 1000, 3);

    // Set the initial state.
    double initial_state[4] = {0};

    // Set the initial guess of the control input vector.
    double initial_guess_control_input[1] = {0};


    // Initialize the solution of the C/GMRES method.
    nmpc_solver.initSolution(0, initial_state, initial_guess_control_input, 1.0e-06, 50);


    // Perform a numerical simulation.
    nmpcsim::simulation(nmpc_solver, initial_state, 0, 10, 0.001, "example1");


    return 0;
}