#include "nmpc_model.hpp"
#include "control_input_saturation_sequence.hpp"
#include "multiple_shooting_cgmres_with_saturation.hpp"
#include "multiple_shooting_cgmres_with_saturation_simulator.hpp"


int main()
{
    NMPCModel nmpc_model;

    ControlInputSaturationSequence control_input_saturation_seq;
    control_input_saturation_seq.appendControlInputSaturation(0, -3, 3, 1.0e-02, 0);
    control_input_saturation_seq.appendControlInputSaturation(1, -1.5, 1.5, 1.0e-02, 0);
    MultipleShootingCGMRESWithSaturation nmpc_solver(control_input_saturation_seq, 0.5, 1.0, 50, 1.0e-06, 1000, 5);


    // Set the initial state.
    double initial_state[4] = {0};

    // Set the initial guess of the control input vector.

    double initial_guess_control_input[] = {0.1, 0.1};

    // Initialize the solution of the C/GMRES method.
    double initial_guess_lagrange_multiplier[2] = {1.0e-03, 1.0e-03};
    nmpc_solver.initSolution(0, initial_state, initial_guess_control_input, initial_guess_lagrange_multiplier, 1.0e-06, 50);

    // Perform a numerical simulation.
    nmpcsim::simulation(nmpc_solver, initial_state, 0, 5, 0.001, "example");

    return 0;
}
