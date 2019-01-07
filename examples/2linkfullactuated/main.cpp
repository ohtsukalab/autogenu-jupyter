#include "nmpc_model.hpp"
#include "continuation_gmres.hpp"
#include "multiple_shooting_cgmres.hpp"
#include "control_input_saturation_sequence.hpp"
#include "multiple_shooting_with_saturation.hpp"
#include "simulator.hpp"


int main()
{
    // Define the model in NMPC.
    NMPCModel nmpc_model;


    // Define the solver of C/GMRES.

    // ContinuationGMRES cgmres_solver(nmpc_model, 0.5, 1.0, 50, 1.0e-06, 1000, 5);

    // MultipleShootingCGMRES cgmres_solver(nmpc_model, 0.5, 1.0, 50, 1.0e-06, 1000, 5);

    //  If you use MultipleShootingWithSaturation, you have to describe the saturaions on the control input in ControlInputSaturationSequence.
    ControlInputSaturationSequence control_input_saturation_seq;
    control_input_saturation_seq.appendControlInputSaturation(0, -3, 3, 1.0e-05);
    control_input_saturation_seq.appendControlInputSaturation(1, -5, 5, 1.0e-05);
    MultipleShootingWithSaturation cgmres_solver(nmpc_model, control_input_saturation_seq, 0.5, 1.0, 50, 1.0e-06, 1000, 5);


    // Define the simulator.
    Simulator cgmres_simulator(nmpc_model);

    // Set the initial state.
    Eigen::VectorXd initial_state(nmpc_model.dimState());
    initial_state = Eigen::VectorXd::Zero(nmpc_model.dimState());

    // Set the initial guess of the control input vector.
    Eigen::VectorXd initial_guess_control_input(nmpc_model.dimControlInput()+nmpc_model.dimConstraints());
    initial_guess_control_input = Eigen::VectorXd::Zero(nmpc_model.dimControlInput()+nmpc_model.dimConstraints());

    // Initialize the solution of the C/GMRES method.
    Eigen::VectorXd initial_guess_lagrange_multiplier(control_input_saturation_seq.dimSaturation());
    initial_guess_lagrange_multiplier << 1.0e-03, 1.0e-03;
    cgmres_solver.initSolution(0, initial_state, initial_guess_control_input, initial_guess_lagrange_multiplier, 1.0e-06, 50);

    // Perform a numerical simulation.
    cgmres_simulator.simulation(cgmres_solver, initial_state, 0, 5, 0.001, "example");


    return 0;
}
