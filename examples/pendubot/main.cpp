#include "nmpc_model.hpp"
#include "continuation_gmres.hpp"
#include "multiple_shooting_cgmres.hpp"
#include "simulator.hpp"


int main()
{
    // set solver and parameters
    NMPCModel nmpc_model;

    // ContinuationGMRES cgmres_solver(nmpc_model, 1.0, 1.0, 50, 1.0e-06, 1000, 5);
    MultipleShootingCGMRES cgmres_solver(nmpc_model, 1.0, 1.0, 50, 1.0e-06, 1000, 5);

    Simulator cgmres_simulator(nmpc_model);

    // initial state
    Eigen::VectorXd initial_state(nmpc_model.dimState());
    initial_state = Eigen::VectorXd::Zero(nmpc_model.dimState());

    // initial guess of the control input
    Eigen::VectorXd initial_guess_control_input(nmpc_model.dimControlInput()+nmpc_model.dimConstraints());
    initial_guess_control_input = Eigen::VectorXd::Zero(nmpc_model.dimControlInput()+nmpc_model.dimConstraints());

    cgmres_solver.initSolution(0, initial_state, initial_guess_control_input, 1.0e-06, 50);
    cgmres_simulator.simulation(cgmres_solver, initial_state, 0, 10, 0.001, "example");

    return 0;
}