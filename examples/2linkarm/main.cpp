#include "nmpc_model.hpp"
#include "continuation_gmres.hpp"
#include "multiple_shooting_cgmres.hpp"
#include "simulator.hpp"


int main()
{
    // set solver and parameters
    NMPCModel nmpc_model;
    ContinuationGMRES cgmres_solver(nmpc_model, 0.5, 1.0, 50, 1.0e-06, 1000, 5);
    Simulator cgmres_simulator(nmpc_model);

    Eigen::VectorXd initial_state(4);
    initial_state << 0.0, 0.0, 0.0, 0.0;

    // initial guess of the control input
    Eigen::VectorXd initial_guess_control_input(2);
    initial_guess_control_input << 0.0, 0.0;

    cgmres_solver.initSolution(0, initial_state, initial_guess_control_input, 1.0e-06, 50);
    cgmres_simulator.simulation(cgmres_solver, initial_state, 0, 10, 0.001, "example");

    return 0;
}
