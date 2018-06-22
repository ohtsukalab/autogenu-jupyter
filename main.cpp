#include "nmpc_model.hpp"
#include "nmpc_solver.hpp"
#include "simulator.hpp"


int main()
{
    // set solver and parameters
    nmpc_solver solver(solverid::newton_gmres_single, 0.5, 50, 1.0e-01, 1.0e-04, 5);
    simulator sim;

    // initial state
    Eigen::VectorXd x0(2);
    x0 << 2.0, 0.0;

    // run simulation
    sim.simulation(solver, x0, 10, 0.01, "example_");

    return 0;
}
