#include "nmpc_model.hpp"
#include "nmpc_solver.hpp"
#include "simulator.hpp"


int main()
{
    nmpc_solver solver(solverid::newton_gmres_single, 1.0, 50, 5, 1.0e-06, 5);
    simulator sim;

    // initial state
    Eigen::VectorXd x0(4);
    x0 << 0.0, 0.0, 0.0, 0.0;

    sim.simulation(solver, x0, 1, 0.001, "example_");

    return 0;
}
