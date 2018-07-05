#include "nmpc_model.hpp"
#include "newton_gmres.hpp"
#include "continuation_gmres.hpp"
#include "simulator.hpp"


int main()
{
    // set solver and parameters
    continuation_gmres solver(1.0, 1.0, 50, 5, 1.0e-05, 1000);
    simulator sim;

    // initial state
    Eigen::VectorXd x0(4), u(2);
    x0 << 0.0, 0.0, 0.0, 0.0;

    Eigen::VectorXd u0(2);
    u0 << 1.0, 1.0;

    solver.init_cgmres(0, x0, u0, 1.0e-06, 50);
    sim.simulation(solver, x0, 10, 0.001, "example");

    return 0;
}

/*
int main()
{
    // set solver and parameters
    nmpc_solver solver(solverid::newton_gmres_single, 0.5, 10, 1.0e-01, 1.0e-04, 5);
    simulator sim;

    // initial state
    Eigen::VectorXd x0(2);
    x0 << 2.0, 0.0;

    // run simulation
    sim.simulation(solver, x0, 10, 0.01, "example_");

    return 0;
}
*/