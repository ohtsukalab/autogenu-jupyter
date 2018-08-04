#include "nmpc_model.hpp"
#include "continuation_gmres.hpp"
#include "simulator.hpp"


int main()
{
    // set solver and parameters
    continuation_gmres solver(0.5, 1.0, 50, 5, 1.0e-05, 1000);
    simulator sim;

    // initial state
    Eigen::VectorXd x0(4);
    x0 << 0.0, 0.0, 0.0, 0.0;

    Eigen::VectorXd u0(1);
    u0 << 0.01;

    solver.init(0, x0, u0, 1.0e-06, 50);
    sim.simulation(solver, x0, 10, 0.001, "example");

    return 0;
}
