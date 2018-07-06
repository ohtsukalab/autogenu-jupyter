#include "nmpc_model.hpp"
#include "newton_gmres.hpp"
#include "continuation_gmres.hpp"
#include "simulator.hpp"


int main()
{
    // set solver and parameters
    newton_gmres solver(0.5, 10, 1.0e-01, 1.0e-04, 5);
    simulator sim;

    // initial state
    Eigen::VectorXd x0(2);
    x0 << 2.0, 0.0;

    // run simulation
    sim.simulation(solver, x0, 10, 0.01, "example_");

    return 0;
}
