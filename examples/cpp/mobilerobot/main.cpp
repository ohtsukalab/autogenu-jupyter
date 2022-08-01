#include "ocp.hpp"
#include "cgmres/zero_horizon_ocp_solver.hpp"
#include "cgmres/multiple_shooting_cgmres_solver.hpp"
#include "cgmres/simulator/simulator.hpp"
#include <string>

int main() {
  // Define the optimal control problem.
  cgmres::OCP_mobilerobot ocp;

  // Define the horizon.
  const double Tf = 1.5;
  cgmres::Horizon horizon(Tf); // fixed-length

  // Define the solver settings.
  cgmres::SolverSettings settings;
  settings.dt = 0.001; // sampling period 
  settings.zeta = 1000;
  settings.finite_difference_epsilon = 1e-08;
  // For initialization.
  settings.max_iter = 50;
  settings.opterr_tol = 1e-06;
  settings.verbose_level = 1;

  // Define the initial time and initial state.
  const double t0 = 0;
  cgmres::Vector<3> x0;
  x0 << 0, 0, 0;

  // Initialize the solution of the C/GMRES method.
  constexpr int kmax_init = 4;
  cgmres::ZeroHorizonOCPSolver<cgmres::OCP_mobilerobot, kmax_init> initializer(ocp, settings);
  cgmres::Vector<4> uc0;
  uc0 << 0.1, 0.1, 0.01, 0.01;
  initializer.set_uc(uc0);
  initializer.solve(t0, x0);

  // Define the C/GMRES solver.
  constexpr int N = 50;
  constexpr int kmax = 10;
  cgmres::MultipleShootingCGMRESSolver<cgmres::OCP_mobilerobot, N, kmax> mpc(ocp, horizon, settings);
  mpc.set_uc(initializer.ucopt());
  mpc.init_x_lmd(t0, x0);
  mpc.init_dummy_mu();


  // Perform a numerical simulation.
  const double tf = 10;
  const double dt = settings.dt;
  const std::string save_dir_name("../simulation_result");
  cgmres::simulation(ocp, mpc, x0, t0, tf, dt, save_dir_name, "mobilerobot");

  std::cout << "\n======================= MPC used in this simulation: =======================" << std::endl;
  std::cout << mpc << std::endl;

  return 0;
}
