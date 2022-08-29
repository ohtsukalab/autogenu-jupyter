#include "ocp.hpp"

#include "cgmres/zero_horizon_ocp_solver.hpp"
#include "cgmres/multiple_shooting_cgmres_solver.hpp"

#include <string>

int main() {
  // Define the optimal control problem.
  cgmres::OCP_mobilerobot ocp;

  // Define the horizon.
  const double Tf = 1.5;
  cgmres::Horizon horizon(Tf); // fixed-length

  // Define the solver settings.
  cgmres::SolverSettings settings;
  settings.sampling_time = 0.001;  
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
  const double tsim = 10.0;
  const double sampling_time = settings.sampling_time;
  const int sim_steps = std::floor(tsim / sampling_time);

  double t = t0;
  cgmres::VectorX x = x0;
  cgmres::VectorX dx = cgmres::VectorX::Zero(x0.size());
  for (int i=0; i<sim_steps; ++i) {
    const auto& u = mpc.uopt()[0]; // const reference to the initial optimal control input 
    dx.setZero();
    ocp.eval_f(t, x, u, dx); // eval the state equation
    const cgmres::VectorX x1 = x + sampling_time * dx;
    mpc.update(t, x);
    x = x1;
    t = t + sampling_time;
    std::cout << "t: " << t << ", x: " << x.transpose() << std::endl;
  }

  std::cout << "\n======================= MPC used in this simulation: =======================" << std::endl;
  std::cout << mpc << std::endl;

  return 0;
}
