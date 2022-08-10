#include "ocp.hpp"

#include "cgmres/zero_horizon_ocp_solver.hpp"
#include "cgmres/multiple_shooting_cgmres_solver.hpp"

#include <string>
#include <memory>

int main() {
  // Define the optimal control problem.
  cgmres::OCP_cartpoleExternalReference ocp;

  // set the external reference ptr
  auto external_refernce = std::make_shared<cgmres::OCP_cartpoleExternalReference::ExternalReference>();
  external_refernce->cart_position = 0.0;
  ocp.external_refernce = external_refernce;

  // Define the horizon.
  const double Tf = 2.0;
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
  cgmres::Vector<4> x0;
  x0 << 0, 0, 0, 0;

  // Initialize the solution of the C/GMRES method.
  constexpr int kmax_init = 1;
  cgmres::ZeroHorizonOCPSolver<cgmres::OCP_cartpoleExternalReference, kmax_init> initializer(ocp, settings);
  cgmres::Vector<1> uc0;
  uc0 << 0.01;
  initializer.set_uc(uc0);
  initializer.solve(t0, x0);

  // Define the C/GMRES solver.
  constexpr int N = 100;
  constexpr int kmax = 5;
  cgmres::MultipleShootingCGMRESSolver<cgmres::OCP_cartpoleExternalReference, N, kmax> mpc(ocp, horizon, settings);
  mpc.set_uc(initializer.ucopt());
  mpc.init_x_lmd(t0, x0);
  mpc.init_dummy_mu();

  // Perform a numerical simulation.
  const double tsim = 10.0;
  const double dt = settings.dt;
  const int sim_steps = std::floor(tsim / dt);

  double t = t0;
  cgmres::VectorX x = x0;
  cgmres::VectorX dx = cgmres::VectorX::Zero(x0.size());
  for (int i=0; i<sim_steps; ++i) {
    const auto& u = mpc.uopt()[0]; // const reference to the initial optimal control input 
    dx.setZero();
    ocp.eval_f(t, x.data(), u.data(), dx.data()); // eval the state equation
    const cgmres::VectorX x1 = x + dt * dx;
    mpc.update(t, x);
    x = x1;
    t = t + dt;
    std::cout << "t: " << t << ", x: " << x.transpose() << std::endl;

    // switch the reference
    if (t > 5.0) {
      external_refernce->cart_position = 1.0;
    }
  }

  std::cout << "\n======================= MPC used in this simulation: =======================" << std::endl;
  std::cout << mpc << std::endl;

  return 0;
}
