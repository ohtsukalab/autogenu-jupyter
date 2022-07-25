#include "cgmres/examples/cartpole.hpp"
#include "cgmres/zero_horizon_ocp_solver.hpp"
#include "cgmres/single_shooting_ocp_solver.hpp"
#include "cgmres/single_shooting_cgmres_solver.hpp"
#include "cgmres/multiple_shooting_cgmres_solver.hpp"

#include <iostream>

int main() {
  cgmres::CartPoleOCP ocp;
  cgmres::Horizon horizon(0.1, 1.0);
  constexpr int N = 10;

  cgmres::SolverSettings settings;
  settings.verbose_level = 2;

  // initialize solution
  constexpr int kmax_init = cgmres::CartPoleOCP::nuc;
  cgmres::ZeroHorizonOCPSolver<cgmres::CartPoleOCP, kmax_init> initializer(ocp, settings);

  cgmres::Vector<cgmres::CartPoleOCP::nuc> uc0;
  uc0 << 0.01, 10.0, 0.01;
  initializer.setSolution(uc0);

  const double t = 1.0;
  cgmres::Vector<cgmres::CartPoleOCP::nx> x0;
  x0.setZero();
  initializer.solve(t, x0);

  // run MPC
  constexpr int kmax = 10;
  cgmres::SingleShootingOCPSolver<cgmres::CartPoleOCP, N, kmax> solver(ocp, horizon, settings);
  solver.setSolution(initializer.getSolution());
  solver.solve(t, x0);

  cgmres::SingleShootingCGMRESSolver<cgmres::CartPoleOCP, N, kmax> ss_cgmres(ocp, horizon, settings);
  ss_cgmres.setSolution(initializer.getSolution());
  ss_cgmres.update(t, x0);

  cgmres::MultipleShootingCGMRESSolver<cgmres::CartPoleOCP, N, kmax> ms_cgmres(ocp, horizon, settings);
  ms_cgmres.setSolution(initializer.getSolution());
  ms_cgmres.setState(x0);
  ms_cgmres.update(t, x0);

  return 0;
}