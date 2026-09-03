#include "cgmres/horizon.hpp"
#include "cgmres/integrator.hpp"
#include "cgmres/solver_settings.hpp"

#include <cmath>
#include <iostream>

namespace {

class StableLinearOCP {
public:
  template <typename State, typename Control, typename Derivative>
  void eval_f(const cgmres::Scalar, const cgmres::MatrixBase<State>& x,
              const cgmres::MatrixBase<Control>& u,
              const cgmres::MatrixBase<Derivative>& dx) const {
    const_cast<Derivative&>(dx.derived()).coeffRef(0) = -x.coeff(0) + u.coeff(0);
  }
};

} // namespace

int main() {
  const cgmres::Horizon horizon(2.0);
  if (horizon.T(1.0) != 2.0) {
    std::cerr << "fixed horizon returned an unexpected length\n";
    return 1;
  }

  const StableLinearOCP ocp;
  cgmres::Vector<1> x;
  cgmres::Vector<1> u;
  x << 1.0;
  u << 0.0;
  const cgmres::VectorX next = cgmres::RK4(ocp, 0.0, 0.1, x, u);
  if (std::abs(next.coeff(0) - std::exp(-0.1)) > 1.0e-6) {
    std::cerr << "RK4 integration produced an unexpected state\n";
    return 1;
  }

  const cgmres::SolverSettings settings;
  if (settings.sampling_time <= 0.0 || settings.finite_difference_epsilon <= 0.0) {
    std::cerr << "default solver settings are invalid\n";
    return 1;
  }

  return 0;
}
