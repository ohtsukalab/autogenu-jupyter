#ifndef SOLVER_SETTINGS_HPP_
#define SOLVER_SETTINGS_HPP_

#include "types.hpp"

namespace cgmres {

struct SolverSettings {
  size_t max_iter = 100;
  Scalar opt_error_tol = 1.0e-04;
  Scalar finite_diference_epsilon = 1.0e-08;
  size_t verbose_level = 0;
};

} // namespace cgmres

#endif // SOLVER_SETTINGS_HPP_