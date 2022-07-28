#ifndef CGMRES__SOLVER_SETTINGS_HPP_
#define CGMRES__SOLVER_SETTINGS_HPP_

#include <iostream>

#include "cgmres/types.hpp"

namespace cgmres {

struct SolverSettings {
  size_t max_iter = 100;
  Scalar opterr_tol = 1.0e-04;

  Scalar finite_difference_epsilon = 1.0e-08;

  Scalar dt = 0.001; // sampling period
  Scalar zeta = 1000.0; // hint: 1.0/dt 

  size_t verbose_level = 0;

  void disp(std::ostream& os) const {
    os << "Soler settings: " << std::endl; 
    os << "  max iter:                  " << max_iter << std::endl;
    os << "  opterr tol:                " << opterr_tol << std::endl;
    os << "  finite difference epsilon: " << finite_difference_epsilon << std::endl;
    os << "  dt (sampling period):      " << dt << std::endl;
    os << "  zeta:                      " << zeta << std::endl;
    os << "  verbose level:             " << verbose_level << std::flush;
  }

  friend std::ostream& operator<<(std::ostream& os, const SolverSettings& settings) {
    settings.disp(os);
    return os;
  }
};

} // namespace cgmres

#endif // CGMRES__SOLVER_SETTINGS_HPP_