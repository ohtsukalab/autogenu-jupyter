#ifndef CGMRES__SOLVER_SETTINGS_HPP_
#define CGMRES__SOLVER_SETTINGS_HPP_

#include <iostream>

#include "cgmres/types.hpp"

namespace cgmres {

///
/// @class SolverSettings
/// @brief Settings of solvers. 
///
struct SolverSettings {
  ///
  /// @brief Maximum number of iterations of the ZeroHorizonOCPSolver method. 
  /// Has nothing to do with SingleShootingCGMRESSolver or MultipleShootingCGMRESSolver. 
  /// Default value is 100.
  ///
  size_t max_iter = 100;

  ///
  /// @brief Termination criterion of the ZeroHorizonOCPSolver method. 
  /// Has nothing to do with SingleShootingCGMRESSolver or MultipleShootingCGMRESSolver. 
  /// Must be non-negative. Default value is 1.0e-04.
  ///
  Scalar opterr_tol = 1.0e-04;

  ///
  /// @brief Epsilon of the finite difference approximation. Must be positive.
  /// Default value is 1.0e-08.
  ///
  Scalar finite_difference_epsilon = 1.0e-08;

  ///
  /// @brief The sampling time of MPC and used in SingleShootingCGMRESSolver
  /// and MultipleShootingCGMRESSolver. Has nothing to do with ZeroHorizonOCPSolver. 
  /// Must be positive. Default is 0.001.
  ///
  Scalar sampling_time = 0.001; 

  ///
  /// @brief The stabilization parameter of the continuation method. 
  /// Typical value is the reciprocal of SolverSettings::sampling_time (sampling period of MPC).
  /// Used in SingleShootingCGMRESSolver and MultipleShootingCGMRESSolver. 
  /// Has nothing to do with ZeroHorizonOCPSolver. 
  /// Must be positive. Default is 1000.0.
  ///
  Scalar zeta = 1000.0; 

  ///
  /// @brief The minimum value of the dummy inputs. 
  /// Mainly used in MultipleShootingCGMRESSolver. 
  /// In SingleShootingCGMRESSolver and ZeroHorizonOCPSolver, this value is 
  /// used only in SingleShootingCGMRESSolver::init_dummy_mu() and 
  /// ZeroHorizonOCPSolver::init_dummy_mu(). Must be non-negaive. Default is 1.0e-03.
  ///
  Scalar min_dummy = 1.0e-03;

  ///
  /// @brief Verbose level. 0: no printings. 1-2: print some things. Default is 0.
  ///
  size_t verbose_level = 0;

  ///
  /// @brief If true, a solver profile is taken.
  ///
  bool profile_solver = true;

  void disp(std::ostream& os) const {
    os << "Soler settings: " << std::endl; 
    os << "  max iter:                  " << max_iter << std::endl;
    os << "  opterr tol:                " << opterr_tol << std::endl;
    os << "  finite difference epsilon: " << finite_difference_epsilon << std::endl;
    os << "  sampling_time:             " << sampling_time << std::endl;
    os << "  zeta:                      " << zeta << std::endl;
    os << "  min dummy:                 " << min_dummy << std::endl;
    os << "  verbose level:             " << verbose_level << std::endl;
    os << "  profile solver:            " << std::boolalpha << profile_solver << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const SolverSettings& settings) {
    settings.disp(os);
    return os;
  }
};

} // namespace cgmres

#endif // CGMRES__SOLVER_SETTINGS_HPP_