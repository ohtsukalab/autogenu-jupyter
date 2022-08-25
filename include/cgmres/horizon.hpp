#ifndef CGMRES__HORIZON_HPP_
#define CGMRES__HORIZON_HPP_

#include <cmath>
#include <stdexcept>
#include <cassert>
#include <iostream>

#include "cgmres/types.hpp"


namespace cgmres {

///
/// @class Horizon
/// @brief Horizon of MPC. 
///
class Horizon {
public:
  ///
  /// @brief Constructs the horizon. If alpha <= 0.0, then the fixed-length Tf is used. If alpha > 0.0, 
  /// then the time-varying length Tf * (1.0-exp(-alpha * (t-t0))) is used.
  /// @param[in] Tf The parameter of the horizon length. Must be positive.
  /// @param[in] alpha The parameter of the time-varying horizon length. Default is 0.0 (i.e., fixed-length horizon).
  /// @param[in] t0 The parameter of the time-varying horizon length. Default is 0.0.
  ///
  Horizon(const Scalar Tf, const Scalar alpha=0.0, const Scalar t0=0.0)
    : Tf_(Tf), alpha_(alpha), t0_(t0) {
    if (Tf <= 0.0) {
      throw std::invalid_argument("[Horizon]: 'Tf' must be positive!");
    }
    time_varying_length_ = (alpha > 0.0);
  }

  ///
  /// @brief Default constructor.
  ///
  Horizon() = default;

  ///
  /// @brief Default destructor.
  ///
  ~Horizon() = default;

  ///
  /// @brief Gets the length of the horizon.
  /// @param[in] t The initial time of the horizon. If this horizon is time-varying (i.e., alpha > 0.0),
  /// then this value must not be less than t0.
  ///
  Scalar T(const Scalar t) const {
    if (time_varying_length_) {
      if (t < t0_) {
        throw std::invalid_argument("[Horizon]: 't' must be greater than or equal to 't0' (" +  std::to_string(t0_) + ") !");
      }
      return Tf_ * (1.0-std::exp(-alpha_*(t-t0_)));
    }
    else {
      return Tf_;
    }
  }

  ///
  /// @brief Resets the length of the horizon (for time-varying horizon).
  /// @param[in] t0 The parameter of the time-varying horizon length. 
  ///
  void reset(const Scalar t0) {
    t0_ = t0;
  }

  void disp(std::ostream& os) const {
    os << "Horizon: "; 
    if (time_varying_length_) {
      os << "time-varying length" << std::endl;
    }
    else {
      os << "fixed length" << std::endl;
    }
    os << "  Tf:    " << Tf_ << std::endl;
    os << "  alpha: " << alpha_ << std::endl;
    os << "  t0:    " << t0_ << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const Horizon& horizon) {
    horizon.disp(os);
    return os;
  }

private:
  Scalar Tf_, alpha_, t0_;
  bool time_varying_length_;
};

} // namespace cgmres

#endif // CGMRES__HORIZON_HPP_