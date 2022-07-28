#ifndef CGMRES__HORIZON_HPP_
#define CGMRES__HORIZON_HPP_

#include <cmath>
#include <stdexcept>
#include <cassert>
#include <iostream>

#include "cgmres/types.hpp"


namespace cgmres {

class Horizon {
public:
  Horizon(const Scalar Tf, const Scalar alpha=0.0, const Scalar t0=0.0)
    : Tf_(Tf), alpha_(alpha), t0_(t0) {
    if (Tf <= 0.0) {
      throw std::invalid_argument("[Horizon]: 'Tf' must be positive!");
    }
    time_varying_length_ = (alpha > 0.0);
  }

  Horizon() = default;

  ~Horizon() = default;

  inline Scalar T(const Scalar t) const {
    if (time_varying_length_) {
      assert(t >= t0_);
      return Tf_ * (1.0-std::exp(-alpha_*(t-t0_)));
    }
    else {
      return Tf_;
    }
  }

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
    os << "  t0:    " << t0_ << std::flush;
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

#endif // CGMRES__HORIZON_H