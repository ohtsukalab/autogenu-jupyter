#ifndef HORIZON_HPP_
#define HORIZON_HPP_

#include <cmath>
#include <stdexcept>

#include "cgmres/types.hpp"


namespace cgmres {

class Horizon {
public:
  Horizon(const Scalar Tf, const Scalar alpha, const Scalar t0=0.0)
    : Tf_(Tf), alpha_(alpha), t0_(t0) {
    if (Tf <= 0.0) {
      throw std::invalid_argument("[Horizon]: 'Tf' must be positive!");
    }
  }

  ~Horizon() = default;

  inline Scalar T(const Scalar t) const {
    return Tf_ * (1.0-std::exp(-alpha_*(t-t0_)));
  }

  void reset(const Scalar t0) {
    t0_ = t0;
  }

private:
  Scalar Tf_, alpha_, t0_;
};

} // namespace cgmres

#endif // HORIZON_H