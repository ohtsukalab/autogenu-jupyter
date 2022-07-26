#ifndef ZERO_HORIZON_NLP_HPP_
#define ZERO_HORIZON_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"

namespace cgmres {

template <class OCP>
class ZeroHorizonNLP {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int dim = nuc;

  ZeroHorizonNLP(const OCP& ocp) 
    : ocp_(ocp),
      lmd_(Vector<nx>::Zero()) {
    assert(nx > 0);
    assert(nu > 0);
    assert(nc >= 0);
  }

  ZeroHorizonNLP() = default;

  ~ZeroHorizonNLP() = default;

  void eval(const Scalar t, const Vector<nx>& x, const Vector<dim>& solution,
            Vector<dim>& fonc) {
    // Compute the Lagrange multiplier over the horizon  
    ocp_.eval_phix(t, x.data(), lmd_.data());
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x.data(), solution.data(), lmd_.data(), fonc.data());
  }

  const OCP& ocp() const { return ocp_; }

  const Vector<nx>& lmd() const { return lmd_; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  OCP ocp_;
  Vector<nx> lmd_;
};

} // namespace cgmres

#endif // ZERO_HORIZON_NLP_HPP_