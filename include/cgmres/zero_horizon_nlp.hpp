#ifndef CGMRES__ZERO_HORIZON_NLP_HPP_
#define CGMRES__ZERO_HORIZON_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"
#include "cgmres/control_input_bounds.hpp"

namespace cgmres {

template <class OCP>
class ZeroHorizonNLP {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = OCP::nub;
  static constexpr int dim = nuc + 2 * nub;

  ZeroHorizonNLP(const OCP& ocp) 
    : ocp_(ocp),
      lmd_(Vector<nx>::Zero()) {
    assert(nx > 0);
    assert(nu > 0);
    assert(nc >= 0);
  }

  ZeroHorizonNLP() = default;

  ~ZeroHorizonNLP() = default;

  void eval_fonc_hu(const Scalar t, const Vector<nx>& x, const Vector<dim>& solution,
                    Vector<dim>& fonc_hu) {
    // Compute the Lagrange multiplier over the horizon  
    ocp_.eval_phix(t, x.data(), lmd_.data());
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x.data(), solution.data(), lmd_.data(), fonc_hu.data());

    if constexpr(nub > 0) {
      const auto uc = solution.template head<nuc>();
      const auto dummy = solution.template segment<nub>(nu);
      const auto mu = solution.template segment<nub>(nuc+nub);
      auto hu = fonc_hu.template head<nuc>();
      auto hdummy = fonc_hu.template segment<nub>(nuc);
      auto hmu = fonc_hu.template segment<nub>(nuc+nub);
      ubounds::eval_hu(ocp_, uc, dummy, mu, hu);
      ubounds::eval_hdummy(ocp_, uc, dummy, mu, hdummy);
      ubounds::eval_hmu(ocp_, uc, dummy, mu, hmu);
    }
  }

  const OCP& ocp() const { return ocp_; }

  const Vector<nx>& lmd() const { return lmd_; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  OCP ocp_;
  Vector<nx> lmd_;
};

} // namespace cgmres

#endif // CGMRES__ZERO_HORIZON_NLP_HPP_