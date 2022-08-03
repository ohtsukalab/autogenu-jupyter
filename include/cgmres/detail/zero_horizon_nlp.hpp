#ifndef CGMRES__ZERO_HORIZON_NLP_HPP_
#define CGMRES__ZERO_HORIZON_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"

#include "cgmres/detail/control_input_bounds.hpp"

namespace cgmres {
namespace detail {

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
    static_assert(OCP::nx > 0);
    static_assert(OCP::nu > 0);
    static_assert(OCP::nc >= 0);
    static_assert(OCP::nub >= 0);
  }

  ZeroHorizonNLP() = default;

  ~ZeroHorizonNLP() = default;

  template <typename VectorType>
  void eval_fonc_hu(const Scalar t, const MatrixBase<VectorType>& x, const Vector<dim>& solution,
                    Vector<dim>& fonc_hu) {
    // Compute the Lagrange multiplier over the horizon  
    ocp_.eval_phix(t, x.derived().data(), lmd_.data());
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x.derived().data(), solution.data(), lmd_.data(), fonc_hu.data());
    if constexpr (nub > 0) {
      const auto uc    = solution.template head<nuc>();
      const auto dummy = solution.template segment<nub>(nuc);
      const auto mu    = solution.template segment<nub>(nuc+nub);
      auto fonc_huc    = fonc_hu.template head<nuc>();
      auto fonc_hdummy = fonc_hu.template segment<nub>(nuc);
      auto fonc_hmu    = fonc_hu.template segment<nub>(nuc+nub);
      ubounds::eval_hu(ocp_, uc, dummy, mu, fonc_huc);
      ubounds::eval_hdummy(ocp_, uc, dummy, mu, fonc_hdummy);
      ubounds::eval_hmu(ocp_, uc, dummy, mu, fonc_hmu);
    }
  }

  void retrive_dummy(Vector<dim>& solution, Vector<dim>& fonc_hu, const Scalar min_dummy) {
    if constexpr (nub > 0) {
      const auto uc    = solution.template head<nuc>();
      auto dummy = solution.template segment<nub>(nuc);
      const auto mu    = solution.template segment<nub>(nuc+nub);
      auto fonc_hmu    = fonc_hu.template segment<nub>(nuc+nub);
      dummy.setZero();
      ubounds::eval_hmu(ocp_, uc, dummy, mu, fonc_hmu);
      ubounds::clip_dummy(dummy, min_dummy);
      dummy.array() = fonc_hmu.array().abs().sqrt();
    }
  }

  void retrive_mu(Vector<dim>& solution, Vector<dim>& fonc_hu) {
    if constexpr (nub > 0) {
      const auto uc    = solution.template head<nuc>();
      const auto dummy = solution.template segment<nub>(nuc);
      auto mu    = solution.template segment<nub>(nuc+nub);
      auto fonc_hdummy = fonc_hu.template segment<nub>(nuc);
      mu.setZero();
      ubounds::eval_hmu(ocp_, uc, dummy, mu, fonc_hdummy);
      mu.array() = fonc_hdummy.array() / (2.0 * dummy.array());
    }
  }

  void synchronize_ocp() { ocp_.synchronize(); }

  const OCP& ocp() const { return ocp_; }

  const Vector<nx>& lmd() const { return lmd_; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  OCP ocp_;
  Vector<nx> lmd_;
};

} // namespace detail
} // namespace cgmres

#endif // CGMRES__ZERO_HORIZON_NLP_HPP_