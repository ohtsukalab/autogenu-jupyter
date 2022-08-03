#ifndef CGMRES__SINGLE_SHOOTING_NLP_HPP_
#define CGMRES__SINGLE_SHOOTING_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"

#include "cgmres/detail/control_input_bounds.hpp"
#include "cgmres/detail/control_input_bounds_shooting.hpp"

namespace cgmres {
namespace detail {

template <class OCP, int N>
class SingleShootingNLP {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = OCP::nub;
  static constexpr int dim = nuc * N + 2 * N * nub;

  SingleShootingNLP(const OCP& ocp, const Horizon& horizon) 
    : ocp_(ocp),
      horizon_(horizon),
      dx_(Vector<nx>::Zero()) {
    static_assert(OCP::nx > 0);
    static_assert(OCP::nu > 0);
    static_assert(OCP::nc >= 0);
    static_assert(OCP::nub >= 0);
    static_assert(N > 0);
    std::fill(x_.begin(), x_.end(), Vector<nx>::Zero());
    std::fill(lmd_.begin(), lmd_.end(), Vector<nx>::Zero());
  }

  SingleShootingNLP() = default;

  ~SingleShootingNLP() = default;

  template <typename VectorType>
  void eval_fonc_hu(const Scalar t, const MatrixBase<VectorType>& x, const Vector<dim>& solution,
                    Vector<dim>& fonc_hu) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    x_[0] = x;
    // Compute the state trajectory over the horizon  
    ocp_.eval_f(t, x_[0].data(), solution.template head<nuc>().data(), dx_.data());
    x_[1] = x_[0] + dt * dx_;
    for (size_t i=1; i<N; ++i) {
      const int inucb2 = i * (nuc + 2 * nub);
      ocp_.eval_f(t+i*dt, x_[i].data(), solution.template segment<nuc>(inucb2).data(), dx_.data());
      x_[i+1] = x_[i] + dt * dx_;
    }
    // Compute the Lagrange multiplier over the horizon  
    ocp_.eval_phix(t+T, x_[N].data(), lmd_[N].data());
    for (size_t i=N-1; i>=1; --i) {
      const int inucb2 = i * (nuc + 2 * nub);
      ocp_.eval_hx(t+i*dt, x_[i].data(), solution.template segment<nuc>(inucb2).data(),
                   lmd_[i+1].data(), dx_.data());
      lmd_[i] = lmd_[i+1] + dt * dx_;
    }
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x_[0].data(), solution.template head<nuc>().data(), lmd_[1].data(), 
                 fonc_hu.template head<nuc>().data());
    for (size_t i=1; i<N; ++i) {
      const int inucb2 = i * (nuc + 2 * nub);
      ocp_.eval_hu(t+i*dt, x_[i].data(), solution.template segment<nuc>(inucb2).data(),
                   lmd_[i+1].data(), fonc_hu.template segment<nuc>(inucb2).data());
    }
    if constexpr (nub > 0) {
      for (size_t i=0; i<N; ++i) {
        const int inucb2 = i * (nuc + 2 * nub);
        const auto uc    = solution.template segment<nuc>(inucb2);
        const auto dummy = solution.template segment<nub>(inucb2+nuc);
        const auto mu    = solution.template segment<nub>(inucb2+nuc+nub);
        auto fonc_huc    = fonc_hu.template segment<nuc>(inucb2);
        auto fonc_hdummy = fonc_hu.template segment<nub>(inucb2+nuc);
        auto fonc_hmu    = fonc_hu.template segment<nub>(inucb2+nuc+nub);
        ubounds::eval_hu(ocp_, uc, dummy, mu, fonc_huc);
        ubounds::eval_hdummy(ocp_, uc, dummy, mu, fonc_hdummy);
        ubounds::eval_hmu(ocp_, uc, dummy, mu, fonc_hmu);
      }
    }
  }

  void retrive_dummy(Vector<dim>& solution, Vector<dim>& fonc_hu, const Scalar min_dummy) {
    if constexpr (nub > 0) {
      for (size_t i=0; i<N; ++i) {
        const int inucb2 = i * (nuc + 2 * nub);
        const auto uc    = solution.template segment<nuc>(inucb2);
        auto dummy = solution.template segment<nub>(inucb2+nuc);
        const auto mu    = solution.template segment<nub>(inucb2+nuc+nub);
        auto fonc_hmu    = fonc_hu.template segment<nub>(inucb2+nuc+nub);
        ubounds::eval_hmu(ocp_, uc, dummy, mu, fonc_hmu);
        ubounds::clip_dummy(dummy, min_dummy);
        dummy.array() = fonc_hmu.array().abs().sqrt();
      }
    }
  }

  void retrive_mu(Vector<dim>& solution, Vector<dim>& fonc_hu) {
    if constexpr (nub > 0) {
      for (size_t i=0; i<N; ++i) {
        const int inucb2 = i * (nuc + 2 * nub);
        const auto uc    = solution.template segment<nuc>(inucb2);
        const auto dummy = solution.template segment<nub>(inucb2+nuc);
        auto mu    = solution.template segment<nub>(inucb2+nuc+nub);
        auto fonc_hdummy = fonc_hu.template segment<nub>(inucb2+nuc);
        mu.setZero();
        ubounds::eval_hdummy(ocp_, uc, dummy, mu, fonc_hdummy);
        mu.array() = - fonc_hdummy.array() / (2.0 * dummy.array());
      }
    }
  }

  void synchronize_ocp() { ocp_.synchronize(); }

  const OCP& ocp() const { return ocp_; }

  const Horizon& horizon() const { return horizon_; }

  const std::array<Vector<nx>, N+1>& x() const { return x_; }

  const std::array<Vector<nx>, N+1>& lmd() const { return lmd_; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  OCP ocp_;
  Horizon horizon_;
  Vector<nx> dx_;
  std::array<Vector<nx>, N+1> x_, lmd_;
};

} // namespace detail
} // namespace cgmres

#endif // CGMRES__SINGLE_SHOOTING_NLP_HPP_