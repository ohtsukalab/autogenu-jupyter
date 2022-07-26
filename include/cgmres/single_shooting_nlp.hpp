#ifndef CGMRES__SINGLE_SHOOTING_NLP_HPP_
#define CGMRES__SINGLE_SHOOTING_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"

namespace cgmres {

template <class OCP, int N>
class SingleShootingNLP {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int dim = nuc * N;

  SingleShootingNLP(const OCP& ocp, const Horizon& horizon) 
    : ocp_(ocp),
      horizon_(horizon),
      dx_(Vector<nx>::Zero()) {
    assert(nx > 0);
    assert(nu > 0);
    assert(nc >= 0);
    assert(N > 0);
    std::fill(x_.begin(), x_.end(), Vector<nx>::Zero());
    std::fill(lmd_.begin(), lmd_.end(), Vector<nx>::Zero());
  }

  SingleShootingNLP() = default;

  ~SingleShootingNLP() = default;

  void eval(const Scalar t, const Vector<nx>& x, const Vector<dim>& solution,
            Vector<dim>& fonc) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    x_[0] = x;
    // Compute the state trajectory over the horizon  
    ocp_.eval_f(t, x_[0].data(), solution.template head<nuc>().data(), dx_.data());
    x_[1] = x_[0] + dt * dx_;
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_f(t+i*dt, x_[i].data(), solution.template segment<nuc>(nuc*i).data(), dx_.data());
      x_[i+1] = x_[i] + dt * dx_;
    }
    // Compute the Lagrange multiplier over the horizon  
    ocp_.eval_phix(t+T, x_[N].data(), lmd_[N].data());
    for (size_t i=N-1; i>=1; --i) {
      ocp_.eval_hx(t+i*dt, x_[i].data(), solution.template segment<nuc>(nuc*i).data(),
                   lmd_[i+1].data(), dx_.data());
      lmd_[i] = lmd_[i+1] + dt * dx_;
    }
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x_[0].data(), solution.template head<nuc>().data(), lmd_[1].data(), 
                 fonc.template head<nuc>().data());
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_hu(t+i*dt, x_[i].data(), solution.template segment<nuc>(nuc*i).data(),
                   lmd_[i+1].data(), fonc.template segment<nuc>(nuc*i).data());
    }
  }

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

} // namespace cgmres

#endif // CGMRES__SINGLE_SHOOTING_NLP_HPP_