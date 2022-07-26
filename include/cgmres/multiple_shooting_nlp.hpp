#ifndef MULTIPLE_SHOOTING_NLP_HPP_
#define MULTIPLE_SHOOTING_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"

namespace cgmres {

template <class OCP, int N>
class MultipleShootingNLP {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int dim = nuc * N;

  MultipleShootingNLP(const OCP& ocp, const Horizon& horizon) 
    : ocp_(ocp),
      horizon_(horizon) {
    assert(nx > 0);
    assert(nu > 0);
    assert(nc >= 0);
    assert(N > 0);
  }

  ~MultipleShootingNLP() = default;

  void eval_fonc_hu(const Scalar t, const Vector<nx>& x0, const Vector<dim>& solution,
                    const std::array<Vector<nx>, N+1>& x, const std::array<Vector<nx>, N+1>& lmd,
                    Vector<dim>& fonc_hu) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x0.data(), solution.template head<nuc>().data(), lmd[1].data(), 
                 fonc_hu.template head<nuc>().data());
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_hu(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(),
                   lmd[i+1].data(), fonc_hu.template segment<nuc>(nuc*i).data());
    }
  }

  void eval_fonc_f(const Scalar t, const Vector<nx>& x0, const Vector<dim>& solution,
                   const std::array<Vector<nx>, N+1>& x, 
                   std::array<Vector<nx>, N+1>& fonc_f) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute optimality error for state.
    ocp_.eval_f(t, x0.data(), solution.template head<nuc>().data(), dx_.data());
    fonc_f[0] = x[1] - x0 - dt * dx_;
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_f(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(), dx_.data());
      fonc_f[i] = x[i+1] - x[i] - dt * dx_;
    }
  }

  void retrive_x(const Scalar t, const Vector<nx>& x0, const Vector<dim>& solution,
                 std::array<Vector<nx>, N+1>& x,
                 const std::array<Vector<nx>, N+1>& fonc_f) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute optimality error for state.
    ocp_.eval_f(t, x0.data(), solution.template head<nuc>().data(), dx_.data());
    x[1] = x0 + dt * dx_  + fonc_f[0];
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_f(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(), dx_.data());
      x[i+1] = x[i] + dt * dx_ + fonc_f[i];
    }
  }

  void eval_fonc_hx(const Scalar t, const Vector<nx>& x0, const Vector<dim>& solution,
                    const std::array<Vector<nx>, N+1>& x, const std::array<Vector<nx>, N+1>& lmd,
                    std::array<Vector<nx>, N+1>& fonc_hx) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute optimality error for lambda.
    ocp_.eval_phix(t+T, x[N].data(), dx_.data());
    fonc_hx[N] = lmd[N] - dx_;
    for (size_t i=N-1; i>=1; --i) {
      ocp_.eval_hx(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(), 
                   lmd[i+1].data(), dx_.data());
      fonc_hx[i] = lmd[i] - lmd[i+1] - dt * dx_;
    }
  }

  void retrive_lmd(const Scalar t, const Vector<nx>& x0, const Vector<dim>& solution,
                   const std::array<Vector<nx>, N+1>& x, std::array<Vector<nx>, N+1>& lmd,
                   const std::array<Vector<nx>, N+1>& fonc_hx) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute optimality error for state.
    ocp_.eval_phix(t+T, x[N].data(), dx_.data());
    lmd[N] = dx_ + fonc_hx[N];
    for (size_t i=N-1; i>=1; --i) {
      ocp_.eval_hx(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(), 
                   lmd[i+1].data(), dx_.data());
      lmd[i] = lmd[i+1] + dt * dx_ + fonc_hx[i];
    }
  }

  const OCP& ocp() const { return ocp_; }

  const Horizon& horizon() const { return horizon_; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  OCP ocp_;
  Horizon horizon_;
  Vector<nx> dx_;
};

} // namespace cgmres

#endif // MULTIPLE_SHOOTING_NLP_HPP_