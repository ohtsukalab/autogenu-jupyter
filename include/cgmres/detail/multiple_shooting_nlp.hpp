#ifndef CGMRES__MULTIPLE_SHOOTING_NLP_HPP_
#define CGMRES__MULTIPLE_SHOOTING_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"

#include "cgmres/detail/control_input_bounds.hpp"
#include "cgmres/detail/control_input_bounds_shooting.hpp"

namespace cgmres {
namespace detail {

template <class OCP, int N>
class MultipleShootingNLP {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = OCP::nub;
  static constexpr int dim = nuc * N;

  MultipleShootingNLP(const OCP& ocp, const Horizon& horizon) 
    : ocp_(ocp),
      horizon_(horizon) {
    static_assert(OCP::nx > 0);
    static_assert(OCP::nu > 0);
    static_assert(OCP::nc >= 0);
    static_assert(OCP::nub >= 0);
    static_assert(N > 0);
  }

  MultipleShootingNLP() = default;

  ~MultipleShootingNLP() = default;

  template <typename VectorType>
  void eval_fonc_hu(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution,
                    const std::array<Vector<nx>, N+1>& x, const std::array<Vector<nx>, N+1>& lmd,
                    Vector<dim>& fonc_hu) {
    assert(x0.size());
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x0.derived().data(), solution.template head<nuc>().data(), lmd[1].data(), 
                 fonc_hu.template head<nuc>().data());
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_hu(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(),
                   lmd[i+1].data(), fonc_hu.template segment<nuc>(nuc*i).data());
    }
  }

  template <typename VectorType>
  void eval_fonc_f(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution,
                   const std::array<Vector<nx>, N+1>& x, 
                   std::array<Vector<nx>, N+1>& fonc_f) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute optimality error for state.
    ocp_.eval_f(t, x0.derived().data(), solution.template head<nuc>().data(), dx_.data());
    fonc_f[0] = x[1] - x0 - dt * dx_;
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_f(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(), dx_.data());
      fonc_f[i] = x[i+1] - x[i] - dt * dx_;
    }
  }

  template <typename VectorType>
  void retrive_x(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution,
                 std::array<Vector<nx>, N+1>& x,
                 const std::array<Vector<nx>, N+1>& fonc_f) {
    const Scalar T = horizon_.T(t);
    const Scalar dt = T / N;
    assert(T >= 0);
    // Compute optimality error for state.
    ocp_.eval_f(t, x0.derived().data(), solution.template head<nuc>().data(), dx_.data());
    x[1] = x0 + dt * dx_  + fonc_f[0];
    for (size_t i=1; i<N; ++i) {
      ocp_.eval_f(t+i*dt, x[i].data(), solution.template segment<nuc>(nuc*i).data(), dx_.data());
      x[i+1] = x[i] + dt * dx_ + fonc_f[i];
    }
  }

  template <typename VectorType>
  void eval_fonc_hx(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution,
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

  template <typename VectorType>
  void retrive_lmd(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution,
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

  void eval_fonc_hu(const Vector<dim>& solution,
                    const std::array<Vector<nub>, N>& dummy, 
                    const std::array<Vector<nub>, N>& mu,
                    Vector<dim>& fonc_hu) const {
    ubounds::eval_fonc_hu<OCP, N>(ocp_, solution, dummy, mu, fonc_hu);
  }

  void eval_fonc_hdummy(const Vector<dim>& solution,
                        const std::array<Vector<nub>, N>& dummy, 
                        const std::array<Vector<nub>, N>& mu,
                        std::array<Vector<nub>, N>& fonc_hdummy) const {
    ubounds::eval_fonc_hdummy<OCP, N>(ocp_, solution, dummy, mu, fonc_hdummy);
  }

  void eval_fonc_hmu(const Vector<dim>& solution,
                     const std::array<Vector<nub>, N>& dummy, 
                     const std::array<Vector<nub>, N>& mu,
                     std::array<Vector<nub>, N>& fonc_hmu) const {
    ubounds::eval_fonc_hmu<OCP, N>(ocp_, solution, dummy, mu, fonc_hmu);
  }

  static void multiply_hdummy_inv(const std::array<Vector<nub>, N>& dummy, 
                                  const std::array<Vector<nub>, N>& mu,
                                  const std::array<Vector<nub>, N>& fonc_hdummy,
                                  const std::array<Vector<nub>, N>& fonc_hmu,
                                  std::array<Vector<nub>, N>& fonc_hdummy_inv) {
    ubounds::multiply_hdummy_inv<OCP, N>(dummy, mu, fonc_hdummy, fonc_hmu,
                                         fonc_hdummy_inv);
  }

  static void multiply_hmu_inv(const std::array<Vector<nub>, N>& dummy, 
                               const std::array<Vector<nub>, N>& mu,
                               const std::array<Vector<nub>, N>& fonc_hdummy,
                               const std::array<Vector<nub>, N>& fonc_hmu,
                               const std::array<Vector<nub>, N>& fonc_hdummy_inv,
                               std::array<Vector<nub>, N>& fonc_hmu_inv) {
    ubounds::multiply_hmu_inv<OCP, N>(dummy, mu, fonc_hdummy, fonc_hmu,
                                      fonc_hdummy_inv, fonc_hmu_inv);
  }

  void retrive_dummy_update(const Vector<OCP::nuc*N>& solution,
                            const std::array<Vector<OCP::nub>, N>& dummy, 
                            const std::array<Vector<OCP::nub>, N>& mu,
                            const Vector<OCP::nuc*N>& solution_update,
                            std::array<Vector<OCP::nub>, N>& dummy_update) {
    ubounds::retrive_dummy_update<OCP, N>(ocp_, solution, dummy, mu, solution_update, dummy_update);
  }

  void retrive_mu_update(const Vector<OCP::nuc*N>& solution,
                         const std::array<Vector<OCP::nub>, N>& dummy, 
                         const std::array<Vector<OCP::nub>, N>& mu,
                         const Vector<OCP::nuc*N>& solution_update,
                         std::array<Vector<OCP::nub>, N>& mu_update) {
    ubounds::retrive_mu_update<OCP, N>(ocp_, solution, dummy, mu, solution_update, mu_update);
  }

  void clip_dummy(std::array<Vector<OCP::nub>, N>& dummy, const Scalar min) {
    ubounds::clip_dummy<OCP, N>(dummy, min);
  }

  void synchronize_ocp() { ocp_.synchronize(); }

  const OCP& ocp() const { return ocp_; }

  const Horizon& horizon() const { return horizon_; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  OCP ocp_;
  Horizon horizon_;
  Vector<nx> dx_;
};

} // namespace detail
} // namespace cgmres

#endif // CGMRES__MULTIPLE_SHOOTING_NLP_HPP_