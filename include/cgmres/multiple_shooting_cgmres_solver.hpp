#ifndef CGMRES__MULTIPLE_SHOOTING_CGMRES_SOLVER_HPP_
#define CGMRES__MULTIPLE_SHOOTING_CGMRES_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/matrixfree_gmres.hpp"
#include "cgmres/multiple_shooting_nlp.hpp"
#include "cgmres/continuation_gmres_condensing.hpp"
#include "cgmres/solver_settings.hpp"

namespace cgmres {

template <class OCP, int N, int kmax>
class MultipleShootingCGMRESSolver {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int dim = nuc * N;

  using MultipleShootingNLP_ = MultipleShootingNLP<OCP, N>;
  using ContinuationGMRES_ = ContinuationGMRESCondensing<MultipleShootingNLP_>;
  using MatrixFreeGMRES_ = MatrixFreeGMRES<kmax, ContinuationGMRES_>;

  MultipleShootingCGMRESSolver(const OCP& ocp, const Horizon& horizon, 
                               const SolverSettings& settings) 
    : nlp_(ocp, horizon),
      continuation_gmres_(nlp_, settings.finite_difference_epsilon, settings.zeta),
      gmres_(),
      settings_(settings),
      solution_(Vector<dim>::Zero()),
      solution_update_(Vector<dim>::Zero()) {
    std::fill(uopt_.begin(), uopt_.end(), Vector<nu>::Zero());
    std::fill(ucopt_.begin(), ucopt_.end(), Vector<nuc>::Zero());
    std::fill(xopt_.begin(), xopt_.end(), Vector<nx>::Zero());
    std::fill(lmdopt_.begin(), lmdopt_.end(), Vector<nx>::Zero());
  }

  MultipleShootingCGMRESSolver() = default;

  ~MultipleShootingCGMRESSolver() = default;

  void set_u(const Vector<nu>& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_u] u.size() must be " + std::to_string(nu));
    }
    for (size_t i=0; i<N; ++i) {
      uopt_[i] = u;
    }
    for (size_t i=0; i<N; ++i) {
      ucopt_[i].template head<nu>() = u;
      ucopt_[i].template tail<nc>().setZero();
    }
    setInnerSolution();
  }

  void set_uc(const Vector<nuc>& uc) {
    if (uc.size() != nuc) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_uc] uc.size() must be " + std::to_string(nuc));
    }
    for (size_t i=0; i<N; ++i) {
      uopt_[i] = uc.template head<nu>();
    }
    for (size_t i=0; i<N; ++i) {
      ucopt_[i] = uc;
    }
    setInnerSolution();
  }

  void set_x(const Vector<nx>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_x] x.size() must be " + std::to_string(nx));
    }
    for (size_t i=1; i<=N; ++i) {
      xopt_[i] = x;
    }
  }

  void set_lmd(const Vector<nx>& lmd) {
    if (lmd.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_lmd] lmd.size() must be " + std::to_string(nx));
    }
    for (size_t i=1; i<=N; ++i) {
      lmdopt_[i] = lmd;
    }
  }

  const std::array<Vector<nu>, N>& uopt() const { return uopt_; }

  const std::array<Vector<nuc>, N>& ucopt() const { return ucopt_; }

  const std::array<Vector<nx>, N+1>& xopt() const { return xopt_; }

  const std::array<Vector<nx>, N+1>& lmdopt() const { return lmdopt_; }

  Scalar optError() const { return continuation_gmres_.optError(); }

  Scalar optError(const Scalar t, const Vector<nx>& x) {
    continuation_gmres_.eval_fonc(t, x, solution_, xopt_, lmdopt_);
    return optError();
  }

  void update(const Scalar t, const Vector<nx>& x) {
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= update solution with C/GMRES =======================" << std::endl;
    }
    const auto gmres_iter 
        = gmres_.template solve<const Scalar, const Vector<nx>&, const Vector<dim>&,
                                const std::array<Vector<nx>, N+1>&, const std::array<Vector<nx>, N+1>&>(
              continuation_gmres_, t, x, solution_, xopt_, lmdopt_, solution_update_);
    const auto opt_error = continuation_gmres_.optError();

    // verbose
    if (settings_.verbose_level >= 1) {
      std::cout << "opt error: " << opt_error 
                << " (opt tol: " << settings_.opterr_tol << ")" <<  std::endl;
    }
    if (settings_.verbose_level >= 2) {
      std::cout << "number of GMRES iter: " << gmres_iter << " (kmax: " << kmax << ")" << std::endl;
    }

    continuation_gmres_.expansion(t, x, solution_, xopt_, lmdopt_, solution_update_, settings_.dt);
    solution_.noalias() += settings_.dt * solution_update_;

    retriveSolution();
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  MultipleShootingNLP_ nlp_;
  ContinuationGMRES_ continuation_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;

  std::array<Vector<nu>, N> uopt_;
  std::array<Vector<nuc>, N> ucopt_;
  std::array<Vector<nx>, N+1> xopt_;
  std::array<Vector<nx>, N+1> lmdopt_;

  Vector<dim> solution_, solution_update_; 

  void setInnerSolution() {
    for (size_t i=0; i<N; ++i) {
      solution_.template segment<nuc>(i*nuc) = ucopt_[i];
    }
  }

  void retriveSolution() {
    for (size_t i=0; i<N; ++i) {
      uopt_[i] = solution_.template segment<nu>(i*nuc);
      ucopt_[i] = solution_.template segment<nuc>(i*nuc);
    }
  }

};

} // namespace cgmres

#endif // CGMRES__MULTIPLE_SHOOTING_CGMRES_SOLVER_HPP_