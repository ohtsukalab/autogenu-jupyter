#ifndef CGMRES__SINGLE_SHOOTING_CGMRES_SOLVER_HPP_
#define CGMRES__SINGLE_SHOOTING_CGMRES_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/matrixfree_gmres.hpp"
#include "cgmres/single_shooting_nlp.hpp"
#include "cgmres/continuation_gmres.hpp"
#include "cgmres/solver_settings.hpp"

namespace cgmres {

template <class OCP, int N, int kmax>
class SingleShootingCGMRESSolver {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int dim = nuc * N;

  using SingleShootingNLP_ = SingleShootingNLP<OCP, N>;
  using ContinuationGMRES_ = ContinuationGMRES<SingleShootingNLP_>;
  using MatrixFreeGMRES_ = MatrixFreeGMRES<kmax, ContinuationGMRES_>;

  SingleShootingCGMRESSolver(const OCP& ocp, const Horizon& horizon, 
                             const SolverSettings& settings) 
    : nlp_(ocp, horizon),
      continuation_gmres_(nlp_, settings.finite_difference_epsilon, settings.zeta),
      gmres_(),
      settings_(settings),
      solution_(Vector<dim>::Zero()),
      solution_update_(Vector<dim>::Zero()) {
    std::fill(uopt_.begin(), uopt_.end(), Vector<nu>::Zero());
    std::fill(ucopt_.begin(), ucopt_.end(), Vector<nuc>::Zero());
  }

  SingleShootingCGMRESSolver() = default;

  ~SingleShootingCGMRESSolver() = default;

  void set_u(const Vector<nu>& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_u] u.size() must be " + std::to_string(nu));
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
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_uc] uc.size() must be " + std::to_string(nuc));
    }
    for (size_t i=0; i<N; ++i) {
      uopt_[i] = uc.template head<nu>();
    }
    for (size_t i=0; i<N; ++i) {
      ucopt_[i] = uc;
    }
    setInnerSolution();
  }

  const std::array<Vector<nu>, N>& uopt() const { return uopt_; }

  const std::array<Vector<nuc>, N>& ucopt() const { return ucopt_; }

  const std::array<Vector<nx>, N+1>& xopt() const { return continuation_gmres_.x(); }

  const std::array<Vector<nx>, N+1>& lmdopt() const { return continuation_gmres_.lmd(); }

  Scalar optError() const { return continuation_gmres_.optError(); }

  void update(const Scalar t, const Vector<nx>& x) {
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= update solution with C/GMRES =======================" << std::endl;
    }
    const auto gmres_iter 
        = gmres_.template solve<const Scalar, const Vector<nx>&, const Vector<dim>&>(
              continuation_gmres_, t, x, solution_, solution_update_);
    const auto opt_error = continuation_gmres_.optError();

    // verbose
    if (settings_.verbose_level >= 1) {
      std::cout << "opt error: " << opt_error 
                << " (opt tol: " << settings_.opterr_tol << ")" <<  std::endl;
    }
    if (settings_.verbose_level >= 2) {
      std::cout << "number of GMRES iter: " << gmres_iter << " (kmax: " << kmax << ")" << std::endl;
    }

    solution_.noalias() += settings_.dt * solution_update_;

    retriveSolution();
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SingleShootingNLP_ nlp_;
  ContinuationGMRES_ continuation_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;

  std::array<Vector<nu>, N> uopt_;
  std::array<Vector<nuc>, N> ucopt_;

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

#endif // CGMRES__SINGLE_SHOOTING_CGMRES_SOLVER_HPP_