#ifndef SINGLE_SHOOTING_OCP_SOLVER_HPP_
#define SINGLE_SHOOTING_OCP_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/matrixfree_gmres.hpp"
#include "cgmres/newton_gmres.hpp"
#include "cgmres/single_shooting_nlp.hpp"
#include "cgmres/solver_settings.hpp"

namespace cgmres {

template <class OCP, int N, int kmax>
class SingleShootingOCPSolver {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int dim = nuc * N;

  using SingleShootingNLP_ = SingleShootingNLP<OCP, N>;
  using NewtonGMRES_ = NewtonGMRES<SingleShootingNLP_>;
  using MatrixFreeGMRES_ = MatrixFreeGMRES<kmax, NewtonGMRES_>;

  SingleShootingOCPSolver(const OCP& ocp, const Horizon& horizon, 
                          const SolverSettings& settings) 
    : nlp_(ocp, horizon),
      newton_gmres_(nlp_, settings.finite_diference_epsilon),
      gmres_(),
      settings_(settings),
      solution_(Vector<dim>::Zero()),
      solution_update_(Vector<dim>::Zero()) {
    std::fill(uopt_.begin(), uopt_.end(), Vector<nu>::Zero());
    std::fill(ucopt_.begin(), ucopt_.end(), Vector<nuc>::Zero());
  }

  ~SingleShootingOCPSolver() = default;

  template <typename VectorType>
  void set_u(const VectorType& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[SingleShootingOCPSolver::set_u] u.size() must be " + std::to_string(nu));
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

  template <typename VectorType>
  void set_uc(const VectorType& uc) {
    if (uc.size() != nuc) {
      throw std::invalid_argument("[SingleShootingOCPSolver::set_uc] uc.size() must be " + std::to_string(nuc));
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

  const std::array<Vector<nx>, N+1>& xopt() const { return newton_gmres_.x(); }

  const std::array<Vector<nx>, N+1>& lmdopt() const { return newton_gmres_.lmd(); }

  void solve(const Scalar t, const Vector<nx>& x) {
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= solve single-shooting OCP =======================" << std::endl;
    }

    for (size_t iter=0; iter<settings_.max_iter; ++iter) {
      const auto gmres_iter 
          = gmres_.template solve<const Scalar, const Vector<nx>&, const Vector<dim>&>(
                newton_gmres_, t, x, solution_, solution_update_);
      const auto opt_error = newton_gmres_.optError();

      // verbose
      if (settings_.verbose_level >= 1) {
        std::cout << "iter " << iter << ": opt error: " << opt_error 
                  << " (opt tol: " << settings_.opt_error_tol << ")" <<  std::endl;
      }
      if (settings_.verbose_level >= 2) {
        std::cout << "         number of GMRES iter: " << gmres_iter << " (kmax: " << kmax << ")" << std::endl;
      }

      solution_.noalias() += solution_update_;
      if (opt_error < settings_.opt_error_tol) {
        if (settings_.verbose_level >= 1) {
          std::cout << "converged!" << std::endl;
        }
        break;
      }
    }

    retriveSolution();
  }

private:
  SingleShootingNLP_ nlp_;
  NewtonGMRES_ newton_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;

  std::array<Vector<nu>, N+1> uopt_;
  std::array<Vector<nuc>, N+1> ucopt_;

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

#endif // SINGLE_SHOOTING_OCP_SOLVER_HPP_