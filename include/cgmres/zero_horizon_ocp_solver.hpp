#ifndef ZERO_HORIZON_OCP_SOLVER_HPP_
#define ZERO_HORIZON_OCP_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/matrixfree_gmres.hpp"
#include "cgmres/newton_gmres.hpp"
#include "cgmres/zero_horizon_nlp.hpp"
#include "cgmres/solver_settings.hpp"

namespace cgmres {

template <class OCP, int kmax>
class ZeroHorizonOCPSolver {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int dim = nuc;

  using ZeroHorizonNLP_ = ZeroHorizonNLP<OCP>;
  using NewtonGMRES_ = NewtonGMRES<ZeroHorizonNLP_>;
  using MatrixFreeGMRES_ = MatrixFreeGMRES<kmax, NewtonGMRES_>;

  ZeroHorizonOCPSolver(const OCP& ocp, const SolverSettings& settings) 
    : nlp_(ocp),
      newton_gmres_(nlp_, settings.finite_diference_epsilon),
      gmres_(),
      settings_(settings),
      uopt_(Vector<nu>::Zero()),
      ucopt_(Vector<nuc>::Zero()),
      solution_(Vector<dim>::Zero()),
      solution_update_(Vector<dim>::Zero()) {
  }

  ~ZeroHorizonOCPSolver() = default;

  template <typename VectorType>
  void set_u(const VectorType& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_u] u.size() must be " + std::to_string(nu));
    }
    uopt_ = u;
    ucopt_.template head<nu>() = u;
    ucopt_.template tail<nc>().setZero();
    setInnerSolution();
  }

  template <typename VectorType>
  void set_uc(const VectorType& uc) {
    if (uc.size() != nuc) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_uc] uc.size() must be " + std::to_string(nuc));
    }
    uopt_ = uc.template head<nu>();
    ucopt_ = uc;
    setInnerSolution();
  }

  const Vector<nu>& uopt() const { return uopt_; }

  const Vector<nuc>& ucopt() const { return ucopt_; }

  const Vector<nx>& lmdopt() const { return newton_gmres_.lmd(); }

  void solve(const Scalar t, const Vector<nx>& x) {
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= solve zero horizon OCP =======================" << std::endl;
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
        std::cout << "         number of GMRES iter: " << gmres_iter 
                  << " (kmax: " << kmax << ")" << std::endl;
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
  ZeroHorizonNLP_ nlp_;
  NewtonGMRES_ newton_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;

  Vector<nu> uopt_;
  Vector<nuc> ucopt_;

  Vector<dim> solution_, solution_update_; 

  void setInnerSolution() {
    solution_ = ucopt_;
  }

  void retriveSolution() {
    uopt_ = solution_.template head<nu>();
    ucopt_ = solution_;
  }

};

} // namespace cgmres

#endif // ZERO_HORIZON_OCP_SOLVER_HPP_