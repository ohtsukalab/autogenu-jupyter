#ifndef CGMRES__ZERO_HORIZON_OCP_SOLVER_HPP_
#define CGMRES__ZERO_HORIZON_OCP_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/solver_settings.hpp"

#include "cgmres/detail/matrixfree_gmres.hpp"
#include "cgmres/detail/zero_horizon_nlp.hpp"
#include "cgmres/detail/newton_gmres.hpp"

namespace cgmres {

template <class OCP, int kmax>
class ZeroHorizonOCPSolver {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = OCP::nub;
  static constexpr int dim = nuc + 2 * nub;

  using ZeroHorizonNLP_ = detail::ZeroHorizonNLP<OCP>;
  using NewtonGMRES_ = detail::NewtonGMRES<ZeroHorizonNLP_>;
  using MatrixFreeGMRES_ = detail::MatrixFreeGMRES<NewtonGMRES_, kmax>;

  ZeroHorizonOCPSolver(const OCP& ocp, const SolverSettings& settings) 
    : nlp_(ocp),
      newton_gmres_(nlp_, settings.finite_difference_epsilon),
      gmres_(),
      settings_(settings),
      uopt_(Vector<nu>::Zero()),
      ucopt_(Vector<nuc>::Zero()),
      solution_(Vector<dim>::Zero()),
      solution_update_(Vector<dim>::Zero()) {
  }

  ZeroHorizonOCPSolver() = default;

  ~ZeroHorizonOCPSolver() = default;

  template <typename VectorType>
  void set_u(const VectorType& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_u] u.size() must be " + std::to_string(nu));
    }
    uopt_ = u;
    ucopt_.template head<nu>() = u;
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

  template <typename VectorType>
  void set_dummy(const MatrixBase<VectorType>& dummy) {
    if (dummy.size() != nub) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_dummy] dummy.size() must be " + std::to_string(nub));
    }
    dummyopt_ = dummy;
    setInnerSolution();
  }

  template <typename VectorType>
  void set_mu(const MatrixBase<VectorType>& mu) {
    if (mu.size() != nub) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_mu] mu.size() must be " + std::to_string(nub));
    }
    muopt_ = mu;
    setInnerSolution();
  }

  const Vector<nu>& uopt() const { return uopt_; }

  const Vector<nuc>& ucopt() const { return ucopt_; }

  const Vector<nub>& dummyopt() const { return dummyopt_; }

  const Vector<nub>& muopt() const { return muopt_; }

  const Vector<nx>& lmdopt() const { return newton_gmres_.lmd(); }

  Scalar optError() const { return newton_gmres_.optError(); }

  template <typename VectorType>
  Scalar optError(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::optError] x.size() must be " + std::to_string(nx));
    }
    newton_gmres_.eval_fonc(t, x, solution_);
    return optError();
  }

  template <typename VectorType>
  void solve(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::update] x.size() must be " + std::to_string(nx));
    }
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= solve zero horizon OCP =======================" << std::endl;
    }

    for (size_t iter=0; iter<settings_.max_iter; ++iter) {
      const auto gmres_iter 
          = gmres_.template solve<const Scalar, const VectorType&, const Vector<dim>&>(
                newton_gmres_, t, x.derived(), solution_, solution_update_);
      const auto opt_error = newton_gmres_.optError();

      // verbose
      if (settings_.verbose_level >= 1) {
        std::cout << "iter " << iter << ": opt error: " << opt_error 
                  << " (opt tol: " << settings_.opterr_tol << ")" <<  std::endl;
      }
      if (settings_.verbose_level >= 2) {
        std::cout << "         number of GMRES iter: " << gmres_iter 
                  << " (kmax: " << kmax << ")" << std::endl;
      }

      solution_.noalias() += solution_update_;
      if (opt_error < settings_.opterr_tol) {
        if (settings_.verbose_level >= 1) {
          std::cout << "converged!" << std::endl;
        }
        break;
      }

    }
    retriveSolution();
  }

  void init_dummy_mu() {
    newton_gmres_.retrive_dummy(solution_, settings_.min_dummy);
    newton_gmres_.retrive_mu(solution_);
    retriveSolution();
  }

  void disp(std::ostream& os) const {
    os << "Zero horizon OCP solver: " << std::endl;
    os << "  kmax: " << kmax << std::endl;
    os << nlp_.ocp() << std::endl;
    os << settings_ << std::flush;
  }

  friend std::ostream& operator<<(std::ostream& os, const ZeroHorizonOCPSolver& solver) {
    solver.disp(os);
    return os;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ZeroHorizonNLP_ nlp_;
  NewtonGMRES_ newton_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;

  Vector<nu> uopt_;
  Vector<nuc> ucopt_;
  Vector<nub> dummyopt_, muopt_;

  Vector<dim> solution_, solution_update_; 

  void setInnerSolution() {
    solution_.template head<nuc>() = ucopt_;
    if constexpr (nub > 0) {
      solution_.template segment<nub>(nuc) = dummyopt_;
      solution_.template segment<nub>(nuc+nub) = muopt_;
    }
  }

  void retriveSolution() {
    uopt_ = solution_.template head<nu>();
    ucopt_ = solution_.template head<nuc>();
    dummyopt_ = solution_.template segment<nub>(nuc);
    muopt_ = solution_.template segment<nub>(nuc+nub);
  }

};

} // namespace cgmres

#endif // CGMRES__ZERO_HORIZON_OCP_SOLVER_HPP_