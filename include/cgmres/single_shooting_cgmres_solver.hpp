#ifndef CGMRES__SINGLE_SHOOTING_CGMRES_SOLVER_HPP_
#define CGMRES__SINGLE_SHOOTING_CGMRES_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/solver_settings.hpp"

#include "cgmres/detail/matrixfree_gmres.hpp"
#include "cgmres/detail/single_shooting_nlp.hpp"
#include "cgmres/detail/continuation_gmres.hpp"

namespace cgmres {

template <class OCP, int N, int kmax>
class SingleShootingCGMRESSolver {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = OCP::nub;
  static constexpr int dim = nuc * N + 2 * N * nub;

  using SingleShootingNLP_ = detail::SingleShootingNLP<OCP, N>;
  using ContinuationGMRES_ = detail::ContinuationGMRES<SingleShootingNLP_>;
  using MatrixFreeGMRES_ = detail::MatrixFreeGMRES<ContinuationGMRES_, kmax>;

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
    if constexpr (nub > 0) {
      std::fill(dummyopt_.begin(), dummyopt_.end(), Vector<nub>::Zero());
      std::fill(muopt_.begin(), muopt_.end(), Vector<nub>::Zero());
    }
  }

  SingleShootingCGMRESSolver() = default;

  ~SingleShootingCGMRESSolver() = default;

  template <typename VectorType>
  void set_u(const MatrixBase<VectorType>& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_u] u.size() must be " + std::to_string(nu));
    }
    for (size_t i=0; i<N; ++i) {
      uopt_[i] = u;
    }
    for (size_t i=0; i<N; ++i) {
      ucopt_[i].template head<nu>() = u;
    }
    setInnerSolution();
  }

  template <typename VectorType>
  void set_uc(const MatrixBase<VectorType>& uc) {
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

  template <typename VectorType>
  void set_dummy(const MatrixBase<VectorType>& dummy) {
    if (dummy.size() != nub) {
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_dummy] dummy.size() must be " + std::to_string(nub));
    }
    for (size_t i=0; i<N; ++i) {
      dummyopt_[i] = dummy;
    }
    setInnerSolution();
  }

  template <typename VectorType>
  void set_mu(const MatrixBase<VectorType>& mu) {
    if (mu.size() != nub) {
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_mu] mu.size() must be " + std::to_string(nub));
    }
    for (size_t i=0; i<N; ++i) {
      muopt_[i] = mu;
    }
    setInnerSolution();
  }

  template <typename T>
  void set_u_array(const T& u_array) {
    if (u_array.size() != N) { 
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_u_array]: 'u_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (u_array[i].size() != nu) {
        throw std::invalid_argument("[SingleShootingCGMRESSolver::set_u_array] u_array[i].size() must be " + std::to_string(nu));
      }
      uopt_[i] = u_array[i];
      ucopt_[i].template head<nu>() = u_array[i];
    }
    setInnerSolution();
  }

  template <typename T>
  void set_uc_array(const T& uc_array) {
    if (uc_array.size() != N) { 
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_uc_array]: 'uc_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (uc_array[i].size() != nuc) {
        throw std::invalid_argument("[SingleShootingCGMRESSolver::uc_array] uc_array[i].size() must be " + std::to_string(nuc));
      }
      uopt_[i] = uc_array[i].template head<nu>();
      ucopt_[i] = uc_array[i];
    }
    setInnerSolution();
  }

  template <typename T>
  void set_dummy_array(const T& dummy_array) {
    if (dummy_array.size() != N) { 
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_dummy_array]: 'dummy_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (dummy_array[i].size() != nub) {
        throw std::invalid_argument("[SingleShootingCGMRESSolver::set_dummy_array] dummy_array[i].size() must be " + std::to_string(nub));
      }
      dummyopt_[i] = dummy_array[i];
    }
    setInnerSolution();
  }

  template <typename T>
  void set_mu_array(const T& mu_array) {
    if (mu_array.size() != N) { 
      throw std::invalid_argument("[SingleShootingCGMRESSolver::set_mu_array]: 'mu_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (mu_array[i].size() != nub) {
        throw std::invalid_argument("[SingleShootingCGMRESSolver::set_mu_array] mu_array[i].size() must be " + std::to_string(nub));
      }
      muopt_[i] = mu_array[i];
    }
    setInnerSolution();
  }

  void init_dummy_mu() {
    continuation_gmres_.retrive_dummy(solution_, settings_.min_dummy);
    continuation_gmres_.retrive_mu(solution_);
    retriveSolution();
  }

  const std::array<Vector<nu>, N>& uopt() const { return uopt_; }

  const std::array<Vector<nuc>, N>& ucopt() const { return ucopt_; }

  const std::array<Vector<nx>, N+1>& xopt() const { return continuation_gmres_.x(); }

  const std::array<Vector<nx>, N+1>& lmdopt() const { return continuation_gmres_.lmd(); }

  const std::array<Vector<nub>, N>& dummyopt() const { return dummyopt_; }

  const std::array<Vector<nub>, N>& muopt() const { return muopt_; }

  Scalar optError() const { return continuation_gmres_.optError(); }

  template <typename VectorType>
  Scalar optError(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[SingleShootingCGMRESSolver::optError] x.size() must be " + std::to_string(nx));
    }
    continuation_gmres_.eval_fonc(t, x, solution_);
    return optError();
  }

  template <typename VectorType>
  void update(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[SingleShootingCGMRESSolver::update] x.size() must be " + std::to_string(nx));
    }
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= update solution with C/GMRES =======================" << std::endl;
    }
    const auto gmres_iter 
        = gmres_.template solve<const Scalar, const VectorType&, const Vector<dim>&>(
              continuation_gmres_, t, x.derived(), solution_, solution_update_);
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

  void disp(std::ostream& os) const {
    os << "Single shooting CGMRES solver: " << std::endl;
    os << "  N:    " << N << std::endl;
    os << "  kmax: " << kmax << std::endl;
    os << nlp_.ocp() << std::endl;
    os << nlp_.horizon() << std::endl;
    os << settings_ << std::flush;
  }

  friend std::ostream& operator<<(std::ostream& os, const SingleShootingCGMRESSolver& solver) {
    solver.disp(os);
    return os;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SingleShootingNLP_ nlp_;
  ContinuationGMRES_ continuation_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;

  std::array<Vector<nu>, N> uopt_;
  std::array<Vector<nuc>, N> ucopt_;
  std::array<Vector<nub>, N> dummyopt_;
  std::array<Vector<nub>, N> muopt_;

  Vector<dim> solution_, solution_update_; 

  void setInnerSolution() {
    for (size_t i=0; i<N; ++i) {
      const int inucb2 = i * (nuc + 2 * nub);
      solution_.template segment<nuc>(inucb2) = ucopt_[i];
      if constexpr (nub > 0) {
        solution_.template segment<nub>(inucb2+nuc) = dummyopt_[i];
        solution_.template segment<nub>(inucb2+nuc+nub) = muopt_[i];
      }
    }
  }

  void retriveSolution() {
    for (size_t i=0; i<N; ++i) {
      const int inucb2 = i * (nuc + 2 * nub);
      uopt_[i] = solution_.template segment<nu>(inucb2);
      ucopt_[i] = solution_.template segment<nuc>(inucb2);
      if constexpr (nub > 0) {
        dummyopt_[i] = solution_.template segment<nub>(inucb2+nuc);
        muopt_[i] = solution_.template segment<nub>(inucb2+nuc+nub);
      }
    }
  }

};

} // namespace cgmres

#endif // CGMRES__SINGLE_SHOOTING_CGMRES_SOLVER_HPP_