#ifndef CGMRES__ZERO_HORIZON_OCP_SOLVER_HPP_
#define CGMRES__ZERO_HORIZON_OCP_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/solver_settings.hpp"
#include "cgmres/timer.hpp"

#include "cgmres/detail/matrixfree_gmres.hpp"
#include "cgmres/detail/zero_horizon_nlp.hpp"
#include "cgmres/detail/newton_gmres.hpp"

namespace cgmres {

///
/// @class ZeroHorizonOCPSolver
/// @brief Zero-horizon OCP solver to initialize the solution of the MPC solvers. 
/// @tparam OCP A definition of the optimal control problem (OCP).
/// @tparam kmax Maximum number of the GMRES iterations. Must be positive.
///
template <class OCP, int kmax>
class ZeroHorizonOCPSolver {
public:
  ///
  /// @brief Dimension of the state. 
  ///
  static constexpr int nx = OCP::nx;

  ///
  /// @brief Dimension of the control input. 
  ///
  static constexpr int nu = OCP::nu;

  ///
  /// @brief Dimension of the equality constraints. 
  ///
  static constexpr int nc = OCP::nc;

  ///
  /// @brief Dimension of the concatenation of the control input and equality constraints. 
  ///
  static constexpr int nuc = nu + nc;

  ///
  /// @brief Dimension of the bound constraints on the control input. 
  ///
  static constexpr int nub = OCP::nub;

  ///
  /// @brief Dimension of the linear problem solved by the GMRES solver. 
  ///
  static constexpr int dim = nuc + 2 * nub;

  using ZeroHorizonNLP_ = detail::ZeroHorizonNLP<OCP>;
  using NewtonGMRES_ = detail::NewtonGMRES<ZeroHorizonNLP_>;
  using MatrixFreeGMRES_ = detail::MatrixFreeGMRES<NewtonGMRES_, kmax>;

  ///
  /// @brief Constructs the zero-horizon OCP solver.
  /// @param[in] ocp A definition of the optimal control problem (OCP).
  /// @param[in] settings Solver settings.
  ///
  ZeroHorizonOCPSolver(const OCP& ocp, const SolverSettings& settings) 
    : newton_gmres_(ZeroHorizonNLP_(ocp), settings.finite_difference_epsilon),
      gmres_(),
      settings_(settings),
      uopt_(Vector<nu>::Zero()),
      ucopt_(Vector<nuc>::Zero()),
      solution_(Vector<dim>::Zero()),
      solution_update_(Vector<dim>::Zero()) {
  }

  ///
  /// @brief Default constructor.
  ///
  ZeroHorizonOCPSolver() = default;

  ///
  /// @brief Default destructor.
  ///
  ~ZeroHorizonOCPSolver() = default;

  ///
  /// @brief Sets the control input vector.
  /// @param[in] u The control input vector. Size must be ZeroHorizonOCPSolver::nu.
  ///
  template <typename VectorType>
  void set_u(const VectorType& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_u] u.size() must be " + std::to_string(nu));
    }
    uopt_ = u;
    ucopt_.template head<nu>() = u;
    setInnerSolution();
  }

  ///
  /// @brief Sets the control input vector and Lagrange multiplier with respect to the equality constraints.
  /// @param[in] uc Concatenatin of the control input vector and Lagrange multiplier with respect to the equality constraints. 
  /// Size must be ZeroHorizonOCPSolver::nuc.
  ///
  template <typename VectorType>
  void set_uc(const VectorType& uc) {
    if (uc.size() != nuc) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_uc] uc.size() must be " + std::to_string(nuc));
    }
    uopt_ = uc.template head<nu>();
    ucopt_ = uc;
    setInnerSolution();
  }

  ///
  /// @brief Sets the dummy input vector with respect to the control input bounds constraint.
  /// @param[in] dummy The dummy input vector. Size must be ZeroHorizonOCPSolver::nub.
  ///
  template <typename VectorType>
  void set_dummy(const MatrixBase<VectorType>& dummy) {
    if (dummy.size() != nub) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_dummy] dummy.size() must be " + std::to_string(nub));
    }
    dummyopt_ = dummy;
    setInnerSolution();
  }

  ///
  /// @brief Sets the Lagrange multiplier with respect to the control input bounds constraint.
  /// @param[in] mu The Lagrange multiplier. Size must be ZeroHorizonOCPSolver::nub.
  ///
  template <typename VectorType>
  void set_mu(const MatrixBase<VectorType>& mu) {
    if (mu.size() != nub) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::set_mu] mu.size() must be " + std::to_string(nub));
    }
    muopt_ = mu;
    setInnerSolution();
  }

  ///
  /// @brief Initializes the dummy input vectors and Lagrange multipliers with respect to the control input bounds constraint.
  ///
  void init_dummy_mu() {
    newton_gmres_.retrive_dummy(solution_, settings_.min_dummy);
    newton_gmres_.retrive_mu(solution_);
    retriveSolution();
  }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal control input vector.
  ///
  const Vector<nu>& uopt() const { return uopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal concatenatins of the control input vector and Lagrange multiplier with respect to the equality constraints.
  ///
  const Vector<nuc>& ucopt() const { return ucopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal dummy input vectors with respect to the control input bounds constraint.
  ///
  const Vector<nub>& dummyopt() const { return dummyopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the Lagrange multipliers with respect to the control input bounds constraint.
  ///
  const Vector<nub>& muopt() const { return muopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal costate vector.
  ///
  const Vector<nx>& lmdopt() const { return newton_gmres_.lmd(); }

  ///
  /// @brief Gets the l2-norm of the current optimality errors.
  /// @return The l2-norm of the current optimality errors.
  ///
  Scalar optError() const { return newton_gmres_.optError(); }

  ///
  /// @brief Computes and gets the l2-norm of the current optimality errors.
  /// @param[in] t Time.
  /// @param[in] x State. Size must be ZeroHorizonOCPSolver::nx.
  /// @return The l2-norm of the current optimality errors.
  ///
  template <typename VectorType>
  Scalar optError(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::optError] x.size() must be " + std::to_string(nx));
    }
    newton_gmres_.synchronize_ocp(); 
    newton_gmres_.eval_fonc(t, x, solution_);
    return optError();
  }

  ///
  /// @brief Solves the zero-horizon optimal control problem by Newton-GMRES method.
  /// @param[in] t Time.
  /// @param[in] x State. Size must be ZeroHorizonOCPSolver::nx.
  ///
  template <typename VectorType>
  void solve(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[ZeroHorizonOCPSolver::update] x.size() must be " + std::to_string(nx));
    }
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= solve zero horizon OCP =======================" << std::endl;
    }

    newton_gmres_.synchronize_ocp(); 
    for (size_t iter=0; iter<settings_.max_iter; ++iter) {
      if (settings_.profile_solver) timer_.tick();
      const auto gmres_iter 
          = gmres_.template solve<const Scalar, const VectorType&, const Vector<dim>&>(
                newton_gmres_, t, x.derived(), solution_, solution_update_);
      const auto opt_error = newton_gmres_.optError();
      solution_.noalias() += solution_update_;
      if (settings_.profile_solver) timer_.tock();

      // verbose
      if (settings_.verbose_level >= 1) {
        std::cout << "iter " << iter << ": opt error: " << opt_error 
                  << " (opt tol: " << settings_.opterr_tol << ")" <<  std::endl;
      }
      if (settings_.verbose_level >= 2) {
        std::cout << "         number of GMRES iter: " << gmres_iter 
                  << " (kmax: " << kmax << ")" << std::endl;
      }

      // check convergence
      if (opt_error < settings_.opterr_tol) {
        if (settings_.verbose_level >= 1) {
          std::cout << "converged!" << std::endl;
        }
        break;
      }
    }
    retriveSolution();
  }

  ///
  /// @brief Get timing result as TimingProfile.
  /// @return Timing profile.
  ///
  TimingProfile getProfile() const {
    return timer_.getProfile();
  }

  void disp(std::ostream& os) const {
    os << "Zero horizon OCP solver: " << std::endl;
    os << "  kmax: " << kmax << std::endl;
    os << newton_gmres_.get_nlp().ocp() << std::endl;
    os << settings_ << std::endl;
    os << timer_.getProfile() << std::flush;
  }

  friend std::ostream& operator<<(std::ostream& os, const ZeroHorizonOCPSolver& solver) {
    solver.disp(os);
    return os;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  NewtonGMRES_ newton_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;
  Timer timer_;

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