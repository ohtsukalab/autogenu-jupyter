#ifndef CGMRES__MULTIPLE_SHOOTING_CGMRES_SOLVER_HPP_
#define CGMRES__MULTIPLE_SHOOTING_CGMRES_SOLVER_HPP_

#include <array>
#include <stdexcept>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/solver_settings.hpp"
#include "cgmres/timer.hpp"

#include "cgmres/detail/matrixfree_gmres.hpp"
#include "cgmres/detail/multiple_shooting_nlp.hpp"
#include "cgmres/detail/continuation_gmres_condensing.hpp"

namespace cgmres {

///
/// @class MultipleShootingCGMRESSolver
/// @brief Multiple-shooting C/GMRES solver for nonlinear MPC. 
/// @tparam OCP A definition of the optimal control problem (OCP).
/// @tparam N Number of discretizationn grids of the horizon. Must be positive.
/// @tparam kmax Maximum number of the GMRES iterations. Must be positive.
///
template <class OCP, int N, int kmax>
class MultipleShootingCGMRESSolver {
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
  static constexpr int dim = nuc * N;

  using MultipleShootingNLP_ = detail::MultipleShootingNLP<OCP, N>;
  using ContinuationGMRES_ = detail::ContinuationGMRESCondensing<MultipleShootingNLP_>;
  using MatrixFreeGMRES_ = detail::MatrixFreeGMRES<ContinuationGMRES_, kmax>;

  ///
  /// @brief Constructs the multiple-shooting C/GMRES solver.
  /// @param[in] ocp A definition of the optimal control problem (OCP).
  /// @param[in] horizon Prediction horizon of MPC.
  /// @param[in] settings Solver settings.
  ///
  MultipleShootingCGMRESSolver(const OCP& ocp, const Horizon& horizon, 
                               const SolverSettings& settings) 
    : continuation_gmres_(MultipleShootingNLP_(ocp, horizon), settings.finite_difference_epsilon, settings.zeta),
      gmres_(),
      settings_(settings),
      solution_(Vector<dim>::Zero()),
      solution_update_(Vector<dim>::Zero()) {
    std::fill(uopt_.begin(), uopt_.end(), Vector<nu>::Zero());
    std::fill(ucopt_.begin(), ucopt_.end(), Vector<nuc>::Zero());
    std::fill(xopt_.begin(), xopt_.end(), Vector<nx>::Zero());
    std::fill(lmdopt_.begin(), lmdopt_.end(), Vector<nx>::Zero());
    if constexpr (nub > 0) {
      std::fill(dummyopt_.begin(), dummyopt_.end(), Vector<nub>::Zero());
      std::fill(muopt_.begin(), muopt_.end(), Vector<nub>::Zero());
    }

    if (settings.finite_difference_epsilon <= 0.0) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'settings.finite_difference_epsilon' must be positive!");
    }
    if (settings.sampling_time <= 0.0) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'settings.sampling_time' must be positive!");
    }
    if (settings.zeta <= 0.0) {
      throw std::invalid_argument("[ContinuationGMRESCondensing]: 'settings.zeta' must be positive!");
    }
  }

  ///
  /// @brief Default constructor.
  ///
  MultipleShootingCGMRESSolver() = default;

  ///
  /// @brief Default destructor.
  ///
  ~MultipleShootingCGMRESSolver() = default;

  ///
  /// @brief Sets the control input vector.
  /// @param[in] u The control input vector. Size must be MultipleShootingCGMRESSolver::nu.
  ///
  template <typename VectorType>
  void set_u(const MatrixBase<VectorType>& u) {
    if (u.size() != nu) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_u] u.size() must be " + std::to_string(nu));
    }
    for (size_t i=0; i<N; ++i) {
      uopt_[i] = u;
    }
    for (size_t i=0; i<N; ++i) {
      ucopt_[i].template head<nu>() = u;
    }
    setInnerSolution();
  }

  ///
  /// @brief Sets the control input vector and Lagrange multiplier with respect to the equality constraints.
  /// @param[in] uc Concatenatin of the control input vector and Lagrange multiplier with respect to the equality constraints. 
  /// Size must be MultipleShootingCGMRESSolver::nuc.
  ///
  template <typename VectorType>
  void set_uc(const MatrixBase<VectorType>& uc) {
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

  ///
  /// @brief Sets the state vector.
  /// @param[in] x The state vector. Size must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename VectorType>
  void set_x(const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_x] x.size() must be " + std::to_string(nx));
    }
    for (size_t i=1; i<=N; ++i) {
      xopt_[i] = x;
    }
  }

  ///
  /// @brief Sets the costate vector.
  /// @param[in] lmd The costate vector. Size must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename VectorType>
  void set_lmd(const MatrixBase<VectorType>& lmd) {
    if (lmd.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_lmd] lmd.size() must be " + std::to_string(nx));
    }
    for (size_t i=1; i<=N; ++i) {
      lmdopt_[i] = lmd;
    }
  }

  ///
  /// @brief Sets the dummy input vector with respect to the control input bounds constraint.
  /// @param[in] dummy The dummy input vector. Size must be MultipleShootingCGMRESSolver::nub.
  ///
  template <typename VectorType>
  void set_dummy(const MatrixBase<VectorType>& dummy) {
    if (dummy.size() != nub) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_dummy] dummy.size() must be " + std::to_string(nub));
    }
    for (size_t i=0; i<N; ++i) {
      dummyopt_[i] = dummy;
    }
  }

  ///
  /// @brief Sets the Lagrange multiplier with respect to the control input bounds constraint.
  /// @param[in] mu The Lagrange multiplier. Size must be MultipleShootingCGMRESSolver::nub.
  ///
  template <typename VectorType>
  void set_mu(const MatrixBase<VectorType>& mu) {
    if (mu.size() != nub) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_mu] mu.size() must be " + std::to_string(nub));
    }
    for (size_t i=0; i<N; ++i) {
      muopt_[i] = mu;
    }
  }

  ///
  /// @brief Sets the array of control input vector.
  /// @param[in] u_array A container (std::vector, std::array, etc.) of the control input vector. Size must be N and size of each element must be MultipleShootingCGMRESSolver::nu.
  ///
  template <typename T>
  void set_u_array(const T& u_array) {
    if (u_array.size() != N) { 
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_u_array]: 'u_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (u_array[i].size() != nu) {
        throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_u_array] u_array[i].size() must be " + std::to_string(nu));
      }
      uopt_[i] = u_array[i];
      ucopt_[i].template head<nu>() = u_array[i];
    }
    setInnerSolution();
  }

  ///
  /// @brief Sets the array of control input vector and Lagrange multiplier with respect to the equality constraints.
  /// @param[in] uc_array A container (std::vector, std::array, etc.) of the concatenatin of the control input vector and Lagrange multiplier 
  /// with respect to the equality constraints. Size must be N and size of each element must be MultipleShootingCGMRESSolver::nuc.
  ///
  template <typename T>
  void set_uc_array(const T& uc_array) {
    if (uc_array.size() != N) { 
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_uc_array]: 'uc_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (uc_array[i].size() != nuc) {
        throw std::invalid_argument("[MultipleShootingCGMRESSolver::uc_array] uc_array[i].size() must be " + std::to_string(nuc));
      }
      uopt_[i] = uc_array[i].template head<nu>();
      ucopt_[i] = uc_array[i];
    }
    setInnerSolution();
  }

  ///
  /// @brief Sets the array of state vector.
  /// @param[in] x_array A container (std::vector, std::array, etc.) of the state vector. Size must be N+1 and size of each element must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename T>
  void set_x_array(const T& x_array) {
    if (x_array.size() != N+1) { 
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_x_array]: 'x_array.size()' must be "+std::to_string(N+1)); 
    } 
    for (size_t i=1; i<=N; ++i) {
      if (x_array.size() != nx) {
        throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_x_array] x_array[i].size() must be " + std::to_string(nx));
      }
      xopt_[i] = x_array[i];
    }
  }

  ///
  /// @brief Sets the array of costate vector.
  /// @param[in] lmd_array A container (std::vector, std::array, etc.) of the costate vector. Size must be N+1 and size of each element must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename T>
  void set_lmd_array(const T& lmd_array) {
    if (lmd_array.size() != N+1) { 
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_lmd_array]: 'lmd_array.size()' must be "+std::to_string(N+1)); 
    } 
    for (size_t i=1; i<=N; ++i) {
      if (lmd_array.size() != nx) {
        throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_lmd_array] lmd_array[i].size() must be " + std::to_string(nx));
      }
      lmdopt_[i] = lmd_array[i];
    }
  }

  ///
  /// @brief Sets the array of the dummy input vector with respect to the control input bounds constraint.
  /// @param[in] dummy_array A container (std::vector, std::array, etc.) of the dummy input vector. 
  /// Size must be N and size of each element must be MultipleShootingCGMRESSolver::nub.
  ///
  template <typename T>
  void set_dummy_array(const T& dummy_array) {
    if (dummy_array.size() != N) { 
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_dummy_array]: 'dummy_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (dummy_array[i].size() != nub) {
        throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_dummy_array] dummy_array[i].size() must be " + std::to_string(nub));
      }
      dummyopt_[i] = dummy_array[i];
    }
  }

  ///
  /// @brief Sets the array of the Lagrange multiplier with respect to the control input bounds constraint.
  /// @param[in] mu_array A container (std::vector, std::array, etc.) of the Lagrange multiplier. 
  /// Size must be N and size of each element must be MultipleShootingCGMRESSolver::nub.
  ///
  template <typename T>
  void set_mu_array(const T& mu_array) {
    if (mu_array.size() != N) { 
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_mu_array]: 'mu_array.size()' must be "+std::to_string(N)); 
    } 
    for (size_t i=0; i<N; ++i) {
      if (mu_array[i].size() != nub) {
        throw std::invalid_argument("[MultipleShootingCGMRESSolver::set_mu_array] mu_array[i].size() must be " + std::to_string(nub));
      }
      muopt_[i] = mu_array[i];
    }
  }

  ///
  /// @brief Initializes the state vectors by simulating the system dynamics over the horizon.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] x Initial state of the horizon. Size must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename VectorType>
  void init_x(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::init_x] x.size() must be " + std::to_string(nx));
    }
    continuation_gmres_.retrive_x(t, x, solution_, xopt_);
  }

  ///
  /// @brief Initializes the costate vectors by simulating the system costate dynamics over the horizon.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] x Initial state of the horizon. Size must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename VectorType>
  void init_lmd(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::init_lmd] x.size() must be " + std::to_string(nx));
    }
    continuation_gmres_.retrive_lmd(t, x, solution_, xopt_, lmdopt_);
  }

  ///
  /// @brief Initializes the state vectors and costate vectors by simulating the system dynamics and system costate dynamcis over the horizon.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] x Initial state of the horizon. Size must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename VectorType>
  void init_x_lmd(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::init_x_lmd] x.size() must be " + std::to_string(nx));
    }
    continuation_gmres_.retrive_x(t, x, solution_, xopt_);
    continuation_gmres_.retrive_lmd(t, x, solution_, xopt_, lmdopt_);
  }

  ///
  /// @brief Initializes the dummy input vectors and Lagrange multipliers with respect to the control input bounds constraint.
  ///
  void init_dummy_mu() {
    continuation_gmres_.retrive_dummy(solution_, dummyopt_, muopt_, settings_.min_dummy);
    continuation_gmres_.retrive_mu(solution_, dummyopt_, muopt_);
  }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal control input vectors over the horizon.
  ///
  const std::array<Vector<nu>, N>& uopt() const { return uopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal concatenatins of the control input vector and Lagrange multiplier with respect to the equality constraints.
  ///
  const std::array<Vector<nuc>, N>& ucopt() const { return ucopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal state vectors over the horizon.
  ///
  const std::array<Vector<nx>, N+1>& xopt() const { return xopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal costate vectors over the horizon.
  ///
  const std::array<Vector<nx>, N+1>& lmdopt() const { return lmdopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the optimal dummy input vectors with respect to the control input bounds constraint over the horizon.
  ///
  const std::array<Vector<nub>, N>& dummyopt() const { return dummyopt_; }

  ///
  /// @brief Getter of the optimal solution.
  /// @return const reference to the Lagrange multipliers with respect to the control input bounds constraint over the horizon.
  ///
  const std::array<Vector<nub>, N>& muopt() const { return muopt_; }

  ///
  /// @brief Gets the l2-norm of the current optimality errors.
  /// @return The l2-norm of the current optimality errors.
  ///
  Scalar optError() const { return continuation_gmres_.optError(); }

  ///
  /// @brief Computes and gets the l2-norm of the current optimality errors.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] x Initial state of the horizon. Size must be MultipleShootingCGMRESSolver::nx.
  /// @return The l2-norm of the current optimality errors.
  ///
  template <typename VectorType>
  Scalar optError(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::optError] x.size() must be " + std::to_string(nx));
    }
    continuation_gmres_.synchronize_ocp(); 
    continuation_gmres_.eval_fonc(t, x, solution_, xopt_, lmdopt_, dummyopt_, muopt_);
    return optError();
  }

  ///
  /// @brief Updates the solution by performing C/GMRES method.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] x Initial state of the horizon. Size must be MultipleShootingCGMRESSolver::nx.
  ///
  template <typename VectorType>
  void update(const Scalar t, const MatrixBase<VectorType>& x) {
    if (x.size() != nx) {
      throw std::invalid_argument("[MultipleShootingCGMRESSolver::update] x.size() must be " + std::to_string(nx));
    }
    if (settings_.verbose_level >= 1) {
      std::cout << "\n======================= update solution with C/GMRES =======================" << std::endl;
    }

    if (settings_.profile_solver) timer_.tick();
    continuation_gmres_.synchronize_ocp(); 
    const auto gmres_iter 
        = gmres_.template solve<const Scalar, const MatrixBase<VectorType>&, const Vector<dim>&,
                                const std::array<Vector<nx>, N+1>&, const std::array<Vector<nx>, N+1>&,
                                const std::array<Vector<nub>, N>&, const std::array<Vector<nub>, N>&>(
              continuation_gmres_, t, x.derived(), solution_, xopt_, lmdopt_, dummyopt_, muopt_, solution_update_);
    const auto opt_error = continuation_gmres_.optError();
    continuation_gmres_.expansion(t, x, solution_, xopt_, lmdopt_, dummyopt_, muopt_, 
                                  solution_update_, settings_.sampling_time, settings_.min_dummy);
    solution_.noalias() += settings_.sampling_time * solution_update_;
    retriveSolution();
    if (settings_.profile_solver) timer_.tock();

    // verbose
    if (settings_.verbose_level >= 1) {
      std::cout << "opt error: " << opt_error << std::endl;
    }
    if (settings_.verbose_level >= 2) {
      std::cout << "number of GMRES iter: " << gmres_iter << " (kmax: " << kmax << ")" << std::endl;
    }
  }

  ///
  /// @brief Get timing result as TimingProfile.
  /// @return Timing profile.
  ///
  TimingProfile getProfile() const {
    return timer_.getProfile();
  }

  void disp(std::ostream& os) const {
    os << "Multiple shooting CGMRES solver: " << std::endl;
    os << "  N:    " << N << std::endl;
    os << "  kmax: " << kmax << "\n" << std::endl;
    os << continuation_gmres_.get_nlp().ocp() << std::endl;
    os << continuation_gmres_.get_nlp().horizon() << std::endl;
    os << settings_ << std::endl;
    os << timer_.getProfile() << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const MultipleShootingCGMRESSolver& solver) {
    solver.disp(os);
    return os;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ContinuationGMRES_ continuation_gmres_;
  MatrixFreeGMRES_ gmres_;
  SolverSettings settings_;
  Timer timer_;

  std::array<Vector<nu>, N> uopt_;
  std::array<Vector<nuc>, N> ucopt_;
  std::array<Vector<nx>, N+1> xopt_;
  std::array<Vector<nx>, N+1> lmdopt_;
  std::array<Vector<nub>, N> dummyopt_;
  std::array<Vector<nub>, N> muopt_;

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