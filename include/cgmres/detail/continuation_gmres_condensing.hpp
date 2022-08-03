#ifndef CGMRES__CONTINUATION_GMRES_CONDENSING_HPP_
#define CGMRES__CONTINUATION_GMRES_CONDENSING_HPP_

#include <stdexcept>

#include "cgmres/types.hpp"

#include "cgmres/detail/macros.hpp"


namespace cgmres {
namespace detail {

template <class NLP>
class ContinuationGMRESCondensing {
public:
  static constexpr int nx = NLP::nx;
  static constexpr int nu = NLP::nu;
  static constexpr int nc = NLP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = NLP::nub;
  static constexpr int dim = NLP::dim;
  static constexpr int N = dim / nuc;

  ContinuationGMRESCondensing(const NLP& nlp, 
                              const Scalar finite_difference_epsilon, 
                              const Scalar zeta) 
    : nlp_(nlp), 
      finite_difference_epsilon_(finite_difference_epsilon),
      zeta_(zeta),
      updated_solution_(Vector<dim>::Zero()), 
      fonc_hu_(Vector<dim>::Zero()), 
      fonc_hu_1_(Vector<dim>::Zero()), 
      fonc_hu_2_(Vector<dim>::Zero()), 
      fonc_hu_3_(Vector<dim>::Zero()), 
      x0_1_(Vector<nx>::Zero()),
      dx_(Vector<nx>::Zero()) {
    std::fill(x_1_.begin(), x_1_.end(), Vector<nx>::Zero());
    std::fill(lmd_1_.begin(), lmd_1_.end(), Vector<nx>::Zero());
    std::fill(fonc_f_.begin(), fonc_f_.end(), Vector<nx>::Zero());
    std::fill(fonc_hx_.begin(), fonc_hx_.end(), Vector<nx>::Zero());
    std::fill(fonc_f_1_.begin(), fonc_f_1_.end(), Vector<nx>::Zero());
    std::fill(fonc_hx_1_.begin(), fonc_hx_1_.end(), Vector<nx>::Zero());
    if constexpr (nub > 0) {
      std::fill(dummy_1_.begin(), dummy_1_.end(), Vector<nub>::Zero());
      std::fill(mu_1_.begin(), mu_1_.end(), Vector<nub>::Zero());
      std::fill(fonc_hdummy_.begin(), fonc_hdummy_.end(), Vector<nub>::Zero());
      std::fill(fonc_hmu_.begin(), fonc_hmu_.end(), Vector<nub>::Zero());
      std::fill(fonc_hdummy_1_.begin(), fonc_hdummy_1_.end(), Vector<nub>::Zero());
      std::fill(fonc_hmu_1_.begin(), fonc_hmu_1_.end(), Vector<nub>::Zero());
      std::fill(dummy_update_.begin(), dummy_update_.end(), Vector<nub>::Zero());
      std::fill(mu_update_.begin(), mu_update_.end(), Vector<nub>::Zero());
    }

    if (finite_difference_epsilon <= 0.0) {
      throw std::invalid_argument("[ContinuationGMRESCondensing]: 'finite_difference_epsilon' must be positive!");
    }
    if (zeta <= 0.0) {
      throw std::invalid_argument("[ContinuationGMRESCondensing]: 'zeta' must be positive!");
    }
  }

  ContinuationGMRESCondensing() = default;

  ~ContinuationGMRESCondensing() = default;

  Scalar optError() const {
    Scalar squared_error = fonc_hu_.squaredNorm();
    if constexpr (nub > 0) {
      for (const auto& e : fonc_hdummy_) {
        squared_error += e.squaredNorm();
      }
      for (const auto& e : fonc_hmu_) {
        squared_error += e.squaredNorm();
      }
    }
    for (const auto& e : fonc_f_) {
      squared_error += e.squaredNorm();
    }
    for (const auto& e : fonc_hx_) {
      squared_error += e.squaredNorm();
    }
    return std::sqrt(squared_error);
  }

  template <typename VectorType>
  void eval_fonc(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution,
                 const std::array<Vector<nx>, N+1>& x, const std::array<Vector<nx>, N+1>& lmd,
                 const std::array<Vector<nub>, N>& dummy, const std::array<Vector<nub>, N>& mu) {
    assert(x0.size() == nx);
    nlp_.eval_fonc_hu(t, x0, solution, x, lmd, fonc_hu_);
    if constexpr (nub > 0) {
      nlp_.eval_fonc_hu(solution, dummy, mu, fonc_hu_);
      nlp_.eval_fonc_hdummy(solution, dummy, mu, fonc_hdummy_);
      nlp_.eval_fonc_hmu(solution, dummy, mu, fonc_hmu_);
    }
    nlp_.eval_fonc_f(t, x0, solution, x, fonc_f_);
    nlp_.eval_fonc_hx(t, x0, solution, x, lmd, fonc_hx_);
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
  void eval_b(const Scalar t, const MatrixBase<VectorType1>& x0, 
              const MatrixBase<VectorType2>& solution, 
              const std::array<Vector<nx>, N+1>& x,
              const std::array<Vector<nx>, N+1>& lmd,
              const std::array<Vector<nub>, N>& dummy,
              const std::array<Vector<nub>, N>& mu,
              const MatrixBase<VectorType3>& solution_update, 
              const MatrixBase<VectorType4>& b_vec) {
    assert(x0.size() == nx);
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(b_vec.size() == dim);

    const Scalar t1 = t + finite_difference_epsilon_;
    nlp_.ocp().eval_f(t, x0.derived().data(), solution.derived().data(), dx_.data());
    x0_1_ = x0 + finite_difference_epsilon_ * dx_; 
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;

    nlp_.eval_fonc_hu(t, x0, solution, x, lmd, fonc_hu_);
    nlp_.eval_fonc_hu(t1, x0_1_, solution, x, lmd, fonc_hu_1_);
    if constexpr (nub > 0) {
      nlp_.eval_fonc_hu(solution, dummy, mu, fonc_hu_);
      nlp_.eval_fonc_hu(solution, dummy, mu, fonc_hu_1_);
    }

    // condensing of x and lmd
    nlp_.eval_fonc_f(t, x0, solution, x, fonc_f_);
    nlp_.eval_fonc_hx(t, x0, solution, x, lmd, fonc_hx_);
    for (size_t i=0; i<=N; ++i) {
      fonc_f_1_[i] = (1.0 - finite_difference_epsilon_*zeta_) * fonc_f_[i];
    }
    for (size_t i=0; i<=N; ++i) {
      fonc_hx_1_[i] = (1.0 - finite_difference_epsilon_*zeta_) * fonc_hx_[i];
    }
    nlp_.retrive_x(t1, x0_1_, solution, x_1_, fonc_f_1_);
    nlp_.retrive_lmd(t1, x0_1_, solution, x_1_, lmd_1_, fonc_hx_1_);

    // condensing of dummy and mu
    if constexpr (nub > 0) {
      nlp_.eval_fonc_hdummy(solution, dummy, mu, fonc_hdummy_);
      nlp_.eval_fonc_hmu(solution, dummy, mu, fonc_hmu_);
      for (size_t i=0; i<N; ++i) {
        dummy_1_[i] = - zeta_ * fonc_hdummy_[i];
      }
      for (size_t i=0; i<N; ++i) {
        mu_1_[i] = - zeta_ * fonc_hmu_[i];
      }
      nlp_.multiply_hdummy_inv(dummy, mu, dummy_1_, mu_1_, fonc_hdummy_1_);
      nlp_.multiply_hmu_inv(dummy, mu, dummy_1_, mu_1_, fonc_hdummy_1_, fonc_hmu_1_);
      for (size_t i=0; i<N; ++i) {
        mu_1_[i] = mu[i] + finite_difference_epsilon_ * fonc_hmu_1_[i];
      }
    }

    nlp_.eval_fonc_hu(t1, x0_1_, solution, x_1_, lmd_1_, fonc_hu_3_);
    if constexpr (nub > 0) {
      nlp_.eval_fonc_hu(solution, dummy_1_, mu_1_, fonc_hu_3_);
    }

    nlp_.eval_fonc_f(t1, x0_1_, solution, x, fonc_f_1_);
    nlp_.eval_fonc_hx(t1, x0_1_, solution, x, lmd, fonc_hx_1_);

    nlp_.retrive_x(t1, x0_1_, updated_solution_, x_1_, fonc_f_1_);
    nlp_.retrive_lmd(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hx_1_);
    if constexpr (nub > 0) {
      nlp_.retrive_mu_update(solution, dummy, mu, solution_update, mu_update_);
      for (size_t i=0; i<N; ++i) {
        mu_1_[i] = mu[i] - finite_difference_epsilon_ * mu_update_[i];
      }
    }

    nlp_.eval_fonc_hu(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hu_2_);
    if constexpr (nub > 0) {
      nlp_.eval_fonc_hu(updated_solution_, dummy_1_, mu_1_, fonc_hu_2_);
    }
    CGMRES_EIGEN_CONST_CAST(VectorType4, b_vec) = (1.0/finite_difference_epsilon_ - zeta_) * fonc_hu_ 
                                                  - (fonc_hu_3_ + fonc_hu_2_ - fonc_hu_1_) / finite_difference_epsilon_;
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
  void eval_Ax(const Scalar t, const MatrixBase<VectorType1>& x0, 
               const MatrixBase<VectorType2>& solution, 
               const std::array<Vector<nx>, N+1>& x,
               const std::array<Vector<nx>, N+1>& lmd,
               const std::array<Vector<nub>, N>& dummy,
               const std::array<Vector<nub>, N>& mu,
               const MatrixBase<VectorType3>& solution_update, 
               const MatrixBase<VectorType4>& ax_vec) {
    assert(x0.size() == nx);
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(ax_vec.size() == dim);

    const Scalar t1 = t + finite_difference_epsilon_;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;

    nlp_.retrive_x(t1, x0_1_, updated_solution_, x_1_, fonc_f_1_);
    nlp_.retrive_lmd(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hx_1_);
    if constexpr (nub > 0) {
      nlp_.retrive_mu_update(solution, dummy, mu, solution_update, mu_update_);
      for (size_t i=0; i<N; ++i) {
        mu_1_[i] = mu[i] - finite_difference_epsilon_ * mu_update_[i];
      }
    }

    nlp_.eval_fonc_hu(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hu_2_);
    if constexpr (nub > 0) {
      nlp_.eval_fonc_hu(updated_solution_, dummy_1_, mu_1_, fonc_hu_2_);
    }
    CGMRES_EIGEN_CONST_CAST(VectorType4, ax_vec) = (fonc_hu_2_ - fonc_hu_1_) / finite_difference_epsilon_;
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void expansion(const Scalar t, const MatrixBase<VectorType1>& x0, 
                 const MatrixBase<VectorType2>& solution, 
                 std::array<Vector<nx>, N+1>& x,
                 std::array<Vector<nx>, N+1>& lmd,
                 std::array<Vector<nub>, N>& dummy,
                 std::array<Vector<nub>, N>& mu,
                 const MatrixBase<VectorType3>& solution_update, 
                 const Scalar dt, const Scalar min_dummy) {
    assert(x0.size() == nx);
    const Scalar t1 = t + finite_difference_epsilon_;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;
    for (size_t i=0; i<N+1; ++i) {
      fonc_f_1_[i]  = (1.0 - finite_difference_epsilon_*zeta_) * fonc_f_[i];
    }
    for (size_t i=0; i<N+1; ++i) {
      fonc_hx_1_[i]  = (1.0 - finite_difference_epsilon_*zeta_) * fonc_hx_[i];
    }
    nlp_.retrive_x(t1, x0_1_, updated_solution_, x_1_, fonc_f_1_);
    nlp_.retrive_lmd(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hx_1_);  
    for (size_t i=0; i<N+1; ++i) {
      fonc_f_[i] = x_1_[i] - x[i];
      x[i].noalias() += (dt/finite_difference_epsilon_) * fonc_f_[i];
    }
    for (size_t i=0; i<N+1; ++i) {
      fonc_hx_[i] = lmd_1_[i] - lmd[i];
      lmd[i].noalias() += (dt/finite_difference_epsilon_) * fonc_hx_[i];
    }
    if constexpr (nub > 0) {
      nlp_.retrive_dummy_update(solution, dummy, mu, solution_update, dummy_update_);
      nlp_.retrive_mu_update(solution, dummy, mu, solution_update, mu_update_);
      for (size_t i=0; i<N; ++i) {
        dummy[i].noalias() += dt * (fonc_hdummy_1_[i] - dummy_update_[i]);
      }
      for (size_t i=0; i<N; ++i) {
        mu[i].noalias() += dt * (fonc_hmu_1_[i] - mu_update_[i]);
      }
      nlp_.clip_dummy(dummy, min_dummy);
    }
  }

  template <typename VectorType>
  void retrive_x(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution, 
                 std::array<Vector<nx>, N+1>& x) {
    assert(x0.size() == nx);
    std::fill(fonc_f_1_.begin(), fonc_f_1_.end(), Vector<nx>::Zero());
    nlp_.retrive_x(t, x0, solution, x, fonc_f_1_);
  }

  template <typename VectorType>
  void retrive_lmd(const Scalar t, const MatrixBase<VectorType>& x0, const Vector<dim>& solution, 
                   const std::array<Vector<nx>, N+1>& x,
                   std::array<Vector<nx>, N+1>& lmd) {
    assert(x0.size() == nx);
    std::fill(fonc_hx_1_.begin(), fonc_hx_1_.end(), Vector<nx>::Zero());
    nlp_.retrive_lmd(t, x0, solution, x, lmd, fonc_hx_1_);
  }

  void retrive_dummy(const Vector<dim>& solution, 
                     std::array<Vector<nub>, N>& dummy,
                     const std::array<Vector<nub>, N>& mu,
                     const Scalar min_dummy) {
    if constexpr (nub > 0) {
      std::fill(fonc_hmu_1_.begin(), fonc_hmu_1_.end(), Vector<nub>::Zero());
      std::fill(dummy.begin(), dummy.end(), Vector<nub>::Zero());
      nlp_.eval_fonc_hmu(solution, dummy, mu, fonc_hmu_1_);
      for (size_t i=0; i<N; ++i) {
        dummy[i].array() = fonc_hmu_1_[i].array().abs().sqrt();
      }
      nlp_.clip_dummy(dummy, min_dummy);
    }
  }

  void retrive_mu(const Vector<dim>& solution, 
                  const std::array<Vector<nub>, N>& dummy,
                  std::array<Vector<nub>, N>& mu) {
    if constexpr (nub > 0) {
      std::fill(fonc_hdummy_1_.begin(), fonc_hdummy_1_.end(), Vector<nub>::Zero());
      std::fill(mu.begin(), mu.end(), Vector<nub>::Zero());
      nlp_.eval_fonc_hdummy(solution, dummy, mu, fonc_hdummy_1_);
      for (size_t i=0; i<N; ++i) {
        mu[i].array() = - fonc_hdummy_1_[i].array() / (2.0 * dummy[i].array());
      }
    }
  }

  const NLP& get_nlp() const { return nlp_; }

  void synchronize_ocp() { nlp_.synchronize_ocp(); }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  NLP nlp_;
  Scalar finite_difference_epsilon_, zeta_; 
  Vector<dim> updated_solution_, fonc_hu_, fonc_hu_1_, fonc_hu_2_, fonc_hu_3_;
  std::array<Vector<nx>, N+1> x_1_, lmd_1_, fonc_f_, fonc_hx_, fonc_f_1_, fonc_hx_1_;
  std::array<Vector<nub>, N> dummy_1_, mu_1_, fonc_hdummy_, fonc_hmu_, 
                             fonc_hdummy_1_, fonc_hmu_1_, dummy_update_, mu_update_;
  Vector<nx> x0_1_, dx_;
};

} // namespace detail
} // namespace cgmres

#endif // CGMRES__CONTINUATION_GMRES_CONDENSING_HPP_