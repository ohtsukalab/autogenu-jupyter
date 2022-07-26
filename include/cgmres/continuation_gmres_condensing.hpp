#ifndef CONTINUATION_GMRES_CONDENSING_HPP_
#define CONTINUATION_GMRES_CONDENSING_HPP_

#include <stdexcept>

#include "cgmres/types.hpp"
#include "cgmres/macros.hpp"


namespace cgmres {

template <class NLP>
class ContinuationGMRESCondensing {
public:
  static constexpr int nx = NLP::nx;
  static constexpr int nu = NLP::nu;
  static constexpr int nc = NLP::nc;
  static constexpr int dim = NLP::dim;
  static constexpr int N = dim / (nu + nc);

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
      x0_1_(Vector<nx>::Zero()) {
    std::fill(x_1_.begin(), x_1_.end(), Vector<nx>::Zero());
    std::fill(lmd_1_.begin(), lmd_1_.end(), Vector<nx>::Zero());
    std::fill(fonc_f_.begin(), fonc_f_.end(), Vector<nx>::Zero());
    std::fill(fonc_hx_.begin(), fonc_hx_.end(), Vector<nx>::Zero());
    std::fill(fonc_f_1_.begin(), fonc_f_1_.end(), Vector<nx>::Zero());
    std::fill(fonc_hx_1_.begin(), fonc_hx_1_.end(), Vector<nx>::Zero());
    if (finite_difference_epsilon <= 0.0) {
      throw std::invalid_argument("[ContinuationGMRESCondensing]: 'finite_difference_epsilon' must be positive!");
    }
    if (zeta <= 0.0) {
      throw std::invalid_argument("[ContinuationGMRESCondensing]: 'zeta' must be positive!");
    }
  }

  ~ContinuationGMRESCondensing() = default;

  Scalar optError() const {
    Scalar squared_error = fonc_hu_.squaredNorm();
    for (const auto& e : fonc_f_) {
      squared_error += e.squaredNorm();
    }
    for (const auto& e : fonc_hx_) {
      squared_error += e.squaredNorm();
    }
    return std::sqrt(squared_error);
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void eval_b(const Scalar t, const Vector<nx>& x0, 
              const MatrixBase<VectorType1>& solution, 
              const std::array<Vector<nx>, N+1>& x,
              const std::array<Vector<nx>, N+1>& lmd,
              const MatrixBase<VectorType2>& solution_update, 
              const MatrixBase<VectorType3>& b_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(b_vec.size() == dim);

    const Scalar t1 = t + finite_difference_epsilon_;
    nlp_.ocp().eval_f(t, x0.data(), solution.derived().data(), x0_1_.data());
    x0_1_.array() *= finite_difference_epsilon_; 
    x0_1_.noalias() += x0;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;

    nlp_.eval_fonc_hu(t, x0, solution, x, lmd, fonc_hu_);
    nlp_.eval_fonc_hu(t1, x0_1_, solution, x, lmd, fonc_hu_1_);

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

    nlp_.eval_fonc_hu(t1, x0_1_, solution, x_1_, lmd_1_, fonc_hu_3_);
    nlp_.eval_fonc_f(t1, x0_1_, solution, x, fonc_f_1_);
    nlp_.eval_fonc_hx(t1, x0_1_, solution, x, lmd, fonc_hx_1_);

    nlp_.retrive_x(t1, x0_1_, updated_solution_, x_1_, fonc_f_1_);
    nlp_.retrive_lmd(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hx_1_);
    nlp_.eval_fonc_hu(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hu_2_);
    EIGEN_CONST_CAST(VectorType3, b_vec) = (1.0/finite_difference_epsilon_ - zeta_) * fonc_hu_ 
                                           - fonc_hu_3_ / finite_difference_epsilon_
                                           - (fonc_hu_2_ - fonc_hu_1_) / finite_difference_epsilon_;
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void eval_Ax(const Scalar t, const Vector<nx>& x0, 
               const MatrixBase<VectorType1>& solution, 
               const std::array<Vector<nx>, N+1>& x,
               const std::array<Vector<nx>, N+1>& lmd,
               const MatrixBase<VectorType2>& solution_update, 
               const MatrixBase<VectorType3>& ax_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(ax_vec.size() == dim);

    const Scalar t1 = t + finite_difference_epsilon_;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;

    nlp_.retrive_x(t1, x0_1_, updated_solution_, x_1_, fonc_f_1_);
    nlp_.retrive_lmd(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hx_1_);
    nlp_.eval_fonc_hu(t1, x0_1_, updated_solution_, x_1_, lmd_1_, fonc_hu_2_);
    EIGEN_CONST_CAST(VectorType3, ax_vec) = (fonc_hu_2_ - fonc_hu_1_) / finite_difference_epsilon_;
  }

  template <typename VectorType1, typename VectorType2>
  void expansion(const Scalar t, const Vector<nx>& x0, 
                 const MatrixBase<VectorType1>& solution, 
                 std::array<Vector<nx>, N+1>& x,
                 std::array<Vector<nx>, N+1>& lmd,
                 const MatrixBase<VectorType2>& solution_update, 
                 const Scalar dt) {
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
  }

  const decltype(auto) x() const { return nlp_.x(); }

  const decltype(auto) lmd() const { return nlp_.lmd(); }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  NLP nlp_;
  Scalar finite_difference_epsilon_, zeta_; 
  Vector<dim> updated_solution_, fonc_hu_, fonc_hu_1_, fonc_hu_2_, fonc_hu_3_;
  std::array<Vector<nx>, N+1> x_1_, lmd_1_, fonc_f_, fonc_hx_, fonc_f_1_, fonc_hx_1_;
  Vector<nx> x0_1_;
};

} // namespace cgmres

#endif // CONTINUATION_GMRES_CONDENSING_HPP_