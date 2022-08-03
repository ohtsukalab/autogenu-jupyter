#ifndef CGMRES__CONTINUATION_GMRES_HPP_
#define CGMRES__CONTINUATION_GMRES_HPP_

#include <stdexcept>

#include "cgmres/types.hpp"

#include "cgmres/detail/macros.hpp"


namespace cgmres {
namespace detail {

template <class NLP>
class ContinuationGMRES {
public:
  static constexpr int nx = NLP::nx;
  static constexpr int nu = NLP::nu;
  static constexpr int nc = NLP::nc;
  static constexpr int dim = NLP::dim;

  ContinuationGMRES(const NLP& nlp, const Scalar finite_difference_epsilon, 
                    const Scalar zeta) 
    : nlp_(nlp), 
      finite_difference_epsilon_(finite_difference_epsilon),
      zeta_(zeta),
      updated_solution_(Vector<dim>::Zero()), 
      fonc_(Vector<dim>::Zero()), 
      fonc_1_(Vector<dim>::Zero()),
      fonc_2_(Vector<dim>::Zero()),
      x_1_(Vector<nx>::Zero()),
      dx_(Vector<nx>::Zero()) {
    if (finite_difference_epsilon <= 0.0) {
      throw std::invalid_argument("[ContinuationGMRES]: 'finite_difference_epsilon' must be positive!");
    }
    if (zeta <= 0.0) {
      throw std::invalid_argument("[ContinuationGMRES]: 'zeta' must be positive!");
    }
  }

  ContinuationGMRES() = default;

  ~ContinuationGMRES() = default;

  Scalar optError() const {
    return fonc_.template lpNorm<2>();
  }

  template <typename VectorType>
  void eval_fonc(const Scalar t, const MatrixBase<VectorType>& x, const Vector<dim>& solution) {
    nlp_.eval_fonc_hu(t, x, solution, fonc_);
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
  void eval_b(const Scalar t, const MatrixBase<VectorType1>& x, 
              const MatrixBase<VectorType2>& solution, 
              const MatrixBase<VectorType3>& solution_update, 
              const MatrixBase<VectorType4>& b_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(b_vec.size() == dim);

    const Scalar t1 = t + finite_difference_epsilon_;
    nlp_.ocp().eval_f(t, x.derived().data(), solution.derived().data(), dx_.data());
    x_1_ = x + finite_difference_epsilon_ * dx_;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;

    nlp_.eval_fonc_hu(t, x, solution, fonc_);
    nlp_.eval_fonc_hu(t1, x_1_, solution, fonc_1_);
    nlp_.eval_fonc_hu(t1, x_1_, updated_solution_, fonc_2_);

    CGMRES_EIGEN_CONST_CAST(VectorType4, b_vec) = (1/finite_difference_epsilon_ - zeta_) * fonc_ 
                                                    - fonc_2_ / finite_difference_epsilon_;
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
  void eval_Ax(const Scalar t, const MatrixBase<VectorType1>& x,
               const MatrixBase<VectorType2>& solution, 
               const MatrixBase<VectorType3>& solution_update, 
               const MatrixBase<VectorType4>& ax_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(ax_vec.size() == dim);
    const Scalar t1 = t + finite_difference_epsilon_;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;
    nlp_.eval_fonc_hu(t1, x_1_, updated_solution_, fonc_2_);
    CGMRES_EIGEN_CONST_CAST(VectorType4, ax_vec) = (fonc_2_ - fonc_1_) / finite_difference_epsilon_;
  }

  void retrive_dummy(Vector<dim>& solution, const Scalar min_dummy) {
    fonc_1_.setZero();
    nlp_.retrive_dummy(solution, fonc_1_, min_dummy);
  }

  void retrive_mu(Vector<dim>& solution) {
    fonc_1_.setZero();
    nlp_.retrive_mu(solution, fonc_1_);
  }

  decltype(auto) x() const { return nlp_.x(); }

  decltype(auto) lmd() const { return nlp_.lmd(); }

  const NLP& get_nlp() const { return nlp_; }

  void synchronize_ocp() { nlp_.synchronize_ocp(); }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  NLP nlp_;
  Scalar finite_difference_epsilon_, zeta_; 
  Vector<dim> updated_solution_, fonc_, fonc_1_, fonc_2_;
  Vector<nx> x_1_, dx_;
};

} // namespace detail
} // namespace cgmres

#endif // CGMRES__CONTINUATION_GMRES_HPP_