#ifndef CGMRES__CONTINUATION_GMRES_HPP_
#define CGMRES__CONTINUATION_GMRES_HPP_

#include <stdexcept>

#include "cgmres/types.hpp"
#include "cgmres/macros.hpp"


namespace cgmres {

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

  void eval_fonc(const Scalar t, const Vector<nx>& x, const Vector<dim>& solution) {
    nlp_.eval(t, x, solution, fonc_);
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void eval_b(const Scalar t, const Vector<nx>& x, 
              const MatrixBase<VectorType1>& solution, 
              const MatrixBase<VectorType2>& solution_update, 
              const MatrixBase<VectorType3>& b_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(b_vec.size() == dim);

    const Scalar t1 = t + finite_difference_epsilon_;
    nlp_.ocp().eval_f(t, x.data(), solution.derived().data(), dx_.data());
    x_1_ = x + finite_difference_epsilon_ * dx_;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;

    nlp_.eval(t, x, solution, fonc_);
    nlp_.eval(t1, x_1_, solution, fonc_1_);
    nlp_.eval(t1, x_1_, updated_solution_, fonc_2_);

    CGMRES_EIGEN_CONST_CAST(VectorType3, b_vec) = (1/finite_difference_epsilon_ - zeta_) * fonc_ 
                                                    - fonc_2_ / finite_difference_epsilon_;
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void eval_Ax(const Scalar t, const Vector<nx>& x, 
               const MatrixBase<VectorType1>& solution, 
               const MatrixBase<VectorType2>& solution_update, 
               const MatrixBase<VectorType3>& ax_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(ax_vec.size() == dim);
    const Scalar t1 = t + finite_difference_epsilon_;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;
    nlp_.eval(t1, x_1_, updated_solution_, fonc_2_);
    CGMRES_EIGEN_CONST_CAST(VectorType3, ax_vec) = (fonc_2_ - fonc_1_) / finite_difference_epsilon_;
  }

  decltype(auto) x() const { return nlp_.x(); }

  decltype(auto) lmd() const { return nlp_.lmd(); }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  NLP nlp_;
  Scalar finite_difference_epsilon_, zeta_; 
  Vector<dim> updated_solution_, fonc_, fonc_1_, fonc_2_;
  Vector<nx> x_1_, dx_;
};

} // namespace cgmres

#endif // CGMRES__CONTINUATION_GMRES_HPP_