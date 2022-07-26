#ifndef NEWTON_GMRES_HPP_
#define NEWTON_GMRES_HPP_

#include <stdexcept>

#include "cgmres/types.hpp"
#include "cgmres/macros.hpp"


namespace cgmres {

template <class NLP>
class NewtonGMRES {
public:
  static constexpr int nx = NLP::nx;
  static constexpr int nu = NLP::nu;
  static constexpr int nc = NLP::nc;
  static constexpr int dim = NLP::dim;

  NewtonGMRES(const NLP& nlp, const Scalar finite_difference_epsilon) 
    : nlp_(nlp), 
      finite_difference_epsilon_(finite_difference_epsilon),
      updated_solution_(Vector<dim>::Zero()), 
      fonc_(Vector<dim>::Zero()), 
      fonc_1_(Vector<dim>::Zero()) {
    if (finite_difference_epsilon <= 0.0) {
      throw std::invalid_argument("[NewtonGMRES]: 'finite_difference_epsilon' must be positive!");
    }
  }

  ~NewtonGMRES() = default;

  Scalar optError() const {
    return fonc_.template lpNorm<2>();
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void eval_b(const Scalar t, const Vector<nx>& x, 
              const MatrixBase<VectorType1>& solution, 
              const MatrixBase<VectorType2>& solution_update, 
              const MatrixBase<VectorType3>& b_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(b_vec.size() == dim);
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;
    nlp_.eval(t, x, solution, fonc_);
    nlp_.eval(t, x, updated_solution_, fonc_1_);
    EIGEN_CONST_CAST(VectorType3, b_vec) = - fonc_ - (fonc_1_ - fonc_) / finite_difference_epsilon_;
  }

  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void eval_Ax(const Scalar t, const Vector<nx>& x, 
               const MatrixBase<VectorType1>& solution, 
               const MatrixBase<VectorType2>& solution_update, 
               const MatrixBase<VectorType3>& ax_vec) {
    assert(solution.size() == dim);
    assert(solution_update.size() == dim);
    assert(ax_vec.size() == dim);
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;
    nlp_.eval(t, x, updated_solution_, fonc_1_);
    EIGEN_CONST_CAST(VectorType3, ax_vec) = (fonc_1_ - fonc_) / finite_difference_epsilon_;
  }

  const decltype(auto) x() const { return nlp_.x(); }

  const decltype(auto) lmd() const { return nlp_.lmd(); }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  NLP nlp_;
  Scalar finite_difference_epsilon_; 
  Vector<dim> updated_solution_, fonc_, fonc_1_;
};

} // namespace cgmres

#endif // NEWTON_GMRES_HPP_