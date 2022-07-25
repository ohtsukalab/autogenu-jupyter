#ifndef CONTINUATION_GMRES_HPP_
#define CONTINUATION_GMRES_HPP_

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
      fonc1_(Vector<dim>::Zero()),
      fonc2_(Vector<dim>::Zero()),
      x1_(Vector<nx>::Zero()) {
  }

  ~ContinuationGMRES() = default;

  Scalar OptError() const {
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

    const Scalar t1 = t + finite_difference_epsilon_;
    nlp_.ocp().eval_f(t, x.data(), solution.derived().data(), x1_.data());
    x1_.array() *= finite_difference_epsilon_; 
    x1_.noalias() += x;
    updated_solution_ = solution + finite_difference_epsilon_ * solution_update;

    nlp_.eval(t, x, solution, fonc_);
    nlp_.eval(t1, x1_, solution, fonc1_);
    nlp_.eval(t1, x1_, updated_solution_, fonc2_);
    EIGEN_CONST_CAST(VectorType3, b_vec) = - (1.0/finite_difference_epsilon_ - zeta_) * fonc_ 
                                           - fonc2_ / finite_difference_epsilon_;
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
    nlp_.eval(t1, x1_, updated_solution_, fonc2_);
    EIGEN_CONST_CAST(VectorType3, ax_vec) = (fonc2_ - fonc1_) / finite_difference_epsilon_;
  }

private:
  NLP nlp_;
  Scalar finite_difference_epsilon_, zeta_; 
  Vector<dim> updated_solution_, fonc_, fonc1_, fonc2_;
  Vector<nx> x1_;
};

} // namespace cgmres

#endif // CONTINUATION_GMRES_HPP_