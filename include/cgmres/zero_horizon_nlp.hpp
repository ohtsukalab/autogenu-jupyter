#ifndef CGMRES__ZERO_HORIZON_NLP_HPP_
#define CGMRES__ZERO_HORIZON_NLP_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"
#include "cgmres/control_input_bounds.hpp"

namespace cgmres {

template <class OCP>
class ZeroHorizonNLP {
public:
  static constexpr int nx = OCP::nx;
  static constexpr int nu = OCP::nu;
  static constexpr int nc = OCP::nc;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = OCP::nub;
  static constexpr int dim = nuc;

  ZeroHorizonNLP(const OCP& ocp) 
    : ocp_(ocp),
      lmd_(Vector<nx>::Zero()) {
    assert(nx > 0);
    assert(nu > 0);
    assert(nc >= 0);
  }

  ZeroHorizonNLP() = default;

  ~ZeroHorizonNLP() = default;

  void eval_fonc_hu(const Scalar t, const Vector<nx>& x, const Vector<dim>& solution,
                    Vector<dim>& fonc_hu) {
    // Compute the Lagrange multiplier over the horizon  
    ocp_.eval_phix(t, x.data(), lmd_.data());
    // Compute the erros in the first order necessary conditions (FONC)
    ocp_.eval_hu(t, x.data(), solution.data(), lmd_.data(), fonc_hu.data());
  }

  void eval_fonc_hu(const Vector<dim>& solution, const Vector<nub>& dummy, 
                    const Vector<nub>& mu, Vector<dim>& fonc_hu) const {
    ubounds::eval_hu(ocp_, solution, dummy, mu, fonc_hu);
  }

  void eval_fonc_hdummy(const Vector<dim>& solution, Vector<nub>& dummy, 
                       Vector<nub>& mu, Vector<nub>& fonc_hdummy) const {
    ubounds::eval_hdummy(ocp_, solution, dummy, mu, fonc_hdummy);
  }

  void eval_fonc_hmu(const Vector<dim>& solution,
                     const Vector<nub>& dummy, 
                     const Vector<nub>& mu,
                     Vector<nub>& fonc_hmu) const {
    ubounds::eval_hmu(ocp_, solution, dummy, mu, fonc_hmu);
  }

  static void multiply_hdummy_inv(const Vector<nub>& dummy, 
                                  const Vector<nub>& mu,
                                  const Vector<nub>& fonc_hdummy,
                                  const Vector<nub>& fonc_hmu,
                                  Vector<nub>& fonc_hdummy_inv) {
    ubounds::multiply_hdummy_inv(dummy, mu, fonc_hdummy, fonc_hmu, fonc_hdummy_inv);
  }

  static void multiply_hmu_inv(const Vector<nub>& dummy, 
                               const Vector<nub>& mu,
                               const Vector<nub>& fonc_hdummy,
                               const Vector<nub>& fonc_hmu,
                               const Vector<nub>& fonc_hdummy_inv,
                               Vector<nub>& fonc_hmu_inv) {
    ubounds::multiply_hmu_inv(dummy, mu, fonc_hdummy, fonc_hmu, fonc_hdummy_inv, fonc_hmu_inv);
  }

  void retrive_dummy_update(const Vector<OCP::nuc>& solution,
                            const Vector<OCP::nub>& dummy, 
                            const Vector<OCP::nub>& mu,
                            const Vector<OCP::nuc>& solution_update,
                            Vector<OCP::nub>& dummy_update) {
    ubounds::retrive_dummy_update(ocp_, solution, dummy, mu, solution_update, dummy_update);
  }

  void retrive_mu_update(const Vector<OCP::nuc>& solution,
                         const Vector<OCP::nub>& dummy, 
                         const Vector<OCP::nub>& mu,
                         const Vector<OCP::nuc>& solution_update,
                         Vector<OCP::nub>& mu_update) {
    ubounds::retrive_mu_update(ocp_, solution, dummy, mu, solution_update, mu_update);
  }

  const OCP& ocp() const { return ocp_; }

  const Vector<nx>& lmd() const { return lmd_; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  OCP ocp_;
  Vector<nx> lmd_;
};

} // namespace cgmres

#endif // CGMRES__ZERO_HORIZON_NLP_HPP_