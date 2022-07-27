#ifndef CGMRES__CONTROL_INPUT_BOUNDS_HPP_
#define CGMRES__CONTROL_INPUT_BOUNDS_HPP_

#include <array>

#include "cgmres/types.hpp"
#include "cgmres/macros.hpp"

namespace cgmres {
namespace ubounds {

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
void eval_hu(const OCP& ocp, const MatrixBase<VectorType1>& u, 
             const MatrixBase<VectorType2>& dummy, 
             const MatrixBase<VectorType3>& mu, 
             const MatrixBase<VectorType4>& hu) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(mu.size() == nub);
    assert(hu.size() == u.size());
    for (int i=0; i<nub; ++i) {
      const auto ui = OCP::ubound_indices[i];
      CGMRES_EIGEN_CONST_CAST(VectorType4, hu).coeffRef(ui)
          += mu.coeff(i) * (2.0*u.coeff(ui) - ocp.umin[i] - ocp.umax[i]);
    }
  }
}

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
void eval_hdummy(const OCP& ocp, const MatrixBase<VectorType1>& u, 
                   const MatrixBase<VectorType2>& dummy, 
                   const MatrixBase<VectorType3>& mu, 
                   const MatrixBase<VectorType4>& hdummy) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(mu.size() == nub);
    assert(hdummy.size() == nub);
    CGMRES_EIGEN_CONST_CAST(VectorType4, hdummy).array() 
        = 2.0 * mu.array() * dummy.array() - Map<const Vector<nub>>(ocp.dummy_weight).array();
  }
}

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
void eval_hmu(const OCP& ocp, const MatrixBase<VectorType1>& u, 
                const MatrixBase<VectorType2>& dummy, 
                const MatrixBase<VectorType3>& mu, 
                const MatrixBase<VectorType4>& hu_bounds) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(mu.size() == nub);
    assert(hu_bounds.size() == nub);
    for (int i=0; i<nub; ++i) {
      const auto ui = OCP::ubound_indices[i];
      CGMRES_EIGEN_CONST_CAST(VectorType4, hu_bounds).coeffRef(i)
          = u.coeff(ui) * (u.coeff(ui) - ocp.umin[i] - ocp.umax[i]) 
              + ocp.umin[i] * ocp.umax[i] + dummy.coeff(i) * dummy.coeff(i);
    }
  }
}

template <typename VectorType1, typename VectorType2, typename VectorType3, 
          typename VectorType4, typename VectorType5, typename VectorType6>
void eval_hdummy_inv(const MatrixBase<VectorType1>& dummy, 
                     const MatrixBase<VectorType2>& mu, 
                     const MatrixBase<VectorType3>& hdummy, 
                     const MatrixBase<VectorType4>& hmu, 
                     const MatrixBase<VectorType5>& hdummy_inv, 
                     const MatrixBase<VectorType6>& hmu_inv) {
  CGMRES_EIGEN_CONST_CAST(VectorType5, hdummy_inv).array() = hmu.array() / (2.0 * dummy.array());
}

template <typename VectorType1, typename VectorType2, typename VectorType3, 
          typename VectorType4, typename VectorType5, typename VectorType6>
void eval_hmu_inv(const MatrixBase<VectorType1>& dummy, 
                  const MatrixBase<VectorType2>& mu, 
                  const MatrixBase<VectorType3>& hdummy, 
                  const MatrixBase<VectorType4>& hmu, 
                  const MatrixBase<VectorType5>& hdummy_inv, 
                  const MatrixBase<VectorType6>& hmu_inv) {
  CGMRES_EIGEN_CONST_CAST(VectorType6, hmu_inv).array() = hdummy.array() / (2.0 * dummy.array())
                                                          - mu.array() * hdummy_inv.array() / dummy.array();
}

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4, typename VectorType5>
void retrive_dummy_update(const OCP& ocp, 
                          const MatrixBase<VectorType1>& u, 
                          const MatrixBase<VectorType2>& dummy, 
                          const MatrixBase<VectorType3>& mu, 
                          const MatrixBase<VectorType4>& u_update, 
                          const MatrixBase<VectorType5>& dummy_update) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(dummy_update.size() == nub);
    for (int i=0; i<nub; ++i) {
      const auto ui = OCP::ubound_indices[i];
      CGMRES_EIGEN_CONST_CAST(VectorType5, dummy_update).coeffRef(ui) 
          = (2.0*u.coeff(ui) - ocp.umin[i] - ocp.umax[i]) * u_update.coeff(ui) / (2.0 * dummy.coeff(i));
    }
  }
}

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4, typename VectorType5>
void retrive_mu_update(const OCP& ocp, 
                       const MatrixBase<VectorType1>& u, 
                       const MatrixBase<VectorType2>& dummy, 
                       const MatrixBase<VectorType3>& mu, 
                       const MatrixBase<VectorType4>& u_update, 
                       const MatrixBase<VectorType5>& mu_udpate) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(mu_udpate.size() == nub);
    for (int i=0; i<nub; ++i) {
      const auto ui = OCP::ubound_indices[i];
      CGMRES_EIGEN_CONST_CAST(VectorType5, mu_udpate).coeffRef(ui) 
          = - mu.coeff(i) * (2.0*u.coeff(ui) - ocp.umin[i] - ocp.umax[i]) * u_update.coeff(ui) 
                          / (dummy.coeff(i) * dummy.coeff(i));
    }
  }
}

} // namespace ubounds
} // namespace cgmres

#endif // CGMRES__CONTROL_INPUT_BOUNDS_HPP_