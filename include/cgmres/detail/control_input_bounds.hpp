#ifndef CGMRES__CONTROL_INPUT_BOUNDS_HPP_
#define CGMRES__CONTROL_INPUT_BOUNDS_HPP_

#include <array>

#include "cgmres/types.hpp"

#include "cgmres/detail/macros.hpp"

namespace cgmres {
namespace detail {
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
        = 2.0 * mu.array() * dummy.array() - Map<const Vector<nub>>(ocp.dummy_weight.data()).array();
  }
}

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
void eval_hmu(const OCP& ocp, const MatrixBase<VectorType1>& u, 
                const MatrixBase<VectorType2>& dummy, 
                const MatrixBase<VectorType3>& mu, 
                const MatrixBase<VectorType4>& hmu) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(mu.size() == nub);
    assert(hmu.size() == nub);
    for (int i=0; i<nub; ++i) {
      const auto ui = OCP::ubound_indices[i];
      CGMRES_EIGEN_CONST_CAST(VectorType4, hmu).coeffRef(i)
          = u.coeff(ui) * (u.coeff(ui) - ocp.umin[i] - ocp.umax[i]) 
              + ocp.umin[i] * ocp.umax[i] + dummy.coeff(i) * dummy.coeff(i);
    }
  }
}

template <typename VectorType1, typename VectorType2, typename VectorType3, 
          typename VectorType4, typename VectorType5>
void multiply_hdummy_inv(const MatrixBase<VectorType1>& dummy, 
                         const MatrixBase<VectorType2>& mu, 
                         const MatrixBase<VectorType3>& hdummy, 
                         const MatrixBase<VectorType4>& hmu, 
                         const MatrixBase<VectorType5>& hdummy_multiplied) {
  CGMRES_EIGEN_CONST_CAST(VectorType5, hdummy_multiplied).array() = hmu.array() / (2.0 * dummy.array());
}

template <typename VectorType1, typename VectorType2, typename VectorType3, 
          typename VectorType4, typename VectorType5, typename VectorType6>
void multiply_hmu_inv(const MatrixBase<VectorType1>& dummy, 
                      const MatrixBase<VectorType2>& mu, 
                      const MatrixBase<VectorType3>& hdummy, 
                      const MatrixBase<VectorType4>& hmu, 
                      const MatrixBase<VectorType5>& hdummy_multiplied, 
                      const MatrixBase<VectorType6>& hmu_multiplied) {
  CGMRES_EIGEN_CONST_CAST(VectorType6, hmu_multiplied).array() = hdummy.array() / (2.0 * dummy.array())
                                                                  - mu.array() * hdummy_multiplied.array() / dummy.array();
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
                          / (2.0 * dummy.coeff(i) * dummy.coeff(i));
    }
  }
}

template <typename VectorType>
void clip_dummy(const MatrixBase<VectorType>& dummy, const Scalar min) {
  assert(min >= 0.0);
  const size_t nub = dummy.size();
  for (size_t i=0; i<nub; ++i) {
    CGMRES_EIGEN_CONST_CAST(VectorType, dummy).coeffRef(i) = std::max(dummy.coeff(i), min);
  }
}

} // namespace ubounds
} // namespace detail
} // namespace cgmres

#endif // CGMRES__CONTROL_INPUT_BOUNDS_HPP_