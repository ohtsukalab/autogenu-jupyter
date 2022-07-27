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
      constexpr int ui = OCP::ubound_indices[i];
      CGMRES_EIGEN_CONST_CAST(VectorType4, hu).coeffRef(ui)
          += mu.coeff(i) * (2.0*u.coeff(ui) - ocp.umin[i] - ocp.umax[i]);
    }
  }
}

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
void eval_hu_dummy(const OCP& ocp, const MatrixBase<VectorType1>& u, 
                   const MatrixBase<VectorType2>& dummy, 
                   const MatrixBase<VectorType3>& mu, 
                   const MatrixBase<VectorType4>& hu_dummy) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(mu.size() == nub);
    assert(hu_dummy.size() == nub);
    CGMRES_EIGEN_CONST_CAST(VectorType4, hu_dummy).array() 
        = 2.0 * mu.array() * dummy.array() - Map<Vector<nub>>(ocp.dummy_weight).array();
  }
}

template <typename OCP, typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
void eval_hu_bounds(const OCP& ocp, const MatrixBase<VectorType1>& u, 
                    const MatrixBase<VectorType2>& dummy, 
                    const MatrixBase<VectorType3>& mu, 
                    const MatrixBase<VectorType4>& hu_bounds) {
  constexpr int nub = OCP::nub;
  if constexpr (nub > 0) {
    assert(dummy.size() == nub);
    assert(mu.size() == nub);
    assert(hu_bounds.size() == nub);
    for (int i=0; i<nub; ++i) {
      constexpr int ui = OCP::ubound_indices[i];
      CGMRES_EIGEN_CONST_CAST(VectorType4, hu_bounds).coeffRef(i)
          = u.coeff(ui) * (u.coeff(ui) - ocp.umin[i] - ocp.umax[i]) 
              + ocp.umin[i] * ocp.umax[i] + dummy.coeff(i) * dummy.coeff(i);
    }
  }
}

} // namespace ubounds
} // namespace cgmres

#endif // CGMRES__CONTROL_INPUT_BOUNDS_HPP_