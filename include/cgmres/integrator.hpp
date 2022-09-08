#ifndef CGMRES__INTEGRATOR_HPP_
#define CGMRES__INTEGRATOR_HPP_

#include "cgmres/types.hpp"

namespace cgmres {

///
/// @brief Computes the next state by the forward Euler.
/// @param[in] ocp Optimal control problem.
/// @param[in] t Time.
/// @param[in] dt Time step.
/// @param[in] x State.
/// @param[in] u Control input.
/// @return State at time t+dt.
///
template <typename OCP, typename StateVectorType, typename ControlInputVectorType>
VectorX ForwardEuler(const OCP& ocp, const Scalar t, const Scalar dt, 
                     const MatrixBase<StateVectorType>& x, 
                     const MatrixBase<ControlInputVectorType >& u) {
  VectorX x1;
  x1.setZero(x.size());
  ocp.eval_f(t, x, u, x1);
  x1.array() *= dt;
  x1.noalias() += x;
  return x1;
}


///
/// @brief Computes the next state by the 4th-order Runge-Kutta method.
/// @param[in] ocp Optimal control problem.
/// @param[in] t Time.
/// @param[in] dt Time step.
/// @param[in] x State.
/// @param[in] u Control input.
/// @return State at time t+dt.
///
template <typename OCP, typename StateVectorType, typename ControlInputVectorType>
VectorX RK4(const OCP& ocp, const Scalar t, const Scalar dt, 
            const MatrixBase<StateVectorType>& x, 
            const MatrixBase<ControlInputVectorType >& u) {
  VectorX k1, k2, k3, k4, x1;
  k1.setZero(x.size());
  k2.setZero(x.size());
  k3.setZero(x.size());
  k4.setZero(x.size());
  x1.setZero(x.size());
  ocp.eval_f(t, x, u, k1);
  x1 = x + 0.5 * dt * k1;
  ocp.eval_f(t+0.5*dt, x1, u, k2);
  x1 = x + dt * 0.5 * (std::sqrt(2.0)-1.0) * k1 + dt*(1.0-(1.0/std::sqrt(2.0))) * k2;
  ocp.eval_f(t+0.5*dt, x1, u, k3);
  x1 = x - dt * 0.5 * std::sqrt(2.0) * k2 + dt * (1.0+(1.0/std::sqrt(2.0))) * k3;
  ocp.eval_f(t+dt, x1, u, k4);
  x1 = x + (dt/6.0) * (k1+(2.0-std::sqrt(2.0))*k2 + (2.0+std::sqrt(2.0))*k3+k4);
  return x1;
}

} // namespace cgmres 

#endif // CGMRES__INTEGRATOR_HPP_