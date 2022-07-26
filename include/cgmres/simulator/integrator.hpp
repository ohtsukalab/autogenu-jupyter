#include "cgmres/types.hpp"

namespace cgmres {


template <typename OCP, typename StateVectorType, typename ControlInputVectorType>
StateVectorType euler(const OCP& ocp, const Scalar t, const Scalar dt, 
                      const StateVectorType x, const ControlInputVectorType u) {
  StateVectorType x1;
  if (x1.size() != x.size()) {
    x1.resize(x.size());
  }
  ocp.eval_f(t, x.data(), u.data(), x1.data());
  x1.array() *= dt;
  x1.noalias() += x;
  return x1;
}

template <typename OCP, typename StateVectorType, typename ControlInputVectorType>
StateVectorType RK4(const OCP& ocp, const Scalar t, const Scalar dt, 
                    const StateVectorType x, const ControlInputVectorType u) {
  StateVectorType k1, k2, k3, k4, x1;
  if (x1.size() != x.size()) {
    k1.resize(x.size());
    k2.resize(x.size());
    k3.resize(x.size());
    k4.resize(x.size());
    x1.resize(x.size());
  }
  ocp.eval_f(t, x.data(), u.data(), k1.data());
  x1 = x + 0.5 * dt * k1;
  ocp.eval_f(t+0.5*dt, x1.data(), u.data(), k2.data());
  x1 = x + dt * 0.5 * (std::sqrt(2.0)-1.0) * k1 + dt*(1.0-(1.0/std::sqrt(2.0))) * k2;
  ocp.eval_f(t+0.5*dt, x1.data(), u.data(), k3.data());
  x1 = x - dt * 0.5 * std::sqrt(2.0) * k2 + dt * (1.0+(1.0/std::sqrt(2.0))) * k3;
  ocp.eval_f(t+dt, x1.data(), u.data(), k4.data());
  x1 = x + (dt/6.0) * (k1+(2.0-std::sqrt(2.0))*k2 + (2.0+std::sqrt(2.0))*k3+k4);
  return x1;
}

} // namespace cgmres