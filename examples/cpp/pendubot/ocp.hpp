#ifndef CGMRES__OCP_PENDUBOT_HPP_ 
#define CGMRES__OCP_PENDUBOT_HPP_ 
 
#define _USE_MATH_DEFINES

#include <cmath>
#include <array>
#include <iostream>

#include "cgmres/types.hpp"

namespace cgmres {

/// 
/// @class OCP_pendubot
/// @brief Definition of the optimal control problem (OCP) of pendubot.
/// 
class OCP_pendubot { 
public:
 
  ///
  /// @brief Dimension of the state. 
  ///
  static constexpr int nx = 4;
 
  ///
  /// @brief Dimension of the control input. 
  ///
  static constexpr int nu = 1;
 
  ///
  /// @brief Dimension of the equality constraints. 
  ///
  static constexpr int nc = 0;
 
  ///
  /// @brief Dimension of the Fischer-Burmeister function (already counded in nc). 
  ///
  static constexpr int nh = 0;
 
  ///
  /// @brief Dimension of the concatenation of the control input and equality constraints. 
  ///
  static constexpr int nuc = nu + nc;
 
  ///
  /// @brief Dimension of the bound constraints on the control input. 
  ///
  static constexpr int nub = 1;

  double m1 = 0.2;
  double m2 = 0.7;
  double l1 = 0.3;
  double l2 = 0.3;
  double d1 = 0.15;
  double d2 = 0.257;
  double J1 = 0.006;
  double J2 = 0.051;
  double g = 9.80665;

  std::array<double, 4> q = {1, 1, 0.1, 0.1};
  std::array<double, 4> q_terminal = {1, 1, 0.1, 0.1};
  std::array<double, 4> x_ref = {M_PI, 0, 0, 0};
  std::array<double, 1> r = {0.1};

  static constexpr std::array<int, nub> ubound_indices = {0};
  std::array<double, nub> umin = {-5.0};
  std::array<double, nub> umax = {5.0};
  std::array<double, nub> dummy_weight = {0.1};

  void disp(std::ostream& os) const {
    os << "OCP_pendubot:" << std::endl;
    os << "  nx:  " << nx << std::endl;
    os << "  nu:  " << nu << std::endl;
    os << "  nc:  " << nc << std::endl;
    os << "  nh:  " << nh << std::endl;
    os << "  nuc: " << nuc << std::endl;
    os << "  nub: " << nub << std::endl;
    os << std::endl;
    os << "  m1: " << m1 << std::endl;
    os << "  m2: " << m2 << std::endl;
    os << "  l1: " << l1 << std::endl;
    os << "  l2: " << l2 << std::endl;
    os << "  d1: " << d1 << std::endl;
    os << "  d2: " << d2 << std::endl;
    os << "  J1: " << J1 << std::endl;
    os << "  J2: " << J2 << std::endl;
    os << "  g: " << g << std::endl;
    os << std::endl;
    Eigen::IOFormat fmt(4, 0, ", ", "", "[", "]");
    Eigen::IOFormat intfmt(1, 0, ", ", "", "[", "]");
    os << "  q: " << Map<const VectorX>(q.data(), q.size()).transpose().format(fmt) << std::endl;
    os << "  q_terminal: " << Map<const VectorX>(q_terminal.data(), q_terminal.size()).transpose().format(fmt) << std::endl;
    os << "  x_ref: " << Map<const VectorX>(x_ref.data(), x_ref.size()).transpose().format(fmt) << std::endl;
    os << "  r: " << Map<const VectorX>(r.data(), r.size()).transpose().format(fmt) << std::endl;
    os << std::endl;
    os << "  ubound_indices: " << Map<const VectorXi>(ubound_indices.data(), ubound_indices.size()).transpose().format(intfmt) << std::endl;
    os << "  umin: " << Map<const VectorX>(umin.data(), umin.size()).transpose().format(fmt) << std::endl;
    os << "  umax: " << Map<const VectorX>(umax.data(), umax.size()).transpose().format(fmt) << std::endl;
    os << "  dummy_weight: " << Map<const VectorX>(dummy_weight.data(), dummy_weight.size()).transpose().format(fmt) << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const OCP_pendubot& ocp) { 
    ocp.disp(os);
    return os;
  }


  ///
  /// @brief Computes the state equation dx = f(t, x, u).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[in] u Control input.
  /// @param[out] dx Evaluated value of the state equation.
  ///
  void eval_f(const double t, const double* x, const double* u, 
              double* dx) const {
    const double x0 = pow(l1, 2);
    const double x1 = m2*x0;
    const double x2 = l1*m2;
    const double x3 = d2*x2;
    const double x4 = x3*cos(x[1]);
    const double x5 = pow(d2, 2);
    const double x6 = J2 + m2*x5;
    const double x7 = J1 + pow(d1, 2)*m1;
    const double x8 = 1.0/(x1 + 2.0*x4 + x6 + x7);
    const double x9 = d2*g*m2*sin(x[0] + x[1]);
    const double x10 = d1*m1;
    const double x11 = x10 + x2;
    const double x12 = sin(x[0]);
    const double x13 = sin(x[1]);
    const double x14 = 2.0*x[1];
    const double x15 = 0.5*l1;
    const double x16 = pow(m2, 2)*x5;
    const double x17 = x15*x16;
    const double x18 = pow(x[2], 2.0);
    const double x19 = x[2]*x[3];
    const double x20 = pow(x[3], 2.0);
    dx[0] = x[2];
    dx[1] = x[3];
    dx[2] = x8*(2.0*d2*l1*m2*x13*x[3]*(x[2] + 0.5*x[3]) - g*x11*x12 + u[0] - x9);
    dx[3] = x8*(0.5*d2*g*l1*m2*x11*sin(x[0] - x[1]) + g*x12*(J2*x10 + m2*(J2*l1 + x10*x5) + x17) - g*x17*sin(x14 + x[0]) - u[0]*(x4 + x6) - x0*x16*(x18 + x19 + 0.5*x20)*sin(x14) - x13*x3*(2.0*J2*x19 + J2*x20 + m2*(x0*x18 + x5*pow(x[2] + x[3], 2.0)) + x18*(J2 + x7)) - x9*(0.5*x1 - x10*x15 + x7))/x6;
 
  }

  ///
  /// Computes the partial derivative of terminal cost with respect to state, 
  /// i.e., phix = dphi/dx(t, x).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[out] phix Evaluated value of the partial derivative of terminal cost.
  ///
  void eval_phix(const double t, const double* x, double* phix) const {
    phix[0] = (1.0/2.0)*q_terminal[0]*(2*x[0] - 2*x_ref[0]);
    phix[1] = (1.0/2.0)*q_terminal[1]*(2*x[1] - 2*x_ref[1]);
    phix[2] = (1.0/2.0)*q_terminal[2]*(2*x[2] - 2*x_ref[2]);
    phix[3] = (1.0/2.0)*q_terminal[3]*(2*x[3] - 2*x_ref[3]);
 
  }

  ///
  /// Computes the partial derivative of the Hamiltonian with respect to state, 
  /// i.e., hx = dH/dx(t, x, u, lmd).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[in] u Concatenatin of the control input and Lagrange multiplier with respect to the equality constraints. 
  /// @param[in] lmd Costate. 
  /// @param[out] hx Evaluated value of the partial derivative of the Hamiltonian.
  ///
  void eval_hx(const double t, const double* x, const double* u, 
               const double* lmd, double* hx) const {
    const double x0 = x[0] + x[1];
    const double x1 = d2*g;
    const double x2 = m2*x1;
    const double x3 = x2*cos(x0);
    const double x4 = d1*m1;
    const double x5 = l1*m2;
    const double x6 = x4 + x5;
    const double x7 = cos(x[0]);
    const double x8 = pow(l1, 2);
    const double x9 = m2*x8;
    const double x10 = cos(x[1]);
    const double x11 = d2*x5;
    const double x12 = x10*x11;
    const double x13 = pow(d2, 2);
    const double x14 = m2*x13;
    const double x15 = J2 + x14;
    const double x16 = pow(d1, 2)*m1;
    const double x17 = J1 + x16;
    const double x18 = 1.0/(2.0*x12 + x15 + x17 + x9);
    const double x19 = lmd[2]*x18;
    const double x20 = 0.5*l1;
    const double x21 = pow(m2, 2)*x13;
    const double x22 = x20*x21;
    const double x23 = 2.0*x[1];
    const double x24 = x23 + x[0];
    const double x25 = g*cos(x24);
    const double x26 = x[0] - x[1];
    const double x27 = 0.5*x1*x5*x6*cos(x26);
    const double x28 = 0.5*x9;
    const double x29 = x17 - x20*x4 + x28;
    const double x30 = x29*x3;
    const double x31 = J2*x4 + m2*(J2*l1 + x13*x4) + x22;
    const double x32 = lmd[3]/x15;
    const double x33 = x18*x32;
    const double x34 = x[2] + 0.5*x[3];
    const double x35 = x2*sin(x0);
    const double x36 = sin(x[0]);
    const double x37 = sin(x[1]);
    const double x38 = x11*x37;
    const double x39 = 0.5*x38/pow(0.5*J1 + 0.5*J2 + x12 + 0.5*x14 + 0.5*x16 + x28, 2);
    const double x40 = pow(x[2], 2.0);
    const double x41 = x[2]*x[3];
    const double x42 = pow(x[3], 2.0);
    const double x43 = x40 + x41 + 0.5*x42;
    const double x44 = 2.0*J2;
    const double x45 = J2 + x17;
    const double x46 = x[2] + x[3];
    const double x47 = J2*x42 + m2*(x13*pow(x46, 2.0) + x40*x8) + x40*x45 + x41*x44;
    const double x48 = x21*x8*sin(x23);
    const double x49 = x38*x[3];
    const double x50 = 2.0*pow(x[2], 1.0);
    const double x51 = 2.0*pow(x46, 1.0);
    const double x52 = pow(x[3], 1.0);
    hx[0] = (1.0/2.0)*q[0]*(2*x[0] - 2*x_ref[0]) + x19*(-g*x6*x7 - x3) + x33*(g*x31*x7 - x22*x25 + x27 - x30);
    hx[1] = lmd[2]*x39*(2.0*d2*l1*m2*x34*x37*x[3] - g*x36*x6 + u[0] - x35) + (1.0/2.0)*q[1]*(2*x[1] - 2*x_ref[1]) + x19*(2.0*d2*l1*m2*x10*x34*x[3] - x3) + x32*x39*(0.5*d2*g*l1*m2*x6*sin(x26) - g*x22*sin(x24) + g*x31*x36 - u[0]*(x12 + x15) - x29*x35 - x38*x47 - x43*x48) + x33*(d2*l1*m2*u[0]*x37 - 1.0*l1*x21*x25 - x12*x47 - 2.0*x21*x43*x8*cos(x23) - x27 - x30);
    hx[2] = lmd[0] + (1.0/2.0)*q[2]*(2*x[2] - 2*x_ref[2]) + 2.0*x19*x49 + x33*(-x38*(m2*(x13*x51 + x50*x8) + x44*x[3] + x45*x50) - x48*(x50 + x[3]));
    hx[3] = lmd[1] + (1.0/2.0)*q[3]*(2*x[3] - 2*x_ref[3]) + x19*(2.0*x34*x38 + 1.0*x49) + x33*(-x38*(x14*x51 + x44*x52 + x44*x[2]) - x48*(1.0*x52 + x[2]));
 
  }

  ///
  /// Computes the partial derivative of the Hamiltonian with respect to control input and the equality constraints, 
  /// i.e., hu = dH/du(t, x, u, lmd).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[in] u Concatenatin of the control input and Lagrange multiplier with respect to the equality constraints. 
  /// @param[in] lmd Costate. 
  /// @param[out] hu Evaluated value of the partial derivative of the Hamiltonian.
  ///
  void eval_hu(const double t, const double* x, const double* u, 
               const double* lmd, double* hu) const {
    const double x0 = d2*l1*m2*cos(x[1]);
    const double x1 = J2 + pow(d2, 2)*m2;
    const double x2 = 1.0/(J1 + pow(d1, 2)*m1 + pow(l1, 2)*m2 + 2.0*x0 + x1);
    hu[0] = lmd[2]*x2 + lmd[3]*x2*(-x0 - x1)/x1 + r[0]*u[0];
 
  }
};

} // namespace cgmres

#endif // CGMRES_OCP_HPP_
