#ifndef CGMRES__OCP_CARTPOLE_HPP_ 
#define CGMRES__OCP_CARTPOLE_HPP_ 
 
#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include "cgmres/types.hpp"

namespace cgmres {

// This class defines the optimal control problem (OCP)
class OCP_cartpole { 
public:
  static constexpr int nx = 4;
  static constexpr int nu = 1;
  static constexpr int nc = 0;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = 1;

  double m_c = 2;
  double m_p = 0.2;
  double l = 0.5;
  double g = 9.80665;

  double q[4] = {2.5, 10, 0.01, 0.01};
  double q_terminal[4] = {2.5, 10, 0.01, 0.01};
  double x_ref[4] = {0, M_PI, 0, 0};
  double r[1] = {1};

  static constexpr int ubound_indices[1] = {0};
  double umin[1] = {-15.0};
  double umax[1] = {15.0};
  double dummy_weight[1] = {0.1};

  void disp(std::ostream& os) const {
    os << "OCP_cartpole:" << std::endl;
    os << "  nx:  " << nx << std::endl;
    os << "  nu:  " << nu << std::endl;
    os << "  nc:  " << nc << std::endl;
    os << "  nuc: " << nuc << std::endl;
    os << "  nub: " << nub << std::endl;
    os << std::endl;
    os << "  m_c: " << m_c << std::endl;
    os << "  m_p: " << m_p << std::endl;
    os << "  l: " << l << std::endl;
    os << "  g: " << g << std::endl;
    os << std::endl;
    Eigen::IOFormat fmt(4, 0, ", ", "", "[", "]");
    Eigen::IOFormat intfmt(1, 0, ", ", "", "[", "]");
    os << "  q: " << Map<const VectorX>(q, 4).transpose().format(fmt) << std::endl;
    os << "  q_terminal: " << Map<const VectorX>(q_terminal, 4).transpose().format(fmt) << std::endl;
    os << "  x_ref: " << Map<const VectorX>(x_ref, 4).transpose().format(fmt) << std::endl;
    os << "  r: " << Map<const VectorX>(r, 1).transpose().format(fmt) << std::endl;
    os << std::endl;
    os << "  ubound_indices: " << Map<const VectorXi>(ubound_indices, 1).transpose().format(intfmt) << std::endl;
    os << "  umin: " << Map<const VectorX>(umin, 1).transpose().format(fmt) << std::endl;
    os << "  umax: " << Map<const VectorX>(umax, 1).transpose().format(fmt) << std::endl;
    os << "  dummy_weight: " << Map<const VectorX>(dummy_weight, 1).transpose().format(fmt) << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const OCP_cartpole& ocp) { 
    ocp.disp(os);
    return os;
  }


  // Computes the state equation f(t, x, u).
  // t : time parameter
  // x : state vector
  // u : control input vector
  // f : the value of f(t, x, u)
  void eval_f(const double t, const double* x, const double* u, 
              double* dx) const {
  double x0 = sin(x[1]);
  double x1 = 1.0/(m_c + m_p*pow(x0, 2));
  double x2 = cos(x[1]);
  double x3 = l*pow(x[1], 2);
  double x4 = m_p*x0;
  dx[0] = x[2];
  dx[1] = x[3];
  dx[2] = x1*(u[0] + x4*(g*x2 + x3));
  dx[3] = x1*(-g*x0*(m_c + m_p) - u[0]*x2 - x2*x3*x4)/l;
 
  }

  // Computes the partial derivative of terminal cost with respect to state, 
  // i.e., dphi/dx(t, x).
  // t    : time parameter
  // x    : state vector
  // phix : the value of dphi/dx(t, x)
  void eval_phix(const double t, const double* x, double* phix) const {
  phix[0] = (1.0/2.0)*q_terminal[0]*(2*x[0] - 2*x_ref[0]);
  phix[1] = (1.0/2.0)*q_terminal[1]*(2*x[1] - 2*x_ref[1]);
  phix[2] = (1.0/2.0)*q_terminal[2]*(2*x[2] - 2*x_ref[2]);
  phix[3] = (1.0/2.0)*q_terminal[3]*(2*x[3] - 2*x_ref[3]);
 
  }

  // Computes the partial derivative of the Hamiltonian with respect to state, 
  // i.e., dH/dx(t, x, u, lmd).
  // t   : time parameter
  // x   : state vector
  // u   : control input vector
  // lmd : the Lagrange multiplier for the state equation
  // hx  : the value of dH/dx(t, x, u, lmd)
  void eval_hx(const double t, const double* x, const double* u, 
               const double* lmd, double* hx) const {
  double x0 = 2*x[1];
  double x1 = sin(x[1]);
  double x2 = cos(x[1]);
  double x3 = g*x2;
  double x4 = pow(x[1], 2);
  double x5 = l*x4;
  double x6 = m_p*(x3 + x5);
  double x7 = pow(x1, 2);
  double x8 = m_c + m_p*x7;
  double x9 = m_p*x1;
  double x10 = x2*x9;
  double x11 = 2*x10/pow(x8, 2);
  double x12 = 1.0/x8;
  double x13 = g*x1;
  double x14 = m_c + m_p;
  double x15 = lmd[3]/l;
  hx[0] = (1.0/2.0)*q[0]*(2*x[0] - 2*x_ref[0]);
  hx[1] = -lmd[2]*x11*(u[0] + x1*x6) + lmd[2]*x12*(x2*x6 + x9*(2*l*x[1] - x13)) + (1.0/2.0)*q[1]*(x0 - 2*x_ref[1]) - x11*x15*(-u[0]*x2 - x10*x5 - x13*x14) + x12*x15*(l*m_p*x4*x7 - l*x0*x10 - m_p*pow(x2, 2)*x5 + u[0]*x1 - x14*x3);
  hx[2] = lmd[0] + (1.0/2.0)*q[2]*(2*x[2] - 2*x_ref[2]);
  hx[3] = lmd[1] + (1.0/2.0)*q[3]*(2*x[3] - 2*x_ref[3]);
 
  }

  // Computes the partial derivative of the Hamiltonian with respect to control 
  // input and the constraints, dH/du(t, x, u, lmd).
  // t   : time parameter
  // x   : state vector
  // u   : control input vector
  // lmd : the Lagrange multiplier for the state equation
  // hu  : the value of dH/du(t, x, u, lmd)
  void eval_hu(const double t, const double* x, const double* u, 
               const double* lmd, double* hu) const {
  double x0 = 1.0/(m_c + m_p*pow(sin(x[1]), 2));
  hu[0] = lmd[2]*x0 + r[0]*u[0] - lmd[3]*x0*cos(x[1])/l;
 
  }
};

} // namespace cgmres

#endif // CGMRES_OCP_HPP_
