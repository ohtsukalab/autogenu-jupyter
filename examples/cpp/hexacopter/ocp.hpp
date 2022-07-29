#ifndef CGMRES__OCP_HEXACOPTER_HPP_ 
#define CGMRES__OCP_HEXACOPTER_HPP_ 
 
#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include "cgmres/types.hpp"

namespace cgmres {

// This class defines the optimal control problem (OCP)
class OCP_hexacopter { 
public:
  static constexpr int nx = 12;
  static constexpr int nu = 6;
  static constexpr int nc = 0;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = 6;

  double m = 1.44;
  double l = 0.23;
  double k = 1.6e-09;
  double Ixx = 0.0348;
  double Iyy = 0.0459;
  double Izz = 0.0977;
  double gamma = 0.01;
  double g = 9.80665;
  double z_ref = 5;

  double q[12] = {1, 1, 1, 0.01, 0.01, 0, 0.01, 0.01, 0.01, 0.1, 0.1, 0.001};
  double q_terminal[12] = {1, 1, 1, 0.01, 0.01, 0, 0.01, 0.01, 0.01, 0.1, 0.1, 0.001};
  double r[6] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

  static constexpr int ubound_indices[6] = {0, 1, 2, 3, 4, 5};
  double umin[6] = {0.144, 0.144, 0.144, 0.144, 0.144, 0.144};
  double umax[6] = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
  double dummy_weight[6] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  void disp(std::ostream& os) const {
    os << "OCP_hexacopter:" << std::endl;
    os << "  nx:  " << nx << std::endl;
    os << "  nu:  " << nu << std::endl;
    os << "  nc:  " << nc << std::endl;
    os << "  nuc: " << nuc << std::endl;
    os << "  nub: " << nub << std::endl;
    os << std::endl;
    os << "  m: " << m << std::endl;
    os << "  l: " << l << std::endl;
    os << "  k: " << k << std::endl;
    os << "  Ixx: " << Ixx << std::endl;
    os << "  Iyy: " << Iyy << std::endl;
    os << "  Izz: " << Izz << std::endl;
    os << "  gamma: " << gamma << std::endl;
    os << "  g: " << g << std::endl;
    os << "  z_ref: " << z_ref << std::endl;
    os << std::endl;
    Eigen::IOFormat fmt(4, 0, ", ", "", "[", "]");
    Eigen::IOFormat intfmt(1, 0, ", ", "", "[", "]");
    os << "  q: " << Map<const VectorX>(q, 12).transpose().format(fmt) << std::endl;
    os << "  q_terminal: " << Map<const VectorX>(q_terminal, 12).transpose().format(fmt) << std::endl;
    os << "  r: " << Map<const VectorX>(r, 6).transpose().format(fmt) << std::endl;
    os << std::endl;
    os << "  ubound_indices: " << Map<const VectorXi>(ubound_indices, 6).transpose().format(intfmt) << std::endl;
    os << "  umin: " << Map<const VectorX>(umin, 6).transpose().format(fmt) << std::endl;
    os << "  umax: " << Map<const VectorX>(umax, 6).transpose().format(fmt) << std::endl;
    os << "  dummy_weight: " << Map<const VectorX>(dummy_weight, 6).transpose().format(fmt) << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const OCP_hexacopter& ocp) { 
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
  double x0 = sin(x[3]);
  double x1 = sin(x[5]);
  double x2 = cos(x[5]);
  double x3 = cos(x[3]);
  double x4 = sin(x[4]);
  double x5 = 1.0/m;
  double x6 = u[0] + u[2] + u[4];
  double x7 = u[1] + u[3] + u[5] + x6;
  double x8 = x5*x7;
  double x9 = 1.0/Ixx;
  double x10 = -Izz;
  double x11 = (1.0/2.0)*u[0];
  double x12 = 1.0/Iyy;
  double x13 = sqrt(3);
  double x14 = 1.0/Izz;
  dx[0] = x[6];
  dx[1] = x[7];
  dx[2] = x[8];
  dx[3] = x[9];
  dx[4] = x[10];
  dx[5] = x[11];
  dx[6] = x8*(x0*x1 + x2*x3*x4);
  dx[7] = x8*(-x0*x2 + x1*x3*x4);
  dx[8] = -g + x3*x5*x7*cos(x[4]);
  dx[9] = l*x9*(-u[1] - 1.0/2.0*u[2] + (1.0/2.0)*u[3] + u[4] + (1.0/2.0)*u[5] - x11) + x9*x[10]*x[11]*(Iyy + x10);
  dx[10] = l*x12*((1.0/2.0)*u[2]*x13 + (1.0/2.0)*u[3]*x13 - 1.0/2.0*u[5]*x13 - x11*x13) + x12*x[11]*x[9]*(-Ixx - x10);
  dx[11] = x14*x[10]*x[9]*(Ixx - Iyy) + x14*(-gamma*x[11] + k*(u[1] + u[3] + u[5] - x6));
 
  }

  // Computes the partial derivative of terminal cost with respect to state, 
  // i.e., dphi/dx(t, x).
  // t    : time parameter
  // x    : state vector
  // phix : the value of dphi/dx(t, x)
  void eval_phix(const double t, const double* x, double* phix) const {
  double x0 = 2*t;
  double x1 = sin(x0);
  double x2 = cos(x0);
  phix[0] = (1.0/2.0)*q_terminal[0]*(-2*x1 + 2*x[0]);
  phix[1] = (1.0/2.0)*q_terminal[1]*(2*x2 + 2*x[1] - 2);
  phix[2] = (1.0/2.0)*q_terminal[2]*(2*x[2] - 2*z_ref - 4*sin(t));
  phix[3] = q_terminal[3]*x[3];
  phix[4] = q_terminal[4]*x[4];
  phix[5] = q_terminal[5]*x[5];
  phix[6] = (1.0/2.0)*q_terminal[6]*(-4*x2 + 2*x[6]);
  phix[7] = (1.0/2.0)*q_terminal[7]*(-4*x1 + 2*x[7]);
  phix[8] = (1.0/2.0)*q_terminal[8]*(2*x[8] - 4*cos(t));
  phix[9] = q_terminal[9]*x[9];
  phix[10] = q_terminal[10]*x[10];
  phix[11] = q_terminal[11]*x[11];
 
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
  double x0 = 2*t;
  double x1 = sin(x0);
  double x2 = cos(x0);
  double x3 = sin(x[3]);
  double x4 = cos(x[4]);
  double x5 = (u[0] + u[1] + u[2] + u[3] + u[4] + u[5])/m;
  double x6 = lmd[8]*x5;
  double x7 = sin(x[5]);
  double x8 = cos(x[3]);
  double x9 = sin(x[4]);
  double x10 = cos(x[5]);
  double x11 = x10*x3;
  double x12 = lmd[6]*x5;
  double x13 = x10*x8;
  double x14 = x3*x7;
  double x15 = lmd[7]*x5;
  double x16 = x7*x8;
  double x17 = -Izz;
  double x18 = lmd[10]*(-Ixx - x17)/Iyy;
  double x19 = lmd[11]/Izz;
  double x20 = x19*(Ixx - Iyy);
  double x21 = lmd[9]*(Iyy + x17)/Ixx;
  hx[0] = (1.0/2.0)*q[0]*(-2*x1 + 2*x[0]);
  hx[1] = (1.0/2.0)*q[1]*(2*x2 + 2*x[1] - 2);
  hx[2] = (1.0/2.0)*q[2]*(2*x[2] - 2*z_ref - 4*sin(t));
  hx[3] = q[3]*x[3] + x12*(-x11*x9 + x7*x8) + x15*(-x13 - x14*x9) - x3*x4*x6;
  hx[4] = q[4]*x[4] + x12*x13*x4 + x15*x16*x4 - x6*x8*x9;
  hx[5] = q[5]*x[5] + x12*(x11 - x16*x9) + x15*(x13*x9 + x14);
  hx[6] = lmd[0] + (1.0/2.0)*q[6]*(-4*x2 + 2*x[6]);
  hx[7] = lmd[1] + (1.0/2.0)*q[7]*(-4*x1 + 2*x[7]);
  hx[8] = lmd[2] + (1.0/2.0)*q[8]*(2*x[8] - 4*cos(t));
  hx[9] = lmd[3] + q[9]*x[9] + x18*x[11] + x20*x[10];
  hx[10] = lmd[4] + q[10]*x[10] + x20*x[9] + x21*x[11];
  hx[11] = -gamma*x19 + lmd[5] + q[11]*x[11] + x18*x[9] + x21*x[10];
 
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
  double x0 = (1.0/3.0)*g*m;
  double x1 = (1.0/2.0)*sqrt(3)*l*lmd[10]/Iyy;
  double x2 = -x1;
  double x3 = l*lmd[9]/Ixx;
  double x4 = (1.0/2.0)*x3;
  double x5 = k*lmd[11]/Izz;
  double x6 = 1.0/m;
  double x7 = sin(x[3]);
  double x8 = sin(x[5]);
  double x9 = cos(x[5]);
  double x10 = cos(x[3]);
  double x11 = sin(x[4]);
  double x12 = lmd[6]*x6*(x10*x11*x9 + x7*x8) + lmd[7]*x6*(x10*x11*x8 - x7*x9) + lmd[8]*x10*x6*cos(x[4]);
  double x13 = x12 - x5;
  double x14 = x13 - x4;
  double x15 = x12 + x5;
  double x16 = x15 + x4;
  hu[0] = (1.0/2.0)*r[0]*(2*u[0] - x0) + x14 + x2;
  hu[1] = (1.0/2.0)*r[1]*(2*u[1] - x0) + x15 - x3;
  hu[2] = (1.0/2.0)*r[2]*(2*u[2] - x0) + x1 + x14;
  hu[3] = (1.0/2.0)*r[3]*(2*u[3] - x0) + x1 + x16;
  hu[4] = (1.0/2.0)*r[4]*(2*u[4] - x0) + x13 + x3;
  hu[5] = (1.0/2.0)*r[5]*(2*u[5] - x0) + x16 + x2;
 
  }
};

} // namespace cgmres

#endif // CGMRES_OCP_HPP_
