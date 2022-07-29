#ifndef CGMRES__OCP_MOBILEROBOT_HPP_ 
#define CGMRES__OCP_MOBILEROBOT_HPP_ 
 
#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include "cgmres/types.hpp"

namespace cgmres {

// This class defines the optimal control problem (OCP)
class OCP_mobilerobot { 
public:
  static constexpr int nx = 3;
  static constexpr int nu = 2;
  static constexpr int nc = 2;
  static constexpr int nuc = nu + nc;
  static constexpr int nub = 2;

  double vx_ref = 0.4;
  double v_min = -0.5;
  double v_max = 0.5;
  double w_min = -0.75;
  double w_max = 0.75;
  double X_1 = 1;
  double Y_1 = 0.25;
  double R_1 = 0.5;
  double X_2 = 2;
  double Y_2 = -0.25;
  double R_2 = 0.5;

  double q[3] = {10, 1, 0.01};
  double r[2] = {0.1, 0.1};
  double x_ref[3] = {0, 0, 0};

  static constexpr int ubound_indices[2] = {0, 1};
  double umin[2] = {-1.0, -1.0};
  double umax[2] = {1.0, 1.0};
  double dummy_weight[2] = {0.1, 0.1};

  double fb_eps[2] = {0.01, 0.01};

  void disp(std::ostream& os) const {
    os << "OCP_mobilerobot:" << std::endl;
    os << "  nx:  " << nx << std::endl;
    os << "  nu:  " << nu << std::endl;
    os << "  nc:  " << nc << std::endl;
    os << "  nuc: " << nuc << std::endl;
    os << "  nub: " << nub << std::endl;
    os << std::endl;
    os << "  vx_ref: " << vx_ref << std::endl;
    os << "  v_min: " << v_min << std::endl;
    os << "  v_max: " << v_max << std::endl;
    os << "  w_min: " << w_min << std::endl;
    os << "  w_max: " << w_max << std::endl;
    os << "  X_1: " << X_1 << std::endl;
    os << "  Y_1: " << Y_1 << std::endl;
    os << "  R_1: " << R_1 << std::endl;
    os << "  X_2: " << X_2 << std::endl;
    os << "  Y_2: " << Y_2 << std::endl;
    os << "  R_2: " << R_2 << std::endl;
    os << std::endl;
    Eigen::IOFormat fmt(4, 0, ", ", "", "[", "]");
    Eigen::IOFormat intfmt(1, 0, ", ", "", "[", "]");
    os << "  q: " << Map<const VectorX>(q, 3).transpose().format(fmt) << std::endl;
    os << "  r: " << Map<const VectorX>(r, 2).transpose().format(fmt) << std::endl;
    os << "  x_ref: " << Map<const VectorX>(x_ref, 3).transpose().format(fmt) << std::endl;
    os << std::endl;
    os << "  ubound_indices: " << Map<const VectorXi>(ubound_indices, 2).transpose().format(intfmt) << std::endl;
    os << "  umin: " << Map<const VectorX>(umin, 2).transpose().format(fmt) << std::endl;
    os << "  umax: " << Map<const VectorX>(umax, 2).transpose().format(fmt) << std::endl;
    os << "  dummy_weight: " << Map<const VectorX>(dummy_weight, 2).transpose().format(fmt) << std::endl;
    os << std::endl;
    os << "  fb_eps: " << Map<const VectorX>(fb_eps, 2).transpose().format(fmt) << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const OCP_mobilerobot& ocp) { 
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
  dx[0] = u[0]*cos(x[2]);
  dx[1] = u[0]*sin(x[2]);
  dx[2] = u[1];
 
  }

  // Computes the partial derivative of terminal cost with respect to state, 
  // i.e., dphi/dx(t, x).
  // t    : time parameter
  // x    : state vector
  // phix : the value of dphi/dx(t, x)
  void eval_phix(const double t, const double* x, double* phix) const {
  phix[0] = (1.0/2.0)*q[0]*(-2*t*vx_ref + 2*x[0]);
  phix[1] = q[1]*x[1];
  phix[2] = q[2]*x[2];
 
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
  double x0 = -2*x[0];
  double x1 = -2*x[1];
  double x2 = u[0]*sin(x[2]);
  double x3 = cos(x[2]);
  hx[0] = (1.0/2.0)*q[0]*(-2*t*vx_ref - x0) + u[2]*(2*X_1 + x0) + u[3]*(2*X_2 + x0);
  hx[1] = q[1]*x[1] + u[2]*(2*Y_1 + x1) + u[3]*(2*Y_2 + x1);
  hx[2] = -lmd[0]*x2 + lmd[1]*u[0]*x3 + q[2]*x[2] - r[0]*x2*(u[0]*x3 - vx_ref);
 
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
  double x0 = cos(x[2]);
  double x1 = -x[0];
  double x2 = -x[1];
  double x3 = -pow(R_1, 2) + pow(-X_1 - x1, 2) + pow(-Y_1 - x2, 2);
  double x4 = -pow(R_2, 2) + pow(-X_2 - x1, 2) + pow(-Y_2 - x2, 2);
  hu[0] = lmd[0]*x0 + lmd[1]*sin(x[2]) + r[0]*x0*(u[0]*x0 - vx_ref);
  hu[1] = lmd[2] + r[1]*u[1];
  hu[2] = -u[2] - x3 + sqrt(fb_eps[0] + pow(u[2], 2) + pow(x3, 2));
  hu[3] = -u[3] - x4 + sqrt(fb_eps[1] + pow(u[3], 2) + pow(x4, 2));
 
  }
};

} // namespace cgmres

#endif // CGMRES_OCP_HPP_
