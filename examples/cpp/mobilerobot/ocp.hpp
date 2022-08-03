#ifndef CGMRES__OCP_MOBILEROBOT_HPP_ 
#define CGMRES__OCP_MOBILEROBOT_HPP_ 
 
#define _USE_MATH_DEFINES

#include <cmath>
#include <array>
#include <iostream>

#include "cgmres/types.hpp"

namespace cgmres {

/// 
/// @class OCP_mobilerobot
/// @brief Definition of the optimal control problem (OCP) of mobilerobot.
/// 
class OCP_mobilerobot { 
public:
 
  ///
  /// @brief Dimension of the state. 
  ///
  static constexpr int nx = 3;
 
  ///
  /// @brief Dimension of the control input. 
  ///
  static constexpr int nu = 2;
 
  ///
  /// @brief Dimension of the equality constraints. 
  ///
  static constexpr int nc = 2;
 
  ///
  /// @brief Dimension of the Fischer-Burmeister function (already counded in nc). 
  ///
  static constexpr int nh = 2;
 
  ///
  /// @brief Dimension of the concatenation of the control input and equality constraints. 
  ///
  static constexpr int nuc = nu + nc;
 
  ///
  /// @brief Dimension of the bound constraints on the control input. 
  ///
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

  std::array<double, 3> q = {10, 1, 0.01};
  std::array<double, 2> r = {0.1, 0.1};
  std::array<double, 3> x_ref = {0, 0, 0};

  static constexpr std::array<int, nub> ubound_indices = {0, 1};
  std::array<double, nub> umin = {-1.0, -1.0};
  std::array<double, nub> umax = {1.0, 1.0};
  std::array<double, nub> dummy_weight = {0.1, 0.1};

  std::array<double, nh> fb_eps = {0.01, 0.01};

  void disp(std::ostream& os) const {
    os << "OCP_mobilerobot:" << std::endl;
    os << "  nx:  " << nx << std::endl;
    os << "  nu:  " << nu << std::endl;
    os << "  nc:  " << nc << std::endl;
    os << "  nh:  " << nh << std::endl;
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
    os << "  q: " << Map<const VectorX>(q.data(), q.size()).transpose().format(fmt) << std::endl;
    os << "  r: " << Map<const VectorX>(r.data(), r.size()).transpose().format(fmt) << std::endl;
    os << "  x_ref: " << Map<const VectorX>(x_ref.data(), x_ref.size()).transpose().format(fmt) << std::endl;
    os << std::endl;
    os << "  ubound_indices: " << Map<const VectorXi>(ubound_indices.data(), ubound_indices.size()).transpose().format(intfmt) << std::endl;
    os << "  umin: " << Map<const VectorX>(umin.data(), umin.size()).transpose().format(fmt) << std::endl;
    os << "  umax: " << Map<const VectorX>(umax.data(), umax.size()).transpose().format(fmt) << std::endl;
    os << "  dummy_weight: " << Map<const VectorX>(dummy_weight.data(), dummy_weight.size()).transpose().format(fmt) << std::endl;
    os << std::endl;
    os << "  fb_eps: " << Map<const VectorX>(fb_eps.data(), fb_eps.size()).transpose().format(fmt) << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const OCP_mobilerobot& ocp) { 
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
    dx[0] = u[0]*cos(x[2]);
    dx[1] = u[0]*sin(x[2]);
    dx[2] = u[1];
 
  }

  ///
  /// Computes the partial derivative of terminal cost with respect to state, 
  /// i.e., phix = dphi/dx(t, x).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[out] phix Evaluated value of the partial derivative of terminal cost.
  ///
  void eval_phix(const double t, const double* x, double* phix) const {
    phix[0] = (1.0/2.0)*q[0]*(-2*t*vx_ref + 2*x[0]);
    phix[1] = q[1]*x[1];
    phix[2] = q[2]*x[2];
 
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
    const double x0 = -2*x[0];
    const double x1 = -2*x[1];
    const double x2 = u[0]*sin(x[2]);
    const double x3 = cos(x[2]);
    hx[0] = (1.0/2.0)*q[0]*(-2*t*vx_ref - x0) + u[2]*(2*X_1 + x0) + u[3]*(2*X_2 + x0);
    hx[1] = q[1]*x[1] + u[2]*(2*Y_1 + x1) + u[3]*(2*Y_2 + x1);
    hx[2] = -lmd[0]*x2 + lmd[1]*u[0]*x3 + q[2]*x[2] - r[0]*x2*(u[0]*x3 - vx_ref);
 
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
    const double x0 = cos(x[2]);
    const double x1 = -x[0];
    const double x2 = -x[1];
    const double x3 = -pow(R_1, 2) + pow(-X_1 - x1, 2) + pow(-Y_1 - x2, 2);
    const double x4 = -pow(R_2, 2) + pow(-X_2 - x1, 2) + pow(-Y_2 - x2, 2);
    hu[0] = lmd[0]*x0 + lmd[1]*sin(x[2]) + r[0]*x0*(u[0]*x0 - vx_ref);
    hu[1] = lmd[2] + r[1]*u[1];
    hu[2] = -u[2] - x3 + sqrt(fb_eps[0] + pow(u[2], 2) + pow(x3, 2));
    hu[3] = -u[3] - x4 + sqrt(fb_eps[1] + pow(u[3], 2) + pow(x4, 2));
 
  }
};

} // namespace cgmres

#endif // CGMRES_OCP_HPP_
