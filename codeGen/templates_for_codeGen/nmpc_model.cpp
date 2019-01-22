#include "nmpc_model.hpp"


// State equation f(t, x, u)
// t : time parameter
// x : state vector
// u : control input vector
// f : the value of f(t, x, u)
void NMPCModel::stateFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f)
{

}


// Partial derivative of terminal cost with respect to state, dphi/dx(t, x)
// phi(t, x) = (q_terminal[0]*(x[0]-x_ref[0])^2 + q_terminal[1]*(x[1]-x_ref[1])^2 + q_terminal[2]*(x[2]-x_ref[2])^2 + q_terminal[3]*(x[3]-x_ref[3])^2)/2
// t    : time parameter
// x    : state vector
// phix : the value of dphi/dx(t, x)
void NMPCModel::phixFunc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> phix)
{

}


// Partial derivative of the Hamiltonian with respect to state, dH/dx(t, x, u, lmd)
// H(t, x, u, lmd) = L(t, x, u) + lmd * f(t, x, u)
// L(t, x, u) = (q[0]*(x[0]-x_ref[0])^2 + q[1]*(x[1]-x_ref[1])^2 + q[2]*(x[2]-x_ref[2])^2 + q[3]*(x[3]-x_ref[3])^2 + r[0]*u[0]^2)/2
// t   : time parameter
// x   : state vector
// u   : control input vector
// lmd : the Lagrange multiplier for the state equation
// hx  : the value of dH/dx(t, x, u, lmd)
void NMPCModel::hxFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hx)
{

}


// Partial derivative of the Hamiltonian with respect to control input and constraints, dH/du(t, x, u, lmd)
// H(t, x, u, lmd) = L(t, x, u) + lmd * f(t, x, u)
// L(t, x, u) = (q[0]*(x[0]-x_ref[0])^2 + q[1]*(x[1]-x_ref[1])^2 + q[2]*(x[2]-x_ref[2])^2 + q[3]*(x[3]-x_ref[3])^2 + r[0]*u[0]^2)/2
// t   : time parameter
// x   : state vector
// u   : control input vector
// lmd : the Lagrange multiplier for the state equation
// hu  : the value of dH/du(t, x, u, lmd)
void NMPCModel::huFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hu)
{

}