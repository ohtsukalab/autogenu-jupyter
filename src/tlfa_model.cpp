// Two link, fully actuated model
#include "tlfa_model.hpp"



// State equation f(t, x, u)
// t : time parameter
// x : state vector
// u : control input vector
// f : the value of f(t, x, u)
void NMPCModel::stateFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f)
{
    f[0] = x[2];
    f[1] = x[3];
    f[2] = -(sin(x[0] + x[1]) * d2 * g * m2 + g * (d1 * m1 + l1 * m2) * sin(x[0]) - 0.2e1 * d2 * (x[2] + x[3] / 0.2e1) * l1 * x[3] * m2 * sin(x[1]) - u[0]) / (0.2e1 * d2 * m2 * l1 * cos(x[1]) + d1 * d1 * m1 + d2 * d2 * m2 + l1 * l1 * m2 + J1 + J2);
    f[3] = (g * d2 * l1 * m2 * (d1 * m1 + l1 * m2) * sin(x[0] - x[1]) / 0.2e1 - d2 * d2 * g * l1 * m2 * m2 * sin(x[0] + 0.2e1 * x[1]) / 0.2e1 - (d1 * d1 * m1 - d1 * l1 * m1 / 0.2e1 + l1 * l1 * m2 / 0.2e1 + J1) * m2 * g * d2 * sin(x[0] + x[1]) - l1 * l1 * m2 * m2 * d2 * d2 * (pow(x[2], 0.2e1) + x[2] * x[3] + pow(x[3], 0.2e1) / 0.2e1) * sin(0.2e1 * x[1]) - l1 * m2 * d2 * ((pow(x[2] + x[3], 0.2e1) * d2 * d2 + pow(x[2], 0.2e1) * l1 * l1) * m2 + (d1 * d1 * m1 + J1 + J2) * pow(x[2], 0.2e1) + 0.2e1 * J2 * x[2] * x[3] + J2 * pow(x[3], 0.2e1)) * sin(x[1]) + g * (d2 * d2 * l1 * m2 * m2 / 0.2e1 + (d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * sin(x[0]) + 0.2e1 * l1 * m2 * (u[1] - u[0] / 0.2e1) * d2 * cos(x[1]) + ((0.2e1 * u[1] - 0.2e1 * u[0]) * d2 * d2 / 0.2e1 + l1 * l1 * u[1]) * m2 + (u[1] - u[0]) * J2 + u[1] * (d1 * d1 * m1 + J1)) / (0.2e1 * d2 * m2 * l1 * cos(x[1]) + d1 * d1 * m1 + d2 * d2 * m2 + l1 * l1 * m2 + J1 + J2) / (d2 * d2 * m2 + J2);
}


// Partial derivative of terminal cost with respect to state, dphi/dx(t, x)
// phi(t, x) = (q_terminal[0]*(x[0]-x_ref[0])^2 + q_terminal[1]*(x[1]-x_ref[1])^2 + q_terminal[2]*(x[2]-x_ref[2])^2 + q_terminal[3]*(x[3]-x_ref[3])^2)/2
// t    : time parameter
// x    : state vector
// phix : the value of dphi/dx(t, x)
void NMPCModel::phixFunc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> phix)
{
    phix[0] = (x[0] - x_ref[0]) * q_terminal[0];
    phix[1] = (x[1] - x_ref[1]) * q_terminal[1];
    phix[2] = (x[2] - x_ref[2]) * q_terminal[2];
    phix[3] = (x[3] - x_ref[3]) * q_terminal[3];
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
    hx[0] = (g * lmd[3] * d2 * l1 * m2 * (d1 * m1 + l1 * m2) * cos(x[0] - x[1]) - cos(x[0] + 0.2e1 * x[1]) * d2 * d2 * g * l1 * m2 * m2 * lmd[3] - 0.2e1 * m2 * g * ((lmd[2] * d2 * d2 + lmd[3] * l1 * l1 / 0.2e1) * m2 + lmd[2] * J2 + lmd[3] * (d1 * d1 * m1 - d1 * l1 * m1 / 0.2e1 + J1)) * d2 * cos(x[0] + x[1]) - 0.2e1 * (d2 * d2 * l1 * (lmd[2] - lmd[3] / 0.2e1) * m2 * m2 + (d1 * d2 * d2 * m1 + J2 * l1) * (lmd[2] - lmd[3]) * m2 + J2 * d1 * m1 * (lmd[2] - lmd[3])) * g * cos(x[0]) + 0.4e1 * (x[0] - x_ref[0]) * q[0] * (d2 * m2 * l1 * cos(x[1]) + (d2 * d2 / 0.2e1 + l1 * l1 / 0.2e1) * m2 + d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 + J2 / 0.2e1) * (d2 * d2 * m2 + J2)) / (0.2e1 * d2 * m2 * l1 * cos(x[1]) + d1 * d1 * m1 + d2 * d2 * m2 + l1 * l1 * m2 + J1 + J2) / (d2 * d2 * m2 + J2) / 0.2e1;
    hx[1] = (-0.2e1 * m2 * (d2 * m2 * l1 * cos(x[1]) + (d2 * d2 / 0.2e1 + l1 * l1 / 0.2e1) * m2 + d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 + J2 / 0.2e1) * lmd[3] * l1 * (d1 * m1 + l1 * m2) * g * d2 * cos(x[0] - x[1]) - 0.4e1 * m2 * m2 * (d2 * m2 * l1 * cos(x[1]) + (d2 * d2 / 0.2e1 + l1 * l1 / 0.2e1) * m2 + d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 + J2 / 0.2e1) * lmd[3] * l1 * g * d2 * d2 * cos(x[0] + 0.2e1 * x[1]) + 0.2e1 * g * lmd[3] * sin(x[1]) * d2 * d2 * l1 * l1 * m2 * m2 * (d1 * m1 + l1 * m2) * sin(x[0] - x[1]) - 0.2e1 * sin(x[0] + 0.2e1 * x[1]) * sin(x[1]) * pow(d2, 0.3e1) * g * l1 * l1 * pow(m2, 0.3e1) * lmd[3] - 0.4e1 * m2 * (d2 * m2 * l1 * cos(x[1]) + (d2 * d2 / 0.2e1 + l1 * l1 / 0.2e1) * m2 + d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 + J2 / 0.2e1) * ((lmd[2] * d2 * d2 + lmd[3] * l1 * l1 / 0.2e1) * m2 - lmd[3] * d1 * l1 * m1 / 0.2e1 + (d1 * d1 * m1 + J1) * lmd[3] + lmd[2] * J2) * g * d2 * cos(x[0] + x[1]) - 0.8e1 * m2 * m2 * (d2 * m2 * l1 * cos(x[1]) + (d2 * d2 / 0.2e1 + l1 * l1 / 0.2e1) * m2 + d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 + J2 / 0.2e1) * lmd[3] * l1 * l1 * (pow(x[2], 0.2e1) + x[2] * x[3] + pow(x[3], 0.2e1) / 0.2e1) * d2 * d2 * cos(0.2e1 * x[1]) - 0.4e1 * m2 * m2 * ((lmd[2] * d2 * d2 + lmd[3] * l1 * l1 / 0.2e1) * m2 - lmd[3] * d1 * l1 * m1 / 0.2e1 + (d1 * d1 * m1 + J1) * lmd[3] + lmd[2] * J2) * l1 * sin(x[1]) * g * d2 * d2 * sin(x[0] + x[1]) - 0.4e1 * pow(m2, 0.3e1) * lmd[3] * pow(l1, 0.3e1) * (pow(x[2], 0.2e1) + x[2] * x[3] + pow(x[3], 0.2e1) / 0.2e1) * sin(x[1]) * pow(d2, 0.3e1) * sin(0.2e1 * x[1]) + 0.8e1 * d2 * d2 * l1 * l1 * m2 * m2 * q[1] * (d2 * d2 * m2 + J2) * (x[1] - x_ref[1]) * pow(cos(x[1]), 0.2e1) + 0.8e1 * m2 * l1 * (((-pow(x[2] + x[3], 0.2e1) * lmd[3] / 0.4e1 + q[1] * (x[1] - x_ref[1]) + (x[2] + x[3] / 0.2e1) * x[3] * lmd[2] / 0.2e1) * d2 * d2 - pow(x[2], 0.2e1) * lmd[3] * l1 * l1 / 0.4e1) * m2 + (-pow(x[2] + x[3], 0.2e1) * J2 / 0.4e1 - pow(x[2], 0.2e1) * (d1 * d1 * m1 + J1) / 0.4e1) * lmd[3] + (q[1] * (x[1] - x_ref[1]) + (x[2] + x[3] / 0.2e1) * x[3] * lmd[2] / 0.2e1) * J2) * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * d2 * cos(x[1]) - 0.4e1 * m2 * ((d2 * d2 * l1 * (lmd[2] - lmd[3] / 0.2e1) * m2 * m2 + (d1 * d2 * d2 * m1 + J2 * l1) * (lmd[2] - lmd[3]) * m2 + J2 * d1 * m1 * (lmd[2] - lmd[3])) * g * sin(x[0]) - u[0] * (((lmd[2] - lmd[3] / 0.2e1) * d2 * d2 + lmd[3] * l1 * l1 / 0.2e1) * m2 + (d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * lmd[3] + lmd[2] * J2)) * l1 * d2 * sin(x[1]) + 0.2e1 * ((x[1] - x_ref[1]) * q[1] * pow(d2, 0.4e1) + 0.2e1 * l1 * l1 * (-pow(x[2] + x[3], 0.2e1) * lmd[3] + q[1] * (x[1] - x_ref[1]) + 0.2e1 * (x[2] + x[3] / 0.2e1) * x[3] * lmd[2]) * d2 * d2 + pow(l1, 0.4e1) * (-0.2e1 * pow(x[2], 0.2e1) * lmd[3] + q[1] * (x[1] - x_ref[1]))) * d2 * d2 * pow(m2, 0.3e1) + (0.4e1 * (d1 * d1 * m1 + J1 + 0.3e1 / 0.2e1 * J2) * q[1] * (x[1] - x_ref[1]) * pow(d2, 0.4e1) + 0.4e1 * ((-pow(x[2] + x[3], 0.2e1) * J2 - pow(x[2], 0.2e1) * (d1 * d1 * m1 + J1)) * lmd[3] + ((0.2e1 * x[1] - 0.2e1 * x_ref[1]) * q[1] + 0.2e1 * (x[2] + x[3] / 0.2e1) * x[3] * lmd[2]) * J2 + q[1] * (d1 * d1 * m1 + J1) * (x[1] - x_ref[1])) * l1 * l1 * d2 * d2 + 0.2e1 * J2 * pow(l1, 0.4e1) * q[1] * (x[1] - x_ref[1])) * m2 * m2 + 0.2e1 * ((d1 * d1 * m1 + J1 + 0.3e1 * J2) * d2 * d2 + 0.2e1 * l1 * l1 * J2) * q[1] * (d1 * d1 * m1 + J1 + J2) * (x[1] - x_ref[1]) * m2 + 0.2e1 * J2 * q[1] * pow(d1 * d1 * m1 + J1 + J2, 0.2e1) * (x[1] - x_ref[1])) / (d2 * d2 * m2 + J2) * pow(0.2e1 * d2 * m2 * l1 * cos(x[1]) + d1 * d1 * m1 + d2 * d2 * m2 + l1 * l1 * m2 + J1 + J2, -0.2e1) / 0.2e1;
    hx[2] = (-0.2e1 * m2 * (0.2e1 * lmd[3] * m2 * (x[2] + x[3] / 0.2e1) * d2 * l1 * cos(x[1]) + ((x[2] * lmd[3] - x[3] * (lmd[2] - lmd[3])) * d2 * d2 + x[2] * lmd[3] * l1 * l1) * m2 + (x[2] * lmd[3] - x[3] * (lmd[2] - lmd[3])) * J2 + x[2] * lmd[3] * (d1 * d1 * m1 + J1)) * d2 * l1 * sin(x[1]) + (0.2e1 * d2 * m2 * l1 * cos(x[1]) + (d2 * d2 + l1 * l1) * m2 + J2 + d1 * d1 * m1 + J1) * (d2 * d2 * m2 + J2) * (q[2] * (x[2] - x_ref[2]) + lmd[0])) / (0.2e1 * d2 * m2 * l1 * cos(x[1]) + (d2 * d2 + l1 * l1) * m2 + J2 + d1 * d1 * m1 + J1) / (d2 * d2 * m2 + J2);
    hx[3] = (0.2e1 * (-lmd[3] * d2 * l1 * m2 * cos(x[1]) + (d2 * d2 * m2 + J2) * (lmd[2] - lmd[3])) * (x[2] + x[3]) * m2 * d2 * l1 * sin(x[1]) + (0.2e1 * d2 * m2 * l1 * cos(x[1]) + (d2 * d2 + l1 * l1) * m2 + J2 + d1 * d1 * m1 + J1) * (q[3] * (x[3] - x_ref[3]) + lmd[1]) * (d2 * d2 * m2 + J2)) / (0.2e1 * d2 * m2 * l1 * cos(x[1]) + (d2 * d2 + l1 * l1) * m2 + J2 + d1 * d1 * m1 + J1) / (d2 * d2 * m2 + J2);
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
    hu[0] = (- d2 * l1 *  m2 *  (-2 * d2 * d2 * m2 * r[0] * u[0] - 2 * J2 * r[0] * u[0] + lmd[3]) * cos(x[1]) + ( r[0] *  u[0] * ( (d2 * d2) + l1 * l1) *  m2 +  (u[0] * (d1 * d1 * m1 + J1 + J2) * r[0]) +  lmd[2] -  lmd[3]) *  (d2 * d2 * m2 + J2)) /  (d2 * d2 * m2 + J2) / (0.2e1 *  d2 *  m2 * l1 * cos(x[1]) + ( (d2 * d2) + l1 * l1) *  m2 +  J2 +  (d1 * d1 * m1) +  J1);
    hu[1] = (r[1] * u[1] * (d2 * d2 * m2 + J2) + lmd[3]) / (d2 * d2 * m2 + J2);
}



int NMPCModel::dimState() const
{
    return dim_state_;
}


int NMPCModel::dimControlInput() const
{
    return dim_control_input_;
}


int NMPCModel::dimConstraints() const
{
    return dim_constraints_;
}