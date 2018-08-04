#include "nmpc_model.hpp"


// state function f(t, x, u)
void nmpc_model::statefunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f)
{
	f[0] = x[2];
	f[1] = x[3];
	f[2] = (-sin(x[0] + x[1]) * cos(x[1]) * d2 * d2 * g * l1 * m2 * m2 - d2 * l1 * m2 * (m2 * l1 * d2 * sin(x[1]) * pow(x[2], 0.2e1) - u[0]) * cos(x[1]) + J2 * (-d2 * l1 * m2 * pow(x[2] + x[3], 0.2e1) * sin(x[1]) + (d1 * m1 + l1 * m2) * g * sin(x[0]) + u[0])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1));
	f[3] = (g * d2 * m2 * (m2 * l1 * d2 * cos(x[1]) + l1 * l1 * m2 + J1) * sin(x[0] + x[1]) - d2 * m2 * (-0.2e1 * d2 * (pow(x[2], 0.2e1) + x[2] * x[3] + pow(x[3], 0.2e1) / 0.2e1) * m2 * l1 * sin(x[1]) + (d1 * m1 + l1 * m2) * g * sin(x[0]) + 0.2e1 * u[0]) * l1 * cos(x[1]) + d2 * m2 * l1 * (pow(x[2], 0.2e1) * l1 * l1 * m2 + pow(x[2] + x[3], 0.2e1) * J2 + pow(x[2], 0.2e1) * J1) * sin(x[1]) - g * J2 * (d1 * m1 + l1 * m2) * sin(x[0]) - u[0] * (l1 * l1 * m2 + J1 + J2)) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1));
}


// partial derivative of terminal cost  dphi/dx(t, x, u)
void nmpc_model::phixfunc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> lmd)
{
	lmd[0] = (x[0] - xf[0]) * sf[0];
	lmd[1] = (x[1] - xf[1]) * sf[1];
	lmd[2] = (x[2] - xf[2]) * sf[2];
	lmd[3] = (x[3] - xf[3]) * sf[3];
}


//  partial derivative of the Hamiltonian  dh/dx(t, x, u)
void nmpc_model::hxfunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hx)
{
	hx[0] = (x[0] - xf[0]) * q[0] + lmd[2] * (-cos(x[0] + x[1]) * cos(x[1]) * d2 * d2 * g * l1 * m2 * m2 + J2 * (d1 * m1 + l1 * m2) * g * cos(x[0])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1)) + lmd[3] * (g * d2 * m2 * (m2 * l1 * d2 * cos(x[1]) + l1 * l1 * m2 + J1) * cos(x[0] + x[1]) - d2 * m2 * (d1 * m1 + l1 * m2) * g * cos(x[0]) * l1 * cos(x[1]) - J2 * (d1 * m1 + l1 * m2) * g * cos(x[0])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1));
	hx[1] = (x[1] - xf[1]) * q[1] + lmd[2] * (-cos(x[0] + x[1]) * cos(x[1]) * d2 * d2 * g * l1 * m2 * m2 + sin(x[0] + x[1]) * sin(x[1]) * d2 * d2 * g * l1 * m2 * m2 - d2 * d2 * l1 * l1 * m2 * m2 * pow(cos(x[1]), 0.2e1) * pow(x[2], 0.2e1) + d2 * l1 * m2 * (m2 * l1 * d2 * sin(x[1]) * pow(x[2], 0.2e1) - u[0]) * sin(x[1]) - J2 * d2 * l1 * m2 * pow(x[2] + x[3], 0.2e1) * cos(x[1])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1)) + 0.2e1 * lmd[2] * (-sin(x[0] + x[1]) * cos(x[1]) * d2 * d2 * g * l1 * m2 * m2 - d2 * l1 * m2 * (m2 * l1 * d2 * sin(x[1]) * pow(x[2], 0.2e1) - u[0]) * cos(x[1]) + J2 * (-d2 * l1 * m2 * pow(x[2] + x[3], 0.2e1) * sin(x[1]) + (d1 * m1 + l1 * m2) * g * sin(x[0]) + u[0])) * pow(pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1), -0.2e1) * cos(x[1]) * d2 * d2 * l1 * l1 * m2 * m2 * sin(x[1]) + lmd[3] * (-sin(x[0] + x[1]) * sin(x[1]) * d2 * d2 * g * l1 * m2 * m2 + g * d2 * m2 * (m2 * l1 * d2 * cos(x[1]) + l1 * l1 * m2 + J1) * cos(x[0] + x[1]) + 0.2e1 * d2 * d2 * m2 * m2 * (pow(x[2], 0.2e1) + x[2] * x[3] + pow(x[3], 0.2e1) / 0.2e1) * l1 * l1 * pow(cos(x[1]), 0.2e1) + d2 * m2 * (-0.2e1 * d2 * (pow(x[2], 0.2e1) + x[2] * x[3] + pow(x[3], 0.2e1) / 0.2e1) * m2 * l1 * sin(x[1]) + (d1 * m1 + l1 * m2) * g * sin(x[0]) + 0.2e1 * u[0]) * l1 * sin(x[1]) + d2 * m2 * l1 * (pow(x[2], 0.2e1) * l1 * l1 * m2 + pow(x[2] + x[3], 0.2e1) * J2 + pow(x[2], 0.2e1) * J1) * cos(x[1])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1)) + 0.2e1 * lmd[3] * (g * d2 * m2 * (m2 * l1 * d2 * cos(x[1]) + l1 * l1 * m2 + J1) * sin(x[0] + x[1]) - d2 * m2 * (-0.2e1 * d2 * (pow(x[2], 0.2e1) + x[2] * x[3] + pow(x[3], 0.2e1) / 0.2e1) * m2 * l1 * sin(x[1]) + (d1 * m1 + l1 * m2) * g * sin(x[0]) + 0.2e1 * u[0]) * l1 * cos(x[1]) + d2 * m2 * l1 * (pow(x[2], 0.2e1) * l1 * l1 * m2 + pow(x[2] + x[3], 0.2e1) * J2 + pow(x[2], 0.2e1) * J1) * sin(x[1]) - g * J2 * (d1 * m1 + l1 * m2) * sin(x[0]) - u[0] * (l1 * l1 * m2 + J1 + J2)) * pow(pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1), -0.2e1) * cos(x[1]) * d2 * d2 * l1 * l1 * m2 * m2 * sin(x[1]);
	hx[2] = (x[2] - xf[2]) * q[2] + lmd[0] + lmd[2] * (-0.2e1 * d2 * d2 * l1 * l1 * m2 * m2 * sin(x[1]) * x[2] * cos(x[1]) - 0.2e1 * J2 * d2 * l1 * m2 * (x[2] + x[3]) * sin(x[1])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1)) + lmd[3] * (0.2e1 * d2 * d2 * m2 * m2 * (0.2e1 * x[2] + x[3]) * l1 * l1 * sin(x[1]) * cos(x[1]) + d2 * m2 * l1 * (0.2e1 * x[2] * l1 * l1 * m2 + 0.2e1 * (x[2] + x[3]) * J2 + 0.2e1 * x[2] * J1) * sin(x[1])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1));
	hx[3] = (x[3] - xf[3]) * q[3] + lmd[1] - 0.2e1 * lmd[2] * J2 * d2 * l1 * m2 * (x[2] + x[3]) * sin(x[1]) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1)) + lmd[3] * (0.2e1 * d2 * d2 * m2 * m2 * l1 * l1 * (x[2] + x[3]) * sin(x[1]) * cos(x[1]) + 0.2e1 * J2 * d2 * l1 * m2 * (x[2] + x[3]) * sin(x[1])) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1));
}


//  partial derivative of the Hamiltonian dh/du(t, x, u)
void nmpc_model::hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hu)
{
	hu[0] = r[0] * u[0] + lmd[2] * (m2 * l1 * d2 * cos(x[1]) + J2) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1)) + lmd[3] * (-0.2e1 * m2 * l1 * d2 * cos(x[1]) - l1 * l1 * m2 - J1 - J2) / (pow(cos(x[1]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 - J2 * (l1 * l1 * m2 + J1));
}