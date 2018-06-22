#include "nmpc_model.hpp"



void nmpc_model::phix(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> lmd)
{
	lmd[0] = (x[0] - xf[0]) * sf[0];
	lmd[1] = (x[1] - xf[1]) * sf[1];
}

void nmpc_model::statefunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f)
{
	f[0] = x[1];
	f[1] = b * u[0] * x[1] + a * x[0];

}


void nmpc_model::hxfunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hx)
{
	hx[0] = a * lmd[1] + q[0] * x[0];
	hx[1] = b * lmd[1] * u[0] + q[1] * x[1] + lmd[0];
}

void nmpc_model::hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hu)
{
	hu[0] = r[0] * u[0] + lmd[1] * b * x[1];	
}


double nmpc_model::stagecost(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	return 0.5 * (q[0] * (x[0] - xf[0]) * (x[0] - xf[0]) + q[1] * (x[1] - xf[1]) * (x[1] - xf[1]) + r[0] * u[0] * u[0]);
}

double nmpc_model::terminalcost(const double t, const Eigen::VectorXd& x)
{
	return 0.5 * (sf[0] * (x[0] - xf[0]) * (x[0] - xf[0]) + sf[1] * (x[1] - xf[1]) * (x[1] - xf[1]));
}

void nmpc_model::Lufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> Ju)
{
	Ju[0] = r[0] * u[0];
}