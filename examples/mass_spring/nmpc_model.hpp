#ifndef NMPC_MODEL_H
#define NMPC_MODEL_H


#define _USE_MATH_DEFINES

#include <eigen3/Eigen/Core>
#include <cmath>

class nmpc_model{
private:
	/* declare parameters of the model here */
	/* dimx: dimension of the state vector, dimu: dimension of the control input vector, dimc: dimenstion of constraints */
	static constexpr int dim_x = 2;
	static constexpr int dim_u = 1;
	static constexpr int dim_c = 0;

	/* parametern in the model */
	static constexpr double a = -1;
	static constexpr double b = -1;

	/* parameters in cost the function */
	double q[dim_x] = {10, 10};
	double r[dim_u] = {1};
	double sf[dim_x] = {10, 10};
	double xf[dim_x] = {0.0, 0.0};


public:
	int dimx, dimu, dimc;
	nmpc_model(){
		dimx = dim_x;
		dimu = dim_u;
		dimc = dim_c;
	}
	void phix(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> lmd);
	void statefunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f);
    void hxfunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hx);
	void hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hu);
	double stagecost(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u);
	double terminalcost(const double t, const Eigen::VectorXd& x);
	void Lufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> Ju);
};


#endif