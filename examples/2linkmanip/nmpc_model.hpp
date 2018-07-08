#ifndef NMPC_MODEL_H
#define NMPC_MODEL_H


#define _USE_MATH_DEFINES

#include <eigen3/Eigen/Core>
#include <cmath>

class nmpc_model{
private:
	/* declare parameters of the model here */
	/* dimx: dimension of the state vector, dimu: dimension of the control input vector, dimc: dimenstion of constraints */
	static constexpr int dim_x = 4;
	static constexpr int dim_u = 2;
	static constexpr int dim_c = 0;

	/* parametern in the model */
	static constexpr double m1 = 0.2;
	static constexpr double m2 = 0.7;
	static constexpr double l1 = 0.3;
	static constexpr double l2 = 0.3;
	static constexpr double d1 = 0.15;
	static constexpr double d2 = 0.257;
	static constexpr double J1 = 0.006;
	static constexpr double J2 = 0.051;
	static constexpr double g = 9.80665;

	/* parameters in cost the function */
	double q[dim_x] = {1, 1, 0.01, 0.01};
	double r[dim_u] = {0.01, 0.01};
	double sf[dim_x] = {1, 1, 0.01, 0.01};
	double xf[dim_x] = {M_PI, 0.0, 0.0, 0.0};


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
};


#endif