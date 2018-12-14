#ifndef NMPC_MODEL_H
#define NMPC_MODEL_H

#define _USE_MATH_DEFINES

#include <eigen3/Eigen/Core>
#include <cmath>

class NMPCModel{
private:
    // declare parameters of the model here 
    static constexpr int dim_state_ = 4;
    static constexpr int dim_control_input_ = 2;
    static constexpr int dim_constraints_ = 0;

    // parametern in the model
    static constexpr double m1 = 0.2;
    static constexpr double m2 = 0.7;
    static constexpr double l1 = 0.3;
    static constexpr double l2 = 0.3;
    static constexpr double d1 = 0.15;
    static constexpr double d2 = 0.257;
    static constexpr double J1 = 0.006;
    static constexpr double J2 = 0.051;
    static constexpr double g = 9.80665;

    // parameters in the cost function 
	double q[dim_state_] = {1, 1, 0.1, 0.1};
	double r[dim_control_input_] = {0.1, 0.1};
	double q_terminal[dim_state_] = {1, 1, 0.1, 0.1};
	double x_ref[dim_state_] = {M_PI, 0.0, 0.0, 0.0};

public:
    int dimState() const;
    int dimControlInput() const;
    int dimConstraints() const;

    void stateFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f);
    void phixFunc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> lmd);
    void hxFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hx);
    void huFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hu);
};


#endif