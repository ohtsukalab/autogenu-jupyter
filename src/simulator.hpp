#ifndef SIMULATOR_H
#define SIMULATOR_H


#include <iostream>
#include <fstream>
#include <chrono>
#include <eigen3/Eigen/Core>
#include "nmpc_model.hpp"
#include "numerical_integrator.hpp"
#include "continuation_gmres.hpp"


class simulator : virtual public numerical_integrator{
private:
    double tsim, ht;
    nmpc_model model;
    void savedata(std::ofstream& x_data, std::ofstream& u_data, std::ofstream& e_data, const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double err);
public:
    void simulation(continuation_gmres solver, const Eigen::VectorXd& x0, const double sim_time, const double sample_ht, const std::string file_name);
};

#endif