#include <iostream>
#include <fstream>
#include <chrono>
#include <eigen3/Eigen/Core>
#include "numerical_integration.h"
#include "newton_gmres.h"
#include "continuation_gmres.h"


enum NMPC_solver{
    newton_gmres = 0;
    continuation_gmres = 1;
}


class simulator : public numerical_integration{
private:
    std::strings fname;
    double tsim, ht;
    NMPC_solver s;
public:
    simulation(NMPC_solver solver, model m);
    ~simulation();
    void fout(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u);
    void simlation(double sim_time, const double sample_ht, const std::string file_name);
    virtual void statefunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& x1) = 0;
};


simulator::simulation(NMPC_solver solver, model m)
{

    std::ofstream x_data(file_name + "x.dat");
    std::ofstream u_data(file_name + "u.dat");
    std::ofstream c_data(file_name + "c.dat");
}

simulator::~simulation()
{
    c_data << fname << "\n";
    c_data << "tsim = " << tsim << "\n";
    c_data << "ht = " << ht << "\n";

    x_data.close();
    u_data.close();
    c_data.close();
}

void simulator::fout(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
    x_data << x << "\n";
    u_data << u << "\n";
}

void simulator::simlation(double sim_time, const double sample_ht, const std::string file_name)
{
    int i, isim;
    double step_time, total_time, t;
    Eigen::VectorXd x(dimx), u(dimu);
    nmpc_model model;

    std::chrono::system_clock::time_point  start, end;

    switch(s) {
        case newton_gmres:
            newton_gmres solver(model, const int division_num, const double h_dir, const int newtonitr_max, const double rho, const int k_max, const double T);
            break;
        case continuation_gmres:
            continuation_gmres solver(model, const int division_num, const double h_dir, const int k_max, const double T_f, const double alpha, const double zeta);
            solver.init(;
    }


    isim = tsim/ht;

    std::cout << " Start" << std::endl;

    total_time = 0;
    for(t=0, i=0; i<isim; i++, t+= ht){
        start = std::chrono::system_clock::now();
        solver.solver(t, x, u);
        end = std::chrono::system_clock::now();
        step_time = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
        total_time += step_time;
        fout(t, x, u);

        adams(statefunc, t, x, u, x1);
        x = x1;
    }

    std::cout << " End" << std::endl;
    std::cout << "CPU time: " << total_time << " [sec]" << std::endl;
}





