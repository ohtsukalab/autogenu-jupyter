#include "simulator.hpp"


void simulator::savedata(std::ofstream& x_data, std::ofstream& u_data, std::ofstream& e_data, const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double err)
{
    int i;

    for(i=0; i<model.dimx; i++)
        x_data << x(i) << " ";
    x_data << "\n";

    for(i=0; i<model.dimu; i++)
        u_data << u(i) << " ";
    u_data << "\n";

    e_data << err << "\n";
}


void simulator::simulation(continuation_gmres solver, const Eigen::VectorXd& x0, const double sim_time, const double sample_ht, const std::string file_name)
{
    int i, isim;
    double step_time, total_time, t, err;
    Eigen::VectorXd x(model.dimx), x1(model.dimx), u(model.dimu);
    std::chrono::system_clock::time_point start, end;

    std::ofstream x_data(file_name + "x.dat");
    std::ofstream u_data(file_name + "u.dat");
    std::ofstream e_data(file_name + "e.dat");
    std::ofstream c_data(file_name + "c.dat");

    x = x0;
    isim = sim_time/sample_ht;

    std::cout << " Start simulation" << std::endl;
    total_time = 0.0;
    for(t=0, i=0; i<isim; i++, t+= sample_ht){
        err = solver.geterr();
        savedata(x_data, u_data, e_data, t, x, u, err);
        runge_kutta_gill(t, x, u, sample_ht, x1);
        start = std::chrono::system_clock::now();
        solver.solvenmpc(t, sample_ht, x, u);
        end = std::chrono::system_clock::now();
        step_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        step_time = step_time * 0.001;
        total_time += step_time;
        x = x1;
    }

    std::cout << " End" << std::endl;
    std::cout << "CPU time: " << total_time << " [sec]" << std::endl;
    
    c_data << file_name << "\n";
    c_data << "simulation time: " << sim_time << " [sec]\n";
    c_data << "CPU time (total): " << total_time << " [sec]\n";
    c_data << "sampling time: " << sample_ht << " [sec]\n";
    c_data << "CPU time (1step): " << total_time/isim << " [sec]\n";

    x_data.close();
    u_data.close();
    e_data.close();
    c_data.close();
}