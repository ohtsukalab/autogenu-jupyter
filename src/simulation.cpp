#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Core>


class simulation : public numerical_integration{
private:
    std::strings fname;
    double tsim, dt;
public:
    simulation( const double sim_time, const double sample_dt, const std::string file_name, );
    ~simulation();
    void fout()

    virtual void statefunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& x1) = 0;
};


simulation::simulation(double sim_time, const double sample_dt, const std::string file_name)
{
    tsim = sim_time;
    dt = sample_dt;

}

simulation::fout()

void simulation::impl()
{
    for(i=0; i<){
        x = x1;
        solver(t, x, u);
        adams(statefunc, t, x, u, x1);
        fout(t, x, u, F);
    }
}





