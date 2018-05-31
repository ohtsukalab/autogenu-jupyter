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