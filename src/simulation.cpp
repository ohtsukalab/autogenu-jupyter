#include <iostream>
#include <eigen3/Eigen/Core>

class simulation : public model {
private:
    std::strings fname;
    double tsim, dt;
public:
    simulation(double tsim, double dt);
    simulation();
    void adams();
    void euler();
};
