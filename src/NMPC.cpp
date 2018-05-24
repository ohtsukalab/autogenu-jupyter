#include <iostream>
#include <eigen3/Eigen/Core>

enum NMPC_solver
{
    
}

class NMPC : public newton_gmres{
private:
    int dim_x, dim_y, dv;
    double T_h, h, 
    void (*func());
public:
    NMPC(dimx, dimy, dv);
    NMPC~();
    void solver(const double t, const Eigen::VectorXd& x, Eigen::VectorXd& u);
};


NMPC::solver


