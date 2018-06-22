#include "nmpc_solver.hpp"


nmpc_solver::nmpc_solver(const solverid s, const double length_horizon, const int dvision_num, const double conv_radius, const double h_diff, const int k_max) : newton_gmres(length_horizon, dvision_num, conv_radius, h_diff, k_max)
{
    if(s != solverid::newton_gmres_single){
        std::cout << "errors in augments of nmpc_solver(): " << std::endl;
        std::exit(0);
    }
    solveridx = s;
}

solverid nmpc_solver::getsolver()
{
    return solveridx;
}