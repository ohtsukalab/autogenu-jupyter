#ifndef NMPC_SOLVER_H
#define NMPC_SOLVER_H


#include "newton_gmres.hpp"

enum struct solverid : unsigned char
{
    newton_gmres_single,
    cgmres_single,
};



class nmpc_solver : public newton_gmres{
private:
    solverid solveridx;

public:
    nmpc_solver(const solverid s, const double length_horizon, const int dvision_num, const int itr_max, const double h_diff, const int k_max);
    solverid getsolver();
};

#endif