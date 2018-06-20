Numerical solvers for nonlinear model predictive control(NMPC)
# Introduction
This project provides numerical solvers for nonlinear model predictive control(NMPC). This project supports
- Newton GMRES method (single-shoooting)


# Requirement
- C++11
- Cmake
- Eigen 3

# Usage
Firstly, write system's equation in nmpc_model.cpp and parameters in nmpc_model.hpp. In main.cpp, you have to select which solver to use, set parameters for solvers and declare solver(). After that, you set initial state x0 and implement simulation with setting simulation time, sampling time, and file name of simulation data.

# References
1. C. T. Kelly, Iterative methods for linear and nonlinear equations, Frontiers in Apllied Mathematics, SIAM (1995)
