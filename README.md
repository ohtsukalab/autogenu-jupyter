# Introduction
This project provides Continuation/GMRES method (C/GMRES), a fast computation algorithm of NMPC, proposed by T. Ohtsuka. 

# Requirement
- C++11
- Cmake
- Eigen 3

# Usage
Firstly, write system's equation in nmpc_model.cpp and parameters in nmpc_model.hpp. In main.cpp, you have to select which solver to use, set parameters for solvers and declare solver(). After that, you set initial state x0 and implement simulation with setting simulation time, sampling time, and file name of simulation data.

# References
1. T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-–574 (2004)
2. C. T. Kelly, Iterative methods for linear and nonlinear equations, Frontiers in Apllied Mathematics, SIAM (1995)
