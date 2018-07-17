# Introduction
This project provides Continuation/GMRES method (C/GMRES), a fast computation algorithm of NMPC.

# Requirement
- C++11
- Cmake
- Eigen 3

# Usage
Write equations of nmpc in nmpc_model.cpp and parameters in nmpc_model.hpp. In main.cpp, you have to set parameters for C/GMRES and set initial state and initial guess solution. 

# References
1. T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-–574 (2004)
2. C. T. Kelly, Iterative methods for linear and nonlinear equations, Frontiers in Apllied Mathematics, SIAM (1995)
