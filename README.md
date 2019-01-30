# Introduction
This project provides the continuation/GMRES method (C/GMRES method) based solvers for nonlinear model predictive control (NMPC) and an automatic code generator for NMPC, called AutoGenU.

The following solvers are provided: 
- The original C/GMRES method (single shooting)
- The multiple shooting based C/GMRES method
- The multiple shooting based C/GMRES method with condensing of variables with respect to the constraints on the saturation function on the control input


# Requirement
- C++11
- Cmake
- Eigen 3
- Python 3.6.x, Jupyter Notebook, SymPy (to generate nmpc_model.hpp, nmpc_model.cpp, main.cpp by AutoGenU.ipynb)
- Python 3.6.x, NumPy, seaborn (to plot simulation data on AutoGenU.ipynb)


# Usage
### AutoGenU
AutoGenU generates following source files under your setting state equation and cost function: 
- nmpc_model.hpp
- nmpc_model.cpp
- main.cpp
- CMakeLists.txt
You can also build source files for numerical simulation, execute numerical simulation, and plot or save simulation result on AutoGenU.


### C/GMRES based solvers
The C/GMRES based solvers in src/solver directory can be used independently of AutoGenU. First, you have to set the following files:  
- nmpc_model.hpp: write parameters in your model  
- nmpc_model.cpp: write equations of your model  
- main.cpp: write parameters of solvers 
AutoGenU.ipynb generates these source files automatically. After setting these files, set CMakeLists.txt and build source files. 


# Demo
### The multiple shooting based C/GMRES method for a pendubot
![pendubot_multiple_shooting_gif](https://raw.githubusercontent.com/wiki/mayataka/CGMRES/images/fps=20.gif)
![pendubot_multiple_shooting_png](https://raw.github.com/wiki/mayataka/CGMRES/images/pendubot_multiple_shooting.png)

# License
MIT

# References
1. [T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)](https://doi.org/10.1016/j.automatica.2003.11.005)
2. [C. T. Kelly, Iterative methods for linear and nonlinear equations, Frontiers in Apllied Mathematics, SIAM (1995)](https://doi.org/10.1137/1.9781611970944)
3. [Y. Shimizu, T. Ohtsuka, M. Diehl, A real‐time algorithm for nonlinear receding horizon control using multiple shooting and continuation/Krylov method, International Journal of Robust and Nonlinear Control, Vol. 19, No. 8, pp. 919-936 (2008)](https://doi.org/10.1002/rnc.1363)