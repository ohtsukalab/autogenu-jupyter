This directory contains files to run closed-loop simulations in C++:
- `ocp.hpp` : A definition of the OCP generated by `autogenu-jupyter`
- `main.cpp` : Executable of the closed-loop simuation.
- `CMakeLists.txt` : CMake script to find `cgmres` C++ library and build the executable.