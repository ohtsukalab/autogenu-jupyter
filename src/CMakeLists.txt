add_library(nmpcsolver STATIC nmpc_solver.cpp newton_gmres.cpp matrixfree_gmres.cpp)
add_library(simulator STATIC simulator.cpp numerical_integrator.cpp)