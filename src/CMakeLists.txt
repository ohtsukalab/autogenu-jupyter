add_library(nmpcsolver STATIC matrixfree_gmres.cpp newton_gmres.cpp continuation_gmres.cpp)
add_library(simulator STATIC simulator.cpp numerical_integrator.cpp)