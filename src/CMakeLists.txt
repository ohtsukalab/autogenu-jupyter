# The original C/GMRES
add_library(cgmres_solver STATIC matrixfree_gmres.cpp init_cgmres.cpp continuation_gmres.cpp)

# The multiple shooting based C/GMRES without condensing of saturations on the control input
add_library(multiple_shooting_cgmres STATIC matrixfree_gmres.cpp init_cgmres.cpp multiple_shooting_cgmres.cpp)

# The multiple shooting based C/GMRES with condensing of saturations on the control input
add_library(multiple_shooting_cgmres_with_saturation STATIC matrixfree_gmres.cpp control_input_saturation.cpp control_input_saturation_sequence.cpp init_cgmres_with_saturation.cpp multiple_shooting_cgmres_with_saturation.cpp)

# Simulator
add_library(simulator STATIC simulator.cpp numerical_integrator.cpp)