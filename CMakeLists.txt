cmake_minimum_required(VERSION 3.1)
project(cgmres_simulator CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")


set(MODEL_DIR ${CMAKE_SOURCE_DIR}/models/hexacopter)
set(SOLVER_DIR ${CMAKE_SOURCE_DIR}/src/solver)
set(SIMULATOR_DIR ${CMAKE_SOURCE_DIR}/src/simulator)

include_directories(${MODEL_DIR})
include_directories(${SOLVER_DIR})
include_directories(${SIMULATOR_DIR})


add_subdirectory(${MODEL_DIR})

add_library(
    multiple_shooting_cgmres
    STATIC
    ${SOLVER_DIR}/multiple_shooting_cgmres.cpp
    ${SOLVER_DIR}/init_cgmres.cpp
    ${SOLVER_DIR}/matrixfree_gmres.cpp
    ${SOLVER_DIR}/linear_funcs.cpp
)
target_include_directories(
    multiple_shooting_cgmres
    PRIVATE
    ${MODEL_DIR}
    ${SOLVER_DIR}
)

add_library(
    multiple_shooting_cgmres_simulator
    STATIC
    ${SIMULATOR_DIR}/save_simulation_data.cpp
    ${SIMULATOR_DIR}/numerical_integrator.cpp
    ${SIMULATOR_DIR}/multiple_shooting_cgmres_simulator.cpp
)
target_include_directories(
    multiple_shooting_cgmres_simulator
    PRIVATE
    ${MODEL_DIR}
    ${SOLVER_DIR}
    ${SIMULATOR_DIR}
)
add_executable(a.out ${MODEL_DIR}/main.cpp)
target_link_libraries(a.out
    PRIVATE
    multiple_shooting_cgmres_simulator
    multiple_shooting_cgmres
    nmpcmodel
)
target_compile_options(
    a.out
    PRIVATE
    -O3
)
