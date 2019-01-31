cmake_minimum_required(VERSION 3.1)
project(cgmres_simulator CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")

set(SOLVER_DIR ${CMAKE_SOURCE_DIR}/src/solver)
set(SIMULATOR_DIR ${CMAKE_SOURCE_DIR}/src/simulator)
set(MODEL_DIR ${CMAKE_SOURCE_DIR}/models/2linkwithsat)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${SOLVER_DIR})
include_directories(${SIMULATOR_DIR})
include_directories(${MODEL_DIR})


add_subdirectory(${MODEL_DIR})

add_library(
    multiple_shooting_cgmres_with_saturation
    STATIC
    ${SOLVER_DIR}/multiple_shooting_cgmres_with_saturation.cpp
    ${SOLVER_DIR}/init_cgmres_with_saturation.cpp
    ${SOLVER_DIR}/control_input_saturation.cpp
    ${SOLVER_DIR}/control_input_saturation_sequence.cpp
    ${SOLVER_DIR}/matrixfree_gmres.cpp
)

add_library(
    multiple_shooting_cgmres_with_saturation_simulator
    STATIC
    ${SIMULATOR_DIR}/multiple_shooting_cgmres_with_saturation_simulator.cpp
    ${SIMULATOR_DIR}/save_simulation_data.cpp
    ${SIMULATOR_DIR}/numerical_integrator.cpp
)


add_executable(a.out ${MODEL_DIR}/main.cpp)
target_link_libraries(a.out multiple_shooting_cgmres_with_saturation)
target_link_libraries(a.out multiple_shooting_cgmres_with_saturation_simulator)
target_link_libraries(a.out nmpcmodel)
