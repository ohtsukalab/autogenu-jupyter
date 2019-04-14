cmake_minimum_required(VERSION 3.1)
project(nmpc_solver CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")

set(MODEL_DIR ${CMAKE_SOURCE_DIR}/examples/pendubot)

set(SOLVER_DIR ${CMAKE_SOURCE_DIR}/src/solver)
set(SIMULATOR_DIR ${CMAKE_SOURCE_DIR}/src/simulator)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${SOLVER_DIR})
include_directories(${SIMULATOR_DIR})
include_directories(${MODEL_DIR})

add_subdirectory(${MODEL_DIR})

add_library(
    cgmres 
    STATIC 
    ${SOLVER_DIR}/linear_funcs.cpp
    ${SOLVER_DIR}/multiple_shooting_cgmres.cpp
    ${SOLVER_DIR}/init_cgmres.cpp
    ${SOLVER_DIR}/matrixfree_gmres.cpp
)

add_library(
    cgmres_simulator 
    STATIC 
    ${SIMULATOR_DIR}/multiple_shooting_cgmres_simulator.cpp
    ${SIMULATOR_DIR}/save_simulation_data.cpp
    ${SIMULATOR_DIR}/numerical_integrator.cpp
)


add_executable(a.out ${MODEL_DIR}/main.cpp)

target_link_libraries(a.out model)
target_link_libraries(a.out cgmres)
target_link_libraries(a.out cgmres_simulator)