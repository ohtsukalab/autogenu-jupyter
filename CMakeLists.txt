cmake_minimum_required(VERSION 3.1)
project(nmpc_solver CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")

set(MODEL_DIR ${CMAKE_SOURCE_DIR}/examples/2linkfullactuated)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(src)
include_directories(${MODEL_DIR})

add_subdirectory(${MODEL_DIR})
add_subdirectory(src)
add_executable(a.out ${MODEL_DIR}/main.cpp)

target_link_libraries(a.out simulator)
target_link_libraries(a.out cgmres_solver) # The original C/GMRES
target_link_libraries(a.out multiple_shooting_cgmres) # The multiple shooting based C/GMRES without condensing of saturations on the control input
target_link_libraries(a.out multiple_shooting_cgmres_with_saturation) # The multiple shooting based C/GMRES with condensing of saturations on the control input
target_link_libraries(a.out model)
