cmake_minimum_required(VERSION 2.8)
project(nmpc_solver CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")

include_directories(${CMAKE_SOURCE_DIR})
include_directories(src)

add_subdirectory(src)
add_executable(a.out main.cpp nmpc_model.cpp)

target_link_libraries(a.out nmpcsolver)
target_link_libraries(a.out simulator)