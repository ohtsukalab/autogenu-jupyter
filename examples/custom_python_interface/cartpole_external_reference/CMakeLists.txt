cmake_minimum_required(VERSION 3.1)
project(cartpole_external_reference CXX)

set(CMAKE_CXX_STANDARD 17)

if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
option(VECTORIZE "Enable -march=native" ON)

find_package(cgmres REQUIRED)

add_subdirectory(python/common)
add_subdirectory(python/${PROJECT_NAME})
