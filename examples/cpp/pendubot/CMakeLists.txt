
cmake_minimum_required(VERSION 3.1)
project(pendubot CXX)

set(CMAKE_CXX_STANDARD 17)

if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
option(VECTORIZE "Enable -march=native" ON)

find_package(cgmres REQUIRED)

add_executable(
  ${PROJECT_NAME}
  main.cpp
)
target_include_directories(
  ${PROJECT_NAME}
  PRIVATE
  ${CGMRES_INCLUDE_DIR}
)
if (VECTORIZE)
  target_compile_options(
    ${PROJECT_NAME}
    PRIVATE
    -march=native
  )
endif()