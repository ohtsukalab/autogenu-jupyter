
cmake_minimum_required(VERSION 3.15...4.2)
project(cartpole CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE
  cgmres::cgmres
)
if (VECTORIZE)
  target_compile_options(
    ${PROJECT_NAME}
    PRIVATE
    -march=native
  )
endif()
