cmake_minimum_required(VERSION 3.15...4.2)
project(cgmres CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CGMRES_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(CGMRES_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

function(cgmres_enable_strict_warnings target scope)
  if (MSVC)
    # Generated OCP callbacks keep a stable signature even when a symbolic
    # expression does not use every argument.
    target_compile_options(${target} ${scope} /W4 /wd4100)
    if (CGMRES_WARNINGS_AS_ERRORS)
      target_compile_options(${target} ${scope} /WX)
    endif()
  else()
    target_compile_options(
      ${target} ${scope} -Wall -Wextra -Wpedantic -Wno-unused-parameter
    )
    if (CGMRES_WARNINGS_AS_ERRORS)
      target_compile_options(${target} ${scope} -Werror)
    endif()
  endif()
endfunction()

function(cgmres_enable_sanitizers target scope)
  if (NOT CGMRES_ENABLE_SANITIZERS)
    return()
  endif()
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
      ${target} ${scope} -fsanitize=address,undefined -fno-omit-frame-pointer
    )
    target_link_options(${target} ${scope} -fsanitize=address,undefined)
  else()
    message(FATAL_ERROR "CGMRES_ENABLE_SANITIZERS requires GCC or Clang")
  endif()
endfunction()

# installation directories
include(GNUInstallDirs)
set(CGMRES_INSTALL_CMAKE_DIR "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake" CACHE STRING "The installation cmake directory")

# the interface library
add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
target_include_directories(
  ${PROJECT_NAME} INTERFACE 
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
cgmres_enable_strict_warnings(${PROJECT_NAME} INTERFACE)
cgmres_enable_sanitizers(${PROJECT_NAME} INTERFACE)

# cmake configs
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  VERSION 0.1
  COMPATIBILITY AnyNewerVersion
)
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake.in
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  INSTALL_DESTINATION ${CGMRES_INSTALL_CMAKE_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  NO_SET_AND_CHECK_MACRO
)

# install files
install(
  TARGETS ${PROJECT_NAME} 
  EXPORT ${PROJECT_NAME}-targets 
)
install(
  FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
        ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  DESTINATION ${CGMRES_INSTALL_CMAKE_DIR}
)
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/include/cgmres
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
# export
export(EXPORT ${PROJECT_NAME}-targets
  FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake
  NAMESPACE cgmres::
)
install(EXPORT ${PROJECT_NAME}-targets
  NAMESPACE cgmres::
  DESTINATION ${CGMRES_INSTALL_CMAKE_DIR}
)
