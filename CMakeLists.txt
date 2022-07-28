cmake_minimum_required(VERSION 3.1)
project(cgmres CXX)

set(CMAKE_CXX_STANDARD 17)

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
