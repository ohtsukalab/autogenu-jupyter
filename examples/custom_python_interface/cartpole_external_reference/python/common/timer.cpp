
// This file was automatically generated by autogenu-jupyter (https://github.com/ohtsukalab/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/timer.hpp"
#include "cgmres/python/timer.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

DEFINE_PYBIND11_MODULE_TIMER()

} // namespace python
} // namespace cgmres
