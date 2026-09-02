#pragma once

namespace cgmres {
namespace python {

inline void bind_horizon(pybind11::module_& m) {
  pybind11::class_<Horizon>(m, "Horizon")
    .def(pybind11::init<const Scalar, const Scalar, const Scalar>(),
          pybind11::arg("Tf"), pybind11::arg("alpha")=0.0, pybind11::arg("t0")=0.0)
    .def(pybind11::init<>())
    .def("clone", [](const Horizon& self) {
       auto copy = self;
       return copy;
     })
    .def("T", &Horizon::T,
          pybind11::arg("t"))
    .def("__str__", [](const Horizon& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace cgmres
