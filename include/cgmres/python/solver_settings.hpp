#pragma once

namespace cgmres {
namespace python {

inline void bind_solver_settings(pybind11::module_& m) {
  pybind11::class_<SolverSettings>(m, "SolverSettings")
    .def(pybind11::init<>())
    .def("clone", [](const SolverSettings& self) {
       auto copy = self;
       return copy;
     })
    .def_readwrite("max_iter", &SolverSettings::max_iter)
    .def_readwrite("opterr_tol", &SolverSettings::opterr_tol)
    .def_readwrite("finite_difference_epsilon", &SolverSettings::finite_difference_epsilon)
    .def_readwrite("sampling_time", &SolverSettings::sampling_time)
    .def_readwrite("zeta", &SolverSettings::zeta)
    .def_readwrite("min_dummy", &SolverSettings::min_dummy)
    .def_readwrite("verbose_level", &SolverSettings::verbose_level)
    .def("__str__", [](const SolverSettings& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace cgmres
