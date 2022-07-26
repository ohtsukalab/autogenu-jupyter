#define DEFINE_PYBIND11_MODULE_SOLVER_SETTINGS() \
PYBIND11_MODULE(solver_settings, m) { \
  py::class_<SolverSettings>(m, "SolverSettings") \
    .def(py::init<>()) \
    .def_readwrite("max_iter", &SolverSettings::max_iter) \
    .def_readwrite("opt_error_tol", &SolverSettings::opt_error_tol) \
    .def_readwrite("finite_difference_epsilon", &SolverSettings::finite_difference_epsilon) \
    .def_readwrite("dt", &SolverSettings::dt) \
    .def_readwrite("zeta", &SolverSettings::zeta) \
    .def_readwrite("verbose_level", &SolverSettings::verbose_level); \
} 