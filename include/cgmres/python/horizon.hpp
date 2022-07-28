#define DEFINE_PYBIND11_MODULE_HORIZON() \
PYBIND11_MODULE(horizon, m) { \
  py::class_<Horizon>(m, "Horizon") \
    .def(py::init<const Scalar, const Scalar, const Scalar>(), \
          py::arg("Tf"), py::arg("alpha")=0.0, py::arg("t0")=0.0) \
    .def(py::init<>()) \ 
    .def("clone", [](const Horizon& self) { \
       auto copy = self; \
       return copy; \
     }) \
    .def("T", &Horizon::T, \
          py::arg("t")) \
    .def("__str__", [](const Horizon& self) { \
        std::stringstream ss; \
        ss << self; \ 
        return ss.str(); \
      }); \
}