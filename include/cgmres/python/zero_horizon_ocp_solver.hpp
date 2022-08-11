#define DEFINE_PYBIND11_MODULE_ZERO_HORIZON_OCP_SOLVER(OCP, KMAX) \
using ZeroHorizonOCPSolver_ = ZeroHorizonOCPSolver<OCP, KMAX>; \
PYBIND11_MODULE(zero_horizon_ocp_solver, m) { \
  py::class_<ZeroHorizonOCPSolver_>(m, "ZeroHorizonOCPSolver") \
    .def(py::init<OCP, SolverSettings>(), \ 
          py::arg("ocp"), py::arg("settings")) \
    .def(py::init<>()) \ 
    .def("clone", [](const ZeroHorizonOCPSolver_& self) { \
       auto copy = self; \
       return copy; \
     }) \
    .def("set_u", [](ZeroHorizonOCPSolver_& self, const VectorX& u) { \
        self.set_u(u); \ 
     }, py::arg("u")) \
    .def("set_uc", [](ZeroHorizonOCPSolver_& self, const VectorX& uc) { \
        self.set_uc(uc); \ 
     }, py::arg("uc")) \
    .def("set_dummy", [](ZeroHorizonOCPSolver_& self, const VectorX& dummy) { \
        self.set_dummy(dummy); \ 
     }, py::arg("dummy")) \
    .def("set_mu", [](ZeroHorizonOCPSolver_& self, const VectorX& mu) { \
        self.set_mu(mu); \ 
     }, py::arg("mu")) \
    .def_property_readonly("uopt", &ZeroHorizonOCPSolver_::uopt) \
    .def_property_readonly("ucopt", &ZeroHorizonOCPSolver_::ucopt) \
    .def_property_readonly("lmdopt", &ZeroHorizonOCPSolver_::lmdopt) \
    .def_property_readonly("dummyopt", &ZeroHorizonOCPSolver_::dummyopt) \
    .def_property_readonly("muopt", &ZeroHorizonOCPSolver_::muopt) \
    .def("opt_error", [](ZeroHorizonOCPSolver_& self, const Scalar t, const VectorX& x) { \
        return self.optError(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("opt_error", static_cast<Scalar (ZeroHorizonOCPSolver_::*)() const>(&ZeroHorizonOCPSolver_::optError)) \
    .def("solve", [](ZeroHorizonOCPSolver_& self, const Scalar t, const VectorX& x) { \
        self.solve(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("get_profile", &ZeroHorizonOCPSolver_::getProfile) \
    .def("__str__", [](const ZeroHorizonOCPSolver_& self) { \
        std::stringstream ss; \
        ss << self; \ 
        return ss.str(); \
      }); \
}