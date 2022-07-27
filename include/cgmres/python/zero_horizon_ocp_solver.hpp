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
        if (u.size() != OCP::nu) { \ 
          throw std::invalid_argument("[ZeroHorizonOCPSolver]: 'u.size()' must be "+std::to_string(OCP::nu)); \ 
        } \ 
        self.set_u(u); \ 
     }, py::arg("u")) \
    .def("set_uc", [](ZeroHorizonOCPSolver_& self, const VectorX& uc) { \
        if (uc.size() != OCP::nuc) { \ 
          throw std::invalid_argument("[ZeroHorizonOCPSolver]: 'uc.size()' must be "+std::to_string(OCP::nuc)); \ 
        } \ 
        self.set_uc(uc); \ 
     }, py::arg("uc")) \
    .def("uopt", &ZeroHorizonOCPSolver_::uopt) \
    .def("ucopt", &ZeroHorizonOCPSolver_::ucopt) \
    .def("lmdopt", &ZeroHorizonOCPSolver_::lmdopt) \
    .def("opt_error", [](ZeroHorizonOCPSolver_& self, const Scalar t, const VectorX& x) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[ZeroHorizonOCPSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        return self.optError(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("opt_error", static_cast<Scalar (ZeroHorizonOCPSolver_::*)() const>(&ZeroHorizonOCPSolver_::optError)) \
    .def("solve", [](ZeroHorizonOCPSolver_& self, const Scalar t, const VectorX& x) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[ZeroHorizonOCPSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        self.solve(t, x); \
    }, py::arg("t"), py::arg("x")); \
}