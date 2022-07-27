#define DEFINE_PYBIND11_MODULE_SINGLE_SHOOTING_CGMRES_SOLVER(OCP, N, KMAX) \
using SingleShootingCGMRESSolver_ = SingleShootingCGMRESSolver<OCP, N, KMAX>; \
PYBIND11_MODULE(single_shooting_cgmres_solver, m) { \
  py::class_<SingleShootingCGMRESSolver_>(m, "SingleShootingCGMRESSolver") \
    .def(py::init<OCP, Horizon, SolverSettings>(), \ 
          py::arg("ocp"), py::arg("horizon"), py::arg("settings")) \
    .def(py::init<>()) \ 
    .def("set_u", [](SingleShootingCGMRESSolver_& self, const VectorX& u) { \
        if (u.size() != OCP::nu) { \ 
          throw std::invalid_argument("[SingleShootingCGMRESSolver]: 'u.size()' must be "+std::to_string(OCP::nu)); \ 
        } \ 
        self.set_u(u); \ 
     }, py::arg("u")) \
    .def("set_uc", [](SingleShootingCGMRESSolver_& self, const VectorX& uc) { \
        if (uc.size() != OCP::nuc) { \ 
          throw std::invalid_argument("[SingleShootingCGMRESSolver]: 'uc.size()' must be "+std::to_string(OCP::nuc)); \ 
        } \ 
        self.set_uc(uc); \ 
     }, py::arg("uc")) \
    .def("uopt", &SingleShootingCGMRESSolver_::uopt) \
    .def("ucopt", &SingleShootingCGMRESSolver_::ucopt) \
    .def("xopt", &SingleShootingCGMRESSolver_::xopt) \
    .def("lmdopt", &SingleShootingCGMRESSolver_::lmdopt) \
    .def("opt_error", [](SingleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[SingleShootingCGMRESSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        return self.optError(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("opt_error", static_cast<Scalar (SingleShootingCGMRESSolver_::*)() const>(&SingleShootingCGMRESSolver_::optError)) \
    .def("update", [](SingleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[ZeroHorizonOCPSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        self.update(t, x); \
    }, py::arg("t"), py::arg("x")); \
}