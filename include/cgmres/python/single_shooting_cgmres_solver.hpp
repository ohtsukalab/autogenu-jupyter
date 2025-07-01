#define DEFINE_PYBIND11_MODULE_SINGLE_SHOOTING_CGMRES_SOLVER(OCP, N, KMAX) \
using SingleShootingCGMRESSolver_ = SingleShootingCGMRESSolver<OCP, N, KMAX>; \
PYBIND11_MODULE(single_shooting_cgmres_solver, m) { \
  py::class_<SingleShootingCGMRESSolver_>(m, "SingleShootingCGMRESSolver") \
    .def(py::init<OCP, Horizon, SolverSettings>(), \ 
          py::arg("ocp"), py::arg("horizon"), py::arg("settings")) \
    .def(py::init<>()) \ 
    .def("clone", [](const SingleShootingCGMRESSolver_& self) { \
       auto copy = self; \
       return copy; \
     }) \
    .def("set_u", [](SingleShootingCGMRESSolver_& self, const VectorX& u) { \
        self.set_u(u); \ 
     }, py::arg("u")) \
    .def("set_uc", [](SingleShootingCGMRESSolver_& self, const VectorX& uc) { \
        self.set_uc(uc); \ 
     }, py::arg("uc")) \
    .def("set_u_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& u_array) { \
        self.set_u_array(u_array); \ 
     }, py::arg("u_array")) \
    .def("set_uc_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& uc_array) { \
        self.set_uc_array(uc_array); \ 
     }, py::arg("uc_array")) \
    .def("set_dummy_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& dummy_array) { \
        self.set_dummy_array(dummy_array); \ 
     }, py::arg("dummy_array")) \
    .def("set_mu_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& mu_array) { \
        self.set_mu_array(mu_array); \ 
     }, py::arg("mu_array")) \
    .def_property_readonly("uopt", &SingleShootingCGMRESSolver_::uopt) \
    .def_property_readonly("ucopt", &SingleShootingCGMRESSolver_::ucopt) \
    .def_property_readonly("xopt", &SingleShootingCGMRESSolver_::xopt) \
    .def_property_readonly("lmdopt", &SingleShootingCGMRESSolver_::lmdopt) \
    .def("opt_error", [](SingleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) { \
        return self.optError(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("opt_error", static_cast<Scalar (SingleShootingCGMRESSolver_::*)() const>(&SingleShootingCGMRESSolver_::optError)) \
    .def("update", [](SingleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) { \
        self.update(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("init_dummy_mu", &SingleShootingCGMRESSolver_::init_dummy_mu) \
    .def("get_profile", &SingleShootingCGMRESSolver_::getProfile) \
    .def("__str__", [](const SingleShootingCGMRESSolver_& self) { \
        std::stringstream ss; \
        ss << self; \ 
        return ss.str(); \
      }); \
}