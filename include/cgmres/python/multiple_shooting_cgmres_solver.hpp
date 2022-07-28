#define DEFINE_PYBIND11_MODULE_MULTIPLE_SHOOTING_CGMRES_SOLVER(OCP, N, KMAX) \
using MultipleShootingCGMRESSolver_ = MultipleShootingCGMRESSolver<OCP, N, KMAX>; \
PYBIND11_MODULE(multiple_shooting_cgmres_solver, m) { \
  py::class_<MultipleShootingCGMRESSolver_>(m, "MultipleShootingCGMRESSolver") \
    .def(py::init<OCP, Horizon, SolverSettings>(), \ 
          py::arg("ocp"), py::arg("horizon"), py::arg("settings")) \
    .def(py::init<>()) \ 
    .def("clone", [](const MultipleShootingCGMRESSolver_& self) { \
       auto copy = self; \
       return copy; \
     }) \
    .def("set_u", [](MultipleShootingCGMRESSolver_& self, const VectorX& u) { \
        if (u.size() != OCP::nu) { \ 
          throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'u.size()' must be "+std::to_string(OCP::nu)); \ 
        } \ 
        self.set_u(u); \ 
     }, py::arg("u")) \
    .def("set_uc", [](MultipleShootingCGMRESSolver_& self, const VectorX& uc) { \
        if (uc.size() != OCP::nuc) { \ 
          throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'uc.size()' must be "+std::to_string(OCP::nuc)); \ 
        } \ 
        self.set_uc(uc); \ 
     }, py::arg("uc")) \
    .def("set_x", [](MultipleShootingCGMRESSolver_& self, const VectorX& x) { \
        if (x.size() != OCP::nx) { \ 
          throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \ 
        } \ 
        self.set_x(x); \ 
     }, py::arg("x")) \
    .def("set_lmd", [](MultipleShootingCGMRESSolver_& self, const VectorX& lmd) { \
        if (lmd.size() != OCP::nx) { \ 
          throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'lmd.size()' must be "+std::to_string(OCP::nx)); \ 
        } \ 
        self.set_lmd(lmd); \ 
     }, py::arg("lmd")) \
    .def("uopt", &MultipleShootingCGMRESSolver_::uopt) \
    .def("ucopt", &MultipleShootingCGMRESSolver_::ucopt) \
    .def("xopt", &MultipleShootingCGMRESSolver_::xopt) \
    .def("lmdopt", &MultipleShootingCGMRESSolver_::lmdopt) \
    .def("opt_error", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        return self.optError(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("opt_error", static_cast<Scalar (MultipleShootingCGMRESSolver_::*)() const>(&MultipleShootingCGMRESSolver_::optError)) \
    .def("update", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        self.update(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("init_x_lmd", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[MultipleShootingCGMRESSolver]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        self.init_x_lmd(t, x); \
    }, py::arg("t"), py::arg("x")) \
    .def("__str__", [](const MultipleShootingCGMRESSolver_& self) { \
        std::stringstream ss; \
        ss << self; \ 
        return ss.str(); \
      }); \
}