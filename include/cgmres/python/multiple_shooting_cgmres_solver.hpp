#pragma once

namespace cgmres {
namespace python {

template <typename OCP, int N, int KMAX>
inline void bind_multiple_shooting_cgmres_solver(pybind11::module_& m) {
  using MultipleShootingCGMRESSolver_ = MultipleShootingCGMRESSolver<OCP, N, KMAX>;
  pybind11::class_<MultipleShootingCGMRESSolver_>(m, "MultipleShootingCGMRESSolver")
    .def(pybind11::init<OCP, Horizon, SolverSettings>(),
          pybind11::arg("ocp"), pybind11::arg("horizon"), pybind11::arg("settings"))
    .def(pybind11::init<>())
    .def("clone", [](const MultipleShootingCGMRESSolver_& self) {
       auto copy = self;
       return copy;
     })
    .def("set_u", [](MultipleShootingCGMRESSolver_& self, const VectorX& u) {
        self.set_u(u);
     }, pybind11::arg("u"))
    .def("set_uc", [](MultipleShootingCGMRESSolver_& self, const VectorX& uc) {
        self.set_uc(uc);
     }, pybind11::arg("uc"))
    .def("set_x", [](MultipleShootingCGMRESSolver_& self, const VectorX& x) {
        self.set_x(x);
     }, pybind11::arg("x"))
    .def("set_lmd", [](MultipleShootingCGMRESSolver_& self, const VectorX& lmd) {
        self.set_lmd(lmd);
     }, pybind11::arg("lmd"))
    .def("set_u_array", [](MultipleShootingCGMRESSolver_& self, const std::vector<VectorX>& u_array) {
        self.set_u_array(u_array);
     }, pybind11::arg("u_array"))
    .def("set_uc_array", [](MultipleShootingCGMRESSolver_& self, const std::vector<VectorX>& uc_array) {
        self.set_uc_array(uc_array);
     }, pybind11::arg("uc_array"))
    .def("set_x_array", [](MultipleShootingCGMRESSolver_& self, const std::vector<VectorX>& x_array) {
        self.set_x_array(x_array);
     }, pybind11::arg("x_array"))
    .def("set_lmd_array", [](MultipleShootingCGMRESSolver_& self, const std::vector<VectorX>& lmd_array) {
        self.set_lmd_array(lmd_array);
     }, pybind11::arg("lmd_array"))
    .def("set_dummy_array", [](MultipleShootingCGMRESSolver_& self, const std::vector<VectorX>& dummy_array) {
        self.set_dummy_array(dummy_array);
     }, pybind11::arg("dummy_array"))
    .def("set_mu_array", [](MultipleShootingCGMRESSolver_& self, const std::vector<VectorX>& mu_array) {
        self.set_mu_array(mu_array);
     }, pybind11::arg("mu_array"))
    .def_property_readonly("uopt", &MultipleShootingCGMRESSolver_::uopt)
    .def_property_readonly("ucopt", &MultipleShootingCGMRESSolver_::ucopt)
    .def_property_readonly("xopt", &MultipleShootingCGMRESSolver_::xopt)
    .def_property_readonly("lmdopt", &MultipleShootingCGMRESSolver_::lmdopt)
    .def_property_readonly("dummyopt", &MultipleShootingCGMRESSolver_::dummyopt)
    .def_property_readonly("muopt", &MultipleShootingCGMRESSolver_::muopt)
    .def("opt_error", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) {
        return self.optError(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("opt_error", static_cast<Scalar (MultipleShootingCGMRESSolver_::*)() const>(&MultipleShootingCGMRESSolver_::optError))
    .def("update", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) {
        self.update(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("init_x", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) {
        self.init_x(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("init_lmd", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) {
        self.init_lmd(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("init_x_lmd", [](MultipleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) {
        self.init_x_lmd(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("init_dummy_mu", &MultipleShootingCGMRESSolver_::init_dummy_mu)
    .def("get_profile", &MultipleShootingCGMRESSolver_::getProfile)
    .def("__str__", [](const MultipleShootingCGMRESSolver_& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace cgmres
