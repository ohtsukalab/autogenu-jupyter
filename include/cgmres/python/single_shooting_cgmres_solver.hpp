#pragma once

namespace cgmres {
namespace python {

template <typename OCP, int N, int KMAX>
inline void bind_single_shooting_cgmres_solver(pybind11::module_& m) {
  using SingleShootingCGMRESSolver_ = SingleShootingCGMRESSolver<OCP, N, KMAX>;
  pybind11::class_<SingleShootingCGMRESSolver_>(m, "SingleShootingCGMRESSolver")
    .def(pybind11::init<OCP, Horizon, SolverSettings>(),
          pybind11::arg("ocp"), pybind11::arg("horizon"), pybind11::arg("settings"))
    .def(pybind11::init<>())
    .def("clone", [](const SingleShootingCGMRESSolver_& self) {
       auto copy = self;
       return copy;
     })
    .def("set_u", [](SingleShootingCGMRESSolver_& self, const VectorX& u) {
        self.set_u(u);
     }, pybind11::arg("u"))
    .def("set_uc", [](SingleShootingCGMRESSolver_& self, const VectorX& uc) {
        self.set_uc(uc);
     }, pybind11::arg("uc"))
    .def("set_u_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& u_array) {
        self.set_u_array(u_array);
     }, pybind11::arg("u_array"))
    .def("set_uc_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& uc_array) {
        self.set_uc_array(uc_array);
     }, pybind11::arg("uc_array"))
    .def("set_dummy_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& dummy_array) {
        self.set_dummy_array(dummy_array);
     }, pybind11::arg("dummy_array"))
    .def("set_mu_array", [](SingleShootingCGMRESSolver_& self, const std::vector<VectorX>& mu_array) {
        self.set_mu_array(mu_array);
     }, pybind11::arg("mu_array"))
    .def_property_readonly("uopt", &SingleShootingCGMRESSolver_::uopt)
    .def_property_readonly("ucopt", &SingleShootingCGMRESSolver_::ucopt)
    .def_property_readonly("xopt", &SingleShootingCGMRESSolver_::xopt)
    .def_property_readonly("lmdopt", &SingleShootingCGMRESSolver_::lmdopt)
    .def("opt_error", [](SingleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) {
        return self.optError(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("opt_error", static_cast<Scalar (SingleShootingCGMRESSolver_::*)() const>(&SingleShootingCGMRESSolver_::optError))
    .def("update", [](SingleShootingCGMRESSolver_& self, const Scalar t, const VectorX& x) {
        self.update(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("init_dummy_mu", &SingleShootingCGMRESSolver_::init_dummy_mu)
    .def("get_profile", &SingleShootingCGMRESSolver_::getProfile)
    .def("__str__", [](const SingleShootingCGMRESSolver_& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace cgmres
