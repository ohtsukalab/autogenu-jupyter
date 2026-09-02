#pragma once

namespace cgmres {
namespace python {

template <typename OCP, int KMAX>
inline void bind_zero_horizon_ocp_solver(pybind11::module_& m) {
  using ZeroHorizonOCPSolver_ = ZeroHorizonOCPSolver<OCP, KMAX>;
  pybind11::class_<ZeroHorizonOCPSolver_>(m, "ZeroHorizonOCPSolver")
    .def(pybind11::init<OCP, SolverSettings>(),
          pybind11::arg("ocp"), pybind11::arg("settings"))
    .def(pybind11::init<>())
    .def("clone", [](const ZeroHorizonOCPSolver_& self) {
       auto copy = self;
       return copy;
     })
    .def("set_u", [](ZeroHorizonOCPSolver_& self, const VectorX& u) {
        self.set_u(u);
     }, pybind11::arg("u"))
    .def("set_uc", [](ZeroHorizonOCPSolver_& self, const VectorX& uc) {
        self.set_uc(uc);
     }, pybind11::arg("uc"))
    .def("set_dummy", [](ZeroHorizonOCPSolver_& self, const VectorX& dummy) {
        self.set_dummy(dummy);
     }, pybind11::arg("dummy"))
    .def("set_mu", [](ZeroHorizonOCPSolver_& self, const VectorX& mu) {
        self.set_mu(mu);
     }, pybind11::arg("mu"))
    .def_property_readonly("uopt", &ZeroHorizonOCPSolver_::uopt)
    .def_property_readonly("ucopt", &ZeroHorizonOCPSolver_::ucopt)
    .def_property_readonly("lmdopt", &ZeroHorizonOCPSolver_::lmdopt)
    .def_property_readonly("dummyopt", &ZeroHorizonOCPSolver_::dummyopt)
    .def_property_readonly("muopt", &ZeroHorizonOCPSolver_::muopt)
    .def("opt_error", [](ZeroHorizonOCPSolver_& self, const Scalar t, const VectorX& x) {
        return self.optError(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("opt_error", static_cast<Scalar (ZeroHorizonOCPSolver_::*)() const>(&ZeroHorizonOCPSolver_::optError))
    .def("solve", [](ZeroHorizonOCPSolver_& self, const Scalar t, const VectorX& x) {
        self.solve(t, x);
    }, pybind11::arg("t"), pybind11::arg("x"))
    .def("get_profile", &ZeroHorizonOCPSolver_::getProfile)
    .def("__str__", [](const ZeroHorizonOCPSolver_& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace cgmres
