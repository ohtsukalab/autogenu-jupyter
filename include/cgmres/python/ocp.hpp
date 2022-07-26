#define DEFINE_PYBIND11_MODULE_OCP(OCP) \
PYBIND11_MODULE(ocp, m) { \
  py::class_<OCP>(m, "OCP") \
    .def(py::init<>()) \
    .def("eval_f", [](const OCP& self, const Scalar t, \
                      const VectorX& x, const VectorX& u, VectorX& dx) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \
        if (u.size() != OCP::nu) { \ 
          throw std::invalid_argument("[OCP]: 'u.size()' must be "+std::to_string(OCP::nu)); \ 
        } \ 
        if (dx.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'dx.size()' must be "+std::to_string(OCP::nx)); \
        } \
        self.eval_f(t, x.data(), u.data(), dx.data()); \ 
     }, py::arg("t"), py::arg("x"), py::arg("u"), py::arg("dx")) \
    .def("eval_phix", [](const OCP& self, const Scalar t, \
                      const VectorX& x, VectorX& phix) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \ 
        if (phix.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'phix.size()' must be "+std::to_string(OCP::nx)); \
        } \
        self.eval_phix(t, x.data(), phix.data()); \
     }, py::arg("t"), py::arg("x"), py::arg("phix")) \
    .def("eval_hx", [](const OCP& self, const Scalar t, \ 
                       const VectorX& x, const VectorX& u, const VectorX& lmd, VectorX& hx) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \
        if (u.size() != OCP::nuc) { \
          throw std::invalid_argument("[OCP]: 'u.size()' must be "+std::to_string(OCP::nuc)); \
        } \
        if (lmd.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'lmd.size()' must be "+std::to_string(OCP::nx)); \
        } \
        if (hx.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'hx.size()' must be "+std::to_string(OCP::nx)); \
        } \
        self.eval_hx(t, x.data(), u.data(), lmd.data(), hx.data()); \
     }, py::arg("t"), py::arg("x"), py::arg("u"), py::arg("lmd"), py::arg("hx")) \
    .def("eval_hu", [](const OCP& self, const Scalar t, \
                       const VectorX& x, const VectorX& u, const VectorX& lmd, VectorX& hu) { \
        if (x.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'x.size()' must be "+std::to_string(OCP::nx)); \
        } \
        if (u.size() != OCP::nuc) { \
          throw std::invalid_argument("[OCP]: 'u.size()' must be "+std::to_string(OCP::nuc)); \
        } \
        if (lmd.size() != OCP::nx) { \
          throw std::invalid_argument("[OCP]: 'lmd.size()' must be "+std::to_string(OCP::nx));\ 
        } \
        if (hu.size() != OCP::nuc) { \
          throw std::invalid_argument("[OCP]: 'hx.size()' must be "+std::to_string(OCP::nuc)); \ 
        } \
        self.eval_hu(t, x.data(), u.data(), lmd.data(), hu.data()); \
     }, py::arg("t"), py::arg("x"), py::arg("u"), py::arg("lmd"), py::arg("hu")); \
}