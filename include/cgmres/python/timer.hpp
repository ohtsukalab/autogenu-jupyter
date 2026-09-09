#pragma once

namespace cgmres {
namespace python {

inline void bind_timer(pybind11::module_& m) {
  pybind11::class_<TimingProfile>(m, "TimingProfile")
    .def(pybind11::init<>())
    .def("clone", [](const TimingProfile& self) {
       auto copy = self;
       return copy;
     })
    .def_readwrite("average_time_ms", &TimingProfile::average_time_ms)
    .def_readwrite("max_time_ms", &TimingProfile::max_time_ms)
    .def_readwrite("counts", &TimingProfile::counts)
    .def("__str__", [](const TimingProfile& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
  pybind11::class_<Timer>(m, "Timer")
    .def(pybind11::init<>())
    .def("clone", [](const Timer& self) {
       auto copy = self;
       return copy;
     })
    .def("reset", &Timer::reset)
    .def("tick", &Timer::tick)
    .def("tock", &Timer::tock)
    .def("get_profile", &Timer::getProfile);
}

} // namespace python
} // namespace cgmres
