#define DEFINE_PYBIND11_MODULE_TIMER() \
PYBIND11_MODULE(timer, m) { \
  py::class_<TimingProfile>(m, "TimingProfile") \
    .def(py::init<>()) \ 
    .def("clone", [](const TimingProfile& self) { \
       auto copy = self; \
       return copy; \
     }) \
    .def_readwrite("average_time_ms", &TimingProfile::average_time_ms) \
    .def_readwrite("max_time_ms", &TimingProfile::max_time_ms) \
    .def_readwrite("counts", &TimingProfile::counts) \
    .def("__str__", [](const TimingProfile& self) { \
        std::stringstream ss; \
        ss << self; \ 
        return ss.str(); \
      }); \
  py::class_<Timer>(m, "Timer") \
    .def(py::init<>()) \ 
    .def("clone", [](const Timer& self) { \
       auto copy = self; \
       return copy; \
     }) \
    .def("reset", &Timer::reset) \
    .def("tick", &Timer::tick) \
    .def("tock", &Timer::tock) \
    .def("get_profile", &Timer::getProfile); \
}