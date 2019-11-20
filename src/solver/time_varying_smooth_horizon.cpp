#include "time_varying_smooth_horizon.hpp"

TimeVaryingSmoothHorizon::TimeVaryingSmoothHorizon(const double T_f, 
                                                   const double alpha, 
                                                   const double initial_time) 
  : T_f_(T_f),
    alpha_(alpha),
    initial_time_(initial_time) {
}

TimeVaryingSmoothHorizon::TimeVaryingSmoothHorizon(const double T_f, 
                                                   const double alpha)
  : T_f_(T_f),
    alpha_(alpha),
    initial_time_(0.0) {
}

TimeVaryingSmoothHorizon::~TimeVaryingSmoothHorizon() {
}

void TimeVaryingSmoothHorizon::resetLength(const double initial_time) {
  initial_time_ = initial_time;
}

void TimeVaryingSmoothHorizon::resetLength(const double T_f, const double alpha, 
                                           const double initial_time) {
  T_f_ = T_f;
  alpha_ = alpha;
  initial_time_ = initial_time;
}