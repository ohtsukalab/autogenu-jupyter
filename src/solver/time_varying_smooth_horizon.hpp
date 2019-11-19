#ifndef TIME_VARYING_SMOOTH_HORIZON_H
#define TIME_VARYING_SMOOTH_HORIZON_H

#include <cmath>

class TimeVaryingSmoothHorizon {
public:
  TimeVaryingSmoothHorizon(const double T_f, const double alpha, 
                           const double initial_time);

  TimeVaryingSmoothHorizon(const double T_f, const double alpha);

  ~TimeVaryingSmoothHorizon();

  inline double getLength(const double time) const {
    return T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
  }

  void resetLength(const double initial_time);
  void resetLength(const double T_f, const double alpha, 
                   const double initial_time);

  TimeVaryingSmoothHorizon(const TimeVaryingSmoothHorizon&) = delete;
  TimeVaryingSmoothHorizon& operator=(const TimeVaryingSmoothHorizon&) = delete;

private:
  double T_f_, alpha_, initial_time_;

};

#endif // TIME_VARYING_SMOOTH_HORIZON_H