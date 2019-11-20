// This class provides time varying smooth horizon. The length of the horizon 
// is given by T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_))).
#ifndef TIME_VARYING_SMOOTH_HORIZON_H
#define TIME_VARYING_SMOOTH_HORIZON_H

#include <cmath>

// Provides time varying smooth horizon. The length of the horizon is given by
// T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_))).
class TimeVaryingSmoothHorizon {
public:
  // Sets the parameters of the horizon.
  TimeVaryingSmoothHorizon(const double T_f, const double alpha, 
                           const double initial_time);

  // Sets the parameters of the horizon. The initial_time is set by zeto, e.g., 
  // the length of the horizon is set by T_f_ * (1.0-std::exp(-alpha_*time)).
  TimeVaryingSmoothHorizon(const double T_f, const double alpha);

  ~TimeVaryingSmoothHorizon();

  inline double getLength(const double time) const {
    return T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
  }

  // Resets the parameters of the horizon.
  void resetLength(const double initial_time);

  // Resets the parameters of the horizon.
  void resetLength(const double T_f, const double alpha, 
                   const double initial_time);

  // Prohibits copy.
  TimeVaryingSmoothHorizon(const TimeVaryingSmoothHorizon&) = delete;
  TimeVaryingSmoothHorizon& operator=(const TimeVaryingSmoothHorizon&) = delete;

private:
  double T_f_, alpha_, initial_time_;
};

#endif // TIME_VARYING_SMOOTH_HORIZON_H