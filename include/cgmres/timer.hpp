#ifndef CGMRES__TIMER_HPP_
#define CGMRES__TIMER_HPP_

#include <chrono>
#include <iostream>

#include "cgmres/types.hpp"


namespace cgmres {

///
/// @class TimingProfile  
/// @brief A profile of the timing benchmark. 
///
struct TimingProfile {
  ///
  /// @brief Average computational time in milliseconds.
  ///
  Scalar average_time_ms = 0;

  ///
  /// @brief Maximum computational time in milliseconds.
  ///
  Scalar max_time_ms = 0;

  ///
  /// @brief Number of timing counts.
  ///
  unsigned long counts = 0;

  void disp(std::ostream& os) const {
    os << "TimingProfile: " << std::endl; 
    os << "  average time: " << average_time_ms << " [ms]" << std::endl;
    os << "  max time:     " << max_time_ms << " [ms]" << std::endl;
    os << "  counts:       " << counts << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const TimingProfile& profile) {
    profile.disp(os);
    return os;
  }
};


///
/// @class Timer
/// @brief A timer for benchmarks. 
///
class Timer {
public:
  ///
  /// @brief Default constructor.
  ///
  Timer() { reset(); }

  ///
  /// @brief Default destructor.
  ///
  ~Timer() = default;

  void reset() {
    counts_ = 0;
    total_elapsed_time_ = 0;
    max_elapsed_time_ = 0;
  }

  ///
  /// @brief Start timer (tick).
  ///
  void tick() {
    time_point_ = std::chrono::high_resolution_clock::now();
  }

  ///
  /// @brief Stop the timer (tock).
  ///
  void tock() {
    const auto now = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<Scalar, std::milli> elapsed_time = now - time_point_;
    total_elapsed_time_ += elapsed_time.count();
    max_elapsed_time_ = std::max(max_elapsed_time_, elapsed_time.count());
    ++counts_;
  }

  ///
  /// @brief Get timing result as TimingProfile.
  /// @return Timing profile.
  ///
  TimingProfile getProfile() const {
    TimingProfile profile;
    if (counts_ > 0) {
      profile.average_time_ms = total_elapsed_time_ / counts_;
      profile.max_time_ms = max_elapsed_time_;
      profile.counts = counts_;
    }
    return profile;
  }

private:
  unsigned long counts_;
  Scalar total_elapsed_time_, max_elapsed_time_;
  std::chrono::high_resolution_clock::time_point time_point_;
};

} // namespace cgmres

#endif // CGMRES__TIMER_HPP_