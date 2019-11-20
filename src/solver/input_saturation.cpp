#include "input_saturation.hpp"

InputSaturation::InputSaturation() 
  : index_(0), 
    min_(0), 
    max_(0), 
    dummy_weight_(0), 
    quadratic_weight_(0) {
}

InputSaturation::InputSaturation(const int index, const double min, 
                                 const double max, const double dummy_weight, 
                                 const double quadratic_weight) 
  : index_(index), 
    min_(min), 
    max_(max), 
    dummy_weight_(dummy_weight), 
    quadratic_weight_(quadratic_weight) {
}

InputSaturation::InputSaturation(const int index, const double min, 
                                 const double max, const double dummy_weight)
  : index_(index), 
    min_(min), 
    max_(max), 
    dummy_weight_(dummy_weight), 
    quadratic_weight_(0) {
}

InputSaturation::InputSaturation(const InputSaturation& other) {
  setParameters(other.index(), other.min(), other.max(), other.dummy_weight(), 
                other.quadratic_weight());
}

InputSaturation& InputSaturation::operator=(const InputSaturation& other) {
  if (this!=&other) {
    setParameters(other.index(), other.min(), other.max(), other.dummy_weight(), 
                  other.quadratic_weight());
  }
  return *this;
}

bool InputSaturation::operator<(const InputSaturation& other) const {
  return index_ < other.index();
}

void InputSaturation::setParameters(const int index, const double min, 
                                    const double max, const double dummy_weight, 
                                    const double quadratic_weight) {
  index_ = index;
  min_ = min;
  max_ = max;
  dummy_weight_ = dummy_weight;
  quadratic_weight_ = quadratic_weight;
}

void InputSaturation::setParameters(const int index, const double min, 
                                    const double max, 
                                    const double dummy_weight) {
  index_ = index;
  min_ = min;
  max_ = max;
  dummy_weight_ = dummy_weight;
  quadratic_weight_ = 0;
}