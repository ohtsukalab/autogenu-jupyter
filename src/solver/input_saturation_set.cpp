#include "input_saturation_set.hpp"

InputSaturationSet::InputSaturationSet()
  : input_saturation_set_() {
}

InputSaturationSet::InputSaturationSet(const InputSaturation& input_saturation)
  : input_saturation_set_({input_saturation}) {
  sort();
}

InputSaturationSet::InputSaturationSet(
    const std::vector<InputSaturation>& input_saturation_set)
  : input_saturation_set_(input_saturation_set) {
  sort();
}

InputSaturationSet::InputSaturationSet(
    const InputSaturationSet& input_saturation_set)
  : input_saturation_set_() {
  if (input_saturation_set.dimSaturation() > 0) {
    for (int i=0; i<input_saturation_set.dimSaturation(); ++i) {
      appendInputSaturation(
          input_saturation_set.getInputSaturation(i));
    }
    sort();
  }
}

InputSaturationSet& InputSaturationSet::operator=(
    const InputSaturationSet& other) {
  if(this!=&other) {
    this->input_saturation_set_.clear();
    if (other.dimSaturation() > 0) {
      for (int i=0; i<other.dimSaturation(); ++i) {
        appendInputSaturation(other.getInputSaturation(i));
      }
    }
    sort();
  }
  return *this;
}

void InputSaturationSet::appendInputSaturation(
    const InputSaturation& input_saturation) {
  if (findSameIndex(input_saturation.index())<input_saturation_set_.size()) {
    input_saturation_set_[findSameIndex(input_saturation.index())] 
        = input_saturation;
  }
  else {
    input_saturation_set_.push_back(input_saturation);
  }
  sort();
}

void InputSaturationSet::appendInputSaturation(const int index, 
                                               const double min, 
                                               const double max, 
                                               const double dummy_weight, 
                                               const double quadratic_weight) {
  if (findSameIndex(index) < input_saturation_set_.size()) {
    input_saturation_set_[findSameIndex(index)] = 
        InputSaturation(index, min, max, dummy_weight, quadratic_weight);
  }
  else {
    input_saturation_set_.push_back(
        InputSaturation(index, min, max, dummy_weight, quadratic_weight));
  }
  sort();
}

void InputSaturationSet::appendInputSaturation(const int index, 
                                               const double min, 
                                               const double max, 
                                               const double dummy_weight) {
  if (findSameIndex(index) < input_saturation_set_.size()) {
    input_saturation_set_[findSameIndex(index)] = 
        InputSaturation(index, min, max, dummy_weight, 0.0);
  }
  else {
    input_saturation_set_.push_back(
        InputSaturation(index, min, max, dummy_weight, 0.0));
  }
  sort();
}

InputSaturation InputSaturationSet::getInputSaturation(const int num) const {
  return input_saturation_set_[num];
}

void InputSaturationSet::generateArray(int* index_array, double* min_array, 
                                       double* max_array, 
                                       double* dummy_weight_array, 
                                       double* quadratic_weight_array) {
  int size = input_saturation_set_.size();
  for (int i=0; i<size; ++i) {
    index_array[i] = index(i);
  }
  for (int i=0; i<size; ++i) {
    min_array[i] = min(i);
  }
  for (int i=0; i<size; ++i) {
    max_array[i] = max(i);
  }
  for (int i=0; i<size; ++i) {
    dummy_weight_array[i] = dummy_weight(i);
  }
  for (int i=0; i<size; ++i) {
    quadratic_weight_array[i] = quadratic_weight(i);
  }
}

void InputSaturationSet::print() {
  std::cout << "index : min : max : dummy_weight : quadratic_weight" 
      << std::endl;
  for (int i=0; i<input_saturation_set_.size(); ++i) {
    std::cout << index(i) << " ";
    std::cout << min(i) << " ";
    std::cout << max(i) << " ";
    std::cout << dummy_weight(i) << " ";
    std::cout << quadratic_weight(i) << std::endl;
  }
}

int InputSaturationSet::findSameIndex(const int index) {
  for (int i=0; i<input_saturation_set_.size(); i++) {
    if (index == input_saturation_set_[i].index()) {
      return i;
    }
  }
  return input_saturation_set_.size();
}

void InputSaturationSet::sort() {
  std::sort(input_saturation_set_.begin(), input_saturation_set_.end());
}