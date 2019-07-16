#include "control_input_saturation_sequence.hpp"


ControlInputSaturationSequence::ControlInputSaturationSequence()
  : control_input_saturation_seq_() {
}

ControlInputSaturationSequence::ControlInputSaturationSequence(
    ControlInputSaturation& control_input_saturation)
  : control_input_saturation_seq_({control_input_saturation}) {
}

ControlInputSaturationSequence::ControlInputSaturationSequence(
    std::vector<ControlInputSaturation>& control_input_saturation_seq)
  : control_input_saturation_seq_(control_input_saturation_seq) {
}

void ControlInputSaturationSequence::appendControlInputSaturation(
    const int index, const double min, const double max, 
    const double dummy_weight, const double quadratic_weight) {
  if (findSameIndex(index) < control_input_saturation_seq_.size()) {
    control_input_saturation_seq_[findSameIndex(index)] = 
        ControlInputSaturation(index, min, max, dummy_weight, quadratic_weight);
  }
  else {
    control_input_saturation_seq_.push_back(
        ControlInputSaturation(index, min, max, dummy_weight, quadratic_weight));
  }
}

void ControlInputSaturationSequence::appendControlInputSaturation(
    const int index, const double min, const double max, 
    const double dummy_weight) {
  if (findSameIndex(index) < control_input_saturation_seq_.size()) {
    control_input_saturation_seq_[findSameIndex(index)] = 
        ControlInputSaturation(index, min, max, dummy_weight, 0);
  }
  else {
    control_input_saturation_seq_.push_back(
        ControlInputSaturation(index, min, max, dummy_weight, 0));
  }
}

int ControlInputSaturationSequence::findSameIndex(const int index) {
  for (int i=0; i<control_input_saturation_seq_.size(); i++) {
    if (index == control_input_saturation_seq_[i].index()) {
      return i;
    }
  }
  return control_input_saturation_seq_.size();
}