#ifndef CONTROL_INPUT_SATURATION_SEQUENCE_H
#define CONTROL_INPUT_SATURATION_SEQUENCE_H

#include <vector>
#include "control_input_saturation.hpp"

// Stores parameters representing the saturation on the control input. 
// ControlInputSaturationSequence is composed of several 
// ControlInputSaturation, which is a class for representing a saturation
// on a control input .Member variables of ControlInputSaturation:
// index_: the index of the element of the control input that is constrained.
// min_: the minimum value of the constrained element of the control input.
// max_: the maximum value of the constrained element of the control input.
// dummy_weight_: the weight parameter on the corresponding dummy input in the 
//   cost function.
// quadratic_weiht_: the weight parameter on the quadratic term for the 
// corresponding dummy input in the cost function.
// Each member variables of i-th ControlInputSaturation can be obtained by 
// access functions, e.g., by index(i), min(i), max(i), dummy_weight(i), 
// and quadratic_weight(i).
class ControlInputSaturationSequence{
public:
  // Constructs ControlInputSaturationSequence with empty.
  ControlInputSaturationSequence();
  // Constructs ControlInputSaturationSequence with a ControlInputSaturation.
  ControlInputSaturationSequence(
      ControlInputSaturation& control_input_saturation);
  // Constructs ControlInputSaturationSequence with several 
  // ControlInputSaturations using std::vector.
  ControlInputSaturationSequence(
      std::vector<ControlInputSaturation>& control_input_saturation_seq);

  // Appends a ControlInputSaturation that has parameters in arguments to 
  // control_input_saturation_seq_. If there is the saturation that has the 
  // same index in arguments, overwrite its min_, max_, dummy_weight_, 
  // and quadratic_weight.
  void appendControlInputSaturation(const int index, const double min, 
                                    const double max, 
                                    const double dummy_weight, 
                                    const double quadratic_weight);
  // Appends a ControlInputSaturation that has parameters in arguments to 
  // control_input_saturation_seq_. Note that quadratic_weight_ is set by 0. 
  // If there is the saturation that has the same index in arguments, 
  // overwrite its min_, max_, dummy_weight_.
  void appendControlInputSaturation(const int index, const double min, 
                                    const double max, 
                                    const double dummy_weight);

  // Returns the number of the ControlInputSaturation, that is, the number 
  // of the elements of the control input.
  inline int dimSaturation() const {
    return control_input_saturation_seq_.size();
  }
  // Returns the index of the element of the control input that is constrained
  // by ControlInputSaturation whose stored index is saturation_index.
  inline int index(const int saturation_index) const {
    return control_input_saturation_seq_[saturation_index].index();
  }
  // Returns the minimum value of the constrained element of the control input
  // that is constrained by ControlInputSaturation whose stored index is 
  // saturation_index.
  inline double min(const int saturation_index) const {
    return control_input_saturation_seq_[saturation_index].min();
  }
  // Returns the maximum value of the constrained element of the control input
  // that is constrained by ControlInputSaturation whose stored index is 
  // saturation_index.
  inline double max(const int saturation_index) const {
    return control_input_saturation_seq_[saturation_index].max();
  }
  // Returns the weight parameter on the linear term of the dummy input 
  // corresponting ControlInpuSaturation whose stored index is saturatio_index.
  inline double dummy_weight(const int saturation_index) const {
    return control_input_saturation_seq_[saturation_index].dummy_weight();
  }
  // Returns the weight parameter on the quadratic term of the dummy input 
  // corresponting ControlInpuSaturation whose stored index is saturatio_index.
  inline double quadratic_weight(const int saturation_index) const {
    return control_input_saturation_seq_[saturation_index].quadratic_weight();
  }

private:
  std::vector<ControlInputSaturation> control_input_saturation_seq_;

  // Returns the number of the saturation that has the same index if it exists 
  // in control_input_saturation_seq_, and returns the total number of 
  // saturations otherwise.
  int findSameIndex(const int index);
};

#endif // CONTROL_INPUT_SATURATION_SEQUENCE_H