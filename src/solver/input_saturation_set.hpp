#ifndef INPUT_SATURATION_SET_H
#define INPUT_SATURATION_SET_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "input_saturation.hpp"

// Stores parameters representing the saturation on the control input. 
// InputSaturationSet is composed of several 
// InputSaturation, which is a class for representing a saturation
// on a control input .Member variables of InputSaturation:
// index_: the index of the element of the control input that is constrained.
// min_: the minimum value of the constrained element of the control input.
// max_: the maximum value of the constrained element of the control input.
// dummy_weight_: the weight parameter on the corresponding dummy input in the 
//   cost function.
// quadratic_weiht_: the weight parameter on the quadratic term for the 
// corresponding dummy input in the cost function.
// Each member variables of i-th InputSaturation can be obtained by 
// access functions, e.g., by index(i), min(i), max(i), dummy_weight(i), 
// and quadratic_weight(i).
class InputSaturationSet {
public:
  // Constructs InputSaturationSet with empty.
  InputSaturationSet();

  // Constructs InputSaturationSet with a InputSaturation.
  InputSaturationSet(const InputSaturation& input_saturation);

  // Constructs InputSaturationSet with several 
  // InputSaturations using std::vector.
  InputSaturationSet(const std::vector<InputSaturation>& input_saturation_seq);

  // Copy constructs InputSaturationSet. Copies 
  // the private member data whose type is std::vector<InputSaturations>.
  InputSaturationSet(const InputSaturationSet& input_saturation_set);

  // Assignment operator copies the private member data whose type is 
  // std::vector<InputSaturations>.
  InputSaturationSet& operator=(const InputSaturationSet& other);

  // Appends a InputSaturation that has parameters in arguments to 
  // input_saturation_set_. Note that quadratic_weight_ is set by 0. 
  // If there is the saturation that has the same index in arguments, 
  // overwrite its min_, max_, dummy_weight_.
  void appendInputSaturation(const InputSaturation& input_saturation);

  // Appends a InputSaturation that has parameters in arguments to 
  // input_saturation_set_. If there is the saturation that has the 
  // same index in arguments, overwrite its min_, max_, dummy_weight_, 
  // and quadratic_weight.
  void appendInputSaturation(const int index, const double min, 
                             const double max, const double dummy_weight, 
                             const double quadratic_weight);

  // Appends a InputSaturation that has parameters in arguments to 
  // input_saturation_set_. Note that quadratic_weight_ is set by 0. 
  // If there is the saturation that has the same index in arguments, 
  // overwrite its min_, max_, dummy_weight_.
  void appendInputSaturation(const int index, const double min, 
                             const double max, const double dummy_weight);

  InputSaturation getInputSaturation(const int num) const;

  void generateArray(int* index_array, double* min_array, 
                     double* max_array, double* dummy_weight_array,
                     double* quadratic_weight_array);

  // Returns the number of the InputSaturation, that is, the number 
  // of the elements of the control input.
  inline int dimSaturation() const {
    return input_saturation_set_.size();
  }
  // Returns the index of the element of the control input that is constrained
  // by InputSaturation whose stored index is saturation_index.
  inline int index(const int saturation_index) const {
    return input_saturation_set_[saturation_index].index();
  }
  // Returns the minimum value of the constrained element of the control input
  // that is constrained by InputSaturation whose stored index is 
  // saturation_index.
  inline double min(const int saturation_index) const {
    return input_saturation_set_[saturation_index].min();
  }
  // Returns the maximum value of the constrained element of the control input
  // that is constrained by InputSaturation whose stored index is 
  // saturation_index.
  inline double max(const int saturation_index) const {
    return input_saturation_set_[saturation_index].max();
  }
  // Returns the weight parameter on the linear term of the dummy input 
  // corresponting ControlInpuSaturation whose stored index is saturatio_index.
  inline double dummy_weight(const int saturation_index) const {
    return input_saturation_set_[saturation_index].dummy_weight();
  }
  // Returns the weight parameter on the quadratic term of the dummy input 
  // corresponting ControlInpuSaturation whose stored index is saturatio_index.
  inline double quadratic_weight(const int saturation_index) const {
    return input_saturation_set_[saturation_index].quadratic_weight();
  }

  // Prints all elements of this object.
  void print();

private:
  std::vector<InputSaturation> input_saturation_set_;

  // Returns the number of the saturation that has the same index if it exists 
  // in input_saturation_set_, and returns the total number of 
  // saturations otherwise.
  int findSameIndex(const int index);

  void sort();
};

#endif // INPUT_SATURATION_SET_H