#ifndef INPUT_SATURATION_H
#define INPUT_SATURATION_H

// Stores parameters representing the saturation on the control input. 
// Member variables of InputSaturation:
// index_: the index of the element of the control input that is constrained.
// min_: the minimum value of the constrained element of the control input.
// max_: the maximum value of the constrained element of the control input.
// dummy_weight_: the weight parameter on the corresponding dummy input in the 
//   cost function.
// quadratic_weiht_: the weight parameter on the quadratic term for the 
// corresponding dummy input in the cost function.
// Each member variables can be obtained by access functions, e.g., by index(), 
// min(), max(), dummy_weight(), and quadratic_weight().
class InputSaturation {
public:
  // Constructs InputSaturation with setting member variables as
  // index_=0, min_=0, max_=0, dummy_weight_=0, and quadratic_weight=0.
  InputSaturation();
  // Constructs InputSaturation with setting member variables 
  // index_, min_, max_, dummy_weight_, and quadratic_weight with values in 
  // arguments.
  InputSaturation(const int index, const double min, const double max, 
                  const double dummy_weight, const double quadratic_weight);
  // Constructs InputSaturation with setting member variables 
  // index_, min_, max_, dummy_weight_, and with values in arguments and 
  // setting quadratic_weight with 0.
  InputSaturation(const int index, const double min, const double max, 
                  const double dummy_weight);

  // Copy constructer of InputSaturation. 
  // Copy all private member variables, index_, min_, max_, dummy_weight_, 
  // and quadratic_weight_.
  InputSaturation(const InputSaturation& other);

  // Copy all private member variables, index_, min_, max_, dummy_weight_, 
  // and quadratic_weight_.
  InputSaturation& operator=(const InputSaturation& other);

  // Operator override for sort.
  bool operator<(const InputSaturation& other) const;

  // Sets index_, min_, max_, dummy_weight_, and quadratic_weight with values
  // in arguments.
  void setParameters(const int index, const double min, const double max, 
                     const double dummy_weight, const double quadratic_weight);

  // Sets index_, min_, max_, and dummy_weight_, with values in arguments and
  // sets quadratic_weight with 0.
  void setParameters(const int index, const double min, const double max, 
                     const double dummy_weight);

  // Returns the index of the element of the control input that is constrained.
  inline int index() const {
    return index_;
  }
  // Returns the minimum value of the constrained element of the control input.
  inline double min() const {
    return min_;
  }
  // Returns the maximum value of the constrained element of the control input.
  inline double max() const {
    return max_;
  }
  // Returns the weight parameter on the corresponding dummy input in the cost
  // function.
  inline double dummy_weight() const {
    return dummy_weight_;
  }
  // Returns the weight parameter on the quadratic term for the corresponding 
  // dummy input in the cost function.
  inline double quadratic_weight() const {
    return quadratic_weight_;
  }

private:
  int index_;
  double min_, max_, dummy_weight_, quadratic_weight_;
};

#endif // INPUT_SATURATION_H