//
// Stores parameters representing the saturation on the control input: the index of the element of the control input that is constrained, the minimum value of the constrained element of the control input, the maximum value of the constrained element of the control input, and the weight parameter on the corresponding dummy input in the cost function.
//

#ifndef CONTROL_INPUT_SATURATION_H
#define CONTROL_INPUT_SATURATION_H


// Stores parameters representing the saturation on the control input. 
class ControlInputSaturation{
private:
    // index_       : the index of the element of the control input that is constrained.
    // min_         : the minimum value of the constrained element of the control input.
    // max_         : the maximum value of the constrained element of the control input.
    // dummy_weight_: the weight parameter on the corresponding dummy input in the cost function.
    int index_;
    double min_;
    double max_;
    double dummy_weight_;
    double quadratic_weight_;


public:
    // Default constructor: sets as index=0, min=0, max=0, and dummy_weight=0.
    ControlInputSaturation();

    // Sets index_, min_, max_, dummy_weight_, quadratic_weight.
    ControlInputSaturation(const int index, const double min, const double max, const double dummy_weight, const double quadratic_weight);
    ControlInputSaturation(const int index, const double min, const double max, const double dummy_weight);

    // Sets index_, min_, max_, dummy_weight_, quadratic_weight.
    void setParams(const int index, const double min, const double max, const double dummy_weight, const double quadratic_weight);
    void setParams(const int index, const double min, const double max, const double dummy_weight);

    inline int index() const{
        return index_;
    }

    inline double min() const{
        return min_;
    }

    inline double max() const{
        return max_;
    }

    inline double dummy_weight() const{
        return dummy_weight_;
    }

    inline double quadratic_weight() const{
        return quadratic_weight_;
    }
};


#endif