//
// Stores the sequence of parameters representing the saturation on the control input: the index of the element of the control input that is constrained, the minimum value of the constrained element of the control input, the maximum value of the constrained element of the control input, and the weight parameter on the corresponding dummy input in the cost function.
//

#ifndef CONTROL_INPUT_SATURATION_SEQUENCE_H
#define CONTROL_INPUT_SATURATION_SEQUENCE_H


#include <vector>
#include "control_input_saturation.hpp"


// Stores sequence of parameters representing the saturation on the control input. 
class ControlInputSaturationSequence{
private:
    std::vector<ControlInputSaturation> control_input_saturation_seq_;

    // Returns the number of the saturation that has the same index if it exists in control_input_saturation_seq_, and returns the total number of saturations otherwise.
    int findSameIndex(const int index);


public:
    // Appends a saturation that has parameters in arguments to control_input_saturation_seq_. If there is the saturation that has the same index in arguments, overwrite its min_, max_, dummy_weight_, and quadratic_weight.
    void appendControlInputSaturation(const int index, const double min, const double max, const double dummy_weight, const double quadratic_weight);
    void appendControlInputSaturation(const int index, const double min, const double max, const double dummy_weight);


    // Returns the total number of saturations.
    int dimSaturation() const;


    inline int index(const int num_saturation) const{
        return control_input_saturation_seq_[num_saturation].index();
    }

    inline double min(const int num_saturation) const{
        return control_input_saturation_seq_[num_saturation].min();
    }

    inline double max(const int num_saturation) const{
        return control_input_saturation_seq_[num_saturation].max();
    }

    inline double dummy_weight(const int num_saturation) const{
        return control_input_saturation_seq_[num_saturation].dummy_weight();
    }

    inline double quadratic_weight(const int num_saturation) const{
        return control_input_saturation_seq_[num_saturation].quadratic_weight();
    }
};


#endif
