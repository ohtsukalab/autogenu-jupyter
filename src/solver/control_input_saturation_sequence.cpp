#include "control_input_saturation_sequence.hpp"



int ControlInputSaturationSequence::findSameIndex(const int index)
{
    for(int i=0; i<control_input_saturation_seq_.size(); i++){
        if(index == control_input_saturation_seq_[i].index()){
            return i;
        }
    }
    return control_input_saturation_seq_.size();
}


void ControlInputSaturationSequence::appendControlInputSaturation(const int index, const double min, const double max, const double dummy_weight, const double quadratic_weight)
{
    if(findSameIndex(index) < control_input_saturation_seq_.size()){
        control_input_saturation_seq_[findSameIndex(index)] = ControlInputSaturation(index, max, min, dummy_weight, quadratic_weight);
    }
    else {
        control_input_saturation_seq_.push_back(ControlInputSaturation(index, max, min, dummy_weight, quadratic_weight));
    }
}


void ControlInputSaturationSequence::appendControlInputSaturation(const int index, const double min, const double max, const double dummy_weight)
{
    if(findSameIndex(index) < control_input_saturation_seq_.size()){
        control_input_saturation_seq_[findSameIndex(index)] = ControlInputSaturation(index, max, min, dummy_weight, 0);
    }
    else {
        control_input_saturation_seq_.push_back(ControlInputSaturation(index, max, min, dummy_weight, 0));
    }
}



int ControlInputSaturationSequence::dimSaturation() const
{
    return control_input_saturation_seq_.size();
}