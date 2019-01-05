#include "control_input_saturation_sequence.hpp"


int ControlInputSaturationSequence::findSameIndex(const ControlInputSaturation control_input_saturation)
{
    for(int i=0; i<control_input_saturation_seq_.size(); i++){
        if(control_input_saturation == control_input_saturation_seq_[i]){
            return i;
        }
    }
    return control_input_saturation_seq_.size()+1;
}


void ControlInputSaturationSequence::appendControlInputSaturation(const ControlInputSaturation control_input_saturation)
{
    if(findSameIndex(control_input_saturation) <= control_input_saturation_seq_.size()){
        control_input_saturation_seq_[findSameIndex(control_input_saturation)] = control_input_saturation;
    }
    else {
        control_input_saturation_seq_.push_back(control_input_saturation);
    }
}


void ControlInputSaturationSequence::sortIndex()
{
    std::sort(control_input_saturation_seq_.begin(), control_input_saturation_seq_.end());
}


int ControlInputSaturationSequence::dimSaturation() const
{
    return control_input_saturation_seq_.size();
}


int ControlInputSaturationSequence::index(const int num_saturation) const
{
    return control_input_saturation_seq_[num_saturation].index();
}


double ControlInputSaturationSequence::max(const int num_saturation) const
{
    return control_input_saturation_seq_[num_saturation].max();
}


double ControlInputSaturationSequence::min(const int num_saturation) const
{
    return control_input_saturation_seq_[num_saturation].min();
}


double ControlInputSaturationSequence::weight(const int num_saturation) const
{
    return control_input_saturation_seq_[num_saturation].weight();
}