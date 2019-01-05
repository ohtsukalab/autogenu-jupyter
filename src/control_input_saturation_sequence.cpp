#include "control_input_saturation_sequence.hpp"



void ControlInputSaturationSequence::appendControlInputSaturation(const int index, const double max, const double min, const double dummy_weight)
{
    if(index <= control_input_saturation_seq_.size()){
        control_input_saturation_seq_[index] = ControlInputSaturation(index, max, min, dummy_weight);
    }
    else {
        control_input_saturation_seq_.push_back(ControlInputSaturation(index, max, min, dummy_weight));
    }
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