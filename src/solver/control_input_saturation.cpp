#include "control_input_saturation.hpp"


ControlInputSaturation::ControlInputSaturation() :
    index_(0), 
    min_(0), 
    max_(0), 
    dummy_weight_(0), 
    quadratic_weight_(0)
{}


ControlInputSaturation::ControlInputSaturation(const int index, const double min, const double max, const double dummy_weight, const double quadratic_weight) : 
    index_(index), 
    min_(min), 
    max_(max), 
    dummy_weight_(dummy_weight), 
    quadratic_weight_(quadratic_weight)
{}


ControlInputSaturation::ControlInputSaturation(const int index, const double min, const double max, const double dummy_weight) : 
    index_(index), 
    min_(min), 
    max_(max), 
    dummy_weight_(dummy_weight), 
    quadratic_weight_(0)
{}


void ControlInputSaturation::setParams(const int index, const double min, const double max, const double dummy_weight, const double quadratic_weight)
{
    index_ = index;
    min_ = min;
    max_ = max;
    dummy_weight_ = dummy_weight;
    quadratic_weight_ = quadratic_weight;
}


void ControlInputSaturation::setParams(const int index, const double min, const double max, const double dummy_weight)
{
    index_ = index;
    min_ = min;
    max_ = max;
    dummy_weight_ = dummy_weight;
    quadratic_weight_ = 0;
}
