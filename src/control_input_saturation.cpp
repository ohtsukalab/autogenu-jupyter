#include "control_input_saturation.hpp"


ControlInputSaturation::ControlInputSaturation()
{
    index_ = 0;
    min_ = 0;
    max_ = 0;
    dummy_weight_ = 0;
}


ControlInputSaturation::ControlInputSaturation(const int index, const double min, const double max, const double dummy_weight)
{
    index_ = index;
    min_ = min;
    max_ = max;
    dummy_weight_ = dummy_weight;
}


void ControlInputSaturation::setParams(const int index, const double min, const double max, const double dummy_weight)
{
    index_ = index;
    min_ = min;
    max_ = max;
    dummy_weight_ = dummy_weight;
}


int ControlInputSaturation::index() const
{
    return index_;
}


double ControlInputSaturation::min() const
{
    return min_;
}


double ControlInputSaturation::max() const
{
    return max_;
}


double ControlInputSaturation::weight() const
{
    return dummy_weight_;
}
