#include "control_input_saturation.hpp"


ControlInputSaturation::ControlInputSaturation()
{
    index_ = 0;
    max_ = 0;
    min_ = 0;
    weight_ = 0;
}


ControlInputSaturation::ControlInputSaturation(const int index, const double max, const double min, const double weight)
{
    index_ = index;
    max_ = max;
    min_ = min;
    weight_ = weight;
}


void ControlInputSaturation::setParams(const int index, const double max, const double min, const double weight)
{
    index_ = index;
    max_ = max;
    min_ = min;
    weight_ = weight;
}


int ControlInputSaturation::index() const
{
    return index_;
}


double ControlInputSaturation::max() const
{
    return max_;
}


double ControlInputSaturation::min() const
{
    return min_;
}


double ControlInputSaturation::weight() const
{
    return weight_;
}
