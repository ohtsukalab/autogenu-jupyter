#ifndef CONTROL_INPUT_SATURATION_SEQUENCE_H
#define CONTROL_INPUT_SATURATION_SEQUENCE_H


#include <vector>
#include "control_input_saturation.hpp"

// Stores saturation of control input
class ControlInputSaturationSequence{
private:
    std::vector<ControlInputSaturation> control_input_saturation_seq_;

    int findSameIndex(const int index);


public:
    void appendControlInputSaturation(const int index, const double max, const double min, const double dummy_weight);
    int dimSaturation() const;
    int index(const int num_saturation) const;
    double max(const int num_saturation) const;
    double min(const int num_saturation) const;
    double weight(const int num_saturation) const;
};


#endif