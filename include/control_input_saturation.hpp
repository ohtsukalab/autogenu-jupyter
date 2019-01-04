//
// The multiple shooting based continuation GMRES (C/GMRES) method, a fast algorithm of nonlinear model predictive control (NMPC).
// This program is witten with reference to "T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)" and "Y. Shimizu, T. Ohtsuka, M. Diehl, A real‐time algorithm for nonlinear receding horizon control using multiple shooting and continuation/Krylov method, International Journal of Robust and Nonlinear Control, Vol. 19, No. 8, pp. 919-936 (2008)".
//

#ifndef CONTROL_INPUT_SATURATION_H
#define CONTROL_INPUT_SATURATION_H

#include <pick_model.hpp>
#include <eigen3/Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "init_cgmres.hpp"


struct Saturation{
public:
    int index;
    double max;
    double min;
};


// Stores saturation of control input
class ControlInputSaturation{
private:
    std::vector<Saturation> control_input_saturation_;

public:
    ControlInputSaturation();
    void appendSaturation(const Saturation saturation);
};


#endif