//
// Define parameters and equations of your model in this file.
//

#ifndef NMPC_MODEL_H
#define NMPC_MODEL_H

#define _USE_MATH_DEFINES

#include <cmath>


// This class stores parameters of NMPC and equations of NMPC.
class NMPCModel{
private:
    // Define parameters of your model here using "static constexpr".



public:
    // State equation of the model.
    void stateFunc(const double t, const double* x, const double* u, double* f);

    // Partial derivative of the terminal cost with respect to state.
    void phixFunc(const double t, const double* x, double* phix);

    // Partial derivative of the Hamiltonian with respect to state.
    void hxFunc(const double t, const double* x, const double* u, const double* lmd, double* hx);

    // Partial derivative of the Hamiltonian with respect to control input and constraints.
    void huFunc(const double t, const double* x, const double* u, const double* lmd, double* hu);



    int dimState() const{
        return dim_state_;
    }

    int dimControlInput() const{
        return dim_control_input_;
    }

    int dimConstraints() const{
        return dim_constraints_;
    }
};


#endif