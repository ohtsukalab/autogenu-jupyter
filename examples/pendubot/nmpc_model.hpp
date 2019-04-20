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
    static constexpr int dim_state_ = 4;
    static constexpr int dim_control_input_ = 1;
    static constexpr int dim_constraints_ = 0;

    static constexpr double m1 = 0.2;
    static constexpr double m2 = 0.7;
    static constexpr double l1 = 0.3;
    static constexpr double l2 = 0.3;
    static constexpr double d1 = 0.15;
    static constexpr double d2 = 0.257;
    static constexpr double J1 = 0.006;
    static constexpr double J2 = 0.051;
    static constexpr double g = 9.80665;


    // Define parameters in the cost function here.
    double q[dim_state_] = {1, 1, 0.1, 0.1};
    double r[dim_control_input_] = {0.1};
    double q_terminal[dim_state_] = {1, 1, 0.1, 0.1};
    double x_ref[dim_state_] = {M_PI, 0.0, 0.0, 0.0};



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