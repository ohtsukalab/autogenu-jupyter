//
// Define parameters and equations of your model in this file.
//

#ifndef NMPC_MODEL_H
#define NMPC_MODEL_H

#define _USE_MATH_DEFINES

#include <Eigen/Core>
#include <cmath>


// This class stores parameters of NMPC and equations of NMPC.
class NMPCModel{
private:
    // Define parameters of your model here using "static constexpr".



public:
    // State equation of the model.
    void stateFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f);

    // Partial derivative of the terminal cost with respect to state.
    void phixFunc(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> phix);

    // Partial derivative of the Hamiltonian with respect to state.
    void hxFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hx);

    // Partial derivative of the Hamiltonian with respect to control input and constraints.
    void huFunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hu);



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