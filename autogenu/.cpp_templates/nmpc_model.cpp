#include "nmpc_model.hpp"


namespace cgmres {

// Computes the state equation f(t, x, u).
// t : time parameter
// x : state vector
// u : control input vector
// f : the value of f(t, x, u)
void NMPCModel::stateFunc(const double t, const double* x, const double* u, double* f) {

}

// Computes the partial derivative of terminal cost with respect to state, i.e., 
// dphi/dx(t, x).
// t    : time parameter
// x    : state vector
// phix : the value of dphi/dx(t, x)
void NMPCModel::phixFunc(const double t, const double* x, double* phix) {

}

// Computes the partial derivative of the Hamiltonian with respect to state, 
// i.e., dH/dx(t, x, u, lmd).
// t   : time parameter
// x   : state vector
// u   : control input vector
// lmd : the Lagrange multiplier for the state equation
// hx  : the value of dH/dx(t, x, u, lmd)
void NMPCModel::hxFunc(const double t, const double* x, const double* u, const double* lmd, double* hx) {

}

// Computes the partial derivative of the Hamiltonian with respect to control 
// input and the constraints, dH/du(t, x, u, lmd).
// t   : time parameter
// x   : state vector
// u   : control input vector
// lmd : the Lagrange multiplier for the state equation
// hu  : the value of dH/du(t, x, u, lmd)
void NMPCModel::huFunc(const double t, const double* x, const double* u, const double* lmd, double* hu) {

}

// Returns the dimension of the state.
int NMPCModel::dim_state() const {
  return dim_state_;
}

// Returns the dimension of the contorl input.
int NMPCModel::dim_control_input() const {
  return dim_control_input_;
}

// Returns the dimension of the constraints.
int NMPCModel::dim_constraints() const {
  return dim_constraints_;
}

} // namespace cgmres