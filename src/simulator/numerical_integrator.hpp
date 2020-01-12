#ifndef NUMERICAL_INTEGRATOR_H
#define NUMERICAL_INTEGRATOR_H

#include "nmpc_model.hpp"


namespace cgmres {

// Supports numerical integration of the state equation of the system described 
// in nmpc_model.hpp for numerical simnulations.
class NumericalIntegrator {
public:
  NumericalIntegrator();

  // Euler method for the state equation.
  void euler(const double current_time, const double* current_state_vec, 
             const double* control_input_vec, const double integration_length,
             double* integrated_state);

  // The four-step Runge-Kutta-Gill method for the state equation.
  void rungeKuttaGill(const double current_time, 
                      const double* current_state_vec, 
                      const double* control_input_vec, 
                      const double integration_length, 
                      double* integrated_state);

private:
  NMPCModel model_;
};

} // namespace cgmres


#endif // NUMERICAL_INTEGRATOR_H