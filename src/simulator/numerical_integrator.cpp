#include "numerical_integrator.hpp"


NumericalIntegrator::NumericalIntegrator() 
  : model_() {
}

void NumericalIntegrator::euler(const double current_time, 
                                const double* current_state_vec, 
                                const double* control_input_vec, 
                                const double integration_length, 
                                double* integrated_state) {
  double dx_vec_[model_.dimState()];
  model_.stateFunc(current_time, current_state_vec, control_input_vec, dx_vec_);
  for (int i=0; i<model_.dimState(); i++) {
    integrated_state[i] = current_state_vec[i] + integration_length*dx_vec_[i];
  }
}

void NumericalIntegrator::rungeKuttaGill(const double current_time, 
                                         const double* current_state_vec, 
                                         const double* control_input_vec, 
                                         const double integration_length, 
                                         double* integrated_state) {
  double k1_vec[model_.dimState()],  k2_vec[model_.dimState()], 
      k3_vec[model_.dimState()], k4_vec[model_.dimState()], 
      tmp_vec[model_.dimState()];

  model_.stateFunc(current_time, current_state_vec, control_input_vec, k1_vec);
  for (int i=0; i<model_.dimState(); i++) {
      tmp_vec[i] = current_state_vec[i] + 0.5*integration_length*k1_vec[i];
  }

  model_.stateFunc(current_time+0.5*integration_length, tmp_vec, 
                   control_input_vec, k2_vec);
  for (int i=0; i<model_.dimState(); i++) {
    tmp_vec[i] = current_state_vec[i] 
        + integration_length*0.5*(std::sqrt(2)-1)*k1_vec[i] 
        + integration_length*(1-(1/std::sqrt(2)))*k2_vec[i];
  }

  model_.stateFunc(current_time+0.5*integration_length, tmp_vec, 
                   control_input_vec, k3_vec);
  for (int i=0; i<model_.dimState(); i++) {
    tmp_vec[i] = current_state_vec[i] 
        - integration_length*0.5*std::sqrt(2)*k2_vec[i] 
        + integration_length*(1+(1/std::sqrt(2)))*k3_vec[i];
  }

  model_.stateFunc(current_time+integration_length, tmp_vec, 
                   control_input_vec, k4_vec);
  for (int i=0; i<model_.dimState(); i++) {
    integrated_state[i] = current_state_vec[i] 
        + (integration_length/6)
        * (k1_vec[i]+(2-std::sqrt(2))*k2_vec[i]
            +(2+std::sqrt(2))*k3_vec[i]+k4_vec[i]);
  }
}