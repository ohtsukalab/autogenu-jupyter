#include "optimal_control_problem.hpp"

OptimalControlProblem::OptimalControlProblem()
  : model_(),
    dim_state_(model_.dimState()),
    dim_control_input_(model_.dimControlInput()),
    dim_constraints_(model_.dimConstraints()),
    dim_control_input_and_constraints_(
        model_.dimControlInput()+model_.dimConstraints()) {
}

int OptimalControlProblem::dim_state() const {
  return dim_state_;
}

int OptimalControlProblem::dim_control_input() const {
  return dim_control_input_;
}

int OptimalControlProblem::dim_constraints() const {
  return dim_constraints_;
}