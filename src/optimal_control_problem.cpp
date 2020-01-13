#include "optimal_control_problem.hpp"


namespace cgmres {

OptimalControlProblem::OptimalControlProblem()
  : model_(),
    dim_state_(model_.dim_state()),
    dim_control_input_(model_.dim_control_input()),
    dim_constraints_(model_.dim_constraints()),
    dim_control_input_and_constraints_(
        model_.dim_control_input()+model_.dim_constraints()) {
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

} // namespace cgmres