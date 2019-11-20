// Abstruct class for optimal control problems. Loads model of NMPCand define 
// dimensions.

#ifndef OPTIMAL_CONTROL_PROBLEM_H
#define OPTIMAL_CONTROL_PROBLEM_H

#include "nmpc_model.hpp"

// Abstruct class for optimal control problems. This class loads model of NMPC
// and define dimensions.
class OptimalControlProblem {
public:
  // Loads model of NMPC and define dimensions.
  OptimalControlProblem();
  virtual ~OptimalControlProblem() = default;

  // Prohibits copy.
  OptimalControlProblem(const OptimalControlProblem&) = delete;
  OptimalControlProblem& operator=(const OptimalControlProblem&) = delete;
  
  // Returns dimension of the state.
  int dim_state() const;

  // Returns dimension of the control input.
  int dim_control_input() const;

  // Returns dimension of the constraints.
  int dim_constraints() const;

  // Returns dimension of the solution of the optimal control problem.
  virtual int dim_solution() const = 0;

protected:
  NMPCModel model_;
  int dim_state_, dim_control_input_, dim_constraints_, 
      dim_control_input_and_constraints_;
};

#endif // OPTIMAL_CONTROL_PROBLEM_H