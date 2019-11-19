#ifndef OPTIMAL_CONTROL_PROBLEM_H
#define OPTIMAL_CONTROL_PROBLEM_H

#include "nmpc_model.hpp"

class OptimalControlProblem {
public:
  OptimalControlProblem();
  virtual ~OptimalControlProblem() = default;

  OptimalControlProblem(const OptimalControlProblem&) = delete;
  OptimalControlProblem& operator=(const OptimalControlProblem&) = delete;
  
  int dim_state() const;
  int dim_control_input() const;
  int dim_constraints() const;

  virtual int dim_solution() const = 0;

protected:
  NMPCModel model_;
  int dim_state_, dim_control_input_, dim_constraints_, 
      dim_control_input_and_constraints_;

};

#endif // OPTIMAL_CONTROL_PROBLEM_H