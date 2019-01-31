//
// Computes the initial solution of the multiple shooting based C/GMRES method with condensing of the saturations on the control input using the Newton GMRES method.
//

#ifndef INIT_CGMRES_WITH_SATURATION_H
#define INIT_CGMRES_WITH_SATURATION_H


#include <Eigen/Core>
#include "matrixfree_gmres.hpp"
#include "control_input_saturation_sequence.hpp"


// Computes the initial solution of the multiple shooting based C/GMRES method with condensing of the saturations on the control input using the Newton GMRES method.
class InitCGMRESWithSaturation final : public MatrixFreeGMRES{
private:
    NMPCModel model_;
    ControlInputSaturationSequence saturation_seq_;

    int dim_control_input_and_constraints_, dim_saturation_, dim_solution_;

    double difference_increment_;
    Eigen::VectorXd solution_update_vec_, lambda_vec_, error_vec_, error_vec_1_, error_vec_2_;


    // Adds partial derivative of the saturation with respect to the control input
    inline void addHamiltonianDerivativeWithControlInput(const Eigen::VectorXd& control_input_and_constraints_vec, const Eigen::VectorXd& saturation_lagrange_multiplier_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_control_input_and_constraints_vec);

    // Computes the optimality for the dummy input
    inline void computeDummyOptimality(const Eigen::VectorXd& dummy_input_vec, const Eigen::VectorXd& saturation_lagrange_multiplier_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_dummy);

    // Computes the optimality of the saturation
    inline void computeSaturationOptimality(const Eigen::VectorXd& control_input_and_constraint_vec, const Eigen::VectorXd& dummy_input_vec, Eigen::Ref<Eigen::VectorXd> optimality_for_saturation);


    // Computes the optimality error vector under current_solution_vec.
    inline void computeOptimalityErrors(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& solution_vec, Eigen::Ref<Eigen::VectorXd> optimality_vec);

    // Computes a vector correspongin to b in Ax=b
    void bFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, Eigen::Ref<Eigen::VectorXd> b_vec) override;

    // Generates a vector corresponding to Ax in Ax=b with using the forward difference approximation.
    void axFunc(const double time_param, const Eigen::VectorXd& state_vec, const Eigen::VectorXd& current_solution_vec, const Eigen::VectorXd& direction_vec, Eigen::Ref<Eigen::VectorXd> ax_vec) override;


public:
    // Sets parameters and allocates vectors.
    InitCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, const double difference_increment, const int max_dim_krylov);


    // Calls the forwardDifferenceGMRES, solves the GMRES, and obtains the solution 
    void solve0stepNOCP(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& initial_guess_control_input_vec, const Eigen::VectorXd& initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration, Eigen::Ref<Eigen::VectorXd> solution_vec);

    // Returns the optimality error for the control input and constraints under current_solution_vec.
    Eigen::VectorXd getControlInputAndConstraintsError(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec);

    // Returns the optimality error for dummy input under current_solution_vec.
    Eigen::VectorXd getDummyInputError(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec);

    // Returns the optimality error for saturation under current_solution_vec.
    Eigen::VectorXd getControlInputSaturationError(const double initial_time, const Eigen::VectorXd& initial_state_vec, const Eigen::VectorXd& current_solution_vec);
};


#endif