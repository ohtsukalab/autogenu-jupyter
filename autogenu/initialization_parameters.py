class InitializationParameters(object):
    """ Parameters for initialization of the NMPC solvers based on the C/GMRES 
        method. Provide these parameters to generate main.cpp of the simulation.

        Attributes: 
            initial_guess_solution: The initial guess of the solution of the 
                optimal control problem for initialization of the solution of 
                NMPC. 
            newton_residual_torelance: The torelance residual of the solution of 
                the optimal control problem for the initialization of the 
                solution of NMPC. The Newton iteration terminates when the 
                optimality error is less than tol_res.
            max_newton_iteration: The maxmum number of Newton iteration for the 
                initialization of the solution of NMPC.
            initial_Lagrange_multiplier: A optional parameter for 
                MSCGMRESWithInputSaturation. This is a part of the initial 
                guess of the solution, the initial guess of the Lagrange
                multiplier with respect the constraints on the saturation 
                function of the control input.
    """

    def __init__(self, initial_guess_solution, newton_residual_torelance, 
                 max_newton_iteration, initial_Lagrange_multiplier=None):
        """ Inits InitializationParameters. """
        self.initial_guess_solution = initial_guess_solution
        self.newton_residual_torelance = newton_residual_torelance 
        self.max_newton_iteration = max_newton_iteration 
        self.initial_Lagrange_multiplier = initial_Lagrange_multiplier