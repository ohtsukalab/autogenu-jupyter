class InitializationParams(object):
    """ Parameters for initialization of the NMPC solvers based on the C/GMRES 
        method. Provide these parameters to generate main.cpp of the simulation.

        Attributes: 
            initial_guess: The initial guess of the solution of the optimal 
                control problem for initialization of the solution of NMPC. 
            tol_res: The torelance residual of the solution of the optimal 
                control problem for the initialization of the solution of NMPC.
                The Newton iteration terminates when the optimality error is 
                less than tol_res.
            max_itr: The maxmum number of Newton iteration for the 
                initialization of the solution of NMPC.
            Lag_multiplier: A optional parameter for 
                MultipleShootingCGMRESWithSaturation. This is a part of the 
                initial guess of the solution, the initial guess of the Lagrange
                multiplier with respect the constraints on the saturation 
                function of the control input.
    """

    def __init__(self, initial_guess, tol_res, max_itr, Lag_multiplier=None):
        """ Inits InitializationParams. """
        self.initial_guess = initial_guess
        self.tol_res = tol_res
        self.max_itr = max_itr
        self.Lag_multiplier = Lag_multiplier