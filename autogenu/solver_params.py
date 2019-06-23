class SolverParams(object):
    """ Parameters of the NMPC solvers based on the C/GMRES method. Provide 
        these parameters to generate main.cpp for the simulation.

        Attributes: 
            Tf: A parameter about the length of the horizon of NMPC.
            alpha: A parameter about the length of the horizon of NMPC.
            horizon_divs: The number of the grid for the discretization
                of the horizon of NMPC.
            finite_diff_step: The small positive value for forward 
                approximation used in the FD-GMRES. 
            zeta: A stabilization parameter of the C/GMRES method. It may work 
                well if you set as zeta=1/sampling_period.
            kmax: Maximam number of the iteration of the Krylov 
                subspace method for the linear problem. 
    """

    def __init__(self, T_f, alpha, horizon_divs, finite_diff_step, zeta, kmax):
        """ Inits SolverParams. """
        self.T_f = T_f
        self.alpha = alpha
        self.horizon_divs = horizon_divs
        self.finite_diff_step = finite_diff_step
        self.zeta = zeta
        self.kmax = kmax