class SolverParameters(object):
    """ Parameters of the NMPC solvers based on the C/GMRES method. Provide 
        these parameters to generate main.cpp for the simulation.

        Attributes: 
            T_f, alpha: Parameter about the length of the horizon of NMPC.
                The length of the horzion at time t is given by 
                T_f * (1-exp(-alpha*t)).
            N: The number of the grid for the discretization
                of the horizon of NMPC.
            finite_difference_increment: The small positive value for forward 
                approximation used in the FD-GMRES. 
            zeta: A stabilization parameter of the C/GMRES method. It may work 
                well if you set as zeta=1/sampling_period.
            kmax: Maximam number of the iteration of the Krylov 
                subspace method for the linear problem. 
    """

    def __init__(self, T_f, alpha, N, finite_difference_increment, zeta, kmax):
        """ Inits SolverParameters. """
        self.T_f = T_f
        self.alpha = alpha
        self.N = N 
        self.finite_difference_increment = finite_difference_increment
        self.zeta = zeta
        self.kmax = kmax