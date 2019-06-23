class SimulationParams(object):
    """ Parameters of the numerical simulation of NMPC. Provide these parameters
        to generate main.cpp for the simulation.
        Attributes: 
            initial_time: The initial time of the simulation.
            initial_state: The state of the controlled system at the beginning of
                           the simulation.
            simulation_time: The length of the simulation.
            sampling_time: The sampling time of the simulation.
    """

    def __init__(self, initial_time, initial_state, simulation_time, sampling_time):
        """ Inits SimulationParams. """
        self.initial_time = initial_time
        self.initial_state = initial_state
        self.simulation_time = simulation_time
        self.sampling_time = sampling_time