import re


class SimulationConditions(object):
    """ Parameters about the simulations. Provide these parameters to plot the
        graph of the simulation and draw the animation of the simulation.
        Attributes: 
            simulation_time(): Returns the simulation time. 
            sampling_period(): Returns the sampling period.
    """

    def __init__(self, model_name):
        """ Inits SimulationConditions with the simulation results. """
        simulation_conditions = open(model_name + '_conditions' + '.dat')
        lines = simulation_conditions.readlines()
        if re.search(r'.', lines[1]):
            pattern = r'([0-9]+\.?[0-9]*)' 
            self.__simulation_time = float(re.findall(pattern, lines[1])[0])
        else:
            pattern = r'([0-9])' 
            self.__simulation_time = float(re.findall(pattern, lines[1])[0])
        if re.search(r'.', lines[3]):
            pattern = r'([0-9]+\.?[0-9]*)' 
            self.__sampling_period = float(re.findall(pattern, lines[3])[0])
        else:
            pattern = r'([0-9])' 
            self.__sampling_period = float(re.findall(pattern, lines[3])[0])

    def simulation_time(self):
        """ Returns the simulation time.
            Returns: 
                The time of the simulation.
        """
        return self.__simulation_time

    def sampling_period(self):
        """ Returns the sampling period of the simulation.
            Returns: 
                The sampling period of the simulation.
        """
        return self.__sampling_period