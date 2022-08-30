import re


class SimulationConditions(object):
    """ Parameters about the simulations. Provide these parameters to plot the
        graph of the simulation and draw the animation of the simulation.
    """

    def __init__(self, ocp_name: str):
        """ Inits SimulationConditions with the simulation results. """
        simulation_conditions = open(ocp_name + '_conditions' + '.log')
        lines = simulation_conditions.readlines()
        if re.search(r'.', lines[1]):
            pattern = r'([0-9]+\.?[0-9]*)' 
            self.simulation_time = float(re.findall(pattern, lines[1])[0])
        else:
            pattern = r'([0-9])' 
            self.simulation_time = float(re.findall(pattern, lines[1])[0])
        if re.search(r'.', lines[3]):
            pattern = r'([0-9]+\.?[0-9]*)' 
            self.sampling_period = float(re.findall(pattern, lines[3])[0]) * 0.0001
        else:
            pattern = r'([0-9])' 
            self.sampling_period = float(re.findall(pattern, lines[3])[0]) * 0.0001