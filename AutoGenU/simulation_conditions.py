import re


class SimulationConditions:
    def __init__(self, file_header):
        # load simulation conditions
        simulation_conditions = open(file_header + '_conditions' + '.dat')
        lines = simulation_conditions.readlines()

        # load float simlation time
        if re.search(r'.',lines[1]):
            pattern = r'([0-9]+\.?[0-9]*)' 
            self.__simulation_time = float(re.findall(pattern, lines
            [1])[0])
        # load int simulation time
        else:
            pattern = r'([0-9])' 
            self.__simulation_time = float(re.findall(pattern, lines[1])[0])

        # load float sampling time
        if re.search(r'.',lines[3]):
            pattern = r'([0-9]+\.?[0-9]*)' 
            self.__sampling_period = float(re.findall(pattern, lines
            [3])[0])
        # load int sampling time
        else:
            pattern = r'([0-9])' 
            self.__sampling_period = float(re.findall(pattern, lines[3])[0])


    def getSimulationTime(self):
        return self.__simulation_time

    def getSamplingPeriod(self):
        return self.__sampling_period

