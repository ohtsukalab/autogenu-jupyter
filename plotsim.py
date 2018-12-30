import numpy as np
import seaborn as sns
import sys
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



class PlotSimulationData:
    def __init__(self, argv):
        # check the argv
        if len(argv) == 2 or len(argv) == 3 or len(argv) == 4:
            file_header = argv[1]
        else:
            print('Input as  $ python3 plotsim.py header  or  $ python3 plotsim.py header save')
            quit()

        # load data
        self.__time_series_state = np.genfromtxt(file_header + '_state' + '.dat')
        self.__time_series_control_input = np.genfromtxt(file_header + '_control_input' + '.dat')
        self.__time_series_error = np.genfromtxt(file_header + '_error' + '.dat')
        self.__simulation_conditions = SimulationConditions(file_header)
        self.__time_sequence = np.linspace(0, self.__simulation_conditions.getSimulationTime(), self.__time_series_state.shape[0])

        # replace NaN with 0
        self.__time_series_state[np.isnan(self.__time_series_state)] = 0
        self.__time_series_control_input[np.isnan(self.__time_series_control_input)] = 0
        self.__time_series_error[np.isnan(self.__time_series_error)] = 0

        # set dimensions
        self.__dim_state = self.__time_series_state.shape[1]
        if self.__time_series_control_input.shape[0] == self.__time_series_control_input.size:
            self.__dim_control_input = 1
        else:
            self.__dim_control_input = self.__time_series_control_input.shape[1]
        
        # set the layout of the graphs and settings  
        self.__num_plot_variables = self.__dim_state+self.__dim_control_input+1
        self.__num_plot_x = int(np.floor(self.__num_plot_variables/np.sqrt(self.__num_plot_variables)))
        self.__num_plot_y = int(np.ceil(self.__num_plot_variables/self.__num_plot_x))

        # adjust figure size 
        self.__figure_scale = 1
        self.__font_scale = 1
        self.__space_scale = 1

        # set format of seaborn
        sns.set_style("ticks")
        sns.set_palette("deep") 
        sns.set_context("paper")
        sns.mpl.pyplot.rc('mathtext', **{'rm':'serif', 'it':'serif:itelic', 'bf':'serif:bold', 'fontset':'cm'})
        sns.mpl.pyplot.rcParams['xtick.direction'] = 'in'
        sns.mpl.pyplot.rcParams['ytick.direction'] = 'in'
        sns.mpl.pyplot.rcParams['pdf.fonttype'] = 42
        sns.mpl.pyplot.rcParams['ps.fonttype'] = 42

        # set the width of lines and axes
        sns.mpl.pyplot.rcParams['lines.linewidth'] = 1
        sns.mpl.pyplot.rcParams['axes.linewidth'] = 0.5


    def setScales(self, figure_scale, font_scale, space_scale):
        self.__figure_scale = figure_scale
        self.__font_scale = font_scale
        self.__space_scale = space_scale


    def drawGraphs(self, argv):
        # set figure size
        sns.mpl.pyplot.figure(figsize=(2.5*self.__num_plot_x*self.__figure_scale, self.__num_plot_y*self.__figure_scale)) 
        # set font size
        sns.mpl.pyplot.rcParams['font.size'] = self.__font_scale*10/self.__num_plot_variables 
        # set the space between the graphs
        sns.mpl.pyplot.subplots_adjust(wspace= self.__space_scale/self.__num_plot_variables, hspace= 2*self.__space_scale/self.__num_plot_variables) 
    
        for i in range(self.__dim_state):
            sns.mpl.pyplot.subplot(self.__num_plot_y, self.__num_plot_x, i+1)
            sns.mpl.pyplot.plot(self.__time_sequence, self.__time_series_state[:,i])
            sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
            sns.mpl.pyplot.ylabel(r'$x_{}$'.format(i+1))
            sns.mpl.pyplot.xlim(0, self.__simulation_conditions.getSimulationTime())

        if self.__dim_control_input > 1:
            for i in range(self.__dim_control_input):
                sns.mpl.pyplot.subplot(self.__num_plot_y, self.__num_plot_x, i+self.__dim_state+1)
                sns.mpl.pyplot.plot(self.__time_sequence, self.__time_series_control_input[:,i])
                sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
                sns.mpl.pyplot.ylabel(r'$u_{}$'.format(i+1))
                sns.mpl.pyplot.xlim(0, self.__simulation_conditions.getSimulationTime())
        else:
            sns.mpl.pyplot.subplot(self.__num_plot_y, self.__num_plot_x, self.__dim_state+1)
            sns.mpl.pyplot.plot(self.__time_sequence, self.__time_series_control_input)
            sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
            sns.mpl.pyplot.ylabel(r'$u$')
            sns.mpl.pyplot.xlim(0, self.__simulation_conditions.getSimulationTime())

        sns.mpl.pyplot.subplot(self.__num_plot_y, self.__num_plot_x, self.__dim_state+self.__dim_control_input+1)
        sns.mpl.pyplot.plot(self.__time_sequence, self.__time_series_error)
        sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
        sns.mpl.pyplot.ylabel(r'$\| F \|$')
        sns.mpl.pyplot.xlim(0, self.__simulation_conditions.getSimulationTime())

        if(len(argv) == 2):
            sns.mpl.pyplot.show()
        elif(len(argv) == 3):
            if(argv[2] == 'show' or argv[2] == 'SHOW'):
                sns.mpl.pyplot.show()
            if(argv[2] == 'save' or argv[2] == 'SAVE'):
                print('Input save file name')
                ans = input('>> ')
                sns.mpl.pyplot.savefig(ans+'.pdf', bbox_inches="tight", pad_inches=0.1)
        elif(len(argv) == 4):
            sns.mpl.pyplot.savefig(argv[3]+'.pdf', bbox_inches="tight", pad_inches=0.1)



if __name__ == "__main__":
    plotsim = PlotSimulationData(sys.argv)
    plotsim.setScales(2, 2, 1.5)
    plotsim.drawGraphs(sys.argv)
