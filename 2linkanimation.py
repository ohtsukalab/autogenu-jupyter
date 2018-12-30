import numpy as np
import seaborn as sns
from matplotlib.animation import FuncAnimation
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


class TwoLinkAnimation:
    def __init__(self, argv):
        # check the argv
        if len(argv) == 2 or len(argv) == 3 or len(argv) == 4:
            file_header = argv[1]
        else:
            print('Input as  $ python3 plotsim.py header  or  $ python3 plotsim.py header save')
            quit()

        # load data
        self.__time_series_state = np.genfromtxt(file_header + '_state' + '.dat')
        self.__simulation_conditions = SimulationConditions(file_header)
        self.__time_sequence = np.arange(0, self.__simulation_conditions.getSimulationTime(), self.__simulation_conditions.getSamplingPeriod())

        # replace NaN with 0
        self.__time_series_state[np.isnan(self.__time_series_state)] = 0

        # check dimensions
        self.__dim_state = self.__time_series_state.shape[1]
        if self.__dim_state != 4:
            print('Dimension of the state is not 4! This may not be data for simulation of 2link arm')
            quit()
        
        # set drawing range
        self.__x_min = -2.5
        self.__x_max = 2.5
        self.__y_min = -2.5
        self.__y_max = 2.5

        # set for drawing animation
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__time_series_state.shape[0]/self.__skip_frames)


    def __updateAnimation(self, i):
        frame = self.__skip_frames * i
        self.__x1 = np.sin(self.__time_series_state[frame,0])
        self.__y1 = - np.cos(self.__time_series_state[frame,0])
        self.__x2 = self.__x1 + np.sin(self.__time_series_state[frame,0]+self.__time_series_state[frame,1])
        self.__y2 = self.__y1 - np.cos(self.__time_series_state[frame,0]+self.__time_series_state[frame,1])
        self.__link1.set_data((0,self.__x1),(0,self.__y1))
        self.__link2.set_data((self.__x1,self.__x2),(self.__y1,self.__y2))
        self.__time_text.set_text('{0:.1f} [s]'.format(self.__simulation_conditions.getSamplingPeriod() * frame))
        return self.__link1, self.__link2, self.__time_text


    def setSkipFrames(self, skip_frames):
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__time_series_state.shape[0]/skip_frames)


    def generateAnimation(self, argv):
        self.__fig = sns.mpl.pyplot.figure()
        self.__ax = sns.mpl.pyplot.axes(xlim=(self.__x_min, self.__x_max), ylim=(self.__y_min, self.__y_max))

        self.__link1, = self.__ax.plot([], [], color = '#0063B1', linewidth = 2)
        self.__link2, = self.__ax.plot([], [], color = '#0063B1', linewidth = 2)

        self.__ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text(0.85, 0.05, '', transform=self.__ax.transAxes, fontsize = 12)

        anime = FuncAnimation(self.__fig, self.__updateAnimation, interval=self.__simulation_conditions.getSamplingPeriod()*1000, frames=self.__total_frames, blit=True)

        if(len(argv) == 2):
            sns.mpl.pyplot.show()
        elif(len(argv) == 3):
            if(argv[2] == 'show' or argv[2] == 'SHOW'):
                sns.mpl.pyplot.show()
            if(argv[2] == 'save' or argv[2] == 'SAVE'):
                print('Input save file name')
                ans = input('>> ')
                anime.save(ans+'.mp4', writer='ffmpeg', fps = int(1/(self.__simulation_conditions.getSamplingPeriod()*self.__skip_frames)))
        elif(len(argv) == 4):
            anime.save(argv[3]+'.mp4', writer='ffmpeg', fps = int(1/(self.__simulation_conditions.getSamplingPeriod()*self.__skip_frames)))



if __name__ == "__main__":
    animation = TwoLinkAnimation(sys.argv)
    animation.setSkipFrames(1)
    animation.generateAnimation(sys.argv)



