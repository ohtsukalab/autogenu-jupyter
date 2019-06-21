from AutoGenU import simulation_conditions as simcon
import numpy as np
import seaborn as sns
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import sys
import re


class TwoLinkArm(object):
    def __init__(self, model_name):
        # load data
        self.file_header = 'models/' + model_name + '/simulation_result/' + model_name
        self.__time_series_state = np.genfromtxt(self.file_header + '_state' + '.dat')
        self.__simulation_conditions = simcon.SimulationConditions(self.file_header)
        self.__time_sequence = np.linspace(0, self.__simulation_conditions.getSimulationTime(), self.__time_series_state.shape[0])
        # replace NaN with 0
        self.__time_series_state[np.isnan(self.__time_series_state)] = 0
        # check dimensions
        self.__dim_state = self.__time_series_state.shape[1]
        if self.__dim_state != 4:
            print('Dimension of the state is not 4! This may not be data for simulation of 2link arm')
            quit()
        # set drawing range
        self.__x_min = -4
        self.__x_max = 4
        self.__y_min = -3
        self.__y_max = 3
        self.__length = 1.2
        # set for drawing animation
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__time_series_state.shape[0]/self.__skip_frames)

    def __updateAnimation(self, i):
        frame = self.__skip_frames * i
        self.__x1 = self.__length*np.sin(self.__time_series_state[frame,0])
        self.__y1 = - self.__length*np.cos(self.__time_series_state[frame,0])
        self.__x2 = self.__x1 + self.__length*np.sin(self.__time_series_state[frame,0]+self.__time_series_state[frame,1])
        self.__y2 = self.__y1 - self.__length*np.cos(self.__time_series_state[frame,0]+self.__time_series_state[frame,1])
        self.__link1.set_data((0,self.__x1),(0,self.__y1))
        self.__link2.set_data((self.__x1,self.__x2),(self.__y1,self.__y2))
        self.__time_text.set_text('{0:.1f} [s]'.format(self.__simulation_conditions.getSamplingPeriod() * frame))
        return self.__link1, self.__link2, self.__time_text

    def setSkipFrames(self, skip_frames):
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__time_series_state.shape[0]/skip_frames)

    def generateAnimation(self):
        self.__fig = sns.mpl.pyplot.figure(figsize=(8,6))
        self.__ax = sns.mpl.pyplot.axes(xlim=(self.__x_min, self.__x_max), ylim=(self.__y_min, self.__y_max))
        self.__link1, = self.__ax.plot([], [], color = '#0063B1', linewidth = 2)
        self.__link2, = self.__ax.plot([], [], color = '#0063B1', linewidth = 2)
        self.__ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text(0.85, 0.05, '', transform=self.__ax.transAxes, fontsize = 14)
        # generate an animation
        anime = FuncAnimation(self.__fig, self.__updateAnimation, interval=self.__simulation_conditions.getSamplingPeriod()*1000, frames=self.__total_frames, blit=True)
        anime.save(self.file_header+'.mp4', writer='ffmpeg', fps = int(1/(self.__simulation_conditions.getSamplingPeriod()*self.__skip_frames)))
        print('The animation of the simlation results is generated at ' + self.file_header + '.mp4\n')


class CartPole(object):
    def __init__(self, model_name):
        # load data
        self.file_header = 'models/' + model_name + '/simulation_result/' + model_name
        self.__time_series_state = np.genfromtxt(self.file_header + '_state' + '.dat')
        self.__simulation_conditions = simcon.SimulationConditions(self.file_header)
        self.__time_sequence = np.linspace(0, self.__simulation_conditions.getSimulationTime(), self.__time_series_state.shape[0])
        # replace NaN with 0
        self.__time_series_state[np.isnan(self.__time_series_state)] = 0
        # check dimensions
        self.__dim_state = self.__time_series_state.shape[1]
        if self.__dim_state != 4:
            print('Dimension of the state is not 4! This may not be data for simulation of cartpole')
            quit()
        # set drawing range
        xabsmax = max(abs(np.amin(self.__time_series_state[:,0])), np.amax(self.__time_series_state[:,0]))
        self.__x_min = - xabsmax -2.5 
        self.__x_max =  xabsmax + 2.5
        self.__y_min = -(self.__x_max-self.__x_min)*0.3*6/8
        self.__y_max = (self.__x_max-self.__x_min)*0.7*6/8
        self.__cart_width = 1.5 
        self.__cart_height = 0.75 
        self.__pole_length = 1.5
        # set for drawing animation
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__time_series_state.shape[0]/self.__skip_frames)

    def __updateAnimation(self, i):
        frame = self.__skip_frames * i
        self.__xc = self.__time_series_state[frame,0]
        self.__yc = 0
        self.__xp = self.__xc + self.__pole_length*np.sin(self.__time_series_state[frame,1])
        self.__yp = 0.5*self.__cart_height - self.__pole_length*np.cos(self.__time_series_state[frame,1])
        self.__ground.set_data((self.__x_min, self.__x_max), (0, 0))
        self.__cartt.set_data((self.__xc-0.5*self.__cart_width, self.__xc+0.5*self.__cart_width), (self.__cart_height, self.__cart_height))
        self.__cartb.set_data((self.__xc-0.5*self.__cart_width, self.__xc+0.5*self.__cart_width), (0, 0))
        self.__cartr.set_data((self.__xc+0.5*self.__cart_width, self.__xc+0.5*self.__cart_width), (0, self.__cart_height))
        self.__cartl.set_data((self.__xc-0.5*self.__cart_width, self.__xc-0.5*self.__cart_width), (0, self.__cart_height))
        self.__pole.set_data((self.__xc, self.__xp),(0.5*self.__cart_height,self.__yp))
        self.__time_text.set_text('{0:.1f} [s]'.format(self.__simulation_conditions.getSamplingPeriod() * frame))
        return  self.__ground, self.__cartt, self.__cartb, self.__cartr, self.__cartl, self.__pole, self.__time_text

    def setSkipFrames(self, skip_frames):
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__time_series_state.shape[0]/skip_frames)

    def generateAnimation(self):
        self.__fig = sns.mpl.pyplot.figure(figsize=(8,6))
        self.__ax = sns.mpl.pyplot.axes(xlim=(self.__x_min, self.__x_max), ylim=(self.__y_min, self.__y_max))
        self.__ground, = self.__ax.plot([], [], color = '#0063B1', linewidth = 0.5)
        self.__cartt, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1)
        self.__cartr, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1)
        self.__cartb, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1)
        self.__cartl, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1)
        self.__pole, = self.__ax.plot([], [], color = '#0063B1', linewidth = 2.5)
        self.__ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text(0.85, 0.05, '', transform=self.__ax.transAxes, fontsize = 14)
        # generate an animation
        anime = FuncAnimation(self.__fig, self.__updateAnimation, interval=self.__simulation_conditions.getSamplingPeriod()*1000, frames=self.__total_frames, blit=True)
        anime.save(self.file_header+'.mp4', writer='ffmpeg', fps = int(1/(self.__simulation_conditions.getSamplingPeriod()*self.__skip_frames)))
        print('The animation of the simlation results is generated at ' + self.file_header + '.mp4\n')


class Hexacopter(object):
    def __init__(self, model_name):
        # load data
        self.file_header = 'models/' + model_name + '/simulation_result/' + model_name
        self.__time_series_state = np.genfromtxt(self.file_header + '_state' + '.dat')
        self.__simulation_conditions = simcon.SimulationConditions(self.file_header)
        self.__time_sequence = np.linspace(0, self.__simulation_conditions.getSimulationTime(), self.__time_series_state.shape[0])
        # replace NaN with 0
        self.__time_series_state[np.isnan(self.__time_series_state)] = 0
        # check dimensions
        self.__dim_state = self.__time_series_state.shape[1]
        if self.__dim_state != 4:
            print('Dimension of the state is not 4! This may not be data for simulation of 2link arm')
            quit()
        # set drawing range
        self.__x_min = -4
        self.__x_max = 4
        self.__y_min = -3
        self.__y_max = 3
        self.__length = 1.2
        # set for drawing animation
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__time_series_state.shape[0]/self.__skip_frames)

    def __hexagonal_disk_world(self, x):
        # Configurations in the body frame
        X_b = [x[0]+np.cos((2/3)*np.pi*i) for i in range(6)]
        Y_b = [x[1]+np.sin((2/3)*np.pi*i) for i in range(6)]
        Z_b = [x[2] for i in range(6)]
        # Configurations in the world frame
        X_w = [X_b[i]*np.cos(x[3])*np.cos(x[4]) + Y_b[i]*(np.cos(x[3])*np.sin(x[4])*np.sin(x[5])-np.sin(x[3])*np.cos(x[5])) + Z_b[i]*(np.cos(x[3])*np.sin(x[4])*np.cos(x[5])+np.sin(x[3])*np.sin(x[5])) for i in range(6)]
        Y_w = [X_b[i]*np.sin(x[3])*np.cos(x[4]) + Y_b[i]*(np.sin(x[3])*np.sin(x[4])*np.sin(x[5])+np.cos(x[3])*np.cos(x[4])) + Z_b[i]*(np.cos(x[3])*np.sin(x[4])*np.sin(x[5])-np.cos(x[3])*np.sin(x[5])) for i in range(6)]
        Z_w = [-X_b[i]*np.sin(x[4]) + Y_b[i]*np.cos(x[4])*np.sin(x[5]) + Z_b[i]*np.cos(x[4])*np.cos(x[5]) for i in range(6)]
        return X_w, Y_w, Z_w
        

    def __updateAnimation(self, i):
        frame = self.__skip_frames * i
        self.__x1 = self.__length*np.sin(self.__time_series_state[frame,0])
        self.__y1 = - self.__length*np.cos(self.__time_series_state[frame,0])
        self.__x2 = self.__x1 + self.__length*np.sin(self.__time_series_state[frame,0]+self.__time_series_state[frame,1])
        self.__y2 = self.__y1 - self.__length*np.cos(self.__time_series_state[frame,0]+self.__time_series_state[frame,1])
        self.__link1.set_data((0,self.__x1),(0,self.__y1))
        self.__link2.set_data((self.__x1,self.__x2),(self.__y1,self.__y2))
        self.__time_text.set_text('{0:.1f} [s]'.format(self.__simulation_conditions.getSamplingPeriod() * frame))

        poly = list(zip(self.__hexagonal_disk_world(self.__time_series_state[frame,:])))
        ax.add_collection3d(art3d.Poly3DCollection([poly],color='m'))
        return plt.Polygon(, )

    def setSkipFrames(self, skip_frames):
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__time_series_state.shape[0]/skip_frames)

    def generateAnimation(self):
        self.__fig = plt.figure(figsize=(8,6))
        self.__ax = plt.add_subplot(111, projection='3d')
        self.__x = []
        self.__y = []
        self.__z = []

        self.__link1, = self.__ax.plot([], [], color = '#0063B1', linewidth = 2)
        self.__link2, = self.__ax.plot([], [], color = '#0063B1', linewidth = 2)
        self.__ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text(0.85, 0.05, '', transform=self.__ax.transAxes, fontsize = 14)
        # generate an animation
        anime = FuncAnimation(self.__fig, self.__updateAnimation, interval=self.__simulation_conditions.getSamplingPeriod()*1000, frames=self.__total_frames, blit=True)
        anime.save(self.file_header+'.mp4', writer='ffmpeg', fps = int(1/(self.__simulation_conditions.getSamplingPeriod()*self.__skip_frames)))
        print('The animation of the simlation results is generated at ' + self.file_header + '.mp4\n')

