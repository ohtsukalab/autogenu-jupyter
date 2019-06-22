from AutoGenU import simulation_conditions as simcon
import numpy as np
import seaborn as sns
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import mpl_toolkits.mplot3d as a3
import mpl_toolkits.mplot3d.art3d as art3d 
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
            sys.exit() 
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
        self.__link1, = self.__ax.plot([], [], color = '#0063B1', linewidth = 3)
        self.__link2, = self.__ax.plot([], [], color = '#0063B1', linewidth = 3)
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
            sys.exit() 
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
        self.__cartt, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1.5)
        self.__cartr, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1.5)
        self.__cartb, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1.5)
        self.__cartl, = self.__ax.plot([], [], color = '#0063B1', linewidth = 1.5)
        self.__pole, = self.__ax.plot([], [], color = '#0063B1', linewidth = 3.5)
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
        if self.__dim_state != 12:
            print('Dimension of the state is not 12! This may not be data for simulation of hexacopter')
            sys.exit() 
        # set for drawing animation
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__time_series_state.shape[0]/self.__skip_frames)

    def __hexagon_world(self, x):
        # Configurations in the body frame
        radius = 0.25
        X_b = [radius*np.cos((1/3)*np.pi*i) for i in range(6)]
        Y_b = [radius*np.sin((1/3)*np.pi*i) for i in range(6)]
        Z_b = [0 for i in range(6)]
        # Configurations in the world frame
        X_w = x[0] + [X_b[i]*np.cos(x[3])*np.cos(x[4]) + Y_b[i]*(np.cos(x[3])*np.sin(x[4])*np.sin(x[5])-np.sin(x[3])*np.cos(x[5])) + Z_b[i]*(np.cos(x[3])*np.sin(x[4])*np.cos(x[5])+np.sin(x[3])*np.sin(x[5])) for i in range(6)]
        Y_w = x[1] + [X_b[i]*np.sin(x[3])*np.cos(x[4]) + Y_b[i]*(np.sin(x[3])*np.sin(x[4])*np.sin(x[5])+np.cos(x[3])*np.cos(x[4])) + Z_b[i]*(np.cos(x[3])*np.sin(x[4])*np.sin(x[5])-np.cos(x[3])*np.sin(x[5])) for i in range(6)]
        Z_w = x[2] + [-X_b[i]*np.sin(x[4]) + Y_b[i]*np.cos(x[4])*np.sin(x[5]) + Z_b[i]*np.cos(x[4])*np.cos(x[5]) for i in range(6)]
        return X_w, Y_w, Z_w

    def __updateAnimation(self, i):
        frame = self.__skip_frames * i
        X, Y, Z = self.__hexagon_world(self.__time_series_state[frame,:])
        self.__line1.set_data((X[0],X[1]), (Y[0],Y[1]))
        self.__line1.set_3d_properties([Z[0],Z[1]])
        self.__line2.set_data((X[1],X[2]), (Y[1],Y[2]))
        self.__line2.set_3d_properties([Z[1],Z[2]])
        self.__line3.set_data((X[2],X[3]), (Y[2],Y[3]))
        self.__line3.set_3d_properties([Z[2],Z[3]])
        self.__line4.set_data((X[3],X[4]), (Y[3],Y[4]))
        self.__line4.set_3d_properties([Z[3],Z[4]])
        self.__line5.set_data((X[4],X[5]), (Y[4],Y[5]))
        self.__line5.set_3d_properties([Z[4],Z[5]])
        self.__line6.set_data((X[5],X[0]), (Y[5],Y[0]))
        self.__line6.set_3d_properties([Z[5],Z[0]])
        self.__time_text.set_text('{0:.1f} [s]'.format(self.__simulation_conditions.getSamplingPeriod() * frame))
        return self.__line1, self.__line2, self.__line3, self.__line4, self.__line5, self.__line6, self.__time_text

    def setSkipFrames(self, skip_frames):
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__time_series_state.shape[0]/skip_frames)

    def generateAnimation(self):
        self.__fig = plt.figure(figsize=(10,10))
        self.__ax = self.__fig.add_subplot(111, projection='3d')
        self.__ax.set_xlabel('x')
        self.__ax.set_ylabel('y')
        self.__ax.set_zlabel('z')
        self.__ax.set_xlim(-4.0, 4.0)
        self.__ax.set_ylim(-4.0, 4.0)
        self.__ax.set_zlim(0.0, 8.0)
        self.__line1, = self.__ax.plot([], [], [], color = '#0063B1', linewidth=5)
        self.__line2, = self.__ax.plot([], [], [], color = '#0063B1', linewidth=5)
        self.__line3, = self.__ax.plot([], [], [], color = '#0063B1', linewidth=5)
        self.__line4, = self.__ax.plot([], [], [], color = '#0063B1', linewidth=5)
        self.__line5, = self.__ax.plot([], [], [], color = '#0063B1', linewidth=5)
        self.__line6, = self.__ax.plot([], [], [], color = '#0063B1', linewidth=5)

        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text2D(0.85, 0, '', transform=self.__ax.transAxes, fontsize = 14)
        # generate an animation
        anime = FuncAnimation(self.__fig, self.__updateAnimation, interval=self.__simulation_conditions.getSamplingPeriod()*1000, frames=self.__total_frames, blit=True)
        anime.save(self.file_header+'.mp4', writer='ffmpeg', fps = int(1/(self.__simulation_conditions.getSamplingPeriod()*self.__skip_frames)))
        print('The animation of the simlation results is generated at ' + self.file_header + '.mp4\n')

