import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import os
import sys


class TwoLinkArm(object):
    """ Generates the animation of the simulation results of a 2link arm.

        Attributes: 
            set_skip_frames(skip_frames): Sets how many frames you want to 
                skip in generating the animation. In the default settings, 
                skip_frames = 1.
            generate_animation(): Draws an animation of the simulation reult
                and saves it as a .mp4 files.
    """

    def __init__(self, log_dir, log_name: str):
        """ Inits TwoLinkArm with loading the simulation results. """
        # Loads the simulation data.
        self.__log_dir = log_dir
        self.__log_name = log_name
        self.__t_data = np.genfromtxt(os.path.join(log_dir, log_name+"_t.log"))
        self.__x_data = np.genfromtxt(os.path.join(log_dir, log_name+"_x.log"))
        self.__sampling_time = self.__t_data[1] - self.__t_data[0] 
        # Replaces NaN with 0.
        self.__x_data[np.isnan(self.__x_data)] = 0
        # Checks the dimension of the state.
        self.__dim_x = self.__x_data.shape[1]
        if self.__dim_x != 4:
            print(
                'Dimension of the state is not 4!\n'
                'This may not be data for simulation of a 2link arm\n'
            )
            sys.exit() 
        # Sets drawing range.
        self.__x_min = -4
        self.__x_max = 4
        self.__y_min = -3
        self.__y_max = 3
        # Sets the length of a link. 
        self.__length = 1.2
        # Sets frames for drawing the animation.
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__x_data.shape[0]/self.__skip_frames)

    def set_skip_frames(self, skip_frames):
        """ Set how many frames you want to skip in generating the animation.

            Args:
                skip_frames: A number of frames to skip.
        """
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__x_data.shape[0]/skip_frames)

    def generate_animation(self):
        """ Generates the animation and saves it as a .mp4 file. """
        self.__fig = plt.figure(figsize=(8, 6))
        self.__ax = plt.axes(
            xlim=(self.__x_min, self.__x_max), 
            ylim=(self.__y_min, self.__y_max)
        )
        self.__link1, = self.__ax.plot([], [], color='#0063B1', linewidth=3)
        self.__link2, = self.__ax.plot([], [], color='#0063B1', linewidth=3)
        self.__ax.tick_params(
            labelbottom=False, 
            labelleft=False, 
            labelright=False, 
            labeltop=False
        )
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text(
            0.85, 
            0.05, 
            '', 
            transform=self.__ax.transAxes, 
            fontsize=14
        )
        # Generates an animation.
        anime = FuncAnimation(
            self.__fig, 
            self.__update_animation, 
            interval=self.__sampling_time*1000*self.__skip_frames, 
            frames=self.__total_frames, 
            blit=True
        )
        anime.save(
            os.path.join(self.__log_dir, self.__log_name+'.mp4'),
            writer='ffmpeg', 
            fps=int(
                1/(self.__sampling_time*self.__skip_frames)
            )
        )
        print(
            'The animation of the simlation results is generated at '
            +self.__log_dir
        )

    def __update_animation(self, i):
        frame = self.__skip_frames * i
        state = self.__x_data[frame, :]
        self.__x1 = self.__length * np.sin(state[0])
        self.__y1 = - self.__length * np.cos(state[0])
        self.__x2 = self.__x1 + self.__length * np.sin(state[0]+state[1])
        self.__y2 = self.__y1 - self.__length * np.cos(state[0]+state[1])
        self.__link1.set_data((0, self.__x1), (0, self.__y1))
        self.__link2.set_data((self.__x1, self.__x2), (self.__y1, self.__y2))
        self.__time_text.set_text(
            '{0:.1f} [s]'.format(self.__sampling_time*frame)
        )
        return self.__link1, self.__link2, self.__time_text


class CartPole(object):
    """ Generates the animation of the simulation results of a cart pole.

        Attributes: 
            set_skip_frames(skip_frames): Sets how many frames you want to 
                skip in generating the animation. In the default settings, 
                skip_frames = 1.
            generate_animation(): Draws an animation of the simulation reult
                and saves it as a .mp4 files.
    """

    def __init__(self, log_dir, log_name: str):
        """ Inits CartPole with loading the simulation results. """
        # Loads the simulation data.
        self.__log_dir = log_dir
        self.__log_name = log_name
        self.__t_data = np.genfromtxt(os.path.join(log_dir, log_name+"_t.log"))
        self.__x_data = np.genfromtxt(os.path.join(log_dir, log_name+"_x.log"))
        self.__sampling_time = self.__t_data[1] - self.__t_data[0] 
        # Replaces NaN with 0.
        self.__x_data[np.isnan(self.__x_data)] = 0
        # Checks the dimension of the state.
        self.__dim_x = self.__x_data.shape[1]
        if self.__dim_x != 4:
            print(
                'Dimension of the state is not 4!\n'
                'This may not be data for simulation of a cartpole\n'
            )
            sys.exit() 
        # Sets the drawing range.
        xabsmax = max(
            abs(np.amin(self.__x_data[:, 0])), 
            abs(np.amax(self.__x_data[:, 0]))
        )
        self.__x_min = - xabsmax - 2.5 
        self.__x_max = xabsmax + 2.5
        self.__y_min = - (self.__x_max - self.__x_min) * 0.3 * (6/8)
        self.__y_max = (self.__x_max - self.__x_min) * 0.7 * (6/8)
        self.__cart_width = 1.5 
        self.__cart_height = 0.75 
        self.__pole_length = 1.5
        # Sets frames for drawing the animation.
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__x_data.shape[0]/self.__skip_frames)

    def set_skip_frames(self, skip_frames):
        """ Set how many frames you want to skip in generating the animation.

            Args:
                skip_frames: A number of frames to skip.
        """
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__x_data.shape[0]/skip_frames)

    def generate_animation(self):
        """ Generates the animation and saves it as a .mp4 file. """
        self.__fig = plt.figure(figsize=(8, 6))
        self.__ax = plt.axes(
            xlim=(self.__x_min, self.__x_max), 
            ylim=(self.__y_min, self.__y_max)
        )
        self.__ground, = self.__ax.plot([], [], color='#0063B1', linewidth=0.5)
        self.__cartt, = self.__ax.plot([], [], color='#0063B1', linewidth=1.5)
        self.__cartr, = self.__ax.plot([], [], color='#0063B1', linewidth=1.5)
        self.__cartb, = self.__ax.plot([], [], color='#0063B1', linewidth=1.5)
        self.__cartl, = self.__ax.plot([], [], color='#0063B1', linewidth=1.5)
        self.__pole, = self.__ax.plot([], [], color='#0063B1', linewidth=3.5)
        self.__ax.tick_params(
            labelbottom=False, 
            labelleft=False, 
            labelright=False, 
            labeltop=False
        )
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text(
            0.85, 0.05, 
            '', 
            transform=self.__ax.transAxes, 
            fontsize=14
        )
        # Generates an animation.
        anime = FuncAnimation(
            self.__fig, 
            self.__update_animation, 
            interval=self.__sampling_time*1000*self.__skip_frames, 
            frames=self.__total_frames, 
            blit=True
        )
        anime.save(
            os.path.join(self.__log_dir, self.__log_name+'.mp4'),
            writer='ffmpeg', 
            fps=int(
                1/(self.__sampling_time*self.__skip_frames)
            )
        )
        print(
            'The animation of the simlation results is generated at '
            +self.__log_dir
        )

    def __update_animation(self, i):
        frame = self.__skip_frames * i
        state = self.__x_data[frame, :]
        self.__xc = state[0]
        self.__yc = 0
        self.__xp = self.__xc + self.__pole_length * np.sin(state[1])
        self.__yp = 0.5 * self.__cart_height - self.__pole_length*np.cos(state[1])
        self.__ground.set_data((self.__x_min, self.__x_max), (0, 0))
        self.__cartt.set_data(
            (self.__xc-0.5*self.__cart_width, self.__xc+0.5*self.__cart_width), 
            (self.__cart_height, self.__cart_height)
        )
        self.__cartb.set_data(
            (self.__xc-0.5*self.__cart_width, self.__xc+0.5*self.__cart_width), 
            (0, 0)
        )
        self.__cartr.set_data(
            (self.__xc+0.5*self.__cart_width, self.__xc+0.5*self.__cart_width), 
            (0, self.__cart_height)
        )
        self.__cartl.set_data(
            (self.__xc-0.5*self.__cart_width, self.__xc-0.5*self.__cart_width), 
            (0, self.__cart_height)
        )
        self.__pole.set_data(
            (self.__xc, self.__xp), 
            (0.5*self.__cart_height, self.__yp)
        )
        self.__time_text.set_text(
            '{0:.1f} [s]'.format(self.__sampling_time*frame)
        )
        return (
            self.__ground, self.__cartt, self.__cartb, self.__cartr, self.__cartl, 
            self.__pole, self.__time_text
        )


class Hexacopter(object):
    """ Generates the animation of the simulation results of a hexacopter.

        Attributes: 
            set_skip_frames(skip_frames): Sets how many frames you want to 
                skip in generating the animation. In the default settings, 
                skip_frames = 1.
            generate_animation(): Draws an animation of the simulation reult
                and saves it as a .mp4 files.
    """

    def __init__(self, log_dir, log_name: str):
        """ Inits Hexacopter with loading the simulation results. """
        # Loads the simulation data.
        self.__log_dir = log_dir
        self.__log_name = log_name
        self.__t_data = np.genfromtxt(os.path.join(log_dir, log_name+"_t.log"))
        self.__x_data = np.genfromtxt(os.path.join(log_dir, log_name+"_x.log"))
        self.__sampling_time = self.__t_data[1] - self.__t_data[0] 
        # Replaces NaN with 0.
        self.__x_data[np.isnan(self.__x_data)] = 0
        # Checks the dimension of the state.
        self.__dim_x = self.__x_data.shape[1]
        if self.__dim_x != 12:
            print(
                'Dimension of the state is not 12!\n'
                'This may not be data for simulation of hexacopter\n'
            )
            sys.exit() 
        self.__radius = 0.25
        # Sets frames for drawing the animation.
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__x_data.shape[0]/self.__skip_frames)

    def set_skip_frames(self, skip_frames):
        """ Set how many frames you want to skip in generating the animation.

            Args:
                skip_frames: A number of frames to skip.
        """
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__x_data.shape[0]/skip_frames)

    def generate_animation(self):
        """ Generates the animation and saves it as a .mp4 file. """
        self.__fig = plt.figure(figsize=(10, 10))
        self.__ax = self.__fig.add_subplot(111, projection='3d')
        self.__ax.set_xlabel('x')
        self.__ax.set_ylabel('y')
        self.__ax.set_zlabel('z')
        self.__ax.set_xlim(-4.0, 4.0)
        self.__ax.set_ylim(-4.0, 4.0)
        self.__ax.set_zlim(0.0, 8.0)
        self.__line1, = self.__ax.plot([], [], [], color='#0063B1', linewidth=5)
        self.__line2, = self.__ax.plot([], [], [], color='#0063B1', linewidth=5)
        self.__line3, = self.__ax.plot([], [], [], color='#0063B1', linewidth=5)
        self.__line4, = self.__ax.plot([], [], [], color='#0063B1', linewidth=5)
        self.__line5, = self.__ax.plot([], [], [], color='#0063B1', linewidth=5)
        self.__line6, = self.__ax.plot([], [], [], color='#0063B1', linewidth=5)
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text2D(
            0.85, 
            0, 
            '', 
            transform=self.__ax.transAxes, 
            fontsize=14
        )
        # Generates an animation.
        anime = FuncAnimation(
            self.__fig, 
            self.__update_animation, 
            interval=self.__sampling_time*1000*self.__skip_frames, 
            frames=self.__total_frames, 
            blit=True
        )
        anime.save(
            os.path.join(self.__log_dir, self.__log_name+'.mp4'),
            writer='ffmpeg', 
            fps=int(
                1/(self.__sampling_time*self.__skip_frames)
            )
        )
        print(
            'The animation of the simlation results is generated at '
            +self.__log_dir
        )

    def __update_animation(self, i):
        frame = self.__skip_frames * i
        X, Y, Z = self.__hexagon_world(self.__x_data[frame, :])
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
        self.__time_text.set_text(
            '{0:.1f} [s]'.format(self.__sampling_time*frame)
        )
        return (
            self.__line1, self.__line2, self.__line3, self.__line4, 
            self.__line5, self.__line6, self.__time_text
        )

    def __hexagon_world(self, x):
        # Configurations in the body frame
        X_b = [self.__radius*np.cos((1/3)*np.pi*i) for i in range(6)]
        Y_b = [self.__radius*np.sin((1/3)*np.pi*i) for i in range(6)]
        Z_b = [0 for i in range(6)]
        # Configurations in the world frame
        X_w = x[0] + [
            X_b[i]*np.cos(x[3])*np.cos(x[4])
            +Y_b[i]*(
                np.cos(x[3])*np.sin(x[4])*np.sin(x[5])-np.sin(x[3])*np.cos(x[5])
            )
            +Z_b[i]*(
                np.cos(x[3])*np.sin(x[4])*np.cos(x[5])+np.sin(x[3])*np.sin(x[5])
            ) 
            for i in range(6)
        ]
        Y_w = x[1] + [
            X_b[i]*np.sin(x[3])*np.cos(x[4])
            +Y_b[i]*(
                np.sin(x[3])*np.sin(x[4])*np.sin(x[5])+np.cos(x[3])*np.cos(x[4])
            )
            +Z_b[i]*(
                np.cos(x[3])*np.sin(x[4])*np.sin(x[5])-np.cos(x[3])*np.sin(x[5])
            ) 
            for i in range(6)
        ]
        Z_w = x[2] + [
            -X_b[i]*np.sin(x[4]) 
            +Y_b[i]*np.cos(x[4])*np.sin(x[5])
            +Z_b[i]*np.cos(x[4])*np.cos(x[5]) 
            for i in range(6)
        ]
        return X_w, Y_w, Z_w


class MobileRobot(object):
    """ Generates the animation of the simulation results of a cart pole.

        Attributes: 
            set_skip_frames(skip_frames): Sets how many frames you want to 
                skip in generating the animation. In the default settings, 
                skip_frames = 1.
            generate_animation(): Draws an animation of the simulation reult
                and saves it as a .mp4 files.
    """

    def __init__(self, log_dir, log_name: str, vx_ref, X1, Y1, R1, X2, Y2, R2):
        """ Inits CartPole with loading the simulation results. """
        # Loads the simulation data.
        self.__log_dir = log_dir
        self.__log_name = log_name
        self.__t_data = np.genfromtxt(os.path.join(log_dir, log_name+"_t.log"))
        self.__x_data = np.genfromtxt(os.path.join(log_dir, log_name+"_x.log"))
        self.__sampling_time = self.__t_data[1] - self.__t_data[0] 
        # Replaces NaN with 0.
        self.__x_data[np.isnan(self.__x_data)] = 0
        # Checks the dimension of the state.
        self.__dim_x = self.__x_data.shape[1]
        if self.__dim_x != 3:
            print(
                'Dimension of the state is not 3!\n'
                'This may not be data for simulation of a mobile robot\n'
            )
            sys.exit() 
        # Sets the drawing range.
        self.__robot_length = 0.2
        self.__robot_width = 0.15
        hlength = 0.5 * self.__robot_length
        hwidth = 0.5 * self.__robot_width
        margin = np.sqrt(hlength**2 + hwidth**2)
        self.__x_min = np.min(self.__x_data[:, 0]) - margin
        self.__x_max = np.max(self.__x_data[:, 0]) + margin
        self.__y_min = np.min(self.__x_data[:, 1]) - margin
        self.__y_max = np.max(self.__x_data[:, 1]) + margin
        # Sets reference trajectory.
        self.__vx_ref = vx_ref
        # Sets obstacles.
        self.__X1, self.__Y1, self.__R1 = X1, Y1, R1-margin
        self.__X2, self.__Y2, self.__R2 = X2, Y2, R2-margin
        # Sets frames for drawing the animation.
        self.__skip_frames = 1
        self.__total_frames = (int)(self.__x_data.shape[0]/self.__skip_frames)

    def set_skip_frames(self, skip_frames):
        """ Set how many frames you want to skip in generating the animation.

            Args:
                skip_frames: A number of frames to skip.
        """
        self.__skip_frames = skip_frames
        self.__total_frames = (int)(self.__x_data.shape[0]/skip_frames)

    def generate_animation(self):
        """ Generates the animation and saves it as a .mp4 file. """
        margin = max(self.__robot_length, self.__robot_width)
        self.__fig = plt.figure(figsize=(10, 5))
        xrange = self.__x_max - self.__x_min
        yrange = self.__y_max - self.__y_min
        max_range = max(xrange, yrange) + margin
        self.__ax = plt.axes(
            xlim=(self.__x_min-margin, self.__x_min+max_range), 
            ylim=(-0.5*(5/10)*max_range, 0.5*(5/10)*max_range)
        )
        obs1 = patches.Circle(xy=(self.__X1, self.__Y1), radius=self.__R1, 
                                  fc='w', ec='gray', linewidth=3)
        obs2 = patches.Circle(xy=(self.__X2, self.__Y2), radius=self.__R2, 
                                  fc='w', ec='gray', linewidth=3)
        self.__ax.add_patch(obs1)
        self.__ax.add_patch(obs2)
        self.__line1, = self.__ax.plot([], [], color='#0063B1', linewidth=5)
        self.__line2, = self.__ax.plot([], [], color='#0063B1', linewidth=5)
        self.__line3, = self.__ax.plot([], [], color='#0063B1', linewidth=5)
        self.__line4, = self.__ax.plot([], [], color='#0063B1', linewidth=5)
        self.__ref, = self.__ax.plot([], color='orange', marker='.',
                                     markersize=20)
        self.__ax.tick_params(
            labelbottom=False, 
            labelleft=False, 
            labelright=False, 
            labeltop=False
        )
        self.__ax.tick_params(color='white')
        self.__time_text = self.__ax.text(
            0.85, 0.05, 
            '', 
            transform=self.__ax.transAxes, 
            fontsize=14
        )
        # Generates an animation.
        anime = FuncAnimation(
            self.__fig, 
            self.__update_animation, 
            interval=self.__sampling_time*1000*self.__skip_frames, 
            frames=self.__total_frames, 
            blit=True
        )
        anime.save(
            os.path.join(self.__log_dir, self.__log_name+'.mp4'),
            writer='ffmpeg', 
            fps=int(
                1/(self.__sampling_time*self.__skip_frames)
            )
        )
        print(
            'The animation of the simlation results is generated at '
            +self.__log_dir
        )

    def __update_animation(self, i):
        frame = self.__skip_frames * i
        p_fl, p_fr, p_bl, p_br = self.__robot_world(self.__x_data[frame, :])
        self.__line1.set_data((p_fl[0], p_fr[0]), (p_fl[1], p_fr[1]))
        self.__line2.set_data((p_fl[0], p_bl[0]), (p_fl[1], p_bl[1]))
        self.__line3.set_data((p_fr[0], p_br[0]), (p_fr[1], p_br[1]))
        self.__line4.set_data((p_bl[0], p_br[0]), (p_bl[1], p_br[1]))
        self.__ref.set_data(self.__get_time(i)*self.__vx_ref, 0)
        self.__time_text.set_text(
            '{0:.1f} [s]'.format(self.__sampling_time*frame)
        )
        return (
            self.__line1, self.__line2, self.__line3, self.__line4, self.__ref,
            self.__time_text
        )

    def __robot_world(self, x):
        hlength = 0.5*self.__robot_length
        hwidth = 0.5*self.__robot_width
        # Configurations in the world frame
        p_fl = (
            x[0] + hlength * np.cos(x[2]) - hwidth * np.sin(x[2]),
            x[1] + hlength * np.sin(x[2]) + hwidth * np.cos(x[2])
        )
        p_fr = (
            x[0] + hlength * np.cos(x[2]) + hwidth * np.sin(x[2]),
            x[1] + hlength * np.sin(x[2]) - hwidth * np.cos(x[2])
        )
        p_bl = (
            x[0] - hlength * np.cos(x[2]) - hwidth * np.sin(x[2]),
            x[1] - hlength * np.sin(x[2]) + hwidth * np.cos(x[2])
        )
        p_br = (
            x[0] - hlength * np.cos(x[2]) + hwidth * np.sin(x[2]),
            x[1] - hlength * np.sin(x[2]) - hwidth * np.cos(x[2])
        )
        return p_fl, p_fr, p_bl, p_br

    def __get_time(self, i):
        return i*self.__t_data[-1]/self.__total_frames