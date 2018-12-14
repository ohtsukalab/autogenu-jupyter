import numpy as np
import seaborn as sns
import sys
import re

def plotsim(argv):
    ############################ load data ############################
    time_series_state = np.genfromtxt(argv[1] + '_state' + '.dat')
    time_series_control_input = np.genfromtxt(argv[1] + '_control_input' + '.dat')
    time_series_error = np.genfromtxt(argv[1] + '_error' + '.dat')
    simulation_conditions = open(argv[1] + '_conditions' + '.dat')

    lines = simulation_conditions.readlines()[1]
    pattern = r'([0-9]+\.?[0-9]*)'
    lists = re.findall(pattern,lines)
    simulation_time = float(lists[0])
    simulation_conditions.close()

    # replace NaN with 0 in simulation data
    time_series_state[np.isnan(time_series_state)] = 0
    time_series_control_input[np.isnan(time_series_control_input)] = 0
    time_series_error[np.isnan(time_series_error)] = 0

    # set dimensions
    dim_state = time_series_state.shape[1]
    num_time_steps = time_series_state.shape[0]
    if time_series_control_input.shape[0] == time_series_control_input.size:
        dim_control_input = 1
    else:
        dim_control_input = time_series_control_input.shape[1]

    # set time sequence
    time_step = simulation_time/num_time_steps
    time_sequence = np.arange(0, simulation_time, time_step)

    # set the layout of the graphs 
    num_plot_variables = dim_state+dim_control_input+1
    num_plot_x = int(np.floor(num_plot_variables/np.sqrt(num_plot_variables)))
    num_plot_y = int(np.ceil(num_plot_variables/num_plot_x))


    ############################ plot data ############################
    scale = 1.5 ## adjust figsize
    sns.set()
    sns.set_style("ticks")
    sns.set_palette("deep") 
    sns.set_context("paper")
    sns.mpl.pyplot.figure(figsize=(num_plot_x*3*scale, num_plot_y*scale*1.2)) #graph size
    sns.mpl.pyplot.rcParams['font.size'] = scale*10/num_plot_variables #font size
    sns.mpl.pyplot.subplots_adjust(wspace= 0.3*scale, hspace= 0.5*scale) #space between graphs
    sns.mpl.pyplot.rcParams['lines.linewidth'] = 1.2 #linewidth
    sns.mpl.pyplot.rc('mathtext', **{'rm':'serif', 'it':'serif:itelic', 'bf':'serif:bold', 'fontset':'cm'})
    sns.mpl.pyplot.rcParams['axes.linewidth'] = 0.5
    sns.mpl.pyplot.rcParams['xtick.direction'] = 'in'
    sns.mpl.pyplot.rcParams['ytick.direction'] = 'in'
    sns.mpl.pyplot.rcParams['pdf.fonttype'] = 42
    sns.mpl.pyplot.rcParams['ps.fonttype'] = 42

    for i in range(dim_state):
        sns.mpl.pyplot.subplot(num_plot_y, num_plot_x, i+1)
        sns.mpl.pyplot.plot(time_sequence, time_series_state[:,i])
        sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
        sns.mpl.pyplot.ylabel(r'$x_{}$'.format(i+1))
        sns.mpl.pyplot.xlim(0, simulation_time)

    if dim_control_input > 1:
        for i in range(dim_control_input):
            sns.mpl.pyplot.subplot(num_plot_y, num_plot_x, i+dim_state+1)
            sns.mpl.pyplot.plot(time_sequence, time_series_control_input[:,i])
            sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
            sns.mpl.pyplot.ylabel(r'$u_{}$'.format(i+1))
            sns.mpl.pyplot.xlim(0, simulation_time)
    else:
        sns.mpl.pyplot.subplot(num_plot_y, num_plot_x, dim_state+1)
        sns.mpl.pyplot.plot(time_sequence, time_series_control_input)
        sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
        sns.mpl.pyplot.ylabel(r'$u$')
        sns.mpl.pyplot.xlim(0, simulation_time)

    sns.mpl.pyplot.subplot(num_plot_y, num_plot_x, dim_state+dim_control_input+1)
    sns.mpl.pyplot.plot(time_sequence, time_series_error)
    sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
    sns.mpl.pyplot.ylabel(r'$\| F \|$')
    sns.mpl.pyplot.xlim(0, simulation_time)


    ############################ show or save data ############################
    if(len(argv) == 2):
        sns.mpl.pyplot.show()
    elif(len(argv) == 3):
        if(argv[2] == 'show' or argv[2] == 'SHOW'):
            sns.mpl.pyplot.show()
        if(argv[2] == 'save' or argv[2] == 'SAVE'):
            print('Input file name')
            ans = input('>> ')
            sns.mpl.pyplot.savefig(ans+'.pdf', bbox_inches="tight", pad_inches=0.1)
    elif(len(argv) == 4):
        sns.mpl.pyplot.savefig(argv[3]+'.pdf', bbox_inches="tight", pad_inches=0.1)


if __name__ == "__main__":
    argc = len(sys.argv)

    if argc != 2 and argc != 3 and argc != 4:
        print('Input as  $ python3 plotsim.py header  or  $ python3 plotsim.py header save')
        quit()
    
    plotsim(sys.argv)


