import numpy as np
import seaborn as sns
import sys
import re

def plotsim(argv):
    ############################ load data ############################
    xs = np.genfromtxt(argv[1] + 'x' + '.dat')
    us = np.genfromtxt(argv[1] + 'u' + '.dat')
    es = np.genfromtxt(argv[1] + 'e' + '.dat')

    xs[np.isnan(xs)] = 0
    us[np.isnan(us)] = 0
    es[np.isnan(es)] = 0


    ld = open(argv[1] + 'c' + '.dat')
    lines = ld.readlines()[1]
    pattern = r'([0-9]+\.?[0-9]*)'
    lists = re.findall(pattern,lines)
    tsim = float(lists[0])

    ld.close()

    dimx = xs.shape[1]
    if us.shape[0] == us.size:
        dimu = 1
    else:
        dimu = us.shape[1]

    ns = xs.shape[0]

    dt = tsim/ns
    ts = np.arange(0,tsim,dt)

    plotno = dimx+dimu+1
    plotnox = int(np.floor(plotno/3))
    plotnoy = int(np.ceil(plotno/plotnox))


    ############################ plot data ############################
    scale = 1.5 ## adjust figsize

    sns.set()
    sns.set_style("ticks")
    sns.set_palette("deep") 
    sns.set_context("paper")
    sns.mpl.pyplot.figure(figsize=(plotnox*3*scale,plotnoy*scale*1.2)) #graph size
    sns.mpl.pyplot.rcParams['font.size'] = scale*10/plotno #font size
    sns.mpl.pyplot.subplots_adjust(wspace= 0.3*scale, hspace= 0.5*scale) #space between graphs
    sns.mpl.pyplot.rcParams['lines.linewidth'] = 1.2 #linewidth
    sns.mpl.pyplot.rc('mathtext', **{'rm':'serif', 'it':'serif:itelic', 'bf':'serif:bold', 'fontset':'cm'})
    sns.mpl.pyplot.rcParams['axes.linewidth'] = 0.5
    sns.mpl.pyplot.rcParams['xtick.direction'] = 'in'
    sns.mpl.pyplot.rcParams['ytick.direction'] = 'in'
    sns.mpl.pyplot.rcParams['pdf.fonttype'] = 42
    sns.mpl.pyplot.rcParams['ps.fonttype'] = 42


    for i in range(dimx):
        sns.mpl.pyplot.subplot(plotnoy, plotnox, i+1)
        sns.mpl.pyplot.plot(ts, xs[:,i])
        sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
        sns.mpl.pyplot.ylabel(r'$x_{}$'.format(i+1))
        sns.mpl.pyplot.xlim(0,tsim)


    if dimu > 1:
        for i in range(dimu):
            sns.mpl.pyplot.subplot(plotnoy, plotnox, i+dimx+1)
            sns.mpl.pyplot.plot(ts, us[:,i])
            sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
            sns.mpl.pyplot.ylabel(r'$u_{}$'.format(i+1))
            sns.mpl.pyplot.xlim(0,tsim)

    else:
        sns.mpl.pyplot.subplot(plotnoy, plotnox, dimx+1)
        sns.mpl.pyplot.plot(ts, us)
        sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
        sns.mpl.pyplot.ylabel(r'$u$')
        sns.mpl.pyplot.xlim(0,tsim)


    F = np.zeros(ns)
    for i in range(ns):
        F[i] = np.sum(es[i])
        F[i] = np.sqrt(F[i])

    sns.mpl.pyplot.subplot(plotnoy, plotnox, dimx+dimu+1)
    sns.mpl.pyplot.plot(ts,F)
    sns.mpl.pyplot.xlabel(r'${\rm Time}$ $[s]$')
    sns.mpl.pyplot.ylabel(r'$\| F \|$')
    sns.mpl.pyplot.xlim(0,tsim)


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


