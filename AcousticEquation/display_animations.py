import sys

sys.path.insert(1, '../Lib/')

from discretisation import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML
import time






def get_anim(parameter_toplot, history, limits, z, propagation) : 
    list_anim = []

    fig,ax = plt.subplots(figsize=(10,3))
    plt.rc('font', size=14)
    ax.set_xlabel("Range (km)")
    ax.set_ylabel(parameter_toplot)
    ax.plot(z[index_source]/1000, 0, 'xr')
    ax.plot(z[index_receivers]/1000, 0, 'xg')
    ax.axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax.grid()
    ax.set_xlim(( z[0]/1000, z[-1]/1000))
    ax.set_ylim(limits)

    line0, = ax.plot([], [], lw=2)

    # initialization function: plot the background of each frame
    def init():
        line0.set_data([], [])
        return (line0,)

    # animation function. This is called sequentially
    def animate(i,z,history, parameter_toplot):
        x = z / 1000
        if parameter_toplot == "Density" :
            y = history[i].rho
        elif parameter_toplot == "Velocity" :
            y = history[i].v  
            x = (z[1:] + z[:-1])/2 / 1000
        elif parameter_toplot == "Pressure" :
            y = history[i].p 
        line0.set_data(x, y)
        #line0.set_color("g")
        return (line0,)

    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(z,history,parameter_toplot), frames=len(history), interval=20, blit=True)
    #writervideo = animation.FFMpegWriter(fps=60)
    anim.save(local_path+'/Animations/animation_'+parameter_toplot+'_'+propagation+".mp4")#, writer=writervideo)
    plt.close()


history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)


#get_anim("Velocity", history_obs[::20], (-0.002,0.002), z, 'forward')
get_anim("Pressure", history_obs[::20], (-0.6,0.8), z, 'forward')
