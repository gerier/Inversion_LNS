#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
sys.path.insert(0, '../Lib/')


from parameters import *
from discretisation import *
from linearised_navier_stokes_acoustic_adjoint_splitchi import *


import numpy as np
import time


# INITIALISATION
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 

# CHOICE
which_method = "forward"

#%% BACKWARD RESOLUTION
if which_method == "forward" : 
    is_reverse = False
    t_end, U_end, history_forwardapriori = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, is_reverse)
    history_reverse = [history_forwardapriori[t] for t in range(len(history_forwardapriori)-1,-1,-1)]
    U_start = history_reverse[-1]

elif which_method == "backward" :
    is_reverse = True 
    t_start, U_start, history_reverse = time_scheme(get_minus_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, reversed_source, is_reverse, history_forwardapriori)

elif which_method  == "checkpointing" : 
    is_reverse = False
    t_start, U_start, history_reverse = time_scheme(get_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse, history_forwardapriori)


# SAVE

np.save(local_path+"/BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), np.array(history_reverse))





#%% RESULTS BACKWARD

U_start.plot(z,"After Tmax = 0")

if display_anim : 
    get_anim("Density", history_reverse[::10], (-0.02,0.02), z, 'backward')
    get_anim("Velocity", history_reverse[::10], (-4,4), z, 'backward')
    get_anim("Pressure", history_reverse[::10], (-2000,2000), z, 'backward')

#%% SUPERIMPOSITION OF BACKWARD AND FORWARD SOLUTION

# load observations     
history_obs = np.load(local_path+"/BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)

# comparison
n = len(history_reverse)
t_ax = np.append(t_ax, [t_ax[-1]+dt])
t_ax = t_ax[::2]
for i in range(0,n,int(n/4)):
    fig,ax = plt.subplots(3,2, figsize=(10,7))
    # forward
    ax[0,0].plot(z,history_obs[i].rho, 'b', label='Forward')
    ax[1,0].plot(z,history_obs[i].p, 'b')
    ax[2,0].plot( (z[1:] + z[:-1])/2,history_obs[i].v, 'b')
    # backward
    ax[0,0].plot(z,history_reverse[n-1-i].rho, 'r', label='Backward')
    ax[1,0].plot(z,history_reverse[n-1-i].p, 'r')
    ax[2,0].plot( (z[1:] + z[:-1])/2,history_reverse[n-1-i].v, 'r')
    # difference
    ax[0,1].plot(z,abs(history_reverse[n-1-i].rho - history_obs[i].rho), 'k', label='Backward')
    ax[1,1].plot(z,abs(history_reverse[n-1-i].p - history_obs[i].p), 'k')
    ax[2,1].plot( (z[1:] + z[:-1])/2,abs(history_reverse[n-1-i].v - history_obs[i].v), 'k')
    # color where the model is different 
    ax[0,1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax[1,1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax[2,1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    # element to display 
    ax[0,0].legend()
    ax[0,0].set_xlabel("Distance on axis x (km)")
    ax[1,0].set_xlabel("Distance on axis x (km)")
    ax[2,0].set_xlabel("Distance on axis x (km)")
    ax[0,0].set_ylabel("Density")
    ax[1,0].set_ylabel("Pressure")
    ax[2,0].set_ylabel("Velocity")
    ax[0,1].set_xlabel("Distance on axis x (km)")
    ax[1,1].set_xlabel("Distance on axis x (km)")
    ax[2,1].set_xlabel("Distance on axis x (km)")
    ax[0,1].set_ylabel("Difference Density")
    ax[1,1].set_ylabel("Difference Pressure")
    ax[2,1].set_ylabel("Difference Velocity")
    ax[0,1].grid()
    ax[1,1].grid()
    ax[2,0].grid()
    fig.suptitle("Comparison of forward and backward solution at time (forward) t = "+str(t_ax[i]))
    plt.show()
   

