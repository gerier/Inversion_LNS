#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
from threading import local
sys.path.insert(1, './Lib/')

from discretisation import *
from linearised_navier_stokes_acoustic_adjoint_splitchi import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
import time

def get_kernels(hist_adjoint, hist_backprop, U0, T0, kappa, gamma, Cv, dt, dz, order_DF=4):
    # initialize
    K = np.zeros((3, len(hist_adjoint[0].rho)))
    if order_DF == 4:
        DF_C = DF_C_4
    else :
        DF_C = DF_C_2
    
    for t in range(len(hist_adjoint)-10):
        # # create variable to be more lisible
                               # at time n
        # create variable to be more lisible
        rho_a = hist_adjoint[t+10].rho     #  at time n+1/2
        rho_p = hist_backprop[t].rho    #  at time n+1/2

        p_a = hist_adjoint[t+10].p         #  at time n+1/2
        p_p = hist_backprop[t].p        #  at time n+1/2
        if t == 0:
            v_a = hist_adjoint[t+10].v
            v_p = hist_backprop[t].v
        
            dtvp = hist_backprop[t].v / dt
        else : 
            v_a = (hist_adjoint[t].v + hist_adjoint[t+10].v)/2              # mean at time n
            v_p = (hist_backprop[t-1].v + hist_backprop[t].v)/2            # mean at time n

            dtvp = (hist_backprop[t].v - hist_backprop[t-1].v) / dt 

        K += dt * get_kernels_centered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dtvp, dz, DF_C)

    return K

# LOAD OBSERVATIONS AND SOLUTION WITH A PRIORI MODEL (in reversed time basis)
history_obs = np.load(local_path+"/BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_reverse = np.load(local_path+"/BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_adjoint = np.load(local_path+"/BackUps/adjoint_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)


# COMPUTE KERNELS
K = get_kernels(history_adjoint, history_reverse, U0, T0, 0, gamma, Cv, dt, h) 


# PLOT
fig,ax = plt.subplots(2,1, figsize=(10,7))
ax[0].plot(z, K[0,:])
ax[1].plot((z[1:] + z[:-1])/2, K[1,:-1])
ax[0].set_xlabel("Altitude")
ax[1].set_xlabel("Altitude")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Wave propagation")
ax[0].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
ax[1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
ax[0].plot(z[index_source],0,'xr')
ax[0].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[1].plot(z[index_source],0,'xr')
ax[1].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[0].grid()
ax[1].grid()
fig.suptitle("Kernels")