#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys

sys.path.insert(1, '../Lib/')

from discretisation import *
from linearised_navier_stokes_acoustic_adjoint import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
import time


def get_kernels_centered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dtvp, dz, DF_C):
    # initialization
    K = np.zeros((3, len(rho_a)))
    
    # main boundary conditions
    BC_vp = [v_p[-2], v_p[-1], v_p[0], v_p[1]]

    # kernel in rho0
    K[0,:] -= interpolation(dtvp * v_a)    
    K[0,:] += DF_C(v_p, dz, BC_vp) * p_a / U0.rho
    
    return K
    


def get_kernels(hist_adjoint, hist_backprop, U0, T0, kappa, gamma, Cv, dt, dz, order_DF=4):
    # initialize
    K = np.zeros((3, len(hist_adjoint[0].rho)))
    if order_DF == 4:
        DF_C = DF_C_4
    else :
        DF_C = DF_C_2

    for t in range(len(hist_adjoint)-1):
        # create variable to be more lisible 
        rho_a = (hist_adjoint[t].rho + hist_adjoint[t+1].rho)/2    # mean at time n+1/2
        rho_p = (hist_backprop[t].rho + hist_backprop[t+1].rho)/2  # mean at time n+1/2
        
        p_a = (hist_adjoint[t].p + hist_adjoint[t+1].p)/2          # mean at time n+1/2
        p_p = (hist_backprop[t].p + hist_backprop[t + 1].p)/2         # mean at time n+1/2 
        
        v_a = hist_adjoint[t].v                                   # at time n
        v_p = hist_backprop[t].v                                  # at time n
        dtvp = (hist_backprop[t +1].v - v_p) / dt
        
        K += dt * get_kernels_centered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dtvp, dz, DF_C)
       
    return K


# LOAD OBSERVATIONS AND SOLUTION WITH A PRIORI MODEL (in reversed time basis)
history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_reverse = np.load("./BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_adjoint = np.load("./BackUps/adjoint_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)


# COMPUTE KERNELS
K = get_kernels(history_adjoint, history_reverse, U0, T0, 0, gamma, Cv, dt, h) 


# PLOT
fig,ax = plt.subplots(3,1, figsize=(10,7))
ax[0].plot(z,K[0,:])
ax[1].plot(z,K[2,:])
ax[2].plot( (z[1:] + z[:-1])/2,K[1,:-1])
ax[0].set_xlabel("Altitude")
ax[1].set_xlabel("Altitude")
ax[2].set_xlabel("Altitude")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
ax[1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
ax[2].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
ax[0].plot(z[index_source],0,'xr')
ax[0].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[1].plot(z[index_source],0,'xr')
ax[1].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[2].plot(z[index_source],0,'xr')
ax[2].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Kernels")

# SAVE RESULT
np.save("./BackUps/kernel_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K)
