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





# LOAD OBSERVATIONS AND SOLUTION WITH A PRIORI MODEL (in reversed time basis)
history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_reverse = np.load("./BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_adjoint = np.load("./BackUps/adjoint_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)


# COMPUTE KERNELS
K = get_kernels(history_adjoint, history_reverse, U0, T0, 0, gamma, Cv, dt, h) 


# PLOT
fig,ax = plt.subplots(2,1, figsize=(10,7))
ax[0].plot(z,K[0,:])
ax[1].plot( (z[1:] + z[:-1])/2,K[1,:-1])
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

# SAVE RESULT
np.save("./BackUps/kernel_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K)
