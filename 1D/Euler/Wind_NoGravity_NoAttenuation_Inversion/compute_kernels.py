#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys

sys.path.insert(1, '../Lib/')

from discretisation import *
from linearised_navier_stokes import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
import time



# LOAD OBSERVATIONS AND SOLUTION WITH A PRIORI MODEL (in reversed time basis)
history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+".npy", allow_pickle=True)
history_reverse = np.load("./BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+".npy", allow_pickle=True)
history_adjoint = np.load("./BackUps/adjoint_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+".npy", allow_pickle=True)


# COMPUTE KERNELS
K = get_kernels(history_adjoint, history_reverse, M, dt, h) 


# PLOT
n = int(len(z) / 2 )
fig,ax = plt.subplots(4,1, figsize=(10,7))
ax[0].plot(z[:n],K[0,:n])
ax[1].plot(z[:n], K[1,:n])
ax[2].plot(z[:n],K[2,:n])
ax[3].plot(z[:n], K[3,:n])
ax[0].set_xlabel("Altitude")
ax[1].set_xlabel("Altitude")
ax[2].set_xlabel("Altitude")
ax[3].set_xlabel("Altitude")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Wind")
ax[2].set_ylabel("Pressure")
ax[3].set_ylabel("Gamma")
ax[0].plot(z[index_source],0,'xr')
ax[0].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[1].plot(z[index_source],0,'xr')
ax[1].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[2].plot(z[index_source],0,'xr')
ax[2].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[3].plot(z[index_source],0,'xr')
ax[3].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
ax[0].grid()
ax[1].grid()
ax[2].grid()
ax[3].grid()
fig.suptitle("Kernels")

# SAVE RESULT
np.save("./BackUps/kernel_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_", K)

np.save("./BackUps/python_version", K[:,50:250])