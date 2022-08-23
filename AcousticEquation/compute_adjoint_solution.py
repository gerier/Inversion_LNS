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



# LOAD OBSERVATIONS AND SOLUTION WITH A PRIORI MODEL (in reversed time basis)
history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_reverse = np.load("./BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)


# PLOT THE PERTURBATION SOURCE OF THE ADJOINT EQUATIONS
recepteur = np.zeros((len(history_obs), len(index_receivers)))

max_obs_rho = 1
max_obs_v = [ max([abs(history_obs[t].v[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs_p = [ max([abs(history_obs[t].p[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs = [max_obs_rho, max_obs_v, max_obs_p]

plt.figure()
for t in range(len(history_reverse)):
    recepteur[t,:] = dchi(history_reverse, history_obs, t*dt, dt, index_receivers, max_obs)[2,index_receivers]
plt.plot(recepteur[:,0])
plt.grid()
plt.show()


# RESOLUTION ADJOINT EQUATIONS                      
U_start = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
t_start, U_start, history_adjoint = time_scheme(get_adjoint_RHS, U_start, T_init, Tmax, z, "adjoint", U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, history_reverse, history_obs, index_receivers, max_obs)


# SAVE RESULT
np.save("./BackUps/adjoint_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), np.array(history_adjoint))


# PLOT RESULT
plot_adjoint_final_state = False
if plot_adjoint_final_state : 
    fig,ax = plt.subplots(3,1, figsize=(10,7))
    ax[0].plot(z,U_start.rho)
    ax[1].plot(z,U_start.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_start.v)
    ax[0].set_xlabel("Altitude")
    ax[1].set_xlabel("Altitude")
    ax[2].set_xlabel("Altitude")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[0].plot(z[index_receivers],np.zeros(len(index_receivers)), 'xr')
    ax[1].plot(z[index_receivers],np.zeros(len(index_receivers)), 'xr')
    ax[2].plot(z[index_receivers],np.zeros(len(index_receivers)), 'xr')
    ax[0].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax[1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax[2].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax[0].plot(z[30],0,'xr')
    ax[0].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
    ax[1].plot(z[30],0,'xr')
    ax[1].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
    ax[2].plot(z[30],0,'xr')
    ax[2].plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    fig.suptitle("Result of adjoint equations (Tmax to T0) ")

