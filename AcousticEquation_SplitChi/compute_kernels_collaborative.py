#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
sys.path.insert(1, './Lib/')

from discretisation import *
from linearised_navier_stokes_acoustic_adjoint_splitchi import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
import time





order_collaboration = [["velocity", "pressure"], ["pressure", "velocity"]]
collaboration = order_collaboration[0]

step = 0.01

model_aprior = deepcopy(U0)


# LOAD OBSERVATIONS 
history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)

max_obs_rho = 1
max_obs_v = [ max([abs(history_obs[t].v[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs_p = [ max([abs(history_obs[t].p[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs = [max_obs_rho, max_obs_v, max_obs_p]


for contrib in [0,1]:
    which_chi = collaboration[contrib]

    # GET SOLUTION OF ACOUSTIC EQUATIONS IN REVERSED TIME BASIS
    is_reverse = False
    U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
    t_end, U_end, history_forwardapriori = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", model_aprior, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, is_reverse)
    history_reverse = [history_forwardapriori[t] for t in range(len(history_forwardapriori)-1,-1,-1)]

    # RESOLUTION ADJOINT EQUATIONS                      
    U_start = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
    t_start, U_start, history_adjoint = time_scheme(get_adjoint_RHS, U_start, T_init, Tmax, z, "adjoint", model_aprior, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, history_reverse, history_obs, index_receivers, max_obs, which_chi)

    # COMPUTE KERNELS
    K = get_kernels(history_adjoint, history_reverse, model_aprior, T0, 0, gamma, Cv, dt, h) 

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

    # UPDATE MODEL
    model_aprior.rho = model_aprior.rho + step * K[0,:] 
    #model_aprior.v = model_aprior.v + step * K[1,:] 
    #model_aprior.p = model_aprior.p + step * K[2,:]

    # PLOT MODEL
    model_aprior.plot(z, "Model after "+str(contrib)+"e contribution")

    # SAVE RESULT
    np.save("./BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), np.array(history_reverse))
    np.save("./BackUps/adjoint_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), np.array(history_adjoint))
    np.save("./BackUps/kernel_contri"+which_chi+str(contrib)+"_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), np.array(history_adjoint))
