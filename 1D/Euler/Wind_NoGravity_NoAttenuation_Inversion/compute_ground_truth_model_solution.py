#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""


from parameters import *
from discretisation import *
from linearised_navier_stokes import *

import time
import numpy as np

import cProfile
cp = cProfile.Profile()
cp.enable()

# INITIALISATION
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z)))


# RESOLUTION
start_time = time.time()
t_end, U_end, history_obs = time_scheme(
    get_RHS, U1, T_init, Tmax, z, "forward", GT_M, h, dt, source)
print(" Total time : ", time.time() - start_time)

np.save("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_" +
        str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_", np.array(history_obs))

np.save("./BackUps/python_version", np.array(history_obs))
np.save("./BackUps/python_version", history_obs[-1].p[50:250])

cp.disable()


# juste les observations
# comparison
n = len(history_obs)
t_ax = np.arange(T_init, Tmax+dt, dt/2)
t_ax = np.append(t_ax, [t_ax[-1]+dt])
t_ax = t_ax[::2]
for i in range(0, n, int(n/4)):
    fig, ax = plt.subplots(3, 1, figsize=(10, 7))
    # forward
    l = int(len(z)/2)
    ax[0].plot(z[:l]/1000, history_obs[i].rho[:l], 'b', label='Modèle réel')
    ax[1].plot(z[:l]/1000, history_obs[i].p[:l], 'b')
    ax[2].plot((z[1:l] + z[:l-1])/2/1000, history_obs[i].v[:l-1], 'b')

    # element to display
    ax[0].legend()
    ax[2].set_xlabel("Distance on axis x (km)")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")

    ax[0].grid()
    ax[1].grid()
    ax[2].set_ylabel("Velocity")
    ax[2].grid()
    fig.suptitle(
        "Comparison of forward and backward solution at time (forward) t = "+str(t_ax[i]))
    plt.show()