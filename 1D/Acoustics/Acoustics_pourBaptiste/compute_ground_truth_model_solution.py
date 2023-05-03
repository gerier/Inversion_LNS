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
U1 = LNS_Variable(np.zeros(len(z)-1), np.zeros(len(z)))


# RESOLUTION
start_time = time.time()
t_end, U_end, history_obs = time_scheme(
    get_RHS, U1, T_init, Tmax, z, "forward", GT_M, h, dt, source)
print(" Total time : ", time.time() - start_time)

np.save("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_" +
        str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_", np.array(history_obs))

cp.disable()


# plot les observations
# comparison
n = len(history_obs)
t_ax = np.arange(T_init, Tmax+2*dt, dt/2)
t_ax = t_ax[::2]
for i in range(0, n, int(n/4)):
    fig, ax = plt.subplots(2, 1, figsize=(10, 7))
    # forward
    l = int(len(z))
    ax[0].plot(z[:l]/1000, history_obs[i].p[:l], 'b', label='Modèle réel')
    ax[1].plot((z[1:l] + z[:l-1])/2/1000, history_obs[i].v[:l-1], 'b')

    # element to display
    ax[0].legend()
    ax[1].set_xlabel("Distance on axis x (km)")
    ax[0].set_ylabel("Pressure")

    ax[1].set_ylabel("Velocity")
    ax[0].grid()
    ax[1].grid()
    fig.suptitle(
        "Comparison of forward and backward solution at time (forward) t = "+str(t_ax[i]))
    plt.show()
