#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""


from parameters import *
from discretisation import *
from linearised_navier_stokes import *


import numpy as np
import time


# INITIALISATION
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z)))

# CHOICE
which_method = "forward"

# %% BACKWARD RESOLUTION
if which_method == "forward":
    is_reverse = False
    t_end, U_end, history_forwardapriori = time_scheme(
        get_RHS, U1, T_init, Tmax, z, "forward", M, h, dt, source)
    history_reverse = history_forwardapriori[::-1]#[history_forwardapriori[t]
                      # for t in range(len(history_forwardapriori)-1, -1, -1)]
    U_start = history_reverse[-1]

# SAVE

np.save("./BackUps/python_version", np.array(history_reverse[0].v[50:250]))

np.save("./BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax) +
        "_"+str(dt)+"_"+str(z[index_source])+"_", np.array(history_reverse))


# %% RESULTS BACKWARD

U_start.plot(z, "After Tmax = 0")

if display_anim:
    get_anim("Density", history_reverse[::10], (-0.02, 0.02), z, 'backward')
    get_anim("Velocity", history_reverse[::10], (-4, 4), z, 'backward')
    get_anim("Pressure", history_reverse[::10], (-2000, 2000), z, 'backward')

# %% SUPERIMPOSITION OF BACKWARD AND FORWARD SOLUTION

# load observations
history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h) +
                      "_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+".npy", allow_pickle=True)

# comparison
n = len(history_reverse)
t_ax = np.append(t_ax, [t_ax[-1]+dt])
t_ax = t_ax[::2]
for i in range(0, n, int(n/4)):
    fig, ax = plt.subplots(3, 2, figsize=(10, 7))
    # forward
    ax[0, 0].plot(z/1000, history_obs[i].rho, 'b', label='Modèle réel')
    ax[1, 0].plot(z/1000, history_obs[i].p, 'b')
    ax[2, 0].plot((z[1:] + z[:-1])/2/1000, history_obs[i].v, 'b')
    # backward
    ax[0, 0].plot(z/1000, history_reverse[n-1-i].rho,
                  'r', label='Modèle a priori')
    ax[1, 0].plot(z/1000, history_reverse[n-1-i].p, 'r')
    ax[2, 0].plot((z[1:] + z[:-1])/2/1000, history_reverse[n-1-i].v, 'r')
    # difference
    ax[0, 1].plot(z/1000, abs(history_reverse[n-1-i].rho -
                  history_obs[i].rho), 'k', label='Backward')
    ax[1, 1].plot(
        z/1000, abs(history_reverse[n-1-i].p - history_obs[i].p), 'k')
    ax[2, 1].plot((z[1:] + z[:-1])/2/1000,
                  abs(history_reverse[n-1-i].v - history_obs[i].v), 'k')
    # element to display
    ax[0, 0].legend()
    ax[0, 0].set_xlabel("Distance on axis x (km)")
    ax[1, 0].set_xlabel("Distance on axis x (km)")
    ax[2, 0].set_xlabel("Distance on axis x (km)")
    ax[0, 0].set_ylabel("Density")
    ax[1, 0].set_ylabel("Pressure")
    ax[2, 0].set_ylabel("Velocity")
    ax[0, 1].set_xlabel("Distance on axis x (km)")
    ax[1, 1].set_xlabel("Distance on axis x (km)")
    ax[2, 1].set_xlabel("Distance on axis x (km)")
    ax[0, 1].set_ylabel("Difference Density")
    ax[1, 1].set_ylabel("Difference Pressure")
    ax[2, 1].set_ylabel("Difference Velocity")
    ax[0, 1].grid()
    ax[1, 1].grid()
    ax[2, 0].grid()
    fig.suptitle(
        "Comparison of forward and backward solution at time (forward) t = "+str(t_ax[i]))
    plt.show()

# %%
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
# %%
a = max(abs(history_obs[-1].rho[100:200]))
b = max(abs(history_obs[-1].rho[250:300]))
c = max(abs(history_obs[-1].rho[0:100]))

print(a)
print(max(b, c))
print(a/max(b, c))
