#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

from linearised_navier_stokes import *
from discretisation import *

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from copy import deepcopy
import pandas as pd
from scipy.interpolate import CubicSpline
from copy import deepcopy



# definition of parameters

# source
f0 = 0.1
# space
z0 = 0
zmax = 20e3
# time 
T_init = 0
Tmax = 20

# user paramter
display_anim = False
case2 = False

#%% MODEL PARAMETERS

rho0 = 1.04898036 
v0 = -4.33914709 #0.00745178154
p0 =88714.4844 
T0 = 294.375305 
c = 344.108887 
g = 0*9.81084824 
mu = 1.27685234e-05
kappa = 0.0258137602 
gamma = 1.40011787 
eta = 2.72986326e-05 
cv = 20.7801247 
M = 28.965 # masse molaire de l'air https://fr.wikipedia.org/wiki/Air

Cv = cv/M
l = eta - (2/3)*mu

param_toset_onmesh = [rho0, v0, p0, T0, c, g, mu, kappa, gamma, eta, Cv]


#%% CONSIDERATION ON MESH

min_v = abs(c) #, abs(v0))
max_v = max(abs(c), abs(v0))
max_dx = min_v / (10 * 2.5 * f0) 

max_dt = max_dx / max_v

print("With this source, you must choose a dx < ", str(max_dx) )

print("If taking the dx max, finaly, you will have to choose a dt < ", str(max_dt))


#%% MESH PARAMETERS

# define the space
h = 100      # size of the mesh
z = np.arange(z0,zmax+h,h)


# define receivers
d_receivers = 3.6e3 # distance between 2 receivers, 
                  # must be a multiple of h
r = np.arange(z0+d_receivers, zmax, d_receivers)

# define Time
t = 0
dt = 0.005

#%% DEFINE MODEL ON MESH

for i,param in enumerate(param_toset_onmesh):
    param_toset_onmesh[i] = param * np.ones(len(z))

[rho0, v0, p0, T0, c, g, mu, kappa, gamma, eta, Cv] = param_toset_onmesh


#%% DEFINE SOURCE ON MESH

# set the source
iz_source = 100
t_ax = np.arange(T_init,Tmax+dt+dt,dt)
source = get_source(t_ax, f0) #scipy.stats.norm.pdf(t_ax,10,1.5)

dist_z = abs(z - z[iz_source])
factor_source = np.exp(-dist_z/1000)

plt.plot(t_ax, source)
plt.title("Form of the source applied on density on point z = "+str(z[iz_source]/1000)+"km")
plt.xlabel('Time')
plt.ylabel("Intensity")
plt.show()

# plot the source in time and space
total_source = np.zeros((len(z), len(t_ax)))
for t in t_ax : 
    for i in range(len(z)):
        total_source[i,int(t/dt)] = source[int(t/dt)] * factor_source[i]


fig, ax = plt.subplots()
plt.imshow(total_source, aspect='auto')
plt.colorbar(label="Intensity")
plt.title("Form of the source applied on density on point z="+ str(z[iz_source]/1000)+"km")
plt.xticks(range(0,len(t_ax),1000), t_ax[0:len(t_ax):1000])
plt.yticks(range(0,len(z),20), z[0:len(z):20]/1000)
plt.xlabel("Time")
plt.ylabel("Distance (km)")

source = [source, factor_source]
reversed_source = [np.flip(source[0])[1:], source[1]]
#%% INTIALISATION

# define vectors of the system
v0_demi = (v0[1:] + v0[:-1])/2
U0 = LNS_Variable(rho0, v0_demi, p0) 
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 

#%% 

# Plot the model
fig,ax = plt.subplots(3,1, figsize=(10,7))
ax[0].plot(z,U0.rho)
ax[1].plot(z,U0.p)
ax[2].plot( (z[1:] + z[:-1])/2,U0.v)
ax[0].set_xlabel("Distance on axis x (km)")
ax[1].set_xlabel("Distance on axis x (km)")
ax[2].set_xlabel("Distance on axis x (km)")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("The background")


#%% RESOLUTION

is_reverse = False
# Resolution in a 1d case
t_end, U_end, history_obs = EE(get_RHS, U1, T_init, Tmax, z, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse)


#%% RESULTS

fig,ax = plt.subplots(3,1, figsize=(10,7))
fig
ax[0].plot(z,U_end.rho)
ax[1].plot(z,U_end.p)

ax[2].plot( (z[1:] + z[:-1])/2,U_end.v)
ax[0].set_xlabel("Distance on axis x (km)")
ax[1].set_xlabel("Distance on axis x (km)")
ax[2].set_xlabel("Distance on axis x (km)")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Perturbation after Tmax = "+ str(Tmax))

#get_ipython().run_line_magic('matplotlib', 'inline')
from IPython.display import HTML
from matplotlib import animation, rc

def get_anim(parameter_toplot, history, limits, z, propagation) : 
    list_anim = []

    fig,ax = plt.subplots(figsize=(10,3))
    ax.set_xlabel("Altitude")
    ax.set_ylabel(parameter_toplot)

    ax.grid()
    ax.set_xlim(( z[0], z[-1]))
    ax.set_ylim(limits)

    line0, = ax.plot([], [], lw=2)

    # initialization function: plot the background of each frame
    def init():
        line0.set_data([], [])
        return (line0,)

    # animation function. This is called sequentially
    def animate(i,z,history, parameter_toplot):
        x = z
        if parameter_toplot == "Density" :
            y = history[i].rho
        elif parameter_toplot == "Velocity" :
            y = history[i].v  
            x = (z[1:] + z[:-1])/2
        elif parameter_toplot == "Pressure" :
            y = history[i].p 
        line0.set_data(x, y)
        return (line0,)

    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(z,history,parameter_toplot), frames=len(history), interval=20, blit=True)
    #writervideo = animation.FFMpegWriter(fps=60)
    anim.save('./animation_'+parameter_toplot+'_'+propagation+".mp4")#, writer=writervideo)
    plt.close()

disp_anim = True
if disp_anim : 
    get_anim("Density", history_obs[::10], (-0.02,0.02), z, 'forward')
    get_anim("Velocity", history_obs[::10], (-4,4), z, 'forward')
    get_anim("Pressure", history_obs[::10], (-2000,2000), z, 'forward')

#%% BACKWARD RESOLUTION
    
# Resolution in a 1d case
U_tmax = deepcopy(U_end)
t_start, U_start, history_reverse = EE(get_minus_RHS, U_end, T_init, Tmax, z, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, reversed_source, is_reverse)


#%% RESULTS BACKWARD

fig,ax = plt.subplots(3,1, figsize=(10,7)
                     )
ax[0].plot(z,U_start.rho)
ax[1].plot(z,U_start.p)
ax[2].plot( (z[1:] + z[:-1])/2,U_start.v)
ax[0].set_xlabel("Distance on axis x (km)")
ax[1].set_xlabel("Distance on axis x (km)")
ax[2].set_xlabel("Distance on axis x (km)")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("After Tmax = 0")


get_anim("Density", history_reverse[::10], (-0.02,0.02), z, 'backward')
get_anim("Velocity", history_reverse[::10], (-4,4), z, 'backward')
get_anim("Pressure", history_reverse[::10], (-2000,2000), z, 'backward')

#%% SUPERIMPOSITION OF BACKWARD AND FORWARD SOLUTION

n = len(history_reverse)
for i in range(0,n,int(n/5)):
    fig,ax = plt.subplots(3,1, figsize=(10,7))
    # forward
    ax[0].plot(z,history_obs[i].rho, 'b', label='Forward')
    ax[1].plot(z,history_obs[i].p, 'b')
    ax[2].plot( (z[1:] + z[:-1])/2,history_obs[i].v, 'b')
    # backward
    ax[0].plot(z,history_reverse[n-1-i].rho, 'r', label='Backward')
    ax[1].plot(z,history_reverse[n-1-i].p, 'r')
    ax[2].plot( (z[1:] + z[:-1])/2,history_reverse[n-1-i].v, 'r')
    # element to display 
    ax[0].legend()
    ax[0].set_xlabel("Distance on axis x (km)")
    ax[1].set_xlabel("Distance on axis x (km)")
    ax[2].set_xlabel("Distance on axis x (km)")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    fig.suptitle("Comparison of forward and backward solution at time (forward) t = "+str(t_ax[i]))
    plt.show()
    
    
