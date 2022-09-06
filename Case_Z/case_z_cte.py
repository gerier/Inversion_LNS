#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""


import sys
sys.path.insert(1, '../Lib/')

from discretisation import *
from linearised_navier_stokes_testdc import *


import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from copy import deepcopy
import pandas as pd
from scipy.interpolate import CubicSpline
from copy import deepcopy

import sys
import time

# definition of parameters

# source
f0 = 0.1
# space
z0 = 0
zmax = 20e3
# time 
T_init = 0
Tmax = 50 #100*0.005

# user paramter
display_anim = False
case2 = False
time_scheme = EE

#%% MODEL PARAMETERS

rho0 = 1.04898036 
v0 = 0*100 #-4.33914709 #0.00745178154
p0 = 88714.4844 
T0 = 294.375305 
c = 344.108887 
g = 0*9.81084824 
mu = 0*1.27685234e-05
kappa = 0.0258137602 
gamma = 1.40011787 
eta = 0*2.72986326e-05 
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
t_ax = np.arange(T_init,Tmax+dt,dt/2)
source = get_source(t_ax, f0) #scipy.stats.norm.pdf(t_ax,10,1.5)

#source = np.linspace(0, len(source)-1, len(source))


dist_z = abs(z - z[iz_source])
factor_source = np.exp(-dist_z/1000)

# plot the source in time and space
total_source = np.zeros((len(z), len(t_ax)))
for t in t_ax : 
    for i in range(len(z)):
        total_source[i,int(round(2*t/dt))] = source[int(round(2*t/dt))] * factor_source[i]


fig, ax = plt.subplots()
plt.imshow(total_source, aspect='auto')
plt.colorbar(label="Intensity")
plt.title("Form of the source applied on density on point z="+ str(z[iz_source]/1000)+"km")
plt.xticks(range(0,len(t_ax),1000), t_ax[0:len(t_ax):1000])
plt.yticks(range(0,len(z),20), z[0:len(z):20]/1000)
plt.xlabel("Time")
plt.ylabel("Distance (km)")
plt.close()

source = [source, factor_source]
reversed_source = [np.flip(source[0][:-1]), source[1]] 

plt.figure()
plt.plot(t_ax, source[0], label="source")
#plt.plot(np.flip(t_ax), reversed_source[0], label="reverse source")
plt.title("Form of the source applied on density on point z = "+str(z[iz_source]/1000)+"km")
plt.xlabel('Time')
plt.ylabel("Intensity")
plt.show()
plt.close()

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
plt.close()

#%% RESOLUTION

is_reverse = False
# Resolution in a 1d case
start_time = time.time()
t_end, U_end, history_obs = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse)
print(" Total time : ", time.time() - start_time)

#fn1 = get_RHS(history_obs[-iterat], t, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse)
#    fn = get_RHS(history_obs[-iterat-1], t, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse)
#diff = history_obs[-iterat].rho - dt * fn1.rho - history_obs[-iterat-1].rho
#    diff2 = fn1.rho - fn.rho
#    print(diff2)
#    trouve = np.all(np.abs(diff) > 1e-12)
#    iterat += 1  
#print(max(diff))    
#print(trouve, iterat-1, diff)


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
plt.close()

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


if display_anim : 
    get_anim("Density", history_obs[::10], (-0.0006,0.0006), z, 'forward')
    get_anim("Velocity", history_obs[::10], (-0.01,0.02), z, 'forward')
    get_anim("Pressure", history_obs[::10], (-6,6), z, 'forward')


#%% BACKWARD RESOLUTION
    
# Resolution in a 1d case
U_tmax = deepcopy(U_end)
is_reverse = True
start_time = time.time() 
#t_start, U_start, history_reverse = time_scheme(get_minus_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, reversed_source, is_reverse, history_obs)
t_start, U_start, history_reverse = time_scheme(get_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse, history_obs)
print(" Total time : ", time.time() - start_time)

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
plt.close()

if display_anim : 
    get_anim("Density", history_reverse[::10], (-0.02,0.02), z, 'backward')
    get_anim("Velocity", history_reverse[::10], (-4,4), z, 'backward')
    get_anim("Pressure", history_reverse[::10], (-2000,2000), z, 'backward')

#%% SUPERIMPOSITION OF BACKWARD AND FORWARD SOLUTION

n = len(history_reverse)
t_ax = np.append(t_ax, [t_ax[-1]+dt])
t_ax = t_ax[::2]
for i in range(0,n,int(n/4)):
    fig,ax = plt.subplots(3,2, figsize=(10,7))
    # forward
    ax[0,0].plot(z,history_obs[i].rho, 'b', label='Forward')
    ax[1,0].plot(z,history_obs[i].p, 'b')
    ax[2,0].plot( (z[1:] + z[:-1])/2,history_obs[i].v, 'b')
    # backward
    ax[0,0].plot(z,history_reverse[n-1-i].rho, 'r', label='Backward')
    ax[1,0].plot(z,history_reverse[n-1-i].p, 'r')
    ax[2,0].plot( (z[1:] + z[:-1])/2,history_reverse[n-1-i].v, 'r')
    # difference
    ax[0,1].plot(z,abs(history_reverse[n-1-i].rho - history_obs[i].rho), 'k', label='Backward')
    ax[1,1].plot(z,abs(history_reverse[n-1-i].p - history_obs[i].p), 'k')
    ax[2,1].plot( (z[1:] + z[:-1])/2,abs(history_reverse[n-1-i].v - history_obs[i].v), 'k')
    # element to display 
    ax[0,0].legend()
    ax[0,0].set_xlabel("Distance on axis x (km)")
    ax[1,0].set_xlabel("Distance on axis x (km)")
    ax[2,0].set_xlabel("Distance on axis x (km)")
    ax[0,0].set_ylabel("Density")
    ax[1,0].set_ylabel("Pressure")
    ax[2,0].set_ylabel("Velocity")
    ax[0,1].set_xlabel("Distance on axis x (km)")
    ax[1,1].set_xlabel("Distance on axis x (km)")
    ax[2,1].set_xlabel("Distance on axis x (km)")
    ax[0,1].set_ylabel("Difference Density")
    ax[1,1].set_ylabel("Difference Pressure")
    ax[2,1].set_ylabel("Difference Velocity")
    ax[0,1].grid()
    ax[1,1].grid()
    ax[2,0].grid()
    fig.suptitle("Comparison of forward and backward solution at time (forward) t = "+str(t_ax[i]))
    plt.show()
    
 
