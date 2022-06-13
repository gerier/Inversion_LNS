#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 12:30:00 2022

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


# Definition of temporal scheme - Euler explicit
def EE(f, U_t, T_init, Tmax, z, *args):
    fig,ax = plt.subplots(3,1)
    ax[0].plot(z,U_t.rho, label=str(T_init))
    ax[1].plot(z,U_t.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    # load the dt
    dt = args[3]

    # create a vector to save state at each time
    history = [deepcopy(U_t)]
    print()
    for t in np.arange(T_init,Tmax,dt):
        U_t = U_t +  f(U_t, t, *args) * dt
        #U_t = apply_sponge(U_t)
        if (t+dt) % 5 < 1e-4 : 
            ax[0].plot(z,U_t.rho, label='t = '+str(t+dt))
            ax[1].plot(z,U_t.p)
            ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)

        history += [deepcopy(U_t)]
        
    # make a plot to vizualise the transformation of the threee quantities    
    ax[0].set_xlabel("Altitude")
    ax[1].set_xlabel("Altitude")
    ax[2].set_xlabel("Altitude")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    fig.legend()
    return t, U_t, np.array(history)


# sponge layer
def apply_sponge(U, index_z0 = 20, alpha = 13):
    # apply a sponge layer for z < z0 and z > zmax - z0
    for i in range(len(U.rho)):
        """
        if i < index_z0 +1:
            U.rho[i] *=  np.exp( - alpha * abs( i - index_z0) / index_z0)
            U.p[i] *=  np.exp( - alpha * abs( i - index_z0) / index_z0)
            # we have to consider that v is define on the center of each cell
            if (i < len(U.rho)-1) and (i < index_z0 - 1) :
                U.v[i] *=  np.exp( - alpha * abs( i - index_z0 + 1/2) / index_z0)
        elif i > len(U.rho) - index_z0 : 
            U.rho[i] *=  np.exp( - alpha * abs( i - len(U.rho) + index_z0) / index_z0)
            U.p[i] *=  np.exp( - alpha * abs( i - len(U.rho) + index_z0) / index_z0)            
            # we have to consider that v is define on the center of each cell
            if (i < len(U.rho)-1) and (i > len(U.rho) - index_z0 + 1) :
                U.v[i] *=  np.exp( - alpha * abs( i - len(U.rho) + index_z0 + 1/2) / index_z0)    
        """
        if i < index_z0 +1:
            sponge = 1 - (1- np.exp(alpha *  (( i - index_z0)/index_z0)**2))/(1-np.exp(alpha))
            if i < len(U.rho)-1:
                U.v[i] *= sponge
            U.p[i] *= sponge
        elif i > len(U.rho) - index_z0 :
            sponge = 1 - (1- np.exp(alpha *  (( i - len(U.rho) + index_z0 + 1)/index_z0)**2))/(1-np.exp(alpha))
            if i < len(U.v):
                U.v[i] *= 1 - (1- np.exp(alpha *  (( i - len(U.v) + index_z0 + 1)/index_z0)**2))/(1-np.exp(alpha)) 
            U.p[i] *= sponge          
    return U

def F(U_t, it, source, n, gamma):
    # define source term
    f = np.zeros((3,n))
    source_time = source[0]
    source_spatial = source[1]
    # compute the contribution on density
    for i in range(len(f[0])):
        f[0,i] = source_time[it]/100 * source_spatial[i]
    # compute the contribution on quantity of mvt
    f[1,:-1] =  U_t.v * (f[0,1:]+f[0,:-1])/2 
    f[2,:] =   gamma * (U_t.p ) / U_t.rho * f[0,:] 
    #f[0,:] = 0 # todo change correctly , here because the source is only on velocity
    return LNS_Variable(0*f[0,:], f[1,:-1], f[2,:])

# Definition of resolution spatial
def get_acoustic_RHS(U_t, t, l, gamma, dz, dt, source, BC, order_DF=4):
    if order_DF == 4:
        DF_DC_backward = DF_DC_4_backward
        DF_DC_forward = DF_DC_4_forward
        DF_C = DF_C_4
    else :
        DF_DC_backward = DF_DC_2_backward
        DF_DC_forward = DF_DC_2_forward
        DF_C = DF_C_2
        
    demi_rho = (U_t.rho[1:] + U_t.rho[:-1]) / 2
    U_t1  = np.zeros((3,len(U_t.rho)))
    
    BC_v = BC[0] * np.ones(4)
    BC_p = BC[1] * np.ones(4)
    
    source = F(U_t, int(t/dt), source, len(U_t1[0,:]), gamma)
    
    U_t1[0,:] = 0
    U_t1[1,:-1] = -DF_C(U_t.p, dz, BC_p, False)
    U_t1[1,:-1] += source.v 
    U_t1[2,:] = - l * DF_C(U_t.v,dz,BC_v)
    U_t1[2,:] += source.p
    
    return LNS_Variable(U_t1[0,:], demi_rho * U_t1[1,:-1] , U_t1[2,:])

#%% Test case

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

# 
rho = 1.04898036 
c = 344.108887
mu = 1.27685234e-05
eta = 2.72986326e-05 
gamma =  1.40011787 

v = -4.33914709 
p =88714.4844 

l = eta - (2/3)*mu

param_toset_onmesh = [rho, v, p, l, gamma]

#%% Mesh parametrisation
min_v = abs(c) #, abs(v0))
max_v = max(abs(c), abs(v))
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

[rho, v, p, l, gamma] = param_toset_onmesh

#%% SOURCE

iz_source = 100
t_ax = np.arange(T_init,Tmax+dt+dt,dt)
source = get_source(t_ax, f0) #scipy.stats.norm.pdf(t_ax,10,1.5)

dist_z = abs(z - z[iz_source])
factor_source = np.exp(-dist_z/1000)

plt.plot(t_ax, source)
plt.title("Form of the source applied on density on point z=12km")
plt.show()

# plot the source in time and space
total_source = np.zeros((len(z), len(t_ax)))
for t in t_ax : 
    for i in range(len(z)):
        total_source[i,int(t/dt)] = source[int(t/dt)] * factor_source[i]


fig, ax = plt.subplots()
plt.imshow(total_source, aspect='auto')
plt.colorbar()
plt.title("Form of the source applied on density on point z="+ str(z[iz_source])+"km")
plt.xticks(range(0,len(t_ax),1000), t_ax[0:len(t_ax):1000])
plt.yticks(range(0,len(z),20), z[0:len(z):20]/1000)
plt.xlabel("Time")
plt.ylabel("Altitude")

source = [source, factor_source]
reversed_source = [np.flip(source[0])[1:], source[1]]


#%% INTIALISATION

# define vectors of the system
v_demi = (v[1:] + v[:-1])/2
U_0 = LNS_Variable(rho, v_demi, p)

 
# Plot the model
fig,ax = plt.subplots(3,1, figsize=(10,7))
ax[0].plot(z,U_0.rho)
ax[1].plot(z,U_0.p)
ax[2].plot( (z[1:] + z[:-1])/2,U_0.v)
ax[0].set_xlabel("Distance on axis x (km)")
ax[1].set_xlabel("Distance on axis x (km)")
ax[2].set_xlabel("Distance on axis x (km)")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Initial state")


#%% RESOLUTION
U_t = deepcopy(U_0)
t_end, U_end, history_obs = EE(get_acoustic_RHS, U_t, T_init, Tmax, z, l, gamma, h, dt, source, [v[0],p[0]])

