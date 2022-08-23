#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

from discretisation import *
from linearised_navier_stokes_testdc import *
#from linearised_navier_stokes_acoustic_adjoint import *


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
zmax = 70e3
# time 
T_init = 0
Tmax = 50 # *0.005 #100*0.005

# user paramter
display_anim = False
case2 = False
time_scheme = EE

#%% MODEL PARAMETERS (A PRIORI MODEL)

rho0 = 1.04898036 
v0 = 0*100 #-4.33914709 #0.00745178154
p0 = 88714.4844 
T0 = 294.375305 
c = 344.108887 
g = 0*9.81084824 
mu = 0*1.27685234e-05
kappa = 0*0.0258137602 
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
index_receivers = [250] #np.arange(z0+d_receivers, zmax, d_receivers)

# define Time
t = 0
dt = 0.005/10

#%% DEFINE A PRIORI MODEL ON MESH

for i,param in enumerate(param_toset_onmesh):
    param_toset_onmesh[i] = param * np.ones(len(z))

[rho0, v0, p0, T0, c, g, mu, kappa, gamma, eta, Cv] = param_toset_onmesh

## DEFINE THE "REAL" MODEL 

[GT_rho0, GT_v0, GT_p0, GT_T0,GT_c, GT_g, GT_mu, GT_kappa, GT_gamma, GT_eta, GT_Cv] = deepcopy(param_toset_onmesh)

factor = 10
if True : # False : 
    for param in GT_v0, GT_rho0:#, GT_p0:
        param[80:120] = factor * param[80:120]
        #param[80:90] = [  (param[90]*factor-param[80]) / 10 * i - 8*factor*param[90]+9*param[80]  for i in range(80,90)]
        #param[110:120] = [  -(param[110]*factor-param[120]) / 10 * i + 12*factor*param[110]-11*param[120]  for i in range(110,120)]
        #param[90:110] *= factor
        #param[120:160] = [  (param[159]-1.1*param[120]) / 40 * i + 4.4*param[120]-3*param[159]  for i in range(120,160)]
# param[80:120]*2 #                        


#%% DEFINE SOURCE ON MESH

# set the source
iz_source =  200#30
t_ax = np.arange(T_init,Tmax+dt,dt/2)
source = get_source(t_ax, f0) #scipy.stats.norm.pdf(t_ax,10,1.5)

dist_z = abs(z - z[iz_source])
dist_z = np.roll(dist_z,-170)
factor_source = np.exp(-2*dist_z/1000)
 

#source = 10*np.ones(len(source)) #np.linspace(0, len(source)-1, len(source))
#factor_source = np.zeros(len(z))
#factor_source[iz_source] = 1

if False:
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
    #plt.close()

source = [source, factor_source]
#reversed_source = [np.flip(source[0][:-1]), source[1]] 

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
GT_v0_demi = (GT_v0[1:] + GT_v0[:-1])/2 
U0 = LNS_Variable(rho0, v0_demi, p0) 
GT_U0 = LNS_Variable(GT_rho0, GT_v0_demi, GT_p0) 
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
fig.suptitle("The  a priori background")



# Plot the model
fig,ax = plt.subplots(3,1, figsize=(10,7))
ax[0].plot(z,GT_U0.rho)
ax[1].plot(z,GT_U0.p)
ax[2].plot( (z[1:] + z[:-1])/2,GT_U0.v)
ax[0].set_xlabel("Distance on axis x (km)")
ax[1].set_xlabel("Distance on axis x (km)")
ax[2].set_xlabel("Distance on axis x (km)")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("The ground truth background")


#%% RESOLUTION

is_reverse = False
# Resolution in a 1d case
start_time = time.time()
t_end, U_end, history_obs = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", GT_U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse)
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


a = history_obs[-1]

#%% BACKWARD RESOLUTION
    
# Resolution in a 1d case
 
 # first do the forward pb with the a priori model
is_reverse = False
# Resolution in a 1d case
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
t_end, U_end, history_forwardapriori = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse)


# make the checkpointing op
U_tmax = deepcopy(U_end)
is_reverse = True
start_time = time.time() 
#t_start, U_start, history_reverse = time_scheme(get_minus_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, reversed_source, is_reverse, history_obs)
#t_start, U_start, history_reverse = time_scheme(get_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse, history_forwardapriori)
history_reverse = [history_forwardapriori[t] for t in range(len(history_obs)-1,-1,-1)]
U_start = history_forwardapriori[-1]
print(" Total time : ", time.time() - start_time)



#%% RESULTS BACKWARD

fig,ax = plt.subplots(3,1, figsize=(10,7))
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
    # color where the model is different 
    ax[0,1].axvspan(z[80], z[120], alpha=0.1, color='grey')
    ax[1,1].axvspan(z[80], z[120], alpha=0.1, color='grey')
    ax[2,1].axvspan(z[80], z[120], alpha=0.1, color='grey')
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
    
 
#%% RESOLUTION ADJOINT EQUATIONS

                                   
max_obs_rho = 1
max_obs_v = [ max([history_obs[t].v[k] for t in range(len(history_obs))]) for k in index_receivers]
max_obs_p = [ max([history_obs[t].p[k] for t in range(len(history_obs))]) for k in index_receivers]

max_obs = [max_obs_rho, max_obs_v, max_obs_v]

U_tmax = deepcopy(U_end)
t_start, U_start, history_adjoint = time_scheme(get_adjoint_RHS, U_end*0, Tmax, T_init,  z, "adjoint", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, history_reverse, history_obs, index_receivers, max_obs)

#%%
# Plot the result
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
ax[0].vlines([z[80], z[120]], *ax[0].get_ylim(), 'r')
ax[1].vlines([z[80], z[120]], *ax[1].get_ylim(), 'r')
ax[2].vlines([z[80], z[120]], *ax[2].get_ylim(), 'r')
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Result of adjoint equations (Tmax to T0) ")


#%% KERNELS

K = get_kernels(history_adjoint, history_reverse, U0, T0, kappa, gamma, Cv, dt, h) 

#%%

# Plot the result
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
ax[0].vlines([z[80], z[120]], *ax[0].get_ylim(), 'r')
ax[1].vlines([z[80], z[120]], *ax[1].get_ylim(), 'r')
ax[2].vlines([z[80], z[120]], *ax[2].get_ylim(), 'r')
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Kernels")

#%% Attendu 

which_param = 'density'

# initialisation of dx
dx_rho = np.zeros(len(history_obs[0].rho))
dx_v = np.zeros(len(history_obs[0].v))
dx_p = np.zeros(len(history_obs[0].p))

# get norm 2
for t in range(len(history_obs)):
    dx_rho += np.power(history_obs[t].rho - history_forwardapriori[t].rho,2)
    dx_v += np.power(history_obs[t].v - history_forwardapriori[t].v,2)
    dx_p += np.power(history_obs[t].p - history_forwardapriori[t].p,2)

# get inf norm to normalization
max_obs_rho = 1
max_obs_v = [ max([history_obs[t].v[k] for t in range(len(history_obs))]) for k in range(len(history_obs[0].v))]
max_obs_p = 1#[ max([history_obs[t].p[k] for t in range(len(history_obs))]) for k in range(len(history_obs[0].p))]

dx = dx_rho / max_obs_rho + dx_p / max_obs_p + interpolation(dx_v / max_obs_v)   

if which_param == 'density': 
    dm = U0.rho - GT_U0.rho
elif which_param == 'velocity': 
    dm = U0.v - GT_U0.v
elif which_param == 'pressure':
    dm = U0.p - GT_U0.p

index_dm = dm != 0

plt.plot(z[index_dm], dx[index_dm]/dm[index_dm])

