#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
sys.path.insert(1, '../Lib/')

from discretisation import *
#from linearised_navier_stokes_testdc import *
from linearised_navier_stokes_acoustic_adjoint import *


import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy



import time



use_backups = False
# definition of parameters

# source
f0 = 0.1
# space
z0 = -25e3
nb_index_neg = int(- z0 / 100)

zmax = 45e3
# time 
T_init = 0
Tmax = 50

# user paramter
display_anim = False
case2 = False
time_scheme = EE

#%% MODEL PARAMETERS (A PRIORI MODEL)

rho0 = 1.04898036 
v0 = 0
p0 = 88714.4844 
T0 = 294.375305 
c = 344.108887 
g = 0
gamma = 1.40011787  
cv = 20.7801247 
M = 28.965 # masse molaire de l'air https://fr.wikipedia.org/wiki/Air

Cv = cv/M

param_toset_onmesh = [rho0, v0, p0, T0, c, g, gamma, Cv]


#%% CONSIDERATION ON MESH

to_check_CFL = False
if to_check_CFL : 
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
index_receivers = [150+nb_index_neg] #np.arange(z0+d_receivers, zmax, d_receivers)

# define Time
t = 0
dt = 0.005/10

#%% DEFINE A PRIORI MODEL ON MESH

for i,param in enumerate(param_toset_onmesh):
    param_toset_onmesh[i] = param * np.ones(len(z))

[rho0, v0, p0, T0, c, g, gamma, Cv] = param_toset_onmesh

## DEFINE THE "REAL" MODEL 

[GT_rho0, GT_v0, GT_p0, GT_T0,GT_c, GT_g, GT_gamma, GT_Cv] = deepcopy(param_toset_onmesh)

factor = 2
obs_start = 90 + nb_index_neg
obs_end = 120 + nb_index_neg

if True : # False : 
    for param in GT_v0, GT_rho0:#, GT_p0:
        param[obs_start:obs_end] = factor * param[obs_start:obs_end]


#%% DEFINE SOURCE ON MESH

# set the source
iz_source =  190#30
t_ax = np.arange(T_init,Tmax+dt,dt/2)
source = get_source(t_ax, f0) #scipy.stats.norm.pdf(t_ax,10,1.5)

dist_z = abs(z - z[iz_source])
dist_z = np.roll(dist_z,-170)

index_source = 50+nb_index_neg
factor_source = (z == z[index_source]) #np.exp(-2*dist_z/1000)
 

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


plt.figure()
plt.plot(t_ax, source[0], label="source")
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
U0.plot(z,"The  a priori background")
# Plot the model
GT_U0.plot(z,"The ground truth background")


#%% RESOLUTION

is_reverse = False
# Resolution in a 1d case
if not use_backups : 
    start_time = time.time()
    t_end, U_end, history_obs = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", GT_U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, is_reverse)
    print(" Total time : ", time.time() - start_time)

    np.save("../BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), np.array(history_obs))

else : 
    history_obs = np.load("../BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]+".npy"))
    U_end = history_obs[-1]

#%% RESULTS

U_end.plot(z,"Perturbation after Tmax = "+ str(Tmax))
plt.close()

#%% BACKWARD RESOLUTION
    
# Resolution in a 1d case
 
 # first do the forward pb with the a priori model
is_reverse = False
# Resolution in a 1d case
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
t_end, U_end, history_forwardapriori = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, is_reverse)


# make the checkpointing op
U_tmax = deepcopy(U_end)
is_reverse = True
start_time = time.time() 
# bakcward : 
#t_start, U_start, history_reverse = time_scheme(get_minus_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, reversed_source, is_reverse, history_obs)
# checkpointing : 
#t_start, U_start, history_reverse = time_scheme(get_RHS, U_end, T_init, Tmax, z,"checkpointing", U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, is_reverse, history_forwardapriori)
history_reverse = [history_forwardapriori[t] for t in range(len(history_obs)-1,-1,-1)]
U_start = history_forwardapriori[-1]
print(" Total time : ", time.time() - start_time)


#%% RESULTS BACKWARD

U_start.plot(z,"After Tmax = 0")

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
    ax[0,1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax[1,1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
    ax[2,1].axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
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
    
 
#%%

# Check the perturbation/source of the adjoint equations
recepteur = np.zeros((len(history_obs), len(index_receivers)))

max_obs_rho = 1
max_obs_v = [ max([history_obs[t].v[k] for t in range(len(history_obs))]) for k in index_receivers]
max_obs_p = [ max([history_obs[t].p[k] for t in range(len(history_obs))]) for k in index_receivers]
max_obs = [max_obs_rho, max_obs_v, max_obs_p]


plt.figure()
for t in range(len(history_reverse)):
    recepteur[t,:] = dchi(history_reverse, history_obs, t*dt, dt, index_receivers, max_obs)[2,index_receivers]
plt.plot(recepteur[:,0])
plt.grid()
plt.show()


#%% RESOLUTION ADJOINT EQUATIONS
                           
U_tmax = deepcopy(U_end)
t_start, U_start, history_adjoint = time_scheme(get_adjoint_RHS, U_end*0, T_init, Tmax, z, "adjoint", U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, history_reverse, history_obs, index_receivers, max_obs)

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


#%% KERNELS

K = get_kernels(history_adjoint, history_reverse, U0, T0, 0, gamma, Cv, dt, h) 

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

#%% Attendu 

def chi(U0, history_obs, index_receivers, max_obs, dt, name_file):
    # init
    CHI = np.zeros(3)
    # get Ucalc
    U1 = deepcopy(U0) * 0 
    t_end, U_end, history_calc = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, is_reverse)

    # Ucalc - Uobs
    for t in range(len(history_obs)):
        #CHI[0] += sum((history_calc[t].rho[index_receivers] - history_obs[t].rho[index_receivers])**2)
        CHI[2] += sum((history_calc[t].p[index_receivers] - history_obs[t].p[index_receivers])**2)
        CHI[1] += sum((history_calc[t].v[index_receivers] - history_obs[t].v[index_receivers])**2)
    CHI *= dt

    # normalisation
    #CHI[0] = CHI[0] / max_obs[0]
    CHI[2] = CHI[2] /  max_obs[2] 
    CHI[1] = CHI[1] / max_obs[1]

    np.save(name_file, CHI[1:])

    return CHI


start_k = 20+nb_index_neg
end_k = 180+nb_index_neg
drho = 0.5 * U0.rho[0]
chi_rho = chi(U0, history_obs, index_receivers, max_obs, dt, './Save_Kernel_info/model_m0')
K_i = np.zeros((3,len(np.arange(start_k,end_k+1,5))))
for it,i in enumerate(range(start_k, end_k+1,5)):
    # add a delta on one component
    U0_rho_di = deepcopy(U0)
    U0_rho_di.rho[i:i+5] += drho
    # compute chi for rho_i
    chi_rho_i = chi(U0_rho_di, history_obs, index_receivers, max_obs, dt, './Save_Kernel_info/model_m_i_'+str(i))
    # compute dX
    dX = chi_rho_i - chi_rho
    # compute dX/drho_i
    K_i[:,it] = dX / drho



# %%
plt.plot(np.arange(start_k,end_k+1,5),K_i[1,:])
#plt.plot(np.arange(start_k,end_k+1,5),K_i[2,:])
# %%
