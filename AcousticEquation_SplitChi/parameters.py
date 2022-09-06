#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
sys.path.insert(1, '../Lib/')

from discretisation import *
from linearised_navier_stokes_acoustic_adjoint_splitchi import *

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

# DEFINITION

####################
# SOURCE PARAMETERS
####################
f0 = 0.1


####################
# DOMAIN PARAMETERS
####################
z0 = -25e3
nb_index_neg = int(- z0 / 100)
zmax = 45e3


################## 
# TIME PARAMETERS
##################
T_init = 0
Tmax = 50


##################
# USER PARAMETERS
##################
display_anim = False
time_scheme = EE
local_path = '/home/deos/s.gerier/PROJECTS/SIMULATIONS/DF_1D/Inversion_LNS/AcousticEquation_SplitChi/'
which_chi = "pressure"

####################################
#  MODEL PARAMETERS (A PRIORI MODEL)
####################################
rho0  = 1.04898036      # density
v0    = 0               # wind
p0    = 88714.4844      # pressure
T0    = 294.375305      # temperature
c     = 344.108887      # sound speed
g     = 0               # gravity
gamma = 1.40011787      # gamma 
cv    = 20.7801247      # calorific capacity
M     = 28.965          # masse molaire de l'air https://fr.wikipedia.org/wiki/Air
Cv    = cv/M
param_toset_onmesh = [rho0, v0, p0, T0, c, g, gamma, Cv]

########################
# CONSIDERATION ON MESH
########################
to_check_CFL = False
if to_check_CFL : 
    min_v = abs(c) #, abs(v0))
    max_v = max(abs(c), abs(v0))
    max_dx = min_v / (10 * 2.5 * f0) 

    max_dt = max_dx / max_v

    print("With this source, you must choose a dx < ", str(max_dx) )
    print("If taking the dx max, finaly, you will have to choose a dt < ", str(max_dt))


##################
# MESH PARAMETERS
#################
# define the space
h = 100      # size of the mesh
z = np.arange(z0,zmax+h,h)
 
# define Time
t = 0
dt = 0.005/10


#######################
# RECEIVERS PARAMETERS
#######################
index_receivers = [150+nb_index_neg]


################################
# DEFINE A PRIORI MODEL ON MESH
################################
for i,param in enumerate(param_toset_onmesh):
    param_toset_onmesh[i] = param * np.ones(len(z))

[rho0, v0, p0, T0, c, g, gamma, Cv] = param_toset_onmesh


###########################
## DEFINE THE "REAL" MODEL 
###########################
[GT_rho0, GT_v0, GT_p0, GT_T0, GT_c, GT_g, GT_gamma, GT_Cv] = deepcopy(param_toset_onmesh)

factor = 1.2
obs_start = 90 + nb_index_neg
obs_end = 120 + nb_index_neg

for param in [GT_c]: #, GT_rho0:
    param[obs_start:obs_end] = factor * param[obs_start:obs_end]


########################
# DEFINE SOURCE ON MESH
########################
# set the source
index_source = 50+nb_index_neg
t_ax = np.arange(T_init,Tmax+dt,dt/2)
source = get_source(t_ax, f0) 
factor_source = (z == z[index_source]) #np.exp(-2*dist_z/1000)
 
source = [source, factor_source]


##############
# PLOT SOURCE
##############
if False:
    # plot the source in time and space
    total_source = np.zeros((len(z), len(t_ax)))
    for t in t_ax : 
        for i in range(len(z)):
            total_source[i,int(round(2*t/dt))] = source[0][int(round(2*t/dt))] * source[1][i]

    fig, ax = plt.subplots()
    plt.imshow(total_source, aspect='auto')
    plt.colorbar(label="Intensity")
    plt.title("Form of the source applied on density on point z="+ str(z[index__source]/1000)+"km")
    plt.xticks(range(0,len(t_ax),1000), t_ax[0:len(t_ax):1000])
    plt.yticks(range(0,len(z),20), z[0:len(z):20]/1000)
    plt.xlabel("Time")
    plt.ylabel("Distance (km)")


plt.figure()
plt.plot(t_ax, source[0], label="source")
plt.title("Form of the source applied on density on point z = "+str(z[index_source]/1000)+"km")
plt.xlabel('Time')
plt.ylabel("Intensity")
plt.show()
plt.close()


##########################
# INTIALISATION OF MODELS
##########################

# define vectors of the system
c_demi = (c[1:] + c[:-1])/2
GT_c_demi = (GT_c[1:] + GT_c[:-1])/2 
U0 = LNS_Model(rho0, c_demi) 
GT_U0 = LNS_Model(GT_rho0, GT_c_demi) 

# Plot the model
U0.plot(z,"The  a priori background")
# Plot the model
GT_U0.plot(z,"The ground truth background")


U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
