#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from Lib_NavierStokesEq.linearised_navier_stokes import *
from Lib_NavierStokesEq.discretisation import *
import sys

from Lib_NavierStokesEq.get_model_msis import get_model
sys.path.insert(1, '../Lib/')

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

N_frames = 3200

def get_true_model(choice, GT_model, arg):
    GT_model_out = deepcopy(GT_model)
    if choice == 1:  # un creneau pour la densité
        factor = arg[0]
        obs_start = arg[1]
        obs_end = arg[2]
        GT_model_out[0][obs_start:obs_end] = factor * \
            GT_model[0][obs_start:obs_end]
    elif choice == 2:  # un créneau pour la vitesse de propagation des ondes
        factor = arg[0]
        obs_start = arg[1]
        obs_end = arg[2]
        GT_model_out[1][obs_start:obs_end] = factor * \
            GT_model[1][obs_start:obs_end]
    elif choice == 3:  # un créneau pour la densité et la propagation des ondes
        factor = arg[0]
        obs_start = arg[1]
        obs_end = arg[2]
        for i in range(len(GT_model_out)):
            GT_model_out[i][obs_start:obs_end] = factor[i] * \
                GT_model[i][obs_start:obs_end]
    elif choice == 4:  # plusieurs créneaux pour la densité et la vitesse de propagation des ondes
        for obs in arg:
            factor = obs[0]
            obs_start = obs[1]
            obs_end = obs[2]
            print(obs_start)
            GT_model_out[0][obs_start:obs_end] = factor * \
                GT_model[0][obs_start:obs_end]
    elif choice == 5:  # plusieurs créneaux pour la densité et la vitesse de propagation des ondes
        gamma, p = arg[0][:]
        for obs in arg[1]:
            factor = obs[0]
            obs_start = obs[1]
            obs_end = obs[2]
            print(obs_start)
            GT_model_out[0][obs_start:obs_end] = factor * \
                GT_model[0][obs_start:obs_end]
            GT_model_out[1][obs_start:obs_end] = np.sqrt(
                gamma * p / GT_model[0][obs_start:obs_end])
    elif choice == 6:  # model given by MSIS file
        [path_file, z, no_gravity, no_wind] = arg
        GT_model_out = get_model(path_file, z, no_gravity, no_wind)
    return GT_model_out

# DEFINITION


####################
# SOURCE PARAMETERS
####################
source_f0 = 0.1


####################
# DOMAIN PARAMETERS
####################
z0 = 0  # -25e3
nb_index_neg = int(- z0 / 100)
zmax = 25e3 #45e3  # 45e3


##################
# TIME PARAMETERS
##################
T_init = 0
Tmax = 110


##################
# USER PARAMETERS
##################
display_anim = False
time_scheme = EE
local_path = "/home/deos/s.gerier/PROJECTS/SIMULATIONS/DF_1D/Inversion_LNS/AcousticEquation_Cleaned/"


####################################
#  MODEL PARAMETERS (A PRIORI MODEL)
####################################
v0 = 0               # wind
c = 344.108887      # sound speed
g = 9.81

########################
# CONSIDERATION ON MESH
########################
to_check_CFL = False
if to_check_CFL:
    min_v = abs(c)  # , abs(v0))
    max_v = max(abs(c), abs(v0))
    max_dx = min_v / (10 * 2.5 * source_f0)

    max_dt = max_dx / max_v

    print("With this source, you must choose a dx < ", str(max_dx))
    print("If taking the dx max, finaly, you will have to choose a dt < ", str(max_dt))


##################
# MESH PARAMETERS
#################
# define the space
h = 100      # size of the mesh
z_aux = np.arange(z0, zmax+h, h)
z = np.zeros(2*len(z_aux))
z[:len(z_aux)] = z_aux
z[len(z_aux):] = z_aux + z_aux[-1]

# define Time
t = 0
dt = 0.05#25


#######################
# RECEIVERS PARAMETERS
#######################
# (np.linspace(200,500,6) + nb_index_neg).astype('int')# [400+nb_index_neg] #(np.linspace(0,180,181) + nb_index_neg).astype('int') #[100+nb_index_neg]
index_receivers = [int(32*1000/h)] # (np.linspace(50,1600,1551) + nb_index_neg).astype('int') #[int(20*1000/h)+nb_index_neg]
index_receivers = (np.arange(6,23,1) *1000/h).astype('int') # [int(20*1000/h), int(25*1000/h), int(32*1000/h)]#


aux_model = get_true_model(6, [], ["./Flores_atmosphere_500km_86p2_4.dat", z, 0, 0])
[aux_rho0, aux_T0, aux_c, aux_p0, aux_g, aux_kappa, aux_mu, aux_eta,aux_v0, aux_cv, aux_gamma, aux_l, aux_Cv] = deepcopy(aux_model)


factor_rho = 0.1
factor_c = -0.015
obs_start = int(15*1000/h + nb_index_neg) #int(20*1000/h + nb_index_neg)
obs_end = int(23*1000/h + nb_index_neg) #int(30/h*1000 + nb_index_neg)

AP_model = deepcopy(aux_model)
z_model = 100.
index_z_model = np.where(z == z_model)[0]
for i in range(len(AP_model)):
    AP_model[i][:] = aux_model[i] #[index_z_model[0]]
[rho0, T0, c, p0, g, kappa, mu, eta, v0, cv, gamma, l, Cv] = AP_model

GT_rho0 = deepcopy(rho0)
GT_rho0[obs_start:obs_end] *= 1+ factor_rho * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  200/(z[obs_end]-z[obs_start])   )

c2 = c**2
GT_c2 = deepcopy(c2)
GT_c2[obs_start:obs_end] *= 1+ factor_c * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  200/(z[obs_end]-z[obs_start])   )

factor_rho = 0.2
factor_c = -0.015
obs_start = int(28*1000/h + nb_index_neg) #int(20*1000/h + nb_index_neg)
obs_end = int(30*1000/h + nb_index_neg) #int(30/h*1000 + nb_index_neg)

GT_rho0[obs_start:obs_end] *= 1+ factor_rho * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  50/(z[obs_end]-z[obs_start])   )
GT_c2[obs_start:obs_end] *= 1+ factor_c * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  50/(z[obs_end]-z[obs_start])   )



v0 =  0*v0#+30
GT_v0 = 0*deepcopy(v0) #+100

gamma = 1.4
GT_gamma = gamma

p0 = c2 * rho0 / gamma
GT_p0 = GT_c2 * GT_rho0 / GT_gamma



########################
# DEFINE SOURCE ON MESH
########################
# set the source fr plots
index_source = int(10*1000/h+nb_index_neg)
t_ax = np.arange(T_init, Tmax+dt, dt/2)

source = [z, source_f0, index_source]


##########################
# INTIALISATION OF MODELS
##########################

# define vectors of the system
v0 = v0[:-1]
GT_v0 = GT_v0[:-1]

M = LNS_Model(rho0, v0, p0, gamma)
GT_M = LNS_Model(GT_rho0, GT_v0, GT_p0, GT_gamma)
z == z[index_source]
# Plot the model
M.plot(z, "The  a priori background")
# Plot the model
GT_M.plot(z, "The ground truth background")
plt.show()

fig, ax = plt.subplots(3,1, figsize=(10,10))
ax[0].plot(z/1000, GT_M.rho)
ax[0].set_ylabel("Density")
ax[1].plot(z/1000, GT_M.p)
ax[1].set_ylabel("Pressure")
ax[2].plot(z/1000, np.sqrt(GT_c2))
ax[2].set_ylabel("Celerity")
for i in range(3):
    ax[i].grid()
    ax[i].set_xlim(5,35)



##########################
# INTIAL CONDITION
##########################
zero_as_initialcondition = True


##########################
# PARAMETRISATION
##########################
# cf. main 