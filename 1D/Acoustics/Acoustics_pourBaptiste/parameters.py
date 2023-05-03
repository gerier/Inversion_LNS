#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from linearised_navier_stokes import *
from discretisation import *
import sys

from get_model_msis import get_model
sys.path.insert(1, '../Lib/')

plt.rcParams.update({'font.size': 13})


def get_true_model(choice, GT_model, arg):
    GT_model_out = deepcopy(GT_model)
    if choice == 1:  # un creneau pour la densité
        factor = arg[0]
        obs_start = arg[1]
        obs_end = arg[2]
        GT_model_out[obs_start:obs_end] = factor * GT_model[obs_start:obs_end]
    elif choice == 2:  # un créneau pour la vitesse de propagation des ondes
        pass
    elif choice == 3:  # un créneau pour la densité et la propagation des ondes
        pass
    elif choice == 4:  # plusieurs créneaux pour la densité et la vitesse de propagation des ondes
        pass
    elif choice == 5:  # plusieurs créneaux pour la densité et la vitesse de propagation des ondes
        pass
    elif choice == 6:  # model given by MSIS file
        [path_file, z, no_gravity, no_wind] = arg
        GT_model_out = get_model(path_file, z, no_gravity, no_wind)
    return GT_model_out

# DEFINITION


####################
# SOURCE PARAMETERS
####################
f0 = 0.1


####################
# DOMAIN PARAMETERS
####################
z0 = -20.e3  # -25e3
nb_index_neg = int(- z0 / 100)
zmax = 20.e3  # 45e3


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
local_path = "/home/deos/s.gerier/PROJECTS/SIMULATIONS/DF_1D/Inversion_LNS/AcousticEquation_Cleaned/"


####################################
#  MODEL PARAMETERS (A PRIORI MODEL)
####################################


########################
# CONSIDERATION ON MESH
########################
to_check_CFL = False
if to_check_CFL:
    min_v = abs(c)  # , abs(v0))
    max_v = max(abs(c), abs(v0))
    max_dx = min_v / (10 * 2.5 * f0)

    max_dt = max_dx / max_v

    print("With this source, you must choose a dx < ", str(max_dx))
    print("If taking the dx max, finaly, you will have to choose a dt < ", str(max_dt))


##################
# MESH PARAMETERS
#################
# define the space
h = 100      # size of the mesh
z = np.arange(z0, zmax+h, h)

# define Time
t = 0
dt = 0.05


#######################
# RECEIVERS PARAMETERS
#######################
# (np.linspace(200,500,6) + nb_index_neg).astype('int')# [400+nb_index_neg] #(np.linspace(0,180,181) + nb_index_neg).astype('int') #[100+nb_index_neg]
index_receivers = [int(15*1000/h)+nb_index_neg]


model_MSIS = False
if model_MSIS:
    ###########################
    # DEFINE THE "REAL" MODEL
    ###########################
    GT_model = get_true_model(6, [], ["", z, 0, 0])
    [GT_rho0, GT_T0, GT_c, GT_p0, GT_g, GT_kappa, GT_mu, GT_eta,
        GT_v0, GT_cv, GT_gamma, GT_l, GT_Cv] = deepcopy(GT_model)

    ###########################
    # DEFINE THE "APRIORI" MODEL
    ###########################
    AP_model = deepcopy(GT_model)
    for i in range(len(AP_model)):
        AP_model[i] = np.mean(AP_model[i]) * np.ones(len(AP_model[i]))
        #AP_model[i] = AP_model[i] *1.1
    [rho0, T0, c, p0, g, kappa, mu, eta, v0, cv, gamma, l, Cv] = AP_model


else:
    aux_model = get_true_model(6, [], ["", z, 0, 0])
    [aux_rho0, aux_T0, aux_c, aux_p0, aux_g, aux_kappa, aux_mu, aux_eta,
        aux_v0, aux_cv, aux_gamma, aux_l, aux_Cv] = deepcopy(aux_model)

    factor = 1.#1
    factor_c = 1.#2
    obs_start = int(16*1000/h + nb_index_neg)
    obs_end = int(17/h*1000 + nb_index_neg)

    # take ana altitude of the model and use this altitude to define an homogeneous model
    AP_model = deepcopy(aux_model)
    z_model = 100.
    index_z_model = np.where(z == z_model)[0]
    for i in range(len(AP_model)):
        AP_model[i][:] = aux_model[i][index_z_model[0]]
    [rho0, T0, c, p0, g, kappa, mu, eta, v0, cv, gamma, l, Cv] = AP_model
    c2 = (c)**2

    # add variation/discontinuities in the model
    GT_rho0 = get_true_model(1, rho0, [factor, obs_start, obs_end])
    GT_c2 = get_true_model(1, c2, [factor_c, obs_start, obs_end])

    GT_c2 += (50000*np.exp(-((z-5000)/5000)**2) / (np.sqrt(2 * np.pi * 5000)))**2

    # cancel the wind
    v0 = 0*v0
    GT_v0 = 0*v0
    # Modify the value of gamma
    gamma = gamma*0+1.4
    GT_gamma = gamma*0+1.4
    # compute the pressure since the model has been modified
    p0 = c2 * rho0 / gamma
    GT_p0 = GT_c2 * GT_rho0 / GT_gamma


########################
# DEFINE SOURCE ON MESH
########################
# set the source fr plots
index_source = int(0*1000/h+nb_index_neg)
t_ax = np.arange(T_init, Tmax+dt, dt/2)

source = [z, f0, index_source]


##########################
# INTIALISATION OF MODELS
##########################

# define vectors of the system
M = LNS_Model(rho0, v0[:-1], p0, gamma, g)
GT_M = LNS_Model(GT_rho0, GT_v0[:-1], GT_p0, GT_gamma, g)

# Plot the model
M.plot(z, "The  a priori background")
# Plot the model
GT_M.plot(z, "The ground truth background")

# plot rho, p and c
fig, ax = plt.subplots(2, 1, figsize=(10, 10))
ax[0].plot(z/1000, GT_M.p)
ax[0].set_ylabel("Pressure")
ax[1].plot(z/1000, np.sqrt(GT_c2))
ax[1].set_ylabel("Celerity")
for i in range(2):
    ax[i].grid()
    #ax[i].set_xlim(5, 35)
