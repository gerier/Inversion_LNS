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
source_f0 = 0.05


####################
# DOMAIN PARAMETERS
####################
z0 = 0  # -25e3
nb_index_neg = int(- z0 / 100)
zmax = 60e3 #45e3  # 45e3


##################
# TIME PARAMETERS
##################
T_init = 0
Tmax = 160


##################
# USER PARAMETERS
##################
display_anim = False
time_scheme = EE



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
h =150      # size of the mesh
z = np.arange(z0, zmax+h, h)

# define Time
t = 0
dt = 0.1#25


########################
# DEFINE SOURCE ON MESH
########################
# set the source fr plots
index_source1 = int(12*1000/h+nb_index_neg)
source_set1 = [z, source_f0, Tmax, index_source1]

index_source2 = int(40*1000/h+nb_index_neg)
source_set2 = [z, source_f0, Tmax, index_source2]

#######################
# RECEIVERS PARAMETERS
#######################
index_receivers = (np.arange(8,46,8) *1000/h).astype('int') # [int(20*1000/h), int(25*1000/h), int(32*1000/h)]#
#index_receivers = [150-45,150+45]
#index_receivers = (np.arange(9,52,5) *1000/h).astype('int') 
index_receivers = (np.arange(9,52,5) *1000/h).astype('int') 


aux_model = get_true_model(6, [], ["./Flores_atmosphere_500km_86p2_4.dat", z, 0, 0])
[aux_rho0, aux_T0, aux_c, aux_p0, aux_g, aux_kappa, aux_mu, aux_eta,aux_v0, aux_cv, aux_gamma, aux_l, aux_Cv] = deepcopy(aux_model)


factor_rho = 0#0.1
factor_c = -0.015
obs_start = int(15*1000/h + nb_index_neg) #int(20*1000/h + nb_index_neg)
obs_end = int(23*1000/h + nb_index_neg) #int(30/h*1000 + nb_index_neg)

AP_model = deepcopy(aux_model)
z_model = 200.
index_z_model = np.where(z == z_model)[0]
for i in range(len(AP_model)):
    AP_model[i][:] = aux_model[i] #[index_z_model[0]]
[rho0, T0, c, p0, g, kappa, mu, eta, v0, cv, gamma, l, Cv] = AP_model

#rho0 = rho0[1]*np.ones(len(rho0))
#c = c[1]*np.ones(len(rho0))

GT_rho0 = deepcopy(rho0)
#GT_rho0[obs_start:obs_end] *= 1.+ factor_rho * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  200/(z[obs_end]-z[obs_start])   )

c2 = c**2
GT_c2 = deepcopy(c2)
#GT_c2[obs_start:obs_end] *= 1.+ factor_c * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  200/(z[obs_end]-z[obs_start])   )

#km_max = 20
#start = 5
#c2[int(start*1e3/h):int((start+km_max)*1e3/h)] *= 1- 0.1* (0.5-0.5*np.cos(2*np.pi*(z[int(start*1e3/h):int((km_max+start)*1e3/h)]-start*1e3)/((km_max)*1e3)))
#km_max = 20
#start = 25
#c2[int(start*1e3/h):int((start+km_max)*1e3/h)] *= 1+ 0.1* (0.5-0.5*np.cos(2*np.pi*(z[int(start*1e3/h):int((km_max+start)*1e3/h)]-start*1e3)/((km_max)*1e3)))
c2[:] = (np.sqrt(GT_c2[index_source1])*1.02 )**2 

factor_rho = 0.05
factor_c = 0.015
#obs_start = int(28*1000/h + nb_index_neg) #int(20*1000/h + nb_index_neg)
#obs_end = int(32*1000/h + nb_index_neg) #int(30/h*1000 + nb_index_neg)

#GT_rho0[obs_start:obs_end] *= 1 + factor_rho * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  50/(z[obs_end]-z[obs_start])   )
#GT_c2[obs_start:obs_end] *= 1 + factor_c * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  50/(z[obs_end]-z[obs_start])   )
#GT_c2 *= 1.05


v0 = 0*v0#+30
GT_v0 = deepcopy(v0) #+100

#GT_v0[obs_start:obs_end] = 5 * np.exp(  - (z[obs_start:obs_end] - ( z[obs_end]+z[obs_start])/2)**2 /  400/(z[obs_end]-z[obs_start])   )


#gamma = 1.4*np.ones(len(rho0))
GT_gamma = gamma #*np.ones(len(rho0))

g = 9.81*np.ones(len(rho0))
GT_g = g#*np.ones(len(rho0))

p0 = c2 * rho0 / gamma
GT_p0 = GT_c2 * GT_rho0 / GT_gamma




##########################
# INTIALISATION OF MODELS
##########################

# define vectors of the system
v0 = v0[:-1]
GT_v0 = GT_v0[:-1]

M = LNS_Model(rho0, v0, p0, gamma, g)
GT_M = LNS_Model(GT_rho0, GT_v0, GT_p0, GT_gamma,GT_g)
# Plot the model
M.plot(z, "The  a priori background")
# Plot the model
GT_M.plot(z, "The ground truth background")
plt.show()

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})
interp_z =( z[1:] + z[:-1] )/2
fig, ax = plt.subplots(4,1, figsize=(10,8))
ax[0].semilogy(z/1000, GT_M.rho, 'red')
ax[0].semilogy(z/1000, M.rho, 'blue')
ax[0].set_ylabel("Density (kg/m3)")
ax[1].semilogy(z/1000, GT_M.p, 'red')
ax[1].semilogy(z/1000, M.p, 'blue')
ax[1].set_ylabel("Pressure (Pa)")
ax[2].plot(z/1000, np.sqrt(GT_c2), 'red')
ax[2].plot(z/1000, np.sqrt(c2), 'blue')
ax[2].set_ylabel("Celerity (m/s)")
ax[3].plot(interp_z/1000, GT_v0, 'red')
ax[3].plot(interp_z/1000, v0, 'blue')
ax[3].set_ylabel("Wind (m/s)")
ax[3].set_xlabel('Altitude (km)')
for i in range(4):
    ax[i].grid()
    #ax[i].set_xlim(0,35)
    yl = ax[i].get_ylim()
    ax[i].vlines(z[index_source1]/1000,*yl,colors='k',linestyles='--')
    for ir in index_receivers : 
        ax[i].vlines(z[ir]/1000,*yl,colors='slategrey',linestyles='dotted')
    ax[i].set_ylim(yl)


##########################
# INTIAL CONDITION
##########################
zero_as_initialcondition = True


##########################
# PARAMETRISATION
##########################
# cf. main 