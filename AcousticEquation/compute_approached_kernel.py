#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
sys.path.insert(1, '../Lib/')

from discretisation import *
from linearised_navier_stokes_acoustic_adjoint import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


#%% Attendu 

def chi(U0, history_obs, index_receivers, max_obs, dt):
    # init
    CHI = np.zeros(3)
    # get Ucalc
    U1 = deepcopy(U0) * 0 
    _, _, history_calc = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, False)

    # Ucalc - Uobs
    for t in range(len(history_obs)):
        #CHI[0] += sum((history_calc[t].rho[index_receivers] - history_obs[t].rho[index_receivers])**2)
        CHI[2] += sum((history_calc[t].p[index_receivers] - history_obs[t].p[index_receivers])**2)
        CHI[1] += sum((history_calc[t].v[index_receivers] - history_obs[t].v[index_receivers])**2)
    CHI *= dt

    if False :
        # plots to test
        plt.figure()
        plt.plot([ history_calc[t].p[index_receivers] for t in range(len(history_obs))], label="calc p")
        plt.plot([ history_obs[t].p[index_receivers] for t in range(len(history_obs))], label="obs p")
        plt.legend()
        plt.figure()
        plt.plot([history_calc[t].p[index_receivers] - history_obs[t].p[index_receivers] for t in range(len(history_obs))], label="diff")
        plt.figure()
        plt.plot([(history_calc[t].p[index_receivers] - history_obs[t].p[index_receivers])**2 for t in range(len(history_obs))], label="diff ** 2")
        plt.show()

        plt.figure()
        plt.plot([ history_calc[t].v[index_receivers] for t in range(len(history_obs))], label="calc v")
        plt.plot([ history_obs[t].v[index_receivers] for t in range(len(history_obs))], label="obs v")
        plt.legend()
        plt.figure()
        plt.plot([history_calc[t].v[index_receivers] - history_obs[t].v[index_receivers] for t in range(len(history_obs))], label="vdiff")
        plt.figure()
        plt.plot([(history_calc[t].v[index_receivers] - history_obs[t].v[index_receivers])**2 for t in range(len(history_obs))], label="vdiff **2")
        plt.show()

        print("For pressure : ", CHI[2], ", for velocity : ", CHI[1])
    
    # normalisation
    #CHI[0] = CHI[0] / max_obs[0]
    CHI[2] = CHI[2] / max_obs[2] 
    CHI[1] = CHI[1] / max_obs[1]

    return CHI/2


history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)

max_obs_rho = 1
max_obs_v = [ max([abs(history_obs[t].v[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs_p = [ max([abs(history_obs[t].p[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs = [max_obs_rho, max_obs_v, max_obs_p]

start_k = -120+nb_index_neg#-100+nb_index_neg #-100
end_k = 320+nb_index_neg #300+nb_index_neg #300

drho = 0.01 * U0.rho[0]
drho_pixel = 1

chi_rho = chi(U0, history_obs, index_receivers, max_obs, dt)

K_i = np.zeros((3,len(np.arange(start_k,end_k+drho_pixel,drho_pixel))))
K_approach = np.zeros(len(K_i[0,:]))

for it,i in enumerate(range(start_k, end_k+drho_pixel,drho_pixel)):
    #print("Didcontinuité à ", i-nb_index_neg)
    # add a delta on one component
    U0_rho_di = deepcopy(U0)
    U0_rho_di.rho[i:i+drho_pixel] += drho
    # compute chi for rho_i
    chi_rho_i = chi(U0_rho_di, history_obs, index_receivers, max_obs, dt)
    # compute dX
    dX = chi_rho_i - chi_rho
    # compute dX/drho_i
    K_i[:,it] = dX / drho

    # make temporary savings
    if (it % 20) <= 1e-4 :
        K_approach = K_i[1,:] + K_i[2,:]
        np.save("./BackUps/kernel_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_approach)
        np.save("./BackUps/kernel_v_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[1,:])
        np.save("./BackUps/kernel_p_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[2,:])



K_approach = K_i[1,:] + K_i[2,:]


plt.figure()
plt.plot(z[np.arange(start_k,end_k+1,drho_pixel)], K_approach)
plt.plot(z[index_source],0)
plt.plot(z[index_receivers],0)
plt.title("Approached kernel")
plt.grid()
plt.show()

if True :
    np.save("./BackUps/kernel_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_approach)
    np.save("./BackUps/kernel_v_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[1,:])
    np.save("./BackUps/kernel_p_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[2,:])

