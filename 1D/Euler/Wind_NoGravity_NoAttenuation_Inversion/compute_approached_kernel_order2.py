#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
sys.path.insert(1, '../Lib/')

from discretisation import *
from linearised_navier_stokes_3eq import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


#%% Attendu 

def chi(U0, history_obs, index_receivers, max_obs, dt):
    # init
    CHI = np.zeros(3)
    # get Ucalc
    U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 

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

start_k = -70+nb_index_neg#-100+nb_index_neg #-100
end_k = 230+nb_index_neg #300+nb_index_neg #300


derivative = "rho"
d_pixel = 1
if derivative == "rho":
    dparam = 0.0001 * U0.rho[0]
elif derivative == "c": 
    dparam = 0.0001 * U0.c[0]

chi_apriori = chi(U0, history_obs, index_receivers, max_obs, dt)

K_i = np.zeros((3,len(np.arange(start_k,end_k+d_pixel,d_pixel))))
K_approach = np.zeros(len(K_i[0,:]))

coef_dparam = [-1, 1, 2]
coef_chi_dparam = [-1/3, 1, -1/6]

for it,i in enumerate(range(start_k, end_k+d_pixel,d_pixel)):
    # add a delta on one component
    dX = - chi_apriori / 2 
    for alpha in range(len(coef_dparam)):

        U0_di = deepcopy(U0)
        if derivative == "rho":
            U0_di.rho[i:i+d_pixel] += coef_dparam[alpha] * dparam
        elif derivative == "c":
            U0_di.c[i:i+d_pixel] += coef_dparam[alpha] * dparam
        # compute chi for parameter i
        chi_i = chi(U0_di, history_obs, index_receivers, max_obs, dt)
        # compute dX
        dX += coef_chi_dparam[alpha] * chi_i

    # compute dX/drho_i
    K_i[:,it] = dX / dparam

    # make temporary savings
    if (it % 20) <= 1e-4 :
        K_approach = K_i[1,:] + K_i[2,:]
        np.save("./BackUps/kernel_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_approach)
        np.save("./BackUps/kernel_v_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[1,:])
        np.save("./BackUps/kernel_p_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[2,:])



K_approach = K_i[1,:] + K_i[2,:]


plt.figure()
plt.plot(z[np.arange(start_k,end_k+1,d_pixel)], K_approach)
plt.plot(z[index_source],0)
plt.plot(z[index_receivers],0)
plt.title("Approached kernel")
plt.grid()
plt.show()

if True :
    np.save("./BackUps/kernel_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_approach)
    np.save("./BackUps/kernel_v_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[1,:])
    np.save("./BackUps/kernel_p_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[2,:])

