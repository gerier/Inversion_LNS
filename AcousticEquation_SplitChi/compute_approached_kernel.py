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

def chi(U0, history_obs, index_receivers, max_obs, dt, which_chi):
    # init
    CHI = np.zeros(3)
    # get Ucalc
    U1 = deepcopy(U0) * 0 
    t_end, U_end, history_calc = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, False)

    # Ucalc - Uobs
    for t in range(len(history_obs)):
        if which_chi == "density" : 
            CHI[0] += sum((history_calc[t].rho[index_receivers] - history_obs[t].rho[index_receivers])**2)
        elif which_chi == "velocity" :
            CHI[2] += sum((history_calc[t].p[index_receivers] - history_obs[t].p[index_receivers])**2)
        elif which_chi == "pressure" :
            CHI[1] += sum((history_calc[t].v[index_receivers] - history_obs[t].v[index_receivers])**2)
    CHI *= dt

    # normalisation
    CHI[0] = CHI[0] / max_obs[0]
    CHI[2] = CHI[2] /  max_obs[2] 
    CHI[1] = CHI[1] / max_obs[1]

    return CHI


history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)

max_obs_rho = 1
max_obs_v = [ max([history_obs[t].v[k] for t in range(len(history_obs))]) for k in index_receivers]
max_obs_p = [ max([history_obs[t].p[k] for t in range(len(history_obs))]) for k in index_receivers]
max_obs = [max_obs_rho, max_obs_v, max_obs_p]

start_k = 0+nb_index_neg
end_k = 200+nb_index_neg



order_collaboration = [["velocity", "pressure"], ["pressure", "velocity"]]
collaboration = order_collaboration[0]

step = 0.01


# define the model 
m_aprior = deepcopy(U0)

for contrib in [0,1]:
    # get the perturbation
    drho = 0.5 * m_aprior.rho[0]
    
    # get the value of chi whith the a piori model
    chi_rho = chi(m_aprior, history_obs, index_receivers, max_obs, dt, which_chi)

    # init kernel
    K_i = np.zeros((3,len(np.arange(start_k,end_k+1,5))))
    
    # define the chi function to use / define the parameter whose sensitivity to the model is to be evaluated
    which_chi = collaboration[contrib]

    for it,i in enumerate(range(start_k, end_k+1,5)):
        # add a delta on one component
        m_aprior_rho_di = deepcopy(m_aprior)
        m_aprior_rho_di.rho[i:i+5] += drho
        # compute chi for rho_i
        chi_rho_i = chi(m_aprior_rho_di, history_obs, index_receivers, max_obs, dt, which_chi)
        # compute dX
        dX = chi_rho_i - chi_rho
        # compute dX/drho_i
        K_i[:,it] = dX / drho

    # get the approached kernel
    K_approach = np.zeros(len(K_i[0,:]))
    K_approach = K_i[1,:] + K_i[2,:]

    #update model
    m_aprior = m_aprior + step * K_approach**(-1) 

    # plot the resulted kernel
    plt.figure()
    plt.plot(z[np.arange(start_k,end_k+1,5)], K_approach)
    plt.plot(z[index_source],0)
    plt.plot(z[index_receivers],0)
    plt.title("Approached kernel")
    plt.grid()
    plt.show()

    # save the resulted kernel
    np.save("./BackUps/kernel_contrib"+which_chi+str(contrib)+"_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_approach)
    np.save("./BackUps/kernel_contrib"+which_chi+str(contrib)+"_v_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[1,:])
    np.save("./BackUps/kernel__contrib"+which_chi+str(contrib)+"_p_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), K_i[2,:])




