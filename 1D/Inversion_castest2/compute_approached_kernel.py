#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys


from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import *
from parameters import *
from parametrisation import * 
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from get_norm import *
from misfit_function import * 

import time
from ApproximationDF.compute_approached_kernel_param import start_k, end_k, chi, saving_dir

which_grad = "rho"


#%% Attendu 


def get_approx_kernel(m0,argf, derivative="rho"): 

    d_pixel = 1

    obs, norm_info, dt, h, size, gamma, sponge_layer, source, index_receivers, Sc,parametrisation = argf

    #for derivative in ["rho", "wind", "pressure", "gamma"]: # 
    if True : 
        print("Treating : [", derivative,"]")

        M = model2paramLNSeq(m0, parametrisation, Sc, gamma, g, size)


        if derivative == "rho":
            dparam = 0.0001 * M.rho[0]
        elif derivative == "wind":
            if M.v[0] == 0:
                dparam = 0.0001
            else :  
                dparam = 0.0001 * M.v[0]
        elif derivative == "pressure": 
            dparam = 0.0001 * M.p[0]
        elif derivative == "gamma": 
            dparam = 0.0001 * M.gamma[0]

        
        chi_apriori = chi(M,obs,norm_info,dt,h,size,gamma,sponge_layer,source,index_receivers, Sc,parametrisation)

        K_i = np.zeros(len(np.arange(start_k,end_k+d_pixel,d_pixel)))

        for it,i in enumerate(range(start_k, end_k+d_pixel,d_pixel)):
            # add a delta on one component
            M_di = deepcopy(M)
            if derivative == "rho":
                M_di.rho[i:i+d_pixel] += dparam
            elif derivative == "wind":
                M_di.v[i:i+d_pixel] += dparam
            elif derivative == "pressure":
                M_di.p[i:i+d_pixel] += dparam
            elif derivative == "gamma":
                M_di.gamma[i:i+d_pixel] += dparam
            # compute chi for pa-rameter i
            chi_i = chi(M_di,obs,norm_info,dt,h,size,gamma,sponge_layer,source,index_receivers, Sc,parametrisation)
            # compute dX
            dX = chi_i - chi_apriori
            # compute dX/drho_i
            K_i[it] = h * dX / dparam


            # make temporary savings
            if (it % 20) <= 1e-4 :
                np.save("./"+saving_dir+"/kernel_"+derivative+"_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(start_k)+"_"+str(end_k), K_i)




        plt.figure()
        plt.plot(z[np.arange(start_k,end_k+1,d_pixel)], K_i)
        plt.plot(z[index_source],0)
        plt.plot(z[index_receivers],np.zeros(len(index_receivers)),'xg')
        plt.title("Approached kernel")
        plt.grid()
        plt.show()

        if True :
            np.save("./"+saving_dir+"/kernel_"+derivative+"_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(start_k)+"_"+str(end_k), K_i)

        return K_i,z[np.arange(start_k,end_k+1,d_pixel)]

#%%

type_parametrisation = 2
type_gradient = 1
type_regul = 0
 
# select parametrisation
if type_parametrisation == 0 :
    parametrisation = ["density", "wind", "velocity"]
    Sc = [1] * len(z) + [1e3] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(M.rho, c0), M.v)
    m = m / Sc
elif type_parametrisation == 1 :
    parametrisation = ["density", "wind", "pressure"]
    Sc = [1] * len(z) + [1e5] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(M.rho, M.p), M.v)
    m = m / Sc
elif type_parametrisation == 2 :
    parametrisation = ["log_density", "wind", "log_pressure"]
    Sc = [1] * len(z) + [1] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(np.log(M.rho), np.log(M.p)), M.v)
    m = m / Sc
elif type_parametrisation == 3 :
    parametrisation = ["log_density", "wind", "log_velocity"]
    Sc = [1] * len(z) + [1] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(np.log(M.rho), np.log(c0)), M.v)
    m = m / Sc


# for the figure
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})


# Define the observation that we want to reproduce
it_end = int(Tmax / dt)
size = len(z)

source = source_set1
index_source = index_source1

start = time.time()
if zero_as_initialcondition:
    U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size))
else:
    U1 = np.load("./initialcondition.npy", allow_pickle=True)[0]
    source[-1] = []

# init sponge layer
sponge_layer_info = get_sponge(size)


history_obs, Uobs = time_scheme(
    get_RHS, U1, 0, it_end, dt, index_receivers, "forward", sponge_layer_info, GT_M, h, source)
print(time.time() - start)

save_for_initial_condition = False
if save_for_initial_condition:
    Uobs.p[:500] = 0
    Uobs.p = np.roll(Uobs.p, -800)
    Uobs.rho[:500] = 0
    Uobs.rho = np.roll(Uobs.rho, -800)
    Uobs.v[:500] = 0
    Uobs.v = np.roll(Uobs.v, -800)(i)
    np.save("./initialcondition.npy", np.array([Uobs]))
    sys.exit()


# define a dictionary to save information on norms
# ord: order of the norm: 1,2,inf
# choice: on time, on receiver (TODO : develop explanation on the options)
# calc: boolean to know if synthetics are normalised by synthetics or by observations
norm_info = {'ord': 2, 'choice': 1, 'calc': False}
norm_obs = get_norm(history_obs, norm_info['ord'], norm_info['choice'])
norm_info["norm_obs"] = norm_obs

argf = [history_obs, norm_info, dt, h, size, gamma, sponge_layer_info, source, index_receivers, Sc,parametrisation]



initial_m = deepcopy(m)
Ki,zi = get_approx_kernel(initial_m,argf, which_grad)



#%%

M = model2paramLNSeq(initial_m, parametrisation, Sc, gamma, g, size)
grad, preconditioner =  compute_gradient(M, 0, it_end, dt, h, size, sponge_layer_info, source, index_receivers, history_obs, norm_info)

if which_grad == "wind":  
    index_grad = 1
elif which_grad == "pressure":
    index_grad = 2
elif which_grad == "rho":
    index_grad = 0   

plt.plot(z/1000,grad[index_grad,:])
plt.plot(zi/1000,Ki)
plt.show()
common_ind = np.zeros(len(z),dtype=bool)
for i in range(len(z)):
    common_ind[i] = z[i] in zi
plt.plot( z[common_ind]/1000, grad[index_grad,common_ind] - Ki)
plt.show()
