#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

from discretisation import *
from linearised_navier_stokes import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from compute_approached_kernel_param import start_k, end_k, chi, saving_dir

#%% Attendu 


history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+".npy", allow_pickle=True)

max_obs_rho = [ max([abs(history_obs[t].rho[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs_v = [ max([abs(history_obs[t].v[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs_p = [ max([abs(history_obs[t].p[k]) for t in range(len(history_obs))]) for k in index_receivers]
max_obs = [ np.max(max_obs_rho)**2, np.max(max_obs_v)**2, np.max(max_obs_p)**2]


derivative = "rho"
d_pixel = 1


for derivative in ["rho", "wind", "pressure", "gamma"]: # 
    print("Treating : [", derivative,"]")

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

    chi_apriori = chi(M, history_obs, index_receivers, max_obs, dt)

    K_i = np.zeros(len(np.arange(start_k,end_k+d_pixel,d_pixel)))
    K_approach = np.zeros(len(K_i))

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
        # compute chi for parameter i
        chi_i = chi(M_di, history_obs, index_receivers, max_obs, dt)
        # compute dX
        dX = chi_i - chi_apriori
        # compute dX/drho_i
        K_i[it] = dX / dparam


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


#%%

fichier = open("parameters.txt", "w")
fichier.write("# Domaine\n")
fichier.write(str(z0)+ ","+ str(zmax)+ ","+ str(h)+'\n')
fichier.write("# Source\n")
fichier.write(str(index_source)+ ','+str(z[index_source])+"\n")
fichier.write(str(f0)+"\n")
fichier.write("# Receivers\n")
fichier.write(str(index_receivers)+ ','+ str(z[index_receivers])+"\n")
fichier.write("# Domain approx\n")
fichier.write(str(start_k)+ ','+ str(end_k)+"\n")
fichier.write("# Time\n")
fichier.write(str(T_init)+ ','+ str(Tmax)+','+str(dt)+"\n")
fichier.close()

np.savetxt("model_density_apriori.csv", M.rho, delimiter=",")
np.savetxt("model_wind_apriori.csv", M.v, delimiter=",")
np.savetxt("model_pressure_apriori.csv", M.p, delimiter=",")
np.savetxt("model_gamma_apriori.csv", M.gamma, delimiter=",")
np.savetxt("model_density_real.csv", GT_M.rho, delimiter=",")
np.savetxt("model_wind_real.csv", GT_M.v, delimiter=",")
np.savetxt("model_pressure_real.csv", GT_M.p, delimiter=",")
np.savetxt("model_gamma_real.csv", GT_M.gamma, delimiter=",")
