import sys
sys.path.insert(1, './Lib/')
sys.path.insert(1, '../Lib/')

from discretisation import *
from linearised_navier_stokes_acoustic_adjoint import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
import time


# LOAD OBSERVATIONS AND SOLUTION
history_obs = np.load("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)
history_reverse = np.load("./BackUps/reverse_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)

history_reverse = history_reverse[::-1]

# DEFINE MODEL A PRIORI
model_aprior = deepcopy(U0)

# DEFINE DIFFERENT CASES TO STUDY 
cases = [ [-10,-5], [38,43], [55,60], [100,105], [147,152], [180,185]]

# DECIDE HICH PARAMETER YOU WANT TO OBSERVE
parameter = "pressure"

# PROCEEED PLOTS TO ANALYSE

for it_cas,cas in enumerate(cases):
    obs_start = cas[0] + nb_index_neg
    obs_end = cas[1] + nb_index_neg

    U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 

    perturb_model = deepcopy(model_aprior)
    perturb_model.rho[obs_start:obs_end] = 1.5 * perturb_model.rho[obs_start:obs_end]

    _,_,history_perturb_m = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", perturb_model, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, False)
    
    t_ax = np.arange(T_init,Tmax+dt,dt/2)
    t_ax = t_ax[::2]
    # PLOT THE SOLUTION, AT DIFFERENT TIME, ON ALL SPACE
    for t in np.array([10,20,30,40,50])/dt:
        if parameter == "velocity" : 
            plt.plot((z[1:]+z[:-1])/2,history_obs[int(t)].v, label="observation")
            plt.plot((z[1:]+z[:-1])/2,history_reverse[int(t)].v, label="solution m")
            plt.plot((z[1:]+z[:-1])/2,history_perturb_m[int(t)].v, label="solution m_i")
        else :
            plt.plot(z,history_obs[int(t)].p, label="observation")
            plt.plot(z,history_reverse[int(t)].p, label="solution m")
            plt.plot(z,history_perturb_m[int(t)].p, label="solution m_i")
        plt.axvspan(z[obs_start], z[obs_end], alpha=0.1, color='grey')
        plt.plot(z[index_source], 0, "xg")
        plt.plot(z[index_receivers], 0, "xr")
        plt.grid()
        plt.legend()
        plt.title("Solution in function of space, at time "+str(t*dt))
        plt.savefig("./Analysis_plots_"+parameter+"/Case"+str(it_cas)+"/solution_function_space_at_time"+str(t*dt)+".png")
        plt.show()

    # PLOT THE SOLUTION AT RECEIVER, ON ALL TIME
    if parameter == "velocity" : 
        plt.plot(t_ax,[history_obs[it].v[index_receivers[0]] for it in range(len(history_obs))], label="observation")
        plt.plot(t_ax,[history_reverse[it].v[index_receivers[0]] for it in range(len(history_obs))], label="solution m")
        plt.plot(t_ax,[history_perturb_m[it].v[index_receivers[0]] for it in range(len(history_obs))], label="solution m_i")
    else : 
        plt.plot(t_ax,[history_obs[it].p[index_receivers[0]] for it in range(len(history_obs))], label="observation")
        plt.plot(t_ax,[history_reverse[it].p[index_receivers[0]] for it in range(len(history_obs))], label="solution m")
        plt.plot(t_ax,[history_perturb_m[it].p[index_receivers[0]] for it in range(len(history_obs))], label="solution m_i")
    plt.grid()
    plt.legend()
    plt.title("Solution in function of time")
    plt.savefig("./Analysis_plots_"+parameter+"/Case"+str(it_cas)+"/solution_function_time.png")
    plt.show()    

    # PLOT DIFFERENCE BETWEEN SOLUTION AND OBSERVATIONS
    if parameter == "velocity" : 
        plt.plot(t_ax,[(history_reverse[it].v[index_receivers[0]] - history_obs[it].v[index_receivers[0]])**2 for it in range(len(history_obs))], label="solution m - obs")
        plt.plot(t_ax,[(history_perturb_m[it].v[index_receivers[0]] - history_obs[it].v[index_receivers[0]])**2 for it in range(len(history_obs))], label="solution m_i - obs")
    else : 
        plt.plot(t_ax,[(history_reverse[it].p[index_receivers[0]] - history_obs[it].p[index_receivers[0]])**2 for it in range(len(history_obs))], label="solution m - obs")
        plt.plot(t_ax,[(history_perturb_m[it].p[index_receivers[0]] - history_obs[it].p[index_receivers[0]])**2 for it in range(len(history_obs))], label="solution m_i - obs")
    plt.title("Difference with the observation")
    plt.legend()
    plt.grid()
    plt.savefig("./Analysis_plots_"+parameter+"/Case"+str(it_cas)+"/diff_function_time.png")
    plt.show()     

    # PLOT DIFFEFERENCE BETWEEN THE BOTH NORMS
    if parameter == "velocity" : 
        plt.plot(t_ax,[(history_reverse[it].v[index_receivers[0]] - history_obs[it].v[index_receivers[0]])**2 - (history_perturb_m[it].v[index_receivers[0]] - history_obs[it].v[index_receivers[0]])**2 for it in range(len(history_obs))])
    else:
        plt.plot(t_ax,[(history_reverse[it].p[index_receivers[0]] - history_obs[it].p[index_receivers[0]])**2 - (history_perturb_m[it].p[index_receivers[0]] - history_obs[it].p[index_receivers[0]])**2 for it in range(len(history_obs))])
    plt.title("Diffference of the norms")
    plt.grid()
    plt.savefig("./Analysis_plots_"+parameter+"/Case"+str(it_cas)+"/diff_norm_function_time.png")
    plt.show()  