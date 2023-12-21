import numpy as np
import matplotlib.pyplot as plt

from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import *
from parameters import *
from Lib_NavierStokesEq.compute_kernels import *
from Lib_NavierStokesEq.checkpointing import *
from get_norm import *
from optimisation_lib_reg import *
from misfit_function import *

import os
import glob
import time

import config
import shutil
import sys

from copy import deepcopy

from Lib_NavierStokesEq.filter import *


type_parametrisation = 2
type_gradient = 5
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


###########
def get_last_test(all_test):
    if all_test == []:
        last_test = 0
    else:
        last_test = int(all_test[0][:-1].split('_')[1])
        for test in all_test:
            current_test = int(test[:-1].split('_')[1])
            if current_test > last_test:
                last_test = current_test
    return last_test


###########

# do the inverse problem for several source frequency
source_freq_cut = 1 / 2**np.arange(2,-1,-1) * source_f0
total_time = np.arange(0,Tmax+dt,dt)
t0 = 1.2/source[1]
pi = 3.141592653589793238462643
signal_source =  - 8 * pi**2 * source[1]**2 * (total_time-t0) * np.exp(-4 * pi**2 * source[1]**2 * (total_time-t0)**2)
source_param = deepcopy(source)

initial_m = deepcopy(m)

for factor_reg in [0.005,0.1, 0.5, 1]:#0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 

    # generate a directory to save results
    if not(os.path.exists("./Resultats/")):
        os.mkdir("./Resultats/")
    all_test = glob.glob("./Resultats/Test_*/")
    last_test = get_last_test(all_test)
    path_test = "./Resultats/Test_"+str(last_test+1)+"_reg_"+str(factor_reg)
    os.mkdir(path_test)

    save_info = [path_test, plot_flatteniter, plot_gradient, plot_hessian, [z, parametrisation, Sc, gamma]]

    # init config parameter
    config.n_iter_alpha = 0
    config.n_fx = 0
    config.n_grad = 0
    config.n_h = 0

    reg = [factor_reg, h, size, deepcopy(initial_m)]


    new_m = deepcopy(initial_m)
    for source_freq in source_freq_cut:
        print("Test with source frequency cutoff:", source_freq_cut)

        # compute observation with this source
        history_obs_filtered = filter_observations(history_obs,source_freq,1/dt)

        # compute the norm of the observations (and add it to the dictionary)
        norm_obs = get_norm(history_obs_filtered, norm_info['ord'], norm_info['choice'])
        norm_info["norm_obs"] = norm_obs

        # filter the source for the synthetics 
        source = source_param + [dt,signal_source,source_freq]

        # define the arguments of the cost function and its derivative
        argf = [history_obs_filtered, norm_info, dt, h, size, gamma, sponge_layer_info, source, index_receivers, Sc,parametrisation]
        argg = [it_end, dt, h, size, sponge_layer_info, source, index_receivers,
                history_obs_filtered, norm_info, gamma, Sc,parametrisation]
        argh = [size]  # [len(history_obs), dt, index_receivers, h, gamma, source]

        # treat the path to save results
        path_test_fr = path_test + "/Step_F0_"+str(source_freq)
        os.mkdir(path_test_fr)
        save_info[0] = path_test_fr

        # compute the model by the inversion
        m = deepcopy(new_m)
        new_m, is_notfound, x_iter = optimisation(m, f0, df0, hessf0, 0.0001, 0.9, argf, argg, argh, alpha_max=1, maxiter=3,
		                                tol_x=1e-8, tol_df=1e-8, true_hessian=False, regularisation=reg, plot_iter=[plot_flatteniter, [z, parametrisation, Sc, gamma]], save_info=save_info, type_gradient=type_gradient,type_regul=type_regul)


        rho, p, c, v = model2allparam(new_m, parametrisation, Sc, gamma, size)
        fig, ax = plt.subplots(4, 1, figsize=(10, 7))
        ax[0].plot(z/1000, rho, label="Result")
        ax[1].plot(z/1000, p, label="Result")
        ax[2].plot(z[:-1]/1000, v, label="Result")
        ax[3].plot(z/1000, c, label="Result")
        ax[0].plot(z/1000, GT_M.rho, label="True")
        ax[1].plot(z/1000, GT_M.p, label="True")
        ax[2].plot(z[:-1]/1000, GT_M.v, label="True")
        ax[3].plot(z/1000, np.sqrt((GT_M.p * GT_M.gamma)/GT_M.rho), label="True")
        ax[0].set_ylabel("Density")
        ax[1].set_ylabel("Pressure")
        ax[2].set_ylabel("Velocity")
        ax[3].set_ylabel("Celerity")
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        ax[3].grid()
        ax[0].legend()
        ax[1].legend()
        ax[2].legend()
        ax[3].legend()
        ax[0].set_title("Comparaison true model and model obtained after inversion")
        ax[3].set_xlabel("Range(km)")
        ax[0].set_xlim(0, 40)
        ax[1].set_xlim(0, 40)
        ax[2].set_xlim(0, 40)
        ax[3].set_xlim(0, 40)
        plt.savefig(path_test_fr+"/solutions.png")
        plt.close()

        solution = np.loadtxt(path_test_fr+"/iterations_informations.txt")
        plt.figure()
        plt.plot(solution[:, 0], solution[:, 3])
        plt.title("Evolution of the cost function")
        plt.savefig(path_test_fr+"/cost_function.png")
        plt.close()

        plt.figure()
        plt.semilogy(solution[:, 0], solution[:, 3])
        plt.title("Evolution of the cost function (log scale)")
        plt.savefig(path_test_fr+"/cost_function_log.png")
        plt.close()

        np.save(path_test_fr+"/solution_density.npy", rho)
        np.save(path_test_fr+"/solution_pressure.npy", p)
        np.save(path_test_fr+"/solution_v.npy", v)

        np.save(path_test_fr+"/aimed_density.npy", GT_M.rho)
        np.save(path_test_fr+"/aimed_wind.npy", GT_M.v)
        np.save(path_test_fr+"/aimed_pressure.npy", GT_M.p)


        shutil.copy2("parameters.py", path_test_fr+"/parameters.py", follow_symlinks=True)
