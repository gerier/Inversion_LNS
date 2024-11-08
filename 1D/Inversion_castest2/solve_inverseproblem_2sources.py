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
from mpi4py import MPI


import os
import glob
import time

import config
import shutil
import sys

from copy import deepcopy

from Lib_NavierStokesEq.filter import *

comm = MPI.COMM_WORLD
rank = 12#comm.Get_rank()

sub_comm = comm.Split(rank < 300,rank % 300)
print(rank)

type_parametrisation = int(rank / 60) #1
type_gradient = (int(rank / 10)) % 6  +1  # 6
type_regul = 0

if (type_gradient > 6 or type_gradient == 0):
    sys.exit() 

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
    Sc = np.array(Sc)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(np.log(M.rho), np.log(M.p)), M.v)
    m = m / Sc
elif type_parametrisation == 3 :
    parametrisation = ["log_density", "wind", "log_velocity"]
    Sc = [1] * len(z) + [1] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(np.log(M.rho), np.log(c0)), M.v)
    m = m / Sc
elif type_parametrisation == 4 :
    parametrisation = ["log_density", "wind", "velocity"]
    Sc = [1] * len(z) + [1e3] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(np.log(M.rho), c0), M.v)
    m = m / Sc
elif type_parametrisation == 6 :
    parametrisation = ["log_pressure", "wind", "log_velocity"]
    Sc = [1] * len(z) + [1] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append( np.log(M.p),np.log(c0)), M.v)
    m = m / Sc
elif type_parametrisation == 5 :
    parametrisation = ["log_pressure", "wind", "velocity"]
    Sc = [1] * len(z) + [1e3] * len(z) + [1] * (len(z)-1)
    c0 = np.sqrt(M.p * gamma / M.rho)
    m = np.append(np.append(np.log(M.p),c0), M.v)
    m = m / Sc
if rank < 360 : 

    # for the figure
    plt.rcParams.update({'font.size': 16})
    plt.rcParams.update({'figure.autolayout': True})


    # Define the observation that we want to reproduce
    it_end = int(Tmax / dt)
    size = len(z)

    start = time.time()

    # init perturbation field
    U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size))
    # init sponge layer
    sponge_layer_info = get_sponge(size)


    history_obs1, Uobs1 = time_scheme(
        get_RHS, U1, 0, it_end, dt, index_receivers, "forward", sponge_layer_info, GT_M, h, source_set1)
    print(time.time() - start)

    plt.figure()
    tt = np.arange(0,Tmax+dt,dt) #- 12
    plt.plot(tt,[history_obs1[t][0]/1000 for t in range(len(history_obs1))])
    plt.plot(tt,[history_obs1[t][-1]/1000 for t in range(len(history_obs1))])

    plt.figure()
    plt.plot(z/1000,Uobs1.p)
    plt.show()

    start = time.time()

    # init perturbation field
    U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size))
    # init sponge layer
    sponge_layer_info = get_sponge(size)

    history_obs2, Uobs2 = time_scheme(
        get_RHS, U1, 0, it_end, dt, index_receivers, "forward", sponge_layer_info, GT_M, h, source_set2)
    print(time.time() - start)

    history_events = [history_obs1]#,history_obs2]

    # define a dictionary to save information on norms
    # ord: order of the norm: 1,2,inf
    # choice: on time, on receiver (TODO : develop explanation on the options)
    # calc: boolean to know if synthetics are normalised by synthetics or by observations
    norm_info = {'ord': 2, 'choice': 2, 'calc': False}


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

    source_param1 = deepcopy(source_set1)
    source_param2 = deepcopy(source_set2)
    source_events = [source_param1]#,source_param2]


    initial_m = deepcopy(m)

    factor_reg = [0.005,0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 5,10]#[0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5]
    reg_rank = factor_reg[rank % 10]
 

    for factor_reg in [reg_rank] : #[0.005,0.1, 0.5, 1]:#0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 

        # generate a directory to save results
        if rank == 0 and not(os.path.exists("./Resultats/")):
            os.mkdir("./Resultats/")
        sub_comm.Barrier()
        all_test = glob.glob("./Resultats/Test_*/")
        last_test = get_last_test(all_test)
        path_test = "./Resultats/Test_"+str(rank)+"_reg_"+str(factor_reg)
        os.mkdir(path_test)

        save_info = [path_test, plot_flatteniter, plot_gradient, [z, parametrisation, Sc, gamma]]

        # init config parameter
        config.n_iter_alpha = 0
        config.n_fx = 0
        config.n_grad = 0
        config.n_h = 0
        config.is_restarted = 0

        reg = [factor_reg, h, size, deepcopy(initial_m)]


        new_m = deepcopy(initial_m)


        # compute the norm of the observations (and add it to the dictionary)
        norm_events = []
        for hist in history_events: 
            print("je suis la")
            norm_obs = get_norm(hist, norm_info['ord'], norm_info['choice'])
            norm_info["norm_obs"] = norm_obs
            norm_events += [deepcopy(norm_info)]

        # define the arguments of the cost function and its derivative
        argf = [history_events, norm_events, dt, h, size, gamma, g, sponge_layer_info, source_events, index_receivers, Sc,parametrisation]
        argg = [it_end, dt, h, size, sponge_layer_info, source_events, index_receivers,
                history_events, norm_events, gamma, g, Sc,parametrisation]
        argh = [size]  # [len(history_obs), dt, index_receivers, h, gamma, source]

        # treat the path to save results
        path_test_fr = path_test 
        save_info[0] = path_test_fr

        target = interpolation(GT_M.v) + np.sqrt(GT_M.gamma * GT_M.p / GT_M.rho)#np.sqrt(GT_c2)

        # compute the model by the inversion
        new_m, is_notfound, x_iter = optimisation(new_m, f0, df0, hessf0, 0.0001, 0.9, argf, argg, argh, alpha_max=1, maxiter=300,
                                        tol_x=1e-8, tol_df=1e-8, true_hessian=False, regularisation=reg, plot_iter=[plot_flatteniter, [z, parametrisation, Sc, gamma]], save_info=save_info, type_gradient=type_gradient,type_regul=type_regul,target=target)


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
        ax[0].legend(loc=2)
        ax[0].set_title("Comparaison true model and model obtained after inversion")
        ax[3].set_xlabel("Range(km)")
        #ax[0].set_xlim(0, 40)
        #ax[1].set_xlim(0, 40)
        #ax[2].set_xlim(0, 40)
        #ax[3].set_xlim(0, 40)
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
