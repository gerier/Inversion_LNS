#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import *
from parameters import *
from get_norm import * 
from parametrisation import *

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


#%% Attendu 


saving_dir = "BackUps"


start_k = int(3/h*1000) 
end_k = int(40/h*1000+nb_index_neg) 



def chi(m0,obs,norm_info,dt,h,size,gamma,sponge_layer,source,index_receivers, Sc,parametrisation):

    # set the model in the parametrisation ["density", "wind","pressure"] to be adapted to implementation of Euler equation
    M = m0

    # init
    CHI = 0
    # get Ucalc
    U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size)) 
    calc,_ = time_scheme(get_RHS, U1, 0, len(obs), dt, index_receivers, "forward", sponge_layer, M, h, source)
    # Ucalc - Uobs / norm
    if norm_info['choice'] == 1 :
        # norm is sum_rec( sum_t (observations(t,r))) or max_rec( max_t (observations(t,r)))
        for ir in range(len(index_receivers)) :
            norm_obs = norm_info["norm_obs"]
            for it in range(len(obs)):
                damping = apodisation(it*dt, 1.2/source[1],source[2])
                CHI += (calc[it][ir] - obs[it][ir])**2 * damping / norm_obs**2

    elif norm_info['choice'] == 2 :
        if norm_info["calc"]:
            # norm is sum_t (signal(t,r))) or max_t (signal(t,r))), where signal is for synthetics or observations
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for ir in range(len(index_receivers)) :
                for it in range(len(obs)):
                    damping = apodisation(it*dt, 1.2/source[1],source[2])
                    CHI += (calc[it][ir]/norm_calc[ir] - obs[it][ir]/norm_obs[ir])**2 
        else : 
            # norm is sum_t (observations(t,r))) or max_t (observations(t,r))),
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for ir in range(len(index_receivers)) :
                for it in range(len(obs)):
                    damping = apodisation(it*dt, 1.2/source[1],source[2])
                    CHI += (calc[it][ir]- obs[it][ir])**2/norm_obs[ir]**2

    elif norm_info['choice'] == 3 :
        if norm_info["calc"]:
            # norm is  (signal(t,r))), where signal is for synthetics or observations
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for ir in range(len(index_receivers)) :
                for it in range(len(obs)):
                    damping = apodisation(it*dt, 1.2/source[1],source[2])
                    CHI += (calc[it][ir]/norm_calc[it][ir] - obs[it][ir]/norm_obs[ir])**2 
        else : 
            # norm is (observations(t,r)))
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for ir in range(len(index_receivers)) :
                for it in range(len(obs)):
                    damping = apodisation(it*dt, 1.2/source[1],source[2])
                    CHI += ((calc[it][ir] - obs[it][ir])/norm_obs[it][ir])**2 

    CHI *= (dt/2)
    return CHI