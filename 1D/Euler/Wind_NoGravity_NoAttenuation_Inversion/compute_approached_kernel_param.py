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


#%% Attendu 

def chi(M, history_obs, index_receivers, max_obs, dt):
    # init
    CHI = 0
    # get Ucalc
    U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
    _, _, history_calc = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", M, h, dt, source)
    # Ucalc - Uobs
    for it,r in enumerate(index_receivers):
        for t in range(len(history_obs)):
            CHI += (history_calc[t].p[r] - history_obs[t].p[r])**2 / max_obs[2] #[it]

    CHI *= (dt/2)
    return CHI


saving_dir = "BackUps"# _rho0times1p2_approx_deltac_d0p0001"


start_k = int(13/h*1000+nb_index_neg) 
end_k = int(22.5/h*1000+nb_index_neg) 



