#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:02:52 2022

@author: s.gerier
"""

import sys
sys.path.insert(0, './Lib/')


from parameters import *
from discretisation import *
from linearised_navier_stokes_acoustic_adjoint import *


import numpy as np
import time


# INITIALISATION
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 


# RESOLUTION
start_time = time.time()
t_end, U_end, history_obs = time_scheme(get_RHS, U1, T_init, Tmax, z, "forward", GT_U0, T0, g, 0, 0, 0, gamma, Cv, h, dt, source, False)
print(" Total time : ", time.time() - start_time)

np.save("./BackUps/observation_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers]), np.array(history_obs))
