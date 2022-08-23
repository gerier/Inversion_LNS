import sys
sys.path.insert(1, './Lib/')
sys.path.insert(1, './AcousticEquation')

from discretisation import *
from linearised_navier_stokes_acoustic_adjoint import *
from parameters import *

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


K = np.load("./BackUps/kernel_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(z[index_receivers])+".npy", allow_pickle=True)


