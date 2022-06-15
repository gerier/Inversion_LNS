#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 08:43:33 2022

@author: s.gerier
"""

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from copy import deepcopy
import pandas as pd
from scipy.interpolate import CubicSpline
from copy import deepcopy
from discretisation import *

#%% VARIABLES 

class LNS_Variable:
    def __init__(self,rho, v, p):
        self.rho = rho
        self.v = v
        self.p = p
    
    def __mul__(self,alpha):
        self.rho *= alpha
        self.v *= alpha
        self.p *= alpha
        return self

    def __add__(self, other):
        self.rho += other.rho
        self.v += other.v
        self.p += other.p
        return self
    
    def __sub__(self, other):
        self.rho -= other.rho
        self.v -= other.v
        self.p -= other.p
        return self
    
    def __neg__(self):
        self.rho = -self.rho
        self.v = -self.v
        self.p = -self.p
        return self
    
    
#%% EQUATIONS

# Definition of function that describes the Linearised Navier-Stokes 

def DF_sigma_v(V, l, mu, dz, BC, DF_C):
    # compute the divergence of the stress tensor stress tensor
    dV = DF_C(V, dz, BC)
    BC = [0,0,0,0]
    return DF_C( (l+2*mu) * dV, dz, BC, False) 

'''
def DF_kappaDT2(T, kappa, gamma, dz, BC):
    # compute the divergene of the temperaure
    dT = DF_DC_2_backward(T, dz, BC)
    return DF_DC_2_forward(kappa * (gamma-1) * dT, dz, BC) 
'''

def DF_kappaDT(T, kappa, gamma, dz, BC, DF_C):
    # compute the divergene of the temperaure
    dT = DF_C(T, dz, BC)
    aux = np.append( np.append([kappa[0] * (gamma[0]-1)], kappa * (gamma-1)), [kappa[-1] * (gamma[-1]-1)])
    BC_aux = BC * np.array([aux[0],aux[0],aux[-1],aux[-1]])
    return DF_C( (aux[1:]+ aux[:-1])/2  * dT, dz, BC_aux, False) 

def DF_Sigma_C_C(U1, U0, dz, l, gamma, is_reverse, DF_C):
    # initialise
    S  = np.zeros((3,len(U1.rho)))
    
    BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2                       
    BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
    BC_rho0 = [U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]  # Points -2, -1, N+1, N+2 
    BC_p0 = [U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
    
    # compute the contribution on density
    # compute the contribution on qunatity of mvt
    S[1,:-1] = DF_C(U1.p, dz, BC1, False)                                                   # grad p1
    # compute the contribution on pressure
    S[2,:] = gamma * U1.p * DF_C(U0.v, dz, BC_v0)                                  # gamma p1 div v0
    #S[2,:] += interpolation(U1.v) * DF_DC(U0.p, dz, BC_p0, U0.v, is_reverse)       # v1 grad p0
    S[2,:] +=  gamma * U0.p * DF_C(U1.v, dz, BC1)                                  # gamma p0 div v1
    #S[2,:] += interpolation(U0.v) * DF_DC(U1.p, dz, BC1, U0.v, is_reverse)         # v0 grad p1
    return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])



def F(U1, U0, it, source, n, gamma):
    # define source term
    f = np.zeros((3,n))
    source_time = source[0]
    source_spatial = source[1]
    # compute the contribution on density
    for i in range(len(f[0])):
        f[0,i] = source_time[it]/100 * source_spatial[i]
    # compute the contribution on quantity of mvt
    f[1,:-1] =  U1.v * (f[0,1:]+f[0,:-1])/2 
    f[2,:] =   gamma * (U0.p + U1.p) / (U0.rho + U1.rho) * f[0,:] 
    #f[0,:] = 0 # todo change correctly , here because the source is only on velocity
    return LNS_Variable(0*f[0,:], 2 * f[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), f[2,:])


# Definition of resolution spatial
def get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
    if order_DF == 4:
        DF_DC_backward = DF_DC_4_backward
        DF_DC_forward = DF_DC_4_forward
        DF_C = DF_C_4
    else :
        DF_DC_backward = DF_DC_2_backward
        DF_DC_forward = DF_DC_2_forward
        DF_C = DF_C_2
        
    RHS_c = F(U1, U0, int(abs(t/dt)), source, len(U1.rho), gamma)  - DF_Sigma_C_C(U1, U0, dz, l, gamma, is_reverse, DF_C)
    
    RHS_dbackward = RHS_c 
    
    RHS_dforward = RHS_dbackward 
    
    return RHS_dforward

def get_minus_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
    return - get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4)
          
#%% SOURCE

def get_source(t,f0):
    t0 = 1.2/f0
    return np.exp(-4 * np.pi**2 * f0**2 * (t-t0)**2)