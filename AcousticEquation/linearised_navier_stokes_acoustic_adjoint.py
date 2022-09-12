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
    
    def __truediv__(self, alpha):
        self.rho = self.rho / alpha
        self.v =  self.v / alpha
        self.p = self.p  / alpha
        return self
    
    def __neg__(self):
        self.rho = -self.rho
        self.v = -self.v
        self.p = -self.p
        return self
    
    def plot(self, abs, title):
        fig,ax = plt.subplots(3,1, figsize=(10,7))
        ax[0].plot(abs,self.rho)
        ax[1].plot(abs,self.p)
        ax[2].plot( (abs[1:] + abs[:-1])/2,self.v)
        ax[0].set_xlabel("Distance on axis x (km)")
        ax[1].set_xlabel("Distance on axis x (km)")
        ax[2].set_xlabel("Distance on axis x (km)")
        ax[0].set_ylabel("Density")
        ax[1].set_ylabel("Pressure")
        ax[2].set_ylabel("Velocity")
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        fig.suptitle(title)


class LNS_Model:
    def __init__(self,rho, c):
        self.rho = rho
        self.c = c
    
    def __mul__(self,alpha):
        self.rho *= alpha
        self.c *= alpha
        return self

    def __add__(self, other):
        self.rho += other.rho
        self.c += other.c
        return self
    
    def __sub__(self, other):
        self.rho -= other.rho
        self.c -= other.c
        return self
    
    def __truediv__(self, alpha):
        self.rho = self.rho / alpha
        self.c =  self.c / alpha
        return self
    
    def __neg__(self):
        self.rho = -self.rho
        self.v = -self.v
        self.c = -self.c
        return self
    
    def plot(self, abs, title):
        fig,ax = plt.subplots(2,1, figsize=(10,7))
        ax[0].plot(abs,self.rho)
        ax[1].plot(abs,self.c)
        ax[0].set_xlabel("Distance on axis x (km)")
        ax[1].set_xlabel("Distance on axis x (km)")
        ax[0].set_ylabel("Density")
        ax[1].set_ylabel("Wave velocity")
        ax[0].grid()
        ax[1].grid()
        fig.suptitle(title)



#%% EQUATIONS

# Definition of function that describes the Linearised Navier-Stokes 

if True : 
    
    def DF_Sigma_C_C(U1, U0, dz, gamma, is_reverse, DF_C):
        # initialise
        S  = np.zeros((3,len(U1.rho)))
        # boundary conditions
        BC_v1 = [U1.v[-2], U1.v[-1], U1.v[0], U1.v[1]]
        BC_p1 = [U1.p[-2], U1.p[-1], U1.p[0], U1.p[1]]
        # compute the contribution on qunatity of mvt
        S[1,:-1] = DF_C(U1.p, dz, BC_p1, False)                                                   # grad p1
        # compute the contribution on pressure
        S[2,:] = U0.rho * U0.c**2 * DF_C(U1.v, dz, BC_v1)                                  # gamma p1 div v0
        return S

    
    
    def F(z, t, f0, index_source):
        # initialisation
        f = np.zeros((3,len(z)))
        # define source term
        source = get_source_t_x(z, t, f0, index_source)
        # compute the contribution on quantity of mvt
        f[1,:-1] = (source[1:]+source[:-1])/2 
        return f

    
    # Definition of resolution spatial
    def get_RHS(U1, t, previous_U, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
        # initialisation
        if order_DF == 4:
            DF_C = DF_C_4
        else :
            DF_C = DF_C_2
        z, f0, index_source = source 
        # computation of the right hand side of the LNS equation
        RHS_c = F(z,t,f0,index_source) \
                - DF_Sigma_C_C(U1, U0, dz, gamma, is_reverse, DF_C)
        # Change the tye of RHS_c to LNS_Variable
        RHS_c = LNS_Variable(RHS_c[0,:], 2 * RHS_c[1,:-1] / (U0.rho[1:] + U0.rho[:-1]), RHS_c[2,:])
        return previous_U + RHS_c * dt

    def get_minus_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
        return - get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4) * dt
          
#%% SOURCE

def get_source(t,f0):
    t0 = 1.2/f0
    return np.exp(-4 * np.pi**2 * f0**2 * (t-t0)**2)

def get_source_t_x(z, t,f0, index_source):
    # define source time function
    time_source = get_source(t,f0) / 100 
    # define the locality of the source
    spatial_source = z == z[index_source]
    # create source from temporal and spatial info
    source = time_source * spatial_source
    return source


#%% ADJOINT EQUAITONS

def dchi(reverse_U, obs_U, t, dt, index_receivers, max_obs):
    n = len(reverse_U[0].rho)
    diff = np.zeros((3, n))

    reverse_U_t = reverse_U[int(np.round(t/dt))]
    obs_U_t = obs_U[len(reverse_U) - 1 - int(round(t/dt))]
    
    diff[0,index_receivers] = (reverse_U_t.rho[index_receivers] - obs_U_t.rho[index_receivers]) / max_obs[0]
    diff[1,index_receivers] = (reverse_U_t.v[index_receivers] - obs_U_t.v[index_receivers]) / max_obs[1]
    diff[2,index_receivers] = (reverse_U_t.p[index_receivers] - obs_U_t.p[index_receivers]) / max_obs[2]
    
    return diff 

 
def get_adjoint_RHS(Ustar, t, previous_Ustar, U0, T0, g, l, mu, kappa, gamma, Cv, dz, dt, source, reverse_U, observation_U, index_receivers, max_obs, order_DF=4):
    if order_DF == 4:
        DF_C = DF_C_4
    else :
        DF_C = DF_C_2
        
    U = np.zeros((3, len(Ustar.rho)))
    
    BC_pstar = [Ustar.p[-2], Ustar.p[-1], Ustar.p[0], Ustar.p[1]]
    BC_vstar = [Ustar.v[-2], Ustar.v[-1], Ustar.v[0], Ustar.v[1]]
    
    Fadjoint = dchi(reverse_U, observation_U, t, dt, index_receivers, max_obs)

    U[1,:-1] += DF_C(Ustar.p, dz, BC_pstar, False)                        # grad ( gamma p* p0)
    # U[1,:-1] += Fadjoint[1,:-1]
    
    U[2,:] += U0.rho * U0.c**2 * DF_C(Ustar.v, dz, BC_vstar)                          # gamma p* div v0
    U[2,:] += Fadjoint[2,:] * (U0.rho * U0.c**2)
    
    return previous_Ustar + LNS_Variable(U[0,:], 2 * U[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), U[2,:]) * dt 


#%% KERNELS

def get_kernels_centered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dtvp, dz, DF_C):
    # initialization
    K = np.zeros((3, len(rho_a)))
    
    # main boundary conditions
    BC_vp = [v_p[-2], v_p[-1], v_p[0], v_p[1]]

    # kernel in rho0
    K[0,:] += interpolation(dtvp * v_a)    
    K[0,:] -= DF_C(v_p, dz, BC_vp) * p_a / U0.rho
    
    # kernel in c0
    K[1,:] = - 2 * DF_C(v_p, dz, BC_vp) * p_a / U0.c

    return K
    


def get_kernels(hist_adjoint, hist_backprop, U0, T0, kappa, gamma, Cv, dt, dz, order_DF=4):
    # initialize
    K = np.zeros((3, len(hist_adjoint[0].rho)))
    if order_DF == 4:
        DF_C = DF_C_4
    else :
        DF_C = DF_C_2

    for t in range(len(hist_adjoint)-1):
        # create variable to be more lisible 
        rho_a = (hist_adjoint[t].rho + hist_adjoint[t+1].rho)/2    # mean at time n+1/2
        rho_p = (hist_backprop[t].rho + hist_backprop[t+1].rho)/2  # mean at time n+1/2
        
        p_a = (hist_adjoint[t].p + hist_adjoint[t+1].p)/2          # mean at time n+1/2
        p_p = (hist_backprop[t].p + hist_backprop[t + 1].p)/2      # mean at time n+1/2 
        
        v_a = hist_adjoint[t].v                                   # at time n
        v_p = hist_backprop[t].v                                  # at time n
        dtvp = (hist_backprop[t +1].v - v_p) / dt

        K += dt * get_kernels_centered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dtvp, dz, DF_C)
       
    return K
