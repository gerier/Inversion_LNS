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
    dT = DF_C(T, dz, BC, False)
    coef = kappa * (gamma-1)
    coef = (coef[1:]+coef[:-1])/2
    
    BC_aux =  np.array([coef[-2:] * dT[:2], coef[:2] * dT[-2:]]) #np.array([aux[0],aux[0],aux[-1],aux[-1]])
    res = DF_C(  coef * dT, dz, BC_aux)

    return res 


decentered_backward_and_forward = True
if not decentered_backward_and_forward : 
    def DF_Sigma_C(U1, U0, dz, gamma, is_reverse, DF_C, DF_DC):
        # initialise
        S  = np.zeros((3,len(U1.rho)))
        
        BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2                       
        BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
        BC_rho0 = [U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]  # Points -2, -1, N+1, N+2 
        BC_p0 = [U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
        
        # compute the contribution on density
        S[0,:] = U0.rho * DF_C(U1.v, dz, BC1)                                              # rho0 div v1
        S[0,:] += interpolation(U1.v) * DF_DC(U0.rho, dz, BC_rho0, U0.v, is_reverse)       # v1 grad rho0
        S[0,:] +=  U1.rho * DF_C(U0.v, dz, BC_v0)                                          # rho1 div v0
        S[0,:] += interpolation(U0.v) * DF_DC(U1.rho, dz, BC1, U0.v, is_reverse)           # v0 grad rho1 
        # compute the contribution on qunatity of mvt
        S[1,:-1] = (U0.rho[1:]+U0.rho[:-1])/2 * U0.v * DF_DC(U1.v, dz, BC1, U0.v, is_reverse)    # roh0 v0 . div v1
        S[1,:-1] += DF_C(U1.p, dz, BC1, False)                                                   # grad p1
        # compute the contribution on pressure
        S[2,:] = gamma * U1.p * DF_C(U0.v, dz, BC_v0)                                  # gamma p1 div v0
        S[2,:] += interpolation(U1.v) * DF_DC(U0.p, dz, BC_p0, U0.v, is_reverse)       # v1 grad p0
        S[2,:] +=  gamma * U0.p * DF_C(U1.v, dz, BC1)                                  # gamma p0 div v1
        S[2,:] += interpolation(U0.v) * DF_DC(U1.p, dz, BC1, U0.v, is_reverse)         # v0 grad p1
        
        return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])
    
    
    def DF_Sigma_D(U1, U0, dz, T0, g, l, mu, kappa, gamma, Cv, DF_C, DF_DC):
        # initiliase
        S  = np.zeros((3,len(U1.rho)))
        
        BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2
        BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
        
        T1 = (U1.p + U0.p) / ((gamma-1) * Cv * (U1.rho + U0.rho)) - T0
        BC_T1 = BC1 #[T1[0], T1[0], T1[-1], T1[-1]]            # Points -2, -1, N+1, N+2 
        
        sigma_V0 = (2*mu + l) * DF_C(U0.v, dz, BC_v0)
        sigma_Vp = (2*mu + l) * DF_C(U1.v, dz, BC1)
        
        # compute the contribution on density
        S[0,:] = 0
        # compute the contribution on quantity of mvt 
        S[1,:-1] = DF_sigma_v( U1.v, l, mu, dz, BC1, DF_C)       # div Sigma_v v1
        # compute the contribution on pressure
        S[2,:] = DF_kappaDT(T1, kappa, gamma, dz, BC1, DF_C)     # div (kappa (gamma-1) grad T)
        S[2,:] += (gamma-1) * sigma_V0 * DF_C(U1.v, dz, BC1)            # (gamma-1) Sigma_v v0: grad v1 
        S[2,:] += (gamma-1) * sigma_Vp * DF_C(U0.v, dz, BC_v0)          # (gamma-1) Sigma_v v1 : grad v0
    
        return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])
    
    
    def G(U1, U0, dz, g, is_reverse, DF_C, DF_DC):
        # initialise
        GU = np.zeros((3,len(U1.rho)))
        rho1_demi = (U1.rho[1:]+U1.rho[:-1])/2
        g_demi = (g[1:]+g[:-1])/2
        rho0_demi = (U0.rho[1:]+U0.rho[:-1])/2
        
        BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
    
        # compute the contribution on density
        GU[0,:] = 0
        # compute the contribution on quantity of mvt
        GU[1,:-1] = rho1_demi * g_demi                                                                 # rho1 g
        GU[1,:-1] -= (rho0_demi * U1.v + rho1_demi * U0.v) * DF_DC(U0.v, dz, BC_v0, U0.v, is_reverse)  # (rho0v1 + rho1v0) div v0
        # compute the contribution on pressure
        GU[2,:] = 0
        
        return LNS_Variable(GU[0,:], GU[1,:-1] / rho0_demi, GU[2,:])
    
    
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
        f[2,:] =   gamma * U0.p / U0.rho * f[0,:] 
        f[0,:] = 0 # todo change correctly , here because the source is only on velocity
        return LNS_Variable(f[0,:], 2 * f[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), f[2,:])
    
    
    # Definition of resolution spatial
    def get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
        if order_DF == 4:
            DF_DC = DF_DC_4
            DF_C = DF_C_4
        else :
            DF_DC = DF_DC_2
            DF_C = DF_C_2
            
        return F(U1, U0, int(round(2*t/dt)), source, len(U1.rho), gamma) + G(U1, U0, dz, g, is_reverse, DF_C, DF_DC) \
            + DF_Sigma_D(U1, U0, dz, T0, g, l, mu, kappa, gamma, R, DF_C, DF_DC) \
            - DF_Sigma_C(U1, U0, dz, gamma, is_reverse, DF_C, DF_DC)
       
    def get_minus_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
        return - get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4)
else : 
    
    def DF_Sigma_C_C(U1, U0, dz, gamma, is_reverse, DF_C):
        # initialise
        S  = np.zeros((3,len(U1.rho)))
        
        #BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2                       
        #BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
        #BC_rho0 = [U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]  # Points -2, -1, N+1, N+2 
        #BC_p0 = [U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
        BC_v1 = [U1.v[-2], U1.v[-1], U1.v[0], U1.v[1]]
        BC_v0 = [U0.v[-2], U0.v[-1], U0.v[0], U0.v[1]]  #[U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
        BC_p1 = [U1.p[-2], U1.p[-1], U1.p[0], U1.p[1]]
        
        # compute the contribution on density
        S[0,:] = U0.rho * DF_C(U1.v, dz, BC_v1)                                              # rho0 div v1
        S[0,:] +=  U1.rho * DF_C(U0.v, dz, BC_v0)                                          # rho1 div v0
        # compute the contribution on qunatity of mvt
        S[1,:-1] = DF_C(U1.p, dz, BC_p1, False)                                                   # grad p1
        # compute the contribution on pressure
        S[2,:] = gamma * U1.p * DF_C(U0.v, dz, BC_v0)                                  # gamma p1 div v0
        S[2,:] +=  gamma * U0.p * DF_C(U1.v, dz, BC_v1)                                  # gamma p0 div v1
        
        return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])

    def DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC):
        # initialise
        S  = np.zeros((3,len(U1.rho)))
        
        #BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2                       
        BC_v1 = [U1.v[-2], U1.v[-1], U1.v[0], U1.v[1]]
        BC_rho0 = [U0.rho[-2], U0.rho[-1], U0.rho[0], U0.rho[1]] #[U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]  # Points -2, -1, N+1, N+2 
        BC_p0 = [U0.p[-2], U0.p[-1], U0.p[0], U0.p[1]]    #[U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
        BC_p1 = [U1.p[-2], U1.p[-1], U1.p[0], U1.p[1]]
        
        # compute the contribution on density
        S[0,:] = interpolation(U1.v) * DF_DC(U0.rho, dz, BC_rho0, U0.v, is_reverse)        # v1 grad rho0
        S[0,:] += interpolation(U0.v) * DF_DC(U1.rho, dz, BC_v1, U0.v, is_reverse)           # v0 grad rho1 
        # compute the contribution on qunatity of mvt
        S[1,:-1] = (U0.rho[1:]+U0.rho[:-1])/2 * U0.v * DF_DC(U1.v, dz, BC_v1, U0.v, is_reverse)    # roh0 v0 . div v1
        # compute the contribution on pressure
        S[2,:] = interpolation(U1.v) * DF_DC(U0.p, dz, BC_p0, U0.v, is_reverse)       # v1 grad p0
        S[2,:] += interpolation(U0.v) * DF_DC(U1.p, dz, BC_p1, U0.v, is_reverse)         # v0 grad p1

        return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])
    
    def DF_Sigma_D_C(U1, U0, dz, T0, g, l, mu, kappa, gamma, Cv, DF_C):
        # initiliase
        S  = np.zeros((3,len(U1.rho)))
        
        BC_v1 = [U1.v[-2], U1.v[-1], U1.v[0], U1.v[1]]  #BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2
        BC_v0 = [U0.v[-2], U0.v[-1], U0.v[0], U0.v[1]]  #BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
        
        T1 = (U1.p + U0.p) / ((gamma-1) * Cv * (U1.rho + U0.rho)) - T0
        BC_T1 = [T1[-2], T1[-1], T1[0], T1[1]]  #BC1 #[T1[0], T1[0], T1[-1], T1[-1]]            # Points -2, -1, N+1, N+2 
        
        sigma_V0 = (2*mu + l) * DF_C(U0.v, dz, BC_v0)
        sigma_Vp = (2*mu + l) * DF_C(U1.v, dz, BC_v1)
        
        # compute the contribution on density
        S[0,:] = 0
        # compute the contribution on quantity of mvt 
        S[1,:-1] = DF_sigma_v( U1.v, l, mu, dz, BC_v1, DF_C)       # div Sigma_v v1
        # compute the contribution on pressure
        S[2,:] = DF_kappaDT(T1, kappa, gamma, dz, BC_T1, DF_C)     # div (kappa (gamma-1) grad T)
        S[2,:] += (gamma-1) * sigma_V0 * DF_C(U1.v, dz, BC_v1)            # (gamma-1) Sigma_v v0: grad v1 
        S[2,:] += (gamma-1) * sigma_Vp * DF_C(U0.v, dz, BC_v0)          # (gamma-1) Sigma_v v1 : grad v0
    
        return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])
    
    
    def G_C(U1, U0, dz, g, is_reverse):
        # initialise
        GU = np.zeros((3,len(U1.rho)))
        rho1_demi = (U1.rho[1:]+U1.rho[:-1])/2
        g_demi = (g[1:]+g[:-1])/2
        rho0_demi = (U0.rho[1:]+U0.rho[:-1])/2
        
        # compute the contribution on density
        GU[0,:] = 0
        # compute the contribution on quantity of mvt
        GU[1,:-1] = rho1_demi * g_demi                                                                 # rho1 g
        # compute the contribution on pressure
        GU[2,:] = 0
        
        return LNS_Variable(GU[0,:], GU[1,:-1] / rho0_demi, GU[2,:])
    
    def G_DC(U1, U0, dz, g, is_reverse, DF_DC):
        # initialise
        GU = np.zeros((3,len(U1.rho)))
        rho1_demi = (U1.rho[1:]+U1.rho[:-1])/2
        rho0_demi = (U0.rho[1:]+U0.rho[:-1])/2
        
        BC_v0 = [U0.v[-2], U0.v[-1], U0.v[0], U0.v[1]]    #[U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
    
        # compute the contribution on density
        GU[0,:] = 0
        # compute the contribution on quantity of mvt
        GU[1,:-1] = - (rho0_demi * U1.v + rho1_demi * U0.v) * DF_DC(U0.v, dz, BC_v0, U0.v, is_reverse)  # (rho0v1 + rho1v0) div v0
        # compute the contribution on pressure
        GU[2,:] = 0
        
        return LNS_Variable(GU[0,:], GU[1,:-1] / rho0_demi, GU[2,:])
    
    
    def F(U1, U0, it, source, n, gamma):
        # define source term
        f = np.zeros((3,n))
        source_time = source[0]
        source_spatial = source[1]
        # compute the contribution on density
        for i in range(len(f[0])):
            f[0,i] = source_time[it]/100 * source_spatial[i]
        # compute the contribution on quantity of mvt
        f[1,:-1] = U1.v * (f[0,1:]+f[0,:-1])/2 
        f[2,:] =  gamma * U0.p / U0.rho * f[0,:] 
        #f[0,:] = 0 # todo change correctly , here because the source is only on velocity
        return LNS_Variable(f[0,:], 2 * f[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), f[2,:])
    
    
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
            
        RHS_c = F(U1, U0, int(np.round(2*t/dt)), source, len(U1.rho), gamma) + G_C(U1, U0, dz, g, is_reverse) \
            + DF_Sigma_D_C(U1, U0, dz, T0, g, l, mu, kappa, gamma, R, DF_C) \
            - DF_Sigma_C_C(U1, U0, dz, gamma, is_reverse, DF_C)

        if not is_reverse : 
            RHS_plusdt = RHS_c + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_4) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_4)) #*0.5
            #RHS_plusdt = RHS_plusdt2 + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_forward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_forward))*0.5         
        else : 
            RHS_plusdt = RHS_c + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_4) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_4)) #*0.5 
            #RHS_plusdt = RHS_plusdt2 + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_backward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_backward))*0.5
        return RHS_plusdt

    def get_minus_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
        return - get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4)
          
#%% SOURCE

def get_source(t,f0):
    t0 = 1.2/f0
    return np.exp(-4 * np.pi**2 * f0**2 * (t-t0)**2)


#%% ADJOINT EQUAITONS

def dchi(reverse_U, obs_U, t, dt, dz, d_receivers, max_obs):
    n = len(reverse_U[0].rho)
    dr = int(d_receivers / dz)
    diff = np.zeros((3, n))
    
    reverse_U_t = reverse_U[len(reverse_U) - 1 - int(np.round(t/dt))]
    obs_U_t = obs_U[int(round(t/dt))]
    diff[0,0:n:dr] = (reverse_U_t.rho[0:n:dr] - obs_U_t.rho[0:n:dr]) / max_obs[0]
    diff[1,0:n-1:dr] = (reverse_U_t.v[0:n:dr] - obs_U_t.v[0:n:dr]) / max_obs[1]
    diff[2,0:n:dr] = (reverse_U_t.p[0:n:dr] - obs_U_t.p[0:n:dr]) / max_obs[2]
    
    return diff
    
def DF_kappaDp_adjoint(p, kappa, gamma, dz, BC, DF_C, DF_DC):
    # compute the divergene of the temperaure
    dp = DF_C(p, dz, BC)
    aux = np.append( np.append([kappa[0] * (gamma[0]-1)], kappa * (gamma-1)), [kappa[-1] * (gamma[-1]-1)])
    BC_auxdp = np.zeros(4)
    return DF_C( (aux[1:]+ aux[:-1])/2  * dp, dz, BC_auxdp, False) # TODO pb BC 


def get_adjoint_RHS(Ustar, t, U0, T0, g, l, mu, kappa, gamma, Cv, dz, dt, source, reverse_U, observation_U, d_receivers, max_obs, order_DF=4):
    if order_DF == 4:
        DF_DC = DF_DC_4
        DF_C = DF_C_4
    else :
        DF_DC = DF_DC_2
        DF_C = DF_C_2
        
    U = np.zeros((3, len(Ustar.rho)))
    
    
    BCstar = [0,0,0,0]
    BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]
    BC_rho0 = [U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]
    BC_p0 = [U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]
    
    ddp = DF_kappaDp(Ustar.p, kappa, gamma, dz, [0,*BCstar,0], DF_C, DF_DC)
    sigma_v0 = DF_C(U0.v, dz, BC_v0)
    
    Frho1 = F(U1, U0, int(t/dt), source, len(Ustar.rho)).rho
    Fadjoint = dchi(reverse_U, observation_U, t, dt, dz, d_receivers, max_obs)
    
    U[0,:] = Ustar.rho * DF_C(U0.v, dz, BC_v0)                                          # rho* div v0
    U[0,:] -= DF_C(  (Ustar.rho[1:]+Ustar.rho[:-1])/2 * U0.v, dz, BCstar)               # div (rho* v0) 
    U[0,:] -= g * interpolation(Ustar.v)                                                # g . m*
    U[0,:] +=  interpolation(U0.v) * DF_C(U0.v, dz, BC_v0) * interpolation(Ustar.v)    # (v0 . div v0) . m*
    U[0,:] -= T0 / U0.rho * ddp                                                         # To/rho0 div( (1-gamma) kappa grad p)
    U[0,:] += Fadjoint[0,:]
    
    U[1,:-1] = (Ustar.rho[1:]+Ustar.rho[:-1])/2 * DF_C(U0.rho, dz, BC_rho0, False)      # rho* grad(rho0)
    U[1,:-1] -= DF_C( U0.rho * Ustar.rho, dz, BCstar, False)                            # div (rho0 rho*)
    #U[1,:-1] -= DF_C( U0.rho * interpolation(U0.v * Ustar.v), dz, BCstar, False)        # div(rho0 v0 otimes m*)
    U[1,:-1] -= (U0.rho[1:] + U0.rho[:-1])/2 * U0.v * DF_DC( Ustar.v, dz, BCstar, U0.v, False)        # div(rho0 v0 otimes m*)
    U[1,:-1] -= DF_sigma_v( Ustar.v, l, mu, dz, BCstar, DF_C, DF_DC)                    # div (Sigma (m*))
    U[1,:-1] += (U0.rho[1:]+U0.rho[:-1])/2 * DF_DC(U0.v, dz, BC_v0, U0.v) * Ustar.v     # rho0 div (v0) m*
    U[1,:-1] -= (Frho1[1:]+Frho1[:-1])/2 *  Ustar.v                                     # F m* 
    U[1,:-1] += (Ustar.p[1:]+Ustar.p[:-1])/2 * DF_C(U0.p, dz, BC_p0, False)             # p* div p0
    U[1,:-1] -= DF_C( gamma * U0.p * Ustar.p, dz, BCstar, False)                        # grad ( gamma p* p0)
    U[1,:-1] += 2 * DF_C(Ustar.p * (gamma-1) * sigma_v0, dz, BCstar, False)             # 2 div ( Sigma_v (v0) (gamma-1) p* )
    U[1,:-1] += Fadjoint[0,:-1]
    
    U[2,:] = - DF_C(Ustar.v, dz, BCstar)                                       # div m*
    U[2,:] += gamma * Ustar.p * DF_C(U0.v, dz, BC_v0)                          # gamma p* div v0
    U[2,:] -= DF_C(  (Ustar.p[1:]+Ustar.p[:-1])/2 * U0.v, dz, BCstar)          # div(p* v0)
    U[2,:] += ddp / ((gamma-1) * Cv * U0.rho)                                  # 1/ ( (gamma-1) cv rho0) div( (1-gamma) kappa grad p*)
    U[2,:] += Fadjoint[2,:]
    
    return LNS_Variable(U[0,:], 2 * U[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), U[2,:]) 


