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
        BC_rho1 = [U1.rho[-2], U1.rho[-1], U1.rho[0], U1.rho[1]] 
        BC_p0 = [U0.p[-2], U0.p[-1], U0.p[0], U0.p[1]]    #[U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
        BC_p1 = [U1.p[-2], U1.p[-1], U1.p[0], U1.p[1]]
        
        # compute the contribution on density
        S[0,:] = interpolation(U1.v) * DF_DC(U0.rho, dz, BC_rho0, U0.v, is_reverse)        # v1 grad rho0
        S[0,:] += interpolation(U0.v) * DF_DC(U1.rho, dz, BC_rho1, U0.v, is_reverse)           # v0 grad rho1 
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
        for i in range(n):
            f[0,i] = source_time[it]/100 * source_spatial[i]
        # compute the contribution on quantity of mvt
        f[1,:-1] = (f[0,1:]+f[0,:-1])/2 # U1.v * 
        f[2,:] =  0*gamma * U0.p / U0.rho * f[0,:] 
        
        return LNS_Variable(0*f[0,:], 2 * f[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), f[2,:])
    
    
    # Definition of resolution spatial
    def get_RHS(U1, t, previous_U, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
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
        

        if False : 
            if not is_reverse : 
                RHS_plusdt = RHS_c + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_4) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_4)) #*0.5
                #RHS_plusdt = RHS_plusdt2 + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_forward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_forward))*0.5         
            else : 
                RHS_plusdt = RHS_c + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_4) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_4)) #*0.5 
                #RHS_plusdt = RHS_plusdt2 + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_backward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_backward))*0.5
            RHS_plusdt = previous_U + RHS_plusdt * dt
        else : 
            if not is_reverse : 
                RHS_demidt = RHS_c + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_backward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_backward)) *0.5
                Udemi  = previous_U +  RHS_demidt * dt
                RHS_plusdt = Udemi + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_forward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_forward))*0.5 * dt      
            else :    
                RHS_demidt = RHS_c + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_forward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_forward)) *0.5 
                Udemi  = previous_U +  RHS_demidt * dt
                RHS_plusdt = Udemi + (G_DC(U1, U0, dz, g, is_reverse, DF_DC_backward) - DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_backward))*0.5*dt
            
        return RHS_plusdt 

    def get_minus_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4):
        return - get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source, is_reverse, order_DF=4) * dt
          
#%% SOURCE

def get_source(t,f0):
    t0 = 1.2/f0
    return np.exp(-4 * np.pi**2 * f0**2 * (t-t0)**2)


#%% ADJOINT EQUAITONS

def dchi(reverse_U, obs_U, t, dt, index_receivers, max_obs):
    n = len(reverse_U[0].rho)
    diff = np.zeros((3, n))

    reverse_U_t = reverse_U[len(reverse_U) - 1 - int(np.round(t/dt))]
    obs_U_t = obs_U[int(round(t/dt))]
    
    diff[0,index_receivers] = (reverse_U_t.rho[index_receivers] - obs_U_t.rho[index_receivers]) / max_obs[0]
    diff[1,index_receivers] = (reverse_U_t.v[index_receivers] - obs_U_t.v[index_receivers]) / max_obs[1]
    diff[2,index_receivers] = (reverse_U_t.p[index_receivers] - obs_U_t.p[index_receivers]) / max_obs[2]
    
    return diff

 
def DF_kappaDp_adjoint(p, kappa, gamma, dz, BC, DF_C):
    # compute the divergene of the temperaure
    dp = DF_C(p, dz, BC, False)
    coef = kappa * (gamma-1)
    coef = (coef[1:]+coef[:-1])/2
    
    BC_aux =  np.array([coef[-2:] * dp[:2], coef[:2] * dp[-2:]]) #np.array([aux[0],aux[0],aux[-1],aux[-1]])
    res = DF_C(  coef * dp, dz, BC_aux)
    return res 


def get_adjoint_RHS(Ustar, t, previous_Ustar, U0, T0, g, l, mu, kappa, gamma, Cv, dz, dt, source, reverse_U, observation_U, index_receivers, max_obs, order_DF=4):
    if order_DF == 4:
        DF_DC = DF_DC_4
        DF_C = DF_C_4
    else :
        DF_DC = DF_DC_2
        DF_C = DF_C_2
        
    U = np.zeros((3, len(Ustar.rho)))
    
    BC_vstar = [Ustar.v[-2], Ustar.v[-1], Ustar.v[0], Ustar.v[1]]
    BC_v0= [U0.v[-2], U0.v[-1], U0.v[0], U0.v[1]]
    BC_rho0 = [U0.rho[-2], U0.rho[-1], U0.rho[0], U0.rho[1]] #[U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]  # Points -2, -1, N+1, N+2 
    BC_rhostar = [Ustar.rho[-2], Ustar.rho[-1], Ustar.rho[0], Ustar.rho[1]]
    BC_p0 = [U0.p[-2], U0.p[-1], U0.p[0], U0.p[1]]    #[U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
    BC_pstar = [Ustar.p[-2], Ustar.p[-1], Ustar.p[0], Ustar.p[1]]
    BC_gamma = [gamma[-2], gamma[-1], gamma[0], gamma[1]]
    
    BC_ddp = [Ustar.p[-2], Ustar.p[-1], Ustar.p[0], Ustar.p[1]]
    ddp = DF_kappaDp_adjoint(Ustar.p, kappa, gamma, dz, BC_ddp, DF_C)
    sigma_v0 = DF_C(U0.v, dz, BC_v0)
    BC_sigmav0 = [sigma_v0[-2], sigma_v0[-1], sigma_v0[0], sigma_v0[1]]
    
    Frho1 = 0*F(Ustar, U0, int(t/dt), source, len(Ustar.rho), gamma).rho
    Fadjoint = dchi(reverse_U, observation_U, t, dt, index_receivers, max_obs)
    
    U[0,:] = Ustar.rho * DF_C(U0.v, dz, BC_v0)                                          # rho* div v0
    U[0,:] -= DF_C(  (Ustar.rho[1:]+Ustar.rho[:-1])/2 * U0.v, dz, [a * b for a, b in zip(BC_rhostar,BC_v0)])               # div (rho* v0) 
    U[0,:] -= g * interpolation(Ustar.v)                                                # g . m*
    U[0,:] +=  interpolation(U0.v) * DF_C(U0.v, dz, BC_v0) * interpolation(Ustar.v)    # (v0 . div v0) . m*
    U[0,:] -= T0 / U0.rho * ddp                                                         # To/rho0 div( (1-gamma) kappa grad p)
    U[0,:] += Fadjoint[0,:]
    
    U[1,:-1] = (Ustar.rho[1:]+Ustar.rho[:-1])/2 * DF_C(U0.rho, dz, BC_rho0, False)      # rho* grad(rho0)
    U[1,:-1] -= DF_C( U0.rho * Ustar.rho, dz, [a * b for a, b in zip(BC_rho0,BC_rhostar)] , False)                            # div (rho0 rho*)
    #U[1,:-1] -= DF_C( U0.rho * interpolation(U0.v * Ustar.v), dz, BCstar, False)        # div(rho0 v0 otimes m*)
    U[1,:-1] -= (U0.rho[1:] + U0.rho[:-1])/2 * U0.v * DF_DC( Ustar.v, dz, BC_vstar, U0.v, False)        # div(rho0 v0 otimes m*)
    U[1,:-1] -= DF_sigma_v( Ustar.v, l, mu, dz, BC_vstar, DF_C)                    # div (Sigma (m*))
    U[1,:-1] += (U0.rho[1:]+U0.rho[:-1])/2 * DF_DC(U0.v, dz, BC_v0, U0.v) * Ustar.v     # rho0 div (v0) m*
    U[1,:-1] -= (Frho1[1:]+Frho1[:-1])/2 *  Ustar.v                                     # F m* 
    U[1,:-1] += (Ustar.p[1:]+Ustar.p[:-1])/2 * DF_C(U0.p, dz, BC_p0, False)             # p* div p0
    U[1,:-1] -= DF_C( gamma * U0.p * Ustar.p, dz, [a * b * c for a,b,c in zip(BC_pstar,BC_p0,BC_gamma)], False)                        # grad ( gamma p* p0)
    U[1,:-1] += 2 * DF_C(Ustar.p * (gamma-1) * sigma_v0, dz, [a * (b-1) * c for a, b, c in zip(BC_pstar, BC_gamma, BC_sigmav0)], False)             # 2 div ( Sigma_v (v0) (gamma-1) p* )
    U[1,:-1] += Fadjoint[1,:-1]
    
    U[2,:] = - DF_C(Ustar.v, dz, BC_vstar)                                       # div m*
    U[2,:] += gamma * Ustar.p * DF_C(U0.v, dz, BC_v0)                          # gamma p* div v0
    U[2,:] -= DF_C(  (Ustar.p[1:]+Ustar.p[:-1])/2 * U0.v, dz, BC_pstar)          # div(p* v0)
    U[2,:] += ddp / ((gamma-1) * Cv * U0.rho)                                  # 1/ ( (gamma-1) cv rho0) div( (1-gamma) kappa grad p*)
    U[2,:] += Fadjoint[2,:]
    
    return previous_Ustar + LNS_Variable(U[0,:], 2 * U[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), U[2,:]) * dt 


#%% KERNELS

def get_kernels_centered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dtvp, dz, DF_C):
    # initialization
    K = np.zeros((3, len(rho_a)))
    
    # main boundary conditions
    BC_rhoa = [rho_a[-2], rho_a[-1], rho_a[0], rho_a[1]]
    BC_rhop = [rho_p[-2], rho_p[-1], rho_p[0], rho_p[1]]
    BC_rho0 = [U0.rho[-2], U0.rho[-1], U0.rho[0], U0.rho[1]]
    BC_pa = [p_a[-2], p_a[-1], p_a[0], p_a[1]]
    BC_pp = [p_p[-2], p_p[-1], p_p[0], p_p[1]]
    BC_v0 = [U0.v[-2], U0.v[-1], U0.v[0], U0.v[1]]            # Points -2, -1, N+1, N+2 
    BC_vp = [v_p[-2], v_p[-1], v_p[0], v_p[1]]
    
    BC_gamma = [gamma[-2], gamma[-1], gamma[0], gamma[1]] 
    BC_ddp = [p_a[-2], p_a[-1], p_a[0], p_a[1]]
    
    BC_gamma_pa_pp = [a * b * c for a, b, c in zip(BC_gamma, BC_pa,BC_pp)]
    BC_rhoa_rhop = [a * b for a, b in zip(BC_rhoa,BC_rhop)]
    
    # other terms
    T_p = (p_p + U0.p) / ((gamma-1) * Cv * (rho_p + U0.rho)) - T0
    sigma_vp = DF_C(v_p, dz, BC_vp)
    BC_sigma_vp = [sigma_vp[-2], sigma_vp[-1], sigma_vp[0], sigma_vp[1]]
    ddp = DF_kappaDp_adjoint(p_a, kappa, gamma, dz, BC_ddp, DF_C)
    
    # kernel in rho0
    #K[0,:] += rho_a * DF_C(v_p, dz, BC_vp)                         # rho* div (v')
    #K[0,:] -= rho_a * DF_C(v_p, dz, BC_vp)                        # rho* div ( v')
    K[0,:] += interpolation(U0.v) * DF_C(v_p, dz, BC_vp) * interpolation(v_a)         # (v0 . nabla) v' m*
    K[0,:] += interpolation(v_p) * DF_C(U0.v, dz, BC_v0) * interpolation(v_a)         # (v' . nabla) v0 m*
    K[0,:] += ((U0.p * rho_p) / (Cv * (gamma-1) * U0.rho**3) - T_p / U0.rho ) * ddp   # ((p0 rho')/(cv (gamma-1)rho0³)-T'/rho0) *  div( (1-gamma) kappa grad p*)
    K[0,:] += interpolation(dtvp * v_a)                                                # dt v' m*
    
    # kernel in v0
    #K[1,:-1] += (rho_a[1:]+rho_a[:-1])/2 * DF_C(rho_p, dz, BC_rhop, False)      # rho* grad(rho')
    #K[1,:-1] -= DF_C( rho_a * rho_p, dz, BC_rhoa_rhop, False)                       # div (rho* rho')
    K[1,:-1] -= (rho_p[1:]+rho_p[:-1])/2 * DF_C(rho_a, dz, BC_rhoa, False)
    K[1,:-1] += (p_a[1:]+p_a[:-1])/2 * DF_C(p_p, dz, BC_pp, False)              # p* grad (p')
    K[1,:-1] -= DF_C( gamma * p_a * p_p, dz, BC_gamma_pa_pp, False)                   # grad (gamma p* p')
    K[1,:-1] -= v_p*v_a * DF_C( U0.rho,dz, BC_rho0, False)
    K[1,:-1] -= U0.v*v_a * DF_C(rho_p, dz, BC_rhop, False)
    K[1,:-1] += 2 * DF_C(p_a * (gamma-1) * sigma_vp, dz, [a * (b-1) * c for a, b, c in zip(BC_pa, BC_gamma,BC_sigma_vp)], False)     # 2 div ( Sigma_v (v') (gamma-1) p* )  
        
    # kernel in p0
    K[2,:] += gamma * p_a * DF_C(v_p, dz, BC_vp)             # gamma p* div(v')
    K[2,:] -= p_a*DF_C(v_p, dz, BC_vp)                        # div (p* v')
    K[2,:] += rho_p / (Cv * (gamma-1) * U0.rho**2) * ddp     # rho'/(cv (gamma-1)rho0³) *  div( (1-gamma) kappa grad p*)
    
    return K
    

def get_kernels_decentered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dz, DF_DC):
    # intialisation
    K = np.zeros((3, len(rho_a)))
    
    BC_v0 = [U0.v[-2], U0.v[-1], U0.v[0], U0.v[1]]            # Points -2, -1, N+1, N+2 
    BC_vp = [v_p[-2], v_p[-1], v_p[0], v_p[1]]
    BC_va = [v_a[-2], v_a[-1], v_a[0], v_a[1]]
    BC_vavp = [a * b for a, b in zip(BC_va,BC_vp)]
    BC_vav0 = [a * b for a, b in zip(BC_va,BC_v0)]
    BC_pa = [p_a[-2], p_a[-1], p_a[0], p_a[1]]
    BC_rhoa = [rho_a[-2], rho_a[-1], rho_a[0], rho_a[1]]
    
    # kernel on rho0
    K[0,:] -= interpolation(v_p) * DF_DC( rho_a, dz, BC_rhoa, U0.v)                        #  v'div (rho*)

    # kernel on v0
    K[1,:-1] -= (U0.rho[1:]+U0.rho[:-1])/2 * DF_DC( v_p*v_a,dz, BC_vavp, U0.v, False)
    K[1,:-1] -= (rho_p[1:]+rho_p[:-1])/2 * DF_DC(U0.v*v_a, dz, BC_vav0, U0.v, False)
    K[1,:-1] += (U0.rho[1:]+U0.rho[:-1])/2 * DF_DC(v_p, dz, BC_vp, U0.v) * v_a    # rho0 grad v' m*
    K[1,:-1] += (rho_p[1:]+rho_p[:-1])/2 * DF_DC(U0.v, dz, BC_v0, U0.v) * v_a    # rho' grad v0 m*
    
    # kernel on p0
    K[2,:] -= interpolation(v_p) * DF_DC( p_a, dz, BC_pa, U0.v)    # div (p* v')
    
    return K


def get_kernels(hist_adjoint, hist_backprop, U0, T0, kappa, gamma, Cv, dt, dz, order_DF=4):
    # initialize
    K = np.zeros((3, len(hist_adjoint[0].rho)))
    if order_DF == 4:
        DF_DC = DF_DC_4
        DF_C = DF_C_4
    else :
        DF_DC = DF_DC_2
        DF_C = DF_C_2

    for t in range(len(hist_adjoint)-1):
        # create variable to be more lisible 
        rho_a = (hist_adjoint[t].rho + hist_adjoint[t+1].rho)/2    # mean at time n+1/2
        rho_p = (hist_backprop[len(hist_backprop) -  1 - t].rho + hist_backprop[len(hist_backprop) -  1 - t].rho)/2  # mean at time n+1/2
        
        p_a = (hist_adjoint[t].p + hist_adjoint[t+1].p)/2          # mean at time n+1/2
        p_p = (hist_backprop[len(hist_backprop) -  1 - t].p + hist_backprop[len(hist_backprop) -  t - 1].p)/2         # mean at time n+1/2 
        
        v_a = hist_adjoint[t].v                                   # at time n
        v_p = hist_backprop[len(hist_backprop) -  1 - t].v                                  # at time n
        
        dtvp = (hist_backprop[len(hist_backprop) - t - 2].v - v_p) / dt
        
        K += dt * get_kernels_centered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dtvp, dz, DF_C)
        K += get_kernels_decentered(rho_a, v_a, p_a, rho_p, v_p, p_p, U0, T0, kappa, gamma, Cv, dz, DF_DC)

    return K