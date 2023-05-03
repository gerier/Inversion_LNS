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
from Lib_NavierStokesEq.discretisation import *

# %% VARIABLES


class LNS_Variable:
    def __init__(self, rho, v, p):
        self.rho = rho
        self.v = v
        self.p = p

    def zeros_like(self, U):
        self.rho = np.zeros_like(U.rho)
        self.v = np.zeros_like(U.v)
        self.p = np.zeros_like(U.p)

    def __mul__(self, alpha):
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
        self.v = self.v / alpha
        self.p = self.p / alpha
        return self

    def __neg__(self):
        self.rho = -self.rho
        self.v = -self.v
        self.p = -self.p
        return self

    def plot(self, abs, title):
        fig, ax = plt.subplots(3, 1, figsize=(10, 7))
        ax[0].plot(abs, self.rho)
        ax[1].plot(abs, self.p)
        ax[2].plot((abs[1:] + abs[:-1])/2, self.v)
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
    def __init__(self, rho, v, p, gamma):
        self.rho = rho
        self.v = v
        self.p = p
        self.gamma = gamma

    def plot(self, abs, title):
        fig, ax = plt.subplots(4, 1, figsize=(10, 7))
        n = int(len(abs)/2)

        ax[0].plot(abs[:n], self.rho[:n])
        ax[1].plot(abs[:n], self.p[:n])
        ax[2].plot((abs[1:n]+abs[:n-1])/2, self.v[:n-1])
        ax[3].plot(abs[:n], np.sqrt(self.p[:n] * self.gamma / self.rho[:n]))
        ax[3].set_xlabel("Range (km)")
        ax[0].set_ylabel("Density")
        ax[1].set_ylabel("Pressure (Pa)")
        ax[2].set_ylabel("Wind (m/s)")
        ax[3].set_ylabel("Celerity (m/s)")
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        ax[3].grid()
        fig.suptitle(title)


# %% SOURCE

def F(z, t, f0, index_source):
    # initialisation
    f = np.zeros((3, len(z)))
    # define source term
    source = get_source_t_x(z, t, f0, index_source)
    # compute the contribution on quantity of mvt
    f[2, :] = source
    return f
    
def get_source(t, f0):
    t0 = 1.2/f0
    pi = 3.141592653589793238462643
    if False : 
        return np.exp(-4 * pi**2 * f0**2 * (t-t0)**2)
    else  :
        return  - 8 * pi**2 * f0**2 * (t-t0) * np.exp(-4 * pi**2 * f0**2 * (t-t0)**2)

def get_source_t_x(z, t, f0, index_source):
    # define source time function
    time_source = get_source(t, f0)
    # define the locality of the source
    spatial_source = z == z[index_source]
    # create source from temporal and spatial info
    source = time_source * spatial_source
    return source


# %% EQUATIONS

# Definition of function that describes the Linearised Navier-Stokes

def get_RHS(t, U1, compute_velocity, M, dz, dt, source, order_DF=4):
    # initialisation
    if order_DF == 4:
        DF_C = DF_C_4
        DF_DC = DF_DC_4_backward
    else:
        DF_C = DF_C_2
        DF_DC = DF_DC_2_backward

    # initialise
    RHS_c = np.zeros((3, len(U1.rho)))
    z, f0, index_source = source

    # boundary conditions
    BC_rho1 = [U1.rho[-2], U1.rho[-1], U1.rho[0], U1.rho[1]]
    BC_v1 = [U1.v[-2],   U1.v[-1],   U1.v[0],   U1.v[1]]
    BC_p1 = [U1.p[-2],   U1.p[-1],   U1.p[0],   U1.p[1]]
    BC_rho0 = [M.rho[-2],  M.rho[-1],  M.rho[0],  M.rho[1]]
    BC_v0 = [M.v[-2],    M.v[-1],    M.v[0],    M.v[1]]
    BC_p0 = [M.p[-2],    M.p[-1],    M.p[0],    M.p[1]]
    BC_vv0rho0 = [a*b*c for a,b,c in zip(BC_rho0, BC_v0,BC_v1)]


    if compute_velocity :
        # auxiliary calculations
        interp_rho0 = (M.rho[1:]+M.rho[:-1])/2
        interp_rho1 = (U1.rho[1:]+U1.rho[:-1])/2

        # compute the contribution on qunatity of mvt
        RHS_c[1, :-1] -= DF_DC(interp_rho0 * M.v * U1.v, dz, BC_vv0rho0, M.v)                 # div( rho0v0 x vp)
        #RHS_c[1, :-1] -= M.v * U1.v * DF_C(M.rho, dz, BC_rho0,False)                 # div( rho0v0 x vp)
        RHS_c[1, :-1] -= DF_C(U1.p, dz, BC_p1, False)                                         # grad p1
        RHS_c[1, :-1] -= (U1.v * interp_rho0 + interp_rho1 * M.v ) * DF_DC(M.v, dz, BC_v0, M.v) # (rhop v0 + rho0 vp) . nabla) v0
        RHS_c[1, :-1] = RHS_c[1, :-1] / interp_rho0

    else : 
        # auxiliary calculations
        interp_v1 = interpolation(U1.v)
        interp_v0 = interpolation(M.v)
        dv0 = DF_C(M.v, dz, BC_v0)
        dv1 = DF_C(U1.v, dz, BC_v1)

        # compute the contribution on density
        RHS_c[0, :] = -M.rho * dv1                                     # rho0 div vp
        RHS_c[0, :] -= U1.rho * dv0                                    # rhop div v0
        RHS_c[0, :] -= interp_v0 * DF_DC(U1.rho, dz, BC_rho1, M.v)     # v0 grad rhop
        RHS_c[0, :] -= interp_v1 * DF_DC(M.rho, dz, BC_rho0, M.v)      # vp grad rho0
        RHS_c[0, :] += F(z, t, f0, index_source)[2,:] * (M.rho / M.gamma / M.p) 

        # compute the contribution on pressure
        RHS_c[2, :] = - interp_v1 * DF_DC(M.p, dz, BC_p0, M.v)         # vp grad p0
        RHS_c[2, :] -= M.gamma * U1.p * dv0                                      # gamma pp div v0
        RHS_c[2, :] -= interp_v0 * DF_DC(U1.p, dz, BC_p1, M.v)                   # v0 grad pp
        RHS_c[2, :] -= M.gamma * M.p * dv1                                       # gamma p0 div vp$
        RHS_c[2, :] += F(z, t, f0, index_source)[2,:] 


    # Change the tye of RHS_c to LNS_Variable
    RHS_c = LNS_Variable(RHS_c[0, :],  RHS_c[1, :-1], RHS_c[2, :])

    return RHS_c * dt 



# %% ADJOINT EQUAITONS

def dchi(reverse_U, obs_U, t, dt, index_receivers, max_obs, n):
    diff = np.zeros((3, n))

    diff[2, index_receivers] = (obs_U - reverse_U) / max_obs

    return diff


def get_adjoint_RHS(t, Ua, compute_velocity, M, dz, dt, index_receivers, reverse_U, observation_U, max_obs, order_DF=4):
    if order_DF == 4:
        DF_C = DF_C_4
        DF_DC = DF_DC_4_backward
    else:
        DF_C = DF_C_2
        DF_DC = DF_DC_2_backward
    
    # initialisation
    RHS = np.zeros((3, len(Ua.rho)))

    # auxiliary calculations
    interp_rhoa = (Ua.rho[1:]+Ua.rho[:-1])/2
    interp_rho0 = (M.rho[1:]+M.rho[:-1])/2
    interp_pa = (Ua.p[1:]+Ua.p[:-1])/2
    rhoa_v0 = interp_rhoa * M.v
    pa_p0 = Ua.p * M.p

    # source
    Fadjoint = dchi(reverse_U, observation_U, t, dt, index_receivers, max_obs, len(Ua.rho))
    # boundary conditions
    BC_pa = [Ua.p[-2], Ua.p[-1],
             Ua.p[0], Ua.p[1]]
    BC_va = [Ua.v[-2], Ua.v[-1],
             Ua.v[0], Ua.v[1]]
    BC_v0 = [M.v[-2], M.v[-1], M.v[0], M.v[1]]
    BC_p0 = [M.p[-2], M.p[-1], M.p[0], M.p[1]]
    BC_rho0 = [M.rho[-2], M.rho[-1], M.rho[0], M.rho[1]]
    BC_rhoa = [Ua.rho[-2], Ua.rho[-1],
               Ua.rho[0], Ua.rho[1]]
    BC_rho0rhoa = [a*b for a,b in zip(BC_rho0, BC_rhoa)]
    BC_rhoav0 = [rhoa_v0[-2], rhoa_v0[-1], rhoa_v0[0], rhoa_v0[1]]
    BC_pap0 = [pa_p0[-2], pa_p0[-1],
                    pa_p0[0], pa_p0[1]]
    BC_pav0 = [a*b for a, b in zip(BC_pa, BC_v0)]

    #     
    dv0 = DF_C(M.v, dz, BC_v0)


    if compute_velocity:
        # contribution on va
        RHS[1, :-1] += DF_C(M.rho * Ua.rho, dz, BC_rho0rhoa, False)      # grad (rhoa rho0)
        RHS[1, :-1] -= interp_rhoa * DF_C(M.rho, dz, BC_rho0, False)     # rhoa grad rho0
        #RHS[1,:-1] += interp_rho0 * DF_C(Ua.rho, dz, BC_rho0, False)
        RHS[1, :-1] -= interp_pa * DF_C(M.p, dz, BC_p0, False)           # pa grad p0
        RHS[1, :-1] += M.gamma * DF_C(pa_p0, dz, BC_pap0, False)        # grad ( gamma pa p0)
        RHS[1, :-1] -= interp_rho0 * DF_DC(M.v, dz, BC_v0, M.v) * Ua.v
        RHS[1, :-1] += interp_rho0 * M.v * DF_DC(Ua.v, dz, BC_va, M.v) 
        #HS[1, :-1] -= Ua.v * M.v * DF_C(M.rho, dz, BC_rho0, False)      # rho0 (grad v0) va
        pass
    else : 
        # contribution on rhoa
        RHS[0, :] = - Ua.rho * dv0         # rhoa div v0
        RHS[0, :] += DF_C(rhoa_v0, dz, BC_rhoav0)      # div (rhoa v0)
        #RHS[0,:] += interpolation(M.v) * DF_DC(Ua.rho, dz, BC_rhoa, M.v)
        RHS[0, :] -= interpolation(M.v) * DF_C(M.v, dz, BC_v0) * interpolation(Ua.v)    # va . (v0 . nabla) v0

        # contribution on pa
        RHS[2, :] += DF_C(Ua.v, dz, BC_va)                      # div(va)
        RHS[2, :] -= M.gamma * Ua.p * DF_C(M.v, dz, BC_v0)      # gamma p* div v0
        RHS[2, :] += DF_C(interp_pa * M.v, dz, BC_pav0)         # div (pa v0)
        RHS[2, :] += Fadjoint[2, :]

    return LNS_Variable(RHS[0, :], RHS[1, :-1] / interp_rho0, RHS[2, :]) * dt


# %% KERNELS

def get_kernels_expr(rho_a, v_a, p_a, rho_p, v_p, p_p, M, dtvp, dz, source, DF_C, DF_DC):
    # initialization
    K = np.zeros((4, len(p_a)))

    # auxiliary calculations
    interp_pa = (p_a[1:]+p_a[:-1])/2
    gamma_pa_pp = M.gamma * p_a * p_p
    interp_pap0 = p_a * M.p
    interp_pap0 = (interp_pap0[1:]+interp_pap0[:-1])/2
    interp_va = interpolation(v_a)
    interp_vp = interpolation(v_p)
    interp_rho0 = (M.rho[1:]+M.rho[:-1])/2
    interp_rhop = (rho_p[1:]+rho_p[:-1])/2
    rho0vpva = interp_rho0 * v_a * v_p
    rhopv0va = interp_rhop * v_a * M.v

    t, z, f0, index_source = source
    f_source = F(z, t, f0, index_source)[2,:] 

    # main boundary conditions
    BC_vp = [v_p[-2], v_p[-1], v_p[0], v_p[1]]
    BC_v0 = [M.v[-2], M.v[-1], M.v[0], M.v[1]]
    BC_va = [v_a[-2], v_a[-1], v_a[0], v_a[1]]
    BC_rhoa = [rho_a[-2], rho_a[-1], rho_a[0], rho_a[1]]
    BC_rhop = [rho_p[-2], rho_p[-1], rho_p[0], rho_p[1]]
    BC_pa = [p_a[-2], p_a[-1], p_a[0], p_a[1]]
    BC_pp = [p_p[-2], p_p[-1], p_p[0], p_p[1]]
    BC_gamma_pa_pp = [gamma_pa_pp[-2], gamma_pa_pp[-1],
                      gamma_pa_pp[0], gamma_pa_pp[1]]
    BC_rhoavp = [a*b for a,b in zip(BC_rhoa, BC_vp)]
    BC_rhoprhoa = [a*b for a,b in zip(BC_rhoa, BC_rhop)]
    BC_v0vp = [a*b for a,b in zip(BC_vp, BC_v0)]
    BC_pavp = [a*b for a,b in zip(BC_vp, BC_pa)]
    BC_rho0vpva = [rho0vpva[-2], rho0vpva[-1], rho0vpva[0], rho0vpva[1]]
    BC_rhopv0va = [rhopv0va[-2], rhopv0va[-1], rhopv0va[0], rhopv0va[1]]

    # kernel in rho0
    K[0, :] -= DF_DC(interpolation(v_p) * rho_a, dz, BC_rhoavp, M.v)    # div(vp rhoa)
    K[0, :] += rho_a * DF_C(v_p, dz, BC_vp)                             # vp grad rhoa
    K[0, :] -= interpolation(dtvp * v_a)                                # va dt vp
    K[0, :] += interpolation(M.v) * interp_va * DF_C(v_p, dz, BC_va)    # v0 (vp . nabla) va
    K[0, :] += interp_vp  * interp_va * DF_C(M.v, dz, BC_v0)            # va (v' . nabla) v0 
    K[0, :] += f_source / (M.gamma * M.p)

    # kernel in v0
    K[1, :-1] = - DF_C(rho_p * rho_a, dz, BC_rhoprhoa, False)
    K[1, :-1] += (rho_a[1:]+rho_a[:-1])/2 * DF_C(rho_p, dz, BC_rhop, False)
    
    K[1, :-1] -= interp_rho0 * v_p * DF_DC(v_a, dz, BC_va, M.v)    # rho0 ma grad vp
    K[1, :-1] -= interp_rhop * v_a * DF_DC(M.v, dz, BC_v0, M.v)    # rho0 ma grad vp
    K[1, :-1] -= DF_DC(rho0vpva,dz, BC_rho0vpva, M.v)
    K[1, :-1] -= DF_DC(rhopv0va,dz, BC_rhopv0va, M.v)

    K[1, :-1] += interp_pa * DF_C(p_p, dz, BC_pp, False)
    K[1, :-1] -= DF_C(M.gamma * p_a * p_p, dz, BC_gamma_pa_pp,False)  # grad (gamma pa pp)
    
    
    # kernel in p0
    K[2, :] = p_a * M.gamma * DF_C(v_p, dz, BC_vp)  # pa gamma div vp
    K[2, :] -= DF_DC(p_a * interpolation(v_p), dz, BC_pavp)   # div(pa vp)
    K[2,:] -= f_source * M.rho / (M.gamma * M.p)

    # kernel in gamma
    #K[3, :] = p_a * p_p * DF_C(M.v, dz, BC_v0)     # pa pp div v0
    #K[3, :] += p_a * M.p * DF_C(v_p, dz, BC_vp)    # pa p0 div vp
    #K[3, :] -= DF_C( interp_pap0 * v_p, dz, [a*b for a,b in zip(BC_vp, BC_interp_pap0)])

    return K

"""
def get_kernels(hist_adjoint, hist_backprop, M, dt, dz, order_DF=4):
    # initialize
    K = np.zeros((4, len(hist_adjoint[0].rho)))
    if order_DF == 4:
        DF_C = DF_C_4
        DF_DC = DF_DC_4_backward
    else:
        DF_C = DF_C_2
        DF_DC = DF_DC_2_backward

    for t in range(1,len(hist_adjoint)):
        # create variable to be more lisible
        rho_a = hist_adjoint[t].rho     
        rho_p = hist_backprop[t].rho       

        p_a = hist_adjoint[t].p           
        p_p = hist_backprop[t].p         

        v_a = (hist_adjoint[t].v + hist_adjoint[t-1].v) /  2      # mean at time n+1/2                               # at time n
        v_p = (hist_backprop[t].v + hist_backprop[t-1].v) / 2     # mean at time n+1/2                              # at time n
        dtvp = -(hist_backprop[t-1].v - hist_backprop[t].v) / dt

        K += dt * get_kernels_expr(rho_a, v_a, p_a,
                                   rho_p, v_p, p_p, M, dtvp, dz, DF_C, DF_DC)


    return K
"""




def get_kernels(adjoint, backprop, adjoint_prec, backprop_prec, M, dt, dz, source, order_DF=4):
    # initialize
    K = np.zeros((4, len(adjoint.rho)))
    if order_DF == 4:
        DF_C = DF_C_4
        DF_DC = DF_DC_4_backward
    else:
        DF_C = DF_C_2
        DF_DC = DF_DC_2_backward


    # create variable to be more lisible
    rho_a = (adjoint.rho + adjoint_prec.rho)/2   # mean at time n+1/2   
    rho_p = (backprop.rho + backprop_prec.rho)/2   # mean at time n+1/2          

    p_a = (adjoint.p + adjoint_prec.p )/2   # mean at time n+1/2          
    p_p = (backprop.p+ backprop_prec.p)/2  # mean at time n+1/2      

    v_a = adjoint.v      # at time n
    v_p = backprop.v     # at time n
    dtvp = -(backprop_prec.v - backprop.v) / dt

    K += dt * get_kernels_expr(rho_a, v_a, p_a,
                                rho_p, v_p, p_p, M, dtvp, dz, source, DF_C, DF_DC)


    return K