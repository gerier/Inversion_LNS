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
import sys
# %% VARIABLES


class LNS_Variable:
    def __init__(self, v, p):
        self.v = v
        self.p = p

    def zeros_like(self, U):
        self.v = np.zeros_like(U.v)
        self.p = np.zeros_like(U.p)

    def __mul__(self, alpha):
        self.v *= alpha
        self.p *= alpha
        return self

    def __add__(self, other):
        self.v += other.v
        self.p += other.p
        return self

    def __sub__(self, other):
        self.v -= other.v
        self.p -= other.p
        return self

    def __truediv__(self, alpha):
        self.v = self.v / alpha
        self.p = self.p / alpha
        return self

    def __neg__(self):
        self.v = -self.v
        self.p = -self.p
        return self

    def plot(self, abs, title):
        fig, ax = plt.subplots(2, 1, figsize=(10, 7))
        ax[0].plot(abs, self.p)
        ax[1].plot((abs[1:] + abs[:-1])/2, self.v)
        ax[0].set_xlabel("Distance on axis x (km)")
        ax[1].set_xlabel("Distance on axis x (km)")
        ax[0].set_ylabel("Pressure")
        ax[1].set_ylabel("Velocity")
        ax[0].grid()
        ax[1].grid()
        fig.suptitle(title)


class LNS_Model:
    def __init__(self, rho, v, p, gamma, g):
        self.rho = rho
        self.v = v
        self.p = p
        self.gamma = gamma
        self.g = g

    def plot(self, abs, title):
        fig, ax = plt.subplots(4, 1, figsize=(10, 7))
        n = int(len(abs))

        ax[0].plot(abs[:n], self.rho[:n])
        ax[1].plot((abs[1:n]+abs[:n-1])/2, self.v[:n-1])
        ax[2].plot(abs[:n], self.p[:n])
        ax[3].plot(abs[:n], self.gamma[:n])
        ax[3].set_xlabel("Range (km)")
        ax[0].set_ylabel("Density")
        ax[1].set_ylabel("Wind")
        ax[2].set_ylabel("Pressure")
        ax[3].set_ylabel(r"$\gamma$")
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        ax[3].grid()
        fig.suptitle(title)


# %% EQUATIONS

# Definition of function that describes the Linearised Navier-Stokes

if True:

    def F(z, t, f0, index_source):
        # initialisation
        f = np.zeros((2, len(z)))
        # define source term
        source = get_source_t_x(z, t, f0, index_source)
        # compute the contribution on quantity of mvt
        f[1, :] = source
        return f
    # Definition of resolution spatial

    def get_RHS(t, previous_U, M, dz, dt, source, order_DF=2):
        # initialisation
        if order_DF == 4:
            DF_C = DF_C_4
            DF_DC = DF_DC_4_backward
        else:
            DF_C = DF_C_2
            DF_DC = DF_DC_2#_backward

        # initialise
        U1 = previous_U
        RHS_c = np.zeros((2, len(U1.p)))
        z, f0, index_source = source

        # boundary conditions
        BC_v1 = [U1.v[-2],   U1.v[-1],   U1.v[0],   U1.v[1]]
        BC_p1 = [U1.p[-2],   U1.p[-1],   U1.p[0],   U1.p[1]]
        BC_p0 = [M.p[-2],    M.p[-1],    M.p[0],    M.p[1]]

        # auxiliary calculations
        interp_rho0 = (M.rho[1:]+M.rho[:-1])/2
        interp_v1 = interpolation(U1.v)
        interp_v0 = interpolation(M.v)
        dv1 = DF_C(U1.v, dz, BC_v1)

        # compute the contribution on qunatity of mvt
        RHS_c[0, :-1] -= DF_C(U1.p, dz, BC_p1, False) / interp_rho0        # grad p1

        # compute the contribution on pressure
        RHS_c[1, :] -= M.gamma * M.p * dv1                                # gamma p0 div vp$
        RHS_c[1, :] += F(z, t, f0, index_source)[1,:]                     # pressure force

        # Change the tye of RHS_c to LNS_Variable
        RHS_c = LNS_Variable(RHS_c[0, :-1], RHS_c[1, :])

        return RHS_c * dt 


# %% SOURCE

def get_source(t, f0):
    t0 = 1.2/f0
    pi = 3.141592653589793238462643
    if True : 
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


# %% ADJOINT EQUAITONS

def dchi(reverse_U, obs_U, t, dt, index_receivers, max_obs):
    pass
    return None


def get_adjoint_RHS(Ustar, t, previous_Ustar, M, dz, dt, reverse_U, observation_U, index_receivers, max_obs, order_DF=4):
    pass
    return None


# %% KERNELS

def get_kernels_expr(rho_a, v_a, p_a, rho_p, v_p, p_p, M, dtvp, dz, DF_C, DF_DC):
    pass
    return None


def get_kernels(hist_adjoint, hist_backprop, M, dt, dz, order_DF=4):
    pass
    return None
