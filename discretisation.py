#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 08:42:13 2022

@author: s.gerier
"""

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from copy import deepcopy
import pandas as pd
from scipy.interpolate import CubicSpline
from copy import deepcopy



#%% SPATIAL SCHEME

# numerical scheme

# ORDER 2
def DF_C_2(U, dz, BC, need_BC=True):
    if need_BC:
        BC = [BC[1],BC[-2]]
    else : 
        BC=  [ [], []]
    new_U = np.append(np.append([BC[0]], U), [BC[1]])
    diff_U = new_U[1:] - new_U[:-1]
    return diff_U / dz
   
def DF_DC_2(U, dz, BC, wind, is_reverse=False):
    diff_U = np.zeros(len(U))
    interp_wind = interpolation(wind)
    if is_reverse:
        interp_wind *= -1
    for i in range(len(U)):
        # depending on the sign of wind, the scheme is forward or backward
        if interp_wind[i] > 0 :
            if i < 1:
                diff_U[i] = (U[i] - BC[0]) /dz   
            else :
                diff_U[i] = (U[i] - U[i-1]) /dz  
        else:
            if i >= len(U) - 1:
                diff_U[i] = (BC[-1] - U[i]) /dz   
            else :
                diff_U[i] = (U[i+1] - U[i]) /dz  
    return diff_U

    
def DF_DC_2_backward(U0, dz, BC):
    U = np.append([BC[1]], U0)
    U = U[1:] - U[:-1]
    return U / dz

def DF_DC_2_forward(U0, dz, BC):
    U = np.append(U0, [BC[-2]])
    U = U[1:] - U[:-1]
    return U / dz

# ORDER 4
def DF_C_4(U, dz, BC, need_BC= True):
    if not need_BC: 
        BC =  BC[1:-1]
    new_U = np.append(np.append(BC[:int(len(BC)/2)], U), BC[-int(len(BC)/2):])
    Umm = new_U[:-3]
    Um = new_U[1:-2]
    Up = new_U[2:-1]
    Upp = new_U[3:]
    diff_U = Umm - 27 * Um + 27 * Up - Upp
    return diff_U / (24 * dz)

def DF_DC_4(U, dz, BC, wind, is_reverse=False):
    diff_U = np.zeros(len(U))
    interp_wind = interpolation(wind)
    print(is_reverse)
    if is_reverse:
        interp_wind *= -1
    for i in range(len(U)):
        # depending on the sign of wind, the scheme is forward or backward
        if interp_wind[i] > 0 :
            if i < 1:
                diff_U[i] = BC[0]/6 + BC[1]/2 - U[0] + U[1]/3 
            elif i == 1 :
                diff_U[i] = BC[1]/6 + U[0]/2 - U[1] + U[2]/3 
            elif i == len(U)-1: 
                diff_U[i] = U[i-2]/6 + U[i-1]/2 - U[i] + BC[-2]/3 
            else :
                diff_U[i] = U[i-2]/6 + U[i-1]/2 - U[i] + U[i+1]/3 
        else:
            if i >= len(U) - 1:
                diff_U[i] = -BC[-1]/3 + BC[-2] - U[i]/2 - U[i-1]/6  
            elif i == len(U) - 2 :
                diff_U[i] = -BC[-2]/3 + U[i+1] - U[i]/2 - U[i-1]/6
            elif i == 0 :
                diff_U[i] = -BC[1]/3 + U[i] - U[i+1]/2 - U[i+2]/6
            else :
                diff_U[i] = -U[i-1]/3 + U[i] - U[i+1]/2 - U[i+2]/6
                
    return diff_U / dz


def DF_DC_4_backward(U0, dz, BC):
    U = np.append(np.append(BC[:2], U0), [BC[-2]])
    Umm = U[:-3]
    Um = U[1:-2]
    Up = U[3:]
    U = Umm/6 + Um/2 - U[2:-1] + Up/3
    return U / dz

def DF_DC_4_forward(U0, dz, BC):
    U = np.append(np.append([BC[1]], U0), BC[-2:])
    Upp = U[3:]
    Up = U[2:-1]
    Um = U[:-3]
    U = - Upp/6 - Up/2 + U[1:-2] - Um/3
    return U / dz


def interpolation(U):
    # compute the value at point i+1/2 for i in [0,N-1]
    U = np.append(np.append([U[-1]], U), [U[0]])
    U = U[1:] + U[:-1]
    return U / 2


#%% BOUNDARY CONDITIONS

# sponge layer
def apply_sponge(U, index_z0 = 20, alpha = 13):
    # apply a sponge layer for z < z0 and z > zmax - z0
    for i in range(len(U.rho)):
        """
        if i < index_z0 +1:
            U.rho[i] *=  np.exp( - alpha * abs( i - index_z0) / index_z0)
            U.p[i] *=  np.exp( - alpha * abs( i - index_z0) / index_z0)
            # we have to consider that v is define on the center of each cell
            if (i < len(U.rho)-1) and (i < index_z0 - 1) :
                U.v[i] *=  np.exp( - alpha * abs( i - index_z0 + 1/2) / index_z0)
        elif i > len(U.rho) - index_z0 : 
            U.rho[i] *=  np.exp( - alpha * abs( i - len(U.rho) + index_z0) / index_z0)
            U.p[i] *=  np.exp( - alpha * abs( i - len(U.rho) + index_z0) / index_z0)            
            # we have to consider that v is define on the center of each cell
            if (i < len(U.rho)-1) and (i > len(U.rho) - index_z0 + 1) :
                U.v[i] *=  np.exp( - alpha * abs( i - len(U.rho) + index_z0 + 1/2) / index_z0)    
        """
        if i < index_z0 +1:
            sponge = 1 - (1- np.exp(alpha *  (( i - index_z0)/index_z0)**2))/(1-np.exp(alpha))
            U.rho[i] *= sponge
            if i < len(U.rho)-1:
                U.v[i] *= sponge
            U.p[i] *= sponge
        elif i > len(U.rho) - index_z0 :
            sponge = 1 - (1- np.exp(alpha *  (( i - len(U.rho) + index_z0 + 1)/index_z0)**2))/(1-np.exp(alpha))
            U.rho[i] *= sponge
            if i < len(U.v):
                U.v[i] *= 1 - (1- np.exp(alpha *  (( i - len(U.v) + index_z0 + 1)/index_z0)**2))/(1-np.exp(alpha)) 
            U.p[i] *= sponge          
    return U


test_sponge = False
if test_sponge : 
    from linearised_navier_stokes import *
    # plot the sponge 
    U = LNS_Variable(np.ones(500), np.ones(499), np.ones(500))
    sponge_U = apply_sponge(U)
    
    fig,ax = plt.subplots(3,1, figsize=(10,7))
    ax[0].vlines(20,0,1, 'k')
    ax[0].vlines(480,0,1,'k')
    ax[1].vlines(20,0,1, 'k')
    ax[1].vlines(480,0,1,'k')
    ax[2].vlines(20,0,1, 'k')
    ax[2].vlines(480,0,1,'k')
    ax[0].plot(sponge_U.rho)
    ax[1].plot(sponge_U.p)
    ax[2].plot(sponge_U.v)
    ax[0].set_xlabel("Altitude")
    ax[1].set_xlabel("Altitude")
    ax[2].set_xlabel("Altitude")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    
    
    
#%% TEMPORAL SCHEME

# Definition of temporal scheme - Runge Kutta of order 4
def RK4(f, U_t, T_init, Tmax, z, *args):
    fig,ax = plt.subplots(3,1, figsize=(10,7))
    ax[0].plot(z,U_t.rho, label=str(T_init))
    ax[1].plot(z,U_t.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    # Load the dt
    dt = args[9]
        
    # create a vector to save state at each time
    history = [deepcopy(U_t)]
    
    for t in np.arange(T_init,Tmax,dt):
        k1 = f(U_t, t, *args)
        aux = k1 * (dt / 2.)
        k2 = f( U_t + aux, t + dt / 2., *args)
        aux = k2 * (dt / 2.)
        k3 = f( U_t + aux, t + dt / 2., *args)
        aux = k3 * (dt)
        k4 = f( U_t + aux, t + dt, *args)
        
        U_t = U_t + (k1 + k2*2 + k3*2 + k4) * (dt/6)
        U_t = apply_sponge(U_t)
        
        # add state to history
        history += [deepcopy(U_t)]

        # make a figure to vizualise the transformation of the three quantities
        if (t+dt) % 5 < 1e-4 : 
            ax[0].plot(z,U_t.rho, label=str(t+dt))
            ax[1].plot(z,U_t.p)
            ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[0].set_xlabel("Altitude")
    ax[1].set_xlabel("Altitude")
    ax[2].set_xlabel("Altitude")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    fig.legend()
    return t, U_t, np.array(history)

# Definition of temporal scheme - Euler explicit
def EE(f, U_t, T_init, Tmax, z, *args):
    fig,ax = plt.subplots(3,1)
    ax[0].plot(z,U_t.rho, label="t = "+str(T_init))
    ax[1].plot(z,U_t.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    # load the dt
    dt = args[9]
    
    # create a vector to save state at each time
    history = [deepcopy(U_t)]
    for t in np.arange(T_init,Tmax,dt):
        U_t = U_t +  f(U_t, t, *args) * dt
        U_t = apply_sponge(U_t)
        if (t+dt) % 5 < 1e-4 : 
            ax[0].plot(z,U_t.rho, label='t = '+str(t+dt))
            ax[1].plot(z,U_t.p)
            ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)

        history += [deepcopy(U_t)]
        
    # make a plot to vizualise the transformation of the threee quantities    
    ax[0].set_xlabel("Altitude")
    ax[1].set_xlabel("Altitude")
    ax[2].set_xlabel("Altitude")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    fig.legend()
    return t, U_t, np.array(history)


# Definition of temporal scheme - Euler explicit with leap-frgo
def EE(f, U_t, T_init, Tmax, z, *args):
    # we suppose here at initialisation that v in given at time t=0, and 
    # rho and p are given at time t=1/2
    fig,ax = plt.subplots(3,1)
    ax[0].plot(z,U_t.rho, label="t = "+str(T_init))
    ax[1].plot(z,U_t.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    # load the dt
    dt = args[9]
    test_f = []
    # create a vector to save state at each time
    history = [deepcopy(U_t)]
    for t in np.arange(T_init,Tmax,dt):
        
        if not args[11] :# need to change the order of treatment in the backward propagation 
            tn = f(U_t, t, *args)
            U_t.v = U_t.v + tn.v * dt
            tndemi = f(U_t, t+dt/2, *args)
            U_t.rho = U_t.rho + tndemi.rho * dt
            U_t.p = U_t.p + tndemi.p * dt
            test_f += [tn.v]
        else :
            tn = f(U_t, t, *args)
            U_t.rho = U_t.rho + tn.rho * dt
            U_t.p = U_t.p + tn.p * dt
            tndemi = f(U_t, t+dt/2, *args)
            U_t.v = U_t.v + tndemi.v * dt
            test_f += [tndemi.v]
#        if True : #not args[11] : 
#            U_t =  U_t + f(U_t, t, *args) * dt 
#            test_f += [f(U_t, t, *args).rho] 
#        else : 
#            U_t =  U_t + f(args[12][len(args[12]) - 2 - int(round(t/dt))], t, *args) * dt   
#            test_f += [f(args[12][len(args[12]) - 2 - int(round(t/dt))], t, *args).rho]  
#        if not args[11] :# need to change the order of treatment in the backward propagation 
#            tn = f(U_t, t, *args)
#            U_t.v = U_t.v + tn.v * dt
#            tndemi = f(U_t, t+dt/2, *args)
#            U_t.rho = U_t.rho + tndemi.rho * dt
#            U_t.p = U_t.p + tndemi.p * dt
#            test_f += [tn.v]
#        else : 
#            aux = deepcopy(args[12][len(args[12]) - 2 - int(round(t/dt))])
#            aux.v = args[12][len(args[12]) - 1 - int(round(t/dt))].v
#            tn = f(aux, t, *args)
#            U_t.rho = U_t.rho + tn.rho * dt
#            U_t.p = U_t.p + tn.p * dt
#            
#            aux = deepcopy(U_t)
#            aux.v = args[12][len(args[12]) - 2 - int(round(t/dt))].v
#            tndemi = f(aux, t+dt/2, *args)
#            U_t.v = U_t.v + tndemi.v * dt
#            test_f += [tndemi.v]
            
#        #U_t = apply_sponge(U_t)
        
        if (t+dt) % 5 < 1e-4 : 
            ax[0].plot(z,U_t.rho, label='t = '+str(t+dt))
            ax[1].plot(z,U_t.p)
            ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)

        history += [deepcopy(U_t)]
        
    # make a plot to vizualise the transformation of the threee quantities    
    ax[0].set_xlabel("Altitude")
    ax[1].set_xlabel("Altitude")
    ax[2].set_xlabel("Altitude")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    fig.legend()
    return t, U_t, np.array(history), test_f
