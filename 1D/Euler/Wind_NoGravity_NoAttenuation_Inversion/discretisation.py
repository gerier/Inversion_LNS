#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 08:42:13 2022

@author: s.gerier
"""

import numpy as np
import matplotlib.pyplot as plt
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
   
def DF_DC_2(U, dz, BC, wind):
    diff_U = np.zeros(len(U))
    interp_wind = interpolation(wind)
    #if is_reverse:
    #    interp_wind *= -1
    for i in range(len(U)):
        # depending on the sign of wind, the scheme is forward or backward
        if interp_wind[i] > 0 :
            if i < 1:
                diff_U[i] = (U[i] - BC[1]) /dz   
            else :
                diff_U[i] = (U[i] - U[i-1]) /dz  
        else:
            if i >= len(U) - 1:
                diff_U[i] = (BC[-2] - U[i]) /dz   
            else :
                diff_U[i] = (U[i+1] - U[i]) /dz  
    return diff_U

    
def DF_DC_2_backward(U0, dz, BC, wind=0):
    U = np.append([BC[1]], U0)
    U = U[1:] - U[:-1]
    return U / dz

def DF_DC_2_forward(U0, dz, BC, wind=0):
    U = np.append(U0, [BC[-2]])
    U = U[1:] - U[:-1]
    return U / dz

# ORDER 4
# def DF_C_4(U, dz, BC, need_BC= True):
#     if not need_BC: 
#         BC =  BC[1:-1]
#     new_U = np.append(np.append(BC[:int(len(BC)/2)], U), BC[-int(len(BC)/2):])
#     Umm = new_U[:-3]
#     Um = new_U[1:-2]
#     Up = new_U[2:-1]
#     Upp = new_U[3:]
#     diff_U = Umm - 27 * Um + 27 * Up - Upp
#     return diff_U / (24 * dz)
def DF_C_4(U, dz, BC, need_BC= True):
    if not need_BC: 
        BC =  BC[1:-1]
    demi_n = int(len(BC)/2) 
    aux_U = np.empty(len(U)+2*demi_n) 
    aux_U[demi_n:-demi_n] = U
    aux_U[:demi_n] = BC[:demi_n]
    aux_U[-demi_n:] = BC[demi_n:]
    DU = aux_U[:-3] - 27 * aux_U[1:-2] + 27 * aux_U[2:-1] - aux_U[3:]
    return DU / (24 * dz)


def DF_DC_4(U, dz, BC, wind):
    diff_U = np.zeros(len(U))
    interp_wind = interpolation(wind)
    
    for i in range(len(U)):
        # depending on the sign of wind, the scheme is forward or backward
        if interp_wind[i] > 0 :
            if i < 1:
                diff_U[i] = BC[0]/6 - BC[1] + U[0]/2 + U[1]/3 
            elif i == 1 :
                diff_U[i] = BC[1]/6 - U[0] + U[1]/2 + U[2]/3 
            elif i == len(U)-1: 
                diff_U[i] = U[i-2]/6 - U[i-1] + U[i]/2 + BC[-2]/3 
            else :
                diff_U[i] = U[i-2]/6 - U[i-1] + U[i]/2 + U[i+1]/3 
        else:
            if i >= len(U) - 1:
                diff_U[i] = -BC[-1]/6 + U[i+1] - U[i]/2 - U[i-1]/3  
            elif i == len(U) - 2 :
                diff_U[i] = -BC[-2]/6 + BC[-1] - U[i]/2 - U[i-1]/3
            elif i == 0 :
                diff_U[i] = -U[i+2]/6 + U[i+1] - U[i]/2 - BC[1]/3
            else :
                diff_U[i] = -U[i+2]/6 + U[i+1] - U[i]/2 - U[i-1]/3
                
    return diff_U / dz


# def DF_DC_4_backward(U0, dz, BC, wind=0):
#     U = np.append(np.append(BC[:2], U0), [BC[-2]])
#     Umm = U[:-3]
#     Um = U[1:-2]
#     Up = U[3:]
#     U = Umm/6 + Um/2 - U[2:-1] + Up/3
#     return U / dz
def DF_DC_4_backward(U0, dz, BC, wind=0):
    aux_U = np.empty(len(U0)+3) 
    aux_U[2:-1] = U0
    aux_U[:2] = BC[:2]
    aux_U[-1] = BC[-2]
    DU = aux_U[:-3]/6 - aux_U[1:-2] + U0/2 + aux_U[3:]/3
    return DU / dz 

# def DF_DC_4_forward(U0, dz, BC, wind=0):
#     U = np.append(np.append([BC[1]], U0), BC[-2:])
#     Upp = U[3:]
#     Up = U[2:-1]
#     Um = U[:-3]
#     U = - Upp/6 - Up/2 + U[1:-2] - Um/3
#     return U / dz
def DF_DC_4_forward(U0, dz, BC, wind=0):
    aux_U = np.empty(len(U0)+3) 
    aux_U[1:-2] = U0
    aux_U[0] = BC[1]
    aux_U[-2:] = BC[-2:]
    DU = - aux_U[3:]/6 + aux_U[2:-1] - U0/2 - aux_U[:-3]/3
    # periodic
    # DU = - np.roll(U0, -2)/6 - np.roll(U0,-1)/2 + U0 - np.roll(U0,1)/3
    return DU / dz

# def interpolation(U):
#     # compute the value at point i+1/2 for i in [0,N-1]
#     U = np.append(np.append([U[-1]], U), [U[0]])
#     U = U[1:] + U[:-1]
#     return U / 2
def interpolation(U):
    # compute the value at point i+1/2 for i in [0,N-1]
    I = np.empty(len(U)+1)
    I[1:-1] = U[1:] + U[:-1]
    I[0] = U[0] + U[-1]
    I[-1] = I[0] 
    return I / 2
    
    
#%% TEMPORAL SCHEME

# Definition of temporal scheme - Euler explicit with leap-frgo
def EE(f_dt, U_t, T_init, Tmax, z, propagation_type, *args):
    # we suppose here at initialisation that v in given at time t=0, and 
    # rho and p are given at time t=1/2
  
    # load the dt
    dt = args[2]
    if Tmax >= T_init : 
        t_tab = np.arange(T_init,Tmax,dt)
    else : 
        t_tab = np.arange(T_init,Tmax,-dt)

    # create a vector to save state at each time
    history = [deepcopy(U_t)]  
    if propagation_type == "adjoint" : 
        #t_tab = t_tab = np.append(t_tab, t_tab[-1]-dt)
        pass#propagation_type = "forward"
     
    for it, t in enumerate(t_tab):

        if propagation_type == "forward" :# need to change the order of treatment in the backward propagation 
            aux = deepcopy(U_t)
            tndemi = f_dt(t, aux, False, *args)
            U_t.p +=   tndemi.p
            U_t.rho += tndemi.rho
            aux = deepcopy(U_t)
            tn= f_dt(t+dt/2, aux, True, *args) 
            U_t.v +=  tn.v
        elif propagation_type == "adjoint" :
            aux = deepcopy(U_t)
            tn = f_dt(t, aux, True, *args)
            U_t.v += tn.v           
            aux = deepcopy(U_t)
            tndemi = f_dt(t, aux, False, *args)
            U_t.p += tndemi.p
            U_t.rho += tndemi.rho
           
        
        history += [deepcopy(U_t)]
        
    #     if (it % 2000) < 0.0004 :
    #        plt.plot(history[-1].v, label="it"+str(it))
    # plt.legend()
    # plt.show()
    return t, U_t, history

