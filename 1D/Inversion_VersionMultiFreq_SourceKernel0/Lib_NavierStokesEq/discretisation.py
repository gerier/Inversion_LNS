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
def EE(f_dt, U_init, it_start, it_end, dt, index_receivers, propagation_type, sponge_layer, *args, compute_h0=False):
    # we suppose here at initialisation that v in given at time t=0, and 
    # rho and p are given at time t=1/2
  
    U_t = deepcopy(U_init)
    # create a vector to save state at each time
    history = [deepcopy(U_t.p[index_receivers])]  
    if propagation_type == "adjoint" : 
            args = args[:2] + (dt, index_receivers) + args[2:]
    else : 
        args = args[:-1] + (dt,) + args[-1:]
    
    h0 = 0
    h0_demi = 0
    U_old = deepcopy(U_init)

    sponge, interp_sponge = sponge_layer

    for it in range(int(it_start), int(it_end)):
        t = it * dt 
        
        if propagation_type == "forward" : 
            aux = deepcopy(U_t)
            tndemi = f_dt(t, aux, False, *args)
            U_t.p +=   tndemi.p
            U_t.rho += tndemi.rho
            U_t.p *= sponge
            U_t.rho *= sponge

            aux = deepcopy(U_t)
            tn= f_dt(t+dt/2, aux, True, *args) 
            U_t.v +=  tn.v
            U_t.v *= interp_sponge 

        elif propagation_type == "adjoint" :
            aux = deepcopy(U_t)
            tn = f_dt(t, aux, True, *args)
            U_t.v += tn.v           
            U_t.v *= interp_sponge

            aux = deepcopy(U_t)
            tndemi = f_dt(t, aux, False, *args)
            U_t.p += tndemi.p
            U_t.rho += tndemi.rho
            U_t.p *= sponge
            U_t.rho *= sponge
              
        if compute_h0 : 
            h0_demi += ((U_t.v - U_old.v) / dt)**2 * dt
            h0 += ((interpolation(U_t.v) - interpolation(U_old.v)) / dt)**2 * dt
            h0[int(1000*15/50)] = 1
            if it == int(it_end)-1 :
                return history, U_t, h0_demi, h0

        history += [deepcopy(U_t.p[index_receivers])]
        if (it*dt) % 1 < 1e-4 and False:
            plt.plot(U_t.v)
            plt.show()

    if False : #it_start == 0 and it_end >= 70/args[-3] : 
        M = args[0]
        z = np.linspace(0,(len(U_t.rho)-1)*args[-3], len(U_t.rho))/1000
        fig, ax = plt.subplots(4,1, figsize=(12,8))
        ax[0].plot(z,M.rho)
        ax[1].plot(z,M.p)
        ax[2].plot(z[:-1],M.v)
        ax[3].plot(z,U_t.p)
        ax[0].set_title("Fin de problème direct")
        plt.show()

    return history, U_t


# %%
def apply_sponge(U, index_z0 = 200, alpha = 13):
    Usponge = deepcopy(U)
    # apply a sponge layer for z < z0 and z > zmax - z0
    for i in range(len(U.rho)):

        if i < 4  or  i > len(U.rho) - 1- 4: 
            sponge = 0
            Usponge.rho[i] *= sponge
            Usponge.p[i] *= sponge

            if i != len(U.rho) - 1:
                Usponge.v[i] *= sponge

        elif i < 4+index_z0 +1 :
            sponge = 1 - (1- np.exp(alpha *  (( i - 4 - index_z0)/index_z0)**2))/(1-np.exp(alpha))
            next_sponge = 1 - (1- np.exp(alpha *  (( i + 1 - 4 - index_z0)/index_z0)**2))/(1-np.exp(alpha))
            sponge_demi = 1/2 * (sponge + next_sponge)
            
            Usponge.rho[i] *= sponge
            Usponge.v[i] *= sponge_demi
            Usponge.p[i] *= sponge
        elif i > len(U.rho) - 1-4-index_z0 : 
            sponge = 1 - (1- np.exp(alpha *  (( (len(U.rho) -4- i - 1) - index_z0)/index_z0)**2))/(1-np.exp(alpha))
            next_sponge = 1 - (1- np.exp(alpha *  (( (len(U.rho) -4 - i  - 2) - index_z0)/index_z0)**2))/(1-np.exp(alpha))
            sponge_demi = 1/2 * (sponge + next_sponge)
                
            Usponge.rho[i] *= sponge
            if i != len(U.rho) - 1 :
                Usponge.v[i-1] *= sponge_demi
            Usponge.p[i] *= sponge   
 
    return Usponge



def get_sponge(size, index_z0 = 200, alpha=13):
    sponge = np.ones(size)
    interp_sponge = np.ones(size-1)
    
    # the first four element are set to 0  
    sponge[:4] = 0 
    sponge[size-4:] = 0
    interp_sponge[:4] = 0 
    interp_sponge[size-4:] = 0
   
    for i in range(4, 4+index_z0 + 1):
        sponge[i] = 1 - (1- np.exp(alpha *  (( i - 4 - index_z0)/index_z0)**2))/(1-np.exp(alpha))
        next_sponge = 1 - (1- np.exp(alpha *  (( i + 1 - 4 - index_z0)/index_z0)**2))/(1-np.exp(alpha))
        interp_sponge[i] = 1/2 * (sponge[i] + next_sponge)
            
    for i in range(size-4-index_z0 , size-4):
        sponge[i] = 1 - (1- np.exp(alpha *  (( (size -4- i - 1) - index_z0)/index_z0)**2))/(1-np.exp(alpha))
        next_sponge = 1 - (1- np.exp(alpha *  (( (size -4 - i  - 2) - index_z0)/index_z0)**2))/(1-np.exp(alpha))
        interp_sponge[i-1] = 1/2 * (sponge[i] + next_sponge)
 
    return sponge, interp_sponge
# %%

# define the laplacian operator that is used for regularisation of the misfit function

def get_laplacian(u, bc, dz, order=2):
    laplacian_u = np.zeros_like(u)

    if order == 2 : 
        for i in range(1, len(u)-1):
            laplacian_u[i] = u[i-1] - 2 * u[i] +  u[i+1]  

        laplacian_u[0] = bc[1] - 2 * u[0] +  u[1]
        laplacian_u[-1] = u[-2] - 2 * u[-1] +  bc[-2]

        laplacian_u = laplacian_u / dz**2  

    elif order == 4 : 
        for i in range(2, len(u)-2):
            laplacian_u[i] = - u[i-2] + 16 * u[i-1] - 30 * u[i] + 16* u[i+1] - u[i+2] 

        laplacian_u[0]  = - bc[0] + 16 * bc[1] - 30 * u[i] + 16* u[i+1]    - u[i+2] 
        laplacian_u[1]  = - bc[1] + 16 * u[0]     - 30 * u[i] + 16* u[i+1]    - u[i+2] 
        laplacian_u[-2] = - u[i-2]   + 16 * u[i-1]   - 30 * u[i] + 16* u[i+1]    - bc[-2]  
        laplacian_u[-1] = - u[i-2]   + 16 * u[i-1]   - 30 * u[i] + 16* bc[-2] - bc[-1] 

        laplacian_u = laplacian_u / 12 / dz**2

    return laplacian_u 


# def get_du_dz(u, bc, dz, order=2):
#     du_dz = np.zeros_like(u)

#     if order == 2 : 
#         for i in range(1, len(u)-1):
#             du_dz = - u[i-1] +  u[i+1]  

#         du_dz[0] = - u[bc[1]]  +  u[i+1]
#         du_dz[-1] = -u[i-1] +  u[bc[-2]]

#         du_dz = du_dz / dz / 2 

#     elif order == 4 : 
#         for i in range(2, len(u)-2):
#             du_dz = u[i-2] - 8 * u[i-1] + 8 * u[i+1] - u[i+2] 

#         du_dz[0]  = u[bc[0]] - 8 * u[bc[1]] + 8 * u[i+1]    - u[i+2] 
#         du_dz[1]  = u[bc[1]] - 8 * u[i-1]   + 8 * u[i+1]    - u[i+2] 
#         du_dz[-2] = u[i-2]   - 8 * u[i-1]   + 8 * u[i+1]    - u[bc[-2]]  
#         du_dz[-1] = u[i-2]   - 8 * u[i-1]   + 8 * u[bc[-2]] - u[bc[-1]] 

#         du_dz = du_dz / 12 / dz

#     return du_dz  

def apply_laplacian_regularisation(m, dz, size):
    param_m = m[2*size:] # apply the regularisation on the velocity perturbation (only)
    laplacian_m = get_laplacian(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    norm2 = sum( laplacian_m**2)
    return dz * norm2 / 2

# def get_derivative_laplacian_regularisation(m, dz, size):
#     param_m = m[2*size:] # apply the regularisation on the velocity perturbation (only)
#     laplacian_m = get_laplacian(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
#     dmdz = get_du_dz(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)

#     dlaplacian = laplacian_m / dmdz

#     return dlaplacian
def get_derivative_laplacian_regularisation(m, dz, size):
    param_m = m[2*size:] # apply the regularisation on the velocity perturbation (only)
    laplacian_m = get_laplacian(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    dlaplacian = get_laplacian(laplacian_m, [laplacian_m[-2], laplacian_m[-1], laplacian_m[0], laplacian_m[1]], dz)
    return dlaplacian


def get_laplacian_matrix(size,order):

    if order == 2:
        laplacian = np.diag( - 2* np.ones(size))
        laplacian += np.diag( np.ones(size-1), 1)
        laplacian += np.diag( np.ones(size-1), -1)

    elif order == 4:
        laplacian += np.diag( 8*np.ones(size-1), 1)
        laplacian += np.diag( -8*np.ones(size-1), -1)
        laplacian += np.diag( -np.ones(size-2), 2)
        laplacian += np.diag( np.ones(size-2), -2)

    return laplacian










def get_gradient(u, bc, dz, order=2):
    laplacian_u = np.zeros_like(u)

    if order == 2 : 
        for i in range(1, len(u)-1):
            laplacian_u[i] = - u[i-1] + u[i+1]  

        laplacian_u[0] = - bc[1] +  u[1]
        laplacian_u[-1] = - u[-2] +  bc[-2]

        laplacian_u = laplacian_u / 2 / dz  

    elif order == 4 : 
        for i in range(2, len(u)-2):
            laplacian_u[i] = - u[i-2] - 8 * u[i-1] + 8 * u[i+1] - u[i+2] 

        laplacian_u[0]  = bc[0] - 8 * bc[1] + 8 * u[i+1]    - u[i+2] 
        laplacian_u[1]  = bc[1] - 8 * u[0]  + 8 * u[i+1]    - u[i+2] 
        laplacian_u[-2] = u[i-2] - 8 * u[i-1]  + 8 * u[i+1]    - bc[-2]  
        laplacian_u[-1] = u[i-2] - 8 * u[i-1]  + 8 * bc[-2] - bc[-1] 

        laplacian_u = laplacian_u / 12 / dz

    return laplacian_u 

def apply_gradient_regularisation(m, dz, size):
    dreg = 0

    param_m = m[2*size:] # apply the regularisation on the velocity perturbation (only)
    gradient_m = get_gradient(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    norm2 = sum(gradient_m**2)
    dreg += dz * norm2 / 2

    param_m = m[size:2*size] # apply the regularisation on the velocity perturbation (only)
    gradient_m = get_gradient(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    norm2 = sum(gradient_m**2)
    dreg+= dz * norm2 / 2
    
    param_m = m[:size] # apply the regularisation on the velocity perturbation (only)
    gradient_m = get_gradient(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    norm2 = sum(gradient_m**2)
    dreg+= dz * norm2 / 2

    return dreg


def get_derivative_gradient_regularisation(m, dz, size):
    dlaplacian = np.zeros_like(m)

    param_m = m[2*size:] # apply the regularisation on the velocity perturbation (only)
    laplacian_m = get_gradient(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    dlaplacian[2*size:] = get_gradient(laplacian_m, [laplacian_m[-2], laplacian_m[-1], laplacian_m[0], laplacian_m[1]], dz)
    
    param_m = m[size:2*size] # apply the regularisation on the velocity perturbation (only)
    laplacian_m = get_gradient(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    dlaplacian[size:2*size] = get_gradient(laplacian_m, [laplacian_m[-2], laplacian_m[-1], laplacian_m[0], laplacian_m[1]], dz)
    
    param_m = m[:size] # apply the regularisation on the velocity perturbation (only)
    laplacian_m = get_gradient(param_m, [param_m[-2], param_m[-1], param_m[0], param_m[1]], dz)
    dlaplacian[:size] = get_gradient(laplacian_m, [laplacian_m[-2], laplacian_m[-1], laplacian_m[0], laplacian_m[1]], dz)
    
    return dlaplacian

