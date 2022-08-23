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


def DF_DC_4_backward(U0, dz, BC, wind=0, is_reverse=0):
    U = np.append(np.append(BC[:2], U0), [BC[-2]])
    Umm = U[:-3]
    Um = U[1:-2]
    Up = U[3:]
    U = Umm/6 + Um/2 - U[2:-1] + Up/3
    return U / dz

def DF_DC_4_forward(U0, dz, BC, wind=0, is_reverse=0):
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
def RK4(f, U_t, T_init, Tmax, z, propagation_type, *args):
    fig,ax = plt.subplots(3,1, figsize=(10,7))
    ax[0].plot(z,U_t.rho, label=str(T_init))
    ax[1].plot(z,U_t.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    # Load the dt
    dt = args[9]
    t_tab = np.arange(T_init,Tmax,dt)        

    # create a vector to save state at each time
    history = [deepcopy(U_t)]
    if propagation_type == "forward" : 
    
        for t in np.arange(T_init,Tmax,dt):
            k1 = f(U_t, t, U_t, *args)
            aux = k1 * (dt / 2.)
            k2 = f( U_t + aux, t + dt / 2., U_t, *args)
            aux = k2 * (dt / 2.)
            k3 = f( U_t + aux, t + dt / 2., U_t, *args)
            aux = k3 * (dt)
            k4 = f( U_t + aux, t + dt, U_t, *args)
            
            U_t = (k1 + k2*2 + k3*2 + k4) /6
            #U_t = apply_sponge(U_t)
            
            # add state to history
            history += [deepcopy(U_t)]

            # make a figure to vizualise the transformation of the three quantities
            if (t+dt) % 5 < 1e-4 : 
                ax[0].plot(z,U_t.rho, label=str(t+dt))
                ax[1].plot(z,U_t.p)
                ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    elif  propagation_type == "checkpointing" :
        index_fr0 = np.arange(0,len(t_tab),1)[:-1:100]
        step_btw_fr0 = index_fr0[1] - index_fr0[0]
        Nfr_1 = 50
        q,r = divmod(step_btw_fr0, Nfr_1+1)
        index_fr1 = np.array([i * q + min(i,r) for i in range(Nfr_1+1)])
        print(index_fr1, index_fr0)
        checkpoint1 = []
        checkpoint0 = deepcopy(args[12][index_fr0])
        
        t_tab = np.append(t_tab, t_tab[-1]+dt)
        history = []
        
        for it, t in enumerate(t_tab):
            # get the nb of the previous frame and the nb of iteration between the previous frame and the actual time
            #which_fr, which_it_infr = divmod(len(t_tab) - 1 - it,step_btw_fr0)
            aux = np.where(len(t_tab) - 1 - it >= index_fr0)
            which_fr = len(aux[0])-1
            which_it_infr = len(t_tab) - 1 - it - index_fr0[aux[0][-1]]
            # get the number of the previous frame of the second chekpointing and the nb of iter
            aux = np.where(which_it_infr >= index_fr1)[0]
            which_fr1 = len(np.where(which_it_infr > index_fr1)[0])-1
            which_it_infr1 = which_it_infr - index_fr1[aux[-1]]
            
            # option 1 : first iteration, need to create the checkpointing1
            if it==0 :
                aux_U_t = deepcopy(checkpoint0[-1])
                checkpoint1 = []
                for it_fr1 in range(index_fr0[-1], len(t_tab)-1):
                    # iteration to record the checkpoint
                    k1 = f(aux_U_t, t_tab[it_fr1], aux_U_t, *args[:-1])
                    aux = k1 * (dt / 2.)
                    k2 = f( aux_U_t + aux, t_tab[it_fr1] + dt / 2., aux_U_t, *args[:-1])
                    aux = k2 * (dt / 2.)
                    k3 = f( aux_U_t + aux, t_tab[it_fr1] + dt / 2., aux_U_t, *args[:-1])
                    aux = k3 * (dt)
                    k4 = f( aux_U_t + aux, t_tab[it_fr1] + dt, aux_U_t, *args[:-1])                  
                    U_t = (k1 + k2*2 + k3*2 + k4) /6
                    # save if it is at a checkpoint time
                    # note that "+1" because we get U^n+1
                    if (it_fr1 - index_fr0[-1] + 1) in index_fr1:
                        checkpoint1 += [deepcopy(aux_U_t)]
                # load U^N+1
                U_t = deepcopy(aux_U_t)
                  
            # option 2 : last iteration, already saved in checkpoint0
            elif  it == len(t_tab) - 1:
                # load U^0
                U_t = deepcopy(checkpoint0[which_fr])  
                              
            # option 3 : actual time has a value in checkpoint0
            elif which_it_infr == 0 :
                # load the value of time t : U^n
                U_t = checkpoint0[which_fr]
                #  need to load new checkpoint1
                aux_U_t = deepcopy(checkpoint0[which_fr-1])
                checkpoint1 = []
                for it_fr1 in range(index_fr0[which_fr-1], index_fr0[which_fr]-1):
                    k1 = f(aux_U_t, t_tab[it_fr1], aux_U_t, *args[:-1])
                    aux = k1 * (dt / 2.)
                    k2 = f( aux_U_t + aux, t_tab[it_fr1] + dt / 2., aux_U_t, *args[:-1])
                    aux = k2 * (dt / 2.)
                    k3 = f( aux_U_t + aux, t_tab[it_fr1] + dt / 2., aux_U_t, *args[:-1])
                    aux = k3 * (dt)
                    k4 = f( aux_U_t + aux, t_tab[it_fr1] + dt, aux_U_t, *args[:-1])
                    U_t = (k1 + k2*2 + k3*2 + k4) /6
                    # save if it is at a checkpoint time 
                    if (it_fr1 - index_fr0[which_fr-1] + 1) in index_fr1:
                        # note that "+1" because we get U^n+1
                        checkpoint1 += [deepcopy(aux_U_t)]
   
            # option 4 :  actual time has a value in checkpoint1
            elif which_it_infr1 == 0:
                U_t = deepcopy(checkpoint1[which_fr1])
                
            # option 5 : no checkpoint exist
            else :
                # load the more recent checkpoint
                if which_fr1 == 0 :
                    aux_U_t = deepcopy(checkpoint0[which_fr])
                else : 
                    aux_U_t = deepcopy(checkpoint1[which_fr1-1])   # -1 because we neglect the checkpoint in 0
                
                for it_fr1 in range(index_fr0[which_fr] + index_fr1[which_fr1], len(t_tab)-it-1):
                    # iterate until finding the actual time
                    k1 = f(aux_U_t, t_tab[it_fr1], aux_U_t, *args[:-1])
                    aux = k1 * (dt / 2.)
                    k2 = f( aux_U_t + aux, t_tab[it_fr1] + dt / 2., aux_U_t, *args[:-1])
                    aux = k2 * (dt / 2.)
                    k3 = f( aux_U_t + aux, t_tab[it_fr1] + dt / 2., aux_U_t, *args[:-1])
                    aux = k3 * (dt)
                    k4 = f( aux_U_t + aux, t_tab[it_fr1] + dt, aux_U_t, *args[:-1])
                    U_t = (k1 + k2*2 + k3*2 + k4) /6
                # save if it is at a checkpoint time
                U_t = deepcopy(aux_U_t)

        
            # add state to history
            history += [deepcopy(U_t)]
            
            # make a figure to vizualise the transformation of the three quantities
            if (t+dt) % 5 < 1e-4 : 
                ax[0].plot(z,U_t.rho, label=str(t+dt))
                ax[1].plot(z,U_t.p)
                ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)

                
                
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

 


#%% 
# Definition of temporal scheme - Euler explicit with leap-frgo
def EE(f_dt, U_t, T_init, Tmax, z, propagation_type, *args):
    # we suppose here at initialisation that v in given at time t=0, and 
    # rho and p are given at time t=1/2
    if False:
        fig,ax = plt.subplots(3,1)
        fig.suptitle("Plot solution of "+propagation_type+ " propagation")
        ax[0].plot(z,U_t.rho, label="t = "+str(T_init))
        ax[1].plot(z,U_t.p)
        ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
        
    # load the dt
    dt = args[9]
    if Tmax >= T_init : 
        t_tab = np.arange(T_init,Tmax,dt)
    else : 
        t_tab = np.arange(T_init,Tmax,-dt)
    
    if propagation_type == "checkpointing":
        index_fr0 = np.arange(0,len(t_tab),1)[:-1:5]
        step_btw_fr0 = index_fr0[1] - index_fr0[0]
        Nfr_1 = 5
        q,r = divmod(step_btw_fr0, Nfr_1+1)
        index_fr1 = np.array([i * q + min(i,r) for i in range(Nfr_1+1)])
        print(index_fr1, index_fr0)
        checkpoint1 = []
        checkpoint0 = deepcopy(args[12][index_fr0])
        
        t_tab = np.append(t_tab, t_tab[-1]+dt)
        history = []
        
    else : 
        # create a vector to save state at each time
        history = [deepcopy(U_t)]  
        if propagation_type == "adjoint" : 
            #t_tab = t_tab = np.append(t_tab, t_tab[-1]-dt)
            propagation_type = "forward"
            
    for it, t in enumerate(t_tab):

        if propagation_type == "forward" :# need to change the order of treatment in the backward propagation 
            tn = f_dt(U_t, t, U_t, *args)
            U_t.v = deepcopy(tn.v)#U_t.v + tn.v
            tndemi = f_dt(U_t, t+dt/2, U_t, *args)
            U_t.rho = deepcopy(tndemi.rho)#U_t.rho + tndemi.rho
            U_t.p = deepcopy(tndemi.p)#U_t.p + tndemi.p 

        elif propagation_type == "backward" :
            tndemi = f_dt(U_t, t, U_t, *args)
            U_t.rho = deepcopy(tndemi.rho)#U_t.rho + tndemi.rho
            U_t.p = deepcopy(tndemi.p)#U_t.p + tndemi.p 
            tn = f_dt(U_t, t+dt/2, U_t, *args)
            U_t.v = deepcopy(tn.v)

        elif propagation_type == "checkpointing" :

            # get the nb of the previous frame and the nb of iteration between the previous frame and the actual time
            aux = np.where(len(t_tab) - 1 - it >= index_fr0)
            which_fr = len(aux[0])-1
            which_it_infr = len(t_tab) - 1 - it - index_fr0[aux[0][-1]]
            # get the number of the previous frame of the second chekpointing and the nb of iter
            aux = np.where(which_it_infr >= index_fr1)[0]
            which_fr1 = len(np.where(which_it_infr > index_fr1)[0])-1
            which_it_infr1 = which_it_infr - index_fr1[aux[-1]]
            
            # option 1 : first iteration, need to create the checkpointing1
            if it==0 :
                aux_U_t = deepcopy(checkpoint0[-1])
                checkpoint1 = []
                for it_fr1 in range(index_fr0[-1], len(t_tab)-1):
                    # iteration to record the checkpoint
                    tn = f_dt(aux_U_t, t_tab[it_fr1], aux_U_t, *args[:-1])
                    aux_U_t.v = deepcopy(tn.v)#aux_U_t.v = aux_U_t.v + tn.v 
                    tndemi = f_dt(aux_U_t, t_tab[it_fr1]+dt/2, aux_U_t, *args[:-1])
                    aux_U_t.rho = deepcopy(tndemi.rho)#aux_U_t.rho = aux_U_t.rho + tndemi.rho 
                    aux_U_t.p = deepcopy(tndemi.p)#aux_U_t.p = aux_U_t.p + tndemi.p
                    # save if it is at a checkpoint time
                    # note that "+1" because we get U^n+1
                    if (it_fr1 - index_fr0[-1] + 1) in index_fr1:
                        checkpoint1 += [deepcopy(aux_U_t)]
                # load U^N+1
                U_t = deepcopy(aux_U_t)
                
                  
            # option 2 : last iteration, already saved in checkpoint0
            elif  it == len(t_tab) - 1:
                # load U^0
                U_t = deepcopy(checkpoint0[which_fr])  
                              
            # option 3 : actual time has a value in checkpoint0
            elif which_it_infr == 0 :
                # load the value of time t : U^n
                U_t = checkpoint0[which_fr]
                #  need to load new checkpoint1
                aux_U_t = deepcopy(checkpoint0[which_fr-1])
                checkpoint1 = []
                for it_fr1 in range(index_fr0[which_fr-1], index_fr0[which_fr]-1):
                    tn = f_dt(aux_U_t, t_tab[it_fr1], aux_U_t, *args[:-1])
                    aux_U_t.v = deepcopy(tn.v)#aux_U_t.v = aux_U_t.v + tn.v 
                    tndemi = f_dt(aux_U_t, t_tab[it_fr1]+dt/2, aux_U_t, *args[:-1])
                    aux_U_t.rho = deepcopy(tndemi.rho)#aux_U_t.rho = aux_U_t.rho + tndemi.rho 
                    aux_U_t.p = deepcopy(tndemi.p)#aux_U_t.p = aux_U_t.p + tndemi.p 
                    # save if it is at a checkpoint time 
                    if (it_fr1 - index_fr0[which_fr-1] + 1) in index_fr1:
                        # note that "+1" because we get U^n+1
                        checkpoint1 += [deepcopy(aux_U_t)]
   
            # option 4 :  actual time has a value in checkpoint1
            elif which_it_infr1 == 0:
                U_t = deepcopy(checkpoint1[which_fr1])
                
            # option 5 : no checkpoint exist
            else :
                # load the more recent checkpoint
                if which_fr1 == 0 :
                    aux_U_t = deepcopy(checkpoint0[which_fr])
                else : 
                    aux_U_t = deepcopy(checkpoint1[which_fr1-1])   # -1 because we neglect the checkpoint in 0
                
                for it_fr1 in range(index_fr0[which_fr] + index_fr1[which_fr1], len(t_tab)-it-1):
                    # iterate until finding the actual time
                    tn = f_dt(aux_U_t, t_tab[it_fr1], aux_U_t, *args[:-1])
                    aux_U_t.v = deepcopy(tn.v)#aux_U_t.v = aux_U_t.v + tn.v 
                    tndemi = f_dt(aux_U_t, t_tab[it_fr1]+dt/2, aux_U_t, *args[:-1])
                    aux_U_t.rho = deepcopy(tndemi.rho)#aux_U_t.rho = aux_U_t.rho + tndemi.rho 
                    aux_U_t.p = deepcopy(tndemi.p)#aux_U_t.p = aux_U_t.p + tndemi.p 
                # save if it is at a checkpoint time
                U_t = deepcopy(aux_U_t)
           
#        if True : #not args[11] : 
#            U_t =  U_t + f(U_t, t, *args) * dt 
#        else : 
#            U_t =  U_t + f(args[12][len(args[12]) - 2 - int(round(t/dt))], t, *args) * dt   
#        if not args[11] :# need to change the order of treatment in the backward propagation 
#            tn = f(U_t, t, *args)
#            U_t.v = U_t.v + tn.v * dt
#            tndemi = f(U_t, t+dt/2, *args)
#            U_t.rho = U_t.rho + tndemi.rho * dt
#            U_t.p = U_t.p + tndemi.p * dt
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
            
#        #U_t = apply_sponge(U_t)
        
        if False:
            if (t+dt) % 5 < 5e-4 : 
                ax[0].plot(z,U_t.rho, label='t = '+str(t+dt))
                ax[1].plot(z,U_t.p)
                ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)

        history += [deepcopy(U_t)]
        
    if False:
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
