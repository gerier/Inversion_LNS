#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Librairies to import
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from copy import deepcopy
import pandas as pd
from scipy.interpolate import CubicSpline


# In[2]:


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

    def __div__(self,alpha):
        self.rho /= alpha
        self.v /= alpha
        self.p /= alpha
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



# In[3]:


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


# In[4]:


# numerical scheme
def DF_C_2(U, dz, BC, need_BC=True):
    if need_BC:
        BC = [BC[1],BC[-2]]
    else : 
        BC=  [ [], []]
    new_U = np.append(np.append([BC[0]], U), [BC[1]])
    diff_U = new_U[1:] - new_U[:-1]
    return diff_U / dz
'''   
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
'''    

def DF_DC_2_backward(U0, dz, BC):
    U = np.append([BC[0]], U0)
    U = U[1:] - U[:-1]
    return U / dz

def DF_DC_2_forward(U0, dz, BC):
    U = np.append(U0, [BC[-1]])
    U = U[1:] - U[:-1]
    return U / dz


def DF_DC_4_backward(U0, dz, BC):
    U = np.append(np.append(BC[:2], U0), [BC[-2]])
    Umm = U[:-3]
    Um = U[1:-2]
    Up = U[3:]
    U = Umm/6 + Um/2 - U[2:-1] + Up/3
    return U / dz

def DF_DC_4_forward(U0, dz, BC):
    U = np.append(np.append([BC[2]], U0), BC[-2:])
    Upp = U[3:]
    Up = U[2:-1]
    Um = U[:-3]
    U = - Upp/6 - Up/2 + U[1:-2] - Um/3
    return U / dz

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
'''
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
'''
def interpolation(U):
    # compute the value at point i+1/2 for i in [0,N-1]
    BC = [U[0],U[-1]]
    U = np.append(np.append([BC[0]], U), [BC[1]])
    U = U[1:] + U[:-1]
    return U / 2


# In[5]:


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


def DF_Sigma_C_C(U1, U0, dz, gamma, is_reverse, DF_C):
    # initialise
    S  = np.zeros((3,len(U1.rho)))
    
    BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2                       
    BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
    BC_rho0 = [U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]  # Points -2, -1, N+1, N+2 
    BC_p0 = [U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
    
    # compute the contribution on density
    S[0,:] = U0.rho * DF_C(U1.v, dz, BC1)                                              # rho0 div v1
    S[0,:] +=  U1.rho * DF_C(U0.v, dz, BC_v0)                                          # rho1 div v0
    # compute the contribution on qunatity of mvt
    S[1,:-1] = DF_C(U1.p, dz, BC1, False)                                                 # grap p1
    # compute the contribution on pressure
    S[2,:] = gamma * U1.p * DF_C(U0.v, dz, BC_v0)                                  # gamma p1 div v0
    S[2,:] +=  gamma * U0.p * DF_C(U1.v, dz, BC1)                                  # gamma p0 div v1
    
    return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])


def DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC):
    # initialise
    S  = np.zeros((3,len(U1.rho)))
    
    BC1 = [0,0,0,0]                                           # Points -2, -1, N+1, N+2                       
    BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
    BC_rho0 = [U0.rho[0], U0.rho[0], U0.rho[-1], U0.rho[-1]]  # Points -2, -1, N+1, N+2 
    BC_p0 = [U0.p[0], U0.p[0], U0.p[-1], U0.p[-1]]            # Points -2, -1, N+1, N+2 
    
    # compute the contribution on density
    S[0,:] = interpolation(U1.v) * DF_DC(U0.rho, dz, BC_rho0)  # v1 grad rho0
    S[0,:] += interpolation(U0.v) * DF_DC(U1.rho, dz, BC1)    # v0 grad rho1 
    # compute the contribution on qunatity of mvt
    S[1,:-1] = (U0.rho[1:]+U0.rho[:-1])/2 * U0.v * DF_DC(U1.v, dz, BC1)    # roh0 v0 . div v1
    # compute the contribution on pressure
    S[2,:] = interpolation(U1.v) * DF_DC(U0.p, dz, BC_p0)  # v0 grad p1               
    S[2,:] += interpolation(U0.v) * DF_DC(U1.p, dz, BC1)  # v1 grad p0
    
    return LNS_Variable(S[0,:], 2 * S[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), S[2,:])


def DF_Sigma_D_C(U1, U0, dz, T0, g, l, mu, kappa, gamma, R, DF_C):
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



def G_C(U1, U0, dz, is_reverse, DF_C):
    # initialise
    GU = np.zeros((3,len(U1.rho)))
    rho1_demi = (U1.rho[1:]+U1.rho[:-1])/2
    g_demi = (g[1:]+g[:-1])/2
    rho0_demi = (U0.rho[1:]+U0.rho[:-1])/2

    # compute the contribution on density
    GU[0,:] = 0
    # compute the contribution on quantity of mvt
    GU[1,:-1] = rho1_demi * g_demi                            # rho1 g
    # compute the contribution on pressure
    GU[2,:] = 0
    
    return LNS_Variable(GU[0,:], GU[1,:-1] / rho0_demi, GU[2,:])


def G_DC(U1, U0, dz, is_reverse, DF_DC):
    # initialise
    GU = np.zeros((3,len(U1.rho)))
    BC_v0 = [U0.v[0], U0.v[0], U0.v[-1], U0.v[-1]]            # Points -2, -1, N+1, N+2 
    rho0_demi = (U0.rho[1:]+U0.rho[:-1])/2
    rho1_demi = (U1.rho[1:]+U1.rho[:-1])/2
    
    # compute the contribution on density
    GU[0,:] = 0
    # compute the contribution on quantity of mvt
    GU[1,:-1] = - (rho0_demi * U1.v + rho1_demi * U0.v) * DF_DC(U0.v, dz, BC_v0)  # (rho0v1 + rho1v0) div v0
    # compute the contribution on pressure
    GU[2,:] = 0
    
    return LNS_Variable(GU[0,:], GU[1,:-1] / rho0_demi, GU[2,:])


def F(U1, U0, it, source, n):
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
    return LNS_Variable(f[0,:], 2 * f[1,:-1] / (U0.rho[1:]+U0.rho[:-1]), f[2,:])


# In[6]:


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
        
    RHS_c = F(U1, U0, int(abs(t/dt)), source, len(U1.rho)) + G_C(U1, U0, dz, is_reverse, DF_C) \
        + DF_Sigma_D_C(U1, U0, dz, T0, g, l, mu, kappa, gamma, R, DF_C) \
        - DF_Sigma_C_C(U1, U0, dz, gamma, is_reverse, DF_C)
    
    RHS_dbackward = RHS_c + (G_DC(U1, U0, dz, is_reverse, DF_DC_backward) + DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_backward))*0.5
    
    RHS_dforward = RHS_dbackward + (G_DC(U1, U0, dz, is_reverse, DF_DC_forward) + DF_Sigma_C_DC(U1, U0, dz, gamma, is_reverse, DF_DC_forward))*0.5 
    
    return RHS_dforward


# In[7]:


# Definition of temporal scheme - Runge Kutta of order 4
def RK4(f, U_t, T_init, Tmax,*args):
    fig,ax = plt.subplots(3,1, figsize=(10,7))
    ax[0].plot(z,U_t.rho, label=str(T_init))
    ax[1].plot(z,U_t.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    # Load the dt
    dt = args[9]
    args = list(args) + [False] # add the value of is_reverse
    if Tmax < T_init:
        dt *= -1
        args[-1] = True 
        
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
    return t, U_t, history

# Definition of temporal scheme - Euler explicit
def EE(f, U_t, T_init, Tmax, *args):
    fig,ax = plt.subplots(3,1)
    ax[0].plot(z,U_t.rho, label=str(T_init))
    ax[1].plot(z,U_t.p)
    ax[2].plot( (z[1:] + z[:-1])/2,U_t.v)
    
    # load the dt
    dt = args[9]
    args = list(args) + [False] # add the value of is_reverse
    if Tmax < T_init:
        dt *= -1
        args[-1] = True 
    
    # create a vector to save state at each time
    history = [deepcopy(U_t)]
    
    for t in np.arange(T_init,Tmax,dt):
        U_t = U_t +  f(U_t, t, *args) * dt
        U_t = apply_sponge(U_t)
        if (t+dt) % 5 < 1e-4 : 
            ax[0].plot(z,U_t.rho, label=str(t+dt))
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
    return t, U_t, history


# # Definition of parameters

# In[8]:


# source
f0 = 0.1
# space
z0 = 0
zmax = 20e3
# time 
T_init = 0
Tmax = 10

# user paramter
display_anim = False
case2 = False

# # Definition de la source


# In[9]:


def get_source(t,f0):
    t0 = 1.2/f0
    return np.exp(-4 * np.pi**2 * f0**2 * (t-t0)**2)


# # Cas 1

# In[10]:


# The parameters

# define model parameters
rho0 = 1.04898036 
v0 = -4.33914709 #0.00745178154
p0 =88714.4844 
T0 = 294.375305 
c = 344.108887 
g = 0*9.81084824 
mu = 1.27685234e-05
kappa = 0.0258137602 
gamma = 1.40011787 
eta = 2.72986326e-05 
cv = 20.7801247 
M = 28.965 # masse molaire de l'air https://fr.wikipedia.org/wiki/Air

Cv = cv/M
l = eta - (2/3)*mu

param_toset_onmesh = [rho0, v0, p0, T0, c, g, mu, kappa, gamma, eta, Cv]


# ### Identify the mesh limits

# In[11]:


min_v = min(abs(c), abs(v0))
max_v = max(abs(c), abs(v0))
max_dx = min_v / (10 * 2.5 * f0) 

max_dt = max_dx / max_v

print("With this source, you must choose a dx < ", str(max_dx) )

print("If taking the dx max, finaly, you will have to choose a dt < ", str(max_dt))


# ### Define the mesh

# In[12]:


# The spatio-temporel domain

# define the space
h = 100      # size of the mesh
z = np.arange(z0,zmax+h,h)


# define receivers
d_receivers = 3.6e3 # distance between 2 receivers, 
                  # must be a multiple of h
r = np.arange(z0+d_receivers, zmax, d_receivers)

# define Time
t = 0
dt = 0.005


# ### Define the model on the mesh

# In[13]:


for i,param in enumerate(param_toset_onmesh):
    param_toset_onmesh[i] = param * np.ones(len(z))

[rho0, v0, p0, T0, c, g, mu, kappa, gamma, eta, Cv] = param_toset_onmesh


# ### Define the source on the temporal mesh

# In[14]:


# set the source
iz_source = 100
t_ax = np.arange(T_init,Tmax+2*dt,dt)
source = get_source(t_ax, f0) #scipy.stats.norm.pdf(t_ax,10,1.5)

dist_z = abs(z - z[iz_source])
factor_source = np.exp(-dist_z/1000)

plt.plot(t_ax, source)
plt.title("Form of the source applied on density on point z=12km")
plt.show()

# plot the source in time and space
total_source = np.zeros((len(z), len(t_ax)))
for t in t_ax : 
    for i in range(len(z)):
        total_source[i,int(t/dt)] = source[int(t/dt)] * factor_source[i]


fig, ax = plt.subplots()
plt.imshow(total_source, aspect='auto')
plt.colorbar()
plt.title("Form of the source applied on density on point z="+ str(z[iz_source])+"km")
plt.xticks(range(0,len(t_ax),1000), t_ax[0:len(t_ax):1000])
plt.yticks(range(0,len(z),20), z[0:len(z):20]/1000)
plt.xlabel("Time")
plt.ylabel("Altitude")

source = [source, factor_source]


# ### Initialisation

# In[15]:


# define vectors of the system
v0_demi = (v0[1:] + v0[:-1])/2
U0 = LNS_Variable(rho0, v0_demi, p0) 
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 


# In[16]:


# Plot the model
fig,ax = plt.subplots(3,1, figsize=(10,7))
ax[0].plot(z,U0.rho)
ax[1].plot(z,U0.p)
ax[2].plot( (z[1:] + z[:-1])/2,U0.v)
ax[0].set_xlabel("Altitude")
ax[1].set_xlabel("Altitude")
ax[2].set_xlabel("Altitude")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("A T = 0, the background")


# #### Resolve the LNS

# In[17]:


# Resolution in a 1d case
t_end, U_end, history_obs = RK4(get_RHS, U1, T_init, Tmax, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source)


# #### Plot the final state

# In[18]:


fig,ax = plt.subplots(3,1, figsize=(10,7))
fig
ax[0].plot(z,U_end.rho)
ax[1].plot(z,U_end.p)
ax[2].plot( (z[1:] + z[:-1])/2,U_end.v)
ax[0].set_xlabel("Altitude")
ax[1].set_xlabel("Altitude")
ax[2].set_xlabel("Altitude")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Perturbation after Tmax = "+ str(Tmax))


# #### Modify the right hand sided to backpropagate and get again the initial state

# In[19]:


# Definition of resolution spatial
#def get_minus_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source):
#    return - get_RHS(U1, t, U0, T0, g, l, mu, kappa, gamma, R, dz, dt, source)


# In[20]:


# Resolution in a 1d case
U_tmax = deepcopy(U_end)
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 
fig,ax = plt.subplots(3,1, figsize=(10,7))
fig
ax[0].plot(z,U1.rho)
ax[1].plot(z,U1.p)
ax[2].plot( (z[1:] + z[:-1])/2,U1.v)
ax[0].set_xlabel("Altitude")
ax[1].set_xlabel("Altitude")
ax[2].set_xlabel("Altitude")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Perturbation after Tmax = "+ str(Tmax))
t_start, U_start, history_reverse = RK4(get_RHS, U1, Tmax, T_init, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source)


# In[21]:


fig,ax = plt.subplots(3,1, figsize=(10,7)
                     )
ax[0].plot(z,U_start.rho)
ax[1].plot(z,U_start.p)
ax[2].plot( (z[1:] + z[:-1])/2,U_start.v)
ax[0].set_xlabel("Altitude")
ax[1].set_xlabel("Altitude")
ax[2].set_xlabel("Altitude")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("After Tmax = 0")


# # Test 2: the model is not constant on z axis

