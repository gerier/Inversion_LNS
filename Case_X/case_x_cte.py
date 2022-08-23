#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:25:07 2022

@author: s.gerier
"""


from linearised_navier_stokes import *
from discretisation import *

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from copy import deepcopy
import pandas as pd
from scipy.interpolate import CubicSpline
from copy import deepcopy



#%% definition of parameters

# source
f0 = 0.1
# space
z0 = 0
zmax = 20e3
# time 
T_init = 0
Tmax = 20

# user paramter
display_anim = False
case2 = False

#%% define the model thanks to a MSISE file

path_file = "/home/deos/s.gerier/PROJECTS/SIMULATIONS/ATMOSPHERIC_MODELS/MSISE/msisehwm_wrapper/OUTPUT/Flores_atmosphere_500km_86p2_4.dat"
path_file = "./Flores_atmosphere_500km_86p2_4.dat"
data = np.genfromtxt(path_file, skip_header=3)#,delimiter=[7,7,7,7,7,7,7,3,3,3,3,5,5,7,7,7])
df = pd.DataFrame(data, columns =["z[m]", "rho[kg/(m^3)]", "T[K]", "c[m/s]", "p[Pa]", "H[m]", "g[m/(s^2)]", "N^2[rad^2/s^2]", "kappa[J/(s.m.K)]", "mu[kg(s.m)]", "mu_vol[kg/(s.m)]", "w_M[m/s]", "w_Z[m/s]", "w_P[m/s]", "c_p[J/(mol.K)]", "c_v[J/(mol.K)]", "gamma"])
v0 = 0

physical_parameters = ["rho[kg/(m^3)]", "T[K]", "c[m/s]", "p[Pa]", "g[m/(s^2)]", "kappa[J/(s.m.K)]", "mu[kg(s.m)]", "mu_vol[kg/(s.m)]", "c_v[J/(mol.K)]", "gamma"]



#%% identify limits of mesh

min_v = min(abs(df["c[m/s]"]))
max_v = max(abs(df["c[m/s]"]))
max_dx = min_v / (10 * 2.5 * f0) 

max_dt = max_dx / max_v

print("With this source, you must choose a dx < ", str(max_dx) )
print("If taking the dx max, finaly, you will have to choose a dt < ", str(max_dt))


#%% MESH PARAMETERS

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

#%% define the mesh

# apply model on the mesh
interp_physical_param = []
for param in physical_parameters:
    if param == "rho[kg/(m^3)]" or param == "p[Pa]" : 
        cs = CubicSpline(df["z[m]"], np.log(df[param]))
        interp_physical_param += [np.exp(cs(z))]
    else : 
        cs = CubicSpline(df["z[m]"], df[param])
        interp_physical_param += [cs(z)]

[rho0, T0, c, p0, g, kappa, mu, eta, cv, gamma] = interp_physical_param
v0 = v0 * np.ones(len(z))
M = 28.965 # masse molaire de l'air https://fr.wikipedia.org/wiki/Air

l = eta - (2/3)*mu
Cv = cv/M

# define vectors of the system
v0_demi = (v0[1:] + v0[:-1])/2
U0 = LNS_Variable(rho0, v0_demi, p0) 
U1 = LNS_Variable(np.zeros(len(z)), np.zeros(len(z)-1), np.zeros(len(z))) 



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
fig.suptitle("The background")
plt.show()
plt.close()
#%% DEFINE SOURCE ON MESH

# set the source
iz_source = 100
t_ax = np.arange(T_init,Tmax+dt+dt,dt)
source = get_source(t_ax, f0) #scipy.stats.norm.pdf(t_ax,10,1.5)

dist_z = abs(z - z[iz_source])
factor_source = np.exp(-dist_z/1000)

source = np.ones(len(source))
source[500:] = 1

plt.figure()
plt.plot(t_ax, source)
plt.title("Form of the source applied on density on point z=12km")
plt.show()
plt.close()

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
plt.close()

source = [source, factor_source]
reversed_source = [np.flip(source[0])[1:], source[1]]

#%% Resolution in a 1d case
t_end, U_end, history_obs = EE(get_RHS, U1, T_init, Tmax, z, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, False)



fig,ax = plt.subplots(3,1, figsize=(10,7))
fig
ax[0].plot(z,U_end.rho)
ax[1].plot(z,U_end.p)

ax[2].plot( (z[1:] + z[:-1])/2,U_end.v)
ax[0].set_xlabel("Distance on axis x (km)")
ax[1].set_xlabel("Distance on axis x (km)")
ax[2].set_xlabel("Distance on axis x (km)")
ax[0].set_ylabel("Density")
ax[1].set_ylabel("Pressure")
ax[2].set_ylabel("Velocity")
ax[0].grid()
ax[1].grid()
ax[2].grid()
fig.suptitle("Perturbation after Tmax = "+ str(Tmax))
plt.close()

#get_ipython().run_line_magic('matplotlib', 'inline')
from IPython.display import HTML
from matplotlib import animation, rc

def get_anim(parameter_toplot, history, limits, z, propagation) : 
    list_anim = []

    fig,ax = plt.subplots(figsize=(10,3))
    ax.set_xlabel("Altitude")
    ax.set_ylabel(parameter_toplot)

    ax.grid()
    ax.set_xlim(( z[0], z[-1]))
    ax.set_ylim(limits)

    line0, = ax.plot([], [], lw=2)

    # initialization function: plot the background of each frame
    def init():
        line0.set_data([], [])
        return (line0,)

    # animation function. This is called sequentially
    def animate(i,z,history, parameter_toplot):
        x = z
        if parameter_toplot == "Density" :
            y = history[i].rho
        elif parameter_toplot == "Velocity" :
            y = history[i].v  
            x = (z[1:] + z[:-1])/2
        elif parameter_toplot == "Pressure" :
            y = history[i].p 
        line0.set_data(x, y)
        return (line0,)

    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(z,history,parameter_toplot), frames=len(history), interval=20, blit=True)
    #writervideo = animation.FFMpegWriter(fps=60)
    anim.save('./animation_'+parameter_toplot+'_'+propagation+".mp4")#, writer=writervideo)
    plt.close()

disp_anim = False
if disp_anim : 
    get_anim("Density", history_obs[::10], (-0.02,0.02), z, 'forward')
    get_anim("Velocity", history_obs[::10], (-4,4), z, 'forward')
    get_anim("Pressure", history_obs[::10], (-2000,2000), z, 'forward')

#%% BACKWARD RESOLUTION


U_tmax = deepcopy(U_end)
t_start, U_start, history_reverse = EE(get_minus_RHS, U_end, T_init, Tmax, z, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, reversed_source, False)


# Plot the model
fig,ax = plt.subplots(3,1, figsize=(10,7))
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
fig.suptitle("A T = 0, the background")
plt.close()

if disp_anim:
    get_anim("Density", history_reverse[::10], (-0.02,0.02), z, 'backward')
    get_anim("Velocity", history_reverse[::10], (-4,4), z, 'backward')
    get_anim("Pressure", history_reverse[::10], (-2000,2000), z, 'backward')
    
    
#%% SUPERIMPOSITION OF BACKWARD AND FORWARD SOLUTION

n = len(history_reverse)
for i in range(0,n,int(n/3)):
    fig,ax = plt.subplots(3,1, figsize=(10,7))
    # forward
    ax[0].plot(z,history_obs[i].rho, 'b', label='Forward')
    ax[1].plot(z,history_obs[i].p, 'b')
    ax[2].plot( (z[1:] + z[:-1])/2,history_obs[i].v, 'b')
    # backward
    ax[0].plot(z,history_reverse[n-1-i].rho, 'r', label='backward')
    ax[1].plot(z,history_reverse[n-1-i].p, 'r')
    ax[2].plot( (z[1:] + z[:-1])/2,history_reverse[n-1-i].v, 'r')
    # element to display 
    ax[0].legend()
    ax[0].set_xlabel("Altitude")
    ax[1].set_xlabel("Altitude")
    ax[2].set_xlabel("Altitude")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    fig.suptitle("Comparison of forward and backward solution at time (forward) t = "+str(t_ax[i]))
    plt.show()
    
    
#%% RESOLUTION ADJOINT EQUATIONS

if False :                                       
	max_obs_rho = max( [max(history_obs[t].rho) for t in range(len(history_obs))])
	max_obs_v = max( [max(history_obs[t].v) for t in range(len(history_obs))])
	max_obs_p = max( [max(history_obs[t].p) for t in range(len(history_obs))])
	
	max_obs = [max_obs_rho, max_obs_v, max_obs_v]
	
	U_end = deepcopy(U_tmax)
	t_start, U_start, history_adjoint = RK4(get_adjoint_RHS, U_end, Tmax, T_init, U0, T0, g, l, mu, kappa, gamma, Cv, h, dt, source, history_reverse, history_obs, d_receivers, max_obs)
	
	# Plot the model
	fig,ax = plt.subplots(3,1, figsize=(10,7))
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
	fig.suptitle("A T = 0, the background")
