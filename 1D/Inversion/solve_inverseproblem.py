import numpy as np
import matplotlib.pyplot as plt

from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import * 
from parameters import *
from Lib_NavierStokesEq.compute_kernels import *
from Lib_NavierStokesEq.checkpointing import *

from optimisation_lib_reg import *

# Define the observation that we want to reproduce
import time

# for the figure
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})


it_end = int(Tmax / dt)
size = len(z)

start = time.time()
U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size))
history_obs, U1 = time_scheme(get_RHS, U1, 0, it_end, dt, index_receivers, "forward", GT_M, h, source)

max_obs_p = [ max([abs(history_obs[t][k]) for t in range(len(history_obs))]) for k in range(len(index_receivers))]
max_obs = [np.max(max_obs_p)**2]


K = compute_gradient(M, 0, it_end, dt, h, size, source, index_receivers, history_obs, max_obs)

print(time.time() - start)
plt.figure()
plt.plot(z/1000,K[0,:])
plt.show()


def f0(m0,obs,max_obs,dt,size,gamma,source,index_receivers, Sc):
    m = Sc * m0
    rho = m[:size]
    p = m[size:2*size]
    v = m[2*size:]
    M = LNS_Model(rho, v, p, gamma)

    # init
    CHI = 0
    # get Ucalc
    U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size)) 
    calc,_ = time_scheme(get_RHS, U1, 0, len(obs), dt, index_receivers, "forward", M, h, source)
    # Ucalc - Uobs
    for r in range(len(index_receivers)) :
        for t in range(len(obs)):
            CHI += (calc[t][r] - obs[t][r])**2 / max_obs

    CHI *= (dt/2)
    return CHI   

def df0(m0, it_end, dt, h, size, source, index_receivers, obs, max_obs, gamma, Sc):
    m = Sc * m0
    rho = m[:size]
    p = m[size:2*size]
    v = m[2*size:]
    M = LNS_Model(rho, v, p, gamma)
    grad =  compute_gradient(M, 0, it_end, dt, h, size, source, index_receivers, obs, max_obs)
    flatten_grad = np.append(np.append(grad[0,:], grad[2,:]), grad[1,:-1])
    return flatten_grad * Sc

def hessf0(m, k, m_prec=None, k_prec=None):
    if True or k_prec is None : 
        h = np.ones_like(k)
    else : 
        h = (k - k_prec + 1e-40) / (m- m_prec + 1e-40)
    return h

def plot_flatteniter(m0, z,title):
    m = Sc * m0
    z = z / 1000 
    size = int((len(m)+ 1)/3)
    rho = m[:size]
    p = m[size:2*size]
    v = m[2*size:]
    M = LNS_Model(rho, v, p, gamma)
    c = np.sqrt(p * 1.4 / rho)
    #M.plot(z, title)
    fig,ax = plt.subplots(4,1, figsize=(10, 7))
    ax[0].plot(z,rho)
    ax[1].plot(z,p)
    ax[2].plot(z[:-1],v)
    ax[3].plot(z,c)
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Pressure")
    ax[2].set_ylabel("Velocity")
    ax[3].set_ylabel("Celerity")
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[3].grid()
    ax[1].axvspan(17, 18, facecolor='k', alpha=0.15)
    ax[3].axvspan(17, 18, facecolor='k', alpha=0.15)
    ax[0].set_title(title)
    ax[3].set_xlabel("Range(km)")
    plt.show()


Sc = [1] * len(z) + [10**5] * len(z) + [10] * (len(z)-1)

argf = [history_obs, max_obs, dt, size, gamma,source,index_receivers, Sc]
argg = [it_end, dt, h, size, source, index_receivers, history_obs, max_obs,gamma, Sc]

m = np.append(np.append(M.rho, M.p), M.v)
m = m / Sc


x, is_notfound, x_iter = optimisation(m, f0, df0, hessf0, 0.0001, 0.9, argf, argg, alpha_max=20, maxiter=100, tol_x=1e-8, tol_df=1e-8, true_hessian=False, regularisation=None, plot_iter=[plot_flatteniter,z])

x = Sc * x

print("_"*16)
print("%10s %10s %10s %10s"%("Is solution","F(m)", "DF(m)", "Nb ITER"))
print("%10s %10.6f %10.6f %10d"%(str(is_notfound[0]),is_notfound[1], is_notfound[2], is_notfound[3]))
print("_"*16)




plt.plot(z/1000, np.sqrt(x[len(z):2*len(z)] * M.gamma / x[:len(z)]))
plt.axvspan(17, 18, facecolor='k', alpha=0.15)
plt.ylabel("Celerity (m/s)")
plt.grid()
plt.show()


##########################7
# Peut penser faire une moyenne glissante sur quelques points pour lissrr le signal 