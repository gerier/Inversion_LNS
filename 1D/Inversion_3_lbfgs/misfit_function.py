from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import * 

from Lib_NavierStokesEq.compute_kernels import *
from Lib_NavierStokesEq.checkpointing import *
from parametrisation import *

from parameters import obs_start,obs_end,h

import config

def f0(m0,obs,norm_info,dt,h,size,gamma,sponge_layer,source,index_receivers, Sc,parametrisation):
    config.n_fx += 1

    # set the model in the parametrisation ["density", "wind","pressure"] to be adapted to implementation of Euler equation
    M = model2paramLNSeq(m0, parametrisation, Sc, gamma, size)

    # init
    CHI = 0
    # get Ucalc
    U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size)) 
    calc,_ = time_scheme(get_RHS, U1, 0, len(obs), dt, index_receivers, "forward", sponge_layer, M, h, source)
    # Ucalc - Uobs / norm
    if norm_info['choice'] == 1 :
        # norm is sum_rec( sum_t (observations(t,r))) or max_rec( max_t (observations(t,r)))
        for r in range(len(index_receivers)) :
            norm_obs = norm_info["norm_obs"]
            for t in range(len(obs)):
                CHI += (calc[t][r] - obs[t][r])**2 / norm_obs**2

    elif norm_info['choice'] == 2 :
        if norm_info["calc"]:
            # norm is sum_t (signal(t,r))) or max_t (signal(t,r))), where signal is for synthetics or observations
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for r in range(len(index_receivers)) :
                for t in range(len(obs)):
                    CHI += (calc[t][r]/norm_calc[r] - obs[t][r]/norm_obs[r])**2 
        else : 
            # norm is sum_t (observations(t,r))) or max_t (observations(t,r))),
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for r in range(len(index_receivers)) :
                for t in range(len(obs)):
                    CHI += (calc[t][r]- obs[t][r])**2/norm_obs[r]**2

    elif norm_info['choice'] == 3 :
        if norm_info["calc"]:
            # norm is  (signal(t,r))), where signal is for synthetics or observations
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for r in range(len(index_receivers)) :
                for t in range(len(obs)):
                    CHI += (calc[t][r]/norm_calc[t][r] - obs[t][r]/norm_obs[r])**2 
        else : 
            # norm is (observations(t,r)))
            norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
            norm_obs = norm_info["norm_obs"]
            for r in range(len(index_receivers)) :
                for t in range(len(obs)):
                    CHI += ((calc[t][r] - obs[t][r])/norm_obs[t][r])**2 

    CHI *= (dt/2)
    return CHI   

def df0(m0, it_end, dt, h, size, sponge_layer, source, index_receivers, obs, norm_info, gamma, Sc,parametrisation):
    config.n_grad += 1

    M = model2paramLNSeq(m0, parametrisation, Sc, gamma, size)
    grad, preconditioner =  compute_gradient(M, 0, it_end, dt, h, size, sponge_layer, source, index_receivers, obs, norm_info)
    flatten_grad = kernel2paramModel(grad, M, parametrisation, gamma)
    return flatten_grad * Sc, preconditioner

# def hessf0(m, k, nstep, dt, index_receivers, dz, gamma, source, m_prec=None, k_prec=None):
#     config.n_h += 1
#     if False and k_prec is None : 
#         h = np.ones_like(k)
#     elif True : 
#         size = int((len(m)+1)/3)
#         sc_m = Sc * m
#         rho = sc_m[:size]
#         p = sc_m[size:2*size]
#         v = sc_m[2*size:]
#         M = LNS_Model(rho, v, p, gamma)

#         U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size)) 
#         _,_,h0_demi, h0 = time_scheme(get_RHS, U1, 0, nstep, dt, index_receivers, "forward", M, dz, source, compute_h0=True)

#         h = np.zeros_like(m)
#         h[:size] = h0
#         h[size:2*size] = h0
#         h[2*size:] = h0_demi
#         #h = h * Sc

        
#     else : 
#         h = (k - k_prec + 1e-40) / (m- m_prec + 1e-40)
#     return h

def hessf0(m, preconditioner,size):
    config.n_h += 1
    
    interp_preconditioner = interpolation(preconditioner)

    h = np.zeros_like(m)
    h[:size] = interp_preconditioner
    h[size:2*size] = interp_preconditioner
    h[2*size:] = preconditioner

    h = 0* h + 1 

    return h

def plot_gradient(grad, z, parametrisation, title):
    size = len(z)

    grad_param1 = grad[:size]
    grad_param2 = grad[size:2*size]
    grad_param3 = grad[2*size:]

    fig, ax = plt.subplots(3,1, figsize=(10, 7))
    ax[0].plot(z,grad_param1)
    ax[1].plot(z,grad_param2)
    ax[2].plot(z[:-1],grad_param3)
    ax[0].set_ylabel(parametrisation[0])
    ax[1].set_ylabel(parametrisation[2])
    ax[2].set_ylabel(parametrisation[1])
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[2].set_xlabel("Range(km)")
    ax[0].set_title(title)

def plot_hessian(h, z, parametrisation, title):
    size = len(z)

    h_param1 = h[:size]
    h_param2 = h[size:2*size]
    h_param3 = h[2*size:]

    fig, ax = plt.subplots(3,1, figsize=(10, 7))
    ax[0].plot(z,h_param1)
    ax[1].plot(z,h_param2)
    ax[2].plot(z[:-1],h_param3)
    ax[0].set_ylabel(parametrisation[0])
    ax[1].set_ylabel(parametrisation[2])
    ax[2].set_ylabel(parametrisation[1])
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[2].set_xlabel("Range(km)")
    ax[0].set_title(title)

def plot_flatteniter(m0, z, parametrisation, Sc, gamma, title):

    z = z / 1000 
    size = int((len(m0)+ 1)/3)

    # set the model parameter depending on the chosen parametrisation 
    rho,p,c,v = model2allparam(m0, parametrisation,Sc, gamma,size)

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
    ax[1].axvspan(obs_start * h /1000, obs_end * h /1000, facecolor='k', alpha=0.15)
    ax[3].axvspan(obs_start * h /1000, obs_end * h /1000, facecolor='k', alpha=0.15)
    ax[0].set_title(title)
    ax[3].set_xlabel("Range(km)")
    ax[0].set_xlim(0,40)
    ax[1].set_xlim(0,40)
    ax[2].set_xlim(0,40)
    ax[3].set_xlim(0,40)


