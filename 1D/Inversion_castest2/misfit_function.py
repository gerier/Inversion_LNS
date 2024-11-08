from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import * 

from Lib_NavierStokesEq.compute_kernels import *
from Lib_NavierStokesEq.checkpointing import *
from parametrisation import *

from parameters import obs_start,obs_end,h

import config




def f0(m0,history_events,norm_events,dt,h,size,gamma,g,sponge_layer,source_events,index_receivers, Sc,parametrisation):
    config.n_fx += 1
    # set the model in the parametrisation ["density", "wind","pressure"] to be adapted to implementation of Euler equation
    M = model2paramLNSeq(m0, parametrisation, Sc, gamma, g, size)

    
    Sum_CHI = 0 
    for event in range(len(history_events)):
        obs = history_events[event]
        source=source_events[event]
        norm_info = norm_events[event]

 
        # init
        CHI = 0
        # get Ucalc
        U1 = LNS_Variable(np.zeros(size), np.zeros(size-1), np.zeros(size)) 
        calc,_ = time_scheme(get_RHS, U1, 0, len(obs), dt, index_receivers, "forward", sponge_layer, M, h, source)
     
        if False : 
            fig,ax = plt.subplots(len(index_receivers),1, figsize=(8,5))
            t = np.arange(0,len(calc)*dt-dt,dt)

            for i in range(len(index_receivers)):
                ax[i].plot(t,[calc[tt][i] for tt in range(len(t))], 'b')
                ax[i].plot(t[:-1],[obs[tt][i] for tt in range(len(t)-1)], 'r')
                ax[i].set_yticks([])
                ax[i].set_ylabel("R"+str(i))
            ax[-1].set_xlabel("Time [s]")
            plt.show()
            fig,ax = plt.subplots(len(index_receivers),1, figsize=(8,5))
            for i in range(len(index_receivers)):
                ax[i].plot(t[:-1],[calc[tt][i] - obs[tt][i] for tt in range(len(t)-1)], 'm')
                ax[i].set_yticks([])
                ax[i].set_ylabel("R"+str(i))
            ax[-1].set_xlabel("Time [s]")
            plt.show()

            sys.exit()

        # Ucalc - Uobs / norm
        if norm_info['choice'] == 1 :
            # norm is sum_rec( sum_t (observations(t,r))) or max_rec( max_t (observations(t,r)))
            for ir in range(len(index_receivers)) :
                norm_obs = norm_info["norm_obs"]
                for it in range(len(obs)):
                    damping = apodisation(it*dt, 1.2/source[1],source[2])
                    CHI += (calc[it][ir] - obs[it][ir])**2 * damping / norm_obs**2

        elif norm_info['choice'] == 2 :
            if norm_info["calc"]:
                # norm is sum_t (signal(t,r))) or max_t (signal(t,r))), where signal is for synthetics or observations
                norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
                norm_obs = norm_info["norm_obs"]
                for ir in range(len(index_receivers)) :
                    for it in range(len(obs)):
                        damping = apodisation(it*dt, 1.2/source[1],source[2])
                        CHI += (calc[it][ir]/norm_calc[ir] - obs[it][ir]/norm_obs[ir])**2 
            else : 
                # norm is sum_t (observations(t,r))) or max_t (observations(t,r))),
                norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
                norm_obs = norm_info["norm_obs"]
                for ir in range(len(index_receivers)) :
                    for it in range(len(obs)):
                        damping = apodisation(it*dt, 1.2/source[1],source[2])
                        CHI += (calc[it][ir]- obs[it][ir])**2/norm_obs[ir]**2

        elif norm_info['choice'] == 3 :
            if norm_info["calc"]:
                # norm is  (signal(t,r))), where signal is for synthetics or observations
                norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
                norm_obs = norm_info["norm_obs"]
                for ir in range(len(index_receivers)) :
                    for it in range(len(obs)):
                        damping = apodisation(it*dt, 1.2/source[1],source[2])
                        CHI += (calc[it][ir]/norm_calc[it][ir] - obs[it][ir]/norm_obs[ir])**2 
            else : 
                # norm is (observations(t,r)))
                norm_calc = get_norm(calc,  norm_info['ord'], norm_info['choice'])
                norm_obs = norm_info["norm_obs"]
                for ir in range(len(index_receivers)) :
                    for it in range(len(obs)):
                        damping = apodisation(it*dt, 1.2/source[1],source[2])
                        CHI += ((calc[it][ir] - obs[it][ir])/norm_obs[it][ir])**2 

        CHI *= (dt/2)
        Sum_CHI += CHI
    return Sum_CHI   

def df0(m0, it_end, dt, h, size, sponge_layer, source_events, index_receivers, history_events, norm_events, gamma, g, Sc,parametrisation):
    config.n_grad += 1
    M = model2paramLNSeq(m0, parametrisation, Sc, gamma, g, size)

    sum_grad_event = np.zeros((3,size))
    for event in range(len(history_events)):

        source = source_events[event]
        obs = history_events[event]
        norm_info = norm_events[event]

        grad, preconditioner =  compute_gradient(M, 0, it_end, dt, h, size, sponge_layer, source, index_receivers, obs, norm_info)
        sum_grad_event += grad 

    flatten_grad = kernel2paramModel(sum_grad_event, M, parametrisation, Sc, gamma, boundary=True)
    return flatten_grad, preconditioner

hessf0 = None

def plot_gradient(grad, z, parametrisation, title):
    size = len(z)
    h = z[1] - z[0]
    z = z /1000

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
    ax[0].axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
    ax[0].axvspan((len(z) - 50) * h /1000, len(z) * h /1000, facecolor='k', alpha=0.15)
    ax[1].axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
    ax[1].axvspan((len(z) - 50) * h /1000, len(z) * h /1000, facecolor='k', alpha=0.15)
    ax[2].axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
    ax[2].axvspan((len(z) - 50) * h /1000, len(z) * h /1000, facecolor='k', alpha=0.15)    
    ax[0].set_xlim(z[0],z[-1])
    ax[1].set_xlim(z[0],z[-1])
    ax[2].set_xlim(z[0],z[-1])
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[2].set_xlabel("Range(km)")
    ax[0].set_title(title)


def plot_flatteniter(m0, z, parametrisation, Sc, gamma, title):

    h = z[1] - z[0]
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
    ax[0].axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
    ax[0].axvspan((len(z) - 50) * h /1000, len(z) * h /1000, facecolor='k', alpha=0.15)
    ax[1].axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
    ax[1].axvspan((len(z) - 50) * h /1000, len(z) * h /1000, facecolor='k', alpha=0.15)
    ax[2].axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
    ax[2].axvspan((len(z) - 50) * h /1000, len(z) * h /1000, facecolor='k', alpha=0.15)
    ax[3].axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
    ax[3].axvspan((len(z) - 50) * h /1000, len(z) * h /1000, facecolor='k', alpha=0.15)
    ax[0].set_xlim(z[0],z[-1])
    ax[1].set_xlim(z[0],z[-1])
    ax[2].set_xlim(z[0],z[-1])
    ax[3].set_xlim(z[0],z[-1])
    ax[0].set_title(title)
    ax[3].set_xlabel("Range(km)")
    #ax[0].set_xlim(0,40)
    #ax[1].set_xlim(0,40)
    #ax[2].set_xlim(0,40)
    #ax[3].set_xlim(0,40)


