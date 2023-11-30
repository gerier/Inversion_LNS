import numpy as np

from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import *
import sys
sys.path.insert(1, './Lib_NavierStokesEq/')
from  checkpointing import *
import math
from scipy import signal



def lissage(signal_brut,L):
    res = np.copy(signal_brut) # duplication des valeurs
    for i in range (1,len(signal_brut)-1): # toutes les valeurs sauf la première et la dernière
        L_g = min(i,L) # nombre de valeurs disponibles à gauche
        L_d = min(len(signal_brut)-i-1,L) # nombre de valeurs disponibles à droite
        Li=min(L_g,L_d)
        res[i]=np.sum(signal_brut[i-Li:i+Li+1])/(2*Li+1)
    return res


def lissage_butter(sig):
    # Fréquence d'échantillonnage
    fe = 1  # Hz
    # Fréquence de nyquist
    f_nyq = fe / 2.  # Hz
    # Fréquence de coupure
    fc = 0.02  # Hz
    # Préparation du filtre de Butterworth en passe-bas
    b, a = signal.butter(4, fc/f_nyq, 'low', analog=False)
    # Application du filtre
    new_K_v = signal.filtfilt(b, a, sig)
    return new_K_v
    
def lissage_mean_exp(sig):
    # Moyenne exponentielle mobile
    # Constante du système
    N = 20
    alpha = 2/(N+1)
    # Préparation de la liste de sortie
    s_m = []
    s_m.append(0)
    # Application du filtre
    for e in sig :
        s_m.append(alpha*e+(1-alpha)*s_m[-1])
    return s_m[1:]



def interp_K_spongelayer(K, index_z0 = 50):
    # fn to replace the values of the kernel in the sponge layer
    # extend the value of the first node that is not is the sponge
    n,m = np.shape(K)
    for k in range(n):
        K[k,:index_z0] = K[k,index_z0]
        K[k,m-index_z0:] = K[k,m-index_z0-1]
    return K



def compute_gradient(M, it_start, it_end, dt, h, x_size, sponge_layer, source, index_receivers, history_obs, norm_info):
    Frame, Ureverse, norm_calc = save_all_frames(N_frames, M, 0, it_end, dt, h, x_size, sponge_layer, source, index_receivers, norm_info["calc"], norm_info)
    norm_info["norm_calc"] = norm_calc
    U_start_adjoint = LNS_Variable(np.zeros(x_size), np.zeros(x_size-1), np.zeros(x_size)) 

    nstep = it_end - it_start
    step_between_frames =  nstep / N_frames - 1

    backprop_prec = deepcopy(Ureverse)
    adjoint_prec = U_start_adjoint

    preconditioner = np.zeros_like(M.v)

    K = np.zeros((4, len(adjoint_prec.rho)))
    for it in range(1,nstep):
        # compute backtracking state
        i_frame = int(math.floor((nstep-it) // (step_between_frames + 1 ))) - 1
        it_frame = (i_frame + 1) * (step_between_frames + 1 ) 
        reverse = compute_from_last_frame(nstep-it, it_frame, Frame[i_frame], M, dt, h, x_size, sponge_layer, source, index_receivers)
        
        # compute adjoint state
        history_reverse = backprop_prec.p[index_receivers]

        _, adjoint = time_scheme(get_adjoint_RHS, U_start_adjoint, it, it+1, dt, index_receivers, "adjoint", sponge_layer, M, h, history_reverse, history_obs[-it], norm_info)
        U_start_adjoint = deepcopy(adjoint)

        # compute kernel
        source_with_time = [it_start + it*dt] + source
        K += get_kernels(adjoint, reverse, adjoint_prec, backprop_prec, M, dt, h, source_with_time) 


        #preconditioner += (reverse.v - backprop_prec.v) * (adjoint.v - adjoint_prec.v) / dt
        preconditioner += (reverse.v)**2 * dt

        adjoint_prec = deepcopy(adjoint)
        backprop_prec = deepcopy(reverse)
 

    #for i in range(3):
    #    K[i,:] = lissage_mean_exp(K[i,:])

    # do not use the information on the sponge layer 
    K = interp_K_spongelayer(K)
    preconditioner = preconditioner.reshape(1,len(preconditioner))
    preconditioner = interp_K_spongelayer(preconditioner)

    return K, preconditioner[0,:]
