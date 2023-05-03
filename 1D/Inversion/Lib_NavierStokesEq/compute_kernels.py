import numpy as np
from parameters import *
from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import *
import sys
sys.path.insert(1, './Lib_NavierStokesEq/')
from  checkpointing import *
import math

def compute_gradient(M, it_start, it_end, dt, h, x_size, source, index_receivers, history_obs, norm_obs):
    Frame, Ureverse = save_all_frames(N_frames, M, 0, it_end, dt, h, x_size, source, index_receivers)
    U_start_adjoint = LNS_Variable(np.zeros(x_size), np.zeros(x_size-1), np.zeros(x_size)) 

    nstep = it_end - it_start
    step_between_frames =  nstep / N_frames - 1

    backprop_prec = deepcopy(Ureverse)
    adjoint_prec = U_start_adjoint

    K = np.zeros((4, len(adjoint_prec.rho)))
    for it in range(1,nstep):
        # compute backtracking state
        i_frame = int(math.floor((nstep-it) // (step_between_frames + 1 ))) - 1
        it_frame = (i_frame + 1) * (step_between_frames + 1 ) 
        reverse = compute_from_last_frame(nstep-it, it_frame, Frame[i_frame], M, dt, h, x_size, source, index_receivers)
        
        # compute adjoint state
        history_reverse = backprop_prec.p[index_receivers]
   
        #test += [history_obs[nstep-it] - history_reverse ]
        _, adjoint = time_scheme(get_adjoint_RHS, U_start_adjoint, it, it+1, dt, index_receivers, "adjoint", M, h, history_reverse, history_obs[-it], norm_obs)
        U_start_adjoint = deepcopy(adjoint)
   
        # compute kernel
        source_with_time = [t] + source
        K += get_kernels(adjoint, reverse, adjoint_prec, backprop_prec, M, dt, h, source_with_time) 

        adjoint_prec = deepcopy(adjoint)
        backprop_prec = deepcopy(reverse)
 
    return K