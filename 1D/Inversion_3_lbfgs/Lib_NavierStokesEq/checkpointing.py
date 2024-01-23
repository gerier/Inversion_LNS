from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import * 
from parameters import zero_as_initialcondition, time_scheme
from copy import deepcopy
from get_norm import *
from parameters import N_frames

def save_all_frames(N_frames, M, it_start, it_end, dt, h, size_x, sponge_layer, source, index_receivers, norm, norm_info=None):

    Frames = []
    Nstep = it_end - it_start
    step_between_frames = Nstep / N_frames - 1


    if zero_as_initialcondition :
        Ut = LNS_Variable(np.zeros(size_x), np.zeros(size_x-1), np.zeros(size_x))
    else : 
        Ut = np.load("./initialcondition.npy", allow_pickle=True)[0]
        source[-1] = []
        
    it_end = 0
    norm_calc = None

    for i in range(N_frames):
        it_start = it_end  
        it_end = it_end + step_between_frames + 1
    
        # _, Ut = EE(get_RHS, Ut, it_start, it_end, dt, index_receivers, "forward", M, h, source)
        history, Ut = time_scheme(get_RHS, Ut, it_start, it_end, dt, index_receivers, "forward", sponge_layer, M, h, source)
       
        if i < N_frames-1:
            Frames += [deepcopy(Ut)]

        # if needs to compute the norm, save all the part of the history in one variable : history_calc
        if norm :
            if i == 0:
                history_calc = deepcopy(history)
            else : 
                history_calc += history
    
    # compute the norm
    if norm:
        norm_calc = get_norm(history_calc, norm_info['ord'], norm_info['choice'])

    return Frames, Ut, norm_calc





def compute_from_last_frame(it_last, it_frame, Frame, M, dt, h, size_x, sponge_layer, source, index_receivers) : 

    if it_last  == it_frame : 
        Ut = Frame
    else : 
        if it_frame == 0 :
            Frame = LNS_Variable(np.zeros(size_x), np.zeros(size_x-1), np.zeros(size_x))

        _, Ut = time_scheme(get_RHS, deepcopy(Frame), it_frame, it_last, dt, index_receivers, "forward", sponge_layer, M, h, source)
    return Ut
