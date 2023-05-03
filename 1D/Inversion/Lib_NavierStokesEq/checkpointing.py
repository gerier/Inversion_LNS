from Lib_NavierStokesEq.discretisation import *
from Lib_NavierStokesEq.linearised_navier_stokes import * 

def save_all_frames(N_frames, M, it_start, it_end, dt, h, size_x, source, index_receivers):

    Frames = []
    Nstep = it_end - it_start
    step_between_frames = Nstep / N_frames - 1

    Ut = LNS_Variable(np.zeros(size_x), np.zeros(size_x-1), np.zeros(size_x))

    it_end = 0
    for i in range(N_frames):
        it_start = it_end  
        it_end = it_end + step_between_frames + 1
    
        _, Ut = EE(get_RHS, Ut, it_start, it_end, dt, index_receivers, "forward", M, h, source)

        if i < N_frames-1:
            Frames += [deepcopy(Ut)]


    return Frames, Ut





def compute_from_last_frame(it_last, it_frame, Frame, M, dt, h, size_x, source, index_receivers) : 

    if it_last  == it_frame : 
        Ut = Frame
    else : 
        if it_frame == 0 :
            Frame = LNS_Variable(np.zeros(size_x), np.zeros(size_x-1), np.zeros(size_x))

        _, Ut = EE(get_RHS, deepcopy(Frame), it_frame, it_last, dt, index_receivers, "forward", M, h, source)
    return Ut