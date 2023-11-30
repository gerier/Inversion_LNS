import numpy as np
from Lib_NavierStokesEq.linearised_navier_stokes import LNS_Model
import sys

def kernel2paramModel(kernel, M, parametrisation, gamma):
    # Input : 
    # K:                array of sensitvity kernel (size 3 x n)
    # M:                current model (in aprametrisation ["density", "wind","pressure"])
    # parametrisation:  list of names of each entry, in which we want to have the kernel, example: ["density", "wind", "velocity"]
    if parametrisation == ["density", "wind", "pressure"]:
        flatten_kernel= np.append(np.append(kernel[0,:], kernel[2,:]), kernel[1,:-1])

    elif parametrisation == ["density", "wind", "velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = M.p / M.rho * kernel[2,:] + kernel[0,:]  
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = 2 * c * M.rho / gamma * kernel[2,:]
        flatten_kernel= np.append(np.append(new_kernel[0,:], new_kernel[2,:]), new_kernel[1,:-1])
    
    elif parametrisation == ["log_density", "wind", "log_pressure"]:
        flatten_kernel= np.append(np.append(kernel[0,:]*M.rho, kernel[2,:]*M.p), kernel[1,:-1])
    
    elif parametrisation == ["log_density", "log_wind", "log_pressure"]: # impossible car vent pas que positif 
        flatten_kernel= np.append(np.append(kernel[0,:]*M.rho, kernel[2,:]*M.p), kernel[1,:-1]*M.v)
    
    elif parametrisation == ["log_density", "wind", "velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = M.p / M.rho * kernel[2,:] + kernel[0,:]  
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = 2 * c * M.rho / gamma * kernel[2,:]
        flatten_kernel= np.append(np.append(new_kernel[0,:]*M.rho, new_kernel[2,:]), new_kernel[1,:-1])

    elif parametrisation == ["log_density", "wind", "log_velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = M.p / M.rho * kernel[2,:] + kernel[0,:]  
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = 2 * c * M.rho / gamma * kernel[2,:]
        flatten_kernel= np.append(np.append(new_kernel[0,:]*M.rho, new_kernel[2,:]*c), new_kernel[1,:-1])

    elif parametrisation == ["log_density", "log_wind", "log_velocity"]:  # impossible car vent pas que positif 
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = M.p / M.rho * kernel[2,:] + kernel[0,:]  
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = 2 * c * M.rho / gamma * kernel[2,:]
        flatten_kernel= np.append(np.append(new_kernel[0,:]*M.rho, new_kernel[2,:]*c), new_kernel[1,:-1]*M.v)
   
    else:
        print("Error: parametrisation unknown")
        sys.exit(-1)
    return flatten_kernel


def model2paramLNSeq(m0, parametrisation, Sc, gamma, size):
    if parametrisation == ["density", "wind", "pressure"] : 
        m = Sc * m0
        rho = m[:size]
        p = m[size:2*size]
        v = m[2*size:]
        M = LNS_Model(rho, v, p, gamma)
    elif parametrisation == ["density", "wind", "velocity"]:
        m = Sc * m0
        rho = m[:size]
        c = m[size:2*size]
        v = m[2*size:]
        p = rho * c**2 / gamma
        M = LNS_Model(rho, v, p, gamma)
    elif parametrisation == ["log_density", "wind", "log_pressure"] : 
        m = Sc * m0
        rho = np.exp(m[:size])
        p = np.exp(m[size:2*size])
        v = m[2*size:]
        M = LNS_Model(rho, v, p, gamma)
    elif parametrisation == ["log_density", "log_wind", "log_pressure"] :  # impossible car vent pas que positif 
        m = Sc * m0
        rho = np.exp(m[:size])
        p = np.exp(m[size:2*size])
        v = np.exp(m[2*size:])
        M = LNS_Model(rho, v, p, gamma)
    elif parametrisation == ["log_density", "wind", "velocity"]:
        m = Sc * m0
        rho = np.exp(m[:size])
        c = m[size:2*size]
        v = m[2*size:]
        p = rho * c**2 / gamma
        M = LNS_Model(rho, v, p, gamma)
    elif parametrisation == ["log_density", "wind", "log_velocity"]:
        m = Sc * m0
        rho = np.exp(m[:size])
        c = np.exp(m[size:2*size])
        v = m[2*size:]
        p = rho * c**2 / gamma
        M = LNS_Model(rho, v, p, gamma)
    elif parametrisation == ["log_density", "log_wind", "log_velocity"]:  # impossible car vent pas que positif 
        m = Sc * m0
        rho = np.exp(m[:size])
        c =  np.exp(m[size:2*size])
        v =  np.exp(m[2*size:])
        p = rho * c**2 / gamma
        M = LNS_Model(rho, v, p, gamma)
    else:
        print("Error: parametrisation unknown")
        sys.exit(-1)
    return M

def model2allparam(m0, parametrisation,Sc, gamma,size):
    if parametrisation == ["density", "wind", "pressure"] : 
        m = Sc * m0
        rho = m[:size]
        p = m[size:2*size]
        v = m[2*size:]
        c = np.sqrt(p * gamma/ rho)
    elif parametrisation == ["density", "wind", "velocity"]:
        m = Sc * m0
        rho = m[:size]
        c = m[size:2*size]
        v = m[2*size:]
        p = rho * c**2 / gamma
    elif parametrisation == ["log_density", "wind", "log_pressure"] : 
        m = Sc * m0
        rho = np.exp(m[:size])
        p = np.exp(m[size:2*size])
        v = m[2*size:]
        c = np.sqrt(p * gamma/ rho)
    elif parametrisation == ["log_density", "log_wind", "log_pressure"] : # impossible car vent pas que positif 
        m = Sc * m0
        rho = np.exp(m[:size])
        p = np.exp(m[size:2*size])
        v = np.exp(m[2*size:])
        c = np.sqrt(p * gamma/ rho)
    elif parametrisation == ["log_density", "wind", "velocity"]:
        m = Sc * m0
        rho = np.exp(m[:size])
        c = m[size:2*size]
        v = m[2*size:]
        p = rho * c**2 / gamma
    elif parametrisation == ["log_density", "wind", "log_velocity"]:
        m = Sc * m0
        rho = np.exp(m[:size])
        c = np.exp(m[size:2*size])
        v = m[2*size:]
        p = rho * c**2 / gamma
    elif parametrisation == ["log_density", "log_wind", "log_velocity"]: # impossible car vent pas que positif 
        m = Sc * m0
        rho = np.exp(m[:size])
        c = np.exp(m[size:2*size])
        v = np.exp(m[2*size:])
        p = rho * c**2 / gamma
    else:
        print("Error: parametrisation unknown")
        sys.exit(-1)
    return rho, p, c, v