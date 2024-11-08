import numpy as np
from Lib_NavierStokesEq.linearised_navier_stokes import LNS_Model
import sys
import matplotlib.pyplot as plt

def interp_K_spongelayer(K, index_z0 = 50):
    # fn to replace the values of the kernel in the sponge layer
    # extend the value of the first node that is not is the sponge
    n,m = np.shape(K)
    for k in range(n):
        K[k,:index_z0] = K[k,index_z0]
        K[k,m-index_z0:] = K[k,m-index_z0-1]
    return K

def kernel2paramModel(kernel, M, parametrisation, Sc, gamma, boundary=True):
    # Input : 
    # K:                array of sensitvity kernel (size 3 x n)
    # M:                current model (in aprametrisation ["density", "wind","pressure"])
    # parametrisation:  list of names of each entry, in which we want to have the kernel, example: ["density", "wind", "velocity"]
    if parametrisation == ["density", "wind", "pressure"]:
        flatten_kernel= np.append(np.append(kernel[0,:], kernel[2,:]), kernel[1,:-1])
        flatten_kernel *= Sc
    elif parametrisation == ["density", "wind", "velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = M.p / M.rho * kernel[2,:] + kernel[0,:]  
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = 2 * c * M.rho / gamma * kernel[2,:]

        new_kernel = interp_K_spongelayer(new_kernel)

        flatten_kernel= np.append(np.append(new_kernel[0,:], new_kernel[2,:]), new_kernel[1,:-1])
        flatten_kernel *= Sc
    elif parametrisation == ["log_density", "wind", "log_pressure"]:

        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = kernel[0,:]*M.rho
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = kernel[2,:]*M.p

        new_kernel = interp_K_spongelayer(new_kernel)

        flatten_kernel= np.append(np.append(new_kernel[0,:], new_kernel[2,:]), new_kernel[1,:-1])
        flatten_kernel *= Sc   

    elif parametrisation == ["log_density", "wind", "velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = M.p / M.rho * kernel[2,:] + kernel[0,:]  
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = 2 * c * M.rho / gamma * kernel[2,:]

        new_kernel = interp_K_spongelayer(new_kernel)

        flatten_kernel= np.append(np.append(new_kernel[0,:]*M.rho, new_kernel[2,:]), new_kernel[1,:-1])
        flatten_kernel *= Sc
    elif parametrisation == ["log_density", "wind", "log_velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = M.p / M.rho * kernel[2,:] + kernel[0,:] 
        new_kernel[0,:] = new_kernel[0,:]*M.rho
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = 2 * c * M.rho / gamma * kernel[2,:]
        new_kernel[2,:] = new_kernel[2,:]*c

        new_kernel = interp_K_spongelayer(new_kernel)

        flatten_kernel= np.append(np.append(new_kernel[0,:], new_kernel[2,:]), new_kernel[1,:-1])
        flatten_kernel *= Sc
    elif parametrisation == ["log_pressure", "wind", "velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = - 2* M.rho / c * kernel[0,:]  
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = kernel[2,:] + M.rho / M.p * kernel[0,:] 

        new_kernel = interp_K_spongelayer(new_kernel)

        flatten_kernel= np.append(np.append(new_kernel[2,:]*M.p, new_kernel[0,:]), new_kernel[1,:-1])
        flatten_kernel *= Sc
    elif parametrisation == ["log_pressure", "wind", "log_velocity"]:
        c = np.sqrt( M.p * gamma / M.rho)
        new_kernel = np.zeros_like(kernel)
        new_kernel[0,:] = - 2* M.rho / c * kernel[0,:] 
        new_kernel[0,:] = new_kernel[0,:]*c
        new_kernel[1,:] = kernel[1,:]
        new_kernel[2,:] = kernel[2,:] + M.rho / M.p * kernel[0,:] 
        new_kernel[2,:] = new_kernel[2,:]*M.p

        new_kernel = interp_K_spongelayer(new_kernel)

        flatten_kernel= np.append(np.append(new_kernel[2,:], new_kernel[0,:]), new_kernel[1,:-1])
        flatten_kernel *= Sc
    else:
        print("Error: parametrisation unknown")
        sys.exit(-1)
    return flatten_kernel


def model2paramLNSeq(m0, parametrisation, Sc, gamma, g, size):
    if parametrisation == ["density", "wind", "pressure"] : 
        m = Sc * m0
        rho = m[:size]
        p = m[size:2*size]
        v = m[2*size:]
        M = LNS_Model(rho, v, p, gamma,g)
    elif parametrisation == ["density", "wind", "velocity"]:
        m = Sc * m0
        rho = m[:size]
        c = m[size:2*size]
        v = m[2*size:]
        p = rho * c**2 / gamma
        M = LNS_Model(rho, v, p, gamma, g)
    elif parametrisation == ["log_density", "wind", "log_pressure"] : 
        m = Sc * m0
        rho = np.exp(m[:size]) 
        p = np.exp(m[size:2*size])
        v = m[2*size:] 
        M = LNS_Model(rho, v, p, gamma, g)
    elif parametrisation == ["log_density", "wind", "velocity"]:
        m = Sc * m0
        rho = np.exp(m[:size])
        c = m[size:2*size]
        v = m[2*size:]
        p = rho * c**2 / gamma
        M = LNS_Model(rho, v, p, gamma, g)
    elif parametrisation == ["log_density", "wind", "log_velocity"]:
        m = Sc * m0
        rho = np.exp(m[:size])
        c = np.exp(m[size:2*size])
        v = m[2*size:]
        p = rho * c**2 / gamma
        M = LNS_Model(rho, v, p, gamma, g)
    elif parametrisation == ["log_pressure", "wind", "velocity"]:
        m = Sc * m0
        p = np.exp(m[:size])
        c = m[size:2*size]
        v = m[2*size:]
        rho = gamma * p / c**2
        M = LNS_Model(rho, v, p, gamma, g)
    elif parametrisation == ["log_pressure", "wind", "log_velocity"]:
        m = Sc * m0
        p = np.exp(m[:size])
        c = np.exp(m[size:2*size])
        v = m[2*size:]
        rho = gamma * p / c**2
        M = LNS_Model(rho, v, p, gamma, g)
    else:
        print("Error: parametrisation unknown")
        sys.exit(-1)
    return M

def model2allparam(m0, parametrisation,Sc, gamma, size):
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
    elif parametrisation == ["log_pressure", "wind", "velocity"]:
        m = Sc * m0
        p = np.exp(m[:size])
        c = m[size:2*size]
        v = m[2*size:]
        rho = p * gamma / c**2 
    elif parametrisation == ["log_pressure", "wind", "log_velocity"]:
        m = Sc * m0
        p = np.exp(m[:size])
        c = np.exp(m[size:2*size])
        v = m[2*size:]
        rho = p * gamma /  c**2 
    else:
        print("Error: parametrisation unknown")
        sys.exit(-1)
    return rho, p, c, v
