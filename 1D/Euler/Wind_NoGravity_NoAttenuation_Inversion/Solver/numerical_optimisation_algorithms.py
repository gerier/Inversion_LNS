import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy



# TODO
# - define compute_kernels
# - define hessian approximation
# - choix de alpha_i dans alpha 1 et alpha max
# - donction phi et phiprime dans linesearch
# - fonction interpolate dans zoom


#%%

def compute_kernel():
    pass
def interpolate():
    pass

#%% Line search algorithm


def phi(alpha, f, xk,p):
    return f(xk + alpha * p)

def phiprime(alpha, xk, p):
    return np.sum(compute_kernel(xk + alpha * p) * p) 


def zoom(low, high, c1, c2, phi, phiprime):
    phizero = phi(0)
    phiprimezero = phiprime(0)

    while True : 
        philow = phi(low)

        alpha = interpolate(low,high)
        phialpha = phi(alpha)
        
        if phialpha > phizero + c1 * alpha * phiprimezero or phialpha >= philow : 
            high = alpha
        else : 
            phiprimealpha = phiprime(alpha)
            if abs(phiprimealpha) <= c2 * phiprimezero : 
                return alpha
            else : 
                high = low
            low = alpha


def line_search(alpha1, alphamax, phi, phiprime, c1=1e-4, c2=0.9):
    alpha = 0
    alpha_prec = 0 
    phizero = phi(0)
    phiprimezero = phiprime(0)

    i = 0
    while True : 
        phialpha = phi(alpha)
        if (phi(alpha) > phizero + c1 * alpha * phiprimezero) or (i > 1 and phialpha >= phialpha_prec):
            alpha = zoom(alpha_prec, alpha)
            return alpha
        phiprimealpha = phiprime(alpha)
        if abs(phiprimealpha) <= - c2 * phiprimezero : 
            return alpha
        if phiprimealpha >= 0 :
            alpha = zoom(alpha, alpha_prec)
            return alpha
        
        alpha_prec = alpha
        phialpha_prec = phialpha 
        alpha = None # choix entre alpha 1 et alpha max
        i = i + 1


def backtracking(rho,c1, phi, phiprime):
    phizero = phi(0)
    phiprimezero = phiprime(0)
    while phi(alpha) <= phizero + c1 * alpha * phiprimezero :
        alpha = rho * alpha
    return alpha 


#%% Algorithm

def steepest_descent_method(x0,itermax,tol):
    xk = deepcopy(x0)
    k = 0
    tol_condition = True
    while (k< itermax and not tol_condition):
        df = compute_kernel(xk)
        alpha = line_search(alpha)
        xk = xk + alpha * df
        k = k+1
        tol_condition = np.norm(df) < tol
    return xk


def nonlin_conjugate_gradient(x0, itermax,tol):
    xk = deepcopy(x0)
    current_df = compute_kernel(xk)
    p = deepcopy(current_df)
    k = 0
    tol_condition = True
    while (k< itermax and not tol_condition):
        alpha = line_search(alpha)
        xk = xk + alpha * p
        current_df = deepcopy(next_df) 
        next_df = compute_kernel(xk)
        beta =  sum( next_df * (next_df - current_df)) / np.norm(current_df)
        p = - next_df + beta * p
        k = k+1
        tol_condition = np.norm(next_df) < tol
    return xk



def hessian_approx(x1,x0, df1, df0):
    pass
    #return Hk

def BFGS(x0,itermax,tol):
    xk = deepcopy(x0)
    current_df = compute_kernel(xk)
    Hk = hessian_approx(xk)
    k = 0   
    tol_condition = True
    while (k< itermax and not tol_condition):
        p = - np.dot(Hk,current_df)
        alpha = line_search(alpha)
        sk = deepcopy(-xk)
        xk = xk + alpha * p
        sk = sk + xk
        yk = deepcopy(-next_df)
        next_df = compute_kernel(xk) 
        yk += next_df
        rhok = 1 / np.sum(yk * sk)
        Q = np.eye(len(yk)) - rhok * np.dot(sk,yk.T)
        Hk = np.dot(np.dot(Q,Hk),Q) + rhok * np.dot(sk,sk.T)
        k = k+1 
        tol_condition = np.norm(next_df) < tol
    return xk
