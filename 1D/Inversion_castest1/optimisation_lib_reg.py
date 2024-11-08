import numpy as np
import numpy.linalg as npl
import sys
import matplotlib.pyplot as plt
import config
from copy import deepcopy
from matplotlib import ticker, cm
from scipy import signal
from Lib_NavierStokesEq.discretisation import interpolation,apply_laplacian_regularisation,get_derivative_laplacian_regularisation, apply_gradient_regularisation, get_derivative_gradient_regularisation
from parametrisation import model2allparam

import config

def wolf1_cond(fx, dfx_r, fx_new, alpha, c1):
    return fx_new <= fx + alpha * c1 * dfx_r

def wolf2strong_cond(r, dfx_r, dfx_new, c2):
    return abs(sum(dfx_new * r)) <= c2 * abs(dfx_r)



#############################################################

def lissage(signal_brut,L):
    signal_ext = np.zeros(len(signal_brut)+8)
    signal_ext[:4] = signal_brut[0]
    signal_ext[-4:] = signal_brut[-1]
    signal_ext[4:-4] = signal_brut
    return signal.medfilt(signal_ext,9)[4:-4]

def update(x, alpha, r):
    size = int((len(x)+1)/3)
    x_new = x + alpha * r

    x_new[0:size] = lissage(x_new[0:size], 5)
    x_new[size:2*size] = lissage(x_new[size:2*size], 5)
    x_new[2*size:] = lissage(x_new[2*size:], 5)

    freeze_around_source = False
    if freeze_around_source : 
        x_new[600-10:600+10] = x[600-10:600+10]
        x_new[600-10+size:600+size+10] = x[size+600-10:size+600+10]
        x_new[2*size+600-10:2*size+600+10] = x[2*size+600-10:2*size+600+10]
    
    return x_new


#############################################################
def best_alpha(array):
    best = 0
    f_best = array[0][2]

    for i in range(1,len(array)):
        if array[i][2] < f_best :
            f_best = array[i][2]
            best = i
    return array[best]

#############################################################
def backtracking(alpha, x, r, fx, dfx_r, f, arg_f, c1, rate = 0.8, maxiter = 100):
    stop = False
    iter = 0
    array_alpha = []

    while not stop  and iter < maxiter:
        config.n_iter_alpha += 1
        #x_new = x+alpha*r
        x_new = update(x, alpha, r)
        fx_new = f(x_new, arg_f)
        if not(np.isnan(fx_new)) :  
                array_alpha += [[alpha, x_new, fx_new]]
        if wolf1_cond(fx, dfx_r, fx_new, alpha, c1):
            stop = True
        else :
            alpha = rate * alpha

        iter +=1

    if iter == maxiter :
        print("ERROR : No acceptable step find in backtracking")
        [alpha, x_new, fx_new] = best_alpha(array_alpha)
        print("Choose alpha = ", alpha)

    return alpha, x_new, fx_new

#############################################################

def quadratic(a, fa, fpa, b, fb):
  D = fa
  C = fpa
  db = b - a
  B = (fb - D - C * db) / (db ** 2)
  xmin = a - C / (2. * B)
  return xmin

def cubic(a,fa,fpa,b,fb,fpb): 
    
    d1 = fpa + fpb - 3 * (fa-fb)/(a-b)
    d2 = np.sign(b-a) * np.sqrt( d1**2 - fpa*fpb)

    xmin = b - (b-a) * (fpb + d2 - d1)/(fpb - fpa + 2*d2)
    return xmin


#############################################################

def zoom(alpha_low, alpha_high, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2):
    stop = False
    iter = 0

    #x_low = x + alpha_low * r
    x_low = update(x, alpha_low, r)
    fx_low = f(x_low, arg_f)
    dfx_low, preconditioner = gradf(x_low, arg_g)

    #x_high = x + alpha_high * r
    x_high = update(x, alpha_high, r)
    fx_high = f(x_high, arg_f)
    dfx_high, preconditioner = gradf(x_high, arg_g)

    if np.isnan(fx_high) or np.isinf(fx_high):
        output = backtracking(alpha_high, x, r, fx, dfx_r, f, arg_f, c1)
        config.is_restarted = 2
        return *output, gradf(output[1], arg_g)[0]

    print(alpha_high, fx_high)


    print("START : alpha_high = ", alpha_high, fx_high)
    array_alpha = []

    if False : #not wolf1_cond(fx, dfx_r, fx_low, alpha_low, c1):
        print("ERROR : alpha low ne respecte pas 1e condition")
        print("Alpha_low = ", alpha_low)
        alpha, x_new, fx_new = backtracking(1, x, r, fx, dfx_r, f, arg_f, c1)
        dfx_new, preconditioner = gradf(x_new, arg_g)
        stop = True

        # sys.exit()

    while not stop:
        print("je suis dans zoom")
        # interpolate an alpha (need to be checked)
        if iter == 0 :
            alpha  = quadratic(alpha_low, fx_low, sum(dfx_low*r), alpha_high, fx_high)
        else :
            if (iter % 10< 1e-4):
                alpha = (alpha_high + alpha_low) / 2
            else : 
                #alpha = cubic(x_low, fx_low, sum(dfx_low*r), x_high, fx_high, x_old, fx_old)
                alpha = cubic(alpha_low, fx_low, sum(dfx_low*r), alpha_high, fx_high, sum(dfx_high*r))

        # check if the new alpha is in the interval defined by alpha_low, alpha_high
        if (alpha < alpha_low or alpha > alpha_high) or alpha == 0 :
            print(alpha, alpha_low, alpha_high)
            alpha = (alpha_high + alpha_low) / 2

        if np.isnan(alpha) :
            alpha = (alpha_high + alpha_low) / 2# backtracking(alpha_high, x, r, fx, dfx_r, f, arg_f, c1)
            #config.is_restarted = 2
            #return *output, gradf(output[1], arg_g)[0]
        
        print( "Zoom step :", alpha, alpha_low, alpha_high)


       
        #x_new = x + alpha * r
        x_new = update(x,alpha,r)
        fx_new = f(x_new, arg_f)
        dfx_new, preconditioner = gradf(x_new, arg_g)
        array_alpha += [[alpha, x_new, fx_new, dfx_new]]

        print(fx_new, fx_low, fx_high)
        print(alpha_low, fx_low, sum(dfx_low*r), alpha_high, fx_high, sum(dfx_high*r))

        if not wolf1_cond(fx, dfx_r, fx_new, alpha, c1) or (fx_new >= fx_low):
            print("je ne fais que modifier alpha high")
            alpha_high = alpha
            x_high = x_new
            fx_high = fx_new

        else   : 
            # add the second part of the condition because can have an intervalle where function only increase but with always a too strong condition
            if wolf2strong_cond(r, dfx_r, dfx_new, c2) or alpha_low == alpha:
                stop = True
            else :
                # need to understand why
                if sum(dfx_new * r) * (alpha_high - alpha_low)  >= 0 :
                    alpha_high = alpha_low
                    x_high = x_low
                    fx_high = fx_low

                alpha_low = alpha
                x_low = x_new
                fx_low = fx_new
                dfx_low = dfx_new
                
                if alpha_low < 1e-16 :
                    return  best_alpha(array_alpha)

        if abs(alpha_high - alpha_low) < 1e-10  :
            return  best_alpha(array_alpha)
                
        x_old = x_new
        fx_old = fx_new
        iter += 1 

    return alpha, x_new, fx_new, dfx_new


##############################################################

def linesearch(alpha_start, alpha_max, x, r ,fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2):
    alpha_prec = 0
    fx_prec = np.Inf
    alpha = alpha_start

    stop = False
    while not stop :
        config.n_iter_alpha += 1
        
        x_new = update(x, alpha, r)
        fx_new = f(x_new, arg_f)

        if (not  wolf1_cond(fx, dfx_r, fx_new, alpha, c1)) or (fx_new >= fx_prec):
            alpha_low = alpha_prec
            alpha_high = alpha
            print( alpha_prec, alpha,  wolf1_cond(fx, dfx_r, fx_new, alpha, c1), fx_new, fx_prec, fx)

            alpha, x_new, fx_new, dfx_new = zoom(alpha_low, alpha_high, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            stop = True
        else :
            dfx_new, preconditioner = gradf(x_new, arg_g)
            if wolf2strong_cond(r, dfx_r, dfx_new, c2):
                stop = True
            else :
                if sum(dfx_new * r) >= 0 :
                    alpha_high = alpha_prec
                    alpha_low = alpha

                    alpha, x_new, fx_new, dfx_new = zoom(alpha_low, alpha_high, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
                    stop = True
                else :
                    # equivalent à une methode de dichotomie 
                    alpha_prec = alpha
                    fx_prec = fx_new
                    alpha = (alpha + alpha_max) / 2
                    print(alpha, alpha_prec, alpha_max)
                    if (abs(alpha - alpha_max)) <1e-6 : 
                        stop = True

    return alpha, x_new, fx_new, dfx_new

def test_linesearch(alpha_start, alpha_max, x, r, fx, dfx, f, gradf, c1, c2, attendu):
    dfx_r = sum(r * dfx)
    new_alpha, x_new, fx_new, dfx_new = linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)

    print("Attendu : ", attendu)
    print("Obtenu : ", new_alpha)


##################################################################


def hessian_approximation(H_old, x_old, x, dfx_old, dfx):
    size = len(x)
    y = dfx - dfx_old
    s = x - x_old
    rho = 1 / sum(y * s)

    s = s.reshape((size,1))
    y = y.reshape((size,1))
    aux_mat = np.eye(size) - rho * np.dot(s, y.T)
    H = np.dot(aux_mat, np.dot(H_old, aux_mat)) + rho * np.dot(s,s.T) # 6.17 nocedal
    return H

def get_direction_descent(option, x, fx,dfx,x_old,fx_old,dfx_old,cond_dfx_old,r_old, iter, h_old=None, true_hessian=False):
    cond_dfx = None
    if option == 1 : # Fletcher-Reeves
        if true_hessian :
            H = hessf(x, dfx, x_old, dfx_old)
            cond_dfx = np.dot(npl.inv(H), dfx)
        else :
            H = np.ones(len(x)) #H = hessf(x, preconditioner, argh)
            index_near_0 = (H < 1e-16)
            cond_dfx = (dfx + 1e-16) / (H + 1e-16)
            cond_dfx[index_near_0] = dfx[index_near_0]
            
        if (iter % len(x) == 0 ) : 
            print("Gradients are far from orthogonal. Restarting gradients") #iter % 100 == 0  :
            beta = 0
        else :
            beta = sum(dfx * cond_dfx) / sum(dfx_old * cond_dfx_old)

        r = - cond_dfx + beta * r_old
  
        h = None

    elif option == 2: # Polak-Ribière-Polyak
        y = dfx - dfx_old
        
        beta =  sum(dfx * y) / sum(dfx_old * dfx_old)

        r =  - dfx + beta * r_old
        h = None
    elif option == 22 : # Polak-Ribière-Polyak with preconditionner
        y = dfx - dfx_old
        s = x - x_old

        if h_old is None : 
            s = x - x_old
            h_old = sum(y*y) / sum(s*y)
        
        h = h_old + y*y/sum(y*s) - (h_old*s)**2/sum(s * (h_old*s))
        rho = sum(y * (h*y)) / sum(y * s)

        M = np.diag(1/ (rho * h))


        beta =  sum(np.diag(M) * dfx * y) / sum(dfx_old * np.diag(M) * dfx_old)

        r =  - np.diag(M)* dfx + beta * r_old
    

    elif option == 3 : # Perry-Shanno
        y = dfx - dfx_old
        s = x - x_old

        beta = ( sum(dfx*y) - 2 * sum(y*y) / sum(y*r_old) * sum(dfx*r_old)) / sum(y * r_old)
        gamma = sum(dfx * r_old) / sum(y * r_old)
        tau = sum(y * s) / sum(y*y) 

        r = tau * (-dfx + gamma * y + beta * r_old)
        h = None
    
    elif option == 4 : # Hager Zhang
        y = dfx - dfx_old
        s = x - x_old

        beta = ( sum(dfx*y) - 2 * sum(y*y) / sum(y*r_old) * sum(dfx*r_old)) / sum(y * r_old)
        eta = -1 / np.sqrt(sum(r_old*r_old)) / min(0.01, np.sqrt(sum(dfx_old * dfx_old)))
        beta_plus = max(beta,eta)

        r = -dfx + beta_plus * r_old
        h = None
    elif option == 44: # Hager Zhang
        y = dfx - dfx_old
        s = x - x_old

        if h_old is None : 
            h_old = sum(y*y) / sum(s*y)

        h = h_old + y*y/sum(y*s) - (h_old*s)**2/sum(s * (h_old*s))
        rho = sum(y * (h*y)) / sum(y * s)

        M = np.diag(1/ (rho * h))

        beta = (sum(dfx*np.diag(M)*y) - 2 *  sum(y*np.diag(M)*y) / sum(y*r_old) * sum(dfx*r_old)) / sum(y*r_old)
        eta = -1 / np.sqrt(sum(r_old*r_old)) / min(0.01, np.sqrt(sum(dfx_old * dfx_old)))
        beta_plus = max(beta,eta)

        r = -np.diag(M)*dfx + beta_plus * r_old

    elif option ==5:
        y = dfx - dfx_old
        s = x - x_old

        eta = 0.5
        beta = sum(dfx*y)/sum(y*r_old) - sum(y*y) / sum(s*y) * sum(dfx*s)/sum(y*r_old)
        beta_plus = max( beta, eta * sum(dfx*r_old) / sum(r_old*r_old))

        r = - dfx + beta_plus * r_old
        h = None

    elif option ==55:
        y = dfx - dfx_old
        s = x - x_old

        if h_old is None : 
            h_old = sum(y*y) / sum(s*y)

        h = h_old + y*y/sum(y*s) - (h_old*s)**2/sum(s * (h_old*s))
        rho = sum(y * (h*y)) / sum(y * s)


        M = np.diag(1/ (rho * h))

        eta = 0.5
        beta = sum(dfx*np.diag(M)*y)/sum(y*r_old) - sum(y*np.diag(M)*y) / sum(s*y) * sum(dfx*s)/sum(y*r_old)
        beta_plus = max( beta, eta * sum(dfx*r_old) / sum(r_old*r_old))

        r = - np.diag(M) * dfx + beta_plus * r_old

    return r, h, cond_dfx


def lbfgs(Y,S,RHO,dfx,dfx_old,x,x_old,m=4):
    q = dfx
    y = dfx - dfx_old # df(x_k) - df(x_{k-1)} = y_{k-1}
    s = x - x_old # x_k - x_{k-1} =  s_{k-1}
    rho = 1 / sum(y * s)

    Y = Y[1:] + [y]
    S = S[1:] + [s]
    RHO = RHO[1:] + [rho]
    
    alp = np.zeros(m)
    for i in range(m):
        alp[m-1-i] = RHO[m-1-i] * sum(S[m-1-i]*q)
        q = q - (alp[m-1-i] * Y[m-1-i])
    
    gamma = sum(S[-1]*Y[-1]) / sum(Y[-1]*Y[-1])
    H0 = gamma * np.identity(len(y))
    r = np.dot(H0,q)
    for i in range(m):
        beta = RHO[i]*sum(Y[i]*r)
        r = r + S[i] *(alp[i] - beta)

    return -r

def get_alpha_max(x,r, parametrisation):
    if parametrisation == ["density", "wind", "pressure"] or parametrisation == ["density", "wind", "velocity"]:
        size = int((len(x)+1)/3)
        index_neg = r[:2*size] < -1e-10 
        if len(index_neg) == 0 :
            alpha_max = 20
        else :  
            alpha_max = np.min(-x[:2*size][index_neg]/r[:2*size][index_neg])

        print(alpha_max)
        return alpha_max
    else : 
        return 1


def optimisation(x0, f0, gradf0, hessf0, c1, c2, argf=[],argg=[], argh=[], alpha_max=20, maxiter=1000, tol_x=1e-8, tol_df=1e-8, true_hessian=False, regularisation=0, plot_iter=None, save_info=None,type_gradient=2,type_regul=0, target=None):
    i_restart = 0 
    if regularisation :

        gamma, dz, size,x_init = regularisation
        arg_f = [dz, size , x_init, *argf]#[x1, *argf]
        arg_g = [dz, size, x_init, *argg]#[x1, *argg]
        
        if type_regul == 0:
            def f(x,argf): 
                size = int(argg[1])
                res = f0(x, *argf[3:]) 

                for i in range(2):
                    norm = sum(argf[2][i*size:(i+1)*size]**2)
                    if (norm < 1e-6):
                        norm = 1
                    res += gamma * sum((x[i*size:(i+1)*size]-argf[2][i*size:(i+1)*size])**2) / norm

                norm = sum(argf[2][2*size:]**2)
                if (norm <1e-6):
                    norm = 1
                res += gamma * sum((x[2*size:]-argf[2][2*size:])**2) / norm

                return res
            

            def f2(x,argf): 
                size = int(argg[1])
                cost_function = f0(x, *argf[3:]) 

                reg = 0
                for i in range(2):
                    norm = sum(argf[2][i*size:(i+1)*size]**2)
                    if (norm < 1e-6): 
                        norm = 1
                    reg += gamma * sum((x[i*size:(i+1)*size]-argf[2][i*size:(i+1)*size])**2) / norm

                norm = sum(argf[2][2*size:]**2)
                if (norm <1e-6):
                    norm = 1
                reg += gamma * sum((x[2*size:]-argf[2][2*size:])**2) / norm

                return [cost_function, reg]
            
            def gradf(x,argg):
                grad = gradf0(x, *argg[3:])
                size = int(argg[1])
                for i in range(2):
                    norm = sum(argg[2][i*size:(i+1)*size]**2)
                    if (norm < 1e-6) :
                        norm= 1
                    grad[0][i*size:(i+1)*size] = grad[0][i*size:(i+1)*size] + 2 * gamma * (x[i*size:(i+1)*size] - argg[2][i*size:(i+1)*size]) / norm
                
                norm = sum(argg[2][2*size:]**2)
                if (norm <1e-6):
                    norm = 1
                grad[0][2*size:] = grad[0][2*size:] + 2 * gamma * (x[2*size:] - argg[2][2*size:]) / norm

                return grad
            
        elif type_regul == 1:
            f = lambda x,argf : f0(x, *argf[3:]) + gamma * apply_gradient_regularisation(x, argf[0], argf[1]) 
            f2 = lambda x,argf : [f0(x, *argf[3:]), apply_gradient_regularisation(x, argf[0], argf[1])] 

            def gradf(x,argg):
                grad = gradf0(x, *argg[3:])
                grad[0][:] += gamma * get_derivative_gradient_regularisation(x, argg[0], argg[1]) 
                return grad

        elif type_regul == 2:
            f = lambda x,argf : f0(x, *argf[3:]) + gamma * apply_laplacian_regularisation(x, argf[0], argf[1]) #npl.norm(x - argf[0])**2 * 0.5
            f2 = lambda x,argf : f0(x, *argf[3:]) , apply_laplacian_regularisation(x, argf[0], argf[1]) #npl.norm(x - argf[0])**2 * 0.5
        
            def gradf(x,argg):
                grad = gradf0(x, *argg[3:])
                grad[0][:] += gamma * get_derivative_laplacian_regularisation(x, argg[0], argg[1]) #(x - argf[0])
                return grad

        #f = lambda x,argf : f0(x, *argf[3:]) + gamma[0] * sum((x[:argf[1]]-argf[2][:argf[1]])**2) \
        #                   + gamma[1] * sum((x[argf[1]:2*argf[1]]-argf[2][argf[1]:2*argf[1]])**2) \
        #                   + gamma[0] * sum((x[2*argf[1]:]-argf[2][2*argf[1]:])**2
        
        # def gradf(x,argg):
        #     grad = gradf0(x, *argg[3:])
        #     grad[0][:argg[1]] = grad[0][:argg[1]] + 2 * gamma[0] * (x[:argg[1]] - argg[2][:argg[1]])
        #     grad[0][argg[1]:2*argg[1]] = grad[0][argg[1]:2*argg[1]] + 2 * gamma[1] * (x[argg[1]:2*argg[1]] - argg[2][argg[1]:2*argg[1]])
        #     grad[0][2*argg[1]:] = grad[0][2*argg[1]:] + 2 * gamma[2] * (x[2*argg[1]:] - argg[2][2*argg[1]:])
        #     return grad
        hessf = None #lambda x, preconditioner, argh : hessf0(x, preconditioner, *argh) #+ gamma * np.eye(len(x0))
    else :
        arg_f = argf
        arg_g = argg
        f = lambda x, argf : f0(x,*argf)
        gradf = lambda x, argg : gradf0(x,*argg)
        hessf = None #lambda  x, preconditioner, argh : hessf0(x, preconditioner, *argh)

    x = deepcopy(x0)
    x_iter = []
    fx = f(x, arg_f)
    dfx, preconditioner = gradf(x,arg_g)

    x_old = 1e6 * x0
    fx_old = None
    dfx_old = None
    h_old = None

    iter = 0
    is_notfound = True

    S = []
    Y = []
    RHO = []
    while iter < maxiter and is_notfound:

        # make steepest descent
        if len(S) < 4 :
            if true_hessian :
                H = hessf(x, dfx, x_old, dfx_old)
                cond_dfx = np.dot(npl.inv(H), dfx)
            else :
                #H = hessf(x, dfx, argh, x_old, dfx_old)
                H = np.ones(len(x)) # hessf(x, preconditioner, argh)
                index_near_0 = (H < 1e-16)
                cond_dfx = dfx # * (1 + 1e-16) / (H + 1e-16) # TODO
                cond_dfx[index_near_0] = dfx[index_near_0]

            r = - cond_dfx
            dfx_r = sum(dfx * r)

            if iter > 0 : 
                S += [ x - x_old]
                Y += [ dfx - dfx_old]
                rho = 1 / sum(Y[-1] * S[-1])
                RHO += [rho]

            if iter == 0 :
                alpha_start = 1
            else :
                alpha_start = min(1,1.01*2*(fx-fx_old)/sum(dfx * r))
            new_alpha, x_new, fx_new =  backtracking(alpha_start, x, r, fx, dfx_r, f, arg_f, c1)#linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            dfx_new, preconditioner = gradf(x_new, arg_g)

        else :
            if (( sum((dfx-dfx_old)*r_old) < 1e-8 and type_gradient in [3,4,5]) or npl.norm(dfx-dfx_old) < 1e-15) :
                # restarts
                i_restart += 1 
                config.is_restarted = 1
                r = - dfx
                h_old = None
                cond_dfx = r
                print( "SUIs passe par la")
            else :
                if type_gradient == 6:
                    r = lbfgs(Y,S,RHO,dfx,dfx_old,x,x_old)
                else : 
                    r, h_old, cond_dfx = get_direction_descent(type_gradient, x, fx,dfx,x_old,fx_old,dfx_old,cond_dfx_old,r_old, iter, h_old=h_old, true_hessian=true_hessian)

            aux = update(x,1,r)
            dfx_r = sum(dfx * r) 
          
            if False :# (np.any(aux[:2*int(len(aux)/3)] < 0 )) or dfx_r > 0 :
              print('r too big for alpha = 1')
              r = -dfx
              h_old = None
              cond_dfx = r
              i_restart += 1 
              config.is_restarted = 1

              dfx_r = sum(dfx * r)

            alpha_max = get_alpha_max(x,r, argf[-1])
            alpha_start = min(alpha_max,1,1.01*2*(fx-fx_old)/sum(dfx * r))
            if alpha_start <=0:
                alpha_start = min(1,alpha_max)
         
            file2 = open(path+"/debug_info.txt","a")
            file2.write("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %1d"%(npl.norm(x-x_old), npl.norm(dfx-dfx_old), fx, fx_old, alpha_start, sum((dfx-dfx_old)*r_old), sum(r_old*r_old), sum((dfx-dfx_old)*(dfx-dfx_old)), config.is_restarted))
            file2.close()


            new_alpha, x_new, fx_new, dfx_new =  linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)



        # update
        iter += 1

        x_old = x
        fx_old = fx
        dfx_old = dfx
        cond_dfx_old = cond_dfx
        r_old = deepcopy(r)
        H_old = H

        x = x_new
        x_iter += [x]
        fx = fx_new
        dfx = dfx_new

        

        print("Current iter:", iter)
        if (fx - fx_old) > 0 :
            print("ERREUR augmente la fonction cout")
            #sys.exit()

        if (iter % 1) <1e-3 and not(plot_iter is None):
            funt_plot, arg_plot = plot_iter
            funt_plot(x, *arg_plot, "Iter n° = "+str(iter))
            #plt.show()

        if (iter % 1) <1e-3 and not(save_info is None):
            path, plot_m, plot_gradient, arg_plot = save_info

            plot_m(x, *arg_plot, "Iter n° = "+str(iter))
            print(iter)
            plt.savefig(path+"/mk_iter_"+str(iter))
            plt.close()

            plot_gradient(dfx_old, *arg_plot[:-2], "Iter n° = "+str(iter))
            plt.savefig(path+"/grad_iter_"+str(iter))
            plt.close()

            misfit, cost_regul = f2(x, arg_f)
            _,_,cc,vv = model2allparam(x, *arg_plot[1:],size)
            actual = cc+interpolation(vv)
            diff_model = np.sum((actual- target)**2) / np.sum((target)**2)
            file = open(path+"/iterations_informations.txt", "a")
            file.write("%10d %4.6e %10d %4.6e %4.6e %4.6e %4.6e %10d %10d %10d %4.6e\n"%(iter,new_alpha, config.n_iter_alpha, fx, misfit, cost_regul, npl.norm(dfx), config.n_fx, config.n_grad, config.is_restarted,diff_model))
            file.close()
            config.n_iter_alpha = 0
            config.n_fx = 0
            config.n_grad = 0
            config.n_h = 0

            size = int((len(x)+ 1)/3)
            _,_,cc,vv = model2allparam(x, *arg_plot[1:],size)
            plt.figure()
            plt.plot(arg_plot[0]/1000, cc+interpolation(vv), label="Current model")
            if target is not None :
                plt.plot(arg_plot[0]/1000, target, label="Target")
            #plt.xlim(5,45)
            plt.grid()
            plt.xlabel("Range (km)")
            plt.ylabel("Effective Celerity (m/s)")
            h = arg_plot[0][2] - arg_plot[0][1] 
            plt.axvspan(0,50 * h /1000, facecolor='k', alpha=0.15)
            plt.axvspan((len(arg_plot[0]) - 50) * h /1000, arg_plot[0][-1] * h /1000, facecolor='k', alpha=0.15)
            plt.xlim(arg_plot[0][0]/1000,arg_plot[0][-1]/1000)
            plt.savefig((path+"/celerity_iter_"+str(iter)))
            plt.close()

        is_notfound = True #npl.norm(x-x_old) >= npl.norm(x_old) * tol_x / 1000000000000000 and npl.norm(dfx) >= tol_df #np.any(abs(dfx) >= tol_df)
        if new_alpha < 1e-18 or np.isnan(fx_new) : 
            is_notfound = False
        print(np.sqrt(npl.norm(dfx)**2 * argf[4]))
        config.is_restarted = 0

    is_notfound = [is_notfound, fx_new, npl.norm(dfx_new), iter]
    return x, is_notfound, x_iter
