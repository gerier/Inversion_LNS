import numpy as np
import numpy.linalg as npl
import sys
import matplotlib.pyplot as plt
import config
from copy import deepcopy
from matplotlib import ticker, cm
from scipy import signal
from Lib_NavierStokesEq.discretisation import apply_laplacian_regularisation,get_derivative_laplacian_regularisation, apply_gradient_regularisation, get_derivative_gradient_regularisation


import config

def wolf1_cond(fx, dfx_r, fx_new, alpha, c1):
    return fx_new <= fx + alpha * c1 * dfx_r

def wolf2strong_cond(r, dfx_r, dfx_new, c2):
    return abs(sum(dfx_new * r)) <= c2 * abs(dfx_r)



#############################################################

def lissage(signal_brut,L):
    res = np.copy(signal_brut) # duplication des valeurs
    for i in range (1,len(signal_brut)-1): # toutes les valeurs sauf la première et la dernière
        L_g = min(i,L) # nombre de valeurs disponibles à gauche
        L_d = min(len(signal_brut)-i-1,L) # nombre de valeurs disponibles à droite
        Li=min(L_g,L_d)
        res[i]= np.median(signal_brut[i-Li:i+Li+1]) #np.sum(signal_brut[i-Li:i+Li+1])/(2*Li+1)
    return res

def lissage(signal_brut,L):
    signal_m = np.roll(signal_brut,1)
    signal_mm = np.roll(signal_brut,2)
    signal_p = np.roll(signal_brut,-1)
    signal_pp = np.roll(signal_brut,-2)
    #res = ( 2*signal_m + 4 * signal_brut + 2* signal_p) / 8
    res = (6*signal_mm + 24*signal_m + 36 * signal_brut + 24* signal_p + 6*signal_pp) / 96

    return res


def lissage(signal_brut,L):
    return signal.medfilt(signal_brut,5)


def update(x, alpha, r):
    size = int((len(x)+1)/3)
    x_new = x + alpha * r

    #x_new[0:size] = lissage(x_new[0:size], 4)
    #x_new[size:2*size] = lissage(x_new[size:2*size], 4)
    #x_new[2*size:] = lissage(x_new[2*size:], 4)

    x_new[0:size] = lissage(x_new[0:size], 4)
    x_new[size:2*size] = lissage(x_new[size:2*size], 4)
    x_new[2*size:] = lissage(x_new[2*size:], 4)

    # plt.figure()
    # plt.plot(x[2*size:])
    # plt.plot(x_new[2*size:])
    # plt.title("alpha = "+str(alpha))
    # plt.show()

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
        #sys.exit()

    return alpha, x_new, fx_new

#############################################################

def quadratic(xa, fxa, dfxa, xb, fxb):
    a = fxb - fxa - dfxa * xb
    b = dfxa * xb**2
    alpha = -b / (2 * a)
    return alpha

def cubic(a, fa, fpa, b, fb, c, fc):
    denom = (b-a)**2 * (c-a)**2 * (b-c)
    d1 = np.zeros((2,2))
    d1[0,0] = (c-a)**2
    d1[0,1] = (a-b)**2
    d1[1,1] = (b-a)**3
    d1[1,0] = (a-c)**3
    d2 = np.zeros((2,1))
    d2[0] = fb - fa - fpa *(b-a)
    d2[1] = fc - fa - fpa * (c-a)

    [A, B] = np.dot(d1, d2)
    A = A / denom
    B = B / denom
    radical = B * B - 3 * A * fpa
    alpha = a + (-B + np.sqrt(radical)) / (3 * A)
    return alpha



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

    array_alpha = []

    if not wolf1_cond(fx, dfx_r, fx_low, alpha_low, c1):
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
            alpha = cubic(x_low, fx_low, sum(dfx_low*r), x_high, fx_high, x_old, fx_old)
  
        # check if the new alpha is in the interval defined by alpha_low, alpha_high
        if (alpha < alpha_low or alpha >= alpha_high) :
            print(alpha, alpha_low, alpha_high)
            alpha = (alpha_high + alpha_low) / 2

        #x_new = x + alpha * r
        x_new = update(x,alpha,r)
        fx_new = f(x_new, arg_f)
        dfx_new, preconditioner = gradf(x_new, arg_g)
        array_alpha += [[alpha, x_new, fx_new, dfx_new]]

        if  wolf1_cond(fx, dfx_r, fx_new, alpha, c1):
            # add the second part of the condition because can have an intervalle where function only increase but with always a too strong condition
            if wolf2strong_cond(r, dfx_r, dfx_new, c2) or alpha_low == alpha:
                stop = True
            else :
                # need to understand why
                if sum(dfx_new * r) * (alpha_high - alpha_low)  >= 0 :
                    alpha_high = alpha_low
                    x_high = x_low
                    fx_high = fx_low
                    dfx_high = dfx_low

                alpha_low = alpha
                x_low = x_new
                fx_low = fx_new
                dfx_low = dfx_new

                if abs(alpha_high - alpha_low) < 1e-15  and sum(dfx_new * r) > 0 :
                    return  best_alpha(array_alpha)

                elif abs(alpha_high - alpha_low) < 1e-15  and sum(dfx_new * r) < 0 :
                    return  best_alpha(array_alpha)

                if alpha_high < alpha_low :
                    aux = alpha_low
                    alpha_low = alpha_high
                    alpha_high =aux
                    aux = fx_low
                    fx_low = fx_high
                    fx_high = aux
                    aux = dfx_low
                    dfx_low = dfx_high
                    dfx_high = aux

                if alpha_low < 1e-16 :
                    return  best_alpha(array_alpha)

        else :
            print("je ne fais que modifier alpha high")
            alpha_high = alpha
            x_high = x_new
            fx_high = fx_new
            dfx_high = dfx_new

            if abs(alpha_low -alpha_high) < 1e-15:
                return  best_alpha(array_alpha)

        x_old = x_new
        fx_old = fx_new


    return alpha, x_new, fx_new, dfx_new



##############################################################

def linesearch(alpha_start, alpha_max, x, r ,fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2):
    alpha_prec = 1e-15
    alpha = alpha_start

    stop = False
    while not stop :
        config.n_iter_alpha += 1

        #x_new  = x + alpha * r
        x_new = update(x, alpha, r)
        fx_new = f(x_new, arg_f)

        if not  wolf1_cond(fx, dfx_r, fx_new, alpha, c1):
            if alpha_prec < alpha :
                alpha_low = alpha_prec
                alpha_high = alpha
            else :
                alpha_high = alpha_prec
                alpha_low = alpha

            alpha, x_new, fx_new, dfx_new = zoom(alpha_low, alpha_high, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            stop = True
        else :
            dfx_new, preconditioner = gradf(x_new, arg_g)
            if wolf2strong_cond(r, dfx_r, dfx_new, c2):
                stop = True
            else :
                if sum(dfx_new * r) >= 0 :
                    if alpha_prec < alpha :
                        alpha_low = alpha_prec
                        alpha_high = alpha
                    else :
                        alpha_high = alpha_prec
                        alpha_low = alpha

                    alpha, x_new, fx_new, dfx_new = zoom(alpha_low, alpha_high, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
                    stop = True
                else :
                    # equivalent à une methode de dichotomie 
                    alpha_prec = alpha
                    alpha = (alpha + alpha_max) / 2
                    print(alpha, alpha_prec, alpha_max)
                    if alpha == alpha_max : 
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
    if option == 1 :
        if true_hessian :
            H = hessf(x, dfx, x_old, dfx_old)
            cond_dfx = np.dot(npl.inv(H), dfx)
        else :
            #H = hessf(x, dfx, argh, x_old, dfx_old)
            H = np.ones(len(x)) #H = hessf(x, preconditioner, argh)
            #H = hessian_approximation(H_old, x_old, x, dfx_old, dfx)
            index_near_0 = (H < 1e-16)
            cond_dfx = (dfx + 1e-16) / (H + 1e-16)
            cond_dfx[index_near_0] = dfx[index_near_0]
            
        if (iter % len(x) == 0 ) : #abs(sum(dfx * cond_dfx_old)) / sum(dfx * cond_dfx) >= 0.1:
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

def optimisation(x0, f0, gradf0, hessf0, c1, c2, argf=[],argg=[], argh=[], alpha_max=20, maxiter=1000, tol_x=1e-8, tol_df=1e-8, true_hessian=False, regularisation=0, plot_iter=None, save_info=None,type_gradient=2,type_regul=0):
    i_restart = 0 
    if regularisation :

        gamma, dz, size,x_init = regularisation
        arg_f = [dz, size , x_init, *argf]#[x1, *argf]
        arg_g = [dz, size, x_init, *argg]#[x1, *argg]
        
        if type_regul == 0:
            def f(x,argf):
                index_source = arg_f[10][2]
                factor_dist_source = 1 + 10 * np.exp(- (np.arange(0,arg_f[1])- index_source)**2/2 /10)
                factor_dist_source_demi = 1 + 10 * np.exp(- (np.arange(0,arg_f[1]-1)+0.5- index_source)**2/2 /10)
                factor_dist_source = np.concatenate([factor_dist_source,factor_dist_source,factor_dist_source_demi])
                return f0(x, *argf[3:]) + gamma * sum(factor_dist_source * (x-argf[2])**2) 
            
            def f2(x,argf):
                index_source = arg_f[10][2]
                factor_dist_source = 1 + 10 * np.exp(- (np.arange(0,arg_f[1])- index_source)**2/2 /10)
                factor_dist_source_demi = 1 + 10 * np.exp(- (np.arange(0,arg_f[1]-1)+0.5- index_source)**2/2 /10)
                factor_dist_source = np.concatenate([factor_dist_source,factor_dist_source,factor_dist_source_demi])
                return [f0(x, *argf[3:]), gamma * sum(factor_dist_source * (x-argf[2])**2)]

            def gradf(x,argg):
                index_source = arg_g[8][2]
                factor_dist_source = 1 + 10 * np.exp(- (np.arange(0,arg_f[1])- index_source)**2/2 /10)
                factor_dist_source_demi = 1 + 10 * np.exp(- (np.arange(0,arg_f[1]-1)+0.5- index_source)**2/2 /10)

                grad = gradf0(x, *argg[3:])
                
                plt.figure()
                plt.plot( grad[0][argg[1]:2*argg[1]]  + 2 * gamma * (x[argg[1]:2*argg[1]] - argg[2][argg[1]:2*argg[1]]))

                grad[0][:argg[1]] = grad[0][:argg[1]] + 2 * gamma * (x[:argg[1]] - argg[2][:argg[1]]) * factor_dist_source
                grad[0][argg[1]:2*argg[1]] = grad[0][argg[1]:2*argg[1]] + 2 * gamma * (x[argg[1]:2*argg[1]] - argg[2][argg[1]:2*argg[1]]) * factor_dist_source
                grad[0][2*argg[1]:] = grad[0][2*argg[1]:] + 2 * gamma * (x[2*argg[1]:] - argg[2][2*argg[1]:]) * factor_dist_source_demi

                #plt.plot(factor_dist_source)
                plt.plot(grad[0][argg[1]:2*argg[1]] )
                plt.show()

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
        hessf = lambda x, preconditioner, argh : hessf0(x, preconditioner, *argh) #+ gamma * np.eye(len(x0))
    else :
        arg_f = argf
        arg_g = argg
        f = lambda x, argf : f0(x,*argf)
        gradf = lambda x, argg : gradf0(x,*argg)
        hessf = lambda  x, preconditioner, argh : hessf0(x, preconditioner, *argh)

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
    while iter < maxiter and is_notfound:

        # make steepest descent
        if iter <= 2 :
            if true_hessian :
                H = hessf(x, dfx, x_old, dfx_old)
                cond_dfx = np.dot(npl.inv(H), dfx)
            else :
                #H = hessf(x, dfx, argh, x_old, dfx_old)
                H = np.ones(len(x)) # hessf(x, preconditioner, argh)
                index_near_0 = (H < 1e-16)
                cond_dfx = dfx * (1 + 1e-16) / (H + 1e-16) # TODO
                cond_dfx[index_near_0] = dfx[index_near_0]

            r = - cond_dfx
            dfx_r = sum(dfx * r)

            if iter == 0 :
                alpha_start = 1
            else :
                alpha_start = min(1,1.01*2*(fx-fx_old)/sum(dfx * r))
            #new_alpha, x_new, fx_new, dfx_new = linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            new_alpha, x_new, fx_new =  backtracking(alpha_start, x, r, fx, dfx_r, f, arg_f, c1)#linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            dfx_new, preconditioner = gradf(x_new, arg_g)

        else :
            if (( sum((dfx-dfx_old)*r_old) < 1e-8 and type_gradient in [4,5,6]) or npl.norm(dfx-dfx_old) < 1e-15) :
                # restarts
                i_restart += 1 
                r = - dfx
                h_old = None
                cond_dfx = r
            else :
                r, h_old, cond_dfx = get_direction_descent(type_gradient, x, fx,dfx,x_old,fx_old,dfx_old,cond_dfx_old,r_old, iter, h_old=h_old, true_hessian=true_hessian)

            aux = update(x,1,r) 
            if (np.any(aux[:2*int(len(aux)/3)] < 0 )):
              print('r too big for alpha = 1')
              r = -dfx
              h_old = None
              cond_dfx = r

            dfx_r = sum(dfx * r)

            alpha_start = min(1,1.01*2*(fx-fx_old)/sum(dfx * r))
            if alpha_start <=0:
                alpha_start = 1
         
            file2 = open(path+"/debug_info.txt","a")
            file2.write("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e"%(npl.norm(x-x_old), npl.norm(dfx-dfx_old), fx, fx_old, alpha_start, sum((dfx-dfx_old)*r_old), sum(r_old*r_old), sum((dfx-dfx_old)*(dfx-dfx_old))))
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
            path, plot_m, plot_gradient, plot_hessian, arg_plot = save_info

            plot_m(x, *arg_plot, "Iter n° = "+str(iter))
            plt.savefig(path+"/mk_iter_"+str(iter))
            plt.close()

            plot_gradient(dfx_old, *arg_plot[:-2], "Iter n° = "+str(iter))
            plt.savefig(path+"/grad_iter_"+str(iter))
            plt.close()

            plot_hessian(H, *arg_plot[:-2], "Iter n° = "+str(iter))
            plt.savefig(path+"/hessian_iter_"+str(iter))
            plt.close()

            misfit, cost_regul = f2(x, arg_f)
            file = open(path+"/iterations_informations.txt", "a")
            file.write("%10d %4.6e %10d %4.6e %4.6e %4.6e %4.6e %10d %10d %10d\n"%(iter,new_alpha, config.n_iter_alpha, fx, misfit, cost_regul, npl.norm(dfx), config.n_fx, config.n_grad, config.n_h))
            file.close()
            config.n_iter_alpha = 0
            config.n_fx = 0
            config.n_grad = 0
            config.n_h = 0

        is_notfound = True #npl.norm(x-x_old) >= npl.norm(x_old) * tol_x / 1000000000000000 and npl.norm(dfx) >= tol_df #np.any(abs(dfx) >= tol_df)

    print("Restart gradient due to small values : ", i_restart)
    is_notfound = [is_notfound, fx_new, npl.norm(dfx_new), iter]
    return x, is_notfound, x_iter
