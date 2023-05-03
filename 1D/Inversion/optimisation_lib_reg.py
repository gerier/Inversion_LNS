import numpy as np
import numpy.linalg as npl
import sys
import matplotlib.pyplot as plt
import config
from copy import deepcopy
from matplotlib import ticker, cm

def wolf1_cond(fx, dfx_r, fx_new, alpha, c1):
    return fx_new <= fx + alpha * c1 * dfx_r

def wolf2strong_cond(r, dfx_r, dfx_new, c2):
    return abs(sum(dfx_new * r)) <= c2 * abs(dfx_r)

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
        x_new = x+alpha*r
        fx_new = f(x_new, arg_f)
        array_alpha += [[alpha, x_new, fx_new]]
        if wolf1_cond(fx, dfx_r, fx_new, alpha, c1):
            stop = True
        else : 
            alpha = rate * alpha
    
        iter +=1

    if iter == maxiter : 
        print("ERROR : No acceptable step find in backtracking")
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

    x_low = x + alpha_low * r
    fx_low = f(x_low, arg_f)
    dfx_low = gradf(x_low, arg_g)
    
    x_high = x + alpha_high * r
    fx_high = f(x_high, arg_f)
    dfx_high = gradf(x_high, arg_g)

    array_alpha = []

    if not wolf1_cond(fx, dfx_r, fx_low, alpha_low, c1):
        print("ERROR : alpha low ne respecte pas 1e condition")
        print("Alpha_low = ", alpha_low)
        alpha, x_new, fx_new = backtracking(1, x, r, fx, dfx_r, f, arg_f, c1)
        dfx_new = gradf(x_new, arg_g)
        stop = True
        # sys.exit()
  
    while not stop:

        # interpolate an alpha (need to be checked)
        if iter == 0 :
            alpha  = quadratic(alpha_low, fx_low, sum(dfx_low*r), alpha_high, fx_high)
        else : 
            alpha = cubic(x_low, fx_low, sum(dfx_low*r), x_high, fx_high, x_old, fx_old)
   
        # check if the new alpha is in the interval defined by alpha_low, alpha_high
        if True or (alpha <= alpha_low or alpha >= alpha_high) :  
            alpha = (alpha_high + alpha_low) / 2

        x_new = x + alpha * r
        fx_new = f(x_new, arg_f)
        dfx_new = gradf(x_new, arg_g)
        array_alpha += [[alpha, x_new, fx_new, dfx_new]]

        if  wolf1_cond(fx, dfx_r, fx_new, alpha, c1):

            if wolf2strong_cond(r, dfx_r, dfx_new, c2):
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

                if abs(alpha_high - alpha_low) < 1e-16  and sum(dfx_new * r) > 0 : 
                    return  best_alpha(array_alpha)
                    
                elif abs(alpha_high - alpha_low) < 1e-16  and sum(dfx_new * r) < 0 : 
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
            alpha_high = alpha
            x_high = x_new
            fx_high = fx_new
            dfx_high = dfx_new 


        x_old = x_new
        fx_old = fx_new


    return alpha, x_new, fx_new, dfx_new



##############################################################

def linesearch(alpha_start, alpha_max, x, r ,fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2):
    alpha_prec = 1e-15
    alpha = alpha_start

    stop = False
    while not stop : 

        x_new  = x + alpha * r
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
            dfx_new = gradf(x_new, arg_g)
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
                    alpha_prec = alpha
                    alpha = (alpha + alpha_max) / 2 

    return alpha, x_new, fx_new, dfx_new

def test_linesearch(alpha_start, alpha_max, x, r, fx, dfx, f, gradf, c1, c2, attendu):
    dfx_r = sum(r * dfx)
    new_alpha, x_new, fx_new, dfx_new = linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)

    print("Attendu : ", attendu)
    print("Obtenu : ", new_alpha)


##################################################################


def optimisation(x0, f0, gradf0, hessf0, c1, c2, argf=[],argg=[], alpha_max=20, maxiter=1000, tol_x=1e-8, tol_df=1e-8, true_hessian=False, regularisation=0, plot_iter=None):
    
    if regularisation : 
        gamma, x1 = regularisation
        arg_f = [x1, *argf]
        arg_g = [x1, *argg]
        f = lambda x,argf : f0(x, *argf[1:]) + gamma * npl.norm(x - argf[0])**2 * 0.5 
        gradf = lambda x,argg : gradf0(x, *argg[1:])[1:] + gamma * (x - argf[0])
        hessf = lambda x, k, x_prec, k_prec : hessf0(x, k, x_prec, k_prec) + gamma * np.eye(2)
    else : 
        arg_f = argf
        arg_g = argg
        f = lambda x, argf : f0(x,*argf)
        gradf = lambda x, argg : gradf0(x,*argg)
        hessf = hessf0
    
    x = deepcopy(x0)
    x_iter = []
    fx = f(x, arg_f)
    dfx = gradf(x,arg_g)


    x_old = 1e6 * x0
    dfx_old = None

    iter = 0
    is_notfound = True
    while iter < maxiter and is_notfound:

        # make steepest descent
        if iter <= 3 :
            if true_hessian : 
                H = hessf(x, dfx, x_old, dfx_old)
                cond_dfx = np.dot(npl.inv(H), dfx)
            else : 
                H = hessf(x, dfx, x_old, dfx_old)
                cond_dfx = (dfx + 1e-40) / (H + 1e-40) 

            r = - cond_dfx
            dfx_r = sum(dfx * r)

            if iter == 0 : 
                alpha_start = 1
            else : 
                alpha_start = min(1,1.01*2*(fx-fx_old)/sum(dfx * r))
            #new_alpha, x_new, fx_new, dfx_new = linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            new_alpha, x_new, fx_new =  backtracking(alpha_start, x, r, fx, dfx_r, f, arg_f, c1)#linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            dfx_new = gradf(x_new, arg_g)

        else :
     
            if true_hessian : 
                H = hessf(x, dfx, x_old, dfx_old)
                cond_dfx = np.dot(npl.inv(H), dfx)
            else : 
                H = hessf(x, dfx, x_old, dfx_old)
                cond_dfx = (dfx + 1e-40) / (H + 1e-40) 

            if (iter % len(x) == 0 ) : #abs(sum(dfx * cond_dfx_old)) / sum(dfx * cond_dfx) >= 0.1:
                print("Gradients are far from orthogonal. Restarting gradients") #iter % 100 == 0  : 
                beta = 0
            else : 
                beta = sum(dfx * cond_dfx) / sum(dfx_old * cond_dfx_old)
            r = - cond_dfx + beta * r_old


            dfx_r = sum(dfx * r)
        
            alpha_start = min(1,1.01*2*(fx-fx_old)/sum(dfx * r))
    
            new_alpha, x_new, fx_new, dfx_new =  linesearch(alpha_start, alpha_max, x, r, fx, dfx_r, f, gradf, arg_f, arg_g, c1, c2)
            #dfx_new = gradf(x_new, arg_g)
     

        # update
        iter += 1
        
        x_old = x
        fx_old = fx
        dfx_old = dfx
        cond_dfx_old = cond_dfx
        r_old = deepcopy(r)

        x = x_new
        x_iter += [x]
        fx = fx_new
        dfx = dfx_new

        if (fx - fx_old) > 0 : 
            print("ERREUR augmente la fonction cout")
            #sys.exit()

        if (iter % 1) <1e-3 and not(plot_iter is None):
            funt_plot, arg_plot = plot_iter
            funt_plot(x, arg_plot, "Iter n° = "+str(iter))

        is_notfound = npl.norm(x-x_old) >= npl.norm(x_old) * tol_x / 1000000000 and npl.norm(dfx) >= tol_df #np.any(abs(dfx) >= tol_df)

    is_notfound = [is_notfound, fx_new, npl.norm(dfx_new), iter]
    return x, is_notfound, x_iter



