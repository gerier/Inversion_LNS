

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline



def get_model(path_file, z, no_gravity, no_wind):
    
    #path_file = "/home/deos/s.gerier/PROJECTS/SIMULATIONS/ATMOSPHERIC_MODELS/MSISE/msisehwm_wrapper/OUTPUT/Flores_atmosphere_500km_86p2_4.dat"
    #path_file = "./Flores_atmosphere_500km_86p2_4.dat"
    data = np.genfromtxt(path_file, skip_header=3)#,delimiter=[7,7,7,7,7,7,7,3,3,3,3,5,5,7,7,7])
    df = pd.DataFrame(data, columns =["z[m]", "rho[kg/(m^3)]", "T[K]", "c[m/s]", "p[Pa]", "H[m]", "g[m/(s^2)]", "N^2[rad^2/s^2]", "kappa[J/(s.m.K)]", "mu[kg(s.m)]", "mu_vol[kg/(s.m)]", "w_M[m/s]", "w_Z[m/s]", "w_P[m/s]", "c_p[J/(mol.K)]", "c_v[J/(mol.K)]", "gamma"])

    physical_parameters = ["rho[kg/(m^3)]", "T[K]", "c[m/s]", "p[Pa]", "g[m/(s^2)]", "kappa[J/(s.m.K)]", "mu[kg(s.m)]", "mu_vol[kg/(s.m)]", "w_P[m/s]", "c_v[J/(mol.K)]", "gamma"]

    interp_physical_param = []
    for param in physical_parameters:
        if param == "rho[kg/(m^3)]" or param == "p[Pa]" : 
            cs = CubicSpline(df["z[m]"], np.log(df[param]))
            interp_physical_param += [np.exp(cs(z))]
        else : 
            cs = CubicSpline(df["z[m]"], df[param])
            interp_physical_param += [cs(z)]

    #n = int(len(z)/2)
    #for param in interp_physical_param :
    #    param[n:] = param[n:0:-1]

    [rho0, T0, c, p0, g, kappa, mu, eta, v0, cv, gamma] = interp_physical_param

    M = 28.965 # masse molaire de l'air https://fr.wikipedia.org/wiki/Air
    l = eta - (2/3)*mu
    Cv = cv/M

    if no_gravity :
        g *= 0
    if no_wind:
        v0 *= 0

    return [rho0, T0, c, p0, g, kappa, mu, eta, v0, cv, gamma, l, Cv]

