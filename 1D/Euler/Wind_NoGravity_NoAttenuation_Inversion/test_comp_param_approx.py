

from readline import set_startup_hook
import numpy as np
import matplotlib.pyplot as plt


config0 = "_Validés/BackUps_ctimes1p2_approx_drho_delta0p0001/Calcul analytique" #"_good_rho0times1p2_deltac_d0p0001"# "_good_rho0times1p2_deltac_d0p0001"#_rho0times1p2"
config1 = "_Validés/BackUps_ctimes1p2_approx_drho_delta0p0001/Approximation"# "_good_rho0times1p2_deltac_d0p0001"# "_good_rho0times1p2_approx_deltarho_d0p0001"#_rho0times1p2_approx_deltarho_d0p0001"
config2 = "_Validés/BackUps_ctimes1p2_approx_dc_delta0p0001/Approximation"#"_good_rho0times1p2_deltac_d0p0001"#_ctimes1p2_approx_rho"# "_good_rho0times1p2_deltac_d0p0001"#_rho0times1p2_approx_deltac_d0p0001"
dir_parametrisation = "AcousticEquation_ParamC_Cleaned/"
local_path = "/home/deos/s.gerier/PROJECTS/SIMULATIONS/DF_1D/Inversion_LNS/AcousticEquation/"

plot_several_graph = False

dt = 0.001
Tmax = 70
z0 = -25e3
h = 100 
zmax = 45e3
source = 30000.0
receivers = "[15000.]"

K = np.load(local_path+dir_parametrisation+"BackUps"+config0+"/kernel_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(source)+"_"+receivers+".npy", allow_pickle=True)
K_approach_rho = np.load(local_path+dir_parametrisation+"BackUps"+config1+"/kernel_p_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(source)+"_"+receivers+".npy", allow_pickle=True) 
K_approach_c = np.load(local_path+dir_parametrisation+"BackUps"+config2+"/kernel_p_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(source)+"_"+receivers+".npy", allow_pickle=True) 
M = np.load(local_path+dir_parametrisation+"BackUps"+config2+"/true_model.npy", allow_pickle=True)

nb_index_neg = int(- z0 / h)
z = np.arange(z0,zmax+h,h)


if zmax == 45e3 : 
    start_k = -70+nb_index_neg
    end_k = 230+nb_index_neg

print("Source at : ", source)
print("Receivers at : ", receivers)


rho0  = 1.04898036 
c     = 344.108887

# make a clear array with position of receivers
if receivers == "" : 
    split_name = config0.split("_")
    start_recpt = int(split_name[3])
    end_recpt = int(split_name[4])
    dr = int(split_name[5][1:])
    index_recpt = np.linspace(start_recpt,end_recpt,dr)
else :
    index_recpt = [float(item) / 1000 for item in receivers[1:-1].split()]


# compute the absolute difference
abs_diff = np.zeros((2,len(K_approach_rho)-1))
abs_diff[0,:] = K[0,start_k:end_k] - K_approach_rho[:-1]
abs_diff[1,:] = K[1,start_k:end_k] - K_approach_c[:-1]

# compute the relative difference
relativ_diff = np.zeros((2,len(K_approach_rho)-1))
for i in range(len(relativ_diff[0,:])):
    if K_approach_rho[i] > 1e-8 : 
        relativ_diff[0,i] = abs_diff[0,i] / K_approach_rho[i]
    else : 
        relativ_diff[0,i] = abs_diff[0,i] / max(K_approach_rho)
    if K_approach_c[i] > 1e-8 : 
        relativ_diff[1,i] = abs_diff[1,i] / K_approach_c[i]
    else : 
        relativ_diff[1,i] = abs_diff[1,i] / max(K_approach_c)


if plot_several_graph : 
    plt.figure()
    plt.plot(z[start_k:end_k]/1000, K[0,start_k:end_k], label="from adjoint")
    plt.plot(z[start_k:end_k]/1000, K_approach_rho[:-1], label="from approximation")
    plt.xlabel("Range (km)")
    plt.ylabel("Kernels (rho)")
    plt.legend()
    plt.grid()
    #plt.show()

    plt.figure()
    plt.plot(z[start_k:end_k]/1000, abs_diff[0,:], 'g', label="abs difference")
    plt.xlabel("Range (km)")
    plt.ylabel("delta Kernels (rho)")
    plt.legend()
    #plt.show()

    plt.figure()
    plt.plot(z[start_k:end_k]/1000, relativ_diff[0,:], 'k', label="relative difference")
    plt.xlabel("Range (km)")
    plt.ylabel("delta Kernels (rho)")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(z[start_k:end_k]/1000, K[1,start_k:end_k], label="from adjoint")
    plt.plot(z[start_k:end_k]/1000, K_approach_c[:-1], label="from approximation")
    plt.xlabel("Range (km)")
    plt.ylabel("Kernels (c)")
    plt.legend()
    plt.grid()
    #plt.show()

    plt.figure()
    plt.plot(z[start_k:end_k]/1000, abs_diff[1,:], 'g', label="absolute difference")
    plt.xlabel("Range (km)")
    plt.ylabel("delta Kernels (c)")
    plt.legend()
    plt.grid()
    #plt.show()

    plt.figure()
    plt.plot(z[start_k:end_k]/1000, relativ_diff[1,:], 'k', label="relative difference")
    plt.xlabel("Range (km)")
    plt.ylabel("delta Kernels (c)")
    plt.legend()
    plt.grid()
    plt.show()


else : 
    fig, ax = plt.subplots(2, figsize=(10,8))
    ax[0].plot(z[start_k:end_k]/1000,M[0].rho[start_k:end_k], label="True model")
    ax[0].set_xlabel("Range (km)")
    ax[0].set_ylabel("Density (kg/m)")
    ax[0].legend()
    ax[0].grid()
    ax[0].plot(source / 1000,0, 'or', markersize=10)
    ax[0].plot(index_recpt, np.zeros(len(index_recpt)), '^k')
    ax[1].plot(z[start_k:end_k]/1000, K[0,start_k:end_k], label="from adjoint")
    ax[1].plot(z[start_k:end_k]/1000, K_approach_rho[:-1], label="from approximation")
    ax[1].plot(z[start_k:end_k]/1000, abs_diff[0,:] * 100, '--', label="absolute difference")
    ax[1].set_xlabel("Range (km)")
    ax[1].set_ylabel("Kernels (rho)")
    ax[1].legend()
    ax[1].grid()
    #plt.show()

    fig, ax = plt.subplots(2, figsize=(10,8))
    ax[0].plot(z[start_k:end_k]/1000,M[0].c[start_k:end_k], label="True model")
    ax[0].set_xlabel("Range (km)")
    ax[0].set_ylabel("Wave velocity (kg/m)")
    ax[0].legend()
    ax[0].grid()
    ax[0].plot(source / 1000,0, 'or', markersize=10)
    ax[0].plot(index_recpt, np.zeros(len(index_recpt)), '^k')
    ax[1].plot(z[start_k:end_k]/1000, K[1,start_k:end_k], label="from adjoint")
    ax[1].plot(z[start_k:end_k]/1000, K_approach_c[:-1], label="from approximation")
    ax[1].plot(z[start_k:end_k]/1000, abs_diff[1,:] * 100, '--', label="absolute difference")
    ax[1].set_xlabel("Range (km)")
    ax[1].set_ylabel("Kernels (rho)")
    ax[1].legend()
    ax[1].grid()
    plt.show()