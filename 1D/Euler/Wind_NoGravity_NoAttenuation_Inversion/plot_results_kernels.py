

import numpy as np
import matplotlib.pyplot as plt

from parameters import z, z0, zmax, Tmax, dt, index_source, index_receivers, M, GT_M, h
from compute_approached_kernel_param import start_k, end_k, saving_dir


directory = "./BackUps/"

file1 = open(directory+'parameters.txt', 'r')
Lines = file1.readlines()

# Strips the newline character
for l in range(len(Lines)):
    if "Domaine" in Lines[l]:
        param_line = Lines[l+1].split(',')
        z0 = int(param_line[0])
        zmax = float(param_line[1])
        h = int(param_line[2])
    elif "Source" in Lines[l]:
        param_line = Lines[l+1].split(',')
        index_source = int(param_line[0])
    elif "Receivers" in Lines[l]:
        param_line = Lines[l+1].split(',')
        index_receivers = np.fromstring(
            param_line[0][1:-1], dtype=int, sep=' ')
    elif "Domain approx" in Lines[l]:
        param_line = Lines[l+1].split(',')
        start_k = int(param_line[0])
        end_k = int(param_line[1])
    elif "Time" in Lines[l]:
        param_line = Lines[l+1].split(',')
        T_init = int(param_line[0])
        Tmax = int(param_line[1])
        dt = float(param_line[2])

z_aux = np.arange(z0, zmax+h, h)
z = np.zeros(2*len(z_aux))
z[:len(z_aux)] = z_aux
z[len(z_aux):] = z_aux + z_aux[-1]


GT_M_rho = np.loadtxt("./"+directory+"/model_density_real.csv")
GT_M_w = np.loadtxt("./"+directory+"/model_wind_real.csv")
GT_M_p = np.loadtxt("./"+directory+"/model_pressure_real.csv")

M_rho = np.loadtxt("./"+directory+"/model_density_apriori.csv")
M_w = np.loadtxt("./"+directory+"/model_wind_apriori.csv")
M_p = np.loadtxt("./"+directory+"/model_pressure_apriori.csv")


Kernels    = np.load("./BackUps/kernel_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_.npy", allow_pickle=True)
K_rho      = np.load("./"+saving_dir+"/kernel_rho_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(start_k)+"_"+str(end_k)+".npy", allow_pickle=True)
K_wind     = np.load("./"+saving_dir+"/kernel_wind_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(start_k)+"_"+str(end_k)+".npy", allow_pickle=True)
K_pressure = np.load("./"+saving_dir+"/kernel_pressure_approx_"+str(z0)+"_"+str(zmax)+"_"+str(h)+"_"+str(Tmax)+"_"+str(dt)+"_"+str(z[index_source])+"_"+str(start_k)+"_"+str(end_k)+".npy", allow_pickle=True)


list_kernel = [K_rho, K_wind, K_pressure]
name_kernel = ["Density", "Wind", "Pressure"]

end_k += 1

# check source and receiver
print("Source at : ", z[index_source])
print("Receivers at : ", z[index_receivers])

# compute the difference
abs_diff = np.zeros((4, end_k-start_k))
for k, kernel in enumerate(list_kernel):
    abs_diff[k, :] = Kernels[k, start_k:end_k] - list_kernel[k]

plt.rcParams['font.size'] = 18

fig, ax = plt.subplots(6, figsize=(12, 16), gridspec_kw={
    'height_ratios': [1, 1, 1, 1, 1, 1]})
ax[0].plot(z[start_k:end_k]/1000, GT_M_rho[start_k:end_k], label="True model")
ax[0].plot(z[start_k:end_k]/1000, M_rho[start_k:end_k], label="Apriori model")
ax[0].set_ylabel("Density\n(kg/m)")
ax[0].legend()
ax[0].grid()
ax[1].plot(z[start_k:end_k]/1000, GT_M_w[start_k:end_k], label="True model")
ax[1].plot(z[start_k:end_k]/1000, M_w[start_k:end_k], label="Apriori model")
ax[1].set_ylabel("Wind\n(m/s)")
ax[1].legend()
ax[1].grid()
ax[2].plot(z[start_k:end_k]/1000, GT_M_p[start_k:end_k], label="True model")
ax[2].plot(z[start_k:end_k]/1000, M_p[start_k:end_k], label="Apriori model")
ax[2].set_ylabel("Pressure\n(Pa)")
ax[2].legend()
ax[2].grid()

ax[0].axvline(x=z[index_source] / 1000, ls='--', color='r')
ax[1].axvline(x=z[index_source] / 1000, ls='--', color='r')
ax[2].axvline(x=z[index_source] / 1000, ls='--', color='r')

for r in index_receivers:
    ax[0].axvline(x=z[r] / 1000, ls='--', color='k')
    ax[1].axvline(x=z[r] / 1000, ls='--', color='k')
    ax[2].axvline(x=z[r] / 1000, ls='--', color='k')

ax[3].plot(z[start_k:end_k]/1000, Kernels[0, start_k:end_k], label="M1")
ax[3].plot(z[start_k:end_k]/1000, list_kernel[0], label="M2")
ax[3].plot(z[start_k:end_k]/1000, abs_diff[0, :] * 10, '--', label=r"10 $\times$ (M1 - M2)")
ax[3].set_ylabel("Kernels\n("+name_kernel[0]+")")
ax[3].legend(loc=4)
ax[3].grid()

ax[4].plot(z[start_k:end_k]/1000, Kernels[1, start_k:end_k], label="M1")
ax[4].plot(z[start_k:end_k]/1000, list_kernel[1], label="M2")
ax[4].plot(z[start_k:end_k]/1000, abs_diff[1, :] *
            10, '--', label=r"10 $\times$ (M1 - M2)")
ax[4].set_xlabel("Range (km)")
ax[4].set_ylabel("Kernels\n("+name_kernel[1]+")")
#ax[4].legend(loc=4)
ax[4].grid()

ax[5].plot(z[start_k:end_k]/1000, Kernels[2, start_k:end_k], label="M1")
ax[5].plot(z[start_k:end_k]/1000, list_kernel[2], label="M2")
ax[5].plot(z[start_k:end_k]/1000, abs_diff[2, :] *
            10, '--', label=r"10 $\times$ (M1 - M2)")
ax[5].set_xlabel("Range (km)")
ax[5].set_ylabel("Kernels\n("+name_kernel[2]+")")
#ax[5].legend(loc=4)
ax[5].grid()

if "plsr" in directory : 
    name_fig = "pls_r"
elif "1r" in directory:
    name_fig = "1r"
else :
    name_fig = directory.split('_')[-1][:-1]
plt.savefig("./"+directory +
            "/"+name_fig+".jpg")
plt.show()




# # compute the difference
# abs_diff = np.zeros((4,len(K_rho)))
# for k, kernel in enumerate([K_rho, K_wind, K_pressure]):
#     abs_diff[k,:] = Kernels[k,start_k:end_k+1] - kernel

# list_kernel = [K_rho, K_wind, K_pressure, K_gamma]
# name_kernel = ["Density", "Wind", "Pressure", "Gamma"]
# for k, kernel in enumerate([K_rho, K_wind, K_pressure, K_gamma]):
#     fig, ax = plt.subplots(5, figsize=(10,40), gridspec_kw={'height_ratios': [1, 1, 1, 1, 2]})
#     ax[0].plot(z[start_k:end_k]/1000,GT_M.rho[start_k:end_k], label="True model")
#     ax[0].plot(z[start_k:end_k]/1000,M.rho[start_k:end_k], label="Apriori model")
#     ax[0].set_xlabel("Range (km)")
#     ax[0].set_ylabel("Density (kg/m)")
#     ax[0].legend()
#     ax[0].grid()
#     ax[1].plot(z[start_k:end_k]/1000,GT_M.v[start_k:end_k], label="True model")
#     ax[1].plot(z[start_k:end_k]/1000,M.v[start_k:end_k], label="Apriori model")
#     ax[1].set_xlabel("Range (km)")
#     ax[1].set_ylabel("Wind (m/s)")
#     ax[1].legend()
#     ax[1].grid()
#     ax[2].plot(z[start_k:end_k]/1000,GT_M.p[start_k:end_k], label="True model")
#     ax[2].plot(z[start_k:end_k]/1000,M.p[start_k:end_k], label="Apriori model")
#     ax[2].set_xlabel("Range (km)")
#     ax[2].set_ylabel("Pressure (Pa)")
#     ax[2].legend()
#     ax[2].grid()
#     ax[3].plot(z[start_k:end_k]/1000,GT_M.gamma[start_k:end_k], label="True model")
#     ax[3].plot(z[start_k:end_k]/1000,M.gamma[start_k:end_k], label="Apriori model")
#     ax[3].set_xlabel("Range (km)")
#     ax[3].set_ylabel("Gamma (?)")
#     ax[3].legend()
#     ax[3].grid()

#     ax[0].axvline(x=z[index_source] / 1000, ls='--', color='r')  #plot(z[index_source] / 1000,0, 'or', markersize=10)
#     ax[0].axvline(x=z[index_receivers] /1000, ls='--', color='k') #ax[0].plot(z[index_receivers] /1000, np.zeros(len(index_receivers)), '^k')
#     ax[1].axvline(x=z[index_source] / 1000, ls='--', color='r') #ax[1].plot(z[index_source] / 1000,0, 'or', markersize=10)
#     ax[1].axvline(x=z[index_receivers] /1000, ls='--', color='k') #ax[1].plot(z[index_receivers] /1000, np.zeros(len(index_receivers)), '^k')
#     ax[2].axvline(x=z[index_source] / 1000, ls='--', color='r') #ax[2].plot(z[index_source] / 1000,0, 'or', markersize=10)
#     ax[2].axvline(x=z[index_receivers] /1000, ls='--', color='k') #ax[2].plot(z[index_receivers]/1000, np.zeros(len(index_receivers)), '^k')
#     ax[3].axvline(x=z[index_source] / 1000, ls='--', color='r')#ax[3].plot(z[index_source] / 1000,0, 'or', markersize=10)
#     ax[3].axvline(x=z[index_receivers] /1000, ls='--', color='k') #ax[3].plot(z[index_receivers]/1000, np.zeros(len(index_receivers)), '^k')
#     ax[4].plot(z[start_k:end_k+1]/1000, Kernels[k,start_k:end_k+1], label="from adjoint")
#     ax[4].plot(z[start_k:end_k+1]/1000, list_kernel[k], label="from approximation")
#     ax[4].plot(z[start_k:end_k+1]/1000, abs_diff[k,:] * 100, '--', label="absolute difference")
#     ax[4].set_xlabel("Range (km)")
#     ax[4].set_ylabel("Kernels ("+name_kernel[k]+")")
#     ax[4].legend()
#     ax[4].grid()

#     plt.savefig("./BackUps/kernel_analytic_autodiff_comp_"+name_kernel[k]+".jpg")

# plt.show()
# fig, ax = plt.subplots(2, figsize=(10,8))
# ax[0].plot(z[start_k:end_k]/1000,M[0].c[start_k:end_k], label="True model")
# ax[0].set_xlabel("Range (km)")
# ax[0].set_ylabel("Wave velocity (kg/m)")
# ax[0].legend()
# ax[0].grid()
# ax[0].plot(source / 1000,0, 'or', markersize=10)
# ax[0].plot(index_recpt, np.zeros(len(index_recpt)), '^k')
# ax[1].plot(z[start_k:end_k]/1000, K[1,start_k:end_k], label="from adjoint")
# ax[1].plot(z[start_k:end_k]/1000, K_approach_c[:-1], label="from approximation")
# ax[1].plot(z[start_k:end_k]/1000, abs_diff[1,:] * 100, '--', label="absolute difference")
# ax[1].set_xlabel("Range (km)")
# ax[1].set_ylabel("Kernels (rho)")
# ax[1].legend()
# ax[1].grid()
 
