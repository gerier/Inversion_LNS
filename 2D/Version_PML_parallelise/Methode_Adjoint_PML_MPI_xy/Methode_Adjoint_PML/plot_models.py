import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})
#plt.rcParams.update({'figure.autolayout': True})

plt.rc('axes', titlesize=16)        # Controls Axes Title
plt.rc('axes', labelsize=16)        # Controls Axes Labels
plt.rc('xtick', labelsize=16)       # Controls x Tick Labels
plt.rc('ytick', labelsize=16)       # Controls y Tick Labels
plt.rc('legend', fontsize=16)       # Controls Legend Font
plt.rc('figure', titlesize=16)      # Controls Figure Title


rho_true_0 = np.loadtxt("MODELS/rho0_true_000000.txt")
p_true_0 = np.loadtxt("MODELS/p0_true_000000.txt")
windx_true_0 = np.loadtxt("MODELS/windx_true_000000.txt")
windy_true_0 = np.loadtxt("MODELS/windy_true_000000.txt")

rho_true_1 = np.loadtxt("MODELS/rho0_true_000001.txt")
p_true_1 =  np.loadtxt("MODELS/p0_true_000001.txt")
windx_true_1 = np.loadtxt("MODELS/windx_true_000001.txt")
windy_true_1 = np.loadtxt("MODELS/windy_true_000001.txt")

rho_true_2 = np.loadtxt("MODELS/rho0_true_000002.txt")
p_true_2 = np.loadtxt("MODELS/p0_true_000002.txt")
windx_true_2 = np.loadtxt("MODELS/windx_true_000002.txt")
windy_true_2 = np.loadtxt("MODELS/windy_true_000002.txt")

rho_true_3 = np.loadtxt("MODELS/rho0_true_000003.txt")
p_true_3 =  np.loadtxt("MODELS/p0_true_000003.txt")
windx_true_3 = np.loadtxt("MODELS/windx_true_000003.txt")
windy_true_3 = np.loadtxt("MODELS/windy_true_000003.txt")

n,m = np.shape(rho_true_3)
# true
rho_true = np.zeros((2*n,2*m))
p_true = np.zeros((2*n,2*m))
windy_true = np.zeros((2*n,2*m))
windx_true = np.zeros((2*n,2*m))

rho_true[:n,:m] = rho_true_0
rho_true[n:,:m] = rho_true_1
rho_true[:n,m:] = rho_true_2
rho_true[n:,m:] = rho_true_3

p_true[:n,:m] = p_true_0
p_true[n:,:m] = p_true_1
p_true[:n,m:] = p_true_2
p_true[n:,m:] = p_true_3

windx_true[:n,:m] = windx_true_0
windx_true[n:,:m] = windx_true_1
windx_true[:n,m:] = windx_true_2
windx_true[n:,m:] = windx_true_3

windy_true[:n,:m] = windy_true_0
windy_true[n:,:m] = windy_true_1
windy_true[:n,m:] = windy_true_2
windy_true[n:,m:] = windy_true_3




h = 100 / 10
Ymax = len(rho_true[:,0]) / h
Xmax = len(rho_true[0,:]) / h




fig, ax = plt.subplots(2, 4, figsize=(16,4),layout="constrained") 
im00 = ax[0,0].imshow(rho_true, origin="lower", extent=[0,Xmax,0,Ymax], cmap="Purples", vmin=1.0, vmax=1.5)
im01 = ax[0,1].imshow(p_true, origin="lower", extent=[0,Xmax,0,Ymax], cmap="Oranges", vmin=98000, vmax=105000)
im02 = ax[0,2].imshow(windy_true, origin="lower", extent=[0,Xmax,0,Ymax], cmap="summer", vmin=340, vmax=360)
im03 = ax[0,3].imshow(windx_true, origin="lower", extent=[0,Xmax,0,Ymax], cmap="summer", vmin=-2, vmax=22)
ax[0,0].set_xlabel("Range (km)", fontsize=18)
ax[0,0].set_ylabel("Alt. (km)", fontsize=18)
ax[0,1].set_xlabel("Range (km)", fontsize=18)
ax[0,1].set_ylabel("Alt. (km)", fontsize=18)
ax[0,2].set_xlabel("Range (km)", fontsize=18)
ax[0,2].set_ylabel("Alt. (km)", fontsize=18)
ax[0,0].set_title("True Density", fontsize=20)
ax[0,1].set_title("True Pressure", fontsize=20)
ax[0,2].set_title("True\nVertical Wind", fontsize=20)
ax[0,3].set_title("True\nhorizontal\nwind", fontsize=20)

plt.show()

# a priori
rho = np.loadtxt("MODELS/rho0_prior.txt").T
p = np.loadtxt("MODELS/p0_prior.txt").T
c = np.loadtxt("MODELS/c_prior.txt").T
windx = np.loadtxt("MODELS/windx_prior.txt").T

im10 = ax[1,0].imshow(rho, origin="lower", extent=[0,Xmax,0,Ymax], cmap="Purples", vmin=1.0, vmax=1.5)
im11 = ax[1,1].imshow(p, origin="lower", extent=[0,Xmax,0,Ymax], cmap="Oranges", vmin=98000, vmax=105000)
im12 = ax[1,2].imshow(windy, origin="lower", extent=[0,Xmax,0,Ymax], cmap="summer", vmin=340, vmax=360)
im13 = ax[1,3].imshow(windx, origin="lower", extent=[0,Xmax,0,Ymax], cmap="summer", vmin=-2, vmax=22)
ax[1,0].set_xlabel("Range (km)", fontsize=18)
ax[1,0].set_ylabel("Alt. (km)", fontsize=18)
ax[1,1].set_xlabel("Range (km)", fontsize=18)
ax[1,1].set_ylabel("Alt. (km)", fontsize=18)
ax[1,2].set_xlabel("Range (km)", fontsize=18)
ax[1,2].set_ylabel("Alt. (km)", fontsize=18)
ax[1,0].set_title("A priori Density", fontsize=20)
ax[1,1].set_title("A priori Pressure", fontsize=20)
ax[1,2].set_title("A priori\nVertical Wind", fontsize=20)
ax[1,3].set_title("A priori\nhorizontal\nwind", fontsize=20)


cbar0 = fig.colorbar(im10, ax=ax[:,0], shrink=0.6)
cbar1 = fig.colorbar(im11, ax=ax[:,1], shrink=0.6)
cbar1.formatter.set_powerlimits((0, 0))
cbar2 = fig.colorbar(im12, ax=ax[:,2], shrink=0.6)
cbar3 = fig.colorbar(im13, ax=ax[:,3], shrink=0.6)
cbar0.set_label('Density (km/m)', fontsize=20)
cbar1.set_label('Pressure (Pa)', fontsize=20)
cbar2.set_label("Windy (m/s)", fontsize=20)
cbar3.set_label("Wind (m/s)", fontsize=20)
plt.show()


