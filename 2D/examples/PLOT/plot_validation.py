import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from matplotlib import ticker

import numpy as np
import glob
import os

dir_adj = "OUTPUT_ADJ/"
dir_df = "OUTPUT_DF/"

NX = 200
NY = 150
NPROCX = 4
NPROCY_ADJ = 1
NPROCY_DF = 2
dx = 100

c_prior = 347.763977
rho0_prior = 1.13837624
gamma = 1.4

win_lim = [37.325, 9.05]

N_PML = 25

def load_kernel(dir_output, NX, NY, NPROCX, NPROCY, dx, kernel="windx", N_PML=N_PML):
	list_kernel_path = glob.glob(f'{dir_output}/K{kernel}*')
	print(f'{dir_output}/K{kernel}*')
	LOC_NX = int(NX / NPROCX)
	LOC_NY = int(NY / NPROCY)
	K = np.zeros((NX,NY))
	for fich in list_kernel_path:

		data = np.loadtxt(fich)

		num_proc = int(fich.split("_")[-1].split(".")[0])
		i = int(num_proc% NPROCX)
		j = int(num_proc/ NPROCX)

		K[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data

	K = K *dx*dx

	K[:,NY-N_PML:] = 0
	K[:N_PML,:] = 0
	K[NX-N_PML,:] = 0

	return K
def load_global_field(param, nproc_x, nproc_y,
                      nx_global, ny_global,
                      dtype=np.float32,
                      folder="."):
    nx_local = nx_global // nproc_x
    ny_local = ny_global // nproc_y
    field = np.zeros((ny_global, nx_global), dtype=dtype)
    for rank in range(nproc_x * nproc_y):
        i_rank = rank % nproc_x
        j_rank = rank // nproc_x
        filename = f"{param}_{rank:06d}.bin"
        filepath = os.path.join(folder, filename)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Manque fichier: {filepath}")
        data = np.fromfile(filepath, dtype=dtype)
        data = data.reshape((ny_local, nx_local))
        #data = np.flipud(data)   # <-- correction ici
        x_start = i_rank * nx_local
        y_start = j_rank * ny_local
        field[y_start:y_start+ny_local,
              x_start:x_start+nx_local] = data
    return field


def convol_gauss_kernel(K, sigma, correction):
    n, m = K.shape
    n_points = int(3*sigma)

    x = np.arange(-n_points, n_points+1)
    X, Y = np.meshgrid(x, x)
    G = np.exp(-(X**2 + Y**2)/(2*sigma**2))

    convol = np.zeros_like(K)
    for i in range(n_points, n-n_points):
        for j in range(n_points, m-n_points):
            sub = K[i-n_points:i+n_points+1,
                    j-n_points:j+n_points+1]
            convol[i,j] = np.sum(sub*G)

    return convol


X = np.linspace(0,(NX-1)*dx/1e3,NX)
Y = np.linspace(0,(NY-1)*dx/1e3,NY) - 0.05

data = np.loadtxt(f"{dir_adj}/pressure_file_obs_001.dat")
time = data[:,0]
obs = data[:,1]

data = np.loadtxt(f"{dir_adj}/pressure_file_001.dat")
time = data[:,0]
synth = data[:,1]

K_wx_adj = load_global_field("Kwindx", NPROCX, NPROCY_ADJ,
                      NX, NY,
                      dtype=np.float32,
                      folder=dir_adj)
K_p_adj = load_global_field("Kp0", NPROCX, NPROCY_ADJ,
                      NX, NY,
                      dtype=np.float32,
                      folder=dir_adj)
K_rho_adj = load_global_field("Krho0", NPROCX, NPROCY_ADJ,
                      NX, NY,
                      dtype=np.float32,
                      folder=dir_adj)

K_wx_df = load_kernel(dir_df+"WINDX/", NX, NY, NPROCX, NPROCY_DF, dx, "")
K_p_df = load_kernel(dir_df+"PRESSURE/", NX, NY, NPROCX, NPROCY_DF, dx, "")
K_rho_df = load_kernel(dir_df+"RHO/", NX, NY, NPROCX, NPROCY_DF, dx, "")

K_wx_adj = convol_gauss_kernel(K_wx_adj, sigma = 1.5, correction=3)
K_p_adj = convol_gauss_kernel(K_p_adj, sigma = 1.5, correction=1)
K_rho_adj = convol_gauss_kernel(K_rho_adj, sigma = 1.5, correction=2)


K_c_adj =  2 * c_prior * rho0_prior / gamma * K_p_adj
K_c_df = 2 * c_prior * rho0_prior / gamma * K_p_df

#K_c_adj = convol_gauss_kernel(K_p_adj, sigma = 1.5)


# ------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------
matplotlib.rcParams.update({'font.size': 16})

labelsize = 18

fig = plt.figure(figsize=(12, 12))
gs = GridSpec(
    nrows=3,
    ncols=2,
    height_ratios=[2.0, 3.0, 3.0],
    #hspace=0.12,
    #wspace=0.08
)

# ------------------------------------------------------------------
# Axes
# ------------------------------------------------------------------
ax0 = fig.add_subplot(gs[0, :])   # top panel
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[1, 1])
ax3 = fig.add_subplot(gs[2, 0])
ax4 = fig.add_subplot(gs[2, 1])

# ------------------------------------------------------------------
# Example plots
# Replace these with your own plotting commands
# ------------------------------------------------------------------

# a)
ax0.plot(time, obs, 'k', alpha=0.5, label="Observation")
ax0.plot(time, synth, 'k', linestyle="--", label="Synthetic")
ax0.set_xlabel("Time [s]")
ax0.set_ylabel("Pressure [Pa]")
ax0.legend()

mask = (time > win_lim[0]) & (time < win_lim[0]+win_lim[1])
ax0.axvline(win_lim[0], color="r", linestyle=":", linewidth=2)
ax0.axvline(win_lim[0]+win_lim[1], color="r", linestyle=":", linewidth=2)

# b)
maxval = np.max(np.max(abs(K_c_adj))/5)
print("Maxval Kc = ", np.max(np.max(abs(K_c_adj))))
print("Quantile 99 Kc = ", np.quantile(abs(K_c_adj), 0.99))
print("Norm RMS de Kc = ", np.sqrt(np.mean(K_c_adj**2)))
print("STD de Kc = ", np.std(K_c_adj))
pcm1 = ax1.pcolormesh(X, Y, K_c_adj, shading='auto',
						cmap=cm.bwr,
						vmin=-maxval, vmax=maxval,
                        rasterized=True)
ax1.set_xlabel("Range [km]")
ax1.set_ylabel("Altitude [km]")

# c)
#maxval = np.max(np.max(abs(K_c_df))/5)
pcm2 = ax2.pcolormesh(X, Y, K_c_df.T, shading='auto',
						cmap=cm.bwr,
						vmin=-maxval, vmax=maxval,
                        rasterized=True)
ax2.set_xlabel("Range [km]")
ax2.set_ylabel("Altitude [km]")

# d)
maxval = np.max(np.max(abs(K_wx_adj))/5)
print("Maxval Kwx = ", np.max(np.max(abs(K_wx_adj))))
print("Quantile 99 Kwx = ", np.quantile(abs(K_wx_adj), 0.99))
print("Norm RMS de Kwx = ", np.sqrt(np.mean(K_wx_adj**2)))
print("STD de Kwx = ", np.std(K_wx_adj))

pcm3 = ax3.pcolormesh(X, Y, K_wx_adj, shading='auto',
						cmap=cm.bwr,
						vmin=-maxval, vmax=maxval,
                        rasterized=True)
ax3.set_xlabel("Range [km]")
ax3.set_ylabel("Altitude [km]")

# e)
#maxval = np.max(np.max(abs(K_wx_df))/5)
pcm4 = ax4.pcolormesh(X, Y, K_wx_df.T, shading='auto',
						cmap=cm.bwr,
						vmin=-maxval, vmax=maxval,
                        rasterized=True)
ax4.set_xlabel("Range [km]")
ax4.set_ylabel("Altitude [km]")

# colorbar
cbformat = ticker.ScalarFormatter()
cbformat.set_scientific('%0.2e')
cbformat.set_powerlimits((0,2))
cbformat.set_useMathText(True)

fig.colorbar(pcm1, ax=ax1, label="Sensitivity [$s^3/m$]", format=cbformat)
fig.colorbar(pcm2, ax=ax2, label="Sensitivity [$s^3/m$]", format=cbformat)
fig.colorbar(pcm3, ax=ax3, label="Sensitivity [$s^3/m$]", format=cbformat)
fig.colorbar(pcm4, ax=ax4, label="Sensitivity [$s^3/m$]", format=cbformat)

# ------------------------------------------------------------------
# Panel labels
# ------------------------------------------------------------------
labels = ['a)', 'b)', 'c)', 'd)', 'e)']
axes = [ax0, ax1, ax2, ax3, ax4]

for lab, ax in zip(labels[1:], axes[1:]):
    ax.text(
        0.0, 1.1,
        lab,
        transform=ax.transAxes,
        fontsize=18,
        fontweight='bold',
        va='top',
        ha='left'
    )
ax0.text(
	0.0, 1.15,
	labels[0],
	transform=ax0.transAxes,
	fontsize=18,
	fontweight='bold',
	va='top',
	ha='left'
)

# ------------------------------------------------------------------
# Tick label size
# ------------------------------------------------------------------
for ax in axes:
    ax.tick_params(axis='both')

plt.tight_layout()


fig.savefig("kernel_validation.pdf")
plt.show()
