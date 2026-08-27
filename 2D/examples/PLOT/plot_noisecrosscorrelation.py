import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from matplotlib import ticker

import numpy as np
import glob

import os

matplotlib.rcParams.update({'font.size': 16})
# =========================
# PARAMETERS
# =========================

model_0p5_path = "model.dat"
model_0p05_path = "model.dat"

dir_0p5_all = "OUTPUT_all/"
dir_0p5_direct = "OUTPUT_direct/"
dir_0p5_refracted = "OUTPUT_refracted/"


dir_0p05 = "./OUTPUT_0p05/"

NX_0p5 = 3740
NY_0p5 = 2250
NPROCX_0p5 = 20
NPROCY_0p5 = 18
dx_0p5 = 20

#NX_0p05 = 750
#NY_0p05 = 450
#NPROCX_0p05 = 2
#NPROCY_0p05 =2
#dx_0p05 = 100
NX_0p05 = 3740
NY_0p05 = 2250
NPROCX_0p05 = 20
NPROCY_0p05 = 18
dx_0p05 = 20

win_lim_0p5_all = [118.2, 20.09]
win_lim_0p5_direct = [118.2, 2.58]
win_lim_0p5_refracted = [134.85, 3.44]
win_lim_0p05 = [106.475, 56.1]


# =========================
# LOAD DATA
# =========================

def load_model(path, NX, NY):
	c = np.zeros(NY)
	wx = np.zeros(NY)

	data = np.loadtxt(path,skiprows=3,delimiter=',')

	z = data[:,0] / 1e3
	c = data[:,3]
	wx = data[:,-4]
	ceff = c + wx

	return z, c, wx, ceff

z_0p5, c_0p5, wx_0p5, ceff_0p5 = load_model(model_0p5_path, NX_0p5, NY_0p5)
z_0p05, c_0p05, wx_0p05, ceff_0p05 = load_model(model_0p05_path, NX_0p05, NY_0p05)

y_model_0p5 = np.linspace(0,NY_0p5*dx_0p5/1e3,NY_0p5) - 0.05
y_model_0p05 = np.linspace(0,NY_0p05*dx_0p05/1e3,NY_0p05) - 0.05

def load_kernel(dir_output, NX, NY, NPROCX, NPROCY, dx, kernel="windx"):
	list_kernel_path = glob.glob(f'{dir_output}/K{kernel}*')
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

X_0p5 = np.linspace(0,(NX_0p5-1)*dx_0p5/1e3,NX_0p5)
Y_0p5 = np.linspace(0,(NY_0p5-1)*dx_0p5/1e3,NY_0p5) - 0.05

X_0p05 = np.linspace(0,(NX_0p05-1)*dx_0p05/1e3,NX_0p05)
Y_0p05 = np.linspace(0,(NY_0p05-1)*dx_0p05/1e3,NY_0p05) - 0.05

data = np.loadtxt(f"{dir_0p5_direct}/pressure_file_001.dat")
t_0p5 = data[:,0]
ts_0p5 = data[:,1]

data = np.loadtxt(f"{dir_0p05}/pressure_file_001.dat")
t_0p05 = data[:,0]
ts_0p05 = data[:,1]


K_0p5_all = load_global_field("Kwindx", NPROCX_0p5, NPROCY_0p5,
                      NX_0p5, NY_0p5,
                      dtype=np.float32,
                      folder=dir_0p5_all)

K_0p5_direct = load_global_field("Kwindx", NPROCX_0p5, NPROCY_0p5,
                      NX_0p5, NY_0p5,
                      dtype=np.float32,
                      folder=dir_0p5_direct)

K_0p5_refracted =  load_global_field("Kwindx", NPROCX_0p5, NPROCY_0p5,
                      NX_0p5, NY_0p5,
                      dtype=np.float32,
                      folder=dir_0p5_refracted)

K_0p05 = load_global_field("Kwindx", NPROCX_0p5, NPROCY_0p5,
                      NX_0p5, NY_0p5,
                      dtype=np.float32,
                      folder=dir_0p05)


# =========================
# FIGURE
# =========================

fig = plt.figure(figsize=(12, 10))
gs = GridSpec(
    nrows=4,
    ncols=6,
    figure=fig,
)

# =========================
# Ligne 1
# =========================
ax_model_x    = fig.add_subplot(gs[0, 0:2])
ax_model_c    = fig.add_subplot(gs[0, 2:4])
ax_model_ceff = fig.add_subplot(gs[0, 4:6])

# plot
ax_model_x.plot(c_0p5, z_0p5, "k")
ax_model_c.plot(wx_0p5, z_0p5, "k")
ax_model_ceff.plot(ceff_0p5, z_0p5, "k")
ax_model_x.plot(c_0p05, z_0p05, "k")
ax_model_c.plot(wx_0p05, z_0p05, "k")
ax_model_ceff.plot(ceff_0p05, z_0p05, "k")
ax_model_x.text(-0.12, 1.2, 'a)', transform=ax_model_x.transAxes,
                fontsize=18, fontweight='bold')
ax_model_c.text(-0.12, 1.2, 'b)', transform=ax_model_c.transAxes,
                fontsize=18, fontweight='bold')
ax_model_ceff.text(-0.12, 1.2, 'c)', transform=ax_model_ceff.transAxes,
                   fontsize=18, fontweight='bold')

ax_model_x.set_ylim(0,NY_0p5*dx_0p5/1e3)
ax_model_c.set_ylim(0,NY_0p5*dx_0p5/1e3)
ax_model_ceff.set_ylim(0,NY_0p5*dx_0p5/1e3)

ax_model_x.set_xlabel("Sound speed [m/s]")
ax_model_c.set_xlabel("Horizontal Wind [m/s]")
ax_model_ceff.set_xlabel("Eff. wave speed [m/s]")

ax_model_x.set_ylabel("Altitude [km]")
ax_model_c.set_ylabel("Altitude [km]")
ax_model_ceff.set_ylabel("Altitude [km]")

ax_model_x.grid()
ax_model_c.grid()
ax_model_ceff.grid()

# =========================
# Ligne 2
# =========================
ax_ts_0p5 = fig.add_subplot(gs[1, 0:3])
ax_K_all  = fig.add_subplot(gs[1, 3:6])

# plot
ax_ts_0p5.plot(t_0p5, ts_0p5, 'k')
ax_ts_0p5.axvline(win_lim_0p5_direct[0], color="r", linestyle=":", linewidth=2)
#ax_ts_0p5.axvline(win_lim_0p5_direct[0]+win_lim_0p5_direct[1], color="r", linestyle=":")
#ax_ts_0p5.axvline(win_lim_0p5_refracted[0], color="r", linestyle="--")
ax_ts_0p5.axvline(win_lim_0p5_refracted[0]+win_lim_0p5_refracted[1], color="r", linestyle="--", linewidth=2)
ax_ts_0p5.set_xlabel("Time [s]")
ax_ts_0p5.set_ylabel("Pressure [Pa]")

mask = (t_0p5 > win_lim_0p5_direct[0]) & (t_0p5 < win_lim_0p5_direct[0]+win_lim_0p5_direct[1])
#ax_ts_0p5.plot(t_0p5[mask], ts_0p5[mask], color='#009E73', solid_capstyle="butt")
ax_ts_0p5.axvspan(xmin=t_0p5[mask][0], xmax=t_0p5[mask][-1], color='#009E73', alpha=0.5)
mask = (t_0p5 > win_lim_0p5_refracted[0]) & (t_0p5 < win_lim_0p5_refracted[0]+win_lim_0p5_refracted[1])
#ax_ts_0p5.plot(t_0p5[mask], ts_0p5[mask], color='#E69F00',solid_capstyle="butt")
ax_ts_0p5.axvspan(xmin=t_0p5[mask][0], xmax=t_0p5[mask][-1], color='#E69F00', alpha=0.5)

ax_ts_0p5.set_xlim(80,170)

ax_ts_0p5.text(0.50, 0.05, r'$W_3$', transform=ax_ts_0p5.transAxes,
              fontsize=18, fontweight='bold', color="r")

ax_ts_0p5.text(0.45, 0.8, r'$W_1$', transform=ax_ts_0p5.transAxes,
              fontsize=18, fontweight='bold', color="#009E73")

ax_ts_0p5.text(0.65, 0.8, r'$W_2$', transform=ax_ts_0p5.transAxes,
              fontsize=18, fontweight='bold', color="#E69F00")


# pcolormesh
maxval = np.max(np.max(abs(K_0p5_all))/100)
pcma = ax_K_all.pcolormesh(X_0p5, Y_0p5, K_0p5_all,
                          shading='auto',
                          cmap=cm.bwr,
                          vmin=-maxval, vmax=maxval,
                          rasterized=True
                          )
cbformat = ticker.ScalarFormatter()
cbformat.set_scientific('%0.2e')
cbformat.set_powerlimits((0,2))
cbformat.set_useMathText(True)
cbar = fig.colorbar(pcma,label="Sensitivity\n[$s^3/m$]", format=cbformat)

ax_K_all.set_xlabel("Range [km]")
ax_K_all.set_ylabel("Altitude [km]")

ax_ts_0p5.text(-0.12, 1.2, 'd)', transform=ax_ts_0p5.transAxes,
               fontsize=18, fontweight='bold')
ax_K_all.text(-0.12, 1.2, 'e)', transform=ax_K_all.transAxes,
              fontsize=18, fontweight='bold')

ax_K_all.text(0.01, 0.8, r'$W_3$', transform=ax_K_all.transAxes,
              fontsize=18, fontweight='bold', color="r")

# =========================
# Ligne 3
# =========================
ax_K_direct = fig.add_subplot(gs[2, 0:3])
ax_K_refrac = fig.add_subplot(gs[2, 3:6])

maxval = np.max(np.max(abs(K_0p5_direct))/100)
pcmd = ax_K_direct.pcolormesh(X_0p5, Y_0p5, K_0p5_direct,
                       shading='auto',
                       cmap=cm.bwr,
                       vmin=-maxval, vmax=maxval,
                       rasterized=True
                       )
cbformat = ticker.ScalarFormatter()
cbformat.set_scientific('%0.2e')
cbformat.set_powerlimits((0,2))
cbformat.set_useMathText(True)
cbar = fig.colorbar(pcmd,label="Sensitivity\n[$s^3/m$]", format=cbformat)

ax_K_direct.set_xlabel("Range [km]")
ax_K_direct.set_ylabel("Altitude [km]")

maxval = np.max(np.max(abs(K_0p5_refracted))/100)
pcmr = ax_K_refrac.pcolormesh(X_0p5, Y_0p5, K_0p5_refracted,
                       shading='auto',
                       cmap=cm.bwr,
                       vmin=-maxval, vmax=maxval,
                       rasterized=True
                       )
cbformat = ticker.ScalarFormatter()
cbformat.set_scientific('%0.2e')
cbformat.set_powerlimits((0,2))
cbformat.set_useMathText(True)
cbar = fig.colorbar(pcmr,label="Sensitivity\n[$s^3/m$]", format=cbformat)


ax_K_refrac.set_xlabel("Range [km]")
ax_K_refrac.set_ylabel("Altitude [km]")


ax_K_direct.text(-0.12, 1.2, 'f)', transform=ax_K_direct.transAxes,
                 fontsize=18, fontweight='bold')
ax_K_refrac.text(-0.12, 1.2, 'g)', transform=ax_K_refrac.transAxes,
                 fontsize=18, fontweight='bold')

ax_K_direct.text(0.01, 0.8, r'$W_1$', transform=ax_K_direct.transAxes,
              fontsize=18, fontweight='bold', color="#009E73")

ax_K_refrac.text(0.01, 0.8, r'$W_2$', transform=ax_K_refrac.transAxes,
              fontsize=18, fontweight='bold', color="#E69F00")

# =========================
# Ligne 4
# =========================
ax_ts_0p05 = fig.add_subplot(gs[3, 0:3])
ax_K_0p05  = fig.add_subplot(gs[3, 3:6])

ax_ts_0p05.plot(t_0p05, ts_0p05, 'k')
ax_ts_0p05.axvline(win_lim_0p05[0], color="r", linestyle=":", linewidth=2)
ax_ts_0p05.axvline(win_lim_0p05[0]+win_lim_0p05[1], color="r", linestyle="--", linewidth=2)
ax_ts_0p05.set_xlabel("Time [s]")
ax_ts_0p05.set_ylabel("Pressure [Pa]")

ax_ts_0p05.set_xlim(80,170)

ax_ts_0p05.text(0.83, 0.8, r'$W_4$', transform=ax_ts_0p05.transAxes,
              fontsize=18, fontweight='bold', color="r")

maxval = np.max(np.max(abs(K_0p05))/100)
pcm = ax_K_0p05.pcolormesh(X_0p05, Y_0p05, K_0p05,
                     shading='auto',
                     cmap=cm.bwr,
                     vmin=-maxval, vmax=maxval,
                     rasterized=True
                     )
cbformat = ticker.ScalarFormatter()
cbformat.set_scientific('%0.2e')
cbformat.set_powerlimits((0,2))
cbformat.set_useMathText(True)
cbar = fig.colorbar(pcm,label="Sensitivity\n[$s^3/m$]", format=cbformat)

ax_K_0p05.set_xlabel("Range [km]")
ax_K_0p05.set_ylabel("Altitude [km]")

ax_ts_0p05.text(-0.12, 1.2, 'h)', transform=ax_ts_0p05.transAxes,
                fontsize=18, fontweight='bold')
ax_K_0p05.text(-0.12, 1.2, 'i)', transform=ax_K_0p05.transAxes,
               fontsize=18, fontweight='bold')

ax_K_0p05.text(0.01, 0.8, r'$W_4$', transform=ax_K_0p05.transAxes,
              fontsize=18, fontweight='bold', color="r")


for ax in [ax_K_0p05, ax_K_refrac, ax_K_direct, ax_K_all]:
    ax.yaxis._set_lim(0,25,auto=False)

plt.tight_layout()
plt.show()
fig.savefig("kernel_crosscorr.pdf")
