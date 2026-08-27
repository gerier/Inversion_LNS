import numpy as np
import pandas as pd
#from geopy.distance import geodesic
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import os
#from pyproj import Geod
from matplotlib import ticker

#sph_proj = Geod(ellps='sphere')

plt.rcParams.update({'font.size': 18})

dir_output = "./OUTPUT/"
list_kernel_path_p0   = glob.glob(dir_output + '/Kp0*')
list_kernel_path_rho0 = glob.glob(dir_output + '/Krho0*')
list_kernel_path_wx   = glob.glob(dir_output + '/Kwindx*')

nb_files = len(glob.glob(dir_output + "/pressure_file_0*.dat"))

atmo_model_prior = "./atmospheric_model_Hukkakero.dat"
atmo_model_true  = "./atmospheric_model_Hukkakero_grav.dat"

NPROCX = 32
NPROCY = 12
NX = 4224
NY = 1512
dx = 100

LOC_NX = int(NX / NPROCX)
LOC_NY = int(NY / NPROCY)

selected_stations = ["I37H6", "I37H4", "I37H2", "I37H0", "I37H3", "I37H9"]
file_path = "location_I37NO.lst"
source = [67.9335, 25.8322]  # (lat, lon)
source_lat, source_lon = source

REC_wr = np.zeros((6,2))
REC_wr[:,0] = [1302.1, 1303.0, 1304.0, 1304.2, 1305.1,1306.1]
REC_wr[:,1] = REC_wr[:,0] + 20
###############################################################################
# Load atmospheric models
###############################################################################

p0       = np.zeros((NX, NY))
p0_prior = np.zeros((NX, NY))
rho0       = np.zeros((NX, NY))
rho0_prior = np.zeros((NX, NY))
c        = np.zeros_like(rho0)
c_prior  = np.zeros_like(rho0_prior)
wx       = np.zeros_like(rho0)
wx_prior = np.zeros_like(rho0_prior)
gamma    = np.zeros((NX, NY))

data = np.loadtxt(atmo_model_prior, skiprows=3, delimiter=',')
for j in range(NY):
    rho0_prior[:, j] = data[j, 1]
    p0_prior[:, j]   = data[j, 4]
    c_prior[:, j]    = data[j, 3]
    gamma[:, j]      = data[j, -1]
    wx_prior[:, j]   = data[j, -4]

data = np.loadtxt(atmo_model_true, skiprows=3, delimiter=',')
for j in range(NY):
    rho0[:, j]  = data[j, 1]
    p0[:, j]    = data[j, 4]
    c[:, j]     = data[j, 3]
    gamma[:, j] = data[j, -1]
    wx[:, j]    = data[j, -4]

###############################################################################
# Load kernels
###############################################################################
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

Kp = np.zeros((NX, NY))
#for fich in list_kernel_path_p0:
#    data    = np.fromfile(fich,dtype=np.float32)
#    data = data.reshape((LOC_NX, LOC_NY))
#    num_proc = int(fich.split("_")[-1].split(".")[0])
#    i = int(num_proc % NPROCX)
#    j = int(num_proc / NPROCX)
#    Kp[i*LOC_NX:(i+1)*LOC_NX, j*LOC_NY:(j+1)*LOC_NY] = data
#Kp = Kp * dx * dx
Kp = load_global_field(
    "Kp0",
    nproc_x=NPROCX,   
    nproc_y=NPROCY,  
    nx_global = NX,
    ny_global = NY,
    folder="./OUTPUT/"
)

Krho = np.zeros((NX, NY))
#for fich in list_kernel_path_rho0:
#    data    = np.fromfile(fich,dtype=np.float32)
#    data = data.reshape((LOC_NX, LOC_NY))
#    num_proc = int(fich.split("_")[-1].split(".")[0])
#    i = int(num_proc % NPROCX)
#    j = int(num_proc / NPROCX)
#    Krho[i*LOC_NX:(i+1)*LOC_NX, j*LOC_NY:(j+1)*LOC_NY] = data
#Krho = Krho * dx * dx
Krho = load_global_field(
    "Krho0",
    nproc_x=NPROCX,   
    nproc_y=NPROCY,  
    nx_global = NX,
    ny_global = NY,
    folder="./OUTPUT/"
)

Kwx = np.zeros((NX, NY))
#for fich in list_kernel_path_wx:
#    data    = np.fromfile(fich,dtype=np.float32)
#    data = data.reshape((LOC_NX, LOC_NY))
#    num_proc = int(fich.split("_")[-1].split(".")[0])
#    i = int(num_proc % NPROCX)
#    j = int(num_proc / NPROCX)
#    Kwx[i*LOC_NX:(i+1)*LOC_NX, j*LOC_NY:(j+1)*LOC_NY] = data
#Kwx = Kwx * dx * dx
Kwx = load_global_field(
    "Kwindx",
    nproc_x=NPROCX,   
    nproc_y=NPROCY,  
    nx_global = NX,
    ny_global = NY,
    folder="./OUTPUT/"
)

Kc = 2 * c_prior.T * rho0_prior.T / gamma.T * Kp

###############################################################################
# Axis helpers
###############################################################################

X = np.linspace(-50, (NX - 1) / 10 - 50, NX)
Y = np.linspace(0, (NY - 1) / 10, NY) - 0.05

cbformat = ticker.ScalarFormatter()
cbformat.set_scientific('%0.2e')
cbformat.set_powerlimits((0, 2))
cbformat.set_useMathText(True)


def add_colorbar(fig, ax, im, label):
    fmt = ticker.ScalarFormatter()
    fmt.set_useMathText(True)
    fmt.set_powerlimits((0, 2))
    cbar = fig.colorbar(im, ax=ax, label=label, format=fmt)
    return cbar


###############################################################################
# Load stations locations
###############################################################################

import math

def haversine(source, lat, lon):
    """
    Calcule la distance en km entre source=(lat, lon)
    et un point (lat, lon).
    """
    lat1, lon1 = source

    R = 6371.0  # rayon moyen de la Terre en km

    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat)
    lon2 = math.radians(lon)

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    )

    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return R * c

df = pd.read_fwf(file_path, widths=[7,11,11], on_bad_lines='skip', skiprows=2, skipfooter=3)
index_to_keep = df['STA'].astype(str).str.contains(
    '-----|STA|SQL',
    case=False,
    na=False
)
#index_to_keep = df.STA.str.contains('|'.join(['-----', 'STA', 'SQL']), case=False)
df = df[~index_to_keep]
df[['LAT', 'LON']] = df[['LAT', 'LON']].astype('float')

#df["dist_km"] = df.apply(lambda row: geodesic(source, (row["LAT"], row["LON"])).km, axis=1)
df["dist_km"] = df.apply(
    lambda row: haversine(source, row["LAT"], row["LON"]),
    axis=1
)

df_sel = df[df["STA"].isin(selected_stations)]
df_sel = df_sel.drop_duplicates()

bary_lat = df_sel["LAT"].mean()
bary_lon = df_sel["LON"].mean()


###############################################################################
# Build the figure  (4 rows)
#
#  Row 1 : 3 subplots — wave speed model (c), wind model (wx),
#                        effective wave speed model (c + wx)
#  Row 2 : Kc 2D  |  Kc 1D (sum over rows)
#  Row 3 : Kwx 2D | Kwx 1D (sum over rows)
#  Row 0 : Lat/Lon map (placeholder) | time-series waveforms
###############################################################################

fig = plt.figure(figsize=(18, 22))

gs = fig.add_gridspec(
    4, 1,
    hspace=0.45,
    height_ratios=[1, 1, 1, 1],
)

# ── Row 0 : map placeholder + time-series ────────────────────────────────────
gs_row4 = gs[0].subgridspec(1, 2, wspace=0.45, width_ratios=[1, 1])

ax_map = fig.add_subplot(gs_row4[0])
ax_ts  = fig.add_subplot(gs_row4[1])

ax_map.scatter(df_sel["LON"], df_sel["LAT"],
               color='#D87373', s=60, zorder=5, label='Stations')

for _, row in df_sel.iterrows():
    ax_map.annotate(row["STA"],
                    xy=(row["LON"], row["LAT"]),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=16, zorder=6)

# Source
#ax_map.scatter(source[1], source[0],
#               marker='*', s=200, color='black', zorder=6, label='Source')

# Barycentre
#ax_map.scatter(bary_lon, bary_lat,
#               marker='+', s=150, color='blue', zorder=6, label='Barycentre')

# Droite passant par source et barycentre, étendue aux bords de l'axe
def line_through_two_points(x1, y1, x2, y2, x_range):
    """Retourne les y correspondant à x_range sur la droite (x1,y1)-(x2,y2)."""
    slope = (y2 - y1) / (x2 - x1) if (x2 - x1) != 0 else np.inf
    return y1 + slope * (x_range - x1)

x_min = min(df_sel["LON"].min(), source[1]) - 1
x_max = max(df_sel["LON"].max(), source[1]) + 1
x_line = np.array([x_min, x_max])
y_line = line_through_two_points(source[1], source[0], bary_lon, bary_lat, x_line)

ax_map.plot(x_line, y_line, 'k--', linewidth=1.5, alpha=0.6, label='Source–barycentre axis')

ax_map.set_xlabel("Longitude [°]")
ax_map.set_ylabel("Latitude [°]")
ax_map.legend(fontsize=9, loc='best')
ax_map.grid(True, alpha=0.3)

ax_map.set_xlim(18.575, 18.645)
ax_map.set_ylim(69.0665, 69.086)

ax_map.text(0, 1.02, "a)", transform=ax_map.transAxes, size=16, weight='bold')

# --- Time-series waveforms ---
station_offsets = []
station_names   = []
for ii in range(nb_files):
    offset = 320 + ii * 2
    station_offsets.append(offset)
    station_names.append(selected_stations[ii])

    file_syn = dir_output + "/pressure_file_%03d.dat" % (ii + 1)
    data_syn = np.loadtxt(file_syn)
    idx_syn  = (data_syn[:, 0] > 1200) & (data_syn[:, 0] < 1400)
    factor_syn = np.max(data_syn[idx_syn, 1])

    file_obs = dir_output + "/pressure_file_obs_%03d.dat" % (ii + 1)
    data_obs = np.loadtxt(file_obs)
    idx_obs  = (data_obs[:, 0] > 1200) & (data_obs[:, 0] < 1400)
    factor_obs = np.max(data_obs[idx_obs, 1])
    
    ax_ts.plot(data_obs[:, 0] / 60, data_obs[:, 1] / factor_obs + offset,
               color='#D87373', linewidth=2, alpha=0.4)
    ax_ts.plot(data_syn[:, 0] / 60, data_syn[:, 1] / factor_syn + offset,
               'k--', linewidth=2, alpha=0.4)

    idx_corr_window = (data_obs[:,0] >= REC_wr[ii,0]) & (data_obs[:,0] < REC_wr[ii,1])
    ax_ts.plot(data_obs[idx_corr_window, 0] / 60, data_obs[idx_corr_window, 1] / factor_obs + offset,
               color='#D87373', linewidth=2)
    ax_ts.plot(data_syn[idx_corr_window, 0] / 60, data_syn[idx_corr_window, 1] / factor_syn + offset,
               'k--', linewidth=2)





# Legend (dummy handles)
ax_ts.plot([], [], color='#D87373', linewidth=2, label='Data')
ax_ts.plot([], [], 'k--', linewidth=2, label='Synthetic')
ax_ts.legend(loc='upper left', framealpha=0.91)
ax_ts.set_xlabel("Time [min]")
#ax_ts.set_ylabel("Distance from source [km]")
ax_ts.set_ylim(318, 332)
ax_ts.set_yticks(station_offsets)
ax_ts.set_yticklabels(station_names)
ax_ts.set_xlim(19.5, 23)
ax_ts.text(0, 1.02, "b)", transform=ax_ts.transAxes, size=16, weight='bold')

# ── Row 1 : atmospheric models ────────────────────────────────────────────────
gs_row1 = gs[1].subgridspec(1, 3, wspace=0.45)

ax_c   = fig.add_subplot(gs_row1[0])
ax_wx  = fig.add_subplot(gs_row1[1])
ax_eff = fig.add_subplot(gs_row1[2])

for ax, profile, xlabel, label in [
    (ax_c,   c_prior[0, :],                    "Sound speed [m/s]",          "c)"),
    (ax_wx,  wx_prior[0, :],                   "Horizontal wind [m/s]",     "d)"),
    (ax_eff, c_prior[0, :] + wx_prior[0, :],   "Effective wave speed [m/s]","e)"),
]:
    ax.plot(profile, Y, '-k', linewidth=2.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Altitude [km]")
    ax.grid(axis='y')
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax.text(0, 1.02, label, transform=ax.transAxes, size=16, weight='bold')

# ── Row 2 : Kc ────────────────────────────────────────────────────────────────
gs_row2 = gs[2].subgridspec(1, 2, wspace=0.15, width_ratios=[3, 1])

ax_Kc2d = fig.add_subplot(gs_row2[0])
ax_Kc1d = fig.add_subplot(gs_row2[1])

maxval_Kc = np.max(np.abs(Kc)) / 25
im_Kc = ax_Kc2d.pcolormesh(
    X, Y, Kc,
    cmap=cm.bwr, vmin=-maxval_Kc, vmax=maxval_Kc, rasterized=True
)
ax_Kc2d.set_xlabel("Range [km]")
ax_Kc2d.set_ylabel("Altitude [km]")
ax_Kc2d.set_xlim(-15, 345)
ax_Kc2d.grid(axis='y')
ax_Kc2d.text(0, 1.02, "f)", transform=ax_Kc2d.transAxes, size=16, weight='bold')
add_colorbar(fig, ax_Kc2d, im_Kc, r"$\bar K^D_c$ - Sensitivity [$s^3/m$]")

Kc_1d = np.sum(Kc, axis=1)
ax_Kc1d.plot(Kc_1d, Y, '-k', linewidth=2.0)
ax_Kc1d.set_xlabel(r"$\mathcal{K}^D_c$ - Sensitivity [$s^3/m$]")
ax_Kc1d.set_ylabel("Altitude [km]")
ax_Kc1d.grid(axis='y')
ax_Kc1d.text(0, 1.02, "g)", transform=ax_Kc1d.transAxes, size=16, weight='bold')

# share y-limits across row 2
yl2 = ax_Kc2d.get_ylim()
ax_Kc1d.set_ylim(yl2)

# ── Row 3 : Kwx ───────────────────────────────────────────────────────────────
gs_row3 = gs[3].subgridspec(1, 2, wspace=0.15, width_ratios=[3, 1])

ax_Kwx2d = fig.add_subplot(gs_row3[0])
ax_Kwx1d = fig.add_subplot(gs_row3[1])

maxval_Kwx = np.max(np.abs(Kwx)) / 25
im_Kwx = ax_Kwx2d.pcolormesh(
    X, Y, Kwx,
    cmap=cm.bwr, vmin=-maxval_Kwx, vmax=maxval_Kwx, rasterized=True
)
ax_Kwx2d.set_xlabel("Range [km]")
ax_Kwx2d.set_ylabel("Altitude [km]")
ax_Kwx2d.set_xlim(-15, 345)
ax_Kwx2d.grid(axis='y')
ax_Kwx2d.text(0, 1.02, "h)", transform=ax_Kwx2d.transAxes, size=16, weight='bold')
add_colorbar(fig, ax_Kwx2d, im_Kwx, r"$\bar K^D_{v_0^1}$ - Sensitivity [$s^3/m$]")

Kwx_1d = np.sum(Kwx, axis=1)
ax_Kwx1d.plot(Kwx_1d, Y, '-k', linewidth=2.0)
ax_Kwx1d.set_xlabel(r"$\mathcal{K}^D_{v_0^1}$ - Sensitivity [$s^3/m$]")
ax_Kwx1d.set_ylabel("Altitude [km]")
ax_Kwx1d.grid(axis='y')
ax_Kwx1d.text(0, 1.02, "i)", transform=ax_Kwx1d.transAxes, size=16, weight='bold')

# share y-limits across row 3
yl3 = ax_Kwx2d.get_ylim()
ax_Kwx1d.set_ylim(yl3)


###############################################################################
# Save & show
###############################################################################

fig.savefig("kernels_summary.pdf", bbox_inches='tight')
plt.show()

