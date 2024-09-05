import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

list_rho_file = glob.glob("./MODELS/rho0*")
list_p_file = glob.glob("./MODELS/p0*")
list_wx_file = glob.glob("./MODELS/windx*")
list_wy_file = glob.glob("./MODELS/windy*")


NPROC_X = 24
NPROC_Y = 4
NX = 7200
demi_NX = int(NX/NPROC_X) 
NY = 1400
demi_NY = int(NY/NPROC_Y)
rho = np.zeros((NX,NY))
p = np.zeros((NX,NY))
wx = np.zeros((NX,NY))
wy = np.zeros((NX,NY))

for fich in list_rho_file:

	data = np.loadtxt(fich)
	
	name = fich.split("_")[-1]
	i_proc = int(name[:-4])
	i_rank = int(i_proc % NPROC_X)
	j_rank = int(i_proc / NPROC_X)
	
	
	rho[demi_NX*i_rank:(i_rank+1)*demi_NX,demi_NY*j_rank:(j_rank+1)*demi_NY] = data
		
plt.figure()
plt.imshow(rho)
plt.colorbar()
plt.title("Density")
plt.show()

for fich in list_p_file:

	data = np.loadtxt(fich)
	
	name = fich.split("_")[-1]
	i_proc = int(name[:-4])
	i_rank = int(i_proc % NPROC_X)
	j_rank = int(i_proc / NPROC_X)
	
	
	p[demi_NX*i_rank:(i_rank+1)*demi_NX,demi_NY*j_rank:(j_rank+1)*demi_NY] = data
		
plt.figure()
plt.imshow(p)
plt.colorbar()
plt.title("Pressure")
plt.show()


plt.figure()
plt.imshow(np.sqrt(p * 1.4 / rho))
plt.colorbar()
plt.title("Celerity")
plt.show()



for fich in list_wx_file:

	data = np.loadtxt(fich)
	
	name = fich.split("_")[-1]
	i_proc = int(name[:-4])
	i_rank = int(i_proc % NPROC_X)
	j_rank = int(i_proc / NPROC_X)
	
	
	wx[demi_NX*i_rank:(i_rank+1)*demi_NX,demi_NY*j_rank:(j_rank+1)*demi_NY] = data
		
plt.figure()
plt.imshow(wx)
plt.colorbar()
plt.title("Horizontal wind")
plt.show()
