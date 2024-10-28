import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

iteration = 60
output =""
list_rho_file = glob.glob("./OUTPUT_INVERSION"+output+"/rho0_%06d_*"%(iteration))
list_p_file = glob.glob("./OUTPUT_INVERSION"+output+"/p0_%06d_*"%(iteration))
list_wx_file = glob.glob("./OUTPUT_INVERSION"+output+"/windx_%06d_*"%(iteration))
list_wy_file = glob.glob("./OUTPUT_INVERSION"+output+"/windy_%06d_*"%(iteration))

list_rho_file_prior = glob.glob("./OUTPUT_INVERSION"+output+"/rho0_000000_*")
list_p_file_prior = glob.glob("./OUTPUT_INVERSION"+output+"/p0_000000_*")
list_wx_file_prior = glob.glob("./OUTPUT_INVERSION"+output+"/windx_000000_*")
list_wy_file_prior = glob.glob("./OUTPUT_INVERSION"+output+"/windy_000000_*")

list_Krho_file = glob.glob("./OUTPUT_INVERSION"+output+"/Krho0_%06d_*"%(iteration))
list_Kp_file = glob.glob("./OUTPUT_INVERSION"+output+"/Kp0_%06d_*"%(iteration))
list_Kwx_file = glob.glob("./OUTPUT_INVERSION"+output+"/Kwindx_%06d_*"%(iteration))
list_Kwy_file = glob.glob("./OUTPUT_INVERSION"+output+"/Kwindy_%06d_*"%(iteration))

list_rho_true_file = glob.glob("./MODELS/rho0_true_*.txt")
list_p_true_file = glob.glob("./MODELS/p0_true_*.txt")
list_wx_true_file = glob.glob("./MODELS/windx_true_*.txt")
list_wy_true_file = glob.glob("./MODELS/windy_true_*.txt")

NX = 504
NPROCX = 12
dx = 100

NY = 396
dy = 100 
NPROCY = 12

rho = np.zeros((NX,NY))
p = np.zeros((NX,NY))
wx = np.zeros((NX,NY))
wy = np.zeros((NX,NY))
Krho = np.zeros((NX,NY))
Kp = np.zeros((NX,NY))
Kwx = np.zeros((NX,NY))
Kwy = np.zeros((NX,NY))
rho_true = np.zeros((NX,NY))
p_true = np.zeros((NX,NY))
wx_true = np.zeros((NX,NY))
wy_true = np.zeros((NX,NY))

LOC_NX = int(NX / NPROCX)
LOC_NY = int(NY / NPROCY)
X = np.arange(0,(NX+1)*dx,dx)/1000
Y = np.arange(0,(NY+1)*dy,dy)/1000

for fich in list_rho_file:
	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	rho[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data

for fich in list_Krho_file:
	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	Krho[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data
		
for fich in list_rho_true_file:
	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	rho_true[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data

fig, ax = plt.subplots(1,2, sharey=True)#, layout='constrained')
vmin = min(np.min(rho_true),np.min(rho))
vmax = max(np.max(rho_true),np.max(rho))

im0 = ax[0].pcolormesh(X,Y,rho_true.T,vmin=vmin,vmax=vmax) 
im1 = ax[1].pcolormesh(X,Y,rho.T,vmin=vmin,vmax=vmax) 
plt.suptitle("Density")
plt.colorbar(im1,ax=ax[1],shrink=0.8)

plt.figure()
plt.imshow(Krho.T)
plt.title("Krho")
plt.colorbar()


plt.figure()
plt.plot(rho[100,:])
plt.plot(rho_true[100,:])
plt.title("Density")
plt.show()



for fich in list_p_file:
	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	p[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data
	
	
for fich in list_Kp_file:
	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	Kp[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data
	
	
for fich in list_p_true_file:
	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	p_true[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data

		
fig, ax = plt.subplots(1,2, sharey=True)#, layout='constrained')
vmin = min(np.min(p_true),np.min(p))
vmax = max(np.max(p_true),np.max(p))
im0 = ax[0].pcolormesh(X,Y,p_true.T,vmin=vmin,vmax=vmax)
im1 = ax[1].pcolormesh(X,Y,p.T,vmin=vmin,vmax=vmax)
plt.suptitle("Pressure")
plt.colorbar(im1,ax=ax[1],shrink=0.8)

plt.figure()
plt.imshow(Kp.T)
plt.title("Kp")
plt.colorbar()

plt.figure()
plt.plot(p[100,:])
plt.plot(p_true[100,:])
plt.title("Pressure")
plt.show()


fig, ax = plt.subplots(1,2, sharey=True)#, layout='constrained')
c_true = np.sqrt(p_true.T * 1.4 / rho_true.T)
c = np.sqrt(p.T * 1.4 / rho.T)
vmin = min(np.min(c_true),np.min(c))
vmax = max(np.max(c_true),np.max(c))
im0 = ax[0].pcolormesh(X,Y,c_true ,vmin=vmin,vmax=vmax,cmap='jet')
im1 = ax[1].pcolormesh(X,Y,c,vmin=vmin,vmax=vmax,cmap='jet')
plt.suptitle("Celerity")
plt.colorbar(im1,ax=ax[1],shrink=0.8)

plt.figure()
plt.plot(np.sqrt(p[100,:] * 1.4 / rho[100,:]))
plt.plot(np.sqrt(p_true[100,:] * 1.4 / rho_true[100,:]))
plt.title("Celerity")
plt.show()


fig, ax = plt.subplots(1,2, sharey=True)#, layout='constrained')
vmin = np.log(min(np.min(rho_true),np.min(rho)))
vmax = np.log(max(np.max(rho_true),np.max(rho)))
im0 = ax[0].pcolormesh(X,Y,np.log(p_true.T),vmin=vmin,vmax=vmax)
im1 = ax[1].pcolormesh(X,Y,np.log(p.T),vmin=vmin,vmax=vmax)
plt.suptitle("Log Pressure")
plt.colorbar(im1,ax=ax[1],shrink=0.8)
plt.show()

for fich in list_wx_file:

	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	wx[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data

for fich in list_Kwx_file:

	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	Kwx[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data


for fich in list_wx_true_file:

	print(fich)
	data = np.loadtxt(fich)
	name = fich.split("_")[-1]
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	wx_true[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data


		
fig, ax = plt.subplots(1,2, sharey=True)#, layout='constrained')

vmin = min(np.min(wx_true),np.min(wx))
vmax = max(np.max(wx_true),np.max(wx))
if vmin == vmax :
	vmin = -0.1
	vmax = 0.1
	
im0 = ax[0].pcolormesh(X,Y,wx_true.T,vmin=vmin,vmax=vmax)
im1 = ax[1].pcolormesh(X,Y,wx.T,vmin=vmin,vmax=vmax)
plt.suptitle("Wind x")
plt.colorbar(im1,ax=ax[1],shrink=0.8)

plt.figure()
plt.imshow(Kwx.T)
plt.title("Kwx")
plt.colorbar()

plt.show()

fig, ax = plt.subplots(1,2, sharey=True)#, layout='constrained')
vmin = min(np.min(c_true)+np.min(abs(wx_true)),np.min(c)+np.min(abs(wx)))
vmax = max(np.max(c_true)+np.max(abs(wx_true)),np.max(c)+np.max(abs(wx)))

im0 = ax[0].pcolormesh(X,Y,c_true + wx_true.T,vmin=vmin,vmax=vmax)
im1 = ax[1].pcolormesh(X,Y,c+wx.T,vmin=vmin,vmax=vmax)
plt.suptitle("Effective wave speed")
plt.colorbar(im1,ax=ax[1],shrink=0.8)

plt.figure()
plt.plot(c.T[100,:] + wx[100,:])
plt.plot(c_true.T[100,:] + wx_true[100,:])
plt.title("Effective wave speed")
plt.show()
