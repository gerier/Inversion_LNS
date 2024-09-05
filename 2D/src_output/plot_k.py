import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

list_kernel_path_p0 = glob.glob('../OUTPUT/Kp0_0000??.txt')
list_model_path_p0 = glob.glob('../MODELS/p0_true*.txt')
list_model_path_p0_prior = glob.glob('MODELS/p0_true*.txt')

nb_files = len(glob.glob("../OUTPUT/pressure_file_*.dat"))


NPROCX = 1
NPROCY = 4
NX = 200
NY = 400
dx = 100
Kp = np.zeros((NX,NY)) 
LOC_NX = int(NX / NPROCX)
LOC_NY = int(NY / NPROCY)
p0 = np.zeros((NX,NY)) 
p0_prior = np.zeros((NX,NY)) 

for fich in list_kernel_path_p0:
	print(fich)
	data = np.loadtxt(fich)
	
	num_proc = int(fich.split("_")[-1].split(".")[0])
	i = int(num_proc% NPROCX)
	j = int(num_proc/ NPROCX)
	
	Kp[i*LOC_NX:(i+1)*LOC_NX,j*LOC_NY:(j+1)*LOC_NY] = data
	
Kp = Kp *dx*dx





lgKp = Kp

for K in [Kp]:#,lgKp]:
	X = np.linspace(0,(len(K[:,0])-1),len(K[:,0]))
	Y = np.linspace(0,(len(K[0,:])-1),len(K[0,:]))

	fig = plt.figure(figsize=(15, 8))
	spec = fig.add_gridspec(2, 7)

	ax0 = fig.add_subplot(spec[0, 0])
	ax0.semilogx(p0[0, :], Y,'-k', linewidth=3.0)
	ax0.semilogx(p0_prior[0, :], Y,'--k', linewidth=3.0)
	ax0.set_xlabel("Pressure [Pa]")
	ax0.set_ylabel("Altitude [km]")

	ax1 = fig.add_subplot(spec[0, 1:-1])
	maxval = np.max(np.max(abs(K))/40)	
	
	
	im1 = ax1.pcolormesh(X,Y,K.T, cmap=cm.bwr, vmin=-maxval, vmax=maxval) 
	ax1.set_xlabel("Range [km]")
	fig.colorbar(im1,label="Sensitivity")

	#ax1.set_xlim(0,NX)
	
	ax2 = fig.add_subplot(spec[0, -1], sharey=ax1)
	int_lgK = np.sum(K.T,axis=1)

	ax2.plot(int_lgK,Y, '-k', linewidth=2.0)
	ax2.set_xlabel(r"Sensitivity [$1$]")
	
	ax3 = fig.add_subplot(spec[1,:])
	

	

plt.show()



