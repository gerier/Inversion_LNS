import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

iteration = 99

list_rho_file = glob.glob("./OUTPUT_INVERSION/rho0_%06d_*"%(iteration))
list_p_file = glob.glob("./OUTPUT_INVERSION/p0_%06d_*"%(iteration))
list_wx_file = glob.glob("./OUTPUT_INVERSION/windx_%06d_*"%(iteration))
list_wy_file = glob.glob("./OUTPUT_INVERSION/windy_%06d_*"%(iteration))

NX = 400
demi_NX = int(NX/2) 
NY = 200
demi_NY = int(NY/2)
rho = np.zeros((NX,NY))
p = np.zeros((NX,NY))
wx = np.zeros((NX,NY))
wy = np.zeros((NX,NY))


		
for fich in list_rho_file:
	print(fich)
	data = np.loadtxt(fich)
	if "0001." in fich:
		rho[demi_NX:, :demi_NY] = data
	elif "0002." in fich:
		rho[:demi_NX,demi_NY:] = data
	elif "0003." in fich:
		rho[demi_NX:,demi_NY:] = data
	else : 
		rho[:demi_NX,:demi_NY] = data

		
plt.figure()
plt.imshow(rho) 
plt.title("Density")
plt.colorbar()
plt.show()



for fich in list_p_file:

	data = np.loadtxt(fich)
	
	if "0001." in fich:
		p[demi_NX:, :demi_NY] = data
	elif "0002." in fich:
		p[:demi_NX,demi_NY:] = data
	elif "0003." in fich:
		p[demi_NX:,demi_NY:] = data
	else : 
		p[:demi_NX,:demi_NY] = data
		
plt.figure()
plt.imshow(p)
plt.title("Pressure")
plt.colorbar()
plt.show()


plt.figure()
plt.imshow(np.sqrt(p * 1.4 / rho), cmap='jet')
plt.title("Celerity")
plt.colorbar()
plt.show()


plt.figure()
plt.imshow(np.log(p))
plt.title("Log Pressure")
plt.colorbar()
plt.show()

for fich in list_wx_file:

	data = np.loadtxt(fich)
	
	if "0001." in fich:
		wx[demi_NX:, :demi_NY] = data
	elif "0002." in fich:
		wx[:demi_NX,demi_NY:] = data
	elif "0003." in fich:
		wx[demi_NX:,demi_NY:] = data
	else : 
		wx[:demi_NX,:demi_NY] = data
		
plt.figure()
plt.imshow(wx)
plt.title("Wind x")
plt.colorbar()
plt.show()



for fich in list_wy_file:

	data = np.loadtxt(fich)
	
	if "0001." in fich:
		wy[demi_NX:, :demi_NY] = data
	elif "0002." in fich:
		wy[:demi_NX,demi_NY:] = data
	elif "0003." in fich:
		wy[demi_NX:,demi_NY:] = data
	else : 
		wy[:demi_NX,:demi_NY] = data
		
plt.figure()
plt.imshow(wy)
plt.title("Wind y")
plt.colorbar()
plt.show()
