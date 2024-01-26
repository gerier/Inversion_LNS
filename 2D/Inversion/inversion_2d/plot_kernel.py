import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

iteration = 1

list_rho_file = glob.glob("./OUTPUT_INVERSION/Krho0_%06d_*"%(iteration))
list_p_file = glob.glob("./OUTPUT_INVERSION/Kp0_%06d_*"%(iteration))
list_wx_file = glob.glob("./OUTPUT_INVERSION/Kwindx_%06d_*"%(iteration))
list_wy_file = glob.glob("./OUTPUT_INVERSION/Kwindy_%06d_*"%(iteration))
#list_rho_file = glob.glob("./OUTPUT/Krho0_true2_*")
#list_p_file = glob.glob("./OUTPUT/Kp0_true2_*")
#list_wx_file = glob.glob("./OUTPUT/Kwindx_true2_*")
#list_wy_file = glob.glob("./OUTPUT/Kwindy_true2_*")

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

print(p[99:102,99:102])

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
