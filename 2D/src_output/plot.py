import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

nb_files = len(glob.glob("./OUTPUT/pressure_file_0*.dat"))


fig1,ax1 = plt.subplots()
fig1,ax2 = plt.subplots()
fig1,ax3 = plt.subplots()

for i in range(0,nb_files):

	file_name = "./OUTPUT/pressure_file_%03d.dat"%(i+1)
	print(file_name)
	data = np.loadtxt(file_name)
	
	file_name = "./OUTPUT/pressure_file_obs_%03d.dat"%(i+1)
	print(file_name)
	dataobs = np.loadtxt(file_name)

	ax1.plot(data[:,0]/60,data[:,1]*10+ 330 + (i-1)*2,'k')
	ax2.plot(dataobs[:,0]/60,dataobs[:,1]*1000+ 330 + (i-1)*2,'k'),
	ax3.plot(dataobs[:,0]/60,(dataobs[:,1]-data[:,1])*1000+ 330 + (i-1)*2,'k'),


	Fs = 1/0.05
	print(Fs, data[:2,0])

start = []
t0 = 1.2/0.1

ax1.set_xlabel("Time (min)")
ax1.set_ylabel("Km from source")
ax1.set_xlim(15,24)
ax1.set_ylim(330,350)

ax2.set_xlabel("Time (min)")
ax2.set_ylabel("Km from source")
#ax2.set_xlim(1,7)
ax2.set_xlim(15,24)
ax2.set_ylim(328,342)

ax3.set_xlabel("Time (min)")
ax3.set_ylabel("Km from source")
#ax3.set_xlim(1,7)
ax3.set_ylim(328,342)
ax3.set_xlim(15,24)
plt.show()
