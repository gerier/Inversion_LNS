import numpy as np
import matplotlib.pyplot as plt


plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.autolayout': True})


file_name = "./OUTPUT/pressure_file_001.dat"
source_name = "./OUTPUT/source_time_function_model.dat"

data = np.loadtxt(file_name)
#source = np.loadtxt(source_name)

plt.plot(data[:,0]/60,data[:,1])
plt.xlabel("Time (min)")
plt.ylabel("Pressure (Pa)")
#plt.figure()
#plt.plot(source[:,0]/60,source[:,1])
plt.show()
