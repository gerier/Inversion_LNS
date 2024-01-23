import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("./OUTPUT/p_checkpoint.txt")
plt.plot(np.linspace(0,1800,1801), data[:,1])


data = np.loadtxt("./OUTPUT/pressure_file_001.dat")
plt.plot(data[::-1,1])


plt.show()
