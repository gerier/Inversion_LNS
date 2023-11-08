import matplotlib.pyplot as plt
import numpy as np


v0 = np.loadtxt("/home/deos/s.gerier/PROJECTS/SIMULATIONS/DF_2D/Roland/Test/p_true.txt")
v1 = np.loadtxt("./OUTPUT/p_true.txt")

plt.imshow(v0)
plt.colorbar()
plt.show()


plt.imshow(v1)
plt.colorbar()
plt.show()

plt.imshow(v0 - v1)
plt.colorbar()
plt.show()
