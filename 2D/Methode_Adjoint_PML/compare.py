import matplotlib.pyplot as plt
import numpy as np


v0 = np.loadtxt("../Methode_Adjoint/OUTPUT/K_windy.txt")
v1 = np.loadtxt("./OUTPUT/K_windy.txt")

plt.imshow(v0)
plt.colorbar()
plt.show()


plt.imshow(v1)
plt.colorbar()
plt.show()

plt.imshow(v0 - v1)
plt.colorbar()
plt.show()
