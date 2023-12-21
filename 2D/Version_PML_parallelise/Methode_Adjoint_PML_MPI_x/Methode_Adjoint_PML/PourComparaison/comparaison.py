import numpy as np
import matplotlib.pyplot as plt


path_parallel = "../OUTPUT/p_true.txt"
path_serial = "./OUTPUT/p_true.txt"


p_parallel = np.loadtxt(path_parallel)
p_serial = np.loadtxt(path_serial)

plt.imshow(p_parallel)
plt.figure()
plt.imshow(p_serial)
plt.figure()
plt.imshow(p_parallel - p_serial)
plt.show()

print(np.max(abs(p_parallel - p_serial)))
