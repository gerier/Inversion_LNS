import numpy as np
import matplotlib.pyplot as plt

path_parallel_0 = "../OUTPUT/Kp0_true_000000.txt"
path_parallel_1 = "../OUTPUT/Kp0_true_000001.txt"
path_parallel_2 = "../OUTPUT/Kp0_true_000002.txt"
path_parallel_3 = "../OUTPUT/Kp0_true_000003.txt"

path_serial = "./OUTPUT/K_p0.txt"


p_parallel_0 = np.loadtxt(path_parallel_0)
p_parallel_1 = np.loadtxt(path_parallel_1)
p_parallel_2 = np.loadtxt(path_parallel_2)
p_parallel_3 = np.loadtxt(path_parallel_3)
p_serial = np.loadtxt(path_serial)

n,m = np.shape(p_serial)
p_parallel = np.zeros((n,m))

if True : 
    p_parallel[:int(n/2),:int(m/2)] = p_parallel_0
    p_parallel[int(n/2):,:int(m/2)] = p_parallel_1
    p_parallel[:int(n/2),int(m/2):] = p_parallel_2
    p_parallel[int(n/2):,int(m/2):] = p_parallel_3

elif False: 
    p_parallel[:int(n/4),:] = p_parallel_0
    p_parallel[int(n/4):2*int(n/4),:] = p_parallel_1
    p_parallel[2*int(n/4):3*int(n/4),:] = p_parallel_2
    p_parallel[3*int(n/4):,:] = p_parallel_3
else : 
    p_parallel[:,:int(m/4)] = p_parallel_0
    p_parallel[:,int(m/4):2*int(m/4)] = p_parallel_1
    p_parallel[:,2*int(m/4):3*int(m/4)] = p_parallel_2
    p_parallel[:,3*int(m/4):] = p_parallel_3


plt.imshow(p_parallel)
plt.figure()
plt.imshow(p_serial)
plt.figure()
plt.imshow(p_parallel - p_serial)
plt.show()

print(np.max(abs(p_parallel - p_serial)))
