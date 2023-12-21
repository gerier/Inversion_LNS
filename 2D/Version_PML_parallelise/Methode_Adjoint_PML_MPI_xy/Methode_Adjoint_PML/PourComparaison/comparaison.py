import numpy as np
import matplotlib.pyplot as plt

path_parallel_0 = "../OUTPUT/p_true_000000.txt"
path_parallel_1 = "../OUTPUT/p_true_000001.txt"
path_parallel_2 = "../OUTPUT/p_true_000002.txt"
path_parallel_3 = "../OUTPUT/p_true_000003.txt"
#path_parallel_0 = "../OUTPUT/p_test_000000.txt"
#path_parallel_1 = "../OUTPUT/p_test_000001.txt"
#path_parallel_2 = "../OUTPUT/p_test_000002.txt"
#path_parallel_3 = "../OUTPUT/p_test_000003.txt"




path_serial = "./OUTPUT/p_true.txt"
#path_serial = "./OUTPUT/p_test.txt"


p_parallel_0 = np.loadtxt(path_parallel_0)
p_parallel_1 = np.loadtxt(path_parallel_1)
p_parallel_2 = np.loadtxt(path_parallel_2)
p_parallel_3 = np.loadtxt(path_parallel_3)
p_serial = np.loadtxt(path_serial)

n,m = np.shape(p_serial)
p_parallel = np.zeros((n,m))

if False : 
    p_parallel = p_parallel_0
elif True :
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

path_parallel_0 = "../OUTPUT/pa_true_000000.txt"
path_parallel_1 = "../OUTPUT/pa_true_000001.txt"
path_parallel_2 = "../OUTPUT/pa_true_000002.txt"
path_parallel_3 = "../OUTPUT/pa_true_000003.txt"
path_parallel_0 = "../OUTPUT/Kwindx_true_000000.txt"
path_parallel_1 = "../OUTPUT/Kwindx_true_000001.txt"
path_parallel_2 = "../OUTPUT/Kwindx_true_000002.txt"
path_parallel_3 = "../OUTPUT/Kwindx_true_000003.txt"

path_serial = "./OUTPUT/pa.txt"
path_serial = "./OUTPUT/K_windx.txt"

p_parallel_0 = np.loadtxt(path_parallel_0)
p_parallel_1 = np.loadtxt(path_parallel_1)
p_parallel_2 = np.loadtxt(path_parallel_2)
p_parallel_3 = np.loadtxt(path_parallel_3)
p_parallel = np.zeros((n,m))
p_serial = np.loadtxt(path_serial)

if False : 
    p_parallel = p_parallel_0
elif True :

    p_parallel[:int(n/2),:int(m/2)] = p_parallel_0
    p_parallel[int(n/2):,:int(m/2)] = p_parallel_1
    p_parallel[:int(n/2),int(m/2):] = p_parallel_2
    p_parallel[int(n/2):,int(m/2):] = p_parallel_3
    
    
plt.imshow(p_parallel)
plt.figure()
plt.imshow(p_serial)
plt.title("p serial")
plt.figure()
plt.imshow(p_parallel - p_serial)
plt.figure()
plt.imshow((p_parallel - p_serial) == 0)
plt.show()

print(np.max(abs(p_parallel - p_serial)))



sis = np.loadtxt("../OUTPUT/pressure_file_001.dat")
sis2 = np.loadtxt("./OUTPUT/pressure_file_001.dat")
plt.figure()

plt.plot(sis[:,0],sis[:,1])
plt.plot(sis2[:,0],sis2[:,1])
plt.show()

