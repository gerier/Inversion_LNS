import numpy as np
import matplotlib.pyplot as plt
# attention, pas nécessairement le meilleur NLCG

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

data_lbfgs = np.loadtxt("iterations_informations.txt")
data_nlcg = np.loadtxt("/home/deos/s.gerier/Bureau/Inversion_LNS/1D/Inversion_3/Resultats/Test_7_reg_0.005/iterations_informations.txt")


fig,ax = plt.subplots(3,1)

ax[0].loglog(data_nlcg[:,4], label="NLCG")
ax[0].loglog(data_lbfgs[:,4], label="LBFGS")
ax[0].set_ylabel("Misfit")


ax[1].plot(np.cumsum(data_nlcg[:,7]), label="NLCG")
ax[1].plot(np.cumsum(data_lbfgs[:,7]), label="LBFGS")
ax[1].set_ylabel("Nb Times in Fx")

ax[2].plot(np.cumsum(data_nlcg[:,8]), label="NLCG")
ax[2].plot(np.cumsum(data_lbfgs[:,8]), label="LBFGS")
ax[2].set_ylabel("Nb Times in DFx")

plt.legend()


plt.show()
