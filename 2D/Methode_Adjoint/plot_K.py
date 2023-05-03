import numpy as np
import matplotlib.pyplot as plt


#K = np.loadtxt("OUTPUT/Kernel.txt")

index_source_x = 100
index_source_y = 100

index_receiver_x = 130
index_receiver_y = 100
 
index_start_x = 80 
index_start_y = 70
index_end_x = 150
index_end_y = 130


K = np.loadtxt("OUTPUT/Kvy.txt")
K = K.T


K_plot = K[index_start_y:index_end_y, index_start_x:index_end_x]
plt.imshow(K_plot, origin="lower")
plt.plot((index_source_x-index_start_x),(index_source_y-index_start_y), 'xr')
plt.plot((index_receiver_x-index_start_x),(index_receiver_y-index_start_y), '^g')
plt.colorbar()
plt.set_cmap('coolwarm')
plt.show()

plt.imshow(K_plot, origin="lower")
plt.plot((index_source_x-index_start_x),(index_source_y-index_start_y), 'xr')
plt.plot((index_receiver_x-index_start_x),(index_receiver_y-index_start_y), '^g')
plt.colorbar()
plt.clim(-1e-7,1e-7)
plt.set_cmap('coolwarm')
plt.show()

K_plot = K_plot[20:-20,20:-10]
plt.imshow(K_plot, origin="lower")
plt.plot((index_source_x-index_start_x-20),(index_source_y-index_start_y-20), 'xr')
plt.plot((index_receiver_x-index_start_x-20),(index_receiver_y-index_start_y-20), '^g')
plt.colorbar()
plt.clim(-1e-7,1e-7)
plt.set_cmap('coolwarm')
plt.show()
