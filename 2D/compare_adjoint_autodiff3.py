import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

directory_adjoint = "Methode_Adjoint (3e copie)"
directory_autodiff = "Methode_Autodifferencation (copie) (copie)"

directory_adjoint = "HomogeneVent_Kp0_vent20ms/Methode_Adjoint_PML_vent"
directory_autodiff = "HomogeneVent_Kp0_vent20ms/Methode_Autodifferenciation_PML_vent"

parameter = "p0"

index_source_x = 40
index_source_y = 100

index_receiver_x = 125
index_receiver_y = 100
 
index_start_x = 80 
index_start_y = 70
index_end_x = 150
index_end_y = 130

   
fig, ax = plt.subplots(1,3, figsize=(15,4))
# Plot Kernel obtained by adjoint
K_adj = np.loadtxt(directory_adjoint+"/OUTPUT/K_"+parameter+".txt")
K_adj = K_adj.T

K_plot = K_adj[index_start_y:index_end_y, index_start_x:index_end_x]
bounds = max( np.min(abs(K_plot)), np.max(abs(K_plot)))
im0 = ax[0].imshow(K_plot, origin="lower", extent=[index_start_x/10,index_end_x/10,index_start_y/10,index_end_y/10])
ax[0].plot((index_source_x)/10,(index_source_y)/10, 'xr')
ax[0].plot((index_receiver_x)/10,(index_receiver_y)/10, '^g')
plt.colorbar(im0, ax=ax[0])
#im0.set_clim(-bounds,bounds)
im0.set_clim(-bounds/100,bounds/100)
im0.set_cmap('bwr')
ax[0].set_title("a) Kernel by\nadjoint method")
ax[0].set_xlabel("Range (km)")
ax[0].set_ylabel("Alt. (km)")

# Plot Kernel obtained by autodifferenciation
K_auto = np.loadtxt(directory_autodiff+"/OUTPUT/Kernel.txt")
K_auto = K_auto.T


K_plot = K_auto[index_start_y:index_end_y, index_start_x:index_end_x]
bounds = max( np.min(abs(K_plot)), np.max(abs(K_plot)))
im1 = ax[1].imshow(K_plot, origin="lower", extent=[index_start_x/10,index_end_x/10,index_start_y/10,index_end_y/10])
ax[1].plot((index_source_x)/10,(index_source_y)/10, 'xr')
ax[1].plot((index_receiver_x)/10,(index_receiver_y)/10, '^g')
plt.colorbar(im1, ax=ax[1])
#im1.set_clim(-bounds,bounds)
im1.set_clim(-bounds/100,bounds/100)
im1.set_cmap('bwr')
ax[1].set_title("b) Kernel by\nautodifferentiation")
ax[1].set_xlabel("Range (km)")
ax[1].set_ylabel("Alt. (km)")

# Plot the relative difference
K_diff = K_auto - K_adj
K_relative_diff = (K_diff + 1e-16) / (np.max(K_auto) + 1e-16)   #(K_adj + 1e-16)
K_plot = K_relative_diff[index_start_y:index_end_y, index_start_x:index_end_x]
bounds = max( np.min(abs(K_plot)), np.max(abs(K_plot)))

im2 = ax[2].imshow(K_plot, origin="lower", extent=[index_start_x/10,index_end_x/10,index_start_y/10,index_end_y/10])
ax[2].plot((index_source_x)/10,(index_source_y)/10, 'xr')
ax[2].plot((index_receiver_x)/10,(index_receiver_y)/10, '^g')
plt.colorbar(im2, ax=ax[2])
#im2.set_clim(-bounds,bounds)
im2.set_clim(-0.1,0.1)
im2.set_cmap('bwr')
ax[2].set_title("c) (Relative) Difference\nbetween a) and b)")
ax[2].set_xlabel("Range (km)")
ax[2].set_ylabel("Alt. (km)")
plt.show()




x = np.arange(index_start_x/10, index_end_x/10, 0.1)
z = np.arange(index_start_y/10, index_end_y/10, 0.1)



tranche_auto = K_auto[index_start_y:index_end_y,index_source_x]
tranche_adj = K_adj[index_start_y:index_end_y,index_source_x]
diff_rs = index_receiver_x - index_source_x

decalage = int(diff_rs *1/ 4)

tranche_auto = K_auto[index_start_y:index_end_y,index_source_x-decalage]
tranche_adj = K_adj[index_start_y:index_end_y,index_source_x-decalage]

fig, ax = plt.subplots(2,3, figsize=(15,8))

ax[0,0].plot(z,tranche_auto, label='autodiff')
ax[0,0].plot(z,tranche_adj, label='adj')
#ax[0,0].plot(z, abs(tranche_auto - tranche_adj) / max(tranche_auto), "--", label="Relative difference")
#ax[0,0].legend()
ax[0,0].grid()
ax[0,0].set_xlabel("Range (km)")
ax[0,0].set_title("Vertical slice,\nbefore source")


ax[0,1].plot(z,tranche_auto, label='autodiff')
ax[0,1].plot(z,tranche_adj, label='adj')
#ax[0,1].legend()
ax[0,1].grid()
ax[0,1].set_xlabel("Range (km)")
ax[0,1].set_title("Vertical slice,\non the source")


decalage = int(diff_rs / 2)

tranche_auto = K_auto[index_start_y:index_end_y,index_source_x+decalage]
tranche_adj = K_adj[index_start_y:index_end_y,index_source_x+decalage]


ax[0,2].plot(z,tranche_auto, label='autodiff')
ax[0,2].plot(z,tranche_adj, label='adj')
#ax[0,2].legend()
ax[0,2].grid()
ax[0,2].set_xlabel("Range (km)")
ax[0,2].set_title("Vertical slice between the\nsource and receiver")



decalage = diff_rs

tranche_auto = K_auto[index_start_y:index_end_y,index_source_x+decalage]
tranche_adj = K_adj[index_start_y:index_end_y,index_source_x+decalage]


ax[1,0].plot(z,tranche_auto, label='autodiff')
ax[1,0].plot(z,tranche_adj, label='adj')
#ax[1,0].legend()
ax[1,0].grid()
ax[1,0].set_xlabel("Range (km)")
ax[1,0].set_title("Vertical slice,\non the receiver")

decalage = int(diff_rs *5/ 4)

tranche_auto = K_auto[index_start_y:index_end_y,index_source_x+decalage]
tranche_adj = K_adj[index_start_y:index_end_y,index_source_x+decalage]


ax[1,1].plot(z,tranche_auto, label='autodiff')
ax[1,1].plot(z,tranche_adj, label='adj')
#ax[1,1].legend()
ax[1,1].grid()
ax[1,1].set_xlabel("Range (km)")
ax[1,1].set_title("Vertical slice,\nafter the receiver")




tranche_auto = K_auto[index_source_y,index_start_x:index_end_x]
tranche_adj = K_adj[index_source_y,index_start_x:index_end_x]


ax[1,2].plot(tranche_auto, label='autodiff')
ax[1,2].plot(tranche_adj, label='adj')
#ax[1,2].legend()
ax[1,2].grid()
ax[1,2].set_title("Horizontal Slice on the\nsource and receiver")
ax[1,2].set_xlabel("Range (km)")

plt.show()
