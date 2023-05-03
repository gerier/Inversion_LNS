



stat_1 = 130
stat_2 = 60
difftime = 0

for stat in [stat_1,stat_2]:
    plt.plot([history_obs[t].rho[stat] for t in range(len(history_obs))])
    tmax = np.argmax([history_obs[t].rho[stat] for t in range(len(history_obs))])
    plt.plot(tmax, history_obs[tmax].rho[stat], 'xr')
    print(tmax * 0.05)
    if stat == stat_1 : 
        difftime -= t_ax[tmax] 
    else: 
        difftime += t_ax[tmax]  

print("Diff de temps : ", difftime)
print("Diff de temps estimé :", (z[stat_2]-z[stat_1])/(c[0]-GT_v0[0]))
print("Diff de temps estimé avec div/2 :", -(z[stat_2]-z[stat_1])/(c[0]-GT_v0[0]/2.))



stat_1 = 220
stat_2 = 170
difftime = 0

for stat in [stat_1,stat_2]:
    plt.plot([history_obs[t].rho[stat] for t in range(len(history_obs))])
    tmax = np.argmax([history_obs[t].rho[stat] for t in range(len(history_obs))])
    plt.plot(tmax, history_obs[tmax].rho[stat], 'xr')
    print(tmax * 0.05)
    if stat == stat_1 : 
        difftime -= t_ax[tmax] 
    else: 
        difftime += t_ax[tmax]  

print("Diff de temps : ", -difftime)
print("Diff de temps estimé :", -(z[stat_2]-z[stat_1])/(c[0]+GT_v0[0]))
print("Diff de temps estimé avec div/2 :", -(z[stat_2]-z[stat_1])/(c[0]+GT_v0[0]/2.))