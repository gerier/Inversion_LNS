import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype="low", analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

fcut = 0.1/4
fs = 1/dt
s_pb = butter_lowpass_filter(a, fcut, fs, order=1)

plt.plot(a, color='silver', label='Signal')
plt.plot(s_pb, color='#cc0000', label='Signal filtré')
plt.grid(True, which='both')
plt.legend(loc="best")
plt.title("Filtre passe-bas du 1er ordre")
plt.show()


#%% 

a = [history_obs[t][10] for t in range(len(history_obs)) ]
import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype="low", analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

fcut = 0.1/4
fs = 1/dt
s_pb = butter_lowpass_filter(a, fcut, fs, order=1)
s_2 =  butter_lowpass_filter(a[550:600], fcut, fs, order=1)
index = np.arange(0,len(a))

plt.figure(figsize=(10,10))
plt.plot(index,a, color='silver', label='Signal')
plt.plot(index[550:600],a[550:600], color='grey', label='Signal sur 50 pas de temps à filtrer')
plt.plot(index,s_pb, color='#cc0000', label='Signal filtré')
plt.plot(index[550:600],s_2, color='blue', label='Signal filtré sur 50 pas de temps')
plt.grid(True, which='both')
plt.legend(loc="best")
plt.title("Filtre passe-bas du 1er ordre")
plt.xlim(500,600)
plt.show()