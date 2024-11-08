import numpy as np
from scipy.signal import butter, lfilter, freqz


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype="low", analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def filter_observations(history_obs,cutoff,fs):
    filtered_obs = np.zeros((len(history_obs[0]),len(history_obs)))
    for r in range(len(history_obs[0])):
        signal = [history_obs[t][r] for t in range(len(history_obs))]
        filtered_obs[r,:] = butter_lowpass_filter(signal, cutoff, fs, order=1)

    filtered_obs = filtered_obs.T.tolist()

    return filtered_obs