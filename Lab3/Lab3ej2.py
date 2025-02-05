import matplotlib.pyplot as plt
import numpy as np

# Exercise 2
# Plot the spectral information

def plot_spectrum(si):
    channel = si.carriers

    freq = np.zeros(len(channel))
    power = np.zeros(len(channel))

    for i in range(0, len(channel)):
        ch_i = si.carriers[i]
        freq[i] = ch_i.frequency/1e12
        power[i] = 10*(np.log10(ch_i.power.signal/0.001))

    plt.plot(freq, power, '.r', label='line 1', linewidth=2)
    plt.ylabel('Power [dBm]')
    plt.xlabel('frequency [THz]')
    plt.show()
