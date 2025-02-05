import gnpy.core.info as gn
import matplotlib.pyplot as plt
import Mylib as ml
import json as js

# Exercise 3
# Create and plot a spectral information representing a comb of 120 WDM channels.
# Each channel has a -3 dB bandwidth of 32 Gbaud, a roll-off of 0.2 and a
# power equal to 0 dBm. This comb is centered around 193 THz and the spacing
# is 45 GHz.

data = {"f center": 194.5e12, "roll off": 0.2, "baud rate": 32e9, "power": 1e-3, "spacing": 45e9}

# write json data into a file
with open("Parameters3.json", "w") as write_file:
    js.dump(data, write_file, indent=4)

with open("Parameters3.json", "r") as read_file:
    data = js.load(read_file)

f_min = data["f center"]-(data["spacing"]*60)
f_max = data["f center"]+(data["spacing"]*60)

si = gn.create_input_spectral_information(f_min, f_max, data["roll off"], data["baud rate"], data["power"], data["spacing"])

ml.plot_spectrum(si)
plt.show()
