import json as js
import gnpy.core.info as gn
import matplotlib.pyplot as plt
import Mylib as ml

# 1. Create a json file with the following parameters (in fundamental units, Hz,
# Baud, W):

data = {"f min": 191.5e12, "f max": 194.5e12, "roll off": 0.2, "baud rate": 32e9, "power": 1e-3, "spacing": 50e9}

# write json data into a file
with open("Parameters.json", "w") as write_file:
    js.dump(data, write_file, indent=4)

# read json data from a file
with open("Parameters.json", "r") as read_file:
    data = js.load(read_file)

# Create spetral information
si = gn.create_input_spectral_information(data['f min'], data['f max'], data['roll off'], data['baud rate'],
                                          data['power'], data['spacing'])

ml.plot_spectrum(si)
plt.show()