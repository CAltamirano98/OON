from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from math import *
from gnpy.core.elements import *


#Exercise 1
#loading fiber data from json file
with open("my_fiber.json" ,"r") as rd_file:
    data = json.load(rd_file)

#loading data fiber
fiber = Fiber(**data["Fiber"])


#instantiating noiseless WDM
with open("eqp.json","r") as read_file:
	data = json.load(read_file)

#cleaning useless values
data["SI"][0].pop("power_range_db")
data["SI"][0].pop('tx_osnr')
data["SI"][0].pop('sys_margins')

p = 10**((data["SI"][0].pop("power_dbm")-30)/10) #dbm to db

WDM_in = create_input_spectral_information(**data["SI"][0],power=p)

#propagate WDM_in throught Fiber

WDM_out = fiber.__call__(WDM_in)

#applying monitor at input and output fiber
mnt1 = Transceiver(uid="receiver")
mnt1._calc_snr(WDM_in)

mnt2 = Transceiver(uid="receiver")
mnt2._calc_snr(WDM_out)

nli_out = [c.power.signal/c.power.nli for c in WDM_out.carriers]#mnt2.osnr_nli

#generating frequency axis
frequency = [ ogg.frequency/(10**12) for ogg in WDM_out[1] ]

plt.ylabel("[dBm] (power of signals)")
plt.xlabel("[THz] channels frequency")
plt.plot(frequency,nli_out,'o')
plt.show()
