from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from math import *
from gnpy.core.elements import *


#Exercise 1

#loading fiber data from json file
with open("my_fiber.json","r") as rd_file:
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

p = 10**((data["SI"][0].pop("power_dbm")-30)/10) #dBm to dB

WDM_in = create_input_spectral_information(**data["SI"][0],power=p)

#propagate WDM_in throught Fiber

WDM_out_fiber = fiber.__call__(WDM_in)

#Calculating loss due to fiber propagation
PowerIn, PowerOut = WDM_in[1][44][4][0], WDM_out_fiber[1][44][4][0]
#Calculating Gain to restore power
FiberLoss = 10*log10(PowerOut/PowerIn)
RecoveringGain = - FiberLoss

print("Setting EDFA gain =",RecoveringGain, " we are able to recover losses")

#Instantiating EDFA
edfa_params = get_edfa_parameters("my_config.json","material/eqp.json")
edfa_params['operational']['gain_target'] = RecoveringGain

edfa = Edfa(**edfa_params)

#propagate output fiber signal throught EDFA
WDM_out_edfa = edfa.__call__(WDM_out_fiber)

PowerIn, PowerOut = WDM_out_fiber[1][44][4][0], WDM_out_edfa[1][44][4][0]
RecoveredSignal =  10*log10(PowerOut/PowerIn)

def MonitorNode( node ):
	transceiver = Transceiver(uid="receiver")
	# refto
	# https://github.com/Telecominfraproject/oopt-gnpy/blob/master/gnpy/core/elements.py
	transceiver._calc_snr(node)
	return transceiver

channel_axis = [ i for i in range(len( WDM_in[1] )) ]
input_power_vect = [ 10*log10(WDM_in[1][i][4][0])+30 for i in range(len( WDM_in[1] )) ]
fiber_power_vect = [ 10*log10(WDM_out_fiber[1][i][4][0])+30 for i in range(len( WDM_in[1] )) ]
output_power_vect = [ 10*log10(WDM_out_edfa[1][i][4][0])+30 for i in range(len( WDM_out_edfa[1] )) ]

plt.ylabel("[dBm] (power of signals)")
plt.xlabel("Channels number")
plt.plot(channel_axis,input_power_vect,'bo')
plt.plot(channel_axis,fiber_power_vect,'+')
plt.plot(channel_axis,output_power_vect,'ro')
plt.show()
