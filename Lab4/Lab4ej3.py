from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt

from gnpy.core.elements import *
import numpy

#Loading EDFA parameters (Optical Amplifier)
edfa_params = get_edfa_parameters("my_config.json","eqp.json")
edfa_instance = Edfa(**edfa_params)


#Instantiating WDM
with open("eqp.json","r") as read_file:
	data = json.load(read_file)


#obtaining parameter json file
n = {}
for el in data["SI"][0]:
	n[el] = data["SI"][0][el]


#setting power interval sweep for test
power_interval = numpy.arange(-2,2,0.5)
WDM_in = []
for pwr in power_interval:
	WDM_in.append(create_input_spectral_information(n["f_min"],n["f_max"], n["roll_off"], n["baud_rate"], 10**((pwr)/10),n["spacing"]))


#propagation of WDM throught EDFA
WDM_out = []
for wdm in WDM_in:
	WDM_out.append( edfa_instance.__call__(wdm) )


# Output monitor
def monitor( wdm ):
	monitor2 = Transceiver(uid="receiver")
	monitor2._calc_snr(wdm)
	osnr_ase  = monitor2.osnr_ase
	return osnr_ase


wdm_ase = [ monitor(wdm)[44] for wdm in WDM_out ]


plt.figure(1)
plt.ylabel("[dB] (OSNR)")
plt.xlabel("[dB] input power")
plt.plot(power_interval,wdm_ase,'*')
plt.show()

sig_pwr = [ 10*log10(wdm.carriers[44].power.signal/0.001) for wdm in WDM_out ]
ase_pwr = [ 10*log10(wdm.carriers[44].power.ase/0.001) for wdm in WDM_out ]

plt.figure(2)
plt.ylabel("[mW] Signal Power, ASE Power")
plt.xlabel("[dB] input power")
plt.plot(power_interval,sig_pwr,'*')
plt.plot(power_interval,ase_pwr,'+')
plt.show()