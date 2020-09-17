from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from gnpy.core.elements import *

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

WDM_instance = create_input_spectral_information(n["f_min"],n["f_max"], n["roll_off"], n["baud_rate"], 10**((n["power_dbm"]-30)/10),n["spacing"])

#Input monitor
monitor1 = Transceiver(uid="receiver")
monitor1._calc_snr(WDM_instance)
osnr_ase  = monitor1.osnr_ase

#propagation of WDM throught EDFA
propagated_WDM = edfa_instance.__call__(WDM_instance)

# Output monitor
monitor2 = Transceiver(uid="receiver")
monitor2._calc_snr(propagated_WDM)
osnr_ase2  = monitor2.osnr_ase

ase_pwr_in = [ c.power.ase for c in WDM_instance.carriers ]
ase_pwr_out = [ c.power.ase for c in propagated_WDM.carriers ]

frequency = [ ogg.frequency/(10**12) for ogg in WDM_instance[1] ]

plt.figure(1)
plt.ylabel("[W] (power of ASE)")
plt.xlabel("[THz] channels frequency")
plt.plot(frequency,ase_pwr_in,'*')
plt.plot(frequency,ase_pwr_out,'+')
plt.show()

plt.figure(2)
plt.ylabel("[dB] (OSNR)")
plt.xlabel("[THz] channels frequency")
plt.plot(frequency,osnr_ase,"o")
plt.plot(frequency,osnr_ase2,"o")
plt.show()
