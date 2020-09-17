from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from math import *
from gnpy.core.elements import Transceiver

#reading json file
with open("eqp.json","r") as read_file:
	data = json.load(read_file)

#obtaining data to compute spectral information
si = data["SI"][0]
#creating spectral information
obj = create_input_spectral_information(si["f_min"],si["f_max"],si["roll_off"],si["baud_rate"],10**((si["power_dbm"]-30)/10),si["spacing"])

#replacing value as requested
amp_power = 10**(3/10)
ase_power = 10**((-40-30)/10)
nli_power = 10**((-43-40)/10)
obj1 = obj._replace(carriers=tuple(c._replace(power = c.power._replace(signal= c.power.signal + amp_power, nli = nli_power, ase= ase_power))
                              for c in obj.carriers))

#LAB 4 - EXERCISE 1

#dati plot
p_,f_ = [],[]
ase_,nli_ = [],[]


for ogg in obj1[1]:
	f_.append(ogg.frequency/(10**12))
	p_.append(10*log10(ogg.power.signal)+30)
	ase_.append(10*log10(ogg.power.ase)+30)
	nli_.append(10*log10(ogg.power.nli)+30)

plt.ylabel("[dBm] (power of signals)")
plt.xlabel("[THz] channels frequency")
plt.plot(f_,p_,'+')
plt.plot(f_,ase_,'*')
plt.plot(f_,nli_,'o')
plt.show()


#LAB 4 - EXERCISE 2
trans = Transceiver(uid="receiver")

trans._calc_snr(obj1)

osnr_ase  = trans.osnr_ase
osnr_nli  = trans.osnr_nli
snr	  = trans.snr

frequency = [ ogg.frequency/(10**12) for ogg in obj1[1] ]

plt.ylabel("[dB] (power of signals)")
plt.xlabel("[THz] channels frequency")
plt.plot(frequency,osnr_ase,'*')
plt.plot(frequency,osnr_nli,'+')
plt.plot(frequency,snr,'o')
plt.show()
