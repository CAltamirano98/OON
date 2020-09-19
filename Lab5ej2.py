from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from math import *
from gnpy.core.elements import *


# Exercise 2

def ReturnPowerWDM(wdm):
    pwr_sig = []
    pwr_ase = []
    pwr_nli = []

    for carrier in WDM_in.carriers:
        pwr_sig.append((carrier.power.signal / 0.001))
        pwr_ase.append((carrier.power.ase / 0.001))
        pwr_nli.append((carrier.power.nli / 0.001))

    return pwr_sig, pwr_ase, pwr_nli


def Monitor(wdm):
    n = Transceiver(uid="receiver")
    n._calc_snr(wdm)
    return n


# loading fiber data from json file
with open("my_fiber.json", "r") as rd_file:
    data = json.load(rd_file)

# loading data fiber
fiber = Fiber(**data["Fiber"])

# instantiating noiseless WDM
with open("eqp.json", "r") as read_file:
    data = json.load(read_file)

# cleaning useless values
data["SI"][0].pop("power_range_db")
data["SI"][0].pop('tx_osnr')
data["SI"][0].pop('sys_margins')

p = 10 ** ((data["SI"][0].pop("power_dbm") - 30) / 10)  # dbm to db

WDM_in = create_input_spectral_information(**data["SI"][0], power=p)

# propagate WDM_in throught Fiber

WDM_out_fiber = fiber.__call__(WDM_in)

# instantiating EDFA
# Loading EDFA parameters (Optical Amplifier)
edfa_params = get_edfa_parameters("my_config.json", "eqp.json")
edfa = Edfa(**edfa_params)

# propagate WDM throught EDFA
WDM_out_edfa = edfa.__call__(WDM_out_fiber)

monitor_2 = Monitor(WDM_out_fiber)
monitor_3 = Monitor(WDM_out_edfa)

ASE = [monitor_2.osnr_ase[44], monitor_3.osnr_ase[44]]
NLI = [monitor_2.osnr_nli[44], monitor_3.osnr_nli[44]]
GSNR = [monitor_2.snr[44], monitor_3.snr[44]]

wdm1_pwr, wdm1_ase, wdm1_nli = ReturnPowerWDM(WDM_in)
wdm2_pwr, wdm2_ase, wdm2_nli = ReturnPowerWDM(WDM_out_fiber)
wdm3_pwr, wdm3_ase, wdm3_nli = ReturnPowerWDM(WDM_out_edfa)
frequency = [ogg.frequency / (10 ** 12) for ogg in WDM_out_edfa[1]]

plt.figure(1)
plt.title("Signal Power, NLI Power, ASE Power")
plt.ylabel("[mW] (power of signals)")
plt.xlabel("[THz] channels frequency")
plt.plot(frequency, wdm1_pwr, 'ro', frequency, wdm1_ase, 'go', frequency, wdm1_nli, 'bo', \
         frequency, wdm2_pwr, 'r', frequency, wdm2_ase, 'g+', frequency, wdm2_nli, 'b+', \
         frequency, wdm3_pwr, 'r*', frequency, wdm3_ase, 'g*', frequency, wdm3_nli, 'b*')
plt.show()

plt.figure(2)
plt.title("GSNR, SNR_NL, OSNR Power")
plt.ylabel("[dB] (power of signals)")
plt.xlabel("[Span element]")
plt.plot([1, 2], GSNR, '*')
plt.plot([1, 2], NLI, '+')
plt.plot([1, 2], ASE, 'o')
plt.show()
