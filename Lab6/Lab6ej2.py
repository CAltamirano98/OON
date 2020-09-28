from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from math import *
from gnpy.core.elements import *

# Exercise 2
_span_ = 10


def MonitorNode(node):
    transceiver = Transceiver(uid="receiver")
    # refto
    # https://github.com/Telecominfraproject/oopt-gnpy/blob/master/gnpy/core/elements.py
    transceiver._calc_snr(node)
    return transceiver


def GenerateSysComponents():
    # loading fiber data from json file
    with open("my_fiber.json", "r") as rd_file:
        data = json.load(rd_file)
    # loading data fiber
    fiber = Fiber(**data["Fiber"])
    # Instantiating EDFA
    edfa_params = get_edfa_parameters("my_config.json", "eqp.json")
    edfa = Edfa(**edfa_params)

    return [fiber, edfa]  # line[0]: fiber, line[1]: edfa


def __Propagate__(wdm_input, line):
    # propagation throught fiber and edfa
    wdm_out = line[1].__call__(line[0].__call__(wdm_input))
    return wdm_out


# instantiating noiseless WDM
with open("eqp.json", "r") as read_file:
    data = json.load(read_file)

# cleaning useless values
data["SI"][0].pop("power_range_db")
data["SI"][0].pop('tx_osnr')
data["SI"][0].pop('sys_margins')

p = 10 ** ((data["SI"][0].pop("power_dbm") - 30) / 10)  # dBm to dB

WDM_in = create_input_spectral_information(**data["SI"][0], power=p)

line = [GenerateSysComponents() for i in range(_span_)]

# propagating wdm throught lines
wdm_n = [WDM_in]
for l in line:
    wdm_n.append(__Propagate__(wdm_n[line.index(l)], l))

# Monitoring propagation at every stage
monitor_n = [MonitorNode(wdm_n[i]) for i in range(len(wdm_n))]

osnr = [monitor_n[i + 1].osnr_ase[44] for i in range(_span_)]
snr = [monitor_n[i + 1].snr[44] for i in range(_span_)]
snr_nl = [monitor_n[i + 1].osnr_nli[44] for i in range(_span_)]

span = [i + 1 for i in range(_span_)]

plt.ylabel("[dB] (power of signals)")
plt.xlabel("Span step")
plt.plot(span, osnr, 'b')
plt.plot(span, snr_nl, 'g')
plt.plot(span, snr, 'ro')
plt.show()
