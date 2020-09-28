from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from gnpy.core.elements import *
import numpy

# Exercise 2
_span_ = 10

power_interval = numpy.arange(-5.0, 2.0, 0.25)


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
data["SI"][0].pop("power_dbm")

# Generating spectral information with a power sweep
WDM_in = []
for power in power_interval:
    p = 10 ** ((power - 30) / 10)  # dBm to W
    WDM_in.append(create_input_spectral_information(**data["SI"][0], power=p))

# Generating components of the line
line = [GenerateSysComponents() for i in range(_span_)]

# Propagating each WDM with different power
WDM_out = []
for wdm in WDM_in:
    wdm_n = [WDM_in[WDM_in.index(wdm)]]
    # propagating wdm throught lines
    for l in line:
        wdm_n.append(__Propagate__(wdm_n[line.index(l)], l))
    WDM_out.append(wdm_n)

# Monitoring propagation at the end of every span
monitor_n = [MonitorNode(WDM_out[i][-1]) for i in range(len(WDM_out))]

# Plot the results compared with power sweep
plt.ylabel("[dB] (SNR_{NL}, SNR_{ASE}, GSNR)")
plt.xlabel("[dBm] (power sweep of signal)")
plt.plot(power_interval, [monitor_n[i].osnr_ase[22] for i in range(len(monitor_n))], 'b+')
plt.plot(power_interval, [monitor_n[i].snr[22] for i in range(len(monitor_n))], 'o')
plt.plot(power_interval, [monitor_n[i].osnr_nli[22] for i in range(len(monitor_n))], 'g+')
plt.show()
