from utilities import *
from gnpy.core.info import *
import json
import matplotlib.pyplot as plt
from math import *
from gnpy.core.elements import *
import numpy
import scipy.constants as sp
from Lab7lib import Ex3

ex3 = Ex3()

__span__ = 10
power_interval = numpy.arange(-5.0,2.0,0.25)

ex3.GenerateLineSystem(__span__)

WDM_out_ex1 = ex3.Ex1_WDM_propagate( power_interval )
WDM_out_ex2,opt_power = ex3.Ex2_WDM_propagate()


# Monitoring propagation at the end of every span
monitor_1 = [ ex3.MonitorNode(WDM_out_ex1[i][-1]) for i in range(len(WDM_out_ex1)) ]
monitor_2 = ex3.MonitorNode(WDM_out_ex2[-1])

#Plot the results compared with power sweep
plt.ylabel("[dB] (SNR_{NL}, SNR_{ASE}, GSNR)")
plt.xlabel("[dBm] (power sweep of signal)")
plt.plot(power_interval, [monitor_1[i].osnr_ase[44] for i in range(len(monitor_1))],'b+')
plt.plot(power_interval, [monitor_1[i].snr[44] for i in range(len(monitor_1))],'+')
plt.plot(power_interval, [monitor_1[i].osnr_nli[44] for i in range(len(monitor_1))],'g+')
plt.plot(opt_power, monitor_2.snr[44], 'o')
plt.plot(opt_power, monitor_2.osnr_ase[44], 'bo')
plt.plot(opt_power, monitor_2.osnr_nli[44], 'go')
plt.show()
