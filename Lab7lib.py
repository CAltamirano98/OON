from utilities import *
from gnpy.core.info import *
import json
from math import *
from gnpy.core.elements import *
import scipy.constants as sp


class Ex3:
    line = None

    def LOGO(self, WDM_data, SyS_Line):
        # ref to slides OSL controller.pdf @ page: 14
        # parameter definition
        Rs = WDM_data["baud_rate"]
        Ks = WDM_data["spacing"] / Rs
        B_opt = WDM_data["f_max"] - WDM_data["f_min"]
        f0 = (WDM_data["f_max"] + WDM_data["f_min"]) / 2

        P_base = sp.h * Rs * f0
        F_0 = 5.5  # obtained from edfa.nf
        ISNRab = 0
        power_v = []
        for node in SyS_Line:
            alpha_i = node[0].lin_loss_coef
            gamma_i = node[0].gamma
            L_eff_i = node[0].effective_length
            beta2_i = node[0].beta2()

            # non-linear efficiency
            arg_i = (0.5 * pi ** 2) * (B_opt / (Ks * Rs)) ** (2 / Ks) * (abs(beta2_i) * Rs * Rs) / alpha_i
            eta_nli_i = 16 / (27 * pi) * log(arg_i) * (alpha_i * (gamma_i * L_eff_i) ** 2) / ((beta2_i) * Rs * Rs)

            F_i = 5.5  # obtained from edfa.nf
            L_i = node[0].fiber_loss
            FLP = (F_i * L_i * P_base) ** 2
            ISNRab += (2 / 3) * (2 * eta_nli_i * FLP) ** (1 / 3)
            # optimal power for each span in line system
            power_v.append(10 * log10(((F_i * L_i * P_base / (2 * eta_nli_i)) ** (1 / 3)) / 0.001))

        return 10 * log10(1 / ISNRab), power_v[0]

    def MonitorNode(self, node):
        transceiver = Transceiver(uid="receiver")
        # refto
        # https://github.com/Telecominfraproject/oopt-gnpy/blob/master/gnpy/core/elements.py
        transceiver._calc_snr(node)
        return transceiver

    def GenerateSysComponents(self):
        # loading fiber data from json file
        with open("my_fiber_ex3.json", "r") as rd_file:
            data = json.load(rd_file)
        # loading data fiber
        fiber = Fiber(**data["Fiber"])
        # Instantiating EDFA
        edfa_params = get_edfa_parameters("my_config.json", "eqp.json")
        edfa = Edfa(**edfa_params)

        return [fiber, edfa]  # line[0]: fiber, line[1]: edfa

    def __Propagate__(self, wdm_input, line):
        # propagation throught fiber and edfa
        wdm_out = line[1].__call__(line[0].__call__(wdm_input))
        return wdm_out

    def InitWDM(self):
        # instantiating noiseless WDM
        with open("eqp.json", "r") as read_file:
            data = json.load(read_file)

        # cleaning notused values
        data["SI"][0].pop("power_range_db")
        data["SI"][0].pop('tx_osnr')
        data["SI"][0].pop('sys_margins')
        data["SI"][0].pop("power_dbm")
        return data

    def GenerateLineSystem(self, __span__):
        self.line = [self.GenerateSysComponents() for i in range(__span__)]

    def Ex1_WDM_propagate(self, power_sweep):
        WDM_in = []
        for power in power_sweep:
            p = 10 ** ((power - 30) / 10)  # dBm to W
            WDM_in.append(create_input_spectral_information(**self.InitWDM()["SI"][0], power=p))
        WDM_out = []
        for wdm in WDM_in:
            wdm_n = [WDM_in[WDM_in.index(wdm)]]
            # propagating wdm throught lines
            for l in self.line:
                wdm_n.append(self.__Propagate__(wdm_n[self.line.index(l)], l))
            WDM_out.append(wdm_n)
        return WDM_out

    def Ex2_WDM_propagate(self):
        snr, OI_power = self.LOGO(self.InitWDM()["SI"][0], self.line)
        p = 10 ** ((OI_power - 30) / 10)  # dBm to W
        WDM_in = create_input_spectral_information(**self.InitWDM()["SI"][0], power=p)
        WDM_out = [WDM_in]
        # Propagating WDM throught lines
        for l in self.line:
            WDM_out.append(self.__Propagate__(WDM_out[self.line.index(l)], l))
        return WDM_out, OI_power
