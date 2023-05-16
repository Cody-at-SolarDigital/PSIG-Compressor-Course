from numpy import sqrt

class CompressorUtils:
    """
    this class contains constant values and functions used in calculation for centrifugal gas compressors

    ###########################################################
    References:
    - GPNO = Gas Pipeline Network Optimization, by Cody Allen
    - GPH = Gas Pipeline Hydraulics, by E. Shashi Menon
    ###########################################################
    """

    @staticmethod
    def comp_head(p_suction: float, p_discharge: float, z_avg: float, mratio: float, t_suction: float,
                  rgas_2: float = 96.3034) -> float:
        """
        actual/standard compressor head calculation

        :param p_suction: suction pressure                                  [psia]
        :param p_discharge: discharge pressure                              [psia]
        :param z_avg: average comrpessibility                               [1]
        :param mratio: (k-1)/k, where k = specific heat ratio               [1]
        :param t_suction: suction temperature                               [degR]
        :param rgas_2: Gas constant, default = 1545/16.043                  [ft*lbf/(lbm*degR)]
        :return: compressor Head                                            [ft*lbf/lbm]
        """
        return z_avg / mratio * t_suction * rgas_2 * ((p_discharge / p_suction) ** mratio - 1)

    @staticmethod
    def calc_comp_consumed_power(eta: float, massflow: float, head: float, mech_eff: float = 1) -> float:
        """
        Calculate the actual power [horsepower] consumed by the compressor for a specific operating point
        convert flow to [lbm/sec] and introduce conversion factor for [ft*lbf/s]
        approximated power: convert flow to [lbm/sec] = 60*60*24 and introduce conversion factor for [ft*lbf/s] = 550

        :param eta: compressor efficiency                                   [1]
        :param massflow: mass flow through compressor                       [lbm/day]
        :param head: compressor head                                        [ft*lbf/lbm]
        :param mech_eff: mechanical train efficiency                        [1]
        :return: power                                                      [horsepower or HP]
        """
        return 1 / (mech_eff * eta * 86400 * 550) * massflow * head
