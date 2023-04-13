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

    @staticmethod
    def get_comp_map_eff_and_speed(h: float, f: float, comp: str = 'c45') -> float:
        """ new method to get compressor speed and efficiency

        :param h: head
        :param f: flow
        :param comp: compressor_model
        :return: (speed,eta) where speed=RPM=NPT, and eta=Eta=istentropic efficiency
        """
        print(f'****************  comp_map_eff{comp}***********************')
        # todo: these values should be loaded via json
        if comp == 'c40':
            speed_coeff = []
            eta_coeff = []
        elif comp == 'c45':
            speed_coeff = []
            eta_coeff = []
        elif comp == 'c65':
            # model coefficients
            speed_coeff = [ 0.,  0.53067103,  0.69780107,  0.10331809, -0.16200128,
               -0.08223556, -0.02178603,  0.02439215,  0.02240949,  0.01975992]
            speed_intercept = 0.016134302866640093
            eta_coeff = [  0.0,   1.64671832,   3.54276976,  -4.46497957,
               -15.45505105, -29.54862438,  45.84447783, -13.97743804,
                41.571843  , -32.512664  ,  -6.20509292,  -3.85618839,
                12.36946311,   2.73892889,  -2.69437958,  -1.52242883,
                 1.56210671,  -6.50088177,  10.14315544,  -5.40818861]
            eta_intercept = 0.8375544267313313

            # standardization parameters: mu=mean, var=variance
            mu = [6.92714286e+03, 8.20365714e+03, 7.72007143e+03, 7.69570000e-01]  # [rpm, Q, Head, Eta]
            var = [3.61042041e+06, 1.22862924e+07, 2.19731243e+07, 1.12604241e-02]  # same order as above
        elif comp == 'c75':
            # model coefficients
            speed_coeff = [ 0.0,  0.62712991,  0.76464215,  0.05535437, -0.15509031,
                -0.08349052, -0.00822576,  0.01047928,  0.01391111,  0.02049904]
            speed_intercept = 0.015590493104049169
            eta_coeff = [  0.0,  -2.47666921,  -1.45398106,   2.84513311,
               -11.7445615 , -23.74548426,  30.25531553, -13.05048949,
                31.45447672, -19.00760562,  -3.19451139,  -7.27199042,
                10.2292593 ,  -8.17331058,  21.34364545, -12.270367  ,
                -4.21022133,  14.88791052, -16.23265985,   5.47532771]
            eta_intercept = 0.4297928285078984

            # standardization parameters: mu=mean, var=variance
            mu = [5.50000000e+03, 1.08585057e+04, 1.49221000e+04, 8.18492857e-01]
            var = [1.00000000e+06, 1.11464157e+07, 5.98813624e+07, 3.81112209e-03]

        # transform flow and head
        transformed_f = (f - mu[1])/sqrt(var[1])
        transformed_h = (h - mu[2])/sqrt(var[2])

        # calculate speed and eta in transformed coordinates
        speed = CompressorUtils.poly_3_2var(speed_coeff, speed_intercept, transformed_f, transformed_h)
        eta = CompressorUtils.poly_3_3var(eta_coeff, eta_intercept, transformed_f, transformed_h, speed)

        # inverse transform back to original coordinates
        predicted_speed = speed*sqrt(var[0]) + mu[0]
        predicted_eta = eta*sqrt(var[3]) + mu[3]
        return predicted_speed, predicted_eta