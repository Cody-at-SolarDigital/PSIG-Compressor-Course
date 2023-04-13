from numpy import array

class FlowUtils():
    """
    This class contains constant values and functions used in calculation for gas flows

    While this class can be implemented directly for calculations, it is safer to use FlowUtilsSpecific, since gas
    properties get set for a specific instance.

    ###########################################################
    References:
    - GPNO = Gas Pipeline Network Optimization, by Cody Allen
    - GPH = Gas Pipeline Hydraulics, by E. Shashi Menon
    ###########################################################
    """

    # set standard conditions and known values
    pb = 14.696                            # standard pressure [psia]
    tb = 520                               # standard temperature [degR]
    kf = 77.54                             # flow equation constant (eq 2.1 GPH, derivation in Appendix 2 of GPNO)
    m_air = 28.97                          # [lbm/lbmol]
    m_methane = 16.04                      # [lbm/lbmol]
    r_universal = 10.73                    # [PSIA*ft^3/(lbmol*degR)]

    """
    Methane (CH4) specific gas constant 
    https://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html
    **** Note that rgas_2 is the same gas constant but expressed in different physical units
    """
    rgas_1 = 0.66895  # [PSIA*ft^3/(lbm*degR)]
    rgas_2 = 96.5625  # [ft*lbf/(lbm*degR)]
    k_sp_heat_ratio = 1.32
    cmf = tb * rgas_1 / pb

    @staticmethod
    def convert_qa_to_mass_flow(q_actual: float, p_suction: float, ksuc: float, cmf: float) -> float:
        """
        convert volumetric flow to mass flow

        :param q_actual: volumetric flow                                    [acfm] = [ft^3/m]
        :param p_suction: suction pressure                                  [psia]
        :param ksuc: constant for real gas conversion, depends on
            suction temperature (see GPNO ch.2):
            ksuc = tb / (pb * t_suction * z)
        :param cmf: constant for converting to mass flow (see GPNO ch.2)
        :return: mass_flow                                                  [lbm/day]
        """
        if abs(cmf) < 1e-8:
            raise RuntimeWarning('cmf very close to 0, check physical units')
        return ksuc * p_suction * q_actual * 60 * 24 / cmf

    @staticmethod
    def convert_qb_to_mass_flow(q_standard: float, rgas_1: float = rgas_1) -> float:
        """
        convert volumetric flow at standard conditions to mass flow.  Note the units: qb must be SCFD and not MMSCFD

        :param q_standard: flow in standard conditions                      [scfd]
        :param rgas: specific gas constant                                  [psia*ft^3/(lbm*degR)]
        :return: mass flow                                                  [lbm/day]
        """
        return q_standard * FlowUtils.pb / (rgas_1 * FlowUtils.tb)

    @staticmethod
    def convert_mass_flow_to_qb(mass_flow: float, rgas_1: float = rgas_1) -> float:
        """
        convert mass flow to volumetric flow at standard conditions

        :param mass_flow: mass flow                                         [lbm/day]
        :param rgas_1: specific gas constant                                [psia*ft^3/(lbm*degR)]
        :return: q_standard                                                 [scfd] or [ft^3/day] @ standard conditions
        """
        return mass_flow * (rgas_1 * FlowUtils.tb) / FlowUtils.pb

    @staticmethod
    def convert_m_to_qa_acfm(mass_flow: float, ksuc: float, p_suction: float, cmf: float = cmf) -> float:
        """
        convert mass flow to actual volume flow in [ACFM] units

        :param mass_flow: mass flow                                         [lbm/day]
        :param ksuc: constant for real gas conversion, depends on
            suction temperature (see GPNO ch.2):
            ksuc = tb / (pb * tf * z)
        :param p_suction: suction pressure                                  [psia]
        :param cmf: constant for converting to mass flow (see GPNO ch.2)
        :return: q_actual                                                   [acfm] or [ft^3/minute]
        """
        return mass_flow * cmf / (ksuc * p_suction) / 1440

    @staticmethod
    def convert_qa_to_qb(q_actual: float, t_suction: float, p_suction: float, z_suction: float) -> float:
        """
        convert actual volume flow qa in [ACFM] to standard flow qb in [SCFM]

        Note that flow values can be slightly off from the InSight performance tag due to our approximation of Z_suction
        *** If abs(% difference) > 3% then adjustment to z_suction may need to be made

        :param q_actual: volumetric flow                                    [acfm] or [ft^3/min]
        :param t_suction: suction temperature                               [degR]
        :param p_suction: suction pressure                                  [psia]
        :param z_suction: suction compressibility factor                    [1]
        :return: q_standard                                                 [scfm] or [ft^3/min] @ standard conditions
        """
        return q_actual * p_suction * FlowUtils.tb / (t_suction * FlowUtils.pb * z_suction)

    @staticmethod
    def calc_pavg(p1: float, p2: float) -> float:
        """
        calculate average pressure in pipe segment.  See GPNO for description of "average"

        :param p1: upstream pressure                    [psia]
        :param p2: downstream pressure                  [psia]
        :return: pavg                                   [psia]
        """
        p1 = array(p1)
        p2 = array(p2)
        return (2 / 3) * (p1 + p2 - p1 * p2 / (p1 + p2))

    @staticmethod
    def calc_z_factor_cnga(sg: float, tavg: float, pavg: float) -> float:
        """
        calculates the compressibility factor z using the california natural gas association (CNGA) method
        (see GPH eq.1.34)

        :param sg: specific gravity                     [1]
        :param tf: average gas temperature              [degF]
        :param pavg: average pressure                   [psia]
        :return: z (compressibility factor)             [1]
        """
        pavg_psig = pavg - 14.7  # put pavg into gauge units
        tavg_abs = array(tavg) + 460  # put temperature into degR
        term = pavg_psig * 344400 * 10 ** (1.785 * sg) / (tavg_abs ** 3.825)
        return 1 / (1 + term)

    @staticmethod
    def general_flow_eq2_2(p1, p2, d, g, tf, l, z, f):
        """ equation 2.2 in GPH """
        return 77.54 * (FlowUtils.tb / FlowUtils.pb) * d ** 2.5 * ((p1 ** 2 - p2 ** 2) / (g * tf * l * z * f)) ** 0.5

    @staticmethod
    def _calc_ksuc(tf, z):
        """
        constant for flow equations

        :param tf: gas flowing temperature              [degR]
        :param z: compressibility factor                [1]
        :return ksuc constant                           [1/PSIA]
        """

        return FlowUtils.tb/(FlowUtils.pb*tf*z)


class FlowUtilsSpecific(FlowUtils):
    """
    Specific instance for changing assumed fuel properties.  Overides functions in FlowUtils to input specific fuel
    properties once instantiated
    """

    def __init__(self, sg: float = 0.6, k_sp_heat_ratio: float = 1.3):
        """ define new gas properties and update constants"""
        self.rgas_1 = FlowUtils.r_universal/(sg*FlowUtils.m_air)
        self.sg = sg
        self.k_sp_heat_ratio = k_sp_heat_ratio
        self.cmf = self._calc_cmf()


    def _calc_cmf(self):
        """ calculate cmf constant based on specified rgas_1 """
        return FlowUtils.tb * self.rgas_1 / FlowUtils.pb

    def convert_qa_to_mass_flow(self, q_actual: float, p_suction: float, ksuc: float) -> float:
        return super().convert_qa_to_mass_flow(q_actual, p_suction, ksuc, self.cmf)

    def convert_qb_to_mass_flow(self, q_standard: float, rgas=None) -> float:
        return super().convert_qb_to_mass_flow(q_standard, self.rgas_1)

    def convert_mass_flow_to_qb(self, mass_flow: float, rgas=None) -> float:
        return super().convert_mass_flow_to_qb(mass_flow, self.rgas_1)

    def convert_m_to_qa_acfm(self, mass_flow: float, ksuc: float, p_suction: float) -> float:
        return super().convert_m_to_qa_acfm(mass_flow, ksuc, p_suction, self.cmf)

