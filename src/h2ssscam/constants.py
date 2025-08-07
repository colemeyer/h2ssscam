import astropy.units as u
import configparser
class Constants:
    """
    Defines physical constants and cutoff parameters for the H₂ fluorescence model.
    Values chosen per McJunkin et al. (2016), Section 4.
    """
    def __init__(self,config_file : str | None= None):
        self.config = self.read_config_files(config_file)
        self.CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
        self.ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1
        # max vibrational (v) and rotational (J) levels for Lyman–Werner bands
        self.VMAX  = float(self.value('VMAX'))
        self.JMAX = float(self.value('JMAX'))
        # model bandpass lambda in [1380,1620] angstroms
        self.BP_MIN = float(self.value('BP_MIN'))* u.AA
        self.BP_MAX = float(self.value('BP_MAX'))* u.AA

        self.LINE_STRENGTH_CUTOFF = float(self.value('LINE_STRENGTH_CUTOFF'))
        self.UNIT = self.value('UNIT')
        self.DLAM = float(self.value('DLAM'))* u.AA
        self.TH2 = float(self.value('TH2'))* u.K 
        self.NH2_TOT = float(self.value('NH2_TOT'))* u.cm**-2 
        self.NH2_CUTOFF = float(self.value('NH2_CUTOFF'))* u.cm**-2 
        self.VELOCITY_DISPERSION = float(self.value('VELOCITY_DISPERSION'))*u.km / u.s 
        self.DOPPLER_SHIFT = float(self.value('DOPPLER_SHIFT'))* u.km / u.s
        self.THI = float(self.value('THI'))* u.K
        self.NHI_TOT = float(self.value('NHI_TOT'))* u.cm**-2
        self.INC_SOURCE = self.value('INC_SOURCE')
    def read_config_files(self,user_config_file):
        pass
    def value(self,parameter_name)
        return self.config['USER'].get(parameter_name, self.config['DEFAULT'][parameter_name])
# ========== FILE: h2_model_constants.py ==========


# ------------------------------------------------------ #
# ----- LINE PARAMETERS -------------------------------- #
# ------------------------------------------------------ #

VMAX, JMAX = 14, 25  # max vibrational (v) and rotational (J) levels for Lyman–Werner bands
BP_MIN, BP_MAX = 1450 * u.AA, 1620 * u.AA  # model bandpass lambda in [1380,1620] angstroms
LINE_STRENGTH_CUTOFF = 0.01  # A_ul/A_tot threshold to include a transition
RESOLVING_POWER = 100000  # instrument resolving power, None = ignore instrumental broadening
UNIT = "CU"  # plotting units; can be 'CU' or 'ERGS'
DLAM = 0.005 * u.AA  # wavelength sampling

# ------------------------------------------------------ #
# ----- H₂ GAS PARAMETERS ------------------------------ #
# ------------------------------------------------------ #
TH2 = 500 * u.K  # kinetic temperature of H2 gas
NH2_TOT = 1e20 * u.cm**-2  # total H2 column density
NH2_CUTOFF = 1e15 * u.cm**-2  # per-level column density cutoff
VELOCITY_DISPERSION = 13 * u.km / u.s  # non-thermal Doppler b-value
DOPPLER_SHIFT = 0 * u.km / u.s  # positive is moving away from us; rho Oph has v_r = -11.4 km/s, zeta Oph has -9 km/s

# ------------------------------------------------------ #
# ----- HI PARAMETERS ---------------------------------- #
# ------------------------------------------------------ #

THI = 3e4 * u.K  # kinetic temperature of HI
NHI_TOT = 1e21 * u.cm**-2  # total HI column density
INC_SOURCE = "BLACKBODY"  # incident source; can be 'BLACKBODY' or 'ISRF'
