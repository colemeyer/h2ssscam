import astropy.units as u
import configparser
import os
from .data_loader import load_config_files


class Constants:
    """
    Defines physical constants and cutoff parameters for the H₂ fluorescence model.
    Values chosen per McJunkin et al. (2016), Section 4.
    """

    def __init__(self, user_config_path: str | None = None):
        self.config = self.read_config_files(user_config_path)
        self.CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
        self.ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1

        # LINE PARAMETERS
        # max vibrational (v) and rotational (J) levels for Lyman–Werner bands
        self.VMAX = self.value("VMAX")
        self.JMAX = self.value("JMAX")
        # model bandpass lambda in [1380,1620] angstroms
        self.BP_MIN = self.value("BP_MIN") * u.AA
        self.BP_MAX = self.value("BP_MAX") * u.AA
        # A_ul/A_tot threshold to include a transition
        self.LINE_STRENGTH_CUTOFF = self.value("LINE_STRENGTH_CUTOFF")
        # instrument resolving power, None = ignore instrumental broadening
        self.RESOLVING_POWER = self.value
        # plotting units; can be 'CU' or 'ERGS'
        self.UNIT = self.value("UNIT", parameter_type=str)
        # wavelength sampling
        self.DLAM = self.value("DLAM") * u.AA

        # H₂ GAS PARAMETERS
        # kinetic temperature of H2 gas
        self.TH2 = self.value("TH2") * u.K
        # total H2 column density
        self.NH2_TOT = self.value("NH2_TOT") * u.cm**-2
        # per-level column density cutoff
        self.NH2_CUTOFF = self.value("NH2_CUTOFF") * u.cm**-2
        # non-thermal Doppler b-value
        self.VELOCITY_DISPERSION = self.value("VELOCITY_DISPERSION") * u.km / u.s
        # positive is moving away from us; rho Oph has v_r = -11.4 km/s, zeta Oph has -9 km/s
        self.DOPPLER_SHIFT = self.value("DOPPLER_SHIFT") * u.km / u.s

        # HI PARAMETERS
        # kinetic temperature of HI
        self.THI = self.value("THI") * u.K
        # total HI column density
        self.NHI_TOT = self.value("NHI_TOT") * u.cm**-2
        # incident source; can be 'BLACKBODY' or 'ISRF'
        self.INC_SOURCE = self.value("INC_SOURCE")

    def read_config_files(self, user_config_path):
        if not os.path.isfile(user_config_path):
            raise ValueError("Incorrect path to user config file")
        if not isinstance(user_config_path, str):
            raise TypeError("Incorrect type of path. String please")
        self.config = load_config_files(user_config_path)

    def value(self, parameter_name, parameter_type=float):
        parameter = self.config["PARAMETERS"].get(parameter_name)
        if parameter is None:
            raise ValueError(f"Missing parameter {parameter} in config files")
        if parameter_type == float:
            return float(parameter)
        return parameter


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
