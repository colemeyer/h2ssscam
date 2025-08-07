import astropy.units as u
import os
from .data_loader import load_config_files
from datetime import datetime


class Constants:
    
    def __init__(self, user_config_path: str | None = None):
        self.config = self.read_config_files(user_config_path)
        self.CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
        self.ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1

        # LINE PARAMETERS
        # max vibrational (v) and rotational (J) levels for Lyman–Werner bands
        self.VMAX = self.value("vmax")
        self.JMAX = self.value("jmax")
        # model bandpass lambda in [1380,1620] angstroms
        self.BP_MIN = self.value("bp_min") * u.AA
        self.BP_MAX = self.value("bp_max") * u.AA
        # A_ul/A_tot threshold to include a transition
        self.LINE_STRENGTH_CUTOFF = self.value("line_strength_cutoff")
        # instrument resolving power, None = ignore instrumental broadening
        self.RESOLVING_POWER = self.value("resolving_power")
        # plotting units; can be 'CU' or 'ERGS'
        self.UNIT = self.value("unit", parameter_type=str)
        # wavelength sampling
        self.DLAM = self.value("dlam") * u.AA

        # H₂ GAS PARAMETERS
        # kinetic temperature of H2 gas
        self.TH2 = self.value("th2") * u.K
        # total H2 column density
        self.NH2_TOT = self.value("nh2_tot") * u.cm**-2
        # per-level column density cutoff
        self.NH2_CUTOFF = self.value("nh2_cutoff") * u.cm**-2
        # non-thermal Doppler b-value
        self.VELOCITY_DISPERSION = self.value("velocity_dispersion") * u.km / u.s
        # positive is moving away from us; rho Oph has v_r = -11.4 km/s, zeta Oph has -9 km/s
        self.DOPPLER_SHIFT = self.value("doppler_shift") * u.km / u.s

        # HI PARAMETERS
        # kinetic temperature of HI
        self.THI = self.value("thi") * u.K
        # total HI column density
        self.NHI_TOT = self.value("nhi_tot") * u.cm**-2
        # incident source; can be 'BLACKBODY' or 'ISRF'
        self.INC_SOURCE = self.value("inc_source", parameter_type=str)

    def value(self, parameter_name, parameter_type=float):
        """
        Load a parameter value from the configparser and transform it to float if required

        Parameters
        ----------
        parameter_name : str
            Name of a parameter
        parameter_type : type, optional
            Transform parameter value into this type; configparser stores everything as strings. Only float is implemented, and by default float

        Returns
        -------
        str | float


        Raises
        ------
        ValueError
            If parameter is not defined in the config file
        """
        parameter = self.config["PARAMETERS"].get(parameter_name)
        if parameter is None:
            raise ValueError(f"Missing parameter {parameter} in config files")
        if parameter_type == float:
            return float(parameter)
        return parameter

    def _set_value(self, parameter_name, value):
        """
        Updates value in configparser

        Parameters
        ----------
        parameter_name : str
            Parameter name to be updated
        value : float | Quantity | str
            Value of the parameter to be updated
        """        
        if type(value) == u.Quantity:
            value = value.value
        if type(value) != str:
            value = str(value)
        self.config["PARAMETERS"][parameter_name] = value

    def read_config_files(self, user_config_path):
        """
        Loads default config file from the package and user's config file if specified

        Parameters
        ----------
        user_config_path : str
            Path to the config file

        Returns
        -------
        configparser
        """        
        config = load_config_files()
        if user_config_path is None:
            return config
        if not os.path.isfile(user_config_path):
            raise ValueError("Incorrect path to user config file")
        if not isinstance(user_config_path, str):
            raise TypeError("Incorrect type of path. String please")
        config.read(user_config_path)

        return config

    def save_config_file(self, output_path):
        """
        Save config file with currently used parameters value

        Parameters
        ----------
        output_path : str
            Output file path
        """
        for key in self.config["PARAMETERS"]:
            self._set_value(key, getattr(self, key.upper()))
        timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
        output_file = os.path.join(output_path, f"config_{timestamp}.ini")

        with open(output_file, "w") as configfile:
            self.config.write(configfile)
