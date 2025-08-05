# ========== FILE: h2_model_constants.py ==========
"""
Defines physical constants and cutoff parameters for the H₂ fluorescence model.
Values chosen per McJunkin et al. (2016), Section 4.
"""
import astropy.units as u

# ------------------------------------------------------ #
# ----- LINE PARAMETERS -------------------------------- #
# ------------------------------------------------------ #

VMAX, JMAX = 14, 25  # max vibrational (v) and rotational (J) levels for Lyman–Werner bands
BP_MIN, BP_MAX = 1450 * u.AA, 1620 * u.AA  # model bandpass lambda in [1380,1620] angstroms
LINE_STRENGTH_CUTOFF = 0.01  # A_ul/A_tot threshold to include a transition
RESOLVING_POWER = 100000  # instrument resolving power, None = ignore instrumental broadening
UNIT = 'CU' # plotting units; can be 'CU' or 'ERGS'
DLAM = 0.005 * u.AA # wavelength sampling

# ------------------------------------------------------ #
# ----- H₂ GAS PARAMETERS ------------------------------ #
# ------------------------------------------------------ #
TH2 = 500 * u.K  # kinetic temperature of H2 gas
NH2_TOT = 1e20 * u.cm**-2  # total H2 column density
NH2_CUTOFF = 1e15 * u.cm**-2  # per-level column density cutoff
VELOCITY_DISPERSION = 13 * u.km / u.s  # non-thermal Doppler b-value
DOPPLER_SHIFT = 0 * u.km / u.s # positive is moving away from us; rho Oph has v_r = -11.4 km/s, zeta Oph has -9 km/s

# ------------------------------------------------------ #
# ----- HI PARAMETERS ---------------------------------- #
# ------------------------------------------------------ #

THI = 3e4 * u.K  # kinetic temperature of HI
NHI_TOT = 1e21 * u.cm**-2  # total HI column density
INC_SOURCE = 'BLACKBODY' # incident source; can be 'BLACKBODY' or 'ISRF'