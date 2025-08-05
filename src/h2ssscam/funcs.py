# # ========== FILE: h2_model_funcs.py ==========
# """
# Functional routines for H₂ fluorescence model.
# Each function implements equations from McJunkin et al. (2016), Section 4.
# """
# from scipy.special import wofz
# import astropy.units as u
# import astropy.constants as c
# from constants import *
# import numpy as np

# CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
# ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1

# # ------------------------------------------------------ #
# # ----- FUNCTION DEFINITIONS --------------------------- #
# # ------------------------------------------------------ #