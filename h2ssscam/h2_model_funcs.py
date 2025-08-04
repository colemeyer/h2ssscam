# ========== FILE: h2_model_funcs.py ==========
"""
Functional routines for H₂ fluorescence model.
Each function implements equations from McJunkin et al. (2016), Section 4.
"""
from scipy.special import wofz
import astropy.units as u
import astropy.constants as c
from h2_model_params import *
import numpy as np

CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1

# ------------------------------------------------------ #
# ----- FUNCTION DEFINITIONS --------------------------- #
# ------------------------------------------------------ #

def calc_dv(T=TH2, b=VELOCITY_DISPERSION, R=RESOLVING_POWER, instr=False):
    """
    Compute total Doppler width dv: thermal + non-thermal [+ instrumental].

    Parameters
    ----------
    T : astropy.units.Quantity
        Kinetic temperature (thermal component).
    b : astropy.units.Quantity
        Non-thermal Doppler b-value.
    R : float or None
        Instrument resolving power.
    instr : bool
        If True and R provided, include instrumental broadening.

    Returns
    -------
    dv : astropy.units.Quantity
        Combined Doppler width (same units as c.c).
    """
    dv_therm = np.sqrt(2 * c.k_B * T / (2 * c.m_p)) # Thermal broadening
    dv_nontherm = b # Non-thermal broadening
    if instr and R: # Instrumental broadening
        dv_instr = c.c / (R * np.sqrt(8 * np.log(2)))
        return np.sqrt(dv_therm**2 + dv_nontherm**2 + dv_instr**2)
    return np.sqrt(dv_therm**2 + dv_nontherm**2)


def boltzmann(Ntot, ju, jl, lam, T):
    """
    Partitioning via Boltzmann distribution.

    Parameters
    ----------
    Ntot : astropy.units.Quantity
        Total column density.
    ju, jl : int or array-like
        Upper and lower rotational quantum numbers.
    lam : astropy.units.Quantity
        Transition wavelength.
    T : astropy.units.Quantity
        Gas temperature.

    Returns
    -------
    N_J : astropy.units.Quantity
        Column density in lower level (cm^-2).
    """
    gu, gl = 2 * ju**2, 2 * jl**2
    pop = Ntot * (gu/gl) * np.exp(-(c.h * c.c / (c.k_B * lam * T)).decompose())
    return pop.to(u.cm**-2)


def calc_flu(ju, jl, lamlu, Aul):
    """
    Calculate oscillator strength f_lu from Einstein A coefficient.

    Parameters
    ----------
    ju, jl : array-like
        Upper and lower rotational levels.
    lamlu : astropy.units.Quantity
        Transition wavelengths.
    Aul : astropy.units.Quantity
        Einstein A coefficients.

    Returns
    -------
    f_ul : array
        Oscillator strengths.

    Notes
    -----
    - Implements Eq. 1 (McJunkin et al. 2016).
    - Original code includes multiplicity (2s + 1); McJunkin et al. 2016 does not.
    """
    gu = 2 * ju + 1
    gl = 2 * jl + 1
    f = c.m_e * c.c / (8 * (np.pi * c.e.esu)**2) * (gu/gl) * lamlu**2 * Aul
    return f.decompose()


def voigt(lam, lam0, gam, dv):
    """
    Compute Voigt profile H(a,y).

    Parameters
    ----------
    lam : astropy.units.Quantity
        Wavelength grid.
    lam0 : astropy.units.Quantity
        Line center wavelength.
    gam : astropy.units.Quantity
        Damping constant Γ.
    dv : astropy.units.Quantity
        Doppler width.

    Returns
    -------
    H : array
        Voigt profile values.

    Notes
    -----
    Implements Eqs. 5-7 (McJunkin et al. 2016).
    """
    nu = c.c / lam
    nu0 = c.c / lam0
    dnu = dv * nu / c.c
    a = gam / (4 * np.pi * dnu)
    y = np.abs(nu - nu0) / dnu
    return np.real(wofz(y + 1j * a))


def calc_siglu(lam, lamlu, Atot, dv, flu):
    """
    Compute absorption cross-section sigma_lu(lambda).

    Parameters
    ----------
    lam : astropy.units.Quantity
        Wavelength grid.
    lamlu : astropy.units.Quantity
        Line center wavelengths.
    Atot : astropy.units.Quantity
        Damping constants.
    dv : astropy.units.Quantity
        Doppler width.
    flu : array
        Oscillator strengths f_lu.

    Returns
    -------
    siglu : astropy.units.Quantity
        sigma_lu(lambda) array in cm^2.

    Notes
    -----
    Implements Eq. 4 (McJunkin et al. 2016).
    """
    siglu = np.zeros((len(lamlu), len(lam))) * u.cm**2
    for i in range(len(lamlu)):
        H_prof = voigt(lam, lamlu[i], Atot[i], dv)
        siglu[i,:] = (np.sqrt(np.pi) * c.e.esu**2 / (c.m_e * c.c * dv)
                      * flu[i] * lamlu[i] * H_prof).to(u.cm**2)
    return siglu


def calc_e(v, j):
    """
    Compute ro-vibrational energy E(v,J).

    Parameters
    ----------
    v : int or array-like
        Vibrational quantum number(s).
    j : int or array-like
        Rotational quantum number(s).

    Returns
    -------
    E : astropy.units.Quantity
        Energy in J (via h c k).

    Notes
    -----
    - All constants are collected and organized by Huber & Herzberg 1979, p. 250. Specific references given below.
    """
    om_e = 4400.39 * u.cm**-1 # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 and residuals could be improved
    x_e = (120.815 * u.cm**-1) / om_e # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 and residuals could be improved
    y_e = (0.7242 * u.cm**-1) / om_e # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 and residuals could be improved
    B_e = 60.864 * u.cm**-1 # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 but residuals are small beyond that too. Good enough for our purposes.
    D_e = 0.0471 * u.cm**-1 # Fink, Wiggins, Rank, JMS 18, 384 (1965).
    H_e = 4.9e-5 * u.cm**-1 # Fink, Wiggins, Rank, JMS 18, 384 (1965).
    a1, a2, a3 = 3.07638 * u.cm**-1, 0.06017 * u.cm**-1, 0.00481 * u.cm**-1 # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 but residuals are small beyond that too. Good enough for our purposes.
    b1, b2 = 0.00274 * u.cm**-1, 0.00040 * u.cm**-1 # Fink, Wiggins, Rank, JMS 18, 384 (1965).
    c1 = 0.5e-5 * u.cm**-1 # Fink, Wiggins, Rank, JMS 18, 384 (1965).

    def single(vi, ji):
        Bv = B_e - a1*(vi+0.5) + a2*(vi+0.5)**2 - a3*(vi+0.5)**3
        Dv = -D_e + b1*(vi+0.5) - b2*(vi+0.5)**2
        Hv = H_e - c1*(vi+0.5)
        Gv = om_e*(vi+0.5) - om_e*x_e*(vi+0.5)**2 + om_e*y_e*(vi+0.5)**3
        FJ = Bv*ji*(ji+1) - Dv*ji**2*(ji+1)**2 + Hv*ji**3*(ji+1)**3
        return (c.h * c.c * (Gv + FJ)).to(u.eV)

    if np.ndim(v) == 1:
        return np.array([[single(vi, jj).value for jj in j] for vi in v]) * u.eV
    else:
        return single(v, j)


def calc_nvj(ntot, T, vmax=VMAX, jmax=JMAX):
    """
    Compute level populations N_vJ for all v,J.

    Parameters
    ----------
    ntot : astropy.units.Quantity
        Total column density.
    T : astropy.units.Quantity
        Temperature.
    vmax, jmax : int
        Max vibrational and rotational levels.

    Returns
    -------
    nvj : array
        Populations per level.

    Notes
    -----
    Implements Eq. 8 (McJunkin et al. 2016).
    """
    vs = np.arange(vmax+1)
    js = np.arange(jmax+1)
    es = calc_e(vs, js)
    nvj = np.zeros((len(vs), len(js)))
    for i in range(nvj.shape[0]):
        for j in range(nvj.shape[1]):
            nvj[i,j] = np.exp(-(es[i,j] / c.k_B / T).decompose())
    nvj = ntot * nvj / np.sum(nvj, axis=None, keepdims=False)
    return nvj


def calc_tau(nvj, siglu):
    """
    Compute optical depth tau(v,J,lambda).

    Parameters
    ----------
    nvj : array
        Level populations (v,J).
    siglu : astropy.units.Quantity
        Absorption cross-sections.

    Returns
    -------
    tau : astropy.units.Quantity
        Optical depth array.

    Notes
    -----
    Implements Eq. 11 (McJunkin et al. 2016).
    """
    return (nvj[:,None] * siglu).decompose()


def calc_abs_rate(I0, tau, tau_all, unit):
    """
    Compute absorbed rate per line.

    Parameters
    ----------
    I0 : array
        Incident continuum intensity.
    tau : array
        Line optical depths.
    tau_all : array
        Total optical depth across lines.

    Returns
    -------
    abs_rate : array
        Absorption rates per line.

    Notes
    -----
    Implements Eq. 12-13 (McJunkin et al. 2016).
    """
    absr = np.zeros_like(tau) * unit
    for i in range(tau.shape[0]):
        tc = tau[i] / tau_all * tau[i]
        absr[i] = I0 * (1 - np.exp(-tc))
    return np.nan_to_num(absr)

def blackbody(lam, temp, unit):
    """
    Compute the photon radiance of a blackbody at a given temperature.

    Parameters
    ----------
    lam : astropy.units.Quantity
        Wavelength array (e.g., in Å or m).
    temp : astropy.units.Quantity
        Blackbody temperature (e.g., in K).

    Returns
    -------
    astropy.units.Quantity
        Photon radiance [ph / (cm2 s sr Å)].
    """

    # Planck spectral radiance B_lambda [W / (m2 sr m)]
    B_lambda = (2 * c.h * c.c**2 / lam**5) / np.expm1((c.h * c.c / (lam * c.k_B * temp)).decompose().value) / u.sr
    B_lambda = B_lambda.to(ERG_UNIT)
    if unit == ERG_UNIT: return B_lambda

    # convert energy radiance to photon radiance:
    N_lambda = B_lambda * (u.ph / ((c.h * c.c) / lam))
    return N_lambda.to(CU_UNIT)


def uv_continuum(lam, unit):
    """
    Empirical UV continuum function.

    Parameters
    ----------
    lam : astropy.units.Quantity
        Wavelength grid.
    G0 : float
        Interstellar UV field strength factor.

    Returns
    -------
    cont : array
        Continuum intensity.

    Notes
    -----
    Original fit comes from Draine (1978).
    """

    ### EXACT FIT FROM DRAINE (1978)
    E = c.h * c.c / lam
    F_E_unit = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.eV**-1
    a1, a2, a3 = 1.658e6 * F_E_unit, -2.152e5 * F_E_unit, 6.919e3 * F_E_unit
    F_E = a1 * (E / u.eV) + a2 * (E / u.eV)**2 + a3 * (E / u.eV)**3

    ### CONVERT TO FUNCTION OF WAVELENGTH
    if unit == CU_UNIT: return (F_E * (c.h*c.c/lam**2)).to(unit)
    else: return (F_E * (c.h*c.c/lam**2) * (c.h * c.c / lam) / u.ph).to(unit)


def dopp_shift(lam, dopp_v):
    """
    Apply non-relativistic Doppler wavelength shift.

    Parameters
    ----------
    lam : astropy.units.Quantity
        Input wavelength array.
    dopp_v : astropy.units.Quantity
        Doppler shift velocity.

    Returns
    -------
    lam_shifted : astropy.units.Quantity
        Shifted wavelength(s): lambda' = lambda * (1 + v/c).
    """
    return lam * (1 + dopp_v / c.c)


def calc_spec(lam, lamlu, Atot, dv, flux_per_trans, source, unit, dopp_v=0):
    """
    Build emergent spectrum from line profiles + continuum.

    Parameters
    ----------
    lam : astropy.units.Quantity
        Wavelength grid.
    lamlu : astropy.units.Quantity
        Line wavelengths.
    Atot : astropy.units.Quantity
        Damping constants.
    dv : astropy.units.Quantity
        Doppler width.
    flux_per_trans : array
        Flux per transition.
    source : array
        Continuum source function.

    Returns
    -------
    spec_norm : array
        Normalized emission-only spectrum.
    spec_tot_norm : array
        Normalized total spectrum including continuum.
    """
    profiles = np.zeros((len(lamlu), len(lam))) * unit
    for i in range(len(lamlu)): profiles[i,:] = flux_per_trans[i] * voigt(lam, lamlu[i], Atot[i], dv) / np.trapz(voigt(lam, lamlu[i], Atot[i], dv), lam).value

    spec = np.sum(profiles, axis=0)
    spec_tot = spec + source
    lam_shifted = dopp_shift(lam, dopp_v)
    
    return lam_shifted, spec, spec_tot