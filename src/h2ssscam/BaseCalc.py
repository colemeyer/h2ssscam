from dataclasses import dataclass
from astropy.units import Quantity
from scipy.special import wofz
import astropy.constants as c
import astropy.units as u
import numpy as np
from .Constants import Constants


@dataclass
class BaseCalc:
    constant: Constants
    _dv_phys: Quantity = None
    _dv_tot: Quantity = None
    _tau: Quantity = None
    _tau_tot: Quantity = None
    _siglu: Quantity = None

    @property
    def dv_phys(self):
        if self._dv_phys is None:
            self._dv_phys = self._calc_dv()  # thermal + non-thermal (Eq. 7)
        return self._dv_phys

    @property
    def dv_tot(self, instr=True):
        if self._dv_tot is None:
            self._dv_tot = self._calc_dv(instr=instr)  # include instrumental
        return self._dv_tot

    def siglu(self, lam=None, hih2_lamlu=None, hih2_Atot=None, hih2_flu=None):
        if self._siglu is None:
            if lam is None or hih2_lamlu is None or hih2_Atot is None or hih2_flu is None:
                raise ValueError("Here's should be an offensive message. To be implemented...")
            self._siglu = self._calc_siglu(lam, hih2_lamlu, hih2_Atot, self.dv_phys, hih2_flu)
        return self._siglu

    def tau(self, hih2_N=None):

        if self._tau is None:
            if self._siglu is None or hih2_N is None:
                raise ValueError("Calculate siglu before...")
            self._tau = self._calc_tau(hih2_N, self._siglu)
        return self._tau

    @property
    def tau_tot(self):
        if self._tau_tot is None:
            if self._tau is None:
                raise ValueError("Calculate tau_tot before...")
            self._tau_tot = self._calc_tau_tot()
        return self._tau_tot
    def reset_parameters(self):
        """Clean calculated and stored values
        """        
        self._dv_phys = None
        self._dv_tot = None
        self._tau = None
        self._tau_tot = None
        self._siglu = None

    def calc_flu(self, ju, jl, lamlu, Aul):
        """Calculate oscillator strength f_lu from Einstein A coefficient.

        Parameters
        ----------
        ju :array-like
            Upper rotation levels
        jl :array-like
            Lower rotation levels
        lamlu : astropy.units.Quantity
            Transition wavelengths.
        Aul : astropy.units.Quantity
            Einstein A coefficients.

        Returns
        -------
        array
            Oscillator strengths.

        Notes
        -----
            - Implements Eq. 1 (McJunkin et al. 2016).
            - Original code includes multiplicity (2s + 1); McJunkin et al. 2016 does not.
        """

        gu = 2 * ju + 1
        gl = 2 * jl + 1
        f = c.m_e * c.c / (8 * (np.pi * c.e.esu) ** 2) * (gu / gl) * lamlu**2 * Aul
        return f.decompose()

    def calc_nvj(self, ntot, T):
        """Compute level populations N_vJ for all v,J.

        Parameters
        ----------
        ntot : astropy.units.Quantity
            Total column density.
        T : astropy.units.Quantity
            Temperature.
        vmax : int, optional
            Max vibrational and rotational levels., by default VMAX
        jmax : int, optional
            Max vibrational and rotational levels., by default JMAX

        Returns
        -------
        array
            Populations per level.

        Notes
        -----
        Implements Eq. 8 (McJunkin et al. 2016).
        """

        vs = np.arange(self.constant.VMAX + 1)
        js = np.arange(self.constant.JMAX + 1)
        es = self._calc_e(vs, js)
        nvj = np.zeros((len(vs), len(js)))
        for i in range(nvj.shape[0]):
            for j in range(nvj.shape[1]):
                nvj[i, j] = np.exp(-(es[i, j] / c.k_B / T).decompose())
        nvj = ntot * nvj / np.sum(nvj, axis=None, keepdims=False)
        return nvj

    def boltzmann(self, Ntot, ju, jl, lam, T):
        """Partitioning via Boltzmann distribution.

        Parameters
        ----------
        Ntot : astropy.units.Quantity
            Total column density.
        ju : int
            Upper rotational quantum numbers.
        jl : int
            Lower rotational quantum numbers.
        lam : astropy.units.Quantity
            Transition wavelength.
        T : astropy.units.Quantity
            Gas temperature.

        Returns
        -------
        astropy.units.Quantity
            Column density in lower level (cm^-2).
        """

        gu, gl = 2 * ju**2, 2 * jl**2
        pop = Ntot * (gu / gl) * np.exp(-(c.h * c.c / (c.k_B * lam * T)).decompose())
        return pop.to(u.cm**-2)

    def blackbody(self, lam, temp, unit):
        """Compute the photon radiance of a blackbody at a given temperature.

        Parameters
        ----------
        lam : astropy.units.Quantity
            Wavelength array (e.g., in Å or m).
        temp : astropy.units.Quantity
            Blackbody temperature (e.g., in K).
        unit : astropy.units.Quantity
            Continuum units or CGS units.

        Returns
        -------
        astropy.units.Quantity
            Photon radiance [ph / (cm2 s sr Å)].
        """

        # Planck spectral radiance B_lambda [W / (m2 sr m)]
        B_lambda = (2 * c.h * c.c**2 / lam**5) / np.expm1((c.h * c.c / (lam * c.k_B * temp)).decompose().value) / u.sr
        B_lambda = B_lambda.to(self.constant.ERG_UNIT)
        if unit == self.constant.ERG_UNIT:
            return B_lambda

        # convert energy radiance to photon radiance:
        N_lambda = B_lambda * (u.ph / ((c.h * c.c) / lam))
        return N_lambda.to(self.constant.CU_UNIT)

    def uv_continuum(self, lam, unit):
        """Empirical UV continuum function.

        Parameters
        ----------
        lam : astropy.units.Quantity
            Wavelength grid.
        unit : astropy.units.Quantity
            Continuum units or CGS units.

        Returns
        -------
        array
            Continuum intensity.

        Notes
        -----
        Original fit comes from Draine (1978).
        """

        ### EXACT FIT FROM DRAINE (1978)
        E = c.h * c.c / lam
        F_E_unit = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.eV**-1
        a1, a2, a3 = 1.658e6 * F_E_unit, -2.152e5 * F_E_unit, 6.919e3 * F_E_unit
        F_E = a1 * (E / u.eV) + a2 * (E / u.eV) ** 2 + a3 * (E / u.eV) ** 3

        ### CONVERT TO FUNCTION OF WAVELENGTH
        if unit == self.constant.CU_UNIT:
            return (F_E * (c.h * c.c / lam**2)).to(unit)
        else:
            return (F_E * (c.h * c.c / lam**2) * (c.h * c.c / lam) / u.ph).to(unit)

    def calc_abs_rate(self, I0, tau, tau_all, unit):
        """Compute absorbed rate per line.

        Parameters
        ----------
        I0 : array
            Incident continuum intensity.
        tau : array
            Line optical depths.
        tau_all : array
            Total optical depth across lines.
        unit : astropy.units.Quantity
            Continuum units or CGS units.

        Returns
        -------
        array
            Absorption rates per line.

        Notes
        ------
        Implements Eq. 12-13 (McJunkin et al. 2016).
        """

        absr = np.zeros_like(tau) * unit
        for i in range(tau.shape[0]):
            tc = tau[i] / tau_all * tau[i]
            absr[i] = I0 * (1 - np.exp(-tc))
        return np.nan_to_num(absr)

    def calc_spec(self, lam, lamlu, Atot, dv, flux_per_trans, source, unit, dopp_v=0):
        """Build emergent spectrum from line profiles + continuum.

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
        unit : astropy.units.Quantity
            Continuum units or CGS units.
        dopp_v : astropy.units.Quantity, optional
            Doppler Velocity.

        Returns
        -------
        array
            Wavelength grid.

        array
            Normalized emission-only spectrum.

        array
            Normalized total spectrum including continuum.
        """

        profiles = np.zeros((len(lamlu), len(lam))) * unit
        for i in range(len(lamlu)):
            profiles[i, :] = (
                flux_per_trans[i]
                * self._voigt(lam, lamlu[i], Atot[i], dv)
                / np.trapezoid(self._voigt(lam, lamlu[i], Atot[i], dv), lam).value
            )

        spec = np.sum(profiles, axis=0)
        spec_tot = spec + source
        lam_shifted = self._dopp_shift(lam, dopp_v)

        return lam_shifted, spec, spec_tot

    def _dopp_shift(self, lam, dopp_v):
        """Apply non-relativistic Doppler wavelength shift.

        Parameters
        ----------
        lam : astropy.units.Quantity
            Input wavelength array.
        dopp_v : astropy.units.Quantity
            Doppler shift velocity.
        Returns
        -------
        astropy.units.Quantity
            Shifted wavelength(s): lambda' = lambda * (1 + v/c).
        """

        return lam * (1 + dopp_v / c.c)

    def _calc_e(self, v, j):
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
        astropy.units.Quantity
            Energy in J (via h c k).

        Notes
        -----
        - All constants are collected and organized by Huber & Herzberg 1979, p. 250. Specific references given below.
        """
        om_e = (
            4400.39 * u.cm**-1
        )  # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 and residuals could be improved
        x_e = (
            120.815 * u.cm**-1
        ) / om_e  # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 and residuals could be improved
        y_e = (
            0.7242 * u.cm**-1
        ) / om_e  # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 and residuals could be improved
        B_e = (
            60.864 * u.cm**-1
        )  # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 but residuals are small beyond that too. Good enough for our purposes.
        D_e = 0.0471 * u.cm**-1  # Fink, Wiggins, Rank, JMS 18, 384 (1965).
        H_e = 4.9e-5 * u.cm**-1  # Fink, Wiggins, Rank, JMS 18, 384 (1965).
        a1, a2, a3 = (
            3.07638 * u.cm**-1,
            0.06017 * u.cm**-1,
            0.00481 * u.cm**-1,
        )  # Herzberg, Howe, CJP 16, 636 (1959); technically only valid for v=0...8 but residuals are small beyond that too. Good enough for our purposes.
        b1, b2 = 0.00274 * u.cm**-1, 0.00040 * u.cm**-1  # Fink, Wiggins, Rank, JMS 18, 384 (1965).
        c1 = 0.5e-5 * u.cm**-1  # Fink, Wiggins, Rank, JMS 18, 384 (1965).

        def single(vi, ji):
            Bv = B_e - a1 * (vi + 0.5) + a2 * (vi + 0.5) ** 2 - a3 * (vi + 0.5) ** 3
            Dv = -D_e + b1 * (vi + 0.5) - b2 * (vi + 0.5) ** 2
            Hv = H_e - c1 * (vi + 0.5)
            Gv = om_e * (vi + 0.5) - om_e * x_e * (vi + 0.5) ** 2 + om_e * y_e * (vi + 0.5) ** 3
            FJ = Bv * ji * (ji + 1) - Dv * ji**2 * (ji + 1) ** 2 + Hv * ji**3 * (ji + 1) ** 3
            return (c.h * c.c * (Gv + FJ)).to(u.eV)

        if np.ndim(v) == 1:
            return np.array([[single(vi, jj).value for jj in j] for vi in v]) * u.eV
        else:
            return single(v, j)

    def _calc_tau(self, nvj, siglu):
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
        astropy.units.Quantity
            Optical depth as a function of wavelength and transition.

        Notes
        -----
        Implements Eq. 11 (McJunkin et al. 2016).
        """
        return (nvj[:, None] * siglu).decompose()

    def _calc_tau_tot(self):
        """
        Sum optical depths over transitions.

        Returns
        -------
        array
            Optical depth as a function of wavelength.
        """
        self._tau_tot = self._tau.sum(axis=0)  # total tau(lambda)
        return self._tau_tot

    def _voigt(self, lam, lam0, gam, dv):
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
        array
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

    def _calc_siglu(self, lam, lamlu, Atot, dv, flu):
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
        astropy.units.Quantity
            sigma_lu(lambda) array in cm^2.

        Notes
        -----
        Implements Eq. 4 (McJunkin et al. 2016).
        """
        siglu = np.zeros((len(lamlu), len(lam))) * u.cm**2
        for i in range(len(lamlu)):
            H_prof = self._voigt(lam, lamlu[i], Atot[i], dv)
            siglu[i, :] = (np.sqrt(np.pi) * c.e.esu**2 / (c.m_e * c.c * dv) * flu[i] * lamlu[i] * H_prof).to(u.cm**2)
        return siglu

    def _calc_dv(self, instr=False):
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
        astropy.units.Quantity
            Combined Doppler width (same units as c.c).
        """
        dv_therm = np.sqrt(2 * c.k_B * self.constant.TH2 / (2 * c.m_p))  # Thermal broadening
        dv_nontherm = self.constant.VELOCITY_DISPERSION  # Non-thermal broadening
        if instr and self.constant.RESOLVING_POWER:  # Instrumental broadening
            dv_instr = c.c / (self.constant.RESOLVING_POWER * np.sqrt(8 * np.log(2)))
            return np.sqrt(dv_therm**2 + dv_nontherm**2 + dv_instr**2)
        return np.sqrt(dv_therm**2 + dv_nontherm**2)
