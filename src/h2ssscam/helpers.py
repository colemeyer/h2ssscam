from constants import TH2,RESOLVING_POWER, VELOCITY_DISPERSION, RESOLVING_POWER
from dataclasses import dataclass
from astropy.units import Quantity
from scipy.special import wofz
import astropy.constants as c
import numpy as np

@dataclass
class BaseCalc:
    dv_phys : Quantity = None
    dv_tot : Quantity = None
    tau : Quantity = None

    def calculate_doppler_widths(self, TH2, RESOLVING_POWER):
        # Compute Doppler widths
        self.dv_phys = self._calc_dv(T=TH2)               # thermal + non-thermal (Eq. 7)
        self.dv_tot = self._calc_dv(T=TH2, R=RESOLVING_POWER, instr=True)      # include instrumental
        return self.dv_phys, self.dv_tot

    def optical_depth(self, lam, hih2_lamlu, hih2_Atot, dv_phys, hih2_flu, hih2_N):
        siglu = self._calc_siglu(lam, hih2_lamlu, hih2_Atot, dv_phys, hih2_flu)  # Eq. 4
        self.tau = self._calc_tau(hih2_N, siglu)                                 # Eq. 11)
        return self.tau

    def total_optical_depth(self):
        self.tau_tot = self.tau.sum(axis=0) # total tau(lambda)
        return self.tau_tot

    def _calc_tau(nvj, siglu):
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

    def _voigt(lam, lam0, gam, dv):
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
        siglu : astropy.units.Quantity
            sigma_lu(lambda) array in cm^2.

        Notes
        -----
        Implements Eq. 4 (McJunkin et al. 2016).
        """
        siglu = np.zeros((len(lamlu), len(lam))) * u.cm**2
        for i in range(len(lamlu)):
            H_prof = self._voigt(lam, lamlu[i], Atot[i], dv)
            siglu[i,:] = (np.sqrt(np.pi) * c.e.esu**2 / (c.m_e * c.c * dv)
                        * flu[i] * lamlu[i] * H_prof).to(u.cm**2)
        return siglu

    def _calc_dv(self,T=TH2, b=VELOCITY_DISPERSION, R=RESOLVING_POWER, instr=False):
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
