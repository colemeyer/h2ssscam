# ========== FILE: h2_model.py ==========
"""
Main script: load line data, compute populations, source function,
absorption rates, and plot emergent H₂ fluorescence spectrum.

To be incorporated:
- Dissociation spectrum
"""
import matplotlib.pyplot as plt
import astropy.units as u
from h2_model_funcs import *
from h2_model_params import *
import numpy as np
# from funkyfresh import set_style
# set_style('AAS', silent=True)

# ------------------------------------------------------ #
# ----- LOAD H2 LINE DATA ------------------------------ #
# ------------------------------------------------------ #

# Abgrall et al. (1993) fluorescence line list
s = np.load('h2fluor_data_Abgrall+1993.npz', allow_pickle=True)
Atot, Auldiss, Aul, lamlu, band, vu, ju, vl, jl = s['Atot'] * u.s**-1, s['Auldiss'] * u.s**-1, s['Aul'] * u.s**-1, s['lamlu'] * u.AA, s['band'], s['vu'], s['ju'], s['vl'], s['jl']

# Filter by v <= VMAX, J <= JMAX
mask_h2 = (vl <= VMAX) & (jl <= JMAX)
Atot, Auldiss, Aul, lamlu, band, vu, ju, vl, jl = Atot[mask_h2], Auldiss[mask_h2], Aul[mask_h2], lamlu[mask_h2], band[mask_h2], vu[mask_h2], ju[mask_h2], vl[mask_h2], jl[mask_h2]


# ------------------------------------------------------ #
# ----- LOAD HI LINE DATA ------------------------------ #
# ------------------------------------------------------ #

# NIST Atomic Spectral Database for HI
s = np.load('hi_data_NIST.npz')
hi_lamlu, hi_jl, hi_ju, hi_Aul, hi_flu = s['lamlu'] * u.AA, s['jl'], s['ju'], s['Aul'] * u.s**-1, s['flu']


# ------------------------------------------------------ #
# ----- PREPARATORY CALCULATIONS ----------------------- #
# ------------------------------------------------------ #

CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1
if UNIT == 'CU': UNIT = CU_UNIT
else: UNIT = ERG_UNIT

# Wavelength grid from 912 to 1800 Å (dlam = 0.01 angstroms)
lam0, lamend, dlam = 912, 1800, 0.1
lam = np.linspace(int(lam0), int(lamend), int((lamend - lam0) / dlam)) * u.AA

# Compute Doppler widths
dv_phys = calc_dv(T=TH2)               # thermal + non-thermal (Eq. 7)
dv_tot = calc_dv(T=TH2, R=RESOLVING_POWER, instr=True)      # include instrumental

# H2 oscillator strengths and level populations
flu = calc_flu(ju, jl, lamlu, Aul)  # Eq. 2
nvj = calc_nvj(NH2_TOT, TH2)        # Eq. 8
sel_levels = np.where(nvj[vl,jl] > NH2_CUTOFF)[0]
Atot_p, Auldiss_p, Aul_p, lamlu_p, band_p, vu_p, ju_p, vl_p, jl_p, flu_p = Atot[sel_levels], Auldiss[sel_levels], Aul[sel_levels], lamlu[sel_levels], band[sel_levels], vu[sel_levels], ju[sel_levels], vl[sel_levels], jl[sel_levels], flu[sel_levels]
nvj_p = nvj[vl_p, jl_p]

# HI calculations
NHI = boltzmann(NHI_TOT, hi_ju, hi_jl, hi_lamlu, THI)
hih2_lamlu, hih2_flu, hih2_Atot, hih2_N = np.append(hi_lamlu, lamlu_p), np.append(hi_flu, flu_p), np.append(hi_Aul, Atot_p), np.append(NHI, nvj_p)


# ------------------------------------------------------ #
# ----- SOURCE FUNCTION -------------------------------- #
# ------------------------------------------------------ #

# Compute absorption cross-sections and optical depths
siglu = calc_siglu(lam, hih2_lamlu, hih2_Atot, dv_phys, hih2_flu)  # Eq. 4
tau = calc_tau(hih2_N, siglu)                                 # Eq. 11
tau_tot = tau.sum(axis=0)                                    # total tau(lambda)

# Incident UV background and attenuated source
if INC_SOURCE == 'BLACKBODY': uv_inc = blackbody(lam, THI, unit=UNIT)
else: uv_inc = uv_continuum(lam, unit=UNIT)                                   # empirical cont.
source = uv_inc * np.exp(-tau_tot)

# Absorption rates for H2 only
tau_h2 = tau[len(hi_lamlu):, :]
abs_rate = calc_abs_rate(uv_inc, tau_h2, tau_tot, unit=UNIT) * dlam     # Eq. 12–13
abs_rate_per_trans = np.sum(abs_rate, axis=1)

# Plot source spectrum
plt.figure()
plt.plot(lam, source, lw=0.5)
plt.xlabel(r'Wavelength (\AA)')
if UNIT == CU_UNIT: plt.ylabel(r'Intensity (arbitrary units)') # CU)')
else: plt.ylabel(r'Intensity ($\mathrm{erg\;cm^{-2}\;s^{-1}\;arcsec^{-1}\;nm^{-1}}$)')
plt.title('Source Spectrum')
plt.show()


# ------------------------------------------------------ #
# ----- EMERGENT SPECTRUM ------------------------------ #
# ------------------------------------------------------ #

# Assemble arrays for relevant emission lines
vljl = []
h2_lamlu = []
h2_Atot = []
flux_per_trans = []
for ui in range(0, len(vu_p)):
    idx_u = np.where((vu == vu_p[ui]) & (ju == ju_p[ui]) & (band == band_p[ui]))[0]
    for idx in idx_u:
        if np.any((Aul[idx] / Atot[idx] < LINE_STRENGTH_CUTOFF) | (lamlu[idx] < BP_MIN) | (lamlu[idx] > BP_MAX)): continue

        vljl.append([vl[idx],jl[idx]])

        h2_lamlu.append(lamlu[idx].value)
        h2_Atot.append(Atot[idx].value)
        flux_per_trans.append((abs_rate_per_trans[ui] * (Aul[idx] / Atot[idx])).value)

h2_lamlu = np.array(h2_lamlu) * u.AA
h2_Atot = np.array(h2_Atot) * u.s**-1
flux_per_trans = np.array(flux_per_trans) * UNIT

lam0, lamend, dlam = 912, 1800, DLAM.to(u.AA).value
lam_highres = np.linspace(int(lam0), int(lamend), int((lamend - lam0) / dlam)) * u.AA
source = np.interp(lam_highres, lam, source)

### Calculate emergent spectrum
lam_shifted, spec, spec_tot = calc_spec(lam_highres, h2_lamlu, h2_Atot, dv_tot, flux_per_trans, source, UNIT, DOPPLER_SHIFT)

### Save emergent spectrum
np.savez_compressed(f'models/h2-fluor-model_R={RESOLVING_POWER}_TH2={int(TH2.value)}_NH2={int(np.log10(NH2_TOT.value))}_THI={int(THI.value)}_NHI={int(np.log10(NHI_TOT.value))}', lam_shifted=lam_shifted, spec=spec.to(UNIT).value, spec_tot=spec_tot.to(UNIT).value)

# Plot emission-only spectrum
plt.plot(lam_shifted, spec, c='k', lw=0.5)
plt.axvline(1608, 0, 1, c='r', lw=0.5, dashes=(8,4))
plt.xlabel(r'Wavelength (\AA)')
if UNIT == CU_UNIT: plt.ylabel(r'Intensity (arbitrary units)') # CU)')
else: plt.ylabel(r'Intensity ($\mathrm{erg\;cm^{-2}\;s^{-1}\;arcsec^{-1}\;nm^{-1}}$)')
plt.xlim([BP_MIN.value, BP_MAX.value])
plt.title('Emergent Spectrum')
plt.show()

# Plot total (emission + continuum) spectrum
plt.plot(lam_shifted, spec_tot, c='k', lw=0.5)
plt.axvline(1608, 0, 1, c='r', lw=0.5, dashes=(8,4))
plt.xlabel(r'Wavelength (\AA)')
if UNIT == CU_UNIT: plt.ylabel(r'Intensity (arbitrary units)') # CU)')
else: plt.ylabel(r'Intensity ($\mathrm{erg\;cm^{-2}\;s^{-1}\;arcsec^{-1}\;nm^{-1}}$)')
plt.title('Emergent Spectrum w/ Continuum')
plt.show()