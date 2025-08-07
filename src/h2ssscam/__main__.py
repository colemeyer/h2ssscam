# ========== FILE: h2_model.py ==========
"""
Main script: load line data, compute populations, source function,
absorption rates, and plot emergent H₂ fluorescence spectrum.

To be incorporated:
- Dissociation spectrum
"""
import os, sys
import pathlib
import astropy.units as u
from h2ssscam.BaseCalc import BaseCalc
from h2ssscam.plotting_funcs import *
import numpy as np
from h2ssscam.data_loader import load_data
from h2ssscam.Constants import Constants

# from funkyfresh import set_style
# set_style('AAS', silent=True)

base_dir = pathlib.Path(os.path.abspath(__file__)).parent
data_dir = base_dir / "data"
models_dir = base_dir / "models"


def main():

    config_file_path = sys.argv[1] if len(sys.argv) == 2 else None
    constant = Constants(config_file_path)
    basecalc = BaseCalc(constant)

    # ------------------------------------------------------ #
    # ----- LOAD H2 LINE DATA ------------------------------ #
    # ------------------------------------------------------ #

    # Abgrall et al. (1993) fluorescence line list
    s = load_data("h2fluor_data_Abgrall+1993")
    Atot, Auldiss, Aul, lamlu, band, vu, ju, vl, jl = (
        s["Atot"] * u.s**-1,
        s["Auldiss"] * u.s**-1,
        s["Aul"] * u.s**-1,
        s["lamlu"] * u.AA,
        s["band"],
        s["vu"],
        s["ju"],
        s["vl"],
        s["jl"],
    )

    # Filter by v <= VMAX, J <= JMAX
    mask_h2 = (vl <= constant.VMAX) & (jl <= constant.JMAX)
    Atot, Auldiss, Aul, lamlu, band, vu, ju, vl, jl = (
        Atot[mask_h2],
        Auldiss[mask_h2],
        Aul[mask_h2],
        lamlu[mask_h2],
        band[mask_h2],
        vu[mask_h2],
        ju[mask_h2],
        vl[mask_h2],
        jl[mask_h2],
    )

    # ------------------------------------------------------ #
    # ----- LOAD HI LINE DATA ------------------------------ #
    # ------------------------------------------------------ #

    # NIST Atomic Spectral Database for HI
    s = load_data("hi_data_NIST")
    hi_lamlu, hi_jl, hi_ju, hi_Aul, hi_flu = s["lamlu"] * u.AA, s["jl"], s["ju"], s["Aul"] * u.s**-1, s["flu"]

    # ------------------------------------------------------ #
    # ----- PREPARATORY CALCULATIONS ----------------------- #
    # ------------------------------------------------------ #

    CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
    ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1
    if constant.UNIT == "CU":
        units = CU_UNIT
    else:
        units = ERG_UNIT

    # Wavelength grid from 912 to 1800 Å (dlam = 0.01 angstroms)
    lam0, lamend, dlam = 912, 1800, 0.1
    lam = np.linspace(int(lam0), int(lamend), int((lamend - lam0) / dlam)) * u.AA

    # Compute Doppler widths
    basecalc.dv_phys
    basecalc.dv_tot

    # H2 oscillator strengths and level populations
    flu = basecalc.calc_flu(ju, jl, lamlu, Aul)  # Eq. 2
    nvj = basecalc.calc_nvj(constant.NH2_TOT, constant.TH2)  # Eq. 8
    sel_levels = np.where(nvj[vl, jl] > constant.NH2_CUTOFF)[0]
    Atot_p, lamlu_p, band_p, vu_p, ju_p, vl_p, jl_p, flu_p = (
        Atot[sel_levels],
        lamlu[sel_levels],
        band[sel_levels],
        vu[sel_levels],
        ju[sel_levels],
        vl[sel_levels],
        jl[sel_levels],
        flu[sel_levels],
    )
    nvj_p = nvj[vl_p, jl_p]

    # HI calculations
    NHI = basecalc.boltzmann(constant.NHI_TOT, hi_ju, hi_jl, hi_lamlu, constant.THI)
    hih2_lamlu, hih2_flu, hih2_Atot, hih2_N = (
        np.append(hi_lamlu, lamlu_p),
        np.append(hi_flu, flu_p),
        np.append(hi_Aul, Atot_p),
        np.append(NHI, nvj_p),
    )

    # ------------------------------------------------------ #
    # ----- SOURCE FUNCTION -------------------------------- #
    # ------------------------------------------------------ #

    # Compute absorption cross-sections and optical depths
    # SPEED TESTING: following 3 lines took ~0.6 seconds to run on 2023 Mac Pro M3 Pro chip
    basecalc._siglu = basecalc._calc_siglu(lam, hih2_lamlu, hih2_Atot, basecalc.dv_phys, hih2_flu)
    tau = basecalc.tau(hih2_N)
    tau_tot = basecalc.tau_tot

    # Incident UV background and attenuated source
    if constant.INC_SOURCE == "BLACKBODY":
        uv_inc = basecalc.blackbody(lam, constant.THI, unit=units)
    else:
        uv_inc = basecalc.uv_continuum(lam, unit=units)  # empirical cont.
    source = uv_inc * np.exp(-tau_tot)

    # Absorption rates for H2 only
    tau_h2 = tau[len(hi_lamlu) :, :]
    abs_rate = basecalc.calc_abs_rate(uv_inc, tau_h2, tau_tot, unit=units) * dlam  # Eq. 12–13
    abs_rate_per_trans = np.sum(abs_rate, axis=1)

    # Plot source spectrum
    plot_spectrum(lam, source, units=units, title=r"Source Spectrum", show=True)

    # ------------------------------------------------------ #
    # ----- EMERGENT SPECTRUM ------------------------------ #
    # ------------------------------------------------------ #

    # Assemble arrays for relevant emission lines
    # SPEED TESTING: following 14 lines took ~1.7 seconds to run on 2023 Mac Pro M3 Pro chip
    vljl = []
    h2_lamlu = []
    h2_Atot = []
    flux_per_trans = []
    for ui in range(0, len(vu_p)):
        idx_u = np.where((vu == vu_p[ui]) & (ju == ju_p[ui]) & (band == band_p[ui]))[0]
        for idx in idx_u:
            if np.any(
                (Aul[idx] / Atot[idx] < constant.LINE_STRENGTH_CUTOFF)
                | (lamlu[idx] < constant.BP_MIN)
                | (lamlu[idx] > constant.BP_MAX)
            ):
                continue

            vljl.append([vl[idx], jl[idx]])

            h2_lamlu.append(lamlu[idx].value)
            h2_Atot.append(Atot[idx].value)
            flux_per_trans.append((abs_rate_per_trans[ui] * (Aul[idx] / Atot[idx])).value)

    h2_lamlu = np.array(h2_lamlu) * u.AA
    h2_Atot = np.array(h2_Atot) * u.s**-1
    flux_per_trans = np.array(flux_per_trans) * units

    lam0, lamend, dlam = 912, 1800, constant.DLAM.to(u.AA).value
    lam_highres = np.linspace(int(lam0), int(lamend), int((lamend - lam0) / dlam)) * u.AA
    source = np.interp(lam_highres, lam, source)

    ### Calculate emergent spectrum
    # SPEED TESTING: following line took ~5.3 seconds to run on 2023 Mac Pro M3 Pro chip
    lam_shifted, spec, spec_tot = basecalc.calc_spec(
        lam_highres, h2_lamlu, h2_Atot, basecalc._dv_tot, flux_per_trans, source, units, constant.DOPPLER_SHIFT
    )

    ### Save emergent spectrum
    np.savez_compressed(
        f"h2-fluor-model_R={constant.RESOLVING_POWER}_TH2={int(constant.TH2.value)}_NH2={int(np.log10(constant.NH2_TOT.value))}_THI={int(constant.THI.value)}_NHI={int(np.log10(constant.NHI_TOT.value))}",
        lam_shifted=lam_shifted,
        spec=spec.to(units).value,
        spec_tot=spec_tot.to(units).value,
    )

    # Plot emission-only spectrum
    plot_spectrum(
        lam_shifted,
        spec,
        xmin=constant.BP_MIN.value,
        xmax=constant.BP_MAX.value,
        ylabel=r"Intensity (arbitrary units)",
        title=r"Emergent Spectrum",
    )
    plt.axvline(1608, 0, 1, c="r", lw=0.5, dashes=(8, 4))
    plt.show()

    # Plot total (emission + continuum) spectrum
    plot_spectrum(
        lam_shifted, spec_tot, ylabel=r"Intensity (arbitrary units)", title=r"Emergent Spectrum w/ Continuum"
    )
    plt.axvline(1608, 0, 1, c="r", lw=0.5, dashes=(8, 4))
    plt.show()


if __name__ == "__main__":
    main()
