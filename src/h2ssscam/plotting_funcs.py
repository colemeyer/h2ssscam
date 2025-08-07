import matplotlib.pyplot as plt
import astropy.units as u

CU_UNIT = u.ph * u.cm**-2 * u.s**-1 * u.sr**-1 * u.AA**-1
ERG_UNIT = u.erg * u.cm**-2 * u.s**-1 * u.arcsec**-2 * u.nm**-1

def plot_spectrum(wavelengths,
                  intensity,
                  xmin=None, xmax=None,
                  ymin=None, ymax=None,
                  units = None,
                  ylabel=None,
                  title=r"Spectrum",
                  show=False):

    plt.figure()
    plt.plot(wavelengths, intensity, lw=0.5)
    plt.title(title)
    plt.xlabel(r"Wavelength (\AA)")

    if (ylabel is not None): print('case 1'); plt.ylabel(ylabel)
    elif (units is not None):
        print('case 2')
        if units == CU_UNIT: unit_string = r"$\mathrm{ph\;cm^{-2}\;s^{-1}\;sr^{-1}\;}\AA^{-1}$"
        elif units == ERG_UNIT: unit_string = r"$\mathrm{erg\;cm^{-2}\;s^{-1}\;arcsec^{-1}\;nm^{-1}}$"
        else: raise ValueError('unknown units')
        plt.ylabel(r"Intensity ("+unit_string+")")
    else: print('case 3'); plt.ylabel(r"Intensity (arbitrary units)")

    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    if show: plt.show()
    

