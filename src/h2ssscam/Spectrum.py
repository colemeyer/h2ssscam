# from dataclasses import dataclass
# from astropy.units import Quantity

# @dataclass
# class Spectrum:
#     wavelength: Quantity
#     yData: Quantity


# @dataclass
# class OpticalDepth(Spectrum):
#     # optical_depth: 
#     _optical_depth: None | Quantity = None
        
#     @property
#     def optical_depth(self):
#         if self._optical_depth is None:
#             self._optical_depth = self.calc_optical_depth()
#         return self._optical_depth
#     def calc_optical_depth(nvj, siglu):
#         """
#         Compute optical depth tau(v,J,lambda).

#         Parameters
#         ----------
#         nvj : array
#             Level populations (v,J).
#         siglu : astropy.units.Quantity
#             Absorption cross-sections.

#         Returns
#         -------
#         tau : astropy.units.Quantity
#             Optical depth array.

#         Notes
#         -----
#         Implements Eq. 11 (McJunkin et al. 2016).
#         """
#         return (nvj[:,None] * siglu).decompose()
