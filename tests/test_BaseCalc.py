"""
Contains the tests for the BaseCalc module.
"""
import numpy as np
import pytest
import astropy.units as u
from astropy.units import Quantity
from astropy.units.core import UnitConversionError
from h2ssscam.BaseCalc import BaseCalc

class TestBaseCalc:
    """
    Contains the tests for the BaseCalc class.
    """
    @pytest.fixture
    @staticmethod
    def base_calc() -> BaseCalc:
        """Return a new BaseCalc instance to be used for a test.

        Returns
        -------
        BaseCalc
            A new BaseCalc instance.
        """
        return BaseCalc()

    @staticmethod
    @pytest.mark.parametrize("lam,dopp_v,expected_dopp_shift", [
        [1500 * u.AA, 10 * u.km / u.s, 1500.0500346142799 * u.AA],
        [np.nan, 10 * u.km / u.s, Quantity(np.nan)],
        [1500 * u.AA, np.nan, None],
        [np.inf, 10 * u.km / u.s, Quantity(np.inf)],
        [1500 * u.AA, np.inf, None]
    ])
    def test_dopp_shift(base_calc, lam: Quantity, dopp_v: Quantity, expected_dopp_shift: Quantity):
        """
        Tests the BaseCalc._dopp_shift method.
        """
        lam_dopp_v_vals = [lam, dopp_v]
        if np.isnan(dopp_v) or np.isinf(dopp_v):
            with pytest.raises(UnitConversionError):
                _ = base_calc._dopp_shift(lam, dopp_v)
        else:
            actual_dopp_shift = base_calc._dopp_shift(lam, dopp_v)
            if np.isnan(lam):
                assert np.isnan(actual_dopp_shift.value) and np.isnan(expected_dopp_shift.value)
            elif np.isinf(lam):
                assert np.isinf(actual_dopp_shift.value) and np.isinf(expected_dopp_shift.value)
            else:
                assert actual_dopp_shift == expected_dopp_shift
