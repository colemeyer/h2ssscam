"""
Contains the tests for the data_loader module.
"""
import os
import pickle
from pathlib import Path
import numpy.testing as npt
import pytest
from h2ssscam.data_loader import load_data

def test_load_data():
    """
    Tests for the correct behavior when the file does not exist.
    """
    with open(Path(os.path.abspath(__file__)).parent / "data" / "load_data_expected.pkl", "rb") as file:
        expected_result = pickle.load(file)

    npt.assert_equal(load_data("h2fluor_data_Abgrall+1993"), expected_result)

def test_load_data_nonexistent_file():
    """
    Tests for the correct behavior when the file does not exist.
    """
    with pytest.raises(FileNotFoundError):
        load_data("nonexistentFile")
