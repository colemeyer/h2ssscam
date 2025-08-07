import numpy as np
import importlib.resources
import configparser
import shutil
from datetime import datetime
import os

def load_data(filename):
    """Returns the NPZ file data in a dict.

    The data must be located in the 'src/h2ssscam/data' directory.

    Args:
        filename: the name of the file to fetch, but without the '.npz' suffix.

    Returns:
        A dict containing the data from the file.

    Raises:
        FileNotFoundError: The file could not be found in the 'src/h2ssscam/data' directory.
    """
    with importlib.resources.files("h2ssscam.data").joinpath(f"{filename}.npz").open("rb") as f:
        with np.load(f, allow_pickle=True) as data:
            return dict(data)


def load_config_files():
    config = configparser.ConfigParser()
    with importlib.resources.files("h2ssscam.data").joinpath(f"config.ini").open("r") as f:
        config.read_file(f)
    return config
