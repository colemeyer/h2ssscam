import numpy as np
import importlib.resources
import configparser
import shutil
from datetime import datetime
import os


def load_data(filename):
    with importlib.resources.files("h2ssscam.data").joinpath(f"{filename}.npz").open("rb") as f:
        with np.load(f, allow_pickle=True) as data:
            return dict(data)


def load_config_files():
    config = configparser.ConfigParser()
    with importlib.resources.files("h2ssscam.data").joinpath(f"config.ini").open("r") as f:
        config.read_file(f)
    return config
