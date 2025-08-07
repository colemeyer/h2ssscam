import numpy as np
import importlib.resources


def load_data(filename):
    with importlib.resources.files("h2ssscam.data").joinpath(f"{filename}.npz").open("rb") as f:
        with np.load(f, allow_pickle=True) as data:
            return dict(data)


def load_config_file():
    pass

def save_config_file():
    pass