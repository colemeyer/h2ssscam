import numpy as np
import importlib.resources

def load_data(filename):
    with importlib.resources.files("h2ssscam.data").joinpath(f"{filename}.npz").open("rb") as f:
            data = np.load(f)
    return data
