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


def create_config_file(output_path: str):
    timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
    with importlib.resources.files("h2ssscam.data").joinpath(f"config.ini").open("r") as f:
        config_content = f.read()

    output_file = os.path.join(output_path, f"config_{timestamp}.ini")
    with open(output_file, "w") as out:
        out.write(config_content)
    print(f"Config file saved as {output_file}")
