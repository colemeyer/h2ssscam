from datetime import datetime
import importlib.resources, os, sys


def create_config_file(output_path: str):
    """Creates config file based on the file with default parameters

    Parameters
    ----------
    output_path : str
        Directory in which file is created
    """
    if not isinstance(output_path, str):
        raise TypeError("Path has to be a string")
    if not os.path.exists(output_path):
        raise ValueError(f'Directory {output_path} does NOT exists')
    timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
    with importlib.resources.files("h2ssscam.data").joinpath(f"config.ini").open("r") as f:
        config_content = f.read()
    output_file = os.path.join(output_path, f"config_{timestamp}.ini")
    with open(output_file, "w") as out:
        out.write(config_content)
    print(f"Config file saved as {output_file}")


if __name__ == "__main__":

    output = sys.argv[1] if len(sys.argv) > 1 else "."
    create_config_file(output)
