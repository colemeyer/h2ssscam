from datetime import datetime
import importlib.resources, os, sys


def create_config_file(output_path: str,file_name:str|None = None):
    """Creates config file based on the file with default parameters

    Parameters
    ----------
    output_path : str
        Directory in which file is created
    file_name : str
        Name of the config file, optional
    """
    if not isinstance(output_path, str):
        raise TypeError("Path has to be a string")
    if not os.path.exists(output_path):
        raise ValueError(f'Directory {output_path} does NOT exists')
    timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
    with importlib.resources.files("h2ssscam.data").joinpath(f"config.ini").open("r") as f:
        config_content = f.read()
    if file_name:
        output_file = os.path.join(output_path, f"{file_name}.ini")
        if os.path.isfile(output_file):
            raise ValueError('File already exists!')
    else:
        output_file = os.path.join(output_path, f"config_{timestamp}.ini")
    with open(output_file, "w") as out:
        out.write(config_content)
    print(f"Config file saved as {output_file}")


if __name__ == "__main__":

    output_path = sys.argv[1] if len(sys.argv) > 1 else "."
    output_filename = sys.argv[2] if len(sys.argv) > 2 else None
    create_config_file(output_path,output_filename)
