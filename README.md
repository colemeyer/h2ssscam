# h2ssscam
Synthetic H2 fluorescence code refactored by Cole Meyer, Daniel Lopez-Sanders, Cassandra S. Cruz, Dominik P. Pacholski during Code/Astro 2025

The code is based on  "EMPIRICALLY ESTIMATED FAR-UV EXTINCTION CURVES FOR CLASSICAL T TAURI STARS", McJunkin, M., France, K., Schindhelm, R., Herczeg, G., Schneider, P. C., & Brown, A. (2016), ApJ, 828, 69. 

![alt text](assets/img/SDG_logo.png "Title")


### Installation instructions
Here are some installation instructions for the average Anaconda user. (Note: in the instructions below, we will assume that you are using a virtual environment named `myenv`.) We've tested this using Python 3.10.
1. Activate your virtual environment:<br>
    `% conda activate myenv`
2. Install the `h2ssscam` package and its dependencies:<br>
    `% pip install h2ssscam`

### Usage instructions

You're now ready to use the `h2ssscam` package! To run the code, execute the following in the terminal:<br>
    `% python -m h2ssscam`

The user may also modify the input parameters for the model. To do so,
1. Open a Jupyter Notebook, Python script, or iPython and load the create_config_file function using<br>
    `from h2ssscam import create_config_file`
2. Create a configuration file by running<br>
    `create_config_file(directory)`,<br>
    where `directory` is the directory in which you would like to create the configuration file.
3. Modify variables in the generated file and save. Run `h2ssscam` as usual using<br>
    `python -m h2ssscam [config_file_path]`.

**Note that package dependencies for `h2ssscam` are managed using the `uv` package manager. Advanced users are encouraged to run `h2ssscam` using `uv` functionalities (e.g., `uv run h2ssscam`).**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) ![A rectangular badge, half black half purple containing the text made at Code Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet.svg)