# h2ssscam
Synthetic H2 fluorescence code refactored by Cole Meyer, Daniel Lopez-Sanders, Cassandra S. Cruz, Dominik P. Pacholski during Code/Astro 2025.

The code is based on  ''Empirically estimated far-UV extinction curves for classical T Tauri stars'', McJunkin, M., France, K., Schindhelm, R., Herczeg, G., Schneider, P. C., & Brown, A. (2016), ApJ, 828, 69.

<img src="https://raw.githubusercontent.com/colemeyer/h2ssscam/main/assets/img/SDG_logo.png" alt="Purple square logo with white border, showing a snake silhouette made of binary code and surrounded by asterisks, with the text H2SSSSCAM in stylized font at the bottom." title="h2ssscam logo" width="400">

### Installation instructions
Here are some installation instructions for the average Anaconda user. (Note: in the instructions below, we will assume that you are using a virtual environment named `myenv`.) We've tested this using Python 3.10.
1. Activate your virtual environment:<br>
    `% conda activate myenv`
2. Install the `h2ssscam` package and its dependencies:<br>
    `% pip install h2ssscam`

### Usage instructions

You're now ready to use the `h2ssscam` package! To run the code, execute the following in the terminal:<br>
    `% python -m h2ssscam`

The user may also modify the input parameters for the model. To do so, navigate into the desired directory and run in the terminal<br>
    `% python -m h2ssscam.create_config_file [directory] [config file name]`<br>
If the user would like to save the config file in the current directory, they should specify `.` in place of `[directory]`. Note that the file extension for the created configuration file will always be `.ini` (even if the user specifies some other extension) to satisfy code requirements. Modify the desired parameters and run the model as usual by executing the following in the terminal:<br>
    `% python -m h2ssscam [directory]/[config file name]`

Complete documentation can be found on Read the Docs: [https://h2ssscam.readthedocs.io/en/latest/index.html](https://h2ssscam.readthedocs.io/en/latest/index.html)<br>
The paper on which this code is based, McJunkin et al. 2016, can be found at [https://iopscience.iop.org/article/10.3847/0004-637X/828/2/69/meta](https://iopscience.iop.org/article/10.3847/0004-637X/828/2/69/meta)<br>
More information about the Code/Astro workshop during which this package was developed can be found at [https://semaphorep.github.io/codeastro/](https://semaphorep.github.io/codeastro/)

**Note that package dependencies for `h2ssscam` are managed using the `uv` package manager. Advanced users are encouraged to run `h2ssscam` using `uv` functionalities (e.g., `uv run h2ssscam`).**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) ![A rectangular badge, half black half purple containing the text made at Code Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet.svg)
