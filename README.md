<img align="center" src="https://raw.githubusercontent.com/carlosborca/Smiles4Psi/main/media/logo/Logo.png" height=260>

Automated conformational search with a genetic algorithm interfaced with a Psi4-input generator for conformers that requires only SMILES as input

| Category | Badges |
|-------------|-------------|
| **Foundation** | [![License](https://img.shields.io/github/license/carlosborca/Smiles4Psi.svg)](https://opensource.org/licenses/LGPL-3.0) [![GitHub Top Languages](https://img.shields.io/github/languages/top/carlosborca/Smiles4Psi)](https://github.com/carlosborca/Smiles4Psi/) |
| **GitHub Info** | [![GitHub Code Size](https://img.shields.io/github/languages/code-size/carlosborca/Smiles4Psi)](https://github.com/carlosborca/Smiles4Psi/) [![GitHub Commits per Month](https://img.shields.io/github/commit-activity/m/carlosborca/Smiles4Psi)](https://github.com/carlosborca/Smiles4Psi/) [![GitHub Last Commit](https://img.shields.io/github/last-commit/carlosborca/Smiles4Psi)](https://github.com/carlosborca/Smiles4Psi/) |

## Overview

Smiles4Psi is a Python 3 program that automates the creation of Psi4 input files for 3D-structures obtained from a genetic algorithm conformational search. The conformers are generated via Open Babel. The only input Smiles4Psi requires is a SMILES string description of the molecule.

## General Information

To run Smiles4Psi, the code requires Python 3.7 and its modules `psi4`, `psi4-rt`, `openbabel`, and `psutil`. So, the instructions to download and install Smiles4Psi and to create a _conda environment_ that includes these modules are presented below.

### Installation

Minimal set of commands to install Smiles4Psi on Linux, MacOS, or Windows (with the Windows Subsystem for Linux). Last tested on 16 August 2021:

#### 1. Install Miniconda:

If you have an installation of _Conda_ in your system, please skip to step 2. Otherwise, _Miniconda_ is required and the installer is available from the Anaconda website. To download the installer from the terminal (in Linux or the Windows Subsystem for Linux):

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

(_Note_) If you use MacOS, replace `Linux` by `MacOSX` in the previous command and the subsequent ones.

Run the installer following on-screen instructions:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

After the installation is complete, close the terminal and start a new shell.

(_Optional_) Disable automatic activation of the _base_ conda environment:

```
conda config --set auto_activate_base false
```

#### 2. Create a _Conda Environment_ for Smiles4Psi

Smiles4Psi requires Python 3 and several modules. Conda offers the possibility of creating an _environment_ that contains all the dependencies required by Smiles4Psi. To download and install the required software tools in a new _S4P_ environment execute the command below and follow on-screen instructions:

```
conda create -n S4P psi4=1.3.2 psi4-rt=1.3 python=3.7 openbabel=3.1.0 psutil=5.8.0 -c psi4 -c conda-forge
```

#### 3. Activate the _S4P_ environment

To use the recently created _S4P_ environment, activate it:

```
conda activate S4P
```

#### 4. Clone Smiles4Psi from its GitHub repository:

In your file system, navigate to the location where you would like to place the root directory of Smiles4Psi and clone it from its corresponding GitHub repository:

```
git clone https://github.com/carlosborca/Smiles4Psi.git
```

### How to run Smiles4Psi

Documentation about the code and runtime options, including examples and defaults, can be found by executing with the `-h` or `--help` option.

```
python Smiles4Psi.py --help
```

The code uses the _current working directory_ as a reference point for execution and will produce new files and directories on that current working directory upon execution. Therefore, one needs to navigate to the folder where the files are meant to be produced and execute the Python 3 code there, within the _S4P_ environment.

```
conda activate S4P
cd /To/Directory/Where/The/Files/Are/Meant/To/Be/Created/
python /Path/to/Smiles4Psi/Smiles4Psi.py
```

The first and third command can be grouped and automated by doing something _horrific_ like creating an alias in `~/.bashrc`, for example:

```
alias Smiles4Psi="conda activate S4P; python ~/Gits/Smiles4Psi/Smiles4Psi.py; conda deactivate"
```

#### Copyright

Copyright (c) 2021, Carlos H. Borca
