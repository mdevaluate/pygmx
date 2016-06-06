# Pygmx

A Cython wrapper for the shared library of GROMACS 5.1.

## Installation

Pygmx requires some python packages and an installation of GROMACS 5.1.

### Python requirements

* Cython
* NumPy

### Gromacs installation

#### Using a module

If gromacs is installed as a module in your system, run

    module load gromacs/5.1

and skip to the section *Installing pygmx*.

#### Manual setup

The shared library of GROMACS 5.1and the corresponding c++ header files need to be present on the system.
If the shared library is installed globally, the environment variable `LD_LIBRARY_PATH` must be set.
The installation process will look for the header files by replacing any `lib` folder
in `LD_LIBRARY_PATH` with `include`.

If no header files are present, simply pull the submodule `gromacs` in this repository.
To build the shared library, follow the [official installation instructions](http://manual.gromacs.org/documentation/5.1.2/install-guide/index.html),
starting with step 4 (creating a build directory).

### Installing pygmx

When the requirements are met, the installation should be easy.
Navigate to the top folder of the repository and run the command

    python setup.py install

This builds the cython modules and installs them into the local python distribution.


### Deploy on intranet

For any installed version of mdevaluate modules do:

    version=dev # or maybe: version=$(python ../mdevaluate/setup.py --version)
    module load gromacs/5.1
    module load mdevaluate/$version
    python setup.py install --prefix=/autohome/niels/modules/mdevaluate-$version
