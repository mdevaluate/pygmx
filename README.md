# Pygmx

A Cython wrapper for the shared library of `GROMACS 5.1`.

## Installation

Pygmx requires some python packages and an installation of `GROMACS 5.1`.

### Python requirements

* Cython
* NumPy

### Gromacs installation

The shared library of `GROMACS 5.1`and the corresponding c++ header files need to be present on the system.
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
