# Pygmx

A Cython wrapper for the shared library of GROMACS 2016

## Installation

Pygmx requires some python packages and an installation of GROMACS 2016.

### Python requirements

Pygmx needs mainly two python packages to be installed, which are available in any major Python distribution.

* Cython
* NumPy

### Gromacs installation

Pygmx requires the shared library and header files of Gromacs 2016 to be installed.

#### Through package manager

On many Unix distributions, Gromacs may be installed through the package manager.
The required packages, which provide the headers, are usually named `gromacs-devel` or `gromcas-dev`.

#### Manual build

To build the shared library manually, follow the [official installation instructions](http://manual.gromacs.org/documentation/2016.4/install-guide/index.html).
Make sure to source the `GMXRC` script, which sets the environment variable `LD_LIBRARY_PATH` to the correct location or set this variable manually.
Don't bother with advanced features like GPU support if the Gromacs installation will only be used for pygmx,
since the file io functions do not use these features.

### Installing pygmx

When the requirements are met, installing pygmx into the local Python installation can be achieved by navigating to the top folder of the repository and running the command

    python setup.py install

This builds the cython modules and installs them into the local python distribution.
To check if everyting works, run the following command from a directory out of the repository

    python -c 'import pygmx'

If this raises an error that refers to `libgromacs.so`, the Gromacs library could not be found.
Make sure it is either in a standard location or `LD_LIBRARY_PATH` is set correctly.
