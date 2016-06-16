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

#### Through package manager

Pygmx requires the shared library and header files of Gromacs 5.1 or higher to be installed.
On many Unix distributions Gromacs may be installed through the package manager.
The required packages, which provide the headers, are usually named `gromacs-devel` or `gromcas-dev`.

#### Manual build

To build the shared library manually, follow the [official installation instructions](http://manual.gromacs.org/documentation/5.1.2/install-guide/index.html).
Make sure to source the `GMXRC` script, which sets the environment variable `LD_LIBRARY_PATH` to the correct location or set this variable manually.
Don't bother with advanced features like GPU support if the Gromacs installation will only be used for pygmx,
since the file io functions do not use these features.

### Installing pygmx

When the requirements are met, the installation should be easy.
Navigate to the top folder of the repository and run the command

    python setup.py install

This builds the cython modules and installs them into the local python distribution.
