# Pygmx

A Cython wrapper for the shared library of GROMACS 5.1.

## Installation

Pygmx requires some python packages and an installation of GROMACS 5.1.

### Python requirements

Pygmx needs mainly two python packages to be installed, which are available in any major Python distribution.

* Cython
* NumPy

### Gromacs installation

Pygmx requires the shared library and header files of Gromacs 5.1 or higher to be installed.

Note that pygmx also supports the Gromacs-2016 version.
Checkout the followng git branch if this version of Gromacs is desired:

    git checkout gromacs-2016

Note that Gromacs is backwards compatible, i. e. files produced by the 5.1 version will also work in the 2016 version.
But Gromacs-2016 is required to load files of the 2016 version.

#### Using a module

If gromacs is installed as a module in your system, run

    module load gromacs/5.1

and skip to the section *Installing pygmx*.
Note that on the AG Vogel intranet pygmx is available through the module `mdevaluate`.

#### Through package manager

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
To check if everyting works run

    python -c 'import pygmx'

If this raises an error that refers to `libgromacs.so`, this means the Gromacs library could not be found.
Make sure it is either in a standard location or `LD_LIBRARY_PATH` is set correctly.
