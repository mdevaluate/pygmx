import os
from setuptools import setup, Extension

from Cython.Build import cythonize
import numpy


def check_header_version(include_path):
    with open(os.path.join(include_path, 'gromacs/version.h')) as f:
        for l in f.readlines():
            if '#define GMX_API_VERSION' in l:
                version = int(l.split()[-1])
                assert version >= 50100, 'Installed gromacs version is too low!'
                return
    print('Gromacs version could not be checked.')

include_dirs = []
library_dirs = []

if 'gromacs' in os.environ.get('LD_LIBRARY_PATH', ''):
    for p in os.environ['LD_LIBRARY_PATH'].split(':'):
        if 'gromacs' in p:
            library_dirs.append(p)
            lib = p
            gmx_root = lib.split('lib')[0]
            include = os.path.join(gmx_root, 'include')
            if os.path.exists(include):
                include_dirs.append(include)
                check_header_version(include)

extensions = [
    Extension(
        'pygmx.gromacs.coordinates',
        ['pygmx/gromacs/coordinates.pyx'],
        include_dirs=include_dirs
    ),
    Extension(
        'pygmx.gromacs.logarithmic',
        ['pygmx/gromacs/logarithmic.pyx'],
        include_dirs=include_dirs
    ),
    Extension(
        'pygmx.tpxio',
        sources=['pygmx/tpxio.pyx'],
        include_dirs=include_dirs,
        libraries=['gromacs'],
        library_dirs=library_dirs,
        language='c++'
    ),
    Extension(
        'pygmx.xtcio',
        sources=['pygmx/xtcio.pyx'],
        include_dirs=include_dirs,
        libraries=['gromacs'],
        library_dirs=library_dirs,
        language='c++'
    ),
    Extension('pygmx.enxio',
              sources=['pygmx/enxio.pyx'],
              include_dirs=include_dirs,
              libraries=['gromacs'],
              library_dirs=library_dirs,
              language='c++'),

#    Extension('pygmx.tngio',
#              sources=['pygmx/tngio.pyx'],
#              include_dirs=include_dirs,
#              libraries=['gromacs'],
#              library_dirs=library_dirs,
#              runtime_library_dirs=library_dirs,
              #language='c++'
#              ),
]

setup(
    name='pygmx',
    description='Python wrapper around the gromacs library for file io.',
    author_email='niels.mueller@physik.tu-darmstadt.de',
    packages=['pygmx', 'pygmx.gromacs'],
    version='0.1',
    requires=['numpy', 'Cython'],
    ext_modules=cythonize(extensions),
)
