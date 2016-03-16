import os
from setuptools import setup, Extension

from Cython.Build import cythonize
import numpy


include_dirs = [numpy.get_include(), 'gromacs/src','gromacs/src/external/tng_io/include',]
library_dirs = []
if 'LD_LIBRARY_PATH' in os.environ:
    lib = os.environ['LD_LIBRARY_PATH'].split(':')[0]
    library_dirs.insert(0, lib)
    include = lib.replace('/lib', '/include')
    if os.path.exists(include):
        include_dirs.insert(0, include)
elif os.path.exists('gromacs/build/lib'):
    library_dirs.insert(0, 'gromacs/build/lib')

if not library_dirs:
    raise OSError("""
        Gromacs library not found.
        Activate a gromacs module or specify environment variable LD_LIBRARY_PATH.
        """)

library_dirs.append('gromacs/src/external/tng_io/build/lib')

library_dirs = [os.path.abspath(p) for p in library_dirs]


extensions = [
    # Extension('mdevaluate.gromacs.coordinates', [
    #           'mdevaluate/gromacs/coordinates.pyx'], include_dirs=include_dirs),
    # Extension('mdevaluate.gromacs.logarithmic', [
    #           'mdevaluate/gromacs/logarithmic.pyx'], include_dirs=include_dirs),
    Extension('pygmx.tpxio',
              sources=['pygmx/tpxio.pyx'],
              include_dirs=include_dirs,
              libraries=['gromacs'],
              library_dirs=library_dirs,
              runtime_library_dirs=library_dirs,
              language='c++'),
    Extension('pygmx.xtcio',
              sources=['pygmx/xtcio.pyx'],
              include_dirs=include_dirs,
              libraries=['gromacs'],
              library_dirs=library_dirs,
              runtime_library_dirs=library_dirs,
              language='c++'),
    Extension('pygmx.tngio',
              sources=['pygmx/tngio.pyx'],
              include_dirs=include_dirs,
              libraries=['gromacs'],
              library_dirs=library_dirs,
              runtime_library_dirs=library_dirs,
              #language='c++'
              ),

]

setup(
    name='pygmx',
    description='Python wrapper for gromacs library.',
    author_email='niels.mueller@physik.tu-darmstadt.de',
    packages=['pygmx', ],
    version='0.0.2',
    requires=['numpy', 'Cython'],
    ext_modules=cythonize(extensions),
)
