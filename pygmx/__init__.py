"""
Python wrapper for the gromacs library.
"""

from .tpxio import TPXReader
from .xtcio import XTCReader
# from .enxio import EDRFile
from .errors import FileTypeError
from .gromacs.reader import index_filename_for_xtc

FILE_EXTENSIONS = {
    'xtc': XTCReader,
    'tpr': TPXReader
}


def open(filename):
    """Open a supported gromacs file with the appropiate reader."""
    ext = filename.split('.')[-1]
    if ext in FILE_EXTENSIONS:
        if ext in ['xtc']:
            indexfile = index_filename_for_xtc(filename)
            return FILE_EXTENSIONS[ext](filename, indexfile)
    else:
        raise FileTypeError('Filetype {} not supported by pygmx.'.format(ext))
