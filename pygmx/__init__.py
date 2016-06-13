"""
Python wrapper for the gromacs library.
"""

import os

from .tpxio import TPXReader
from .xtcio import XTCReader
# from .enxio import EDRFile
from .errors import FileTypeError
from .gromacs.reader import index_filename_for_xtc
from .gromacs.xtcindex import index_xtcfile

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
            if not os.path.exists(indexfile):
                print('Gneerating Index for xtc file. This may take a while...')
                index_xtcfile(filename)
            return FILE_EXTENSIONS[ext](filename, indexfile)
        else:
            return FILE_EXTENSIONS[ext](filename)
    else:
        raise FileTypeError('Filetype {} not supported by pygmx.'.format(ext))
