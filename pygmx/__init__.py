"""
Python wrapper for the gromacs library for their file formats.

Currently this supports xtc and tpr files within the `open` function.
Trajectories in trr format may be read with `.gromacs.reader.TRRReader`, which is experimental.
"""

import os

from .tpxio import TPXReader, make_xtcframe_whole
from .xtcio import XTCReader, read_xtcframe, append_xtcfile
from .trrio import TRRReader
from .enxio import EDRFile
from .errors import FileTypeError
from .gromacs.reader import index_filename_for_xtc
from .gromacs.xtcindex import index_xtcfile

FILE_EXTENSIONS = {
    'xtc': XTCReader,
    'trr': TRRReader,
    'tpr': TPXReader,
    'edr': EDRFile,
}


def open(filename, ignore_index_timestamps=False):
    """Open a supported gromacs file. Currently supported file formats: tpr, xtc."""
    ext = filename.split('.')[-1]
    if ext in FILE_EXTENSIONS:
        if ext in ['xtc', 'trr']:
            indexfile = index_filename_for_xtc(filename)
            if not os.path.exists(indexfile) and ext == 'xtc':
                print('Generating Index for xtc file. This may take a while...')
                index_xtcfile(filename)
            return FILE_EXTENSIONS[ext](filename, indexfile, ignore_timestamps=ignore_index_timestamps)
        else:
            return FILE_EXTENSIONS[ext](filename)
    else:
        raise FileTypeError('Filetype {} not supported by pygmx.'.format(ext))
