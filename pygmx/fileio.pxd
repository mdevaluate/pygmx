
from utility cimport *

cdef extern from "gromacs/fileio/gmxfio.h":
    ctypedef struct t_fileio:
        pass

    int gmx_fio_seek(t_fileio *fio, gmx_off_t fpos)

    gmx_off_t gmx_fio_ftell(t_fileio *fio)

    int xtc_seek_time(t_fileio *fio, real time, int natoms, gmx_bool bSeekForwardOnly)
