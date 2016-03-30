# from libc.stdio cimport FILE

from utility cimport *

#cdef extern from "gromacs/fileio/gmx_system_xdr.h":
#    ctypedef struct XDR:
#        pass


cdef extern from "gromacs/fileio/gmxfio.h":

    ctypedef struct t_fileio:
        pass

    void gmx_fio_rewind(t_fileio *fio)

    int gmx_fio_seek(t_fileio *fio, gmx_off_t fpos)

    gmx_off_t gmx_fio_ftell(t_fileio *fio)

#    int xtc_seek_frame(t_fileio *fio, int frame, int natoms)
    int xtc_seek_time(t_fileio *fio, real time, int natoms, gmx_bool bSeekForwardOnly)

    void gmx_fio_setdebug(t_fileio *fio, gmx_bool bDebug)
