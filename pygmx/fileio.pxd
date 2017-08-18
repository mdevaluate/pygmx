from libc.stdio cimport FILE

from utility cimport *

#cdef extern from "gromacs/fileio/gmx_system_xdr.h":
#    ctypedef struct XDR:
#        pass

cdef extern from "gromacs/fileio/filetypes.h":
    int fn2ftp(const char *fn)

#cdef extern from "gromacs/fileio/gmx_internal_xdr.h":
#    ctypedef struct XDR:
#        pass
#        
#    ctypedef enum xdr_op:
#        XDR_ENCODE = 0
#        XDR_DECODE = 1
#        XDR_FREE   = 2
#
#    void xdrstdio_create (XDR *__xdrs, FILE *__file, xdr_op __xop)

cdef extern from "gromacs/fileio/gmxfio.h":

    ctypedef struct t_fileio:
        pass
        # FILE           *fp               # the file pointer */
        # gmx_bool        bRead             # the file is open for reading */
        # gmx_bool        bDouble           # write doubles instead of floats */
        # gmx_bool        bReadWrite        # the file is open for reading and writing */
        # char        *fn                   # the file name */
        # # XDR         *xdr                  # the xdr data pointer */
        # int  xdrmode              # the xdr mode */
        # int          iFTP                 # the file type identifier */

        # t_fileio    *next, *prev          # next and previous file pointers in the linked list */
        # # tMPI_Lock_t  mtx;                  # content locking mutex. This is a fast lock
    

    void gmx_fio_rewind(t_fileio *fio)

    int gmx_fio_seek(t_fileio *fio, gmx_off_t fpos)

    gmx_off_t gmx_fio_ftell(t_fileio *fio)

#    int xtc_seek_frame(t_fileio *fio, int frame, int natoms)
    int xtc_seek_time(t_fileio *fio, real time, int natoms, gmx_bool bSeekForwardOnly)

    void gmx_fio_setdebug(t_fileio *fio, gmx_bool bDebug)
