# C-API in gromacs/utility

from libc.stdint cimport int64_t


cdef extern from "gromacs/utility/basedefinitions.h":
    ctypedef int gmx_bool
    ctypedef int64_t gmx_int64_t

cdef extern from "gromacs/utility/real.h":
    ctypedef double real

cdef extern from "gromacs/utility/futil.h":
    ctypedef gmx_int64_t    gmx_off_t

cdef extern from "gromacs/utility/smalloc.h":
    void snew(void *ptr, int nelem)
    void sfree(void *ptr)


cdef inline cstr(instr):
    if isinstance(instr, str):
        instr = instr.encode()
    return instr
