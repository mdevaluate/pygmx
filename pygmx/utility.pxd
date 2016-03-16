# C-API in gromacs/utility

#cdef extern from "inttypes.h":
ctypedef int __int64

cdef extern from "gromacs/utility/basedefinitions.h":
    ctypedef int gmx_bool
    ctypedef __int64 gmx_int64_t

cdef extern from "gromacs/utility/real.h":
    ctypedef double real

cdef extern from "gromacs/utility/futil.h":
    ctypedef gmx_int64_t    gmx_off_t
