
from utility cimport real

cdef extern from "gromacs/math/vectypes.h":
    ctypedef real   rvec[3]
    ctypedef real   matrix[3][3]
    ctypedef real   tensor[3][3]
