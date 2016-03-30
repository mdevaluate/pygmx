# Wrapper for xtcio.h: I/O for xtc files.

from cpython.array cimport array
from array import array

import numpy as np
cimport numpy as np

from utility cimport *
from math cimport *
from fileio cimport *


cdef extern from "gromacs/fileio/xtcio.h":
    t_fileio *open_xtc(const char *filename, const char *mode)

    void close_xtc(t_fileio *fio)

    int read_first_xtc(t_fileio *fio,
                       int *natoms, int *step, real *time,
                       matrix box, rvec **x, real *prec, gmx_bool *_bOK)

    int read_next_xtc(t_fileio *fio,
                      int natoms, int *step, real *time,
                      matrix box, rvec *x, real *prec, gmx_bool *_bOK)


if sizeof(real) == 4:
    np_real = np.float32
else:
    np_real = np.float


#cdef array get_xtc_index_by_frames(t_fileio *fio, int length, int natoms):
#    cdef:
        #gmx_bool _bOK
##        int frame = 1
#    cdef array cache = array('L')
#    gmx_fio_rewind(fio)
#    cache.append(gmx_fio_ftell(fio))
#    while xtc_seek_frame(fio, frame*500, natoms) == 0:
#        cache.append(gmx_fio_ftell(fio))
#        #print(frame, cache[-1])
#        frame += 1
#        if frame == length:
#            break
#    return cache

cdef array get_xtc_index(t_fileio *fio):
    cdef:
        gmx_bool _bOK
        int natoms, step, frame, state = 1
        real time, prec
        matrix box
        rvec *x
    cdef array cache = array('L')
    gmx_fio_rewind(fio)
    cache.append(gmx_fio_ftell(fio))
    read_first_xtc(fio, &natoms, &step, &time, box, &x, &prec, &_bOK)
    cache.append(gmx_fio_ftell(fio))
    while read_next_xtc(fio, natoms, &step, &time, box, x, &prec, &_bOK):
        cache.append(gmx_fio_ftell(fio))
    # the last index is invalid
    return cache[:-1]


cdef class XTCReader:
    cdef:
        t_fileio *fio
        int natoms, cur_step
        real start_time, timestep, prec, cur_time
        bint has_cache
        rvec *coords
        array _cache

    @property
    def cache(self):
        return self._cache

    def seek(self, int frame):
        if self.has_cache:
            gmx_fio_seek(self.fio, self.cache[frame])


    def make_cache(self):
        self._cache = get_xtc_index(self.fio)
        self.has_cache = True

    def __cinit__(self):
        self.has_cache = False

    def __init__(self, filename, make_cache=True):
        if isinstance(filename, str):
            filename = filename.encode()

        cdef:
            int step
            matrix box
            gmx_bool _bOK
            real time, prec
            rvec *x

        self.fio = open_xtc(filename, b'r')
        read_first_xtc(self.fio, &self.natoms, &step, &time, box, &x, &prec, &_bOK)

        if make_cache:
            self.make_cache()

    def __len__(self):
        return len(self.cache)

    def __getitem__(self, frame):
        cdef matrix box
        cdef gmx_bool _bOK
        cdef real time
        cdef np.ndarray[real, ndim=2] coords = np.empty((self.natoms, 3), dtype=np_real)
        if frame < len(self):
            self.seek(frame)
            read_next_xtc(self.fio, self.natoms, &self.cur_step, &time, box,
                <rvec *>coords.data, &self.prec, &_bOK)
            if _bOK:
                return coords
            else:
                raise
        else:
            raise IndexError('Frame {} is out of range for trajectory of length {}.'.format(frame, len(self)))
