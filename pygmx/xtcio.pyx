# Wrapper for xtcio.h: I/O for xtc files.

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


cdef class XTCReader:
    cdef:
        t_fileio *fio
        int natoms, cur_step
        real start_time, timestep, prec, cur_time
        dict _cache
        #rvec *coords

    def nearest_cache(self, index):
        nearest = 0
        cur_dist = 99999999
        for ind in self._cache:
            dist = abs(ind - index)
            if dist < cur_dist:
                nearest = ind
                cur_dist = dist
        return nearest


    def seek(self, index):
        if index in self._cache:
            gmx_fio_seek(self.fio, self._cache[index])
        else:
            nearest = self.nearest_cache(index)
            gmx_fio_seek(self.fio, self._cache[nearest])
            self.cur_step = index
            self.cur_time = index * self.timestep + self.start_time
            xtc_seek_time(self.fio, self.cur_time, self.natoms, index > nearest)
            self._cache[index] = gmx_fio_ftell(self.fio)


    def __init__(self, filename):
        if isinstance(filename, str):
            filename = filename.encode()

        cdef:
            matrix box
            gmx_bool _bOK
            real time
        cdef rvec *coords # = np.empty((1,3), dtype=np.float32)

        self._cache = {}
        self.fio = open_xtc(filename, b'r')
        read_first_xtc(self.fio, &self.natoms, &self.cur_step, &self.start_time, box,
            &coords, &self.prec, &_bOK)
        self._cache[self.cur_step] = gmx_fio_ftell(self.fio)

        read_next_xtc(self.fio, self.natoms, &self.cur_step, &time, box,
            coords, &self.prec, &_bOK)

        self._cache[<int>self.cur_step] = gmx_fio_ftell(self.fio)
        self.timestep = time - self.start_time


    def __getitem__(self, item):
        cdef matrix box
        cdef gmx_bool _bOK
        cdef real time
        cdef np.ndarray[real, ndim=2] coords = np.empty((self.natoms, 3), dtype=np_real)

        self.seek(item)

        read_next_xtc(self.fio, self.natoms, &self.cur_step, &time, box,
            <rvec *>coords.data, &self.prec, &_bOK)

        return coords
