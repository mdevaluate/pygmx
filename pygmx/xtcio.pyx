# Wrapper for xtcio.h: I/O for xtc files.

from cpython.array cimport array
from array import array
from xdrlib import Unpacker

import numpy as np
cimport numpy as np
import os

from utility cimport *
from math cimport *
from fileio cimport *
from .gromacs.reader import INDEX_MAGIC, SubscriptableReader, XTCFrame
from .errors import InvalidIndexException, InvalidMagicException, XTCError


cdef extern from "gromacs/fileio/xtcio.h":
    t_fileio *open_xtc(const char *filename, const char *mode)

    void close_xtc(t_fileio *fio)

    int read_first_xtc(t_fileio *fio,
                       int *natoms, gmx_int64_t *step, real *time,
                       matrix box, rvec **x, real *prec, gmx_bool *_bOK)

    int read_next_xtc(t_fileio *fio,
                      int natoms, gmx_int64_t *step, real *time,
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
        int natoms, frame, state = 1
        gmx_int64_t step
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
        bint has_cache, has_times
        rvec *coords
        array _cache, _times
        public str filename

    @property
    def cache(self):
        if self.has_cache:
            return self._cache

    @cache.setter
    def cache(self, indices):
        self._cache = array('L')
        for i in indices:
            self._cache.append(i)
        self.has_cache = True

    @property
    def times(self):
        if self.has_times:
            return self._times

    @times.setter
    def times(self, times):
        self._times = array('f')
        for t in times:
            self._times.append(t)
        self.has_times = True

#    @property
#    def filename(self):
#        return self._filename.decode()

    def seek(self, int frame):
        if self.has_cache:
            gmx_fio_seek(self.fio, self.cache[frame])

    def make_cache(self):
        self._cache = get_xtc_index(self.fio)
        self.has_cache = True

    def load_cache(self, indexfile):

        xtc_stat = os.stat(self.filename)
        c_time = int(xtc_stat.st_ctime)
        m_time = int(xtc_stat.st_mtime)
        size = xtc_stat.st_size

        with open(indexfile, 'rb') as fd:
            unpacker = Unpacker(SubscriptableReader(fd))

            if unpacker.unpack_hyper() != INDEX_MAGIC:
                raise InvalidMagicException
            if unpacker.unpack_hyper() != c_time:
                raise InvalidIndexException
            if unpacker.unpack_hyper() != m_time:
                raise InvalidIndexException
            if unpacker.unpack_hyper() != size:
                raise InvalidIndexException

            self._cache = array('L')
            self._times = array('f')

            try:
                while True:
                    self._cache.append(unpacker.unpack_hyper())
                    self._times.append(unpacker.unpack_float())
            except EOFError:
                pass
        self.has_cache = True

    def __cinit__(self):
        self.has_cache = False
        self.has_times = False

    def __init__(self, filename, indexfile=None, make_cache=False):
        if isinstance(filename, str):
            self.filename = filename
            filename = filename.encode()
        else:
            self.filename = filename.decode()

        cdef:
            gmx_int64_t step
            matrix box
            gmx_bool _bOK
            real time, prec
            rvec *x

        if not os.path.exists(filename):
            raise OSError('File not found: {}'.format(filename))
        if filename.decode().split('.')[-1] != 'xtc':
            raise XTCError('File is not of xtc type: {}'.format(filename))

        self.fio = open_xtc(filename, b'r')
        read_first_xtc(self.fio, &self.natoms, &step, &time, box, &x, &prec, &_bOK)

        if indexfile is not None:
            try:
                self.load_cache(indexfile)
            except InvalidIndexException:
                if make_cache:
                    pass
                else:
                    raise

        if make_cache:
            self.make_cache()

    def __len__(self):
        if self.has_cache:
            return len(self.cache)

    def __getitem__(self, frame):
        cdef matrix box
        cdef gmx_bool _bOK
        cdef real time
        cdef np.ndarray[real, ndim=2] coords = np.empty((self.natoms, 3), dtype=np_real)
        if frame < len(self):
            self.seek(frame)
            read_next_xtc(self.fio, self.natoms, <gmx_int64_t *>self.cur_step, &time, box,
                <rvec *>coords.data, &self.prec, &_bOK)
            if _bOK:
                frame = XTCFrame()
                frame._coordinates = coords
                frame.index = self.cur_step
                frame.time = time
                frame.box = box
                return frame
            else:
                raise
        else:
            raise IndexError('Frame {} is out of range for trajectory of length {}.'.format(frame, len(self)))
