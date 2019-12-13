from cpython.array cimport array
import cython
from array import array
from xdrlib import Unpacker, Packer

import numpy as np
cimport numpy as np
import os

from utility cimport *
from math cimport *
from fileio cimport *
from .gromacs.reader import INDEX_MAGIC, SubscriptableReader, XTCFrame
from .errors import InvalidIndexException, InvalidMagicException, XTCError


cdef extern from "gromacs/fileio/trrio.h":

    ctypedef struct gmx_trr_header_t:
        gmx_bool    bDouble   #/* Double precision?                   */
        int         ir_size   #/* Backward compatibility              */
        int         e_size    #/* Backward compatibility              */
        int         box_size  #/* Non zero if a box is present        */
        int         vir_size  #/* Backward compatibility              */
        int         pres_size #/* Backward compatibility              */
        int         top_size  #/* Backward compatibility              */
        int         sym_size  #/* Backward compatibility              */
        int         x_size    #/* Non zero if coordinates are present */
        int         v_size    #/* Non zero if velocities are present  */
        int         f_size    #/* Non zero if forces are present      */

        int         natoms    #/* The total number of atoms           */
        gmx_int64_t step      #/* Current step number                 */
        int         nre       #/* Backward compatibility              */
        real        t         #/* Current time                        */
        real        lamb    #/* Current value of lambda             */
        int         fep_state #/* Current value of alchemical state   */

    t_fileio *gmx_trr_open(const char *fn, const char *mode)
    void gmx_trr_close(t_fileio *fio)

    gmx_bool gmx_trr_read_frame_header(t_fileio *fio, gmx_trr_header_t *header, gmx_bool *bOK)
    gmx_bool gmx_trr_read_frame(t_fileio *fio, gmx_int64_t *step, real *t, real *l,
                                rvec *box, int *natoms, rvec *x, rvec *v, rvec *f)


if sizeof(real) == 4:
    np_real = np.float32
else:
    np_real = np.float



cdef tuple get_trr_index(t_fileio *fio):
    cdef:
        gmx_bool bOK

        gmx_trr_header_t header

        int natoms, frame, state = 1
        gmx_int64_t step
        real time, prec, lamb
        matrix box
        rvec *x
    cdef array cache = array('L')
    cdef array times = array('f')
    gmx_fio_rewind(fio)
    cache.append(gmx_fio_ftell(fio))
    #gmx_trr_read_frame_header(fio, &header, &bOK)
    #cache.append(gmx_fio_ftell(fio))
    while gmx_trr_read_frame(fio, &step, &time, &lamb, box, &natoms, NULL, NULL, NULL):
        cache.append(gmx_fio_ftell(fio))
        times.append(time)
    # the last index is invalid
    return cache[:-1], times[:-1]


class TRRFrame:

    def __init__(self, index, time, box, positions, velocities=None, forces=None):
        self.index= index
        self.time = time
        self.box = box
        self.positions = positions
        self.velocities = velocities
        self.forces =forces


cdef class TRRReader:
    cdef:
        t_fileio *fio
        public str filename
        gmx_trr_header_t header
        int natoms
        gmx_int64_t cur_step
        real start_time, timestep, prec, cur_time
        bint has_cache, has_times
        array _cache, _times
    
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

    def make_cache(self):
        self._cache, self._times = get_trr_index(self.fio)
        self.has_cache = self.has_times = True

    def seek(self, int frame):
        if self.has_cache:
            gmx_fio_seek(self.fio, self.cache[frame])

    def load_cache(self, indexfile, ignore_time=False):

        xtc_stat = os.stat(self.filename)
        c_time = int(xtc_stat.st_ctime)
        m_time = int(xtc_stat.st_mtime)
        size = xtc_stat.st_size

        with open(indexfile, 'rb') as fd:
            unpacker = Unpacker(fd.read())
            # first 4 * 8 bytes are used for checks, followed by N * (8 + 4) bytes of data
            length = int((len(unpacker.get_buffer()) - 32) / 12)

        if length < 0:
            raise InvalidIndexException

        if unpacker.unpack_hyper() != INDEX_MAGIC:
            raise InvalidMagicException
        if unpacker.unpack_hyper() != c_time and not ignore_time:
            raise InvalidIndexException
        if unpacker.unpack_hyper() != m_time and not ignore_time:
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
        self.has_times = True

    def write_cache(self, index_filename):
        if not self.has_cache:
            raise Exception("Can not write indexfile, no cache loaded.")

        packer = Packer()
        xtc_stat = os.stat(self.filename)
        c_time = int(xtc_stat.st_ctime)
        m_time = int(xtc_stat.st_mtime)
        size = xtc_stat.st_size

        packer.pack_hyper(INDEX_MAGIC)
        # TODO Maybe store only the hashsum of the file for identification?
        packer.pack_hyper(c_time)
        packer.pack_hyper(m_time)
        packer.pack_hyper(size)

        for i, t in zip(self.cache, self.times):
            packer.pack_hyper(i)
            packer.pack_float(t)

        with open(index_filename, 'wb') as idx_fd:
            idx_fd.write(packer.get_buffer())
            packer.reset()
        

    def __cinit__(self):
        self.has_cache = False
        self.has_times = False

    def __init__(self, filename, indexfile=None, make_cache=False, ignore_timestamps=False):

        self.filename = str(filename)
        self.fio = gmx_trr_open(self.filename.encode(), b'r')

        cdef gmx_bool bOK 
        gmx_trr_read_frame_header(self.fio, &self.header, &bOK)
        self.natoms = self.header.natoms

        if make_cache:
            self.make_cache()

        if indexfile is not None:
            if os.path.exists(indexfile):
                try:
                    self.load_cache(indexfile, ignore_time=ignore_timestamps)
                except (InvalidIndexException, InvalidMagicException):
                    self.make_cache()
                    self.write_cache(indexfile)
            else:
                print("Creating Indexfile for TRR file.")
                self.make_cache()
                self.write_cache(indexfile)

    def __len__(self):
        if self.has_cache:
            return len(self.cache)

    def __getitem__(self, frame):
        cdef matrix box
        cdef gmx_bool _bOK
        cdef real time, lamb
        cdef np.ndarray[real, ndim=2] x = np.empty((self.natoms, 3), dtype=np_real)
        cdef np.ndarray[real, ndim=2] v = np.empty((self.natoms, 3), dtype=np_real)
        cdef np.ndarray[real, ndim=2] f = np.empty((self.natoms, 3), dtype=np_real)
        if frame < len(self):
            self.seek(frame)

            _bOK = gmx_trr_read_frame(self.fio, &self.cur_step, &time, &lamb,
                                box, &self.natoms, <rvec *>x.data, <rvec *>v.data, <rvec *>f.data)

            if _bOK:

                frame = TRRFrame(self.cur_step, time, box, x, v, f)
                return frame
            else:
                raise
        else:
            raise IndexError('Frame {} is out of range for trajectory of length {}.'.format(frame, len(self)))



    