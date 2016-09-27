
from cpython.array cimport array
from array import array

cimport numpy as np
import numpy as np

from utility cimport *
from mdtypes cimport *

cdef extern from "gromacs/fileio/enxio.h":
    ctypedef struct gmx_enxnm_t:
        char *name
        char *unit

    ctypedef struct ener_file_t:
        pass

    ctypedef struct t_enxblock:
        pass

    ctypedef struct t_enxframe:
        double          t            # Timestamp of this frame	                     */
        gmx_int64_t     step         # MD step	                             */
        gmx_int64_t     nsteps       # The number of steps between frames            */
        double          dt           # The MD time step                              */
        int             nsum         # The number of terms for the sums in ener      */
        int             nre          # Number of energies			     */
        int             e_size       # Size (in bytes) of energies		     */
        int             e_alloc      # Allocated size (in elements) of ener          */
        t_energy       *ener         # The energies                                  */
        int             nblock       # Number of following energy blocks             */
        t_enxblock     *block        # The blocks                                    */
        int             nblock_alloc # The number of blocks allocated                */

    ener_file_t open_enx(const char *fn, const char *mode)
    void close_enx(ener_file_t ef)

    void do_enxnms(ener_file_t ef, int *nre, gmx_enxnm_t **enms)

    gmx_bool do_enx(ener_file_t ef, t_enxframe *fr)

    void free_enxframe(t_enxframe *ef)


cdef class EDRFile:

    cdef:
        ener_file_t efile
        gmx_enxnm_t *etypes
        int n_etypes
        t_enxframe *frames

    @property
    def types(self):

        types = []
        for i in range(self.n_etypes):
            types.append((self.etypes[i].name.decode(), self.etypes[i].unit.decode()))
        return types

    def read(self, steps=None):
        cdef:
            t_enxframe *frame
            np.ndarray[float, ndim=2] energies
            np.ndarray[float, ndim=1] times
        if steps is None:
            energies = np.empty((1,self.n_etypes), np.float32)
            times = np.empty((1,), np.float32)
        else:
            energies = np.empty((steps,self.n_etypes), np.float32)
            times = np.empty((steps,), np.float32)
        snew(frame, 1)

        i = 0
        while True:
            cont = do_enx(self.efile, frame)

            for j in range(frame.nre):
                energies[i, j] = frame.ener[j].e
            times[i] = frame.t

            if cont:
                i += 1
                if len(energies) <= i:
                    energies = np.resize(energies, (len(energies)+1, self.n_etypes))
                    times = np.resize(times, (len(times)+1,))
            else:
                break

        free_enxframe(frame)
        sfree(frame)

        return times, energies


    def __init__(self, filename):
        filename = cstr(filename)

        self.efile = open_enx(filename, b'r')

        self.etypes = NULL
        do_enxnms(self.efile, &self.n_etypes, &self.etypes)
        print(self.n_etypes)

        #do_enx(self.efile, self.frames)


    def __getitem__(self, item):
        cdef:
            int f = 0
            t_enxframe frame = self.frames[f]

        cdef array values = array('d')

        while frame.nblock > 0:
            if item < frame.nre:
                values.append(frame.ener[item].e)
            i += 1
            frame = self.frames[i]
            print(i, frame.t, frame.nblock)

        return values
