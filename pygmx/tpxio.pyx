# cython: c_string_type=unicode, c_string_encoding=utf8
# Cython wrapper around tpxio.cpp

from libc cimport stdio

import numpy as np
#cimport numpy as np

from utility cimport *
from math cimport *
from mdtypes cimport *
from topology cimport *


cdef extern from "gromacs/fileio/tpxio.h":
    ctypedef struct t_tpxheader:
      int   bIr        # Non zero if input_rec is present		*/
      int   bBox       # Non zero if a box is present			*/
      int   bTop       # Non zero if a topology is present		*/
      int   bX         # Non zero if coordinates are present		*/
      int   bV         # Non zero if velocities are present		*/
      int   bF         # Non zero if forces are present (no longer
                       #    supported, but retained so old .tpr can be read) */
      int   natoms     # The total number of atoms			*/
      int   ngtc       # The number of temperature coupling groups    */
      real  _lambda     # Current value of lambda			*/
      int   fep_state  # Current value of the alchemical state --


    void read_tpxheader(const char *fn, t_tpxheader *tpx, bint TopOnlyOK, int *version, int *generation)
    int read_tpx(const char *fname,
                 t_inputrec *ir,
                 matrix box,
                 int *natoms,
                 rvec *x,
                 rvec *v,
                 gmx_mtop_t *mtop)


def cstr(instr):
    if isinstance(instr, str):
        instr = instr.encode()
    return instr

cdef atoms_from_topology(gmx_mtop_t *topology):
    cdef:
        t_atoms c_atoms
        int moltype, resind

    atoms = []
    residues = []
    res_id = 0
    for i_molblock in range(topology.nmolblock):
        moltype = topology.molblock[i_molblock].type
        c_atoms = topology.moltype[moltype].atoms
        for n in range(topology.molblock[i_molblock].nmol):
            mol_atoms = []
            res_id += 1
            for i_atom in range(c_atoms.nr):
                resind = c_atoms.atom[i_atom].resind
                resname = c_atoms.resinfo[resind].name[0]
                if resname not in residues:
                    residues.append(resname)
                resid = residues.index(resname) + 1
                mol_atoms.append((
                    res_id,
                    resname,
                    c_atoms.atomname[i_atom][0],
                ))
            atoms += mol_atoms
    return np.array(atoms)


ctypedef object (*atom_func)(t_atom)


cdef atom_charge(t_atom atom):
    return atom.q


cdef atom_mass(t_atom atom):
    return atom.m


cdef index_groups_from_topology(gmx_mtop_t *topology):
    # retrieve the index groups from the topology->groups ?
    pass


cdef open_tpx(const char* filename, t_inputrec *ir, matrix box, int *natoms, gmx_mtop_t *top):
    
    # unterdrücke GMX output durch internen stderr buffer...
    cdef char buffer[stdio.BUFSIZ]
    stdio.setbuf(stdio.stderr, buffer)
    cdef rvec *x, *v
    return_code = read_tpx(filename, ir, box, natoms, x, v, top)

    for i in range(stdio.BUFSIZ):
        buffer[i] = 0
    stdio.fflush(stdio.stderr)
    stdio.fseek(stdio.stderr, 0, stdio.SEEK_END)
    stdio.setbuf(stdio.stderr, NULL)
    return return_code


cdef read_ffparams(gmx_ffparams_t *ffparams, gmx_bool bShowNumbers):
    cdef char buffer[stdio.BUFSIZ]
    stdio.setbuf(stdio.stdout, buffer)
    pr_ffparams(stdio.stdout, 0, '', ffparams, bShowNumbers)
    stdio.fflush(stdio.stderr)
    stdio.fseek(stdio.stderr, 0, stdio.SEEK_END)
    stdio.setbuf(stdio.stderr, NULL)
    return buffer


cdef class TPXReader:
    cdef:
        t_tpxheader header
        t_inputrec input_record
        gmx_mtop_t topology
        real box[3][3]
        readonly int  n_atoms, n_tcouple_groups, n_mol
        readonly char *topology_name

    cdef _map_atoms(self, object (*func)(t_atom)):
        cdef t_atoms atoms
        map = []
        for i_molblock in range(self.topology.nmolblock):
            moltype = self.topology.molblock[i_molblock].type
            nmol = self.topology.molblock[i_molblock].nmol
            atoms = self.topology.moltype[moltype].atoms
            mol_map = []
            for i_atom in range(atoms.nr):
                mol_map.append(func(atoms.atom[i_atom]))
            map += mol_map * nmol
        return map

    property atoms:
        """
        Get a list of tuples describing all atoms in the system:

        Returns:
            List of tuples: (atom_name, residue_id, resiude_name)
        """
        def __get__(self):
            return atoms_from_topology(&self.topology)

    property mass2:
        def __get__(self):
            return np.array(self._map_atoms(atom_mass))

    property mass:
        "Get the masses of all atoms in the system."
        def __get__(self):
            cdef t_atoms atoms
            masses = []
            for i_molblock in range(self.topology.nmolblock):
                moltype = self.topology.molblock[i_molblock].type
                nmol = self.topology.molblock[i_molblock].nmol
                atoms = self.topology.moltype[moltype].atoms
                mol_masses = []
                for i_atom in range(atoms.nr):
                    mol_masses.append(atoms.atom[i_atom].m)
                masses += mol_masses * nmol
            return np.array(masses)

    property charge:
        "Get the partial charge of all atoms in the system."
        def __get__(self):
            cdef t_atoms atoms
            charges = []
            for i_molblock in range(self.topology.nmolblock):
                moltype = self.topology.molblock[i_molblock].type
                nmol = self.topology.molblock[i_molblock].nmol
                atoms = self.topology.moltype[moltype].atoms
                mol_charges = []
                for i_atom in range(atoms.nr):
                    mol_charges.append(atoms.atom[i_atom].q)
                charges += mol_charges * nmol
            return np.array(charges)

    property type:
        "Get the type of all atoms in the system."
        def __get__(self):
            cdef t_atoms atoms
            types = []
            for i_molblock in range(self.topology.nmolblock):
                moltype = self.topology.molblock[i_molblock].type
                nmol = self.topology.molblock[i_molblock].nmol
                atoms = self.topology.moltype[moltype].atoms
                mol_type = []
                for i_atom in range(atoms.nr):
                    mol_type.append(atoms.atom[i_atom].type)
                types += mol_type * nmol
            return np.array(types)

    property molecules:
        "Get molecule indices from topology."
        def __get__(self):
            mols = [0] * self.n_atoms
            molid = 0
            for i in range(self.topology.mols.nr):
                molid += 1
                for j in range(self.topology.mols.index[i], self.topology.mols.index[i+1]):
                    mols[j] = molid
            return mols
    # @property
    # def nsteps(self):
    #     return self.input_record.nsteps

    # @property
    # def nstxout(self):
    #     return self.input_record.nstxout

    # @property
    # def nstvout(self):
    #     return self.input_record.nstvout

    # @property
    # def nstfout(self):
    #     return self.input_record.nstfout

    # @property
    # def nstxout_compressed(self):
    #     return self.input_record.nstxout_compressed

    # @property
    # def nstenergy(self):
    #     return self.input_record.nstenergy


    def read_ff(self):
        return read_ffparams(&self.topology.ffparams, True)

    def __cinit__(self, filename):
        filename = cstr(filename)

        # Arrays for coordinates etc can be NULL, if they shouldn't be read
        # cdef int version, generation
        # read_tpxheader(<char *>filename, &self.header, True, &version, &generation)
        # self.n_atoms = self.header.natoms
        # self.n_tcouple_groups = self.header.ngtc
        # cdef np.ndarray[real, ndim=2] coordinates = np.empty((self.n_atoms, 3), dtype=np.float32)
        # cdef np.ndarray[real, ndim=2] velocites = np.empty((self.n_atoms, 3), dtype=np.float32)
        # cdef np.ndarray[real, ndim=2] forces = np.empty((self.n_atoms, 3), dtype=np.float32)
        open_tpx(
            <char *>filename,
            NULL,
            self.box,
            &self.n_atoms,
            &self.topology
        )
        self.topology_name = self.topology.name[0]


@cython.binding(True)
def make_xtcframe_whole(coords, box, <TPXReader>reader):
    cdef t_atoms = gmx_mtop_global_atoms(reader.topology)
    cdef np.ndarray[real, ndim=2] b = np.asarray(box, dtype=np.float32)
    cdef np.ndarray[real, ndim=2] x = np.asarray(coords, dtype=np.float32)
    rm_gropbc(const t_atoms *atoms, <rvec *>x.data, <matrix> b)
    return x