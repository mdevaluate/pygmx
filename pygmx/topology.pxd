
from utility cimport *
from math cimport *
from mdtypes cimport *

cdef extern from "gromacs/topology/atoms.h":
    ctypedef struct t_atom:
        real           m, q;        # Mass and charge                      */
        real           mB, qB;      # Mass and charge for Free Energy calc */
        unsigned short type;        # Atom type                            */
        unsigned short typeB;       # Atom type for Free Energy calc       */
        int            ptype;       # Particle type                        */
        int            resind;      # Index into resinfo (in t_atoms)      */
        int            atomnumber;  # Atomic Number or 0                   */
        char           elem[4];     # Element name                         */

    ctypedef struct t_resinfo:
        char          **name;       # Pointer to the residue name          */
        int             nr;         # Residue number                       */
        unsigned char   ic;         # Code for insertion of residues       */
        int             chainnum;   # Iincremented at TER or new chain id  */
        char            chainid;    # Chain identifier written/read to pdb */
        char          **rtp;        # rtp building block name (optional)   */

    ctypedef struct t_pdbinfo:
        pass

    ctypedef struct t_atoms:
        int            nr          # Nr of atoms                          */
        t_atom        *atom        # Array of atoms (dim: nr)             */
                                    # The following entries will not       */
                                    # always be used (nres==0)             */
        char          ***atomname  # Array of pointers to atom name       */
                                    # use: (*(atomname[i]))                */
        char          ***atomtype  # Array of pointers to atom types      */
                                    # use: (*(atomtype[i]))                */
        char          ***atomtypeB # Array of pointers to B atom types    */
                                    # use: (*(atomtypeB[i]))               */
        int              nres      # The number of resinfo entries        */
        t_resinfo       *resinfo   # Array of residue names and numbers   */
        t_pdbinfo       *pdbinfo   # PDB Information, such as aniso. Bfac */

    ctypedef struct t_atomtypes:
        pass

    ctypedef struct t_grps:
        int   nr;                   # Number of different groups           */
        int  *nm_ind;               # Index in the group names             */

#ctypedef t_atoms *t_atoms_ptr

cdef extern from "gromacs/topology/symtab.h":
  ctypedef struct t_symtab:
    pass

cdef extern from "gromacs/topology/block.h":
    ctypedef struct t_block:
        pass
    ctypedef struct t_blocka:
        pass

cdef extern from "gromacs/topology/idef.h":
  ctypedef struct t_ilist:
    pass

  ctypedef struct gmx_ffparams_t:
    pass

cdef extern from "gromacs/topology/topology.h":
    ctypedef struct gmx_moltype_t:
        char          **name;         # Name of the molecule type            */
        t_atoms         atoms;        # The atoms in this molecule           */
        t_ilist         ilist[0]; # Interaction list with local indices  */
        t_block         cgs;          # The charge groups                    */
        t_blocka        excls;        # The exclusions                       */

    ctypedef struct gmx_molblock_t:
        int            type;        # The molcule type index in mtop.moltype */
        int            nmol;        # The number of molecules in this block  */
        int            natoms_mol;  # The number of atoms in one molecule    */
        int            nposres_xA;  # The number of posres coords for top A  */
        rvec          *posres_xA;   # The posres coords for top A            */
        int            nposres_xB;  # The number of posres coords for top B  */
        rvec          *posres_xB;   # The posres coords for top B            */

    ctypedef struct gmx_groups_t:
        pass
#        t_grps            grps[0]   # Groups of things                     */
#        int               ngrpname      # Number of groupnames                 */
#        char           ***grpname       # Names of the groups                  */
#        int               ngrpnr[0]
#        unsigned char    *grpnr[0]  # Group numbers or NULL                */

    ctypedef struct gmx_mtop_t:
        char           **name       # Name of the topology                 */
        gmx_ffparams_t   ffparams
        int              nmoltype
        gmx_moltype_t   *moltype
        int              nmolblock
        gmx_molblock_t  *molblock
        bint             bIntermolecularInteractions    # Are there intermolecular
                                                        # interactions?            */
        t_ilist         *intermolecular_ilist           # List of intermolecular interactions
                                                        # using system wide atom indices,
                                                        # either NULL or size F_NRE           */
        int              natoms
        int              maxres_renum                 # Parameter for residue numbering      */
        int              maxresnr                     # The maximum residue number in moltype */
        t_atomtypes      atomtypes                    # Atomtype properties                  */
        t_block          mols                         # The molecules                        */
        gmx_groups_t     groups
        t_symtab         symtab                       # The symbol table                     */

# cdef extern from "gromacs/topology/topology.h":
        # generate a t_atoms struct for the system from gmx_mtop_t
        # t_atoms* mtop2atoms(gmx_mtop_t *mtop)
