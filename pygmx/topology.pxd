
from libc.stdio cimport FILE

from utility cimport *
from math cimport *
from mdtypes cimport *

cdef extern from "gromacs/topology/atoms.h":
    ctypedef struct t_atom:
        real           m, q        # Mass and charge                      */
        real           mB, qB      # Mass and charge for Free Energy calc */
        unsigned short type        # Atom type                            */
        unsigned short typeB       # Atom type for Free Energy calc       */
        int            ptype       # Particle type                        */
        int            resind      # Index into resinfo (in t_atoms)      */
        int            atomnumber  # Atomic Number or 0                   */
        char           elem[4]     # Element name                         */

    ctypedef struct t_resinfo:
        char          **name       # Pointer to the residue name          */
        int             nr         # Residue number                       */
        unsigned char   ic         # Code for insertion of residues       */
        int             chainnum   # Iincremented at TER or new chain id  */
        char            chainid    # Chain identifier written/read to pdb */
        char          **rtp        # rtp building block name (optional)   */

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
        int   nr                   # Number of different groups           */
        int  *nm_ind               # Index in the group names             */

#ctypedef t_atoms *t_atoms_ptr

cdef extern from "gromacs/topology/symtab.h":
  ctypedef struct t_symtab:
    pass

cdef extern from "gromacs/topology/block.h":
    ctypedef struct t_block:
        int      nr           #/* The number of blocks          */
        int     *index        #/* Array of indices (dim: nr+1)  */
        int      nalloc_index #/* The allocation size for index */

    ctypedef struct t_blocka:
        pass

cdef extern from "gromacs/topology/idef.h":
    ctypedef struct t_ilist:
        pass

    ctypedef struct t_idef:
        pass

    #cdef enum t_ft_enum:
    #    F_LJ

    ctypedef union t_iparams:
        pass

    ctypedef int t_functype

    ctypedef struct gmx_cmap_t:
        pass

    ctypedef struct gmx_ffparams_t:
        int         ntypes
        int         atnr
        t_functype *functype
        t_iparams  *iparams
        double      reppow    # The repulsion power for VdW: C12*r^-reppow   */
        real        fudgeQQ   # The scaling factor for Coulomb 1-4: f*q1*q2  */
        gmx_cmap_t  cmap_grid # The dihedral correction maps                 */


    void pr_ffparams(FILE *fp, int indent, const char *title,
                     const gmx_ffparams_t *ffparams, gmx_bool bShowNumbers)

cdef extern from "gromacs/topology/topology.h":
    ctypedef struct gmx_moltype_t:
        char          **name         # Name of the molecule type            */
        t_atoms         atoms        # The atoms in this molecule           */
        t_ilist         ilist[0] # Interaction list with local indices  */
        t_block         cgs          # The charge groups                    */
        t_blocka        excls        # The exclusions                       */

    ctypedef struct gmx_molblock_t:
        int     type               #*< The molecule type index in mtop.moltype  */
        int     nmol               #*< The number of molecules in this block    */
        int     nposres_xA         #*< The number of posres coords for top A    */
        rvec   *posres_xA          #*< Position restraint coordinates for top A */
        int     nposres_xB         #*< The number of posres coords for top B    */
        rvec   *posres_xB          #*< Position restraint coordinates for top B */

        # Convenience information, derived from other gmx_mtop_t contents     */
        int     natoms_mol         #*< The number of atoms in one molecule      */
        int     globalAtomStart    #*< Global atom index of the first atom in the block */
        int     globalAtomEnd      #*< Global atom index + 1 of the last atom in the block */
        int     globalResidueStart #*< Global residue index of the first residue in the block */
        int     residueNumberStart #*< Residue numbers start from this value if the number of residues per molecule is <= maxres_renum */

    ctypedef struct gmx_groups_t:
        pass
#        t_grps            grps[0]   # Groups of things                     */
#        int               ngrpname      # Number of groupnames                 */
#        char           ***grpname       # Names of the groups                  */
#        int               ngrpnr[0]
#        unsigned char    *grpnr[0]  # Group numbers or NULL                */

    ctypedef struct gmx_mtop_t:
        char           **name      # Name of the topology                 */
        gmx_ffparams_t   ffparams
        int              nmoltype
        gmx_moltype_t   *moltype
        int              nmolblock
        gmx_molblock_t  *molblock
        gmx_bool         bIntermolecularInteractions # Are there intermolecular
                                                     #  * interactions?            */
        t_ilist         *intermolecular_ilist        # List of intermolecular interactions
                                                      # * using system wide atom indices,
                                                      # * either NULL or size F_NRE           */
        int              natoms
        int              maxres_renum                # Parameter for residue numbering      */
        int              maxresnr                    # The maximum residue number in moltype */
        t_atomtypes      atomtypes                   # Atomtype properties                  */
        t_block          mols                        # The molecules                        */
        gmx_groups_t     groups
        t_symtab         symtab                      # The symbol table                     */

    ctypedef struct t_topology:
        char          **name                       # /* Name of the topology                 */
        t_idef          idef                       # /* The interaction function definition  */
        t_atoms         atoms                      # /* The atoms                            */
        t_atomtypes     atomtypes                  # /* Atomtype properties                  */
        t_block         cgs                        # /* The charge groups                    */
        t_block         mols                       # /* The molecules                        */
        gmx_bool        bIntermolecularInteractions# /* Inter.mol. int. ?   */
        t_blocka        excls                      # /* The exclusions                       */
        t_symtab        symtab                     # /* The symbol table                     */

    ctypedef struct gmx_localtop_t:
        t_idef        idef        # /* The interaction function definition  */
        #t_atomtypes   atomtypes   # /* Atomtype properties                  */
        #t_block       cgs         # /* The charge groups                    */
        #   t_blocka      excls       # /* The exclusions                       */

    void init_mtop(gmx_mtop_t *mtop)
    void done_top(t_topology *top)

# cdef extern from "gromacs/topology/topology.h":
        # generate a t_atoms struct for the system from gmx_mtop_t
        # t_atoms* mtop2atoms(gmx_mtop_t *mtop)
cdef extern from "gromacs/topology/mtop_util.h":
    t_atoms gmx_mtop_global_atoms(const gmx_mtop_t *mtop)
    t_topology gmx_mtop_t_to_t_topology(gmx_mtop_t *mtop, bint freeMTop)
    gmx_localtop_t *gmx_mtop_generate_local_top(const gmx_mtop_t *mtop, bool freeEnergyInteractionsAtEnd)



cdef extern from "gromacs/pbcutil/rmpbc.h":
    ctypedef struct gmx_rmpbc_t:
        pass

    gmx_rmpbc_t gmx_rmpbc_init(const t_idef *idef, int ePBC, int natoms)
    void gmx_rmpbc_done(gmx_rmpbc_t gpbc)
    void rm_gropbc(const t_atoms *atoms, rvec x[], const matrix box)
    void gmx_rmpbc(gmx_rmpbc_t gpbc, int natoms, const matrix box, rvec x[])




