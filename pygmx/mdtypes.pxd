
from utility cimport *
from math cimport *

cdef extern from "gromacs/trajectory/energy.h":
    ctypedef struct t_energy:
        real    e
        double  eav
        double  esum


cdef extern from "gromacs/mdtypes/inputrec.h":
    ctypedef struct t_simtemp:
        pass
    ctypedef struct t_lambda:
        pass
    ctypedef struct t_expanded:
        pass
    ctypedef struct t_rot:
        pass
    ctypedef struct t_IMD:
        pass
    ctypedef struct t_grpopts:
        int       ngtc
        real     *ref_t          # Coupling temperature	per group   */

    ctypedef struct t_cosines:
        pass
    ctypedef struct t_swapcoords:
        pass

    ctypedef struct t_inputrec:
        # t_inputrec()
        # explicit t_inputrec(const t_inputrec &) = delete
        # t_inputrec &operator=(const t_inputrec &) = delete
        # ~t_inputrec()
        int             eI                      # Integration method                 */
        gmx_int64_t     nsteps                  # number of steps to be taken         */
        int             simulation_part         # Used in checkpointing to separate chunks */
        gmx_int64_t     init_step               # start at a stepcount >0 (used w. convert-tpr)    */
        int             nstcalcenergy           # frequency of energy calc. and T/P coupl. upd.   */
        int             cutoff_scheme           # group or verlet cutoffs     */
        int             ns_type                 # which ns method should we use?               */
        int             nstlist                 # number of steps before pairlist is generated    */
        int             ndelta                  # number of cells per rlong           */
        int             nstcomm                 # number of steps after which center of mass  */
                                                 # motion is removed               */
        int             comm_mode               # Center of mass motion removal algorithm      */
        int             nstlog                  # number of steps after which print to logfile    */
        int             nstxout                 # number of steps after which X is output */
        int             nstvout                 # id. for V                   */
        int             nstfout                 # id. for F                   */
        int             nstenergy               # number of steps after which energies printed */
        int             nstxout_compressed      # id. for compressed trj (.xtc,.tng)           */
        double          init_t                  # initial time (ps)              */
        double          delta_t                 # time step (ps)              */
        real            x_compression_precision # precision of x in compressed trajectory file */
        real            fourier_spacing         # requested fourier_spacing, when nk? not set  */
        int             nkx, nky, nkz           # number of k vectors in each spatial dimension*/
                                                 # for fourier methods for long range electrost.*/
        int             pme_order               # interpolation order for PME                  */
        real            ewald_rtol              # Real space tolerance for Ewald, determines   */
                                                 # the real/reciprocal space relative weight    */
        real            ewald_rtol_lj           # Real space tolerance for LJ-Ewald            */
        int             ewald_geometry          # normal/3d ewald, or pseudo-2d LR corrections */
        real            epsilon_surface         # Epsilon for PME dipole correction            */
        int             ljpme_combination_rule  # Type of combination rule in LJ-PME          */
        int             ePBC                    # Type of periodic boundary conditions        */
        int             bPeriodicMols           # Periodic molecules                           */
        gmx_bool        bContinuation           # Continuation run: starting state is correct */
        int             etc                     # temperature coupling               */
        int             nsttcouple              # interval in steps for temperature coupling   */
        gmx_bool        bPrintNHChains          # whether to print nose-hoover chains        */
        int             epc                     # pressure coupling                            */
        int             epct                    # pressure coupling type          */
        int             nstpcouple              # interval in steps for pressure coupling      */
        real            tau_p                   # pressure coupling time (ps)         */
        tensor          ref_p                   # reference pressure (kJ/(mol nm^3))      */
        tensor          compress                # compressability ((mol nm^3)/kJ)        */
        int             refcoord_scaling        # How to scale absolute reference coordinates  */
        rvec            posres_com              # The COM of the posres atoms                  */
        rvec            posres_comB             # The B-state COM of the posres atoms          */
        int             andersen_seed           # Random seed for Andersen thermostat (obsolete) */
        real            verletbuf_tol           # Per atom pair energy drift tolerance (kJ/mol/ps/atom) for list buffer  */
        real            rlist                   # short range pairlist cut-off (nm)       */
        real            rtpi                    # Radius for test particle insertion           */
        int             coulombtype             # Type of electrostatics treatment             */
        int             coulomb_modifier        # Modify the Coulomb interaction              */
        real            rcoulomb_switch         # Coulomb switch range start (nm)     */
        real            rcoulomb                # Coulomb cutoff (nm)                     */
        real            epsilon_r               # relative dielectric constant                 */
        real            epsilon_rf              # relative dielectric constant of the RF       */
        int             implicit_solvent        # No (=explicit water), or GBSA solvent models */
        int             gb_algorithm            # Algorithm to use for calculation Born radii  */
        int             nstgbradii              # Frequency of updating Generalized Born radii */
        real            rgbradii                # Cutoff for GB radii calculation              */
        real            gb_saltconc             # Salt concentration (M) for GBSA models       */
        real            gb_epsilon_solvent      # dielectric coeff. of implicit solvent     */
        real            gb_obc_alpha            # 1st scaling factor for Bashford-Case GB      */
        real            gb_obc_beta             # 2nd scaling factor for Bashford-Case GB      */
        real            gb_obc_gamma            # 3rd scaling factor for Bashford-Case GB      */
        real            gb_dielectric_offset    # Dielectric offset for Still/HCT/OBC     */
        int             sa_algorithm            # Algorithm for SA part of GBSA                */
        real            sa_surface_tension      # Energy factor for SA part of GBSA */
        int             vdwtype                 # Type of Van der Waals treatment              */
        int             vdw_modifier            # Modify the VdW interaction                   */
        real            rvdw_switch             # Van der Waals switch range start (nm)        */
        real            rvdw                    # Van der Waals cutoff (nm)           */
        int             eDispCorr               # Perform Long range dispersion corrections    */
        real            tabext                  # Extension of the table beyond the cut-off,   *
                                                #  * as well as the table length for 1-4 interac. */
        real            shake_tol               # tolerance for shake             */
        int             efep                    # free energy calculations                     */
        t_lambda       *fepvals                 # Data for the FEP state                       */
        gmx_bool        bSimTemp                # Whether to do simulated tempering            */
        t_simtemp      *simtempvals             # Variables for simulated tempering            */
        gmx_bool        bExpanded               # Whether expanded ensembles are used          */
        t_expanded     *expandedvals            # Expanded ensemble parameters              */
        int             eDisre                  # Type of distance restraining                 */
        real            dr_fc                   # force constant for ta_disre         */
        int             eDisreWeighting         # type of weighting of pairs in one restraints    */
        gmx_bool        bDisreMixed             # Use comb of time averaged and instan. viol's    */
        int             nstdisreout             # frequency of writing pair distances to enx   */
        real            dr_tau                  # time constant for memory function in disres    */
        real            orires_fc               # force constant for orientational restraints  */
        real            orires_tau              # time constant for memory function in orires    */
        int             nstorireout             # frequency of writing tr(SD) to enx           */
        real            em_stepsize             # The stepsize for updating           */
        real            em_tol                  # The tolerance               */
        int             niter                   # Number of iterations for convergence of      */
                                                 # steepest descent in relax_shells             */
        real            fc_stepsize             # Stepsize for directional minimization        */
                                                 # in relax_shells                              */
        int             nstcgsteep              # number of steps after which a steepest       */
                                                 # descents step is done while doing cg         */
        int             nbfgscorr               # Number of corrections to the hessian to keep */
        int             eConstrAlg              # Type of constraint algorithm                 */
        int             nProjOrder              # Order of the LINCS Projection Algorithm      */
        real            LincsWarnAngle          # If bond rotates more than %g degrees, warn   */
        int             nLincsIter              # Number of iterations in the final Lincs step */
        gmx_bool        bShakeSOR               # Use successive overrelaxation for shake      */
        real            bd_fric                 # Friction coefficient for BD (amu/ps)         */
        gmx_int64_t     ld_seed                 # Random seed for SD and BD                    */
        int             nwall                   # The number of walls                          */
        int             wall_type               # The type of walls                            */
        real            wall_r_linpot           # The potentail is linear for r<=wall_r_linpot */
        int             wall_atomtype[2]        # The atom type for walls                      */
        real            wall_density[2]         # Number density for walls                     */
        real            wall_ewald_zfac         # Scaling factor for the box for Ewald         */

        # COM pulling data */
        # gmx_bool              bPull             # Do we do COM pulling?                        */
        # struct pull_params_t *pull              # The data for center of mass pulling          */
        # TODO: Remove this by converting pull into a ForceProvider
        # struct pull_t        *pull_work         # The COM pull force calculation data structure */

        # AWH bias data */
        gmx_bool                 bDoAwh    # Use awh biasing for PMF calculations?        */
        # gmx::AwhParams          *awhParams # AWH biasing parameters                       */
        # TODO: Remove this by converting AWH into a ForceProvider
        # gmx::Awh                *awh       # AWH work object */

        # Enforced rotation data */
        gmx_bool                 bRot           # Calculate enforced rotation potential(s)?    */
        t_rot                   *rot            # The data for enforced rotation potentials    */

        int                      eSwapCoords    # Do ion/water position exchanges (CompEL)?    */
        t_swapcoords            *swap

        gmx_bool                 bIMD           # Allow interactive MD sessions for this .tpr? */
        t_IMD                   *imd            # Interactive molecular dynamics               */

        real                     cos_accel      # Acceleration for viscosity calculation       */
        tensor                   deform         # Triclinic deformation velocities (nm/ps)     */
        int                      userint1       # User determined parameters                   */
        int                      userint2
        int                      userint3
        int                      userint4
        real                     userreal1
        real                     userreal2
        real                     userreal3
        real                     userreal4
        t_grpopts                opts          # Group options                */
        gmx_bool                 bQMMM         # QM/MM calculation                            */
        int                      QMconstraints # constraints on QM bonds                      */
        int                      QMMMscheme    # Scheme: ONIOM or normal                      */
        real                     scalefactor   # factor for scaling the MM charges in QM calc.*/

        # Fields for removed features go here (better caching) */
        gmx_bool                 bAdress      # Whether AdResS is enabled - always false if a valid .tpr was read
        gmx_bool                 useTwinRange # Whether twin-range scheme is active - always false if a valid .tpr was read

        # gmx::KeyValueTreeObject *params
