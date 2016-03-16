
from utility cimport *
from math cimport *

cdef extern from "gromacs/legacyheaders/types/inputrec.h":
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
        int             eI                       # Integration method                 */
        int             nsteps                   # number of steps to be taken			*/
        int             simulation_part          # Used in checkpointing to separate chunks */
        int             init_step                # start at a stepcount >0 (used w. convert-tpr)    */
        int             nstcalcenergy            # frequency of energy calc. and T/P coupl. upd.	*/
        int             cutoff_scheme            # group or verlet cutoffs     */
        int             ns_type                  # which ns method should we use?               */
        int             nstlist                  # number of steps before pairlist is generated	*/
        int             ndelta                   # number of cells per rlong			*/
        int             nstcomm                  # number of steps after which center of mass	*/
                                               # motion is removed				*/
        int             comm_mode                # Center of mass motion removal algorithm      */
        int             nstlog                   # number of steps after which print to logfile	*/
        int             nstxout                  # number of steps after which X is output	*/
        int             nstvout                  # id. for V					*/
        int             nstfout                  # id. for F					*/
        int             nstenergy                # number of steps after which energies printed */
        int             nstxout_compressed       # id. for compressed trj (.xtc,.tng)           */
        double          init_t                   # initial time (ps)              */
        double          delta_t                  # time step (ps)				*/
        double          x_compression_precision  # precision of x in compressed trajectory file */
        double          fourier_spacing          # requested fourier_spacing, when nk? not set  */
        int             nkx, nky, nkz            # number of k vectors in each spatial dimension*/
                                               # for fourier methods for long range electrost.*/
        int             pme_order                # interpolation order for PME                  */
        double          ewald_rtol               # double space tolerance for Ewald, determines   */
                                               # the double/reciprocal space relative weight    */
        double          ewald_rtol_lj            # double space tolerance for LJ-Ewald            */
        int             ewald_geometry           # normal/3d ewald, or pseudo-2d LR corrections */
        double          epsilon_surface          # Epsilon for PME dipole correction            */
        int             ljpme_combination_rule   # Type of combination rule in LJ-PME          */
        int             ePBC                     # Type of periodic boundary conditions		*/
        int             bPeriodicMols            # Periodic molecules                           */
        bint            bContinuation            # Continuation run: starting state is correct	*/
        int             etc                      # temperature coupling               */
        int             nsttcouple               # interval in steps for temperature coupling   */
        bint            bPrintNHChains           # whether to print nose-hoover chains        */
        int             epc                      # pressure coupling                            */
        int             epct                     # pressure coupling type			*/
        int             nstpcouple               # interval in steps for pressure coupling      */
        double          tau_p                    # pressure coupling time (ps)			*/
        tensor          ref_p                    # reference pressure (kJ/(mol nm^3))		*/
        tensor          compress                 # compressability ((mol nm^3)/kJ)        */
        int             refcoord_scaling         # How to scale absolute reference coordinates  */
        rvec            posres_com               # The COM of the posres atoms                  */
        rvec            posres_comB              # The B-state COM of the posres atoms          */
        int             andersen_seed            # Random seed for Andersen thermostat (obsolete) */
        double          verletbuf_tol            # Per atom pair energy drift tolerance (kJ/mol/ps/atom) for list buffer  */
        double          rlist                    # short range pairlist cut-off (nm)		*/
        double          rtpi                     # Radius for test particle insertion           */
        int             coulombtype              # Type of electrostatics treatment             */
        int             coulomb_modifier         # Modify the Coulomb interaction              */
        double          rcoulomb_switch          # Coulomb switch range start (nm)		*/
        double          rcoulomb                 # Coulomb cutoff (nm)		                */
        double          epsilon_r                # relative dielectric constant                 */
        double          epsilon_rf               # relative dielectric constant of the RF       */
        int             implicit_solvent         # No (=explicit water), or GBSA solvent models */
        int             gb_algorithm             # Algorithm to use for calculation Born radii  */
        int             nstgbradii               # Frequency of updating Generalized Born radii */
        double          rgbradii                 # Cutoff for GB radii calculation              */
        double          gb_saltconc              # Salt concentration (M) for GBSA models       */
        double          gb_epsilon_solvent       # dielectric coeff. of implicit solvent     */
        double          gb_obc_alpha             # 1st scaling factor for Bashford-Case GB      */
        double          gb_obc_beta              # 2nd scaling factor for Bashford-Case GB      */
        double          gb_obc_gamma             # 3rd scaling factor for Bashford-Case GB      */
        double          gb_dielectric_offset     # Dielectric offset for Still/HCT/OBC     */
        int             sa_algorithm             # Algorithm for SA part of GBSA                */
        double          sa_surface_tension       # Energy factor for SA part of GBSA */
        int             vdwtype                  # Type of Van der Waals treatment              */
        int             vdw_modifier             # Modify the VdW interaction                   */
        double          rvdw_switch              # Van der Waals switch range start (nm)        */
        double          rvdw                     # Van der Waals cutoff (nm)	        */
        int             eDispCorr                # Perform Long range dispersion corrections    */
        double          tabext                   # Extension of the table beyond the cut-off,   *
                                               #  * as well as the table length for 1-4 interac. */
        double          shake_tol                # tolerance for shake				*/
        int             efep                     # free energy calculations                     */
        t_lambda       *fepvals                  # Data for the FEP state                       */
        bint            bSimTemp                 # Whether to do simulated tempering            */
        t_simtemp      *simtempvals              # Variables for simulated tempering            */
        bint            bExpanded                # Whether expanded ensembles are used          */
        t_expanded     *expandedvals             # Expanded ensemble parameters              */
        int             eDisre                   # Type of distance restraining                 */
        double          dr_fc                    # force constant for ta_disre			*/
        int             eDisreWeighting          # type of weighting of pairs in one restraints	*/
        bint            bDisreMixed              # Use comb of time averaged and instan. viol's	*/
        int             nstdisreout              # frequency of writing pair distances to enx   */
        double          dr_tau                   # time constant for memory function in disres    */
        double          orires_fc                # force constant for orientational restraints  */
        double          orires_tau               # time constant for memory function in orires    */
        int             nstorireout              # frequency of writing tr(SD) to enx           */
        double          em_stepsize              # The stepsize for updating			*/
        double          em_tol                   # The tolerance				*/
        int             niter                    # Number of iterations for convergence of      */
                                               # steepest descent in relax_shells             */
        double          fc_stepsize              # Stepsize for directional minimization        */
                                               # in relax_shells                              */
        int             nstcgsteep               # number of steps after which a steepest       */
                                               # descents step is done while doing cg         */
        int             nbfgscorr                # Number of corrections to the hessian to keep */
        int             eConstrAlg               # Type of constraint algorithm                 */
        int             nProjOrder               # Order of the LINCS Projection Algorithm      */
        double          LincsWarnAngle           # If bond rotates more than %g degrees, warn   */
        int             nLincsIter               # Number of iterations in the final Lincs step */
        bint            bShakeSOR                # Use successive overrelaxation for shake      */
        double          bd_fric                  # Friction coefficient for BD (amu/ps)         */
        int             ld_seed                  # Random seed for SD and BD                    */
        int             nwall                    # The number of walls                          */
        int             wall_type                # The type of walls                            */
        double          wall_r_linpot            # The potentail is linear for r<=wall_r_linpot */
        int             wall_atomtype[2]         # The atom type for walls                      */
        double          wall_density[2]          # Number density for walls                     */
        double          wall_ewald_zfac          # Scaling factor for the box for Ewald         */

        # COM pulling data */
        bint                  bPull              # Do we do COM pulling?                        */
        #struct pull_params_t *pull               # The data for center of mass pulling          */
        #struct pull_t        *pull_work          # The COM pull force calculation data structure  TODO this pointer should live somewhere else */

        # Enforced rotation data */
        bint            bRot                     # Calculate enforced rotation potential(s)?    */
        t_rot          *rot                      # The data for enforced rotation potentials    */

        int             eSwapCoords              # Do ion/water position exchanges (CompEL)?    */
        t_swapcoords   *swap

        bint            bIMD                     # Allow interactive MD sessions for this .tpr? */
        t_IMD          *imd                      # Interactive molecular dynamics               */

        double          cos_accel                # Acceleration for viscosity calculation       */
        tensor          deform                   # Triclinic deformation velocities (nm/ps)     */
        int             userint1                 # User determined parameters                   */
        int             userint2
        int             userint3
        int             userint4
        double          userdouble1
        double          userdouble2
        double          userdouble3
        double          userdouble4
        t_grpopts       opts           # Group options				*/
        t_cosines       ex[0]        # Electric field stuff	(spatial part)		*/
        t_cosines       et[0]        # Electric field stuff	(time part)		*/
        bint            bQMMM          # QM/MM calculation                            */
        int             QMconstraints  # constraints on QM bonds                      */
        int             QMMMscheme     # Scheme: ONIOM or normal                      */
        double          scalefactor    # factor for scaling the MM charges in QM calc.*/

        # Fields for removed features go here (better caching) */
        bint            bAdress        # Whether AdResS is enabled - always false if a valid .tpr was read
        bint            useTwinRange   # Whether twin-range scheme is active - always false if a valid .tpr was read
