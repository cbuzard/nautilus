! -----------------------------------------------------------------------
!   #     #     #     #     #  #######  ###  #        #     #   #####   
!   ##    #    # #    #     #     #      #   #        #     #  #     #  
!   # #   #   #   #   #     #     #      #   #        #     #  #        
!   #  #  #  #     #  #     #     #      #   #        #     #   #####   
!   #   # #  #######  #     #     #      #   #        #     #        #  
!   #    ##  #     #  #     #     #      #   #        #     #  #     #  
!   #     #  #     #   #####      #     ###  #######   #####    #####   
! -----------------------------------------------------------------------
!                            ?
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~|^"~~~~~~~~~~~~~~~~~~~~~~~~~o~~~~~~~~~~~
!        o                   |                  o      __o
!         o                  |                 o     |X__>
!       ___o                 |                __o
!     (X___>--             __|__            |X__>     o
!                         |     \                   __o
!                         |      \                |X__>
!  _______________________|_______\________________
! <                                                \____________   _
!  \                                                            \ (_)
!   \    O       O       O                                       >=)
!    \__________________________________________________________/ (_)
!
!                            ___
!                           / o \
!                      __   \   /   _
!                        \__/ | \__/ \
!                       \___//|\\___/\
!                        ___/ | \___
!                             |     \
!                            /
! -----------------------------------------------------------------------
! 
! A fast 1D gas-grain chemical model by FH (2008)
! Based upon the OSU gas-grain chemical model
! Uses the OSU chemical network
! Updated from gg_osu_2006v1d by RTG/VW
! Rate equations from Hasegawa & Herbst (1992)
! Modified rates following Caselli et al. (1998)
! Stiff solver for sparse Jacobians: LSODES/ODEPACK (Hindmarsh 1983)
! Turbulent mixing implemented through first order operator splitting
! 
! April 2011 VW  An appoximative calculation of the X-ray ionization rate 
! have been added according to Glasgold at al. (1999)
!
! INPUT FILES
! 
! All parameters file can have comments, either a full line or the end of a line. The comment character being the '!' character. 
! Blanck lines are ignored
! 
! parameters.in : parameter file of the code, with various flags
! 
! abundances.in : Give initial abundances for a set of species (a reduced number or all. Default minimum values are applied for those
! that do not exist in this file.
! 
! activation_energies.in : Activation energy for endothermic reactions
! 
! element.in : name and mass in AMU of all elemental species existing in the simulation
! 
! gas_reactions.in : Reaction that occurs in gas phase
! 
! gas_species.in : species that are involved in gas phase reactions
! 
! grain_reactions.in : Reactions that occurs on grain surface
! 
! grain_species.in : Species that are involved in grain surface reactions
! 
! surface_parameters.in : various energies and parameters for diffusion and movements on the grain surface
!
! OUTPUT FILES
! *.out files are output files. *.tmp files are file that are generated at each timestep, either to continue a 
! simulation or check if there is a problem
!
! abundances.*.out : writing in binary format the abundances of all species at each timestep (one file per timestep-output
!
! abundances.tmp : writing in ASCII format the abundances of all species at the current timestep-output
!
! rates.*.out : writing in binary format the critical reactions rates
!
! info.out : writing information on the code, in particular about the code version (commit ID and branch)
!
! species.out : Writing the list of species and their corresponding index
!
! elemental_abundances.tmp/out : writing information about prime elements, their total abundances and mass
!
!
! 
! -----------------------------------------------------------------------

PROGRAM Gasgrain

use global_variables
use iso_fortran_env
use shielding
use utilities
use dust_temperature
use nautilus_main

implicit none

! Parameters for DLSODES. RTOL is the RELATIVE_ABUNDANCE parameter in global_variables.f90
integer :: itol = 2 !< ITOL = 1 or 2 according as ATOL (below) is a scalar or array.
integer :: itask = 1 !< ITASK = 1 for normal computation of output values of Y at t = TOUT.
integer :: istate = 1 !< ISTATE = integer flag (input and output).  Set ISTATE = 1.
integer :: iopt = 1 !< IOPT = 0 to indicate no optional inputs used.
integer :: mf = 21 !< method flag.  Standard values are:
!!\n          10  for nonstiff (Adams) method, no Jacobian used
!!\n          121 for stiff (BDF) method, user-supplied sparse Jacobian
!!\n          222 for stiff method, internally generated sparse Jacobian
real(double_precision) :: atol = 1.d-99 !< absolute tolerance parameter (scalar or array).
!!\n          The estimated local error in Y(i) will be controlled so as
!!\n          to be roughly less (in magnitude) than
!!\n             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!!\n             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!!\n          Thus the local error test passes if, in each component,
!!\n          either the absolute error is less than ATOL (or ATOL(i)),
!!\n          or the relative error is less than RTOL.
!!\n          Use RTOL = 0.0 for pure absolute error control, and
!!\n          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!!\n          control.  Caution: actual (global) errors may exceed these
!!\n          local tolerances, so choose them conservatively.
real(double_precision) :: integration_timestep !< Timestep of the present step, starting from current_time [s]

call initialisation()

call initialize_work_arrays()

! The real time loop
do while (current_time.lt.0.9*STOP_TIME)

  ! Log spacing of time outputs when there is no diffusion
  if (current_time.gt.1.d-2) then
    integration_timestep = current_time*(10.**(1.d0/OUTPUT_PER_DECADE)-1.)
  else
    integration_timestep = 1.d0 * YEAR
  endif
  
  call get_structure_properties(time=current_time, & ! Inputs
                              Av=visual_extinction, density=H_number_density, & ! Outputs
                              gas_temperature=gas_temperature, grain_temperature=dust_temperature) ! Outputs
                                
  current_time = current_time + integration_timestep ! Final time for chemistry T -> T + integration_timestep
  timestep = timestep + 1

  write(Output_Unit,'(A,I5,A,1PD10.3,A)') 'Time=',timestep,', TIME=',current_time/YEAR,' yrs'

  call integrate_chemical_scheme(integration_timestep,itol,atol,itask,istate,iopt,mf)

  ! Output of the rates once every 10 chemical outputs
  call write_current_rates()

  if (istate.eq.-3) stop

  if (mod(timestep,wstep).eq.0) then
    call write_current_output()
  endif

enddo

call write_abundances('abundances.tmp')

contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Chemically evolve for a given time delta_t
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine integrate_chemical_scheme(delta_t,itol,atol,itask,istate,iopt,mf)

  use global_variables
  
  implicit none

  ! Inputs
  real(double_precision), intent(in) :: delta_t !<[in] time during which we must integrate
  integer, intent(in) :: itol !<[in] ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
  integer, intent(in) :: itask !<[in] ITASK  = 1 for normal computation of output values of Y at t = TOUT.
  integer, intent(in) :: iopt !<[in] IOPT   = 0 to indicate no optional inputs used.
  integer, intent(in) :: mf !<[in] method flag.  Standard values are:
!!\n          10  for nonstiff (Adams) method, no Jacobian used
!!\n          121 for stiff (BDF) method, user-supplied sparse Jacobian
!!\n          222 for stiff method, internally generated sparse Jacobian
  real(double_precision), intent(in) :: atol !<[in] integrator tolerance

  ! Outputs
  integer, intent(out) :: istate !<[out] ISTATE = 2  if DLSODES was successful, negative otherwise.
!!\n          -1 means excess work done on this call (perhaps wrong MF).
!!\n          -2 means excess accuracy requested (tolerances too small).
!!\n          -3 means illegal input detected (see printed message).
!!\n          -4 means repeated error test failures (check all inputs).
!!\n          -5 means repeated convergence failures (perhaps bad Jacobian
!!\n             supplied or wrong choice of MF or tolerances).
!!\n          -6 means error weight became zero during problem. (Solution
!!\n             component i vanished, and ATOL or ATOL(i) = 0.)
!!\n          -7 means a fatal error return flag came from sparse solver
!!\n             CDRV by way of DPRJS or DSOLSS.  Should never happen.
!!\n          A return with ISTATE = -1, -4, or -5 may result from using
!!\n          an inappropriate sparsity structure, one that is quite
!!\n          different from the initial structure.  Consider calling
!!\n          DLSODES again with ISTATE = 3 to force the structure to be
!!\n          reevaluated. 
 
  ! Locals
  integer :: i
  real(double_precision) :: t !< The local time, starting from 0 to delta_t
  real(double_precision), dimension(nb_species) :: satol !< Array that contain the absolute tolerance, 
  !! one value per species, either a minimum value or a value related to its abundance.
  real(double_precision), dimension(nb_species) :: temp_abundances !< Temporary array that contain the 
  !! abundances for the duration of the timestep, a sort of buffer.

  real(double_precision) :: t_stop_step

  t_stop_step = delta_t
  t = 0.d0

  temp_abundances(1:nb_species) = abundances(1:nb_species)

  do while (t.lt.t_stop_step)

    istate = 1

    ! Adaptive absolute tolerance to avoid too high precision on very abundant species,
    ! H2 for instance. Helps running a bit faster

    do i=1,nb_species
      satol(i) = max(atol, 1.d-16 * temp_abundances(i))
    enddo

    ! Feed IWORK with IA and JA

    call set_work_arrays(Y=temp_abundances)
    

    call dlsodes(get_temporal_derivatives,nb_species,temp_abundances,t,t_stop_step,itol,RELATIVE_TOLERANCE,&
    satol,itask,istate,iopt,rwork,lrw,iwork,liw,get_jacobian,mf)       

    ! Whenever the solver fails converging, print the reason.
    ! cf odpkdmain.f for translation
    if (istate.ne.2) then
      write(*,*)  'ISTATE = ', ISTATE
    endif

    call check_conservation(temp_abundances)

  enddo

  abundances(1:nb_species) = temp_abundances(1:nb_species)

  return 
  end subroutine integrate_chemical_scheme

! ======================================================================
! ======================================================================
END program gasgrain
