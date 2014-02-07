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
use diffusion
use input_output
use model_1D
use ode_solver
use iso_fortran_env
use shielding

implicit none

integer :: lrw, liw

real(double_precision), dimension(:), allocatable :: temp_abundances ! nb_species

integer :: itol = 2
integer :: itask = 1
integer :: istate = 1
integer :: iopt = 1
integer :: mf = 21
real(double_precision) :: atol = 1.d-99
real(double_precision) :: T, TOUT, TIN

integer :: ipts

call initialisation()

allocate(temp_abundances(nb_species))

! Initializing T
! T = local, TIME = global
T = 0.d0
TIME = 0.d0

lrw = 20 + 3 * nb_nonzeros_values*nb_species + 21 * nb_species
liw = 31 + 3 * nb_nonzeros_values*nb_species + 21 * nb_species

! Allocate the JA array
allocate(JA(liw))

! The real time loop
do while (t.lt.0.9*STOP_TIME)

  call ztimestep() ! Determination of the diffusive timestep, modification of the global parameter DIFFUSIVE_TIMESTEP
  TIME = T + DIFFUSIVE_TIMESTEP ! Final time for chemistry T -> T + DIFFUSIVE_TIMESTEP

  timestep = timestep + 1

  ! Store the current time in TIN (for 1D calculation)
  TIN = T

  write(Output_Unit,'(A,I5,A,1PD10.3,A)') 'Time=',timestep,', TIME=',TIME/TYEAR,' yrs'

  do ipts=1,nptmax ! Start of the spatial loop for chemistry

    iptstore = ipts

    ! T being changed in dlsode, needs to be defined again

    T=TIN

    ! Feed 1D physical structure
    TAU=TAU1D(ipts)
    X_IONISATION_RATE=X_IONISATION_RATE1D(ipts)
    ! XNT is the total density of H atoms (both H and H2)
    ! But DENS1D is the real number density
    ! Problem with the sound speed if NH > NH2 (cf phys_1D)
    XNT=2.*DENS1D(ipts)
    TEMP=TEMP1D(ipts)
    DTEMP=DTEMP1D(ipts)

    ! Chemical evolution for each spatial point

    temp_abundances(:nb_species) = ZXN(:,ipts)
    abundances(:nb_species) = ZXN(:,ipts)

    call integrate_chemical_scheme(T,temp_abundances,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

    ! Output of the rates once every 10 chemical outputs
    !      if ((mod(it,wstepr).eq.0).and.(ipts.eq.irateout)) then
    call write_current_rates()
    !      endif

    if (istate.eq.-3) stop

    ZXN(:,ipts) = abundances(:) ! putting back

  enddo ! end of the spatial loop for chemistry 

  ! Generic call to zdiffusion, a subroutine calling the chosen numerical scheme 
  call zdiffusion() ! diffusion 

  if (mod(timestep,wstep).eq.0) then
    call write_current_output()
  endif

enddo

if (nptmax.eq.1) call write_abundances('abundances.tmp')

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine that contain all initialisation that needs to be done in the code
!! before the integration
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine initialisation()
use global_variables

implicit none

! Locals
integer :: i ! For loops

! Read list of prime elements, including their atomic mass (in AMU)
call read_element_in()

! Get various size needed for allocatable arrays
call get_array_sizes()

! Init global allocatable arrays. From now on, we can read data files
call initialize_global_arrays()

! Read simulation parameters
call read_parameters_in()

! Read list of species, either for gas or grain reactions
call read_species()

! Read list of reactions for gas and grains
call read_reactions()

! Overwrite the parameter file to update the syntax, organization and add default values for parameters that did not exist or weren't defined
call write_parameters()

! Set initial abundances. Will define a minimum value for species that are not present in the input file
call read_abundances()

! From the total list of species, determine the exact number of species that are in gas phase, and on the surface of grains. 
! This is different from the list in input files were the lists correspond to species NEEDED for reactions in gas phase or on grains. 
! Gas species can be needed for surface reactions, and vice versa.
call get_gas_surface_species()

! Build spatial mesh 
call mesh()

! Initialization of elemental/chemical quantities
call index_datas()

timestep=0

! Calculate the initial abundances for all elements that compose 
call get_elemental_abundance(all_abundances=abundances, el_abundances=INITIAL_ELEMENTAL_ABUNDANCE)

! Store initial elemental abundances
call write_elemental_abundances(filename='elemental_abundances.out', el_abundances=INITIAL_ELEMENTAL_ABUNDANCE)

! Recompute initial_dtg_mass_ratio to remove He
! In the following, initial_dtg_mass_ratio is used as a H/dust mass ratio
do i=1,NB_PRIME_ELEMENTS
  if (species_name(PRIME_ELEMENT_IDX(I)).EQ.YHE) then
    initial_dtg_mass_ratio = initial_dtg_mass_ratio*(1.d0+4*INITIAL_ELEMENTAL_ABUNDANCE(I))
    ! Mean molecular weight (cgs) 
    ! Approximated here (the exact calculus would require a sume over AWT
    ! Used for the diffusion in disks, not for chemistry
    mean_molecular_weight = 2.d0 + 4.d0*INITIAL_ELEMENTAL_ABUNDANCE(I)/(1.d0+INITIAL_ELEMENTAL_ABUNDANCE(I))
  endif
enddo

! Compute the grain abundance
GTODN=(4.D+0*PI*GRAIN_DENSITY*grain_radius*grain_radius*grain_radius)/(3.D+0*initial_dtg_mass_ratio*AMU)

where(species_name.EQ.YGRAIN) abundances=1.0/GTODN

! Set the electron abundance via conservation===========
! And check at the same time that nls_init has the same elemental abundance
! as nls_control
call check_conservation(abundances(1:nb_species))

! 1D physical structure (nls_phys_1D)
call phys_1D()

! Write species name/index correspondance
call write_species()

! Initialize indices of reactants and products 
call set_chemical_reactants()

! Initialize the arrays that list, for each species, the reactions using it as a reactant
!! max_reactions_same_species is set here. nb_reactions_using_species and relevant_reactions array are set here.
call init_relevant_reactions()

! Initializing ZXN
do ipts=1,nptmax
  ZXN(:,ipts) = abundances(:)
enddo

! Calculate the optimum number for temporary solving-arrays in ODEPACK, based on the number of non-zeros values in 
!! the jacobian
call count_nonzeros()

! Write information about the code at the very end of initialisation
call write_general_infos()

end subroutine initialisation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Read and retrieve information for species. 
!! Save index for prime elements in the full species array
!!
!! Save other interesting index for special species
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine index_datas()
use global_variables

implicit none

! Locals
real(double_precision) :: MSUM
integer :: ILAB, j, k, i, isptemp
integer :: KSUM ! sum of number of primary element composing the species. If equal to 1, the current species is elemental
real(double_precision) :: mass_tmp !< temporary value to exchange two index in the mass array
character(len=11) :: name_tmp !< temporary value to exchange two index in the name array

! Set elements' characteristics=========================================
! --- Find the atomic species associated with a given element

ILAB=1
do J=1,nb_species
  KSUM=0
  ! ------ Calculate species' elemental_mass
  do K=1,NB_PRIME_ELEMENTS
    KSUM=KSUM+species_composition(K,J)
  enddo
  ! ------ Check for atomic species
  if ((KSUM.EQ.1).AND.(ICG(J).EQ.0).AND.&
  (species_name(J)(:1).NE.'J          ').AND.(species_name(J)(:1).NE.'X          ')) then
    if (ILAB.GT.NB_PRIME_ELEMENTS) then
      write(Error_unit, *) '***More fundamental elements than NB_PRIME_ELEMENTS***'
      call exit(3)
    endif       
    ! --------- Save species number
    PRIME_ELEMENT_IDX(ILAB)=J
    ILAB=ILAB+1
  endif

  ! ------ Check for electron species number
  IF (species_name(J).EQ.'e-         ') then
    ISPE=J
  endif
enddo



! --- Re-arrange order of elements to match species_composition columns (reactions file)
do J=1,NB_PRIME_ELEMENTS-1
  if (species_composition(J,PRIME_ELEMENT_IDX(J)).NE.1) then
    do K=J+1,NB_PRIME_ELEMENTS
      if (species_composition(J,PRIME_ELEMENT_IDX(K)).EQ.1) then
        ISPTEMP=PRIME_ELEMENT_IDX(K)
        mass_tmp = elemental_mass(k)
        name_tmp = element_name(k)
        
        elemental_mass(k) = elemental_mass(j)
        element_name(k)  = element_name(j)
        PRIME_ELEMENT_IDX(K)=PRIME_ELEMENT_IDX(J)
        
        elemental_mass(k) = mass_tmp
        element_name(k) = name_tmp
        PRIME_ELEMENT_IDX(J)=ISPTEMP
      endif
    enddo
  endif
enddo

! Set species' characteristics==========================================
! --- Set special species labels
YH     = 'H          '
YJH    = 'JH         '
YH2    = 'H2         '
YJH2   = 'JH2        '
YHE    = 'He         '
YHEP   = 'He+        '
YE     = 'e-         '
YGRAIN = 'GRAIN0     '
YCO    = 'CO         '

! --- Set reference species
do I=1,nb_species 
  ! ------ Calculate elemental_masses
  MSUM=0.d0
  do K=1,NB_PRIME_ELEMENTS 
    MSUM=MSUM+elemental_mass(K)*species_composition(K,I) 
  enddo 
  AWT(I)=MSUM
  if (species_name(I).EQ.YE) AWT(I)=1.D+0/1836.D+0 
  if (species_name(I).EQ.YGRAIN .OR. species_name(I).EQ.'GRAIN-      ')&
  AWT(I)=4.0*PI*grain_radius*grain_radius*grain_radius*GRAIN_DENSITY/3.0/AMU
enddo

! Initialize the Av/NH ratio
! Can be scaled for different dust/gas ratios
! Here initial_dtg_elemental_mass_ratio is the original dust/gas elemental_mass ratio (from nls_control.d)
! initial_dtg_elemental_mass_ratio is changed later into the dust/hydrogen elemental_mass ratio

AV_NH_ratio = 5.34d-22 * (initial_dtg_mass_ratio / 1.d-2)

! Find ITYPE first and last reactions===================================
do I=0,NITYPE-1
  IRXSTA(I)=0
  IRXFIN(I)=0
  do J=1,nb_reactions
    if ((ITYPE(J).EQ.I).AND.(IRXSTA(I).EQ.0)) IRXSTA(I)=J
    if (ITYPE(J).EQ.I) IRXFIN(I)=J
  enddo
enddo

! Find the index of CO and H2
do i=1,nb_species
  if (species_name(i).eq.YH2) INDH2=i
  if (species_name(i).eq.YCO) INDCO=i
  if (species_name(i).eq.YHE) INDHE=i
enddo

! Compute nb_sites_per_grain = number of sites per grain
nb_sites_per_grain = SITE_DENSITY*4.d0*PI*grain_radius**2

! Initialise reaction rates=============================================
call init_reaction_rates()

return 
end subroutine index_datas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Chemically evolve from T to TOUT the given spatial point
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine integrate_chemical_scheme(T,temp_abundances,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

  use global_variables
  
  implicit none
  
  ! Inputs
  integer, intent(in) :: liw
  integer, intent(in) :: lrw
  integer, intent(in) :: itol
  integer, intent(in) :: itask
  integer, intent(in) :: iopt
  integer, intent(in) :: mf
  real(double_precision), intent(in) :: atol

  ! Outputs
  integer, intent(out) :: istate
  real(double_precision), intent(out) :: t
  real(double_precision), intent(out) :: TOUT
  real(double_precision), intent(out), dimension(nb_species) :: temp_abundances
 
  ! Locals
  integer, dimension(liw) :: IWORK
  real(double_precision), dimension(lrw) :: RWORK
  integer :: i
  real(double_precision), dimension(nb_species) :: satol

  real(double_precision) :: TIN

  integer :: NNZ

  ! Initialize work arrays

  iwork(:) = 0
  rwork(:) = 0.d0
  IWORK(5) = 5
  RWORK(6) = 3.154D14
  IWORK(6) = 10000
  IWORK(7) = 2

  if (timestep.eq.1) then
    IWORK(6)=2000
  endif
  
  ! Changing the time to avoid T + DT = T 

  TOUT=TIME

  TIN = T
  TOUT = TOUT - TIN
  T = 0.d0      

  temp_abundances(:) = abundances(:)

  call shieldingsetup()

  ! Don't stop running, Forrest

  do while (t.lt.tout)

    istate = 1

    ! Adaptive absolute tolerance to avoid too high precision on very abundant species,
    ! H2 for instance. Helps running a bit faster

    do i=1,nb_species
      satol(i)=max(atol,1.d-16*temp_abundances(i))
    enddo

    ! set_constant_ratesd is already called in compute IAJA
    !      call set_constant_rates(Y)

    ! Feed IWORK with IA and JA

    call computeIAJA(temp_abundances)
    NNZ=IA(nb_species+1)-1
    iwork(30+1:30+nb_species+1)=IA(1:nb_species+1)
    iwork(31+nb_species+1:31+nb_species+NNZ)=JA(1:NNZ)

    call dlsodes(get_temporal_derivatives,nb_species,temp_abundances,t,tout,itol,RELATIVE_TOLERANCE,&
    satol,itask,istate,iopt,rwork,lrw,iwork,liw,get_jacobian,mf)       

    ! Whenever the solver fails converging, print the reason.
    ! cf odpkdmain.f for translation
    if (istate.ne.2) write(*,*)  'IPTS = ', ipts, 'ISTATE = ', ISTATE

    call check_conservation(temp_abundances)

    ! Stop, Forrest
  enddo

  T = TOUT + TIN 

  abundances(:) = temp_abundances(:)

  return 
  end subroutine integrate_chemical_scheme

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Check if elementary abundances are conserved. If not, display a warning. 
!! The option CONSERVATION_TYPE can override abundances and skip the warning display.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine check_conservation(temp_abundances)
  use global_variables
  
  implicit none
  
  ! Inputs/outputs
  real(double_precision), intent(inout), dimension(nb_species) :: temp_abundances !< [in,out] The list of abundances for all species. 
!! Can be overwritten if CONSERVATION_TYPE is on, and force conservation.
  
  ! Locals
  real(double_precision), dimension(NB_PRIME_ELEMENTS) :: elemental_abundance
  real(double_precision) :: CHASUM

  integer :: i, k

  ! Prevent too low abundances
  do i=1,nb_species
    if (temp_abundances(i).le.1.d-99) then
      temp_abundances(i) = 1.d-99
    endif
  enddo

  ! --- Conserve electrons
  CHASUM=0.d0
  do I=1,nb_species
    if (I.NE.ISPE) CHASUM=CHASUM+ICG(I)*temp_abundances(I)
  enddo
  if (CHASUM.LE.0.d0) CHASUM=MINIMUM_INITIAL_ABUNDANCE
  temp_abundances(ISPE)=CHASUM

  ! --- Conserve other elements if selected
  if (CONSERVATION_TYPE.GT.0) then
    do K=1,CONSERVATION_TYPE
      elemental_abundance(K)=0.d0
    enddo
    do I=1,nb_species
      do K=1,CONSERVATION_TYPE
        if (I.NE.PRIME_ELEMENT_IDX(K)) elemental_abundance(K)=elemental_abundance(K)+species_composition(K,I)*temp_abundances(I)
      enddo
    enddo
    do K=1,CONSERVATION_TYPE
      temp_abundances(PRIME_ELEMENT_IDX(K))=INITIAL_ELEMENTAL_ABUNDANCE(K)-elemental_abundance(K)
      if (temp_abundances(PRIME_ELEMENT_IDX(K)).LE.0.d0) temp_abundances(PRIME_ELEMENT_IDX(K))=MINIMUM_INITIAL_ABUNDANCE
    enddo
  endif

  ! Check for conservation
  call get_elemental_abundance(all_abundances=temp_abundances, el_abundances=elemental_abundance)
  
  call write_elemental_abundances(filename='elemental_abundances.tmp', el_abundances=elemental_abundance)

  do k=1,NB_PRIME_ELEMENTS
    if (abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-elemental_abundance(K))/INITIAL_ELEMENTAL_ABUNDANCE(K).ge.0.01d0) then 
      write(Error_unit,*) 'Caution: Element ', trim(species_name(PRIME_ELEMENT_IDX(K))), 'is not conserved'
      write(Error_unit,*) 'Relative difference: ', abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-elemental_abundance(K)) / &
                           INITIAL_ELEMENTAL_ABUNDANCE(K)
    endif
    if (species_name(PRIME_ELEMENT_IDX(K)).eq.YH) then
      if (abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-temp_abundances(INDH2)*2.D0)/INITIAL_ELEMENTAL_ABUNDANCE(K).ge.0.01d0) then
        write(Error_unit,*) 'H is too depleted on the grains !!!!'
      endif
    endif
    if (species_name(PRIME_ELEMENT_IDX(K)).eq.YHE) then
      if (abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-temp_abundances(INDHE))/INITIAL_ELEMENTAL_ABUNDANCE(K).ge.0.01d0) then
        write(Error_unit,*) 'He is too depleted on the grains !!!!'
      endif
    endif       
  enddo

  ! VW fev 2012 add a test for the helium and H2 abundance in the gas phase
  ! prevent excessive depletion


  return
  end subroutine check_conservation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine that calculate the elemantal abundances from
!! all abundances of all species
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_elemental_abundance(all_abundances, el_abundances)

implicit none

! Inputs
real(double_precision), intent(in), dimension(nb_species) :: all_abundances !< [in] List of abundances for all existing species

! Outputs
real(double_precision), intent(out), dimension(NB_PRIME_ELEMENTS) :: el_abundances !< [out] list of abundances for all fundamental elements

! Locals
integer :: i,j

el_abundances(1:NB_PRIME_ELEMENTS) = 0.0d0

do i=1,nb_species
  do j=1,NB_PRIME_ELEMENTS
    el_abundances(j) = el_abundances(j) + species_composition(j,i) * all_abundances(i)
  enddo
enddo

end subroutine get_elemental_abundance

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Computes the H2 and CO column density for their self-shielding
!! for each ipts, the column density is incremented recursively
!! from the results for NH2 and NCO computed by set_dependant_rates the spatial step before
!
!> @note ZN is shifted with respect to N
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine SHIELDINGSETUP()

  use global_variables
  implicit none
  
  ! Locals
  real(double_precision) :: XNH2,XNCO

  if (iptstore.eq.1) then
    ZNH2(iptstore)=0.d0
    ZNCO(iptstore)=0.d0
  else

    XNH2=ZXN(indH2,iptstore-1)
    XNCO=ZXN(indCO,iptstore-1)

    ZNH2(iptstore)=ZNH2(iptstore-1)+XNT*zstepsize*XNH2
    ZNCO(iptstore)=ZNCO(iptstore-1)+XNT*zstepsize*XNCO
  endif

  return
end subroutine SHIELDINGSETUP


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief compute surface info (thermodynamic, quantum and kinetic data) 
!! from datafiles
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine init_reaction_rates()
    use global_variables
    
    implicit none

    ! Locals
    real(double_precision), dimension(nb_species) :: REA1,REA2,REA3,REA4
    real(double_precision), dimension(nb_reactions) :: REA5
    real(double_precision), dimension(nb_species) :: SMASS
    real(double_precision) :: SMA,REDMAS,STICK,EVFRAC,DHFSUM,SUM1,SUM2
    integer, dimension(nb_reactions) :: INT1
    integer :: NGS,NEA,NPATH,NEVAP,BADFLAG,ATOMS
    character(len=11), dimension(5,nb_reactions) :: GSREAD
    character(len=11), dimension(nb_species) :: GSPEC

    real(double_precision) :: cond
    integer :: i, j,k,l,n4, n5, n6
    
    character(len=80) :: filename !< name of the file to be read
    character(len=80) :: line
    character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
    integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
    integer :: error !< to store the state of a read instruction

    logical :: isDefined

    ! Set accretion rate====================================================

    ! COND is used to calculate R_acc = (sigma_d) * <v_i> * n_i * n_d
    ! Uses 'Mean' Maxwellian speed, rather than RMS (is this justifiable?)

    COND=PI*grain_radius*grain_radius*SQRT(8.0d0*K_B/PI/AMU)

    ! --- Evaluate sticking coeff and accretion rate factor for each species
    STICK=0.d0
    do I=1,nb_species
      if (ICG(I).EQ.0) then
        STICK=sticking_coeff_neutral
      endif 
      if (ICG(I).GT.0) then 
        STICK=sticking_coeff_positive 
      endif
      if (ICG(I).LT.0) then 
        STICK=sticking_coeff_negative 
      endif
      !         if (species_name(I).EQ.YH2)      STICK=0.D+0 
      !         if (species_name(I).EQ.YHE)      STICK=0.D+0 
      !         if (species_name(I).EQ.YH)       STICK=0.D+0 
      if (species_name(I).EQ.YHEP)     STICK=0.D+0 
      if (species_name(I).EQ.'e-         ')      STICK=0.D+0
      if (species_name(I).EQ.'H+         ')     STICK=0.D+0
      if (species_name(I).EQ.YGRAIN)   STICK=0.D+0
      if (species_name(I).EQ.'GRAIN-     ') STICK=0.D+0
      !         if (species_name(I).EQ.'H-')     STICK=0.D+0
      !         if (species_name(I).EQ.'H2+')    STICK=0.D+0

      if (I.GT.nb_gaseous_species) STICK=0.D+0
      CONDSP(I)=COND*STICK/SQRT(AWT(I))
    enddo


    ! Read in molecular information for surface rates=======================
    
    
    ! Reading list of species for gas phase
    filename = 'surface_parameters.in'
    inquire(file=filename, exist=isDefined)
    if (isDefined) then
      call get_linenumber(filename=filename, nb_lines=NGS)

      open(10, file=filename, status='old')
      
      i = 0
      do
        read(10, '(a)', iostat=error) line
        if (error /= 0) exit
          
        ! We get only what is on the left of an eventual comment parameter
          comment_position = index(line, comment_character)
        
        ! if there are comments on the current line, we get rid of them
        if (comment_position.ne.0) then
          line = line(1:comment_position - 1)
        end if
        
        if (line.ne.'') then
          i = i + 1
          read(line, '(A11,I4,F7.0,F6.0,D8.1,27X,F8.2)') GSPEC(I),INT1(I),REA1(I),REA2(I),REA3(I),REA4(I)
        
        end if
      end do
      close(10)
      
    else
      write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
      call exit(1)
    end if
    
    filename = 'activation_energies.in'
    inquire(file=filename, exist=isDefined)
    if (isDefined) then
      call get_linenumber(filename=filename, nb_lines=NEA)

      open(10, file=filename, status='old')
      
      i = 0
      do
        read(10, '(a)', iostat=error) line
        if (error /= 0) exit
          
        ! We get only what is on the left of an eventual comment parameter
          comment_position = index(line, comment_character)
        
        ! if there are comments on the current line, we get rid of them
        if (comment_position.ne.0) then
          line = line(1:comment_position - 1)
        end if
        
        if (line.ne.'') then
          i = i + 1
          read(line, '(5A11,D9.2)') (GSread(L,i),L=1,5),REA5(i)
        
        end if
      end do
      close(10)
      
    else
      write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
      call exit(1)
    end if

    ! --- Transfer from dummies to arrays with correct species numbers
    do I=1,nb_species
      SMASS(I)=0.d0
      ED(I)=0.d0
      EB(I)=0.d0
      DEB(I)=0.d0
      DHF(I)=0.d0
      do J=1,NGS
        if (species_name(I).EQ.GSPEC(J)) then
          SMASS(I)=dble(INT1(J))
          ED(I)=REA1(J)
          EB(I)=REA2(J)
          DEB(I)=REA3(J)
          DHF(I)=REA4(J)
          if ((species_name(I).NE.YJH).AND.(species_name(I).NE.YJH2).AND.&
          (EBFAC.GE.0.d0)) EB(I)=EBFAC*ED(I)
        endif
      enddo
      !IF(species_name(I) == 'JN2O2      ') write(*,*) ED(I)
    enddo

    do I=1,nb_reactions
      EA(I)=0.d0
      do J=1,NEA
        if (SYMBOL(4,I)(:1).EQ.'J') then
          if ((SYMBOL(1,I).EQ.GSread(1,J)).AND.&
          (SYMBOL(2,I).EQ.GSread(2,J)).AND.&
          (SYMBOL(4,I).EQ.GSread(3,J)).AND.&
          (SYMBOL(5,I).EQ.GSread(4,J)).AND.&
          (SYMBOL(6,I).EQ.GSread(5,J))) EA(I)=REA5(J)
        else
          if ((SYMBOL(1,I).EQ.GSread(1,J)).AND.&
          (SYMBOL(2,I).EQ.GSread(2,J)).AND.&
          (SYMBOL(4,I).EQ.GSread(3,J)(2:)).AND.&
          (SYMBOL(5,I).EQ.GSread(4,J)(2:)).AND.&
          (SYMBOL(6,I).EQ.GSread(5,J)(2:))) EA(I)=REA5(J)
        endif
      enddo
      !IF(symbol(4,i) == 'JO2H       ') write(*,*)  symbol(:,i), Ea(i)
    enddo

    ! Set up constants, quantum rate info===================================
    do I=1,nb_species
      CHF(I)=0.d0
      RQ1(I)=0.d0
      RQ2(I)=0.d0
      ! ------ For species which have been assigned surface info, SMASS=/=0
      if (SMASS(I).NE.0) then
        SMA=dble(SMASS(I))
        ! --------- Set characteristic frequency
        CHF(I)=SQRT(2.0d0*K_B/PI/PI/AMU * SITE_DENSITY*ED(I)/SMA)
        ! --------- Set quantum rates
        if (DEB(I).GE.1.0D-38) then
          RQ1(I)=DEB(I)*K_B/4.0d0/H_BARRE/nb_sites_per_grain
        else
          RQ1(I)=0.d0
        endif
        RQ2(I)=CHF(I)/nb_sites_per_grain*&
        EXP(-2.0d0*SITE_SPACING/H_BARRE*SQRT(2.0d0*AMU*SMA*K_B*EB(I)))
      endif
    enddo

    ! === Cycle all reactions
    do J=1,nb_reactions

      ! ------ Initialise all XJ rate factors, and get species 1 & 2
      XJ(J)=1.0d0
      JSP1(J)=0
      JSP2(J)=0
      do I=1,nb_species
        if (SYMBOL(1,J).EQ.species_name(I)) JSP1(J)=I
        if (SYMBOL(2,J).EQ.species_name(I)) JSP2(J)=I
      enddo

      ! === ITYPE 14 - SURFACE REACTIONS
      if (ITYPE(J).EQ.14) then
        NPATH=0

        ! ------ Check for branching
        do K=1,nb_reactions
          if (((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(2,J).EQ.SYMBOL(2,K))).OR.&
          ((SYMBOL(2,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(1,J).EQ.SYMBOL(2,K)))) then
          if (SYMBOL(4,K)(:1).EQ.'J          ') NPATH=NPATH+1
        endif
      enddo

      ! ------ Branching ratio
      if (NPATH.EQ.0) then
        XJ(J)=0.d0
      else
        XJ(J)=XJ(J)/dble(NPATH)
      endif

      ! ------ Factor of 2 for same species reactions
      if (JSP1(J).EQ.JSP2(J)) XJ(J)=XJ(J)/2.0d0

      ! ------ Calculate evaporation fraction
      NEVAP=0
      do K=1,nb_reactions
        if ((SYMBOL(4,J)(:1).EQ.'J          ').AND.(A(K).NE.0.d0)) then
          if ((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(2,J).EQ.SYMBOL(2,K)).AND.&
          (SYMBOL(4,J)(2:).EQ.SYMBOL(4,K)).AND.&
          (SYMBOL(5,J)(2:).EQ.SYMBOL(5,K)).AND.&
          (SYMBOL(6,J)(2:).EQ.SYMBOL(6,K)).AND.&
          (SYMBOL(4,K)(:1).NE.'J          ')) NEVAP=NEVAP+1
        endif
        if ((SYMBOL(4,J)(:1).NE.'J          ').AND.(A(J).NE.0.d0)) then
          if ((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(2,J).EQ.SYMBOL(2,K)).AND.&
          (SYMBOL(4,J).EQ.SYMBOL(4,K)(2:)).AND.&
          (SYMBOL(5,J).EQ.SYMBOL(5,K)(2:)).AND.&
          (SYMBOL(6,J).EQ.SYMBOL(6,K)(2:)).AND.&
          (SYMBOL(4,K)(:1).EQ.'J          ')) NEVAP=NEVAP+1
        endif
      enddo

      N4=0
      N5=0
      N6=0
      do I=nb_gaseous_species+1,nb_species
        if (SYMBOL(4,J)(:1).EQ.'J          ') then
          if (SYMBOL(4,J).EQ.species_name(I)) N4=I
          if (SYMBOL(5,J).EQ.species_name(I)) N5=I
          if (SYMBOL(6,J).EQ.species_name(I)) N6=I
        endif
        if ((SYMBOL(4,J)(:1).NE.'J          ').AND.&
        (SYMBOL(4,J)(:1).NE.'X          ')) then
        if (SYMBOL(4,J).EQ.species_name(I)(2:)) N4=I
        if (SYMBOL(5,J).EQ.species_name(I)(2:)) N5=I
        if (SYMBOL(6,J).EQ.species_name(I)(2:)) N6=I
      endif
    enddo

    DHFSUM=DHF(JSP1(J))+DHF(JSP2(J))-DHF(N4)
    if (N5.NE.0) DHFSUM=DHFSUM-DHF(N5)
    if (N6.NE.0) DHFSUM=DHFSUM-DHF(N6)
    ! ------ Convert from kcal to J, from J to K
    DHFSUM=DHFSUM*4.184D+03/1.38054D-23
    ! ------ Convert from #moles-1 to #reactions-1
    DHFSUM=DHFSUM/AVOGADRO

    DHFSUM=DHFSUM+EA(J)

    SUM1=ED(N4)
    if (N5.NE.0) SUM1=MAX(ED(N4),ED(N5))
    if (N6.NE.0) SUM1=MAX(ED(N4),ED(N5),ED(N6))

    ATOMS=0
    do K=1,NB_PRIME_ELEMENTS
      ATOMS=ATOMS+species_composition(K,N4)
    enddo

    SUM2=1.0d0-(SUM1/DHFSUM)
    if (ATOMS.EQ.2) SUM2=SUM2**(3*ATOMS-5)
    if (ATOMS.GT.2) SUM2=SUM2**(3*ATOMS-6)
    !         SUM2=SUM2**(3*ATOMS-6)
    SUM2=ARRK*SUM2
    EVFRAC=SUM2/(1+SUM2)

    !        V.W. Jul 2006 f=evfrac=0.009 for H2O (Kroes & Anderson 2006) 
    !         if (SYMBOL(4,J).EQ.'H2O     ') then
    !                 EVFRAC=0.009
    !                EVFRAC_H2O=0.009
    !         endif

    BADFLAG=0
    if (DHF(JSP1(J)).LE.-999.0) then
      EVFRAC=0.d0
      BADFLAG=BADFLAG+1
    endif
    if (DHF(JSP2(J)).LE.-999.0) then
      EVFRAC=0.d0
      BADFLAG=BADFLAG+1
    endif
    if (DHF(N4).LE.-999.0) then
      EVFRAC=0.d0
      BADFLAG=BADFLAG+1
    endif
    if (N5.NE.0) then
      EVFRAC=0.d0
      BADFLAG=BADFLAG+1
    endif
    if (N6.NE.0) then
      EVFRAC=0.d0
      BADFLAG=BADFLAG+1
    endif

    if (EVFRAC.GE.1.0d0) EVFRAC=1.0d0
    if (EVFRAC.LE.0.d0) EVFRAC=0.d0
    if (NEVAP.EQ.0) EVFRAC=0.d0
    if (DHFSUM.LE.0.d0) EVFRAC=0.d0

    if (SYMBOL(4,J)(:1).EQ.'J          ') then
      EVFRAC=1.0d0-EVFRAC
    endif

    XJ(J)=XJ(J)*EVFRAC

    ! ------ Calculate quantum activation energy
    REDMAS = SMASS(JSP1(J)) * SMASS(JSP2(J)) / (SMASS(JSP1(J)) + SMASS(JSP2(J)))
    ACT1(J) = 2.0d0 * ACT/H_BARRE * SQRT(2.0d0*AMU*REDMAS*K_B*EA(J))
  endif

  ! === ITYPE 16 - C.R. DESORPTION
  if (ITYPE(J).EQ.16) then
    if (SMASS(JSP1(J)).EQ.0) XJ(J)=0.d0
  endif

  ! === ITYPE 99 - ACCRETION
  if (ITYPE(J).EQ.99) then
    ! ------ Save tag of resultant grain surface species
    do I=1,nb_species
      if (SYMBOL(4,J).EQ.species_name(I)) JSP2(J)=I
    enddo
  endif

enddo

! === Zero dummy H2 formation rxns, if necc.
!      if (IS_GRAIN_REACTIONS.NE.0) then
!         XJ(1)=0.d0
!         XJ(2)=0.d0
!      endif

return
end subroutine init_reaction_rates

! ======================================================================
! ======================================================================
END program gasgrain
