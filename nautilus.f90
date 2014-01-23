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
! Files:
! nautilus.f90          The main program
! nls_control.d         Control parameters
! nls_header_mod.f90    Declarations and constants
! nls_io.f90            Input/output
! nls_odes.f90          Rate equations and jacobian
! nls_phys_1D.f90       1D structure (overwrites some parameters)
! nls_mesh.f90          1D mesh
! nls_diffusion.f90     Turbulent diffusion routines
!
! opkd*.f               Odepack files
! 
! nls_surf_input.d      Surface parameters (desorption energies, etc..)
! gg_react_osu_03_2008  Chemical network
! gg_photodiss*         H2 and CO photodissociation
! nls_init.d            Initial abundances (must be present even if IREAD=0)
!
! nlso_srates.d         Output of some surface quantities (incomplete)
! nlso_tail.d           Chemical compasition at the last timestep (if used in 0D)
! output_1D.******      1D outputs (binaries) 
! rates1D.******        Outputs of the rate coefficients for 1 specific mesh point
! nlso_mod.d            Output of the modified rate coefficients (not checked)
! nlso_spec.d           List of species with their names and numbers
! 
! April 2011 VW  An appoximative calculation of the X-ray ionization rate 
! have been added according to Glasgold at al. (1999)
! 
! -----------------------------------------------------------------------

PROGRAM Gasgrain

use global_variables

use diffusion
use input_output
use model_1D
use ode_solver
use iso_fortran_env

implicit none

integer :: lrw, liw

real(double_precision), dimension(:), allocatable :: Y ! nb_species

integer :: itol, itask, istate, iopt, mf
real(double_precision) :: atol
real(double_precision) :: T, TOUT, TIN

data itol, itask, istate, iopt, mf, atol/2,1,1,1,021,1.d-99/

call get_array_sizes()

allocate(Y(nb_species))

call initialize_global_arrays()

call read_input_files()

call write_parameters()

call read_abundances()

call get_gas_surface_species()

! Dimension of the work arrays for the solver 
! The number of non zero values is checked with the testjac flag
! NJAC should be around the largest printed value

if (testjac.eq.1) then
  write(*,*) '--------------------------------------------'
  write(*,*) 'Dummy run to check'
  write(*,*) 'the number of non zero elements per column !'
  write(*,*) '--------------------------------------------'
  lrw = 20 + 9*nb_species*NPTMAX + 10000*nb_species*NPTMAX
  liw = 31 + nb_species*NPTMAX + 10000*nb_species*NPTMAX
else
  lrw = 20 + 3 * NJAC*nb_species + 21 * nb_species
  liw = 31 + 3 * NJAC*nb_species + 21 * nb_species
endif

! Build spatial mesh 
call mesh()

! Initialization of elemental/chemical quantities

call initial()

timestep=0
call start(TOUT)

! 1D physical structure (nls_phys_1D)
call phys_1D()

! Write species name/index correspondance
call write_species()

! Global variable for the current spatial point
iptstore=1

! Initialize indices of reactants and products 

call chemsetup()

! Initializing ZXN
do ipts=1,nptmax
  ZXN(:,ipts) = XN(:)
enddo

! Initializing T
! T = local, TIME = global
T = 0.d0
TIME = 0.d0

! Allocate the JA array
allocate(JA(1:liw))

! The real time loop
do while (t.lt.0.9*STOP_TIME)

  call ztimestep() ! Determination of the diffusive timestep, modification of the global parameter zdt
  TIME=T+ZDT ! Final time for chemistry T -> T + ZDT

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

    Y(:nb_species) = ZXN(:,ipts)
    XN(:nb_species) = ZXN(:,ipts)

    call evolve (T,Y,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

    ! Output of the rates once every 10 chemical outputs
    !      if ((mod(it,wstepr).eq.0).and.(ipts.eq.irateout)) then
    call write_current_rates()
    !      endif

    if (istate.eq.-3) stop

    ZXN(:,ipts) = XN(:) ! putting back

  enddo ! end of the spatial loop for chemistry 

  ! Generic call to zdiffusion, a subroutine calling the chosen numerical scheme 
  call zdiffusion() ! diffusion 

  if (mod(timestep,wstep).eq.0) then
    call write_current_output()
  endif

enddo

if (nptmax.eq.1) call write_abundances('abundances.tmp')

stop

contains

! ======================================================================
! ======================================================================
subroutine INITIAL()
use global_variables

implicit none

! Locals
real(double_precision), dimension(NEMAX) :: MASS
real(double_precision) :: MSUM
integer :: ILAB, KSUM, j, k, i, isptemp

! Set elements' characteristics=========================================
! --- Find the atomic species associated with a given element

ILAB=1
do J=1,nb_species
  KSUM=0
  ! ------ Calculate species' mass
  do K=1,NEMAX
    KSUM=KSUM+IELM(K,J)
  enddo
  ! ------ Check for atomic species
  if ((KSUM.EQ.1).AND.(ICG(J).EQ.0).AND.&
  (SPEC(J)(:1).NE.'J          ').AND.(SPEC(J)(:1).NE.'X          ')) then
  if (ILAB.GT.NEMAX) then
    STOP '***More elements than NEMAX***'
  endif       
  ! --------- Save species number
  ISPELM(ILAB)=J
  ILAB=ILAB+1
endif

! ------ Check for electron species number
IF (SPEC(J).EQ.'e-         ') ISPE=J
enddo

! --- Re-arrange order of elements to match IELM columns (reactions file)
do J=1,NEMAX-1
  if (IELM(J,ISPELM(J)).NE.1) then
    do K=J+1,NEMAX
      if (IELM(J,ISPELM(K)).EQ.1) then
        ISPTEMP=ISPELM(K)
        ISPELM(K)=ISPELM(J)
        ISPELM(J)=ISPTEMP
      endif
    enddo
  endif
enddo

! --- Set elements' masses
do I=1,NEMAX
  if (SPEC(ISPELM(I)).EQ.'H          ') MASS(I)=1.d0
  if (SPEC(ISPELM(I)).EQ.'D          ') MASS(I)=2.d0
  if (SPEC(ISPELM(I)).EQ.'He         ') MASS(I)=4.d0
  if (SPEC(ISPELM(I)).EQ.'C          ') MASS(I)=12.d0
  if (SPEC(ISPELM(I)).EQ.'N          ') MASS(I)=14.d0
  if (SPEC(ISPELM(I)).EQ.'O          ') MASS(I)=16.d0
  if (SPEC(ISPELM(I)).EQ.'Na         ') MASS(I)=23.d0
  if (SPEC(ISPELM(I)).EQ.'Mg         ') MASS(I)=24.d0
  if (SPEC(ISPELM(I)).EQ.'Si         ') MASS(I)=28.d0
  if (SPEC(ISPELM(I)).EQ.'P          ') MASS(I)=31.d0
  if (SPEC(ISPELM(I)).EQ.'S          ') MASS(I)=32.d0
  if (SPEC(ISPELM(I)).EQ.'Cl         ') MASS(I)=35.d0
  if (SPEC(ISPELM(I)).EQ.'Fe         ') MASS(I)=56.d0
  if (SPEC(ISPELM(I)).EQ.'F          ') MASS(I)=19.d0
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
  ! ------ Calculate masses
  MSUM=0.d0
  do K=1,NEMAX 
    MSUM=MSUM+MASS(K)*IELM(K,I) 
  enddo 
  AWT(I)=MSUM
  if (SPEC(I).EQ.YE) AWT(I)=1.D+0/1836.D+0 
  if (SPEC(I).EQ.YGRAIN .OR. SPEC(I).EQ.'GRAIN-      ')&
  AWT(I)=4.0*PI*grain_radius*grain_radius*grain_radius*GRAIN_DENSITY/3.0/AMU
enddo

! Initialize the Av/NH ratio
! Can be scaled for different dust/gas ratios
! Here initial_dtg_mass_ratio is the original dust/gas mass ratio (from nls_control.d)
! initial_dtg_mass_ratio is changed later into the dust/hydrogen mass ratio

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
  if (SPEC(i).eq.YH2) INDH2=i
  if (SPEC(i).eq.YCO) INDCO=i
  if (SPEC(i).eq.YHE) INDHE=i
enddo

! Compute nb_sites_per_grain = number of sites per grain
nb_sites_per_grain = SITE_DENSITY*4.d0*PI*grain_radius**2

! Initialise reaction rates=============================================
call grainrate()

return 
end subroutine initial

! ======================================================================
! ======================================================================
subroutine start(tout)

use global_variables

implicit none

! Output
real(double_precision), intent(out) :: TOUT

! Locals
character(len=11), dimension(nb_species) :: SREAD
real(double_precision), dimension(nb_species) :: Y
integer :: i, j, k

! Initialise times======================================================
TOUT=0.0d0
TIME=0.0d0

! Compute elemental abundances

do J=1,NEMAX
  ELEMS(J)=0.0D+0
enddo
do I=1,nb_species 
  do J=1,NEMAX
    ELEMS(J)=ELEMS(J)+IELM(J,I)*XN(I)
  enddo
enddo

! Recompute initial_dtg_mass_ratio to remove He
! In the following, initial_dtg_mass_ratio is used as a H/dust mass ratio
do i=1,nemax
  if (SPEC(ISPELM(I)).EQ.YHE) then
    initial_dtg_mass_ratio = initial_dtg_mass_ratio*(1.d0+4*ELEMS(I))
    ! Mean molecular weight (cgs) 
    ! Approximated here (the exact calculus would require a sume over AWT
    ! Used for the diffusion in disks, not for chemistry
    MEANW = 2.d0 + 4.d0*ELEMS(I)/(1.d0+ELEMS(I))
  endif
enddo

! Compute the grain abundance

GTODN=(4.D+0*PI*GRAIN_DENSITY*grain_radius*grain_radius*grain_radius)/(3.D+0*initial_dtg_mass_ratio*AMU)

where(SPEC.EQ.YGRAIN) XN=1.0/GTODN

  ! Set the electron abundance via conservation===========
  ! And check at the same time that nls_init has the same elemental abundance
  ! as nls_control

  Y(:)=XN(:)
  call conserve(Y)
  XN(:)=Y(:) 

  return 
  end subroutine start

  ! ======================================================================
  ! ======================================================================

  subroutine evolve(T,Y,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

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
  real(double_precision), intent(out), dimension(nb_species) :: Y
 
  ! Locals
  integer, dimension(liw) :: IWORK
  real(double_precision), dimension(lrw) :: RWORK
  real(double_precision), dimension(nb_species) :: DUMMYPDJ, DUMMYY
  integer IDUMMY
  integer :: i
  real(double_precision), dimension(nb_species) :: satol

  real(double_precision) :: TIN

  integer :: NNZ

  !EXTERNAL FCHEM,JAC


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

  Y(:) = XN(:)

  ! if testjac is 1, print non zero elements per column of the Jacobian
  ! Done in odes/JACVW

  if (TESTJAC.eq.1) then
    DUMMYY=1.d-5
    call ratcon()
    do IDUMMY=1,nb_species
      call jacvw(dummyy,idummy,dummypdj)
    enddo
    STOP
  endif

  call shieldingsetup()

  ! Don't stop running, Forrest

  do while (t.lt.tout)

    istate = 1

    ! Adaptive absolute tolerance to avoid too high precision on very abundant species,
    ! H2 for instance. Helps running a bit faster

    do i=1,nb_species
      satol(i)=max(atol,1.d-16*y(i))
    enddo

    ! ratcond is already called in compute IAJA
    !      call ratcon(Y)

    ! Feed IWORK with IA and JA

    call computeIAJA(Y)
    NNZ=IA(nb_species+1)-1
    iwork(30+1:30+nb_species+1)=IA(1:nb_species+1)
    iwork(31+nb_species+1:31+nb_species+NNZ)=JA(1:NNZ)

    call dlsodes(fchem,nb_species,y,t,tout,itol,RELATIVE_TOLERANCE,satol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)       

    ! Whenever the solver fails converging, print the reason.
    ! cf odpkdmain.f for translation
    if (istate.ne.2) write(*,*)  'IPTS = ', ipts, 'ISTATE = ', ISTATE

    call conserve(Y)

    ! Stop, Forrest
  enddo

  T = TOUT + TIN 

  XN=Y

  return 
  end subroutine evolve




  ! ======================================================================
  ! ======================================================================
  subroutine conserve(Y)
  use global_variables
  
  implicit none
  
  ! Inputs/outputs
  real(double_precision), intent(inout), dimension(nb_species) :: Y
  
  ! Locals
  real(double_precision), dimension(NEMAX) :: ELMSUM
  real(double_precision) :: CHASUM

  integer :: i, k

  ! Prevent too low abundances
  do i=1,nb_species
    if (Y(i).le.1.d-99) Y(i)=1.d-99
  enddo

  ! --- Conserve electrons
  CHASUM=0.0D+0
  do I=1,nb_species
    if (I.NE.ISPE) CHASUM=CHASUM+ICG(I)*Y(I)
  enddo
  if (CHASUM.LE.0.0D+0) CHASUM=XNMIN
  Y(ISPE)=CHASUM

  ! --- Conserve other elements if selected
  if (ICONS.GT.0) then
    do K=1,ICONS
      ELMSUM(K)=0.0D+0
    enddo
    do I=1,nb_species
      do K=1,ICONS
        if (I.NE.ISPELM(K)) ELMSUM(K)=ELMSUM(K)+IELM(K,I)*Y(I)
      enddo
    enddo
    do K=1,ICONS
      Y(ISPELM(K))=ELEMS(K)-ELMSUM(K)
      if (Y(ISPELM(K)).LE.0.0D+0) Y(ISPELM(K))=XNMIN
    enddo
  endif

  ! Check for conservation
  ELMSUM(:)=0.0D+0
  do I=1,nb_species
    do K=1,NEMAX
      ELMSUM(K)=ELMSUM(K)+IELM(K,I)*Y(I)
    enddo
  enddo

  do k=1,nemax
    if (abs(ELEMS(K)-ELMSUM(K))/ELEMS(K).ge.0.01d0) then 
      write(*,*)  'CAUTION : Element ',SPEC(ISPELM(K)), 'is not conserved'
      write(*,*)  'Relative difference: ', abs(ELEMS(K)-ELMSUM(K))/ELEMS(K)
    endif
    if (SPEC(ISPELM(K)).eq.YH) then
      if (abs(ELEMS(K)-Y(INDH2)*2.D0)/ELEMS(K).ge.0.01d0) write(*,*) 'H is too depleted on the grains !!!!'
    endif
    if (SPEC(ISPELM(K)).eq.YHE) then
      if (abs(ELEMS(K)-Y(INDHE))/ELEMS(K).ge.0.01d0) write(*,*) 'He is too depleted on the grains !!!!'
    endif       
  enddo

  ! VW fev 2012 add a test for the helium and H2 abundance in the gas phase
  ! prevent excessive depletion


  return
  end subroutine CONSERVE

  ! ======================================================================
  ! Dummy jacobians when the code is run with mf=222
  ! Not to use unless, the solver has big problems converging
  ! ======================================================================
  !      subroutine DUMMY 
  !      implicit none
  !      integer :: N,J
  !      real(double_precision) :: T,Y,IAN, JAN, PDJ 

  !      entry jac (N,T,Y,ML,MU,PD,NROWPD)
  !      entry jac (N, T, Y, J, IAN, JAN, PDJ)
  !      return

  !      END



  ! ======================================================================
  ! Computes the H2 and CO column density for their self-shielding
  ! for each ipts, the column density is incremented recursively
  ! from the results for NH2 and NCO computed by RATCON2 the spatial step before
  ! NB: ZN is shifted with respect to N
  ! ======================================================================
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



    ! ======================================================================
    ! ======================================================================
  subroutine GRAINRATE()
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

    ! Set accretion rate====================================================

    ! COND is used to calculate R_acc = (sigma_d) * <v_i> * n_i * n_d
    ! Uses 'Mean' Maxwellian speed, rather than RMS (is this justifiable?)

    COND=PI*grain_radius*grain_radius*SQRT(8.0d0*K_B/PI/AMU)

    ! --- Evaluate sticking coeff and accretion rate factor for each species
    STICK=0.0D+0
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
      !         if (SPEC(I).EQ.YH2)      STICK=0.D+0 
      !         if (SPEC(I).EQ.YHE)      STICK=0.D+0 
      !         if (SPEC(I).EQ.YH)       STICK=0.D+0 
      if (SPEC(I).EQ.YHEP)     STICK=0.D+0 
      if (SPEC(I).EQ.'e-         ')      STICK=0.D+0
      if (SPEC(I).EQ.'H+         ')     STICK=0.D+0
      if (SPEC(I).EQ.YGRAIN)   STICK=0.D+0
      if (SPEC(I).EQ.'GRAIN-     ') STICK=0.D+0
      !         if (SPEC(I).EQ.'H-')     STICK=0.D+0
      !         if (SPEC(I).EQ.'H2+')    STICK=0.D+0

      if (I.GT.nb_gaseous_species) STICK=0.D+0
      CONDSP(I)=COND*STICK/SQRT(AWT(I))
    enddo


    ! Read in molecular information for surface rates=======================
    open(unit=10, file='nls_surf_fev2012.dat', status='OLD')
    read(10,*)

    ! --- Read info into dummy arrays
    NGS=0
    do I=1,nb_species
      700    CONTINUE
      read(10,'(A11,I4,F7.0,F6.0,D8.1,27X,F8.2)') GSPEC(I),INT1(I),REA1(I),REA2(I),REA3(I),REA4(I)
      if (GSPEC(I).EQ.'X          ') GOTO 700
      if (GSPEC(I).EQ.'           ') EXIT
      NGS=NGS+1
    enddo


    ! --- Read activation energies into dummy arrays
    read(10,'(I4)') NEA
    do J=1,NEA
      read(10,'(5A11,D9.2)') (GSread(L,J),L=1,5),REA5(J)
    enddo
    close(10)

    ! --- Transfer from dummies to arrays with correct species numbers
    do I=1,nb_species
      SMASS(I)=0
      ED(I)=0.0D+0
      EB(I)=0.0D+0
      DEB(I)=0.0D+0
      DHF(I)=0.0D+0
      do J=1,NGS
        if (SPEC(I).EQ.GSPEC(J)) then
          SMASS(I)=INT1(J)
          ED(I)=REA1(J)
          EB(I)=REA2(J)
          DEB(I)=REA3(J)
          DHF(I)=REA4(J)
          if ((SPEC(I).NE.YJH).AND.(SPEC(I).NE.YJH2).AND.&
          (EBFAC.GE.0.0D+0)) EB(I)=EBFAC*ED(I)
        endif
      enddo
      !IF(SPEC(I) == 'JN2O2      ') write(*,*) ED(I)
    enddo

    do I=1,nb_reactions
      EA(I)=0.0D+0
      do J=1,NEA
        if (SYMBOL(4,I)(:1).EQ.'J          ') then
          if ((SYMBOL(1,I).EQ.GSread(1,J)).AND.&
          (SYMBOL(2,I).EQ.GSread(2,J)).AND.&
          (SYMBOL(4,I).EQ.GSread(3,J)).AND.&
          (SYMBOL(5,I).EQ.GSread(4,J)).AND.&
          (SYMBOL(6,I).EQ.GSread(5,J))) EA(I)=REA5(J)
        endif
        if (SYMBOL(4,I)(:1).NE.'J          ') then
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
      CHF(I)=0.0D+0
      RQ1(I)=0.0D+0
      RQ2(I)=0.0D+0
      ! ------ For species which have been assigned surface info, SMASS=/=0
      if (SMASS(I).NE.0) then
        SMA=REAL(SMASS(I))
        ! --------- Set characteristic frequency
        CHF(I)=SQRT(2.0D+0*K_B/PI/PI/AMU * SITE_DENSITY*ED(I)/SMA)
        ! --------- Set quantum rates
        if (DEB(I).GE.1.0D-38) then
          RQ1(I)=DEB(I)*K_B/4.0d0/H_BARRE/nb_sites_per_grain
        else
          RQ1(I)=0.0D+0
        endif
        RQ2(I)=CHF(I)/nb_sites_per_grain*&
        EXP(-2.0D+0*SITE_SPACING/H_BARRE*SQRT(2.0D+0*AMU*SMA*K_B*EB(I)))
      endif
    enddo

    ! === Cycle all reactions
    do J=1,nb_reactions

      ! ------ Initialise all XJ rate factors, and get species 1 & 2
      XJ(J)=1.0D+0
      JSP1(J)=0
      JSP2(J)=0
      do I=1,nb_species
        if (SYMBOL(1,J).EQ.SPEC(I)) JSP1(J)=I
        if (SYMBOL(2,J).EQ.SPEC(I)) JSP2(J)=I
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
        XJ(J)=0.0D+0
      else
        XJ(J)=XJ(J)/REAL(NPATH)
      endif

      ! ------ Factor of 2 for same species reactions
      if (JSP1(J).EQ.JSP2(J)) XJ(J)=XJ(J)/2.0D+0
      !        write(*,*) SYMBOL(1,J)
      !       write(*,*) SYMBOL(2,J)
      !       write(*,*) XJ(J)

      !       stop

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
          if (SYMBOL(4,J).EQ.SPEC(I)) N4=I
          if (SYMBOL(5,J).EQ.SPEC(I)) N5=I
          if (SYMBOL(6,J).EQ.SPEC(I)) N6=I
        endif
        if ((SYMBOL(4,J)(:1).NE.'J          ').AND.&
        (SYMBOL(4,J)(:1).NE.'X          ')) then
        if (SYMBOL(4,J).EQ.SPEC(I)(2:)) N4=I
        if (SYMBOL(5,J).EQ.SPEC(I)(2:)) N5=I
        if (SYMBOL(6,J).EQ.SPEC(I)(2:)) N6=I
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
    do K=1,NEMAX
      ATOMS=ATOMS+IELM(K,N4)
    enddo

    SUM2=1.0D+0-(SUM1/DHFSUM)
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
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    endif
    if (DHF(JSP2(J)).LE.-999.0) then
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    endif
    if (DHF(N4).LE.-999.0) then
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    endif
    if (N5.NE.0) then
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    endif
    if (N6.NE.0) then
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    endif

    if (EVFRAC.GE.1.0D+0) EVFRAC=1.0D+0
    if (EVFRAC.LE.0.0D+0) EVFRAC=0.0D+0
    if (NEVAP.EQ.0) EVFRAC=0.0D+0
    if (DHFSUM.LE.0.0D+0) EVFRAC=0.0D+0

    !         EVFRAC=0.0D+0

    if (SYMBOL(4,J)(:1).EQ.'J          ') then
      EVFRAC=1.0D+0-EVFRAC
    endif

    !        V.W. Jul 2006 f=evfrac=0.009 for H2O (Kroes & Anderson 2006) 
    !        if (SYMBOL(4,J).EQ.'JH2O    ') then
    !            EVFRAC=1.0D+0-EVFRAC_H2O
    !         endif

    XJ(J)=XJ(J)*EVFRAC
    !        write(*,*) SYMBOL(1,J)
    !       write(*,*) SYMBOL(2,J)
    !       write(*,*) XJ(J)

    !       stop

    ! ------ Calculate quantum activation energy
    REDMAS=REAL(SMASS(JSP1(J))*SMASS(JSP2(J)))/&
    REAL(SMASS(JSP1(J))+SMASS(JSP2(J)))
    ACT1(J)=2.0D+0 * ACT/H_BARRE * SQRT(2.0D+0*AMU*REDMAS*K_B*EA(J))
  endif

  ! === ITYPE 16 - C.R. DESORPTION
  if (ITYPE(J).EQ.16) then
    if (SMASS(JSP1(J)).EQ.0) XJ(J)=0.0D+0
  endif

  ! === ITYPE 99 - ACCRETION
  if (ITYPE(J).EQ.99) then
    ! ------ Save tag of resultant grain surface species
    do I=1,nb_species
      if (SYMBOL(4,J).EQ.SPEC(I)) JSP2(J)=I
    enddo
  endif

enddo

! === Zero dummy H2 formation rxns, if necc.
!      if (IS_GRAIN_REACTIONS.NE.0) then
!         XJ(1)=0.0D+0
!         XJ(2)=0.0D+0
!      endif

return
end subroutine grainrate

! ======================================================================
! ======================================================================
END program gasgrain
