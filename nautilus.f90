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
use constants
use diffusion
use input_output
use model_1D
use ode_solver

implicit none

integer :: lrw, liw

real(double_precision), dimension(NSMAX) :: Y

integer :: itol, itask, istate, iopt, mf
real(double_precision) :: atol
real(double_precision) :: T, TOUT, TIN

data itol, itask, istate, iopt, mf, atol/2,1,1,1,021,1.d-99/

CALL FILESET
CALL READINPUT

! Dimension of the work arrays for the solver 
! The number of non zero values is checked with the testjac flag
! NJAC should be around the largest printed value

if (testjac.eq.1) then
  write(*,*) '--------------------------------------------'
  write(*,*) 'Dummy run to check'
  write(*,*) 'the number of non zero elements per column !'
  write(*,*) '--------------------------------------------'
  lrw = 20 + 9*NSMAX*NPTMAX + 10000*NSMAX*NPTMAX
  liw = 31 + NSMAX*NPTMAX + 10000*NSMAX*NPTMAX
else
  lrw = 20 + 3 * NJAC*NSMAX + 21 * NSMAX
  liw = 31 + 3 * NJAC*NSMAX + 21 * NSMAX
endif

! Build spatial mesh 
call mesh

! Initialization of elemental/chemical quantities

CALL INITIAL

IT=0
CALL START(TOUT)

! 1D physical structure (nls_phys_1D)
call phys_1D

! Write species name/index correspondance
CALL WRITESPEC

! Global variable for the current spatial point
iptstore=1

! Initialize indices of reactants and products 

call chemsetup

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
do while (t.lt.0.9*tfinal)

  call ztimestep() ! Determination of the diffusive timestep, modification of the global parameter zdt
  TIME=T+ZDT ! Final time for chemistry T -> T + ZDT

  it=it+1

  ! Store the current time in TIN (for 1D calculation)
  TIN = T

  WRITE(NERR,5) 'IT=',IT,', TIME=',TIME/TYEAR,' yrs'
  5    FORMAT (A3,I5,A7,1PD10.3,A4)

  do ipts=1,nptmax ! Start of the spatial loop for chemistry

    iptstore = ipts

    ! T being changed in dlsode, needs to be defined again

    T=TIN

    ! Feed 1D physical structure
    TAU=TAU1D(ipts)
    ZETAX=ZETAX1D(ipts)
    ! XNT is the total density of H atoms (both H and H2)
    ! But DENS1D is the real number density
    ! Problem with the sound speed if NH > NH2 (cf phys_1D)
    XNT=2.*DENS1D(ipts)
    TEMP=TEMP1D(ipts)
    DTEMP=DTEMP1D(ipts)

    ! Chemical evolution for each spatial point

    Y(:nsmax) = ZXN(:,ipts)
    XN(:nsmax) = ZXN(:,ipts)

    call evolve (T,Y,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

    ! Output of the rates once every 10 chemical outputs
    !      if ((mod(it,wstepr).eq.0).and.(ipts.eq.irateout)) then
    call write_rates()
    !      endif

    if (istate.eq.-3) stop

    ZXN(:,ipts) = XN(:) ! putting back

  enddo ! end of the spatial loop for chemistry 

  ! Generic call to zdiffusion, a subroutine calling the chosen numerical scheme 
  call zdiffusion ! diffusion 

  if (mod(it,wstep).eq.0) then
    call write_outputs()
  endif

ENDDO

if (nptmax.eq.1) call writetail

stop

contains

! ======================================================================
! ======================================================================
SUBROUTINE FILESET
use constants
use global_variables
implicit none


OPEN (UNIT=NCON,FILE='nls_control.d',STATUS='OLD')
!~       OPEN (UNIT=NMOD,FILE='nlso_mod.d',STATUS='UNKNOWN')
OPEN (UNIT=NSP, FILE='nlso_spec.d',STATUS='UNKNOWN')
OPEN (UNIT=NGR, FILE='nls_surf_fev2012.dat',STATUS='OLD')
OPEN (UNIT=NJR, FILE='nls_gas_fev2012.dat',STATUS='OLD')
!      OPEN (UNIT=NJR, FILE='nls_gas_update.dat',STATUS='OLD')
OPEN (UNIT=NJR2, FILE='nls_grain_fev2012.dat',STATUS='OLD')
OPEN (UNIT=NTAI,FILE='nlso_tail.d',STATUS='UNKNOWN')
OPEN (UNIT=NINI,FILE='nls_init.d',STATUS='OLD') 
OPEN (UNIT=CODIS,FILE='gg_CO_Photodiss.d',STATUS='OLD')
OPEN (UNIT=H2DIS,FILE='gg_H2_Photodiss.d',STATUS='OLD')

RETURN 
END subroutine FILESET

! ======================================================================
! ======================================================================
SUBROUTINE INITIAL
use global_variables
use constants
implicit none

real(double_precision), dimension(NEMAX) :: MASS
real(double_precision) :: MSUM
integer :: ILAB, KSUM, j, k, i, isptemp

! Set elements' characteristics=========================================
! --- Find the atomic species associated with a given element

ILAB=1
DO J=1,NSMAX
  KSUM=0
  ! ------ Calculate species' mass
  DO K=1,NEMAX
    KSUM=KSUM+IELM(K,J)
  ENDDO
  ! ------ Check for atomic species
  IF ((KSUM.EQ.1).AND.(ICG(J).EQ.0).AND.&
  (SPEC(J)(:1).NE.'J          ').AND.(SPEC(J)(:1).NE.'X          ')) THEN
  IF (ILAB.GT.NEMAX) then
    STOP '***More elements than NEMAX***'
  endif       
  ! --------- Save species number
  ISPELM(ILAB)=J
  ILAB=ILAB+1
ENDIF

! ------ Check for electron species number
IF (SPEC(J).EQ.'e-         ') ISPE=J
ENDDO

! --- Re-arrange order of elements to match IELM columns (reactions file)
DO J=1,NEMAX-1
  IF (IELM(J,ISPELM(J)).NE.1) THEN
    DO K=J+1,NEMAX
      IF (IELM(J,ISPELM(K)).EQ.1) THEN
        ISPTEMP=ISPELM(K)
        ISPELM(K)=ISPELM(J)
        ISPELM(J)=ISPTEMP
      ENDIF
    ENDDO
  ENDIF
ENDDO

! --- Set elements' masses
DO I=1,NEMAX
  IF (SPEC(ISPELM(I)).EQ.'H          ') MASS(I)=1.d0
  IF (SPEC(ISPELM(I)).EQ.'D          ') MASS(I)=2.d0
  IF (SPEC(ISPELM(I)).EQ.'He         ') MASS(I)=4.d0
  IF (SPEC(ISPELM(I)).EQ.'C          ') MASS(I)=12.d0
  IF (SPEC(ISPELM(I)).EQ.'N          ') MASS(I)=14.d0
  IF (SPEC(ISPELM(I)).EQ.'O          ') MASS(I)=16.d0
  IF (SPEC(ISPELM(I)).EQ.'Na         ') MASS(I)=23.d0
  IF (SPEC(ISPELM(I)).EQ.'Mg         ') MASS(I)=24.d0
  IF (SPEC(ISPELM(I)).EQ.'Si         ') MASS(I)=28.d0
  IF (SPEC(ISPELM(I)).EQ.'P          ') MASS(I)=31.d0
  IF (SPEC(ISPELM(I)).EQ.'S          ') MASS(I)=32.d0
  IF (SPEC(ISPELM(I)).EQ.'Cl         ') MASS(I)=35.d0
  IF (SPEC(ISPELM(I)).EQ.'Fe         ') MASS(I)=56.d0
  IF (SPEC(ISPELM(I)).EQ.'F          ') MASS(I)=19.d0
ENDDO

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
DO I=1,NSMAX 
  ! ------ Calculate masses
  MSUM=0.d0
  DO K=1,NEMAX 
    MSUM=MSUM+MASS(K)*IELM(K,I) 
  ENDDO 
  AWT(I)=MSUM
  IF (SPEC(I).EQ.YE) AWT(I)=1.D+0/1836.D+0 
  IF (SPEC(I).EQ.YGRAIN .OR. SPEC(I).EQ.'GRAIN-      ')&
  AWT(I)=4.0*PI*RD*RD*RD*RHOD/3.0/AMU
ENDDO

! Initialize the Av/NH ratio
! Can be scaled for different dust/gas ratios
! Here DTOGM is the original dust/gas mass ratio (from nls_control.d)
! DTOGM is changed later into the dust/hydrogen mass ratio

RAVNH=5.34d-22*(DTOGM/1.d-2)

! Find ITYPE first and last reactions===================================
DO I=0,NITYPE-1
  IRXSTA(I)=0
  IRXFIN(I)=0
  DO J=1,NKMAX
    IF ((ITYPE(J).EQ.I).AND.(IRXSTA(I).EQ.0)) IRXSTA(I)=J
    IF (ITYPE(J).EQ.I) IRXFIN(I)=J
  ENDDO
ENDDO

! Find the index of CO and H2
do i=1,nsmax
  if (SPEC(i).eq.YH2) INDH2=i
  if (SPEC(i).eq.YCO) INDCO=i
  if (SPEC(i).eq.YHE) INDHE=i
enddo

! Compute TNS = number of sites per grain
TNS = SNS*4.d0*PI*RD**2

! Initialise reaction rates=============================================
CALL GRAINRATE

RETURN 
END SUBROUTINE INITIAL

! ======================================================================
! ======================================================================
SUBROUTINE START(TOUT)

use global_variables
use constants
implicit none

character(len=11), dimension(NSMAX) :: SREAD
real(double_precision), dimension(NSMAX) :: Y

integer :: i, j, k
real(double_precision) :: TOUT

! Initialise times======================================================
TOUT=0.0D+0
TIME=0.0D+0

! Set initial abundances================================================
DO I=1,NSMAX
  XN(I)=XNMIN
  DO K=1,NS0
    IF (SPEC(I).EQ.XS0(K)) THEN
      XN(I)=XN0(K)
    ENDIF
  ENDDO
ENDDO



! Compute elemental abundances

DO J=1,NEMAX
  ELEMS(J)=0.0D+0
ENDDO
DO I=1,NSMAX 
  DO J=1,NEMAX
    ELEMS(J)=ELEMS(J)+IELM(J,I)*XN(I)
  ENDDO
ENDDO

! Recompute DTOGM to remove He
! In the following, DTOGM is used as a H/dust mass ratio
do i=1,nemax
  IF (SPEC(ISPELM(I)).EQ.YHE) then
    DTOGM = DTOGM*(1.d0+4*ELEMS(I))
    ! Mean molecular weight (cgs) 
    ! Approximated here (the exact calculus would require a sume over AWT
    ! Used for the diffusion in disks, not for chemistry
    MEANW = 2.d0 + 4.d0*ELEMS(I)/(1.d0+ELEMS(I))
  endif
enddo

! Compute the grain abundance

GTODN=(4.D+0*PI*RHOD*RD*RD*RD)/(3.D+0*DTOGM*AMU)

where(SPEC.EQ.YGRAIN) XN=1.0/GTODN


  ! Read initial abundances if IREAD is switched on=======================
  IF (IREAD.NE.0) THEN

    ! Skip the first line on nls_init.d
    READ (NINI,*) 

    ! ------ Read species and abundances
    READ (NINI,110) (SREAD(I),XN(I),I=1,NSMAX)
    110 FORMAT (5(A11,2X,1PE12.5,2X))
    READ (NINI,*)
    ! ------ Check if species in nls_init.d correspond to the reaction file
    DO I=1,NSMAX
      IF (SREAD(I).NE.SPEC(I)) THEN
        WRITE (NERR,*) 'Input species in init file ',&
        'do not match those in reaction file'
        STOP
      ENDIF
    ENDDO
  ENDIF


  ! Set the electron abundance via conservation===========
  ! And check at the same time that nls_init has the same elemental abundance
  ! as nls_control

  Y(:)=XN(:)
  CALL CONSERVE(Y)
  XN(:)=Y(:) 

  RETURN 
  END subroutine start

  ! ======================================================================
  ! ======================================================================

  SUBROUTINE EVOLVE (T,Y,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

  use global_variables
  use constants
  implicit none

  integer :: liw,lrw

  integer, dimension(liw) :: IWORK
  real(double_precision), dimension(NSMAX) :: Y
  real(double_precision), dimension(lrw) :: RWORK
  real(double_precision), dimension(NSMAX) :: DUMMYPDJ, DUMMYY
  integer IDUMMY
  integer :: i
  integer :: itol, itask, istate, iopt, mf
  real(double_precision) :: atol
  real(double_precision), dimension(nsmax) :: satol

  real(double_precision) :: T, TOUT, TIN

  integer :: NNZ

  !EXTERNAL FCHEM,JAC


  ! Initialize work arrays

  iwork(:) = 0
  rwork(:) = 0.d0
  IWORK(5) = 5
  RWORK(6) = 3.154D14
  IWORK(6) = 10000
  IWORK(7) = 2

  if (it.eq.1) IWORK(6)=2000

  ! Changing the time to avoid T + DT = T 

  TOUT=TIME

  TIN = T
  TOUT = TOUT - TIN
  T = 0.d0      

  Y(:) = XN(:)

  ! If testjac is 1, print non zero elements per column of the Jacobian
  ! Done in odes/JACVW

  if (TESTJAC.eq.1) then
    DUMMYY=1.d-5
    call ratcon()
    do IDUMMY=1,NSMAX
      call JACVW(DUMMYY,IDUMMY,DUMMYPDJ)
    enddo
    STOP
  endif

  call shieldingsetup

  ! Don't stop running, Forrest

  do while (t.lt.tout)

    istate = 1

    ! Adaptive absolute tolerance to avoid too high precision on very abundant species,
    ! H2 for instance. Helps running a bit faster

    do i=1,nsmax
      satol(i)=max(atol,1.d-16*y(i))
    enddo

    ! ratcond is already called in compute IAJA
    !      call ratcon(Y)

    ! Feed IWORK with IA and JA

    call computeIAJA(Y)
    NNZ=IA(NSMAX+1)-1
    iwork(30+1:30+NSMAX+1)=IA(1:NSMAX+1)
    iwork(31+NSMAX+1:31+NSMAX+NNZ)=JA(1:NNZ)

    call dlsodes(fchem,nsmax,y,t,tout,itol,rtol,satol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)       

    ! Whenever the solver fails converging, print the reason.
    ! cf odpkdmain.f for translation
    if (istate.ne.2) write(*,*)  'IPTS = ', ipts, 'ISTATE = ', ISTATE

    CALL CONSERVE(Y)

    ! Stop, Forrest
  enddo

  T = TOUT + TIN 

  XN=Y

  RETURN 
  END SUBROUTINE EVOLVE




  ! ======================================================================
  ! ======================================================================
  SUBROUTINE CONSERVE(Y)
  use global_variables
  use constants
  implicit none

  real(double_precision), dimension(NSMAX) :: Y
  real(double_precision), dimension(NEMAX) :: ELMSUM
  real(double_precision) :: CHASUM

  integer :: i, k

  ! Prevent too low abundances
  do i=1,NSMAX
    if (Y(i).le.1.d-99) Y(i)=1.d-99
  enddo

  ! --- Conserve electrons
  CHASUM=0.0D+0
  DO I=1,NSMAX
    IF (I.NE.ISPE) CHASUM=CHASUM+ICG(I)*Y(I)
  ENDDO
  IF (CHASUM.LE.0.0D+0) CHASUM=XNMIN
  Y(ISPE)=CHASUM

  ! --- Conserve other elements if selected
  IF (ICONS.GT.0) THEN
    DO K=1,ICONS
      ELMSUM(K)=0.0D+0
    ENDDO
    DO I=1,NSMAX
      DO K=1,ICONS
        IF (I.NE.ISPELM(K)) ELMSUM(K)=ELMSUM(K)+IELM(K,I)*Y(I)
      ENDDO
    ENDDO
    DO K=1,ICONS
      Y(ISPELM(K))=ELEMS(K)-ELMSUM(K)
      IF (Y(ISPELM(K)).LE.0.0D+0) Y(ISPELM(K))=XNMIN
    ENDDO
  ENDIF

  ! Check for conservation
  ELMSUM(:)=0.0D+0
  DO I=1,NSMAX
    DO K=1,NEMAX
      ELMSUM(K)=ELMSUM(K)+IELM(K,I)*Y(I)
    ENDDO
  ENDDO

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


  RETURN
  END SUBROUTINE CONSERVE

  ! ======================================================================
  ! Dummy jacobians when the code is run with mf=222
  ! Not to use unless, the solver has big problems converging
  ! ======================================================================
  !      SUBROUTINE DUMMY 
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
  SUBROUTINE SHIELDINGSETUP

  use global_variables
  implicit none
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
  end SUBROUTINE SHIELDINGSETUP



    ! ======================================================================
    ! ======================================================================
    SUBROUTINE GRAINRATE
    use global_variables
    use constants
    implicit none

    real(double_precision), dimension(nsmax) :: REA1,REA2,REA3,REA4
    real(double_precision), dimension(nkmax) :: REA5
    real(double_precision), dimension(nsmax) :: SMASS
    real(double_precision) :: SMA,REDMAS,STICK,EVFRAC,DHFSUM,SUM1,SUM2
    integer, dimension(NKMAX) :: INT1
    integer :: NGS,NEA,NPATH,NEVAP,BADFLAG,ATOMS
    character(len=11), dimension(5,NKMAX) :: GSREAD
    character(len=11), dimension(nsmax) :: GSPEC

    real(double_precision) :: cond
    integer :: i, j,k,l,n4, n5, n6

    ! Set accretion rate====================================================

    ! COND is used to calculate R_acc = (sigma_d) * <v_i> * n_i * n_d
    ! Uses 'Mean' Maxwellian speed, rather than RMS (is this justifiable?)

    COND=PI*RD*RD*SQRT(8.0D+0*BOLTZ/PI/AMU)

    ! --- Evaluate sticking coeff and accretion rate factor for each species
    STICK=0.0D+0
    DO I=1,NSMAX
      IF (ICG(I).EQ.0) THEN
        STICK=STICK0
      ENDIF 
      IF (ICG(I).GT.0) THEN 
        STICK=STICKP 
      ENDIF
      IF (ICG(I).LT.0) THEN 
        STICK=STICKN 
      ENDIF
      !         IF (SPEC(I).EQ.YH2)      STICK=0.D+0 
      !         IF (SPEC(I).EQ.YHE)      STICK=0.D+0 
      !         IF (SPEC(I).EQ.YH)       STICK=0.D+0 
      IF (SPEC(I).EQ.YHEP)     STICK=0.D+0 
      IF (SPEC(I).EQ.'e-         ')      STICK=0.D+0
      IF (SPEC(I).EQ.'H+         ')     STICK=0.D+0
      IF (SPEC(I).EQ.YGRAIN)   STICK=0.D+0
      IF (SPEC(I).EQ.'GRAIN-     ') STICK=0.D+0
      !         IF (SPEC(I).EQ.'H-')     STICK=0.D+0
      !         IF (SPEC(I).EQ.'H2+')    STICK=0.D+0

      IF (I.GT.NSGAS) STICK=0.D+0
      CONDSP(I)=COND*STICK/SQRT(AWT(I))
    ENDDO

    ! Read in molecular information for surface rates=======================
    READ(NGR,*)

    ! --- Read info into dummy arrays
    NGS=0
    DO I=1,NSMAX
      700    CONTINUE
      READ(NGR,705) GSPEC(I),INT1(I),REA1(I),REA2(I),REA3(I),REA4(I)
      705    FORMAT(A11,I4,F7.0,F6.0,D8.1,27X,F8.2)
      IF (GSPEC(I).EQ.'X          ') GOTO 700
      IF (GSPEC(I).EQ.'           ') EXIT
      NGS=NGS+1
    ENDDO


    ! --- Read activation energies into dummy arrays
    READ(NGR,720) NEA
    720 FORMAT(I4)
    DO J=1,NEA
      READ(NGR,730) (GSREAD(L,J),L=1,5),REA5(J)
      730   FORMAT(5A11,D9.2)
    ENDDO

    ! --- Transfer from dummies to arrays with correct species numbers
    DO I=1,NSMAX
      SMASS(I)=0
      ED(I)=0.0D+0
      EB(I)=0.0D+0
      DEB(I)=0.0D+0
      DHF(I)=0.0D+0
      DO J=1,NGS
        IF (SPEC(I).EQ.GSPEC(J)) THEN
          SMASS(I)=INT1(J)
          ED(I)=REA1(J)
          EB(I)=REA2(J)
          DEB(I)=REA3(J)
          DHF(I)=REA4(J)
          IF ((SPEC(I).NE.YJH).AND.(SPEC(I).NE.YJH2).AND.&
          (EBFAC.GE.0.0D+0)) EB(I)=EBFAC*ED(I)
        ENDIF
      ENDDO
      !IF(SPEC(I) == 'JN2O2      ') write(*,*) ED(I)
    ENDDO

    DO I=1,NKMAX
      EA(I)=0.0D+0
      DO J=1,NEA
        IF (SYMBOL(4,I)(:1).EQ.'J          ') THEN
          IF ((SYMBOL(1,I).EQ.GSREAD(1,J)).AND.&
          (SYMBOL(2,I).EQ.GSREAD(2,J)).AND.&
          (SYMBOL(4,I).EQ.GSREAD(3,J)).AND.&
          (SYMBOL(5,I).EQ.GSREAD(4,J)).AND.&
          (SYMBOL(6,I).EQ.GSREAD(5,J))) EA(I)=REA5(J)
        ENDIF
        IF (SYMBOL(4,I)(:1).NE.'J          ') THEN
          IF ((SYMBOL(1,I).EQ.GSREAD(1,J)).AND.&
          (SYMBOL(2,I).EQ.GSREAD(2,J)).AND.&
          (SYMBOL(4,I).EQ.GSREAD(3,J)(2:)).AND.&
          (SYMBOL(5,I).EQ.GSREAD(4,J)(2:)).AND.&
          (SYMBOL(6,I).EQ.GSREAD(5,J)(2:))) EA(I)=REA5(J)
        ENDIF
      ENDDO
      !IF(symbol(4,i) == 'JO2H       ') write(*,*)  symbol(:,i), Ea(i)
    ENDDO

    ! Set up constants, quantum rate info===================================
    DO I=1,NSMAX
      CHF(I)=0.0D+0
      RQ1(I)=0.0D+0
      RQ2(I)=0.0D+0
      ! ------ For species which have been assigned surface info, SMASS=/=0
      IF (SMASS(I).NE.0) THEN
        SMA=REAL(SMASS(I))
        ! --------- Set characteristic frequency
        CHF(I)=SQRT(2.0D+0*BOLTZ/PI/PI/AMU * SNS*ED(I)/SMA)
        ! --------- Set quantum rates
        IF (DEB(I).GE.1.0D-38) THEN
          RQ1(I)=DEB(I)*BOLTZ/4.0D+0/HBAR/TNS
        ELSE
          RQ1(I)=0.0D+0
        ENDIF
        RQ2(I)=CHF(I)/TNS*&
        EXP(-2.0D+0*ACM/HBAR*SQRT(2.0D+0*AMU*SMA*BOLTZ*EB(I)))
      ENDIF
    ENDDO

    ! === Cycle all reactions
    DO J=1,NKMAX

      ! ------ Initialise all XJ rate factors, and get species 1 & 2
      XJ(J)=1.0D+0
      JSP1(J)=0
      JSP2(J)=0
      DO I=1,NSMAX
        IF (SYMBOL(1,J).EQ.SPEC(I)) JSP1(J)=I
        IF (SYMBOL(2,J).EQ.SPEC(I)) JSP2(J)=I
      ENDDO

      ! === ITYPE 14 - SURFACE REACTIONS
      IF (ITYPE(J).EQ.14) THEN
        NPATH=0

        ! ------ Check for branching
        DO K=1,NKMAX
          IF (((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(2,J).EQ.SYMBOL(2,K))).OR.&
          ((SYMBOL(2,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(1,J).EQ.SYMBOL(2,K)))) THEN
          IF (SYMBOL(4,K)(:1).EQ.'J          ') NPATH=NPATH+1
        ENDIF
      ENDDO

      ! ------ Branching ratio
      IF (NPATH.EQ.0) THEN
        XJ(J)=0.0D+0
      ELSE
        XJ(J)=XJ(J)/REAL(NPATH)
      ENDIF

      ! ------ Factor of 2 for same species reactions
      IF (JSP1(J).EQ.JSP2(J)) XJ(J)=XJ(J)/2.0D+0
      !        write(*,*) SYMBOL(1,J)
      !       write(*,*) SYMBOL(2,J)
      !       write(*,*) XJ(J)

      !       stop

      ! ------ Calculate evaporation fraction
      NEVAP=0
      DO K=1,NKMAX
        IF ((SYMBOL(4,J)(:1).EQ.'J          ').AND.(A(K).NE.0.d0)) THEN
          IF ((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(2,J).EQ.SYMBOL(2,K)).AND.&
          (SYMBOL(4,J)(2:).EQ.SYMBOL(4,K)).AND.&
          (SYMBOL(5,J)(2:).EQ.SYMBOL(5,K)).AND.&
          (SYMBOL(6,J)(2:).EQ.SYMBOL(6,K)).AND.&
          (SYMBOL(4,K)(:1).NE.'J          ')) NEVAP=NEVAP+1
        ENDIF
        IF ((SYMBOL(4,J)(:1).NE.'J          ').AND.(A(J).NE.0.d0)) THEN
          IF ((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
          (SYMBOL(2,J).EQ.SYMBOL(2,K)).AND.&
          (SYMBOL(4,J).EQ.SYMBOL(4,K)(2:)).AND.&
          (SYMBOL(5,J).EQ.SYMBOL(5,K)(2:)).AND.&
          (SYMBOL(6,J).EQ.SYMBOL(6,K)(2:)).AND.&
          (SYMBOL(4,K)(:1).EQ.'J          ')) NEVAP=NEVAP+1
        ENDIF
      ENDDO

      N4=0
      N5=0
      N6=0
      DO I=NSGAS+1,NSMAX
        IF (SYMBOL(4,J)(:1).EQ.'J          ') THEN
          IF (SYMBOL(4,J).EQ.SPEC(I)) N4=I
          IF (SYMBOL(5,J).EQ.SPEC(I)) N5=I
          IF (SYMBOL(6,J).EQ.SPEC(I)) N6=I
        ENDIF
        IF ((SYMBOL(4,J)(:1).NE.'J          ').AND.&
        (SYMBOL(4,J)(:1).NE.'X          ')) THEN
        IF (SYMBOL(4,J).EQ.SPEC(I)(2:)) N4=I
        IF (SYMBOL(5,J).EQ.SPEC(I)(2:)) N5=I
        IF (SYMBOL(6,J).EQ.SPEC(I)(2:)) N6=I
      ENDIF
    ENDDO

    DHFSUM=DHF(JSP1(J))+DHF(JSP2(J))-DHF(N4)
    IF (N5.NE.0) DHFSUM=DHFSUM-DHF(N5)
    IF (N6.NE.0) DHFSUM=DHFSUM-DHF(N6)
    ! ------ Convert from kcal to J, from J to K
    DHFSUM=DHFSUM*4.184D+03/1.38054D-23
    ! ------ Convert from #moles-1 to #reactions-1
    DHFSUM=DHFSUM/AVOGADRO

    DHFSUM=DHFSUM+EA(J)

    SUM1=ED(N4)
    IF (N5.NE.0) SUM1=MAX(ED(N4),ED(N5))
    IF (N6.NE.0) SUM1=MAX(ED(N4),ED(N5),ED(N6))

    ATOMS=0
    DO K=1,NEMAX
      ATOMS=ATOMS+IELM(K,N4)
    ENDDO

    SUM2=1.0D+0-(SUM1/DHFSUM)
    IF (ATOMS.EQ.2) SUM2=SUM2**(3*ATOMS-5)
    IF (ATOMS.GT.2) SUM2=SUM2**(3*ATOMS-6)
    !         SUM2=SUM2**(3*ATOMS-6)
    SUM2=ARRK*SUM2
    EVFRAC=SUM2/(1+SUM2)

    !        V.W. Jul 2006 f=evfrac=0.009 for H2O (Kroes & Anderson 2006) 
    !         IF (SYMBOL(4,J).EQ.'H2O     ') THEN
    !                 EVFRAC=0.009
    !                EVFRAC_H2O=0.009
    !         ENDIF

    BADFLAG=0
    IF (DHF(JSP1(J)).LE.-999.0) THEN
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    ENDIF
    IF (DHF(JSP2(J)).LE.-999.0) THEN
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    ENDIF
    IF (DHF(N4).LE.-999.0) THEN
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    ENDIF
    IF (N5.NE.0) THEN
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    ENDIF
    IF (N6.NE.0) THEN
      EVFRAC=0.0D+0
      BADFLAG=BADFLAG+1
    ENDIF

    IF (EVFRAC.GE.1.0D+0) EVFRAC=1.0D+0
    IF (EVFRAC.LE.0.0D+0) EVFRAC=0.0D+0
    IF (NEVAP.EQ.0) EVFRAC=0.0D+0
    IF (DHFSUM.LE.0.0D+0) EVFRAC=0.0D+0

    !         EVFRAC=0.0D+0

    IF (SYMBOL(4,J)(:1).EQ.'J          ') THEN
      EVFRAC=1.0D+0-EVFRAC
    ENDIF

    !        V.W. Jul 2006 f=evfrac=0.009 for H2O (Kroes & Anderson 2006) 
    !        IF (SYMBOL(4,J).EQ.'JH2O    ') THEN
    !            EVFRAC=1.0D+0-EVFRAC_H2O
    !         ENDIF

    XJ(J)=XJ(J)*EVFRAC
    !        write(*,*) SYMBOL(1,J)
    !       write(*,*) SYMBOL(2,J)
    !       write(*,*) XJ(J)

    !       stop

    ! ------ Calculate quantum activation energy
    REDMAS=REAL(SMASS(JSP1(J))*SMASS(JSP2(J)))/&
    REAL(SMASS(JSP1(J))+SMASS(JSP2(J)))
    ACT1(J)=2.0D+0 * ACT/HBAR * SQRT(2.0D+0*AMU*REDMAS*BOLTZ*EA(J))
  ENDIF

  ! === ITYPE 16 - C.R. DESORPTION
  IF (ITYPE(J).EQ.16) THEN
    IF (SMASS(JSP1(J)).EQ.0) XJ(J)=0.0D+0
  ENDIF

  ! === ITYPE 99 - ACCRETION
  IF (ITYPE(J).EQ.99) THEN
    ! ------ Save tag of resultant grain surface species
    DO I=1,NSMAX
      IF (SYMBOL(4,J).EQ.SPEC(I)) JSP2(J)=I
    ENDDO
  ENDIF

ENDDO

! === Zero dummy H2 formation rxns, if necc.
!      IF (IDUST.NE.0) THEN
!         XJ(1)=0.0D+0
!         XJ(2)=0.0D+0
!      ENDIF

RETURN
END SUBROUTINE grainrate

! ======================================================================
! ======================================================================
END program gasgrain
