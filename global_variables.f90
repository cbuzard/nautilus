!******************************************************************************
! MODULE: global_variables
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Contains global variables and parameters
!
!******************************************************************************

module global_variables
use numerical_types
use iso_fortran_env

implicit none
integer, parameter :: nptmax = 1 ! number of space points (TODO : allocatable)
integer :: NKMAX
integer :: NSMAX
integer, parameter :: nb_gaseous_species=485
integer, parameter :: nb_surface_species=199
integer, parameter :: NEMAX=13
integer :: nb_species_for_grain !< number of species on the grain surface
integer :: NK2 !< number of reactions on the grain surface
integer :: nb_species_for_gas !< number of species in gas phase
integer :: NK1 !< number of reactions in gas phase
integer, parameter :: NITYPE=100, NOPMAX=1000
integer, parameter :: NL1=105,NL2=52,NL3=43
real(double_precision), parameter :: RXNMIN=1.0D-99
real(double_precision), parameter :: BOLTZ=1.38054D-16,GRAV=6.668D-8,PI=3.1415926535898d0
real(double_precision), parameter :: PLANCK=6.62565D-27,HBAR=1.05459D-27,AMU=1.66043D-24
real(double_precision), parameter :: ECHARGE=1.60219D-19,AVOGADRO=6.0221415D+23
real(double_precision), parameter :: TYEAR=3.1536D+07
real(double_precision), parameter :: AU=1.5D+13

real(double_precision) :: meanw, Omega2

real(double_precision) :: RTOL

character(len=11), allocatable, dimension(:) :: SPEC

character(len=11), allocatable, dimension(:,:) :: SYMBOL
character(len=11) :: YJH,YJH2,YH,YH2,YHE,YHEP,YE,YGRAIN,YCO
integer :: INDCO, INDH2, INDHE
character(len=11), allocatable, dimension(:) :: XS0

real(double_precision), allocatable, dimension(:) :: XJ,A,B,C,R
real(double_precision), allocatable, dimension(:) :: XK
real(double_precision), allocatable, dimension(:) :: XN,XNI,XN0,CTS
real(double_precision), allocatable, dimension(:) :: DXDT,DXDTP,DXDTN,AWT
real(double_precision), dimension(nemax) :: ELEMS,ESUM,EMERR
real(double_precision), dimension(nopmax) :: DENS,AV,TIMS
real(double_precision), allocatable, dimension(:,:) :: XNOP,CTSOP
real(double_precision), allocatable, dimension(:,:) :: DXDTOP
real(double_precision), dimension(6, nopmax) :: EQP
real(double_precision), allocatable, dimension(:) :: RDIF1,RDIF2,EX1,EX2,EA, Tmin,Tmax
real(double_precision), allocatable, dimension(:) :: ACT1
real(double_precision), allocatable, dimension(:) :: TINDIF,TINACC,TINEVA
real(double_precision), allocatable, dimension(:) :: ED,EB,DEB,DHF
real(double_precision), allocatable, dimension(:) :: CHF,CONDSP,RQ1,RQ2
real(double_precision) :: CHARGE,TATOM,CSUM,TSUM,CERR,TERR,DXDTS,DTOGM,GTODN, RAVNH
real(double_precision) :: RD,RHOD,STICK0,STICKP,STICKN,XNMIN
real(double_precision) :: XNT,XNT0,TEMP
real(double_precision) :: TEMP0,DTEMP,DTEMP0,TAU,TAU0,ZETA,ZETA0,XNTI,ZETAX
real(double_precision) :: UVGAS,UVGRA,LAYERS,ACM,SNS,TNS,ACT,TSMAX,CRT,CRFE,EBFAC
real(double_precision) :: TSTART,TFINAL,TIME,ALPHA,BFAC,NF,A1,B1,C1,A2,B2,C2,ARRK
real(double_precision), dimension(nl1) :: N1H2,T1H2
real(double_precision), dimension(nl2) :: N2CO,T2CO
real(double_precision), dimension(nl3) :: N2H2
real(double_precision), dimension(nl3) :: T2H2,AV2,T2AV

integer, dimension(4,nopmax) :: IEQP
integer, dimension(nopmax) :: IORDTM
integer, allocatable, dimension(:) :: INUM, itype, jsp1, jsp2, FORMULA,NUM
integer, allocatable, dimension(:,:) :: ICRTBL,ICROCC
integer, allocatable, dimension(:) :: ICRNUM,ORDSP, icg
integer, allocatable, dimension(:,:) :: IELM
integer, dimension(nemax) :: ISPELM
integer :: ISPE
integer, dimension(0:nitype-1) :: IRXSTA,IRXFIN
integer :: IT,NT,ITFLAG,OTPD,DP,NDP,NS0,IODR,IREFSP, ISORD

integer :: IDENS,ITEMP,IDUST,IGRQM,ICONS,IMODH,IREAD
integer :: IPOUT,IPMON,IPLOT,IPDET,IPRXN,IPORD,ISABS


! Diffusion and 1D variables
real(double_precision), dimension(nptmax) :: zspace, zaspace ! space variables
real(double_precision), allocatable, dimension(:,:) ::  ZXN
real(double_precision) :: zdt ! diffusive timestep
real(double_precision) :: zstepsize ! Spatial resolution
real(double_precision) :: Hsize ! Size of the computing box
real(double_precision) :: diffty ! Turbulent diffusivity
real(double_precision) :: Mcenter ! Central mass
real(double_precision) :: Distr ! Radial distance
real(double_precision) :: Densmax ! Maximum density of the profile
real(double_precision) :: TAUBC ! Av at the edge of the computing box
integer :: idiff ! Diffusivity flag
real(double_precision), dimension(nptmax) :: TEMP1D, DTEMP1D, DEnb_species_for_gasD, TAU1D, ZETAX1D ! 1D physical structure
real(double_precision), dimension(nptmax) :: DIFF1D ! 1D diffusivity profile
real(double_precision), dimension(nptmax) :: ZNCO,ZNH2 ! 1D column density (for the self shielding)
real(double_precision) :: NCO,NH2 ! column density (for the self shielding)
integer :: istep, wstep, wstepr, irateout
integer :: testjac, njac

integer :: iptstore
integer :: ipts

! For FCHEMVW

character(len=11), allocatable, dimension(:) :: SPEC2
integer, allocatable, dimension(:,:) :: REACT

! For IA and JA

integer, allocatable, dimension(:) :: IA
integer, allocatable :: JA(:)

contains 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to retrieve the number of lines of a given file whose
!! filename is passed as an argument. 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_linenumber(filename, nb_lines)

  implicit none
  
  ! Input
  character(len=*), intent(in) :: filename !< [in] the filename of the file we want the number of lines
  
  ! Output
  integer, intent(out) :: nb_lines !< [out] the number of line of the input file
  
  ! Local
  integer :: error
  logical test
  !------------------------------------------------------------------------------
  nb_lines = 0
  
  ! Read in filenames and check for duplicate filenames
  inquire (file=filename, exist=test)
  if (.not.test) then
    write(Error_Unit,'(a,a,a)') 'Error: the file "',trim(filename),'" does not exist.'
  end if
  open(15, file=filename, status='old')
  
  error = 0
  do while(error.eq.0)
    read(15,*,iostat=error)
    nb_lines = nb_lines + 1
  enddo
  close(15)

end subroutine get_linenumber

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to determine array sizes, namely number of reactions, 
!! of species, for gas, grain and in total. 
!! some global size are set (nb_species_for_gas, NK1, nb_species_for_gas, NK2, NSMAX, NKMAX)\n
!! This routine prepare allocation of global dynamical arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_array_sizes()

implicit none

! We get various sizes
call get_linenumber(filename='gas_species.in', nb_lines=nb_species_for_gas)
call get_linenumber(filename='gas_reactions.in', nb_lines=NK1)

call get_linenumber(filename='grain_species.in', nb_lines=nb_species_for_grain)
call get_linenumber(filename='grain_reactions.in', nb_lines=NK2)

NSMAX = nb_species_for_gas + nb_species_for_grain ! The total number of species, sum of species in gas and grain
NKMAX = NK1 + NK2 ! The total number of reactions, sum of species in gas and grain

end subroutine get_array_sizes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2013
!
! DESCRIPTION: 
!> @brief Routine that allocate global arrays once their sizes are set
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine initialize_global_arrays()

implicit none

allocate(spec(nsmax))
allocate(xs0(nsmax))
allocate(xn(nsmax))
allocate(xni(nsmax))
allocate(xn0(nsmax))
allocate(cts(nsmax))
allocate(dxdt(nsmax))
allocate(dxdtp(nsmax))
allocate(dxdtn(nsmax))
allocate(awt(nsmax))
allocate(xnop(nsmax, nopmax))
allocate(ctsop(nsmax, nopmax))
allocate(dxdtop(nsmax, nopmax))
allocate(tindif(nsmax))
allocate(tinacc(nsmax))
allocate(tineva(nsmax))
allocate(ed(nsmax))
allocate(eb(nsmax))
allocate(deb(nsmax))
allocate(dhf(nsmax))
allocate(chf(nsmax))
allocate(condsp(nsmax))
allocate(rq1(nsmax))
allocate(rq2(nsmax))
allocate(icrtbl(nsmax, nkmax))
allocate(icrocc(nsmax, nkmax))
allocate(icrnum(nsmax))
allocate(ordsp(nsmax))
allocate(icg(nsmax))
allocate(zxn(nsmax, nptmax))
allocate(spec2(nsmax+1))
allocate(ia(nsmax+1))
allocate(ielm(nemax,nsmax))

allocate(xj(nkmax))
allocate(a(nkmax))
allocate(b(nkmax))
allocate(c(nkmax))
allocate(r(nkmax))
allocate(xk(nkmax))
allocate(rdif1(nkmax))
allocate(rdif2(nkmax))
allocate(ex1(nkmax))
allocate(ex2(nkmax))
allocate(ea(nkmax))
allocate(tmin(nkmax))
allocate(tmax(nkmax))
allocate(act1(nkmax))
allocate(inum(nkmax))
allocate(itype(nkmax))
allocate(jsp1(nkmax))
allocate(jsp2(nkmax))
allocate(formula(nkmax))
allocate(num(nkmax))
allocate(react(nkmax, 7))
allocate(symbol(7,nkmax))



end subroutine initialize_global_arrays
end module global_variables
