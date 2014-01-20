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
integer :: nb_reactions !< total number of reactions
integer :: nb_species !< total number of species
integer :: nb_gaseous_species !< number of species that are gaseous
integer :: nb_surface_species !< number of species that are on the surface of grains
integer, parameter :: NEMAX=13
integer :: nb_species_for_grain !< number of species involved in grain surface reactions
integer :: nb_surface_reactions !< number of reactions on the grain surface
integer :: nb_species_for_gas !< number of species involved in gas phase reactions
integer :: nb_gas_phase_reactions !< number of reactions in gas phase
integer, parameter :: NITYPE=100, NOPMAX=1000
integer, parameter :: NL1=105,NL2=52,NL3=43
real(double_precision), parameter :: RXNMIN=1.0D-99
real(double_precision), parameter :: BOLTZ = 1.3806488d-16 !< Boltzmann constant in CGS (cm^2 g s^⁻2 K-1)
real(double_precision), parameter :: GRAV = 6.67384d-8 !< Gravitationnal constant in CGS (cm^3 g-1 s-2)
real(double_precision), parameter :: PI = 3.1415926535898d0
real(double_precision), parameter :: PLANCK = 6.62606957d-27
real(double_precision), parameter :: HBAR = 1.05459D-27
real(double_precision), parameter :: AMU = 1.66043D-24
real(double_precision), parameter :: ECHARGE = 1.60219D-19
real(double_precision), parameter :: AVOGADRO = 6.0221415D+23
real(double_precision), parameter :: TYEAR = 3.1536D+07
real(double_precision), parameter :: AU = 1.5D+13

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
real(double_precision) :: CHARGE
real(double_precision) :: TATOM
real(double_precision) :: CSUM
real(double_precision) :: TSUM
real(double_precision) :: CERR
real(double_precision) :: TERR
real(double_precision) :: DXDTS
real(double_precision) :: initial_dtg_mass_ratio
real(double_precision) :: GTODN
real(double_precision) :: AV_NH_ratio
real(double_precision) :: grain_radius
real(double_precision) :: RHOD
real(double_precision) :: sticking_coeff_neutral
real(double_precision) :: sticking_coeff_positive
real(double_precision) :: sticking_coeff_negative
real(double_precision) :: XNMIN
real(double_precision) :: XNT
real(double_precision) :: initial_gas_density
real(double_precision) :: TEMP
real(double_precision) :: initial_gas_temperature
real(double_precision) :: DTEMP
real(double_precision) :: initial_dust_temperature
real(double_precision) :: TAU
real(double_precision) :: TAU0
real(double_precision) :: ZETA
real(double_precision) :: ZETA0
real(double_precision) :: XNTI
real(double_precision) :: ZETAX
real(double_precision) :: UVGAS
real(double_precision) :: UVGRA
real(double_precision) :: LAYERS
real(double_precision) :: ACM
real(double_precision) :: SNS
real(double_precision) :: nb_sites_per_grain
real(double_precision) :: ACT
real(double_precision) :: TSMAX
real(double_precision) :: CRT
real(double_precision) :: CRFE
real(double_precision) :: EBFAC
real(double_precision) :: TSTART
real(double_precision) :: TFINAL
real(double_precision) :: TIME
real(double_precision) :: ALPHA
real(double_precision) :: BFAC
real(double_precision) :: NF
real(double_precision) :: A1
real(double_precision) :: B1
real(double_precision) :: C1
real(double_precision) :: A2
real(double_precision) :: B2
real(double_precision) :: C2
real(double_precision) :: ARRK
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
integer :: timestep,NT,ITFLAG,OTPD,DP,NDP,NS0,IODR,IREFSP, ISORD

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
real(double_precision) :: Denb_species ! Maximum density of the profile
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
!! some global size are set (nb_species_for_gas, nb_gas_phase_reactions, nb_species_for_gas, nb_surface_reactions, nb_species, nb_reactions)\n
!! This routine prepare allocation of global dynamical arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_array_sizes()

implicit none

! We get various sizes
call get_linenumber(filename='gas_species.in', nb_lines=nb_species_for_gas)
call get_linenumber(filename='gas_reactions.in', nb_lines=nb_gas_phase_reactions)

call get_linenumber(filename='grain_species.in', nb_lines=nb_species_for_grain)
call get_linenumber(filename='grain_reactions.in', nb_lines=nb_surface_reactions)

nb_species = nb_species_for_gas + nb_species_for_grain ! The total number of species, sum of species in gas and grain
nb_reactions = nb_gas_phase_reactions + nb_surface_reactions ! The total number of reactions, sum of species in gas and grain

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

allocate(spec(nb_species))
allocate(xs0(nb_species))
allocate(xn(nb_species))
allocate(xni(nb_species))
allocate(xn0(nb_species))
allocate(cts(nb_species))
allocate(dxdt(nb_species))
allocate(dxdtp(nb_species))
allocate(dxdtn(nb_species))
allocate(awt(nb_species))
allocate(xnop(nb_species, nopmax))
allocate(ctsop(nb_species, nopmax))
allocate(dxdtop(nb_species, nopmax))
allocate(tindif(nb_species))
allocate(tinacc(nb_species))
allocate(tineva(nb_species))
allocate(ed(nb_species))
allocate(eb(nb_species))
allocate(deb(nb_species))
allocate(dhf(nb_species))
allocate(chf(nb_species))
allocate(condsp(nb_species))
allocate(rq1(nb_species))
allocate(rq2(nb_species))
allocate(icrtbl(nb_species, nb_reactions))
allocate(icrocc(nb_species, nb_reactions))
allocate(icrnum(nb_species))
allocate(ordsp(nb_species))
allocate(icg(nb_species))
allocate(zxn(nb_species, nptmax))
allocate(spec2(nb_species+1))
allocate(ia(nb_species+1))
allocate(ielm(nemax,nb_species))

allocate(xj(nb_reactions))
allocate(a(nb_reactions))
allocate(b(nb_reactions))
allocate(c(nb_reactions))
allocate(r(nb_reactions))
allocate(xk(nb_reactions))
allocate(rdif1(nb_reactions))
allocate(rdif2(nb_reactions))
allocate(ex1(nb_reactions))
allocate(ex2(nb_reactions))
allocate(ea(nb_reactions))
allocate(tmin(nb_reactions))
allocate(tmax(nb_reactions))
allocate(act1(nb_reactions))
allocate(inum(nb_reactions))
allocate(itype(nb_reactions))
allocate(jsp1(nb_reactions))
allocate(jsp2(nb_reactions))
allocate(formula(nb_reactions))
allocate(num(nb_reactions))
allocate(react(nb_reactions, 7))
allocate(symbol(7,nb_reactions))



end subroutine initialize_global_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2013
!
! DESCRIPTION: 
!> @brief From the label of each species, determine if this is a gas
!! phase or a surface species. Set the values of nb_gaseous_species and
!! nb_surface_species global parameters that count the total number of species
!! in each category.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_gas_surface_species()

implicit none

! Locals
integer :: i

! We retrieve the total number of gas and surface species (not the ones that are involved in reactions, but the actual position of the species)
nb_gaseous_species = 0
nb_surface_species = 0
do i=1,nb_species
  if (SPEC(i)(1:1).eq.'J') then
    nb_surface_species = nb_surface_species + 1
  else
    nb_gaseous_species = nb_gaseous_species + 1
  endif
end do

end subroutine get_gas_surface_species

end module global_variables
