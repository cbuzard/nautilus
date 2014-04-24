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
use utilities

implicit none
integer :: nb_reactions !< total number of reactions
integer :: nb_species !< total number of species
integer :: nb_gaseous_species !< number of species that are gaseous
integer :: nb_surface_species !< number of species that are on the surface of grains
integer :: NB_PRIME_ELEMENTS !< Number of prime element that compose molecules, such as H, He, C and so on.
integer :: nb_species_for_grain !< number of species involved in grain surface reactions (can be gas or grain phase elements)
integer :: nb_surface_reactions !< number of reactions on the grain surface
integer :: nb_species_for_gas !< number of species involved in gas phase reactions (can be gas or grain phase elements)
integer :: nb_gas_phase_reactions !< number of reactions in gas phase
real(double_precision), parameter :: MINIMUM_RATE_COEFFICIENT=1.0D-99 !< Minimum rate coefficient (Below, coefficients are forced to 0)
real(double_precision), parameter :: K_B = 1.3806488d-16 !< Boltzmann constant in CGS (cm^2 g s^⁻2 K-1)
real(double_precision), parameter :: PI = 3.1415926535898d0 !< The number Pi
real(double_precision), parameter :: H_BARRE = 1.054571628d-27 !< Reduced Planck constant h/2*pi in CGS (g cm2 s-1)
real(double_precision), parameter :: AMU = 1.66053892d-24 !< Atomic mass unit in g
real(double_precision), parameter :: AVOGADRO = 6.02214129d23 !< avogadro number : number of atom in 1 mol
real(double_precision), parameter :: YEAR = 3.15576d7 !< one year in seconds
real(double_precision), parameter :: AU = 1.49597871d13 !< Astronomical unit in cm (mean earth-sun distance)

real(double_precision) :: RELATIVE_TOLERANCE

! Name of key species
character(len=11) :: YH     = 'H          ' !< Gas phase Hydrogen
character(len=11) :: YJH    = 'JH         ' !< Hydrogen on grains
character(len=11) :: YH2    = 'H2         ' !< Gas phase Dihydrogen
character(len=11) :: YJH2   = 'JH2        ' !< Dihydrogen on grains
character(len=11) :: YHE    = 'He         ' !< Gas phase Helium
character(len=11) :: YHEP   = 'He+        ' !< Gas phase Helium+
character(len=11) :: YE     = 'e-         ' !< Gas phase electrons
character(len=11) :: YGRAIN = 'GRAIN0     ' !< Grain
character(len=11) :: YCO    = 'CO         ' !< Gas phase CO

integer :: INDCO !< Index corresponding to CO in nb_species length arrays
integer :: INDH2 !< Index corresponding to H2 in nb_species length arrays
integer :: INDHE !< Index corresponding to He in nb_species length arrays
integer :: INDEL !< Index corresponding to e- in nb_species length arrays


! Arrays about prime elements
real(double_precision), allocatable, dimension(:) :: INITIAL_ELEMENTAL_ABUNDANCE !< dim(NB_PRIME_ELEMENTS) Store abundances for all 
!! elemental species before running the simulation.
real(double_precision), allocatable, dimension(:) :: elemental_mass !< dim(NB_PRIME_ELEMENTS) elemental mass.
character(len=11), allocatable, dimension(:) :: element_name !< dim(NB_PRIME_ELEMENTS) elemental mass.
integer, allocatable, dimension(:) :: PRIME_ELEMENT_IDX ! < dim(NB_PRIME_ELEMENTS) Tell for each prime element its index in the global array of all elements.

! Arrays about species
character(len=11), allocatable, dimension(:) :: species_name !< dim(nb_species)
integer, allocatable, dimension(:,:) :: species_composition !< dim(NB_PRIME_ELEMENTS,nb_species) number of atom of each element composition the given species.
real(double_precision), allocatable, dimension(:) :: abundances !< dim(nb_species) Species abundances
real(double_precision), allocatable, dimension(:) :: AWT !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: TINDIF !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: TINACC !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: TINEVA !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: ED !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: EB !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: DEB !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: DHF !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: CHF !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: CONDSP !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: RQ1 !< dim(nb_species)
real(double_precision), allocatable, dimension(:) :: RQ2 !< dim(nb_species)
integer, allocatable, dimension(:) :: SPECIES_CHARGE !< dim(nb_species) !< electric charge [in e-] for each species, 0 if neutral, positive or negative if ions.

! Arrays about reactions
character(len=11), allocatable, dimension(:,:) :: REACTION_SUBSTANCES_NAMES !< dim(7,nb_reactions)
integer, allocatable, dimension(:,:) :: REACTION_SUBSTANCES_ID !< dim(7, nb_reactions) for all reactions, list for reactants (first 3) and products (last 4).
real(double_precision), allocatable, dimension(:) :: XJ !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: A !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: B !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: C !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: XK !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: RDIF1 !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: RDIF2 !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: EX1 !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: EX2 !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: EA !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: Tmin !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: Tmax !< dim(nb_reactions)
real(double_precision), allocatable, dimension(:) :: ACT1 !< dim(nb_reactions)
integer, allocatable, dimension(:) :: REACTION_TYPE !< dim(nb_reactions) For each reaction, what is its type (cosmic ray evaporation, etc...)
integer, allocatable, dimension(:) :: jsp1 !< dim(nb_reactions)
integer, allocatable, dimension(:) :: jsp2 !< dim(nb_reactions)
integer, allocatable, dimension(:) :: FORMULA !< dim(nb_reactions)
integer, allocatable, dimension(:) :: REACTION_ID !< dim(nb_reactions) index of the reactions (one of the columns of the concerned file, 
!! declaring a given number for each reaction, like a hashtag.

real(double_precision) :: initial_dtg_mass_ratio
real(double_precision) :: GTODN
real(double_precision) :: AV_NH_ratio
real(double_precision) :: grain_radius
real(double_precision) :: GRAIN_DENSITY
real(double_precision) :: sticking_coeff_neutral
real(double_precision) :: sticking_coeff_positive
real(double_precision) :: sticking_coeff_negative
real(double_precision) :: MINIMUM_INITIAL_ABUNDANCE
real(double_precision) :: H_number_density !< Total H number density (both H and H2)
real(double_precision) :: initial_gas_density
real(double_precision) :: gas_temperature !< current gas temperature [K]
real(double_precision) :: initial_gas_temperature !< initial gas temperature [K], simulation parameter
real(double_precision) :: dust_temperature !< current dust temperature [K]
real(double_precision) :: initial_dust_temperature !< initial dust temperature [K], simulation parameter
real(double_precision) :: visual_extinction !< visual extinction of the molecular cloud (or other astronomical object). 
!! It's the magnitude attenuation, difference from the absolute magnitude of the object and its apparent magnitude
real(double_precision) :: INITIAL_VISUAL_EXTINCTION
real(double_precision) :: CR_IONISATION_RATE
real(double_precision) :: UV_FLUX
real(double_precision) :: SITE_SPACING
real(double_precision) :: SITE_DENSITY
real(double_precision) :: nb_sites_per_grain
real(double_precision) :: ACT
real(double_precision) :: TSMAX
real(double_precision) :: CRT
real(double_precision) :: CRFE
real(double_precision) :: EBFAC
real(double_precision) :: START_TIME !< Start time of the simulation [s]
real(double_precision) :: STOP_TIME !< Stop time of the simulation [s]
real(double_precision) :: current_time !< Global current time of the simulation [s]
real(double_precision) :: ARRK

integer, parameter :: MAX_NUMBER_REACTION_TYPE=100 !< Max number of various reaction type
! The following arrays start at 0 because the index correspond to the reaction type as indexed elsewhere, and there is a type 0 for reactions.
integer, dimension(0:MAX_NUMBER_REACTION_TYPE-1) :: type_id_start !< list of id start for each reaction type given their type number
integer, dimension(0:MAX_NUMBER_REACTION_TYPE-1) :: type_id_stop !< list of id stop for each reaction type given their type number

integer :: IS_GRAIN_REACTIONS
integer :: GRAIN_TUNNELING_DIFFUSION
integer :: CONSERVATION_TYPE
integer :: IMODH
integer :: IS_ABSORPTION

! About the optimization so that, for each species, we only check the reactions we know the species is involved.
integer :: max_reactions_same_species !< Maximum number of reactions in which any species will be involved. Used to set the array 'relevant_reactions'
integer, dimension(:,:), allocatable :: relevant_reactions !< dim(max_reactions_same_species, nb_species) For each species, store the list of reactions involving the species. 
!! When 0 are encountered for a given species, this means that no more reactions involve it. dimensions : (max_reactions_same_species, nb_species)
integer, dimension(:), allocatable :: nb_reactions_using_species !< dim(nb_species) For each species, the number of reactions in which it is used

! For LSODES
integer :: lrw, liw
integer, dimension(:), allocatable :: IWORK !< dim(liw)
real(double_precision), dimension(:), allocatable :: RWORK !< dim(lrw)
integer :: nb_nonzeros_values !< number of non-zeros values in the jacobian. This is usefull for ODEPACK, to increase speed

! Diffusion and 1D variables
real(double_precision) :: X_IONISATION_RATE
real(double_precision) :: NCO ! column density [cm-2] (for the self shielding)
real(double_precision) :: NH2 ! column density [cm-2] (for the self shielding)

integer :: timestep
integer :: OUTPUT_PER_DECADE
integer :: wstep
integer :: wstepr

contains 

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

! We get the number of reactions and species
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
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine that allocate global arrays once their sizes are set
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine initialize_global_arrays()

implicit none

allocate(species_name(nb_species))
allocate(abundances(nb_species))
allocate(awt(nb_species))
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
allocate(SPECIES_CHARGE(nb_species))
allocate(nb_reactions_using_species(nb_species))

allocate(xj(nb_reactions))
allocate(a(nb_reactions))
allocate(b(nb_reactions))
allocate(c(nb_reactions))
allocate(xk(nb_reactions))
allocate(rdif1(nb_reactions))
allocate(rdif2(nb_reactions))
allocate(ex1(nb_reactions))
allocate(ex2(nb_reactions))
allocate(ea(nb_reactions))
allocate(tmin(nb_reactions))
allocate(tmax(nb_reactions))
allocate(act1(nb_reactions))
allocate(REACTION_TYPE(nb_reactions))
allocate(jsp1(nb_reactions))
allocate(jsp2(nb_reactions))
allocate(formula(nb_reactions))
allocate(REACTION_ID(nb_reactions))
allocate(REACTION_SUBSTANCES_ID(7,nb_reactions))
allocate(REACTION_SUBSTANCES_NAMES(7,nb_reactions))

allocate(INITIAL_ELEMENTAL_ABUNDANCE(NB_PRIME_ELEMENTS))
allocate(PRIME_ELEMENT_IDX(NB_PRIME_ELEMENTS))
allocate(species_composition(NB_PRIME_ELEMENTS,nb_species))



end subroutine initialize_global_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
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
  if (species_name(i)(1:1).eq.'J') then
    nb_surface_species = nb_surface_species + 1
  else
    nb_gaseous_species = nb_gaseous_species + 1
  endif
end do

end subroutine get_gas_surface_species

end module global_variables
