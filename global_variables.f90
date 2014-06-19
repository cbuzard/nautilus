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

! Theses 3 parameters are only intended to easy the transition when one want to add a reagent or a 
!! product (product are fairly easy, reagent are not). But some things needed to be modified in the 
!! code anyway when you want to add either a reagent or a product
integer, parameter :: MAX_REAGENTS = 3 !< The maximum number of reagents for one reaction. 
!! Do not think that changing this parameter alone is sufficient to allow the code to handle it directly !!
integer, parameter :: MAX_PRODUCTS = 5 !< The maximum number of products for one reaction. 
!! Do not think that changing this parameter alone is sufficient to allow the code to handle it directly !!
integer, parameter :: MAX_COMPOUNDS = MAX_REAGENTS + MAX_PRODUCTS !< Total maximum number of compounds for one reaction (reagents + products)
!! Warning: If this number change, get_jacobian(N, T, Y, J, IAN, JAN, PDJ) must be actualised, since each reagent and product has
!! its own variable, a new one must be created for the new column possible. 

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
real(double_precision), parameter :: PLANCK_CONSTANT = 6.62565d-27 !< Planck constant in CGS (erg.s or g cm2 s-1)
real(double_precision), parameter :: SPEED_OF_LIGHT = 2.99792458d10 !< speed of light in CGS (cm/s)
real(double_precision), parameter :: PI = 3.1415926535898d0 !< The number Pi
real(double_precision), parameter :: H_BARRE = 1.054571628d-27 !< Reduced Planck constant h/2*pi in CGS (g cm2 s-1)
real(double_precision), parameter :: AMU = 1.66053892d-24 !< Atomic mass unit in g
real(double_precision), parameter :: ELECTRON_MASS = 0.000548579909d0 !< Electron mass in AMU (close to 1/1836)
real(double_precision), parameter :: AVOGADRO = 6.02214129d23 !< avogadro number : number of atom in 1 mol
real(double_precision), parameter :: YEAR = 3.15576d7 !< one year in seconds
real(double_precision), parameter :: AU = 1.49597871d13 !< Astronomical unit in cm (mean earth-sun distance)

real(double_precision) :: RELATIVE_TOLERANCE !< relative tolerance parameter (scalar) of DLSODES (of the ODEPACK package to solve ODE's)

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
real(double_precision), allocatable, dimension(:) :: INITIAL_ELEMENTAL_ABUNDANCE !< dim(NB_PRIME_ELEMENTS) Store abundances 
!! (relative to H) for all elemental species before running the simulation [number ratio]
real(double_precision), allocatable, dimension(:) :: elemental_mass !< dim(NB_PRIME_ELEMENTS) elemental mass [a.m.u]
character(len=11), allocatable, dimension(:) :: element_name !< dim(NB_PRIME_ELEMENTS) elemental mass.
integer, allocatable, dimension(:) :: PRIME_ELEMENT_IDX ! < dim(NB_PRIME_ELEMENTS) Tell for each prime element its index in the global array of all elements.

! Arrays about species
character(len=11), allocatable, dimension(:) :: species_name !< dim(nb_species)
integer, allocatable, dimension(:,:) :: species_composition !< dim(NB_PRIME_ELEMENTS,nb_species) number of atom of each element composition the given species.
real(double_precision), allocatable, dimension(:) :: abundances !< dim(nb_species) Species abundances (relative to H) [number ratio]
real(double_precision), allocatable, dimension(:) :: SPECIES_MASS !< dim(nb_species) Species mass [a.m.u]
real(double_precision), allocatable, dimension(:) :: THERMAL_HOPING_RATE !< dim(nb_species) Diffusion rate by thermal hopping [s-1]:
!! 1/time required for an adsorbed particle to sweep over a number of sites equivalent to the whole grain surface 
real(double_precision), allocatable, dimension(:) :: CR_HOPING_RATE !< dim(nb_species) Diffusion rate by thermal hopping due to cosmic rays heating [s-1]:
!! 1/time required for an adsorbed particle to sweep over a number of sites equivalent to the whole grain surface
real(double_precision), allocatable, dimension(:) :: ACCRETION_RATES !< dim(nb_species) Accretion rate for a given species onto the grain surface [s-1]
real(double_precision), allocatable, dimension(:) :: EVAPORATION_RATES !< dim(nb_species) evaporation rate for a given species [s-1]
real(double_precision), allocatable, dimension(:) :: EVAPORATION_RATESCR !< dim(nb_species) evaporation rate  due to cosmic rays for a given species [s-1]
real(double_precision), allocatable, dimension(:) :: BINDING_ENERGY !< dim(nb_species) [K] Binding energy of a species to the surface (specific to each species). Parameter read in the file surface_parameters.in
real(double_precision), allocatable, dimension(:) :: DIFFUSION_BARRIER !< dim(nb_species) [K] Potential energy barrier between adjacent surface potential energy wells. 
!! It is usually a fraction of the binding energy. This parameter is used to compute the rate for thermal hopping and the diffusion 
!! by quantum tunneling using Hasegawa’s formalism. Parameter read in the file surface_parameters.in for some species and computed 
!! in the model for the others. (specific of each species)
real(double_precision), allocatable, dimension(:) :: GAP_ENERGY_BANDS !< dim(nb_species) [K] gap between the energy bands corresponding 
!! to the fundamental and first excited state (Ricca et al., 1969). This parameter is used to compute the diffusion of species on 
!! the surface by tunneling using the formalism by Hollenbach & Salpeter (1970) and Watson (1976). It is the parameter dEb in the 
!! equation 20 of Reboussin et al. (2014). Values of this parameters from the literature are also given in the paper.
real(double_precision), allocatable, dimension(:) :: FORMATION_ENTHALPY !< dim(nb_species) Enthalpy formation (read in kcal/mol and then converted into Kelvin/reaction via DHFSUM)
real(double_precision), allocatable, dimension(:) :: VIBRATION_FREQUENCY !< dim(nb_species) Characteristic vibration frequency [s-1] of the adsorbed species  as from a harmonic oscillator hypothesis (Hasegawa & Herbst 1992)
real(double_precision), allocatable, dimension(:) :: ACC_RATES_PREFACTOR !< dim(nb_species) Interrim calculation variable for ACCRETION_RATES
!! that contain cross section, probability and constant needed for thermal motion
real(double_precision), allocatable, dimension(:) :: TUNNELING_RATE_TYPE_1 !< dim(nb_species) Quantum tunneling diffusion rate [s-1] (Watson 1976) (dEB.BOLTZ) / (4.HBAR.nb_sites_per_grain)
real(double_precision), allocatable, dimension(:) :: TUNNELING_RATE_TYPE_2 !< dim(nb_species) Quantum tunneling diffusion rate [s-1] (Hasegawa & Herbst 1992) VIBRATION_FREQUENCY / nb_sites_per_grain.EXP(-2.SITE_SPACING / HBAR.(2.AMU.SMA.BOLTZ.EB)^1/2)
integer, allocatable, dimension(:) :: SPECIES_CHARGE !< dim(nb_species) !< electric charge [in e-] for each species, 0 if neutral, positive or negative if ions.

! Arrays about reactions
character(len=11), allocatable, dimension(:,:) :: REACTION_COMPOUNDS_NAMES !< dim(MAX_COMPOUNDS,nb_reactions)
integer, allocatable, dimension(:,:) :: REACTION_COMPOUNDS_ID !< dim(MAX_COMPOUNDS, nb_reactions) for all reactions, list for reagents (first 3) and products (last 5).
real(double_precision), allocatable, dimension(:) :: branching_ratio !< dim(nb_reactions) Branching ratio of each reaction
real(double_precision), allocatable, dimension(:) :: RATE_A !< dim(nb_reactions) Coefficient used to compute the reaction rate. Formula (and unit) is different in function of the reaction type.
real(double_precision), allocatable, dimension(:) :: RATE_B !< dim(nb_reactions) Coefficient used to compute the reaction rate. Formula (and unit) is different in function of the reaction type.
real(double_precision), allocatable, dimension(:) :: RATE_C !< dim(nb_reactions) Coefficient used to compute the reaction rate. Formula (and unit) is different in function of the reaction type.
real(double_precision), allocatable, dimension(:) :: reaction_rates !< dim(nb_reactions) reaction rate [unit depend on the reaction]
real(double_precision), allocatable, dimension(:) :: ACTIVATION_ENERGY !< dim(nb_reactions) Activation energy of reactions [K]
real(double_precision), allocatable, dimension(:) :: REACTION_TMIN !< dim(nb_reactions) min temperature boundary of each reactions [K]
real(double_precision), allocatable, dimension(:) :: REACTION_TMAX !< dim(nb_reactions) max temperature boundary of each reactions [K]
real(double_precision), allocatable, dimension(:) :: quantum_activation_energy !< dim(nb_reactions) Quantum activation energy
integer, allocatable, dimension(:) :: REACTION_TYPE !< dim(nb_reactions) For each reaction, what is its type (cosmic ray evaporation, etc...)
integer, allocatable, dimension(:) :: RATE_FORMULA !< dim(nb_reactions) The index tracing the formula used for each specific 
!! reaction, defining its reaction rates in function of temperature and abundances.
integer, allocatable, dimension(:) :: REACTION_ID !< dim(nb_reactions) index of the reactions (one of the columns of the concerned file, 
!! declaring a given number for each reaction, like a hashtag.

! Specific variables for first or second reagent of each reactions
integer, allocatable, dimension(:) :: reagent_1_idx !< dim(nb_reactions) Index of the first reagent species involved in the reaction
integer, allocatable, dimension(:) :: reagent_2_idx !< dim(nb_reactions) Index of the second reagent species involved in the reaction
real(double_precision), allocatable, dimension(:) :: THERMAL_DIFFUSION_RATE_1 !< dim(nb_reactions) [s-1] Diffusion rates used to compute the grain reaction rate for reagent 1
real(double_precision), allocatable, dimension(:) :: THERMAL_DIFFUSION_RATE_2 !< dim(nb_reactions) [s-1] Diffusion rates used to compute the grain reaction rate for reagent 2
real(double_precision), allocatable, dimension(:) :: CR_DIFFUSION_RATE_1 !< dim(nb_reactions) [s-1] Diffusion rates used to compute 
!! the grain reaction rate by cosmic rays heating for reagent 1. Used for reaction_type=21
real(double_precision), allocatable, dimension(:) :: CR_DIFFUSION_RATE_2 !< dim(nb_reactions) [s-1] Diffusion rates used to compute 
!! the grain reaction rate by cosmic rays heating for reagent 2. Used for reaction_type=21
real(double_precision), allocatable, dimension(:) :: EVAP_OVER_ACC_RATIO_1 !< dim(nb_reactions) EVAPORATION_RATES/ACCRETION_RATES for reagent 1
real(double_precision), allocatable, dimension(:) :: EVAP_OVER_ACC_RATIO_2 !< dim(nb_reactions) EVAPORATION_RATES/ACCRETION_RATES for reagent 2

real(double_precision) :: initial_dtg_mass_ratio !< [no unit] initial dust to gas mass ratio
real(double_precision) :: GTODN !< Gas to dust number ratio. 1/GTODN is equivalent to the grain abundance [no unit]
real(double_precision) :: AV_NH_ratio !< Extinction over total hydrogen column density [mag/cm-2]
real(double_precision) :: grain_radius !< Grain radius [cm]
real(double_precision) :: GRAIN_DENSITY !< grain density [g/cm^3]
real(double_precision) :: sticking_coeff_neutral  !< sticking coefficient for neutral  species on grain surface [no unit]
real(double_precision) :: sticking_coeff_positive !< sticking coefficient for positive species on grain surface [no unit]
real(double_precision) :: sticking_coeff_negative !< sticking coefficient for negative species on grain surface [no unit]
real(double_precision) :: MINIMUM_INITIAL_ABUNDANCE !< minimum value of the abundance (relative to H) [number ratio]
real(double_precision) :: H_number_density !< [part/cm^3] Total H number density (both H and H2), representing the total gas density
real(double_precision) :: initial_gas_density !< [part/cm^3] initial gas density of the structure
real(double_precision) :: gas_temperature !< current gas temperature [K]
real(double_precision) :: initial_gas_temperature !< initial gas temperature [K], simulation parameter
real(double_precision) :: dust_temperature !< current dust temperature [K]
real(double_precision) :: initial_dust_temperature !< initial dust temperature [K], simulation parameter
real(double_precision) :: visual_extinction !< visual extinction [mag] of the molecular cloud (or other astronomical object). 
!! It's the magnitude attenuation, difference from the absolute magnitude of the object and its apparent magnitude
real(double_precision) :: INITIAL_VISUAL_EXTINCTION !< initial visual extinction [mag] 
real(double_precision) :: CR_IONISATION_RATE !< cosmic ray ionisation rate [s-1]
real(double_precision) :: UV_FLUX !< Scale factor for the UV flux, in unit of the reference flux (1.=nominal)
real(double_precision) :: SITE_SPACING !< site spacing [cm]
real(double_precision) :: SITE_DENSITY !< site density [cm-2]
real(double_precision) :: nb_sites_per_grain !< Number of site per grain (site density * surface of the grain)
real(double_precision) :: ACTIVATION_BARRIER_WIDTH !< grain reaction activation energy barrier width. [cm]
real(double_precision) :: PEAK_GRAIN_TEMPERATURE !< Peak grain temperature when struck by a cosmic ray [K]
real(double_precision) :: PEAK_DURATION !< Peak duration [s] of PEAK_GRAIN_TEMPERATURE
real(double_precision) :: FE_IONISATION_RATE !< (cosmic) Fe-ion--grain encounter [s-1 grain-1] (for 0.1 micron grain) 
!! For cosmic photo desorptions, only Fe-ions are efficient to heat grains. 
real(double_precision) :: DIFF_DESORP_DEFAULT_RATIO !< [no unit] DIFFUSION_BARRIER/BINDING_ENERGY. Ratio used if DIFFUSION_BARRIER is not known
real(double_precision) :: START_TIME !< Start time of the simulation [s]
real(double_precision) :: STOP_TIME !< Stop time of the simulation [s]
real(double_precision) :: current_time !< Global current time of the simulation [s]
real(double_precision) :: VIB_TO_DISSIP_FREQ_RATIO !< [no unit] For the RRK (Rice Ramsperger-Kessel) desorption mechanism. Ratio of the vibration 
!! frequency (proper energy of a species when it is created on a grain) to the dissipation frequency (energy needed by the 
!! molecule to be evaporated from the grain surface). This ratio help to determine if a species evaporate after its formation 
!! on the grain surface. Since the dissipation frequency is usually unknown, this ratio is a free parameter. A common value is 1%.

integer, parameter :: MAX_NUMBER_REACTION_TYPE=100 !< Max number of various reaction type
! The following arrays start at 0 because the index correspond to the reaction type as indexed elsewhere, and there is a type 0 for reactions.
integer, dimension(0:MAX_NUMBER_REACTION_TYPE-1) :: type_id_start !< list of id start for each reaction type given their type number
integer, dimension(0:MAX_NUMBER_REACTION_TYPE-1) :: type_id_stop !< list of id stop for each reaction type given their type number

character(len=80) :: GRAIN_TEMPERATURE_TYPE = 'fixed' !< ('gas', 'computed', 'table', 'fixed') How the grain temperature is computed in the code
procedure(get_grain_temperature_interface), pointer :: get_grain_temperature !< Pointer toward the routine that will calculate grain temperature

abstract interface 
  subroutine get_grain_temperature_interface(time, gas_temperature, grain_temperature)
  import
  
  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  
  ! Outputs
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  end subroutine get_grain_temperature_interface
end interface


integer :: IS_GRAIN_REACTIONS !< Accretion, grain surface reactions
integer :: GRAIN_TUNNELING_DIFFUSION !< How grain tunneling diffusion is handled
integer :: CONSERVATION_TYPE !< 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc
integer :: MODIFY_RATE_FLAG !< Modify rates flag ; 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only
integer :: IS_ABSORPTION !< H2 AND CO SELF-SHIELDING

! About IS_STRUCTURE_EVOLUTION, describing the evolution of the physical structure properties with time
integer :: IS_STRUCTURE_EVOLUTION = 0 !< if 1, physical structure properties evolve with time. They come from structure_evolution.dat file, containing
!! {time [Myr], Av [mag], number density [part/cm3], gas temperature [K] and possibly grain temperature [K]} for the structure
procedure(get_structure_properties_interface), pointer :: get_structure_properties

abstract interface 
  subroutine get_structure_properties_interface(time, Av, density, gas_temperature, grain_temperature)
    import 
    
    implicit none
  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]
  
  end subroutine get_structure_properties_interface
end interface

! About the optimization so that, for each species, we only check the reactions we know the species is involved.
integer :: max_reactions_same_species !< Maximum number of reactions in which any species will be involved. Used to set the array 'relevant_reactions'
integer, dimension(:,:), allocatable :: relevant_reactions !< dim(max_reactions_same_species, nb_species) For each species, store the list of reactions involving the species. 
!! When 0 are encountered for a given species, this means that no more reactions involve it. dimensions : (max_reactions_same_species, nb_species)
integer, dimension(:), allocatable :: nb_reactions_using_species !< dim(nb_species) For each species, the number of reactions in which it is used

! For LSODES
integer :: lrw !< declared length of RWORK (in user dimension).
integer :: liw !< declared length of IWORK (in user dimension).
integer, dimension(:), allocatable :: IWORK !< dim(liw) integer work array of length at least 30.
real(double_precision), dimension(:), allocatable :: RWORK !< dim(lrw) real work array of length at least:
!!\n             20 + 16*NEQ            for MF = 10,
!!\n             20 + (2 + 1./LENRAT)*NNZ + (11 + 9./LENRAT)*NEQ
!!\n                                    for MF = 121 or 222,
!!\n          where:
!!\n          NNZ    = the number of nonzero elements in the sparse
!!\n                   Jacobian (if this is unknown, use an estimate), and
!!\n          LENRAT = the real to integer wordlength ratio (usually 1 in
!!\n                   single precision and 2 in double precision).
!!\n          In any case, the required size of RWORK cannot generally
!!\n          be predicted in advance if MF = 121 or 222, and the value
!!\n          above is a rough estimate of a crude lower bound.  Some
!!\n          experimentation with this size may be necessary.
!!\n          (When known, the correct required length is an optional
!!\n          output, available in IWORK(17).)
integer :: nb_nonzeros_values !< number of non-zeros values in the jacobian. This is usefull for ODEPACK, to increase speed

! Diffusion and 1D variables
real(double_precision) :: X_IONISATION_RATE !< Ionisation rate due to X-rays [s-1]
real(double_precision) :: NCO ! column density [cm-2] (for the self shielding)
real(double_precision) :: NH2 ! column density [cm-2] (for the self shielding)

logical :: first_step_done = .false. !< do we have currently done the first step of the integrator ?
integer :: NB_OUTPUTS !< Total number of outputs in the simulation
character(len=80) :: OUTPUT_TYPE !< Type of output times sampling. linear, log, table\n
!! linear: Output times are linearly spaced\n 
!! log   : Outputs times are log-spaced 
!! table : Outputs times are read from time_evolution.dat'

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
allocate(SPECIES_MASS(nb_species))
allocate(THERMAL_HOPING_RATE(nb_species))
allocate(CR_HOPING_RATE(nb_species))
allocate(ACCRETION_RATES(nb_species))
allocate(EVAPORATION_RATES(nb_species))
allocate(EVAPORATION_RATEScr(nb_species))
allocate(BINDING_ENERGY(nb_species))
allocate(DIFFUSION_BARRIER(nb_species))
allocate(GAP_ENERGY_BANDS(nb_species))
allocate(FORMATION_ENTHALPY(nb_species))
allocate(VIBRATION_FREQUENCY(nb_species))
allocate(ACC_RATES_PREFACTOR(nb_species))
allocate(TUNNELING_RATE_TYPE_1(nb_species))
allocate(TUNNELING_RATE_TYPE_2(nb_species))
allocate(SPECIES_CHARGE(nb_species))
allocate(nb_reactions_using_species(nb_species))

allocate(branching_ratio(nb_reactions))
allocate(RATE_A(nb_reactions))
allocate(RATE_B(nb_reactions))
allocate(RATE_C(nb_reactions))
allocate(reaction_rates(nb_reactions))
allocate(THERMAL_DIFFUSION_RATE_1(nb_reactions))
allocate(CR_DIFFUSION_RATE_1(nb_reactions))
allocate(THERMAL_DIFFUSION_RATE_2(nb_reactions))
allocate(CR_DIFFUSION_RATE_2(nb_reactions))
allocate(EVAP_OVER_ACC_RATIO_1(nb_reactions))
allocate(EVAP_OVER_ACC_RATIO_2(nb_reactions))
allocate(ACTIVATION_ENERGY(nb_reactions))
allocate(REACTION_TMIN(nb_reactions))
allocate(REACTION_TMAX(nb_reactions))
allocate(quantum_activation_energy(nb_reactions))
allocate(REACTION_TYPE(nb_reactions))
allocate(reagent_1_idx(nb_reactions))
allocate(reagent_2_idx(nb_reactions))
allocate(RATE_FORMULA(nb_reactions))
allocate(REACTION_ID(nb_reactions))
allocate(REACTION_COMPOUNDS_ID(MAX_COMPOUNDS,nb_reactions))
allocate(REACTION_COMPOUNDS_NAMES(MAX_COMPOUNDS,nb_reactions))

allocate(INITIAL_ELEMENTAL_ABUNDANCE(NB_PRIME_ELEMENTS))
allocate(PRIME_ELEMENT_IDX(NB_PRIME_ELEMENTS))
allocate(species_composition(NB_PRIME_ELEMENTS,nb_species))

species_name(1:nb_species) = ''
abundances(1:nb_species) = 0.d0
SPECIES_MASS(1:nb_species) = 0.d0
THERMAL_HOPING_RATE(1:nb_species) = 0.d0
ACCRETION_RATES(1:nb_species) = 0.d0
EVAPORATION_RATES(1:nb_species) = 0.d0
BINDING_ENERGY(1:nb_species) = 0.d0
DIFFUSION_BARRIER(1:nb_species) = 0.d0
GAP_ENERGY_BANDS(1:nb_species) = 0.d0
FORMATION_ENTHALPY(1:nb_species) = 0.d0
VIBRATION_FREQUENCY(1:nb_species) = 0.d0
ACC_RATES_PREFACTOR(1:nb_species) = 0.d0
TUNNELING_RATE_TYPE_1(1:nb_species) = 0.d0
TUNNELING_RATE_TYPE_2(1:nb_species) = 0.d0
SPECIES_CHARGE(1:nb_species) = 0
nb_reactions_using_species(1:nb_species) = 0

branching_ratio(1:nb_reactions) = 0.d0
RATE_A(1:nb_reactions) = 0.d0
RATE_B(1:nb_reactions) = 0.d0
RATE_C(1:nb_reactions) = 0.d0
reaction_rates(1:nb_reactions) = 0.d0
THERMAL_DIFFUSION_RATE_1(1:nb_reactions) = 0.d0
THERMAL_DIFFUSION_RATE_2(1:nb_reactions) = 0.d0
EVAP_OVER_ACC_RATIO_1(1:nb_reactions) = 0.d0
EVAP_OVER_ACC_RATIO_2(1:nb_reactions) = 0.d0
ACTIVATION_ENERGY(1:nb_reactions) = 0.d0
REACTION_TMIN(1:nb_reactions) = 0.d0
REACTION_TMAX(1:nb_reactions) = 0.d0
quantum_activation_energy(1:nb_reactions) = 0.d0
REACTION_TYPE(1:nb_reactions) = 0
reagent_1_idx(1:nb_reactions) = 0
reagent_2_idx(1:nb_reactions) = 0
RATE_FORMULA(1:nb_reactions) = 0
REACTION_ID(1:nb_reactions) = 0
REACTION_COMPOUNDS_ID(1:MAX_COMPOUNDS,1:nb_reactions) = 0
REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS,1:nb_reactions) = ''
INITIAL_ELEMENTAL_ABUNDANCE(1:NB_PRIME_ELEMENTS) = 0.d0
PRIME_ELEMENT_IDX(1:NB_PRIME_ELEMENTS) = 0
species_composition(1:NB_PRIME_ELEMENTS,nb_species) = 0

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
