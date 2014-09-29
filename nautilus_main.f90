!******************************************************************************
! MODULE: nautilus_main
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the main routines of nautilus. Usefull to 
!! use them in several programs such as nautilus and unitary_tests for instance. \n\n
!!
!
!******************************************************************************

module nautilus_main

use iso_fortran_env
use global_variables
use structure
use input_output
use ode_solver
use dust_temperature_module

implicit none

contains

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
    if (I.NE.INDEL) CHASUM=CHASUM+SPECIES_CHARGE(I)*temp_abundances(I)
  enddo
  if (CHASUM.LE.0.d0) CHASUM=MINIMUM_INITIAL_ABUNDANCE
  temp_abundances(INDEL)=CHASUM

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

! VW fev 2012 add a test for the helium and H2 abundance in the gas phase
! warn for excessive depletion

  do k=1,NB_PRIME_ELEMENTS
    if (abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-elemental_abundance(K))/INITIAL_ELEMENTAL_ABUNDANCE(K).ge.0.01d0) then 
      write(Error_unit,*) 'Caution: Element ', trim(species_name(PRIME_ELEMENT_IDX(K))), 'is not conserved'
      write(Error_unit,*) 'Relative difference: ', abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-elemental_abundance(K)) / &
                           INITIAL_ELEMENTAL_ABUNDANCE(K)
    endif
    if (species_name(PRIME_ELEMENT_IDX(K)).eq.YH) then
      if (abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-(temp_abundances(INDH2)*2.D0+temp_abundances(INDH))/ &
            INITIAL_ELEMENTAL_ABUNDANCE(K)).ge.0.01d0) then
        write(Error_unit,'(a,a,a,es10.2e2,a,a,a,a,a,es10.2e2)') 'H is too depleted on the grains: initial &
                            Y(',trim(species_name(K)),') =',INITIAL_ELEMENTAL_ABUNDANCE(K), ' ; Y(',&
                            trim(species_name(INDH2)),')*2 + Y(',trim(species_name(INDH)),') =',&
                                temp_abundances(INDH2)*2.d0+temp_abundances(INDH)
      endif
    endif
    if (species_name(PRIME_ELEMENT_IDX(K)).eq.YHE) then
      if (abs(INITIAL_ELEMENTAL_ABUNDANCE(K)-temp_abundances(INDHE))/INITIAL_ELEMENTAL_ABUNDANCE(K).ge.0.01d0) then
        write(Error_unit,'(2(a,a,a,es10.2e2))') 'He is too depleted on the grains: initial Y(',trim(species_name(K)),') =',&
                             INITIAL_ELEMENTAL_ABUNDANCE(K), ' ; Y(',trim(species_name(INDHE)),')*2 =',temp_abundances(INDHE)
      endif
    endif       
  enddo


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
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief To check coherence of the parameter files, reactions and so on.
!! Writing information in 'info.out', appending to the file created by 'write_general_infos'
! There is a switch in the parameters.in file
! The following tests are done :
!   - check that all species are both produced and destroyed
!   - Check the balance of the reactions for elements and charges
!   - check for reactions with alpha = 0
!   - Check that Tmin < Tmax for all reactions
!   - check for reactions with the same ID (must have the same reactants and products and 
!       have complementary T range)
!   - check that each gas species as a grain equivalent
!   - check that each gas neutral has a depletion reaction
!   - check that each grain species has a desorption reaction for each type
!   - check that each surface species has a binding energy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine preliminary_tests()

implicit none

! Parameters
character(len=80), parameter :: information_file = 'info.out' !< File in which all warnings will be displayed

! Locals
integer :: i, j, species, reaction, compound, element !< For loops

! To check production and destruction reactions for each species
integer :: nb_production_reactions !< Number of reactions that produce a given species.
integer :: nb_destruction_reactions !< Number of reactions that destruct a given species.
integer :: production_reaction_id !< last reaction id for production of a given species
integer :: destruction_reaction_id !< last reaction id for destruction of a given species
character(len=11) :: tmp_name !< Name of one species. Only a local temporary variable used to construct the reaction string line
character(len=80) :: reaction_line !< will store the line of a given reaction, that will be constructed by concatenation

! To check if reactions are equilibrated
integer :: left_sum, right_sum !< The sum of one prime element for each side of a reaction

! To check reactions with identical ID's
integer reaction2 !< Index for the reaction second loop
integer :: ref_id !< reference ID to find twin reactions
real(double_precision) :: range1_max, range2_min !< To test overlap between temperature ranges of reactions

! To check reactions with alpha = 0
logical :: alpha_equal_zero !< If any (at least one) reaction has alpha = 0, then display a warning.

! To check gas and grain species
logical :: no_grain_equivalent !< true if a gas species has no grain equivalent.

! To check that some species are at least reactant of one reaction of a given type
logical :: no_itype !< False if one species is the reactant of at least one reaction of a given type

! To check reaction_ID index start and stop
integer :: id_start, id_stop

!-------------------------------------------------



if (IS_TEST.eq.1) then
  open(12, file=information_file, position='append')

  write(12, '(a)') '!----------------------------'
  write(12, '(a)') '!     Preliminary tests      |'
  write(12, '(a)') '!----------------------------'
  close(12)

  ! Test if all species whose index is 'grain' type, are effectively "on-grain" species
  do i = nb_gaseous_species+1,nb_species
    if(species_name(i)(:1).NE.'J          ') then
      write(Error_unit, *) 'Error: Species number ',i,' (',trim(species_name(i)),' is expected to be a grain species but is not'
      call exit(13)
    endif
  enddo

  ! Check if all species have production AND destruction reactions
  ! Display a warning if there is only one production or destruction reactions
  open(12, file=information_file, position='append')

  write(12, *) ' ### CHECK ### if all species have production AND destruction reactions'
  do species=1,nb_species
    ! This way of doing things is not the fastest, but it is the only convenient way to calculate the number of 
    ! reactions involving one species without counting twice the same reaction when for instance, we have H + H -> H2
    
    nb_production_reactions = 0
    nb_destruction_reactions = 0
    
    do reaction=1,nb_reactions
      if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS, reaction).eq.species)) then
        nb_destruction_reactions = nb_destruction_reactions + 1
        destruction_reaction_id = reaction
      endif
      
      if (any(REACTION_COMPOUNDS_ID(MAX_REACTANTS+1:MAX_COMPOUNDS, reaction).eq.species)) then
        nb_production_reactions = nb_production_reactions + 1
        production_reaction_id = reaction
      endif
    enddo
    
    if (nb_destruction_reactions.eq.0) then
      write(Error_unit,*) 'Error: ',trim(species_name(species)), ' have no destruction reaction.'
      call exit(14)
    else if (nb_destruction_reactions.eq.1) then
      write(12,'(a,a,a,i0,a)') 'Warning: ',trim(species_name(species)), ' have only one destruction reaction (number ',&
                          &REACTION_ID(destruction_reaction_id),'):'
      
      ! We construct the reaction string display
      reaction_line = trim(REACTION_COMPOUNDS_NAMES(1,destruction_reaction_id))
      do compound=2,MAX_REACTANTS
        tmp_name = REACTION_COMPOUNDS_NAMES(compound,destruction_reaction_id)
        if (tmp_name.ne.'') then
          reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
        endif
      enddo
      reaction_line = trim(reaction_line)//" -> "//trim(REACTION_COMPOUNDS_NAMES(MAX_REACTANTS+1,destruction_reaction_id))
      do compound=MAX_REACTANTS+2,MAX_COMPOUNDS
        tmp_name = REACTION_COMPOUNDS_NAMES(compound,destruction_reaction_id)
        if (tmp_name.ne.'') then
          reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
        endif
      enddo
      write(12,*) trim(reaction_line)
    endif
    
    if (nb_production_reactions.eq.0) then
      write(Error_unit,*) 'Error: ',trim(species_name(species)), ' have no production reaction.'
      call exit(14)
    else if (nb_production_reactions.eq.1) then
      write(12,'(a,a,a,i0,a)') 'Warning: ',trim(species_name(species)), ' have only one production reaction (number ',&
                           &REACTION_ID(production_reaction_id),'):'
      
      ! We construct the reaction string display
      reaction_line = trim(REACTION_COMPOUNDS_NAMES(1,production_reaction_id))
      do compound=2,MAX_REACTANTS
        tmp_name = REACTION_COMPOUNDS_NAMES(compound,production_reaction_id)
        if (tmp_name.ne.'') then
          reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
        endif
      enddo
      reaction_line = trim(reaction_line)//" -> "//trim(REACTION_COMPOUNDS_NAMES(MAX_REACTANTS+1,production_reaction_id))
      do compound=MAX_REACTANTS+2,MAX_COMPOUNDS
        tmp_name = REACTION_COMPOUNDS_NAMES(compound,production_reaction_id)
        if (tmp_name.ne.'') then
          reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
        endif
      enddo
      write(12,*) trim(reaction_line)
    endif
  enddo
  close(12)

  ! CHECK that reactions are equilibrated (for prime elements)
  do reaction=1,nb_reactions
    do element=1,NB_PRIME_ELEMENTS
      left_sum = 0
      do compound=1,MAX_REACTANTS
        tmp_name = REACTION_COMPOUNDS_NAMES(compound,reaction)
        if (tmp_name.ne.'') then
          left_sum = left_sum + species_composition(element,REACTION_COMPOUNDS_ID(compound, reaction))
        endif
      enddo
      
      right_sum = 0
      do compound=MAX_REACTANTS+1,MAX_COMPOUNDS
        tmp_name = REACTION_COMPOUNDS_NAMES(compound,reaction)
        if (tmp_name.ne.'') then
          right_sum = right_sum + species_composition(element,REACTION_COMPOUNDS_ID(compound, reaction))
        endif
      enddo
      
      if (left_sum.ne.right_sum) then
        write(Error_Unit,'(a,i0,a,a)') 'Error: The reaction ',REACTION_ID(reaction), ' is not equilibrated in ',&
                                        trim(element_name(element))
        call exit(15)
      endif
    enddo
  enddo

  ! CHECK that reactions are equilibrated (for charge)
  do reaction=1,nb_reactions
    left_sum = 0
    do compound=1,MAX_REACTANTS
      tmp_name = REACTION_COMPOUNDS_NAMES(compound,reaction)
      if (tmp_name.ne.'') then
        left_sum = left_sum + SPECIES_CHARGE(REACTION_COMPOUNDS_ID(compound, reaction))
      endif
    enddo
    
    right_sum = 0
    do compound=MAX_REACTANTS+1,MAX_COMPOUNDS
      tmp_name = REACTION_COMPOUNDS_NAMES(compound,reaction)
      if (tmp_name.ne.'') then
        right_sum = right_sum + SPECIES_CHARGE(REACTION_COMPOUNDS_ID(compound, reaction))
      endif
    enddo
    
    if (left_sum.ne.right_sum) then
      write(Error_Unit,'(a,i0,a,a)') 'Error: The reaction ',REACTION_ID(reaction), ' is not equilibrated in electric charge.'
      call exit(15)
    endif
  enddo

  ! Check reactions with alpha equal 0
  open(12, file=information_file, position='append')
  alpha_equal_zero = .false.
  write(12, *) ' ### CHECK ### reactions with alpha = 0'
  do reaction=1,nb_reactions
    if (RATE_A(reaction).eq.0.d0) then
      alpha_equal_zero = .true.
      write(12,'(a,i0,a,a)') 'Warning: The reaction ',REACTION_ID(reaction), ' has an alpha = 0'
    endif
  enddo
  close(12)
  
  if (alpha_equal_zero) then
    write(Error_Unit,'(a,i0,a,2(es10.3e2,a))') 'Error: Some reactions have alpha = 0. Please change alpha or comment the '
    write(Error_Unit,'(a,i0,a,2(es10.3e2,a))') 'reaction with "!" because this can cause errors in '
    write(Error_Unit,'(a,i0,a,2(es10.3e2,a))') 'branching ratio. Refers to "info.out" for the list of reactions.'
  call exit(16)
  endif
  
  ! Check if tmin < tmax for all reactions (else, trange is not correctly defined)
  do reaction=1,nb_reactions
    if (REACTION_TMIN(reaction).gt.REACTION_TMAX(reaction)) then
      write(Error_Unit,'(a,i0,a,2(es10.3e2,a))') 'Error: The reaction ',REACTION_ID(reaction), ' has Tmin(',&
                       REACTION_TMIN(reaction),') > Tmax(',REACTION_TMAX(reaction),')'
      call exit(15)
    endif
  enddo
  
  ! Check reactions with the same reaction ID. We want them to have the same reactants and products. We also want them to have
  !! complementary temperature ranges
  do reaction=1,nb_reactions-1
    ref_id = REACTION_ID(reaction)
    
    ! For each reaction, we check any other reaction (above it) that have the same ID
    ! If so, we compare compounds that MUST be equal
    
    do reaction2=reaction+1,nb_reactions
      if (REACTION_ID(reaction2).eq.ref_id) then
        ! Check that they have the same reactants and products.
        if (any(REACTION_COMPOUNDS_ID(1:MAX_COMPOUNDS,reaction).ne.REACTION_COMPOUNDS_ID(1:MAX_COMPOUNDS,reaction2))) then
          write(Error_Unit,'(a,i0,a)') 'Error: The reactions with ID=',REACTION_ID(reaction), ' have different compounds.'
          call exit(17)
        endif
        
        ! Check that temperature ranges do not overlap
        if (REACTION_TMIN(reaction).gt.REACTION_TMIN(reaction2)) then
          ! We ensure that range1 is inferior to range2
          range1_max = REACTION_TMAX(reaction2)
          range2_min = REACTION_TMIN(reaction)
        else
          range1_max = REACTION_TMAX(reaction)
          range2_min = REACTION_TMIN(reaction2)
        endif
        
        if (range2_min.lt.range1_max) then
          write(Error_Unit,'(a,i0,a)') 'Error: The reactions with ID=',REACTION_ID(reaction), ' have overlapping temperature ranges'
          write(Error_Unit,'(a,4(es10.3e2,a))') 'range1 = [',REACTION_TMIN(reaction),';',REACTION_TMAX(reaction),&
                           '] ; range2 = [',REACTION_TMIN(reaction2),':',REACTION_TMAX(reaction2),']'
          call exit(17)
        endif
        
      endif
    enddo
  enddo
  
  ! Check gas neutral species
  do species=1,nb_gaseous_species
    ! We are only interested in neutral species
    if (SPECIES_CHARGE(species).ne.0) then
      cycle
    endif
  
    ! Grain0 represent a grain, thus you can't have a grain on the surface of itself, we skip this one.
    if ((species_name(species).eq.YGRAIN).or.species_name(species).eq.'XH') then
      cycle
    endif
    
    ! Check that each gas species has a grain equivalent (J+name)
    tmp_name = 'J'//trim(species_name(species))
    
    no_grain_equivalent = .true.
    do i = nb_gaseous_species+1,nb_species
      if (species_name(i).eq.tmp_name) then
        no_grain_equivalent = .false.
      endif
    enddo
    
    if (no_grain_equivalent) then
      write(Error_Unit,'(4a)') 'Error: The species ',trim(species_name(species)), ' has no grain equivalent : ',trim(tmp_name)
      call exit(18)
    endif
    
    ! Check that each gas neutral species has an ITYPE=99 reaction
    no_itype = .true. !<There is no itype 99 reaction while we do not find a reactant of one reaction witht he correct ID
    do reaction=type_id_start(99),type_id_stop(99)
      if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS,reaction).eq.species)) then
        no_itype = .false.
      endif
    enddo
    
    if (no_itype) then
      write(Error_Unit,'(4a)') 'Error: The species ',trim(species_name(species)), ' has no adsorption reaction (ITYPE=99)'
      call exit(19)
    endif
  enddo
  
  ! Check surface species
  do species=nb_gaseous_species+1,nb_species
    
    ! Check that each surface species possess one desorption reaction (for ITYPE 15, 16, 66, 67)
    no_itype = .true. !< There is no itype 15 reaction while we do not find a reactant of one reaction with the correct ID
    do reaction=type_id_start(15),type_id_stop(15)
      if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS,reaction).eq.species)) then
        no_itype = .false.
      endif
    enddo
    
    if (no_itype) then
      write(Error_Unit,'(4a)') 'Error: The species ',trim(species_name(species)), ' has no desorption reaction (ITYPE=15)'
      call exit(20)
    endif
    
    ! Check that each surface species possess one desorption reaction (for ITYPE 15, 16, 66, 67)
    no_itype = .true. !< There is no itype 16 reaction while we do not find a reactant of one reaction with the correct ID
    do reaction=type_id_start(16),type_id_stop(16)
      if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS,reaction).eq.species)) then
        no_itype = .false.
      endif
    enddo
    
    if (no_itype) then
      write(Error_Unit,'(4a)') 'Error: The species ',trim(species_name(species)), ' has no desorption reaction (ITYPE=16)'
      call exit(20)
    endif
    
    ! Check that each surface species possess one desorption reaction (for ITYPE 15, 16, 66, 67)
    no_itype = .true. !< There is no itype 66 reaction while we do not find a reactant of one reaction with the correct ID
    do reaction=type_id_start(66),type_id_stop(66)
      if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS,reaction).eq.species)) then
        no_itype = .false.
      endif
    enddo
    
    if (no_itype) then
      write(Error_Unit,'(4a)') 'Error: The species ',trim(species_name(species)), ' has no desorption reaction (ITYPE=66)'
      call exit(20)
    endif
    
    ! Check that each surface species possess one desorption reaction (for ITYPE 15, 16, 66, 67)
    no_itype = .true. !< There is no itype 67 reaction while we do not find a reactant of one reaction with the correct ID
    do reaction=type_id_start(67),type_id_stop(67)
      if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS,reaction).eq.species)) then
        no_itype = .false.
      endif
    enddo
    
    if (no_itype) then
      write(Error_Unit,'(4a)') 'Error: The species ',trim(species_name(species)), ' has no desorption reaction (ITYPE=67)'
      call exit(20)
    endif
  enddo
  
  ! Check reactions id indexes (if they overlap, or miss some indexes)
  do i=0,MAX_NUMBER_REACTION_TYPE-1
    if ((type_id_start(i).eq.0).and.(type_id_stop(i).eq.0)) then
      cycle ! Reaction type not defined
    endif
    
    if ((type_id_start(i).gt.type_id_stop(i))) then
      write(Error_Unit,'(3(a,i0),a)') 'Error: For ITYPE=',i, ' index start > index stop (', type_id_start(i),&
                                    ' > ', type_id_stop(i),')'
      call exit(21)
    endif
  enddo
  
  ! To reconstruct the range of indexes only by sub-ranges of each ITYPE. We must then get the total range or reactions ID
  ! We assume that type_id_start(j) < type_id_stop(j)
  id_start = -1 !< We must not initialize to 0 because 0 is expected to be the final value. Thus, to avoid problems, -1 is a safe value
  id_stop = -1
  do i=0,MAX_NUMBER_REACTION_TYPE-1
    if (type_id_stop(i).eq.0) then
      cycle ! Reaction type not defined
    endif
    
    !< We must not initialize to 0 because 0 is expected to be the final value. Thus, to avoid problems, -1 is a safe value
    if (id_start.eq.-1) then
      id_start = type_id_start(i)
    endif
    
    if (id_stop.eq.-1) then
      id_stop = type_id_stop(i)
    endif
    
    do j=0,MAX_NUMBER_REACTION_TYPE-1
      if ((id_stop+1).eq.type_id_start(j)) then
        id_stop = type_id_stop(j)
      else if ((id_start-1).eq.type_id_stop(j)) then
        id_start = type_id_start(j)
      endif
    enddo
  enddo
  
  
  if (id_start.gt.1) then 
    write(Error_Unit,'(a,i0,a)') 'Error: Reaction ID = ',id_start-1, " doesn't belong to any reaction type when looking"
    write(Error_Unit,'(a)')      '       into type_id_start(:) and type_id_stop(:).'
    write(Error_Unit,'(a,2(i0,a))') 'Reaction range indexed : [',id_start,';', id_stop, ']'
    write(Error_Unit,*) REACTION_COMPOUNDS_NAMES(:,id_start-1)
    call exit(21)
  endif
  if (id_stop.ne.maxval(type_id_stop(0:MAX_NUMBER_REACTION_TYPE-1))) then
    write(Error_Unit,'(a,i0,a)') 'Error: Reaction ID = ',id_stop+1, " doesn't belong to any reaction type when looking"
    write(Error_Unit,'(a)')      '       into type_id_start(:) and type_id_stop(:).'
    write(Error_Unit,'(a,2(i0,a))') 'Reaction range indexed: [',id_start,';', id_stop, ']'
    call exit(21)
  endif
  
! Above tests are done only if IS_TEST=1. 
endif

end subroutine preliminary_tests

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
real(double_precision), dimension(:), allocatable :: temp_abundances !< Temporary variable to store abundances summed over all 1D points

! Initializing global time
current_time = 0.d0

! Read list of prime elements, including their atomic mass (in AMU)
call read_element_in()

! Get various size needed for allocatable arrays
call get_array_sizes()

! Read simulation parameters. Need to reed it before initializing arrays because some of the arrays sizes are carved into it.
call read_parameters_in()

! Initialize structure evolution
select case(IS_STRUCTURE_EVOLUTION)
  case(0)
    get_structure_properties => get_structure_properties_fixed

  case(1)
    call init_structure_evolution()
    
    get_structure_properties => get_structure_properties_table
    
  case default
    write(error_unit,*) 'The is_structure_evolution="', IS_STRUCTURE_EVOLUTION,'" cannot be found.'
    write(error_unit,*) 'Values possible : 0: no ; 1: yes'
    write(error_unit, '(a)') 'Error in structure: subroutine init_structure_evolution' 
    call exit(9)
end select

! Initialize grain temperature
select case(GRAIN_TEMPERATURE_TYPE)
  case('fixed') !Tgrain = initial_Tgrain
    get_grain_temperature => get_grain_temperature_fixed

  case('table') ! Tgrain read from a table
    get_grain_temperature => get_grain_temperature_table

  case('gas') ! Tgrain = Tgas
    get_grain_temperature => get_grain_temperature_gas

  case('computed') ! Tgrain computed consistently
    get_grain_temperature => get_grain_temperature_computed
    
  case default
    write(error_unit,*) 'The GRAIN_TEMPERATURE_TYPE="', GRAIN_TEMPERATURE_TYPE,'" cannot be found.'
    write(error_unit,*) 'Values possible : fixed, table, gas, computed'
    write(error_unit, '(a)') 'Error in subroutine initialisation.' 
    call exit(10)
end select

! Init global allocatable arrays. From now on, we can read data files
call initialize_global_arrays()

!0D, 1D_diff, 1D_no_diff
! Initialize structure pointers
select case(STRUCTURE_TYPE)
  case('0D') ! No structure at all. Point towards routines that do almost nothing.
    get_timestep => get_timestep_0D
    structure_diffusion => structure_no_diffusion

  case('1D_diff') ! 1D structure with species diffusion
    call init_1D_static()
    get_timestep => get_timestep_1D_diff
    structure_diffusion => structure_diffusion_1D
  
  case('1D_no_diff') ! 1D structure but without species diffusion
    call init_1D_static()
    get_timestep => get_timestep_0D
    structure_diffusion => structure_no_diffusion
    
  case default
    write(error_unit,*) 'The STRUCTURE_TYPE="', STRUCTURE_TYPE,'" cannot be found.'
    write(error_unit,*) 'Values possible : 0D, 1D_diff, 1D_no_diff'
    write(error_unit, '(a)') 'Error in subroutine initialisation.' 
    call exit(21)
end select


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

! Initialization of elemental/chemical quantities
call index_datas()

! Calculate the initial abundances for all elements that compose 
allocate(temp_abundances(nb_species))
temp_abundances(1:nb_species) = sum(abundances(1:nb_species,1:spatial_resolution), dim=2)
call get_elemental_abundance(all_abundances=temp_abundances(1:nb_species), el_abundances=INITIAL_ELEMENTAL_ABUNDANCE)

! Store initial elemental abundances
call write_elemental_abundances(filename='elemental_abundances.out', el_abundances=INITIAL_ELEMENTAL_ABUNDANCE)

! Recompute initial_dtg_mass_ratio to remove He
! In the following, initial_dtg_mass_ratio is used as a H/dust mass ratio
do i=1,NB_PRIME_ELEMENTS
  if (species_name(PRIME_ELEMENT_IDX(I)).EQ.YHE) then
    initial_dtg_mass_ratio = initial_dtg_mass_ratio*(1.d0+4*INITIAL_ELEMENTAL_ABUNDANCE(I))
  endif
enddo

! Compute the grain abundance
GTODN = (4.d0 * PI * GRAIN_DENSITY * grain_radius * grain_radius * grain_radius) / (3.d0 * initial_dtg_mass_ratio * AMU)

abundances(INDGRAIN,1:spatial_resolution) = 1.0 / GTODN ! TODO do we must divide by spatial_resolution in the 1D case???

! Set the electron abundance via conservation===========
! And check at the same time that nls_init has the same elemental abundance
! as nls_control
! Make comparison for the sum of abundances over 1D dimension
if (spatial_resolution.eq.1) then
  call check_conservation(abundances(1:nb_species, 1))
  ! Make a check only routine for 1D case were it is very complicated to modify any abundances since its spread throughout 1D points

  ! 1D physical structure (nls_phys_1D)
  call get_structure_properties(time=current_time, & ! Inputs
                              Av=visual_extinction(1), density=H_number_density(1), & ! Outputs
                              gas_temperature=gas_temperature(1)) ! Outputs
endif

call get_grain_temperature(time=current_time, gas_temperature=gas_temperature(1), Av=visual_extinction(1), & ! Inputs
      grain_temperature=dust_temperature(1)) ! Outputs

! We can't add input variables in dlsodes called routines, so we must store the values as global variables
actual_gas_temp = gas_temperature(1)
actual_dust_temp = dust_temperature(1)
actual_av = visual_extinction(1)
actual_gas_density = H_number_density(1)

! Write species name/index correspondance
call write_species()

! Initialize indices of reactants and products 
call set_chemical_reactants()

! Initialize the arrays that list, for each species, the reactions using it as a reactant
!! max_reactions_same_species is set here. nb_reactions_using_species and relevant_reactions array are set here.
call init_relevant_reactions()

! Calculate the optimum number for temporary solving-arrays in ODEPACK, based on the number of non-zeros values in 
!! the jacobian
call count_nonzeros()

! Write information about the code at the very end of initialisation
call write_general_infos()

! Do preliminary tests, before starting the integration
! Writing information in 'info.out', appending to the file created by 'write_general_infos'
call preliminary_tests()

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
  if ((KSUM.EQ.1).AND.(SPECIES_CHARGE(J).EQ.0).AND.&
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
  IF (species_name(J).EQ.YE) then
    INDEL=J
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

! --- Set reference species
do I=1,nb_species 
  ! ------ Calculate elemental_masses
  MSUM=0.d0
  do K=1,NB_PRIME_ELEMENTS 
    MSUM=MSUM+elemental_mass(K)*species_composition(K,I) 
  enddo 
  SPECIES_MASS(I)=MSUM
  if (species_name(I).EQ.YE) SPECIES_MASS(I) = ELECTRON_MASS ! electron mass in amu
  if (species_name(I).EQ.YGRAIN .OR. species_name(I).EQ.'GRAIN-      ')&
  SPECIES_MASS(I)=4.0*PI*grain_radius*grain_radius*grain_radius*GRAIN_DENSITY/3.0/AMU
enddo

! Initialize the Av/NH ratio
! Can be scaled for different dust/gas ratios
! Here initial_dtg_elemental_mass_ratio is the original dust/gas elemental_mass ratio (from nls_control.d)
! initial_dtg_elemental_mass_ratio is changed later into the dust/hydrogen elemental_mass ratio

AV_NH_ratio = 5.34d-22 * (initial_dtg_mass_ratio / 1.d-2)

! Find ITYPE first and last reactions===================================
do I=0,MAX_NUMBER_REACTION_TYPE-1
  type_id_start(I)=0
  type_id_stop(I)=0
  do J=1,nb_reactions
    if ((REACTION_TYPE(J).EQ.I).AND.(type_id_start(I).EQ.0)) type_id_start(I)=J
    if (REACTION_TYPE(J).EQ.I) type_id_stop(I)=J
  enddo
enddo

! Find the index of CO, H2, H, He and grain0
do i=1,nb_species
  if (species_name(i).eq.YH2) INDH2=i
  if (species_name(i).eq.YH) INDH=i
  if (species_name(i).eq.YCO) INDCO=i
  if (species_name(i).eq.YHE) INDHE=i
  if (species_name(i).eq.YGRAIN) INDGRAIN=i
enddo

! Compute nb_sites_per_grain = number of sites per grain
nb_sites_per_grain = SURFACE_SITE_DENSITY*4.d0*PI*grain_radius**2

! Initialise reaction rates=============================================
call init_reaction_rates()

return 
end subroutine index_datas

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
    character(len=11), dimension(MAX_COMPOUNDS,nb_reactions) :: GSREAD
    character(len=11), dimension(nb_species) :: GSPEC

    real(double_precision) :: cond
    integer :: i, j,k,l,n4, n5, n6
    
    character(len=80) :: filename !< name of the file to be read
    character(len=80) :: line_format
    character(len=200) :: line
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
      if (SPECIES_CHARGE(I).EQ.0) then
        STICK = sticking_coeff_neutral
      else if (SPECIES_CHARGE(I).GT.0) then 
        STICK = sticking_coeff_positive 
      else ! (SPECIES_CHARGE(I).LT.0)
        STICK = sticking_coeff_negative 
      endif
      !         if (species_name(I).EQ.YH2)      STICK=0.d0 
      !         if (species_name(I).EQ.YHE)      STICK=0.d0 
      !         if (species_name(I).EQ.YH)       STICK=0.d0 
      if (species_name(I).EQ.YHEP)     STICK=0.d0 
      if (species_name(I).EQ.'e-         ')      STICK=0.d0
      if (species_name(I).EQ.'H+         ')     STICK=0.d0
      if (species_name(I).EQ.YGRAIN)   STICK=0.d0
      if (species_name(I).EQ.'GRAIN-     ') STICK=0.d0
      !         if (species_name(I).EQ.'H-')     STICK=0.d0
      !         if (species_name(I).EQ.'H2+')    STICK=0.d0

      if (I.GT.nb_gaseous_species) STICK=0.d0
      ACC_RATES_PREFACTOR(I)=COND*STICK/SQRT(SPECIES_MASS(I))
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
      
      ! Definition of the line format, common to gas_reaction.in and grain_reaction.in
      write(line_format, '(a,i0,a,i0,a)') '(', MAX_REACTANTS, 'A11,4x,', MAX_PRODUCTS, 'A11,D9.2)'

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
!~           write(*,*) "i=",i
!~           write(*,*) line_format
!~           write(*,*) line
          read(line, line_format) (GSread(L,i),L=1,MAX_COMPOUNDS),REA5(i)
        
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
      BINDING_ENERGY(I)=0.d0
      DIFFUSION_BARRIER(I)=0.d0
      GAP_ENERGY_BANDS(I)=0.d0
      FORMATION_ENTHALPY(I)=0.d0
      do J=1,NGS
        if (species_name(I).EQ.GSPEC(J)) then
          SMASS(I)=dble(INT1(J))
          BINDING_ENERGY(I)=REA1(J)
          DIFFUSION_BARRIER(I)=REA2(J)
          GAP_ENERGY_BANDS(I)=REA3(J)
          FORMATION_ENTHALPY(I)=REA4(J)
          if ((species_name(I).NE.YJH).AND.(species_name(I).NE.YJH2).AND.(DIFF_BINDING_RATIO.GE.0.d0)) then
            DIFFUSION_BARRIER(I)=DIFF_BINDING_RATIO*BINDING_ENERGY(I)
          endif
        endif
      enddo
      !IF(species_name(I) == 'JN2O2      ') write(*,*) BINDING_ENERGY(I)
    enddo

    ! For each reaction, we search if there is an activation energy defined for it.
    ! the "all()" function can compare array element by element to ensure that everything is equal one by one. usefull to find 
    ! if we have the good reaction
    do I=1,nb_reactions
      ACTIVATION_ENERGY(I)=0.d0
      do J=1,NEA
        if (REACTION_COMPOUNDS_NAMES(MAX_REACTANTS+1,I)(:1).EQ.'J') then
          if (all(REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS,I).EQ.GSread(1:MAX_COMPOUNDS,J))) then
            ACTIVATION_ENERGY(I) = REA5(J)
          endif
        else
          if (all((REACTION_COMPOUNDS_NAMES(1:MAX_REACTANTS,I).EQ.GSread(1:MAX_REACTANTS,J))).AND.&
          (all(REACTION_COMPOUNDS_NAMES(MAX_REACTANTS+1:MAX_COMPOUNDS,I).EQ.GSread(MAX_REACTANTS+1:MAX_COMPOUNDS,J)(2:)))) then
            ACTIVATION_ENERGY(I) = REA5(J)
          endif
        endif
      enddo
      !IF(REACTION_COMPOUNDS_NAMES(4,i) == 'JO2H       ') write(*,*)  REACTION_COMPOUNDS_NAMES(:,i), ACTIVATION_ENERGY(i)
    enddo

    ! Set up constants, quantum rate info===================================
    do I=1,nb_species
      VIBRATION_FREQUENCY(I)=0.d0
      TUNNELING_RATE_TYPE_1(I)=0.d0
      TUNNELING_RATE_TYPE_2(I)=0.d0
      ! ------ For species which have been assigned surface info, SMASS=/=0
      if (SMASS(I).NE.0) then
        SMA=dble(SMASS(I))
        ! --------- Set characteristic frequency
        VIBRATION_FREQUENCY(I)=SQRT(2.0d0*K_B/PI/PI/AMU * SURFACE_SITE_DENSITY*BINDING_ENERGY(I)/SMA)
        ! --------- Set quantum rates
        if (GAP_ENERGY_BANDS(I).GE.1.0D-38) then
          TUNNELING_RATE_TYPE_1(I)=GAP_ENERGY_BANDS(I)*K_B/4.0d0/H_BARRE/nb_sites_per_grain
        else
          TUNNELING_RATE_TYPE_1(I)=0.d0
        endif
        TUNNELING_RATE_TYPE_2(I) = VIBRATION_FREQUENCY(I) / nb_sites_per_grain * &
                 EXP(-2.0d0*DIFFUSION_BARRIER_THICKNESS/H_BARRE*SQRT(2.0d0*AMU*SMA*K_B*DIFFUSION_BARRIER(I)))
      endif
    enddo

    ! === Cycle all reactions
    do J=1,nb_reactions

      ! ------ Initialise all branching_ratio rate factors, and get species 1 & 2
      branching_ratio(J)=1.0d0
      reactant_1_idx(J)=0
      reactant_2_idx(J)=0
      do I=1,nb_species
        if (REACTION_COMPOUNDS_NAMES(1,J).EQ.species_name(I)) reactant_1_idx(J)=I
        if (REACTION_COMPOUNDS_NAMES(2,J).EQ.species_name(I)) reactant_2_idx(J)=I
      enddo

      ! === ITYPE 14 AND 21 - SURFACE REACTIONS
      if (REACTION_TYPE(J).EQ.14 .OR. REACTION_TYPE(J).EQ.21) then
        NPATH=0

        ! ------ Check for branching
        do K=1,nb_reactions
           if(REACTION_TYPE(K).EQ.REACTION_TYPE(J)) then
             if (((REACTION_COMPOUNDS_NAMES(1,J).EQ.REACTION_COMPOUNDS_NAMES(1,K)).AND.&
                  (REACTION_COMPOUNDS_NAMES(2,J).EQ.REACTION_COMPOUNDS_NAMES(2,K))).OR.&
                 ((REACTION_COMPOUNDS_NAMES(2,J).EQ.REACTION_COMPOUNDS_NAMES(1,K)).AND.&
                  (REACTION_COMPOUNDS_NAMES(1,J).EQ.REACTION_COMPOUNDS_NAMES(2,K)))) then
                if (REACTION_COMPOUNDS_NAMES(4,K)(:1).EQ.'J          ') NPATH=NPATH+1
             endif
           endif
        enddo

      ! ------ Branching ratio
      if (NPATH.EQ.0) then
        branching_ratio(J)=0.d0
      else
        branching_ratio(J)=branching_ratio(J)/dble(NPATH)
      endif

      ! ------ Factor of 2 for same species reactions
      if (reactant_1_idx(J).EQ.reactant_2_idx(J)) branching_ratio(J)=branching_ratio(J)/2.0d0

      ! ------ Calculate evaporation fraction
      NEVAP=0
      do K=1,nb_reactions
        if ((REACTION_COMPOUNDS_NAMES(4,J)(:1).EQ.'J          ').AND.(RATE_A(K).NE.0.d0)) then
          if ((REACTION_COMPOUNDS_NAMES(1,J).EQ.REACTION_COMPOUNDS_NAMES(1,K)).AND.&
          (REACTION_COMPOUNDS_NAMES(2,J).EQ.REACTION_COMPOUNDS_NAMES(2,K)).AND.&
          (REACTION_COMPOUNDS_NAMES(4,J)(2:).EQ.REACTION_COMPOUNDS_NAMES(4,K)).AND.&
          (REACTION_COMPOUNDS_NAMES(5,J)(2:).EQ.REACTION_COMPOUNDS_NAMES(5,K)).AND.&
          (REACTION_COMPOUNDS_NAMES(6,J)(2:).EQ.REACTION_COMPOUNDS_NAMES(6,K)).AND.&
          (REACTION_COMPOUNDS_NAMES(4,K)(:1).NE.'J          ')) NEVAP=NEVAP+1
        endif
        if ((REACTION_COMPOUNDS_NAMES(4,J)(:1).NE.'J          ').AND.(RATE_A(J).NE.0.d0)) then
          if ((REACTION_COMPOUNDS_NAMES(1,J).EQ.REACTION_COMPOUNDS_NAMES(1,K)).AND.&
          (REACTION_COMPOUNDS_NAMES(2,J).EQ.REACTION_COMPOUNDS_NAMES(2,K)).AND.&
          (REACTION_COMPOUNDS_NAMES(4,J).EQ.REACTION_COMPOUNDS_NAMES(4,K)(2:)).AND.&
          (REACTION_COMPOUNDS_NAMES(5,J).EQ.REACTION_COMPOUNDS_NAMES(5,K)(2:)).AND.&
          (REACTION_COMPOUNDS_NAMES(6,J).EQ.REACTION_COMPOUNDS_NAMES(6,K)(2:)).AND.&
          (REACTION_COMPOUNDS_NAMES(4,K)(:1).EQ.'J          ')) NEVAP=NEVAP+1
        endif
      enddo

      N4=0
      N5=0
      N6=0
      do I=nb_gaseous_species+1,nb_species
        if (REACTION_COMPOUNDS_NAMES(4,J)(:1).EQ.'J          ') then
          if (REACTION_COMPOUNDS_NAMES(4,J).EQ.species_name(I)) N4=I
          if (REACTION_COMPOUNDS_NAMES(5,J).EQ.species_name(I)) N5=I
          if (REACTION_COMPOUNDS_NAMES(6,J).EQ.species_name(I)) N6=I
        endif
        if ((REACTION_COMPOUNDS_NAMES(4,J)(:1).NE.'J          ').AND.&
        (REACTION_COMPOUNDS_NAMES(4,J)(:1).NE.'X          ')) then
        if (REACTION_COMPOUNDS_NAMES(4,J).EQ.species_name(I)(2:)) N4=I
        if (REACTION_COMPOUNDS_NAMES(5,J).EQ.species_name(I)(2:)) N5=I
        if (REACTION_COMPOUNDS_NAMES(6,J).EQ.species_name(I)(2:)) N6=I
      endif
    enddo

    DHFSUM=FORMATION_ENTHALPY(reactant_1_idx(J))+FORMATION_ENTHALPY(reactant_2_idx(J))-FORMATION_ENTHALPY(N4)
    if (N5.NE.0) DHFSUM=DHFSUM-FORMATION_ENTHALPY(N5)
    if (N6.NE.0) DHFSUM=DHFSUM-FORMATION_ENTHALPY(N6)
    ! ------ Convert from kcal to J, from J to K
    DHFSUM=DHFSUM*4.184d03/1.38054D-23
    ! ------ Convert from #moles-1 to #reactions-1
    DHFSUM=DHFSUM/AVOGADRO

    DHFSUM=DHFSUM+ACTIVATION_ENERGY(J)

    SUM1=BINDING_ENERGY(N4)
    if (N5.NE.0) SUM1=MAX(BINDING_ENERGY(N4),BINDING_ENERGY(N5))
    if (N6.NE.0) SUM1=MAX(BINDING_ENERGY(N4),BINDING_ENERGY(N5),BINDING_ENERGY(N6))

    ATOMS=0
    do K=1,NB_PRIME_ELEMENTS
      ATOMS=ATOMS+species_composition(K,N4)
    enddo

    SUM2=1.0d0-(SUM1/DHFSUM)
    if (ATOMS.EQ.2) SUM2=SUM2**(3*ATOMS-5)
    if (ATOMS.GT.2) SUM2=SUM2**(3*ATOMS-6)
    !         SUM2=SUM2**(3*ATOMS-6)
    SUM2=VIB_TO_DISSIP_FREQ_RATIO*SUM2
    EVFRAC=SUM2/(1+SUM2)

    !        V.W. Jul 2006 f=evfrac=0.009 for H2O (Kroes & Anderson 2006) 
    !         if (REACTION_COMPOUNDS_NAMES(4,J).EQ.'H2O     ') then
    !                 EVFRAC=0.009
    !                EVFRAC_H2O=0.009
    !         endif

    BADFLAG=0
    if (FORMATION_ENTHALPY(reactant_1_idx(J)).LE.-999.0) then
      EVFRAC=0.d0
      BADFLAG=BADFLAG+1
    endif
    if (FORMATION_ENTHALPY(reactant_2_idx(J)).LE.-999.0) then
      EVFRAC=0.d0
      BADFLAG=BADFLAG+1
    endif
    if (FORMATION_ENTHALPY(N4).LE.-999.0) then
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

    if (REACTION_COMPOUNDS_NAMES(4,J)(:1).EQ.'J          ') then
      EVFRAC=1.0d0-EVFRAC
    endif

    branching_ratio(J)=branching_ratio(J)*EVFRAC

    ! ------ Calculate quantum activation energy
    REDMAS = SMASS(reactant_1_idx(J)) * SMASS(reactant_2_idx(J)) / (SMASS(reactant_1_idx(J)) + SMASS(reactant_2_idx(J)))
    SURF_REACT_PROBA(J) = 2.0d0 * CHEMICAL_BARRIER_THICKNESS/H_BARRE * SQRT(2.0d0*AMU*REDMAS*K_B*ACTIVATION_ENERGY(J))
  endif

  ! === ITYPE 16 - C.R. DESORPTION
  if (REACTION_TYPE(J).EQ.16) then
    if (SMASS(reactant_1_idx(J)).EQ.0) branching_ratio(J)=0.d0
  endif

  ! === ITYPE 99 - ACCRETION
  if (REACTION_TYPE(J).EQ.99) then
    ! ------ Save tag of resultant grain surface species
    do I=1,nb_species
      if (REACTION_COMPOUNDS_NAMES(4,J).EQ.species_name(I)) reactant_2_idx(J)=I
    enddo
  endif

enddo

! === Zero dummy H2 formation rxns, if necc.
!      if (IS_GRAIN_REACTIONS.NE.0) then
!         branching_ratio(1)=0.d0
!         branching_ratio(2)=0.d0
!      endif

return
end subroutine init_reaction_rates

end module nautilus_main
