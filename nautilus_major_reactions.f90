program nautilus_major_reactions

use numerical_types
use iso_fortran_env
use utilities
use nautilus_main

implicit none

! Parameters
real(double_precision), parameter :: PERCENTAGE_THRESHOLD = 1.d0 !< Percentage below which the reaction will not be displayed

! Locals
character(len=80) :: filename_output
integer :: species, output, reaction, i ! index for loops
logical :: isDefined
real(double_precision) :: tmp !< temporary variable
character(len=80) :: output_format

integer :: nb_sample_1D = 1

! /!\ Variable names with _out are variable that already exist in the nautilus code, but here they are arrays, one value per output.

real(double_precision), dimension(:,:,:), allocatable :: abundances_out !< abundances_out over time for each species. (nb_outputs, nb_species)

real(double_precision), dimension(:), allocatable :: time
real(double_precision), dimension(:,:), allocatable :: gas_temperature_out
real(double_precision), dimension(:,:), allocatable :: dust_temperature_out
real(double_precision), dimension(:,:), allocatable :: density
real(double_precision), dimension(:,:), allocatable :: visual_extinction_out !< visual extinction
real(double_precision), dimension(:,:), allocatable :: zeta

! For rates
real(double_precision), allocatable, dimension(:,:) :: reaction_rates_out ! (nb_outputs, nb_reactions)

! User asked values
logical :: change_time = .true. !< If true, ask the user for a value
logical :: change_species = .true. !< If true, ask the user for a value
logical :: change_space = .true. !< If true, ask the user for a value
logical :: print_types = .true. !< print different reaction types the first time
logical :: wrong_species, wrong_output, wrong_1D, wrong_action !< Flags for while loops when asking the user something
integer :: user_action !< Ask the user what he wants to do after the first run
character(len=11) :: user_species !< The species designed by the user
integer :: user_species_id !< corresponding id of the desired species of the user
integer :: output_id !< designed output id by the user
integer :: user_1D_id !< designed spatial id by the user

! For outputs
real(double_precision), allocatable, dimension(:) :: destructions !< production rates for one given species. Reactions where this species is not involved have 0 value
integer, allocatable, dimension(:) :: destructions_id !< corresponding ID sorted from the lest important to the most importants at the end
real(double_precision), allocatable, dimension(:) :: productions !< destruction rates for one given species. Reactions where this species is not involved have 0 value
integer, allocatable, dimension(:) :: productions_id !< corresponding ID sorted from the lest important to the most importants at the end
real(double_precision) :: destructions_sum
real(double_precision) :: productions_sum
real(double_precision) :: percentage
character(len=80) :: reaction_line

! Initialise all variables from the global_variables module. Only some of them are used here.
call initialisation()

! We calculate the total number of outputs by checking for each file if it exist or not.
nb_outputs = 0
isDefined = .true.
do while(isDefined)
  nb_outputs = nb_outputs + 1
  write(filename_output, '(a,i0.6,a)') 'abundances.',nb_outputs,'.out'
  inquire(file=filename_output, exist=isDefined)

enddo
nb_outputs = nb_outputs - 1

write(*,'(a, i0)') 'Spatial resolution: ', nb_sample_1D
write(*,'(a, i0)') 'Number of time outputs: ', nb_outputs
write(*,'(a, i0)') 'Number of species: ', nb_species
write(*,'(a, i0)') 'Number of reactions: ', nb_reactions

! We allocate the output arrays

allocate(time(nb_outputs))
allocate(gas_temperature_out(nb_sample_1D, nb_outputs))
allocate(dust_temperature_out(nb_sample_1D, nb_outputs))
allocate(density(nb_sample_1D, nb_outputs))
allocate(visual_extinction_out(nb_sample_1D, nb_outputs))
allocate(zeta(nb_sample_1D, nb_outputs))

allocate(abundances_out(nb_outputs, nb_species+1, nb_sample_1D)) ! We create an extra species that will always have an abundance of 1

allocate(reaction_rates_out(nb_outputs, nb_reactions))

! Outputs
allocate(destructions(nb_reactions))
allocate(destructions_id(nb_reactions))
allocate(productions(nb_reactions))
allocate(productions_id(nb_reactions))

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Reading unformatted outputs...'
! We read output files
do output=1,nb_outputs
  write(filename_output, '(a,i0.6,a)') 'abundances.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) time(output)
  read(10) gas_temperature_out(1:nb_sample_1D, output), dust_temperature_out(1:nb_sample_1D, output), &
           density(1:nb_sample_1D, output), visual_extinction_out(1:nb_sample_1D, output), zeta(1:nb_sample_1D, output)
  read(10) abundances_out(output,1:nb_species, 1:nb_sample_1D)
  close(10)
  
  write(filename_output, '(a,i0.6,a)') 'rates.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) reaction_rates_out(output,1:nb_reactions)
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'

! For non existing reactants (whose index is 'nb_species+1') in reactions, we create a new species whose abundance is always 1, so that we can calculate the fluxes 
!! more easily.
abundances_out(1:nb_outputs, nb_species+1, 1:nb_sample_1D) = 1.d0

!######################################################
! User asked section
!######################################################

write(*,*) "Warning by Christophe Cossou:"
write(*,*) "Please check gas density. There can be a factor of 2 and I never know which one is to be used"

! What time ?
20 if (change_time) then
  wrong_output = .true.
  do while (wrong_output)
    write(*,'(a,i0,a)') 'Select the output time (from 1 to ', nb_outputs, '):'
    read(*,*) output_ID
    
    
    if (output_ID.gt.nb_outputs) then
      write(*,'(a,i0)') "Choose output between 1 and ", nb_outputs
    else
      wrong_output = .false.
    endif
  enddo
  change_time = .false.
endif

! Which species ?
30 if (change_species) then
  wrong_species = .true.
  do while (wrong_species)
    write(*,*) 'Select one species:'
    read(*,*) user_species
    
    user_species_id = -1
    do i=1, nb_species
      if (user_species.eq.species_name(i)) then
      user_species_id = i
      endif
    enddo
    
    if (user_species_id.eq.-1) then
      write(*,*) "'", trim(user_species), "' doesn't exist"
    else
      wrong_species = .false.
    endif
  enddo
  change_species = .false.
endif

! What spatial point ?
40 if (change_space) then
  if (nb_sample_1D.ne.1) then
    wrong_1D = .true.
    do while (wrong_1D)
      write(*,'(a,i0,a)') 'Select the spatial point (from 1 to ', nb_sample_1D, '):'
      read(*,*) user_1D_id
      
      
      if ((user_1D_id.gt.nb_sample_1D).or.(user_1D_id.lt.0)) then
      write(*,'(a,i0)') "Choose spatial point between 1 and ", nb_sample_1D
      else
      wrong_1D = .false.
      endif
    enddo
  else
    ! Default value if only one point
    write(*,*) "We are in 0D, skipping spatial point choosing"
    user_1D_id = 1
  endif
  change_space = .false.
endif

!######################################################
! End User asked section
!######################################################

destructions(1:nb_reactions) = 0.d0
productions(1:nb_reactions) = 0.d0
destructions_sum = 0.d0
productions_sum = 0.d0


do reaction=1, nb_reactions
  if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS, reaction).eq.user_species_id)) then
    tmp = 1.d0
    do i=1, MAX_REACTANTS
      ! Skip if blanck species
      if (REACTION_COMPOUNDS_ID(i, reaction).ne.nb_species+1) then
        tmp = tmp * abundances_out(output_id, REACTION_COMPOUNDS_ID(i, reaction), 1)
      endif
    enddo
    destructions(reaction) = reaction_rates_out(output_id, reaction) * tmp
  endif
  
  if (any(REACTION_COMPOUNDS_ID(MAX_REACTANTS+1:MAX_COMPOUNDS, reaction).eq.user_species_id)) then
    tmp = 1.d0
    do i=1, MAX_REACTANTS
      ! Skip if blanck species
      if (REACTION_COMPOUNDS_ID(i, reaction).ne.nb_species+1) then
        tmp = tmp * abundances_out(output_id, REACTION_COMPOUNDS_ID(i, reaction), 1)
      endif
    enddo
    productions(reaction) = reaction_rates_out(output_id, reaction) * tmp
  endif
enddo
  
call get_sorted_index(productions(1:nb_reactions), productions_id(1:nb_reactions))
call get_sorted_index(destructions(1:nb_reactions), destructions_id(1:nb_reactions))

productions_sum = sum(productions(1:nb_reactions))
destructions_sum = sum(destructions(1:nb_reactions))

if (print_types) then
  write(*,*) "---------------------------------------- ITYPE ------------------------------------------"
  write(*,*) " 0     : Gas-phase reactions"
  write(*,*) " 1-2   : Photodissoc/ionisation by cosmic rays"
  write(*,*) " 3     : Gas phase photodissociations/ionisations by UV"
  write(*,*) " 4-8   : Bimolecular gas phase reactions"
  write(*,*) " 10-11 : ad-hoc formation on grains"
  write(*,*) " 15    : Thermal evaporation"
  write(*,*) " 14    : Grain surface reactions"
  write(*,*) " 16    : Cosmic-ray evaporation"
  write(*,*) " 17-18 : Photodissociations by Cosmic rays on grain surfaces"
  write(*,*) " 19-20 : Photodissociations by UV photons on grain surfaces"
  write(*,*) " 21    : Cosmic rays diffusion (Laura)"
  write(*,*) " 30    : Direct formation process X + JY -> JXY"
  write(*,*) " 31    : JC-X->JCX and JCH-X->JCHX"
  write(*,*) " 66    : CO photodesorption by external UV"
  write(*,*) " 67    : CO photodesorption by CR generated UV"
  write(*,*) " 98    : Storage of H2S under a refractory form"
  write(*,*) " 99    : Adsorption on grains"
  print_types = .false.
endif

write(output_format,*)'(i2,4x,a,4x,es13.7,4x,f5.1,"%")'


write(*,*) "--------------------------------- PRODUCTION (cm-3 s-1) ----------------------------------    --------"
percentage = 100.d0
i = nb_reactions
do while(percentage.gt.PERCENTAGE_THRESHOLD)
  reaction = productions_id(i)
  percentage = productions(reaction) / productions_sum * 100.d0
  call display_reaction(REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_line)
  
  write(*, output_format) REACTION_TYPE(reaction), trim(reaction_line), productions(reaction), percentage
  
  ! Decrement index
  i = i - 1
  
  ! Just in case, to avoid negative index for arrays
  if (i.lt.0) then
    percentage = 0.d0
  endif
enddo

write(*,*) "--------------------------------- DESTRUCTION (cm-3 s-1) ---------------------------------    --------"
percentage = 100.d0
i = nb_reactions
do while(percentage.gt.PERCENTAGE_THRESHOLD)
  reaction = destructions_id(i)
  percentage = destructions(reaction) / destructions_sum * 100.d0
  call display_reaction(REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_line)
  
  write(*,output_format) REACTION_TYPE(reaction), trim(reaction_line), destructions(reaction), percentage
  
  ! Decrement index
  i = i - 1
  
  ! Just in case, to avoid negative index for arrays
  if (i.lt.0) then
    percentage = 0.d0
  endif
enddo

! Ask user if he want another run
wrong_action = .true.
do while (wrong_action)
	write(*,*) "What do you want to do?"
	write(*,*) "0:quit ; 1: change all ; 2:change time ; 3: change species ; 4: change spatial point"
	read(*,*) user_action
	if ((user_action.ge.0).and.(user_action.le.4)) then
	  wrong_action = .false.
	else
	  write(*,*) "Error: Action must be between 0 and 4"
	endif
	
	if ((user_action.eq.4).and.(nb_sample_1D.eq.1)) then
	  wrong_action = .true.
	  write(*,*) "/!\ You are in 0D !"
	endif
enddo

select case(user_action)
  case(0)
    call exit(0) ! Exiting normally

  case(1) ! change all
    change_time = .true.
    change_species = .true.
    change_space = .true.
    goto 20

  case(2) ! change time
    change_time = .true.
    goto 20

  case(3) ! change species
    change_species = .true.
    goto 30

  case(4) ! change spatial point
    change_space = .true.
    goto 40
    
  case default
    write(error_unit,*) 'The is_structure_evolution="', IS_STRUCTURE_EVOLUTION,'" cannot be found.'
    write(error_unit,*) 'Values possible : 0: no ; 1: yes'
    write(error_unit, '(a)') 'Error in structure: subroutine init_structure_evolution' 
    call exit(9)
end select


contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Return a string to display a reaction
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine display_reaction(names,reaction_line)
implicit none

! Input
character(len=11), dimension(MAX_COMPOUNDS), intent(in) :: names !<[in] REACTION_COMPOUNDS_NAMES for a given reaction

! Output
character(len=80), intent(out) :: reaction_line !<[out] String that represent the reaction

! Locals
character(len=11) :: tmp_name
integer :: compound

reaction_line = trim(names(1))
do compound=2,MAX_REACTANTS
tmp_name = names(compound)
  if (tmp_name.ne.'') then
    reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
  endif
enddo

reaction_line = trim(reaction_line)//" -> "//trim(names(MAX_REACTANTS+1))

do compound=MAX_REACTANTS+2,MAX_COMPOUNDS
tmp_name = names(compound)
  if (tmp_name.ne.'') then
  reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
  endif
enddo

end subroutine display_reaction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine coming from Numerical Recipes index_sp function
!! adapted for this code. The routine return the ordered indexes corresponding
!! to the element of the input arrays, from lowest to highest. \n\n
!! For a given array [3, 1, 2, 4, 5], the function will return [2, 3, 1, 4, 5]
!! because the lowest element is the 2nd, then come the 3rd, the 1st 
!! and then 4-th and 5-th
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_sorted_index(arr,index)
implicit none
real(double_precision), dimension(:), intent(in) :: arr !<[in] the input array for the sorted process
integer, dimension(:), intent(out) :: index !<[out] the index of the elements, from lowest to highest corresponding values
integer, parameter :: nn=15, nstack=50
real(double_precision) :: a
integer :: n,k,i,j,indext,jstack,l,r
integer, dimension(nstack) :: istack
integer :: tmp_index

if (size(index).ne.size(arr)) then
  write(error_unit,*) 'in get_sorted_index. size are different for the two arguments'
  call exit(99)
endif

n = size(index)

! initialize list of index
do i=1,n
  index(i) = i
enddo

jstack=0
l=1
r=n
do
  if (r-l < nn) then
    do j=l+1,r
      indext=index(j)
      a=arr(indext)
      do i=j-1,l,-1
        if (arr(index(i)) <= a) exit
        index(i+1)=index(i)
      end do
      index(i+1)=indext
    end do
    if (jstack == 0) return
    r=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
  else
    k=(l+r)/2
    
    ! swaping indexes
    tmp_index = index(k)
    index(k) = index(l+1)
    index(l+1) = tmp_index
    
    call icomp_xchg(arr,index(l),index(r))
    call icomp_xchg(arr,index(l+1),index(r))
    call icomp_xchg(arr,index(l),index(l+1))
    i=l+1
    j=r
    indext=index(l+1)
    a=arr(indext)
    do
      do
        i=i+1
        if (arr(index(i)) >= a) exit
      end do
      do
        j=j-1
        if (arr(index(j)) <= a) exit
      end do
      if (j < i) exit
      tmp_index = index(i)
      index(i) = index(j)
      index(j) = tmp_index
    end do
    index(l+1)=index(j)
    index(j)=indext
    jstack=jstack+2
    if (jstack > nstack) then
      write(error_unit, *) 'indexx: nstack too small'
      call exit(99)
    endif
    if (r-i+1 >= j-l) then
      istack(jstack)=r
      istack(jstack-1)=i
      r=j-1
    else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    end if
  end if
end do
end subroutine get_sorted_index

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> christophe cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine coming from numerical recipes icomp_xchg function
!! adapted for this code. the routine will swap the two indexes i and j
!! if the corresponding value of the first (i) is bigger than 
!! the value of the second (j).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine icomp_xchg(arr,i,j)
real(double_precision), dimension(:), intent(in) :: arr !<[in] the reference array
integer, intent(inout) :: i,j !<[in,out] index we will swap if necessary
integer :: swp
if (arr(j) < arr(i)) then
  swp=i
  i=j
  j=swp
end if
end subroutine icomp_xchg

end program nautilus_major_reactions