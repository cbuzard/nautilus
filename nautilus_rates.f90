program nautilus_rates

use numerical_types
use iso_fortran_env
use utilities

implicit none

integer, parameter :: MAX_REACTANTS = 3 !< The maximum number of reactants for one reaction. 
!! Do not think that changing this parameter alone is sufficient to allow the code to handle it directly !!
integer, parameter :: MAX_PRODUCTS = 5 !< The maximum number of products for one reaction. 
!! Do not think that changing this parameter alone is sufficient to allow the code to handle it directly !!
integer, parameter :: MAX_COMPOUNDS = MAX_REACTANTS + MAX_PRODUCTS !< Total maximum number of compounds for one reaction (reactants + products)
!! Warning: If this number change, get_jacobian(N, T, Y, J, IAN, JAN, PDJ) must be actualised, since each reactant and product has
!! its own variable, a new one must be created for the new column possible. 

! Locals
character(len=80) :: filename_output
integer :: species, output, reaction ! index for loops
integer :: reactant1, reactant2, reactant3 ! to store indexes of the 3 possible reactants of a given reaction
integer :: error ! to store the state of a read instruction
logical :: isDefined

integer :: nb_sample_1D = 1

character(len=80) :: rate_format, time_format !< string to store specific format used to output datas

character(len=11), dimension(:), allocatable :: species_name !< species name as a string
real(double_precision), dimension(:,:,:), allocatable :: abundances !< abundances over time for each species. (nb_outputs, nb_species)

! Temporary values to store the total number of species and reactions
integer :: nb_species_for_grain !< number of species involved in grain surface reactions
integer :: nb_surface_reactions !< number of reactions on the grain surface
integer :: nb_species_for_gas !< number of species involved in gas phase reactions
integer :: nb_gas_phase_reactions !< number of reactions in gas phase

integer :: nb_outputs !< The total number of outputs
integer :: nb_species !< the total number of species
integer :: nb_reactions !< the total number of reactions

real(double_precision), dimension(:), allocatable :: time
real(double_precision), dimension(:,:), allocatable :: gas_temperature
real(double_precision), dimension(:,:), allocatable :: dust_temperature
real(double_precision), dimension(:,:), allocatable :: density
real(double_precision), dimension(:,:), allocatable :: visual_extinction !< visual extinction
real(double_precision), dimension(:,:), allocatable :: zeta

! For rates
real(double_precision), allocatable, dimension(:,:) :: reaction_rates ! (nb_outputs, nb_reactions)
character(len=11), allocatable, dimension(:,:) :: REACTION_COMPOUNDS_NAMES
integer, allocatable, dimension(:) :: REACTION_ID !< index of the reactions (one of the columns of the concerned file, 
!! declaring a given number for each reaction, like a hashtag.
integer, allocatable, dimension(:,:) :: REACTION_COMPOUNDS_ID

! Output of the code
real(double_precision), allocatable, dimension(:,:) :: reaction_fluxes

! We calculate the total number of outputs by checking for each file if it exist or not.
nb_outputs = 0
isDefined = .true.
do while(isDefined)
  nb_outputs = nb_outputs + 1
  write(filename_output, '(a,i0.6,a)') 'abundances.',nb_outputs,'.out'
  inquire(file=filename_output, exist=isDefined)

enddo
nb_outputs = nb_outputs - 1

! We get the number of reactions and species
call get_linenumber(filename='gas_species.in', nb_lines=nb_species_for_gas)
call get_linenumber(filename='gas_reactions.in', nb_lines=nb_gas_phase_reactions)

call get_linenumber(filename='grain_species.in', nb_lines=nb_species_for_grain)
call get_linenumber(filename='grain_reactions.in', nb_lines=nb_surface_reactions)

nb_species = nb_species_for_gas + nb_species_for_grain ! The total number of species, sum of species in gas and grain
nb_reactions = nb_gas_phase_reactions + nb_surface_reactions ! The total number of reactions, sum of species in gas and grain

write(*,'(a, i5)') 'Spatial resolution: ', nb_sample_1D
write(*,'(a, i5)') 'Number of time outputs: ', nb_outputs
write(*,'(a, i5)') 'Number of species: ', nb_species
write(*,'(a, i6)') 'Number of reactions: ', nb_reactions

! We allocate the output arrays
allocate(species_name(nb_species))

allocate(time(nb_outputs))
allocate(gas_temperature(nb_sample_1D, nb_outputs))
allocate(dust_temperature(nb_sample_1D, nb_outputs))
allocate(density(nb_sample_1D, nb_outputs))
allocate(visual_extinction(nb_sample_1D, nb_outputs))
allocate(zeta(nb_sample_1D, nb_outputs))

allocate(abundances(nb_outputs, nb_species+1, nb_sample_1D)) ! We create an extra species that will always have an abundance of 1

allocate(REACTION_COMPOUNDS_NAMES(MAX_COMPOUNDS,nb_reactions))
allocate(reaction_rates(nb_outputs, nb_reactions))
allocate(REACTION_ID(nb_reactions))
allocate(REACTION_COMPOUNDS_ID(3, nb_reactions))

! Outputs
allocate(reaction_fluxes(nb_outputs, nb_reactions))

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Reading unformatted outputs...'
! We read output files
do output=1,nb_outputs
  write(filename_output, '(a,i0.6,a)') 'abundances.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) time(output)
  read(10) gas_temperature(1:nb_sample_1D, output), dust_temperature(1:nb_sample_1D, output), density(1:nb_sample_1D, output), &
           visual_extinction(1:nb_sample_1D, output), zeta(1:nb_sample_1D, output)
  read(10) abundances(output,1:nb_species, 1:nb_sample_1D)
  close(10)
  
  write(filename_output, '(a,i0.6,a)') 'rates.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) species_name(1:nb_species)
  read(10) REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, 1:nb_reactions)
  read(10) reaction_rates(output,1:nb_reactions)
  read(10) REACTION_ID(1:nb_reactions)
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'

! For non existing reactants in reactions, we create a new species whose abundance is always 1, so that we can calculate the fluxes 
!! more easily.
abundances(1:nb_outputs, nb_species+1, 1:nb_sample_1D) = 1.d0

! We retrieve indexes of each reactants for each reaction
call set_only_reactants()

! We replace blanck species by 'XXX' for the outputs constrains
do reaction=1,nb_reactions
  do species=1,MAX_COMPOUNDS
    if (REACTION_COMPOUNDS_NAMES(species, reaction).eq.'   ') then
      REACTION_COMPOUNDS_NAMES(species, reaction) = 'XXX'
    endif
  enddo
enddo

! We write fluxes for all reactions and all output times
do output=1, nb_outputs
  do reaction=1, nb_reactions
    reactant1 = REACTION_COMPOUNDS_ID(1, reaction)
    reactant2 = REACTION_COMPOUNDS_ID(2, reaction)
    reactant3 = REACTION_COMPOUNDS_ID(3, reaction)

    reaction_fluxes(output, reaction) = reaction_rates(output, reaction) * abundances(output, reactant1, 1) * &
                                        abundances(output, reactant2, 1) * abundances(output, reactant3, 1)
  enddo
enddo

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Writing rates ASCII files...'

write(rate_format, '(a,i0,a,i0,a)') '(',MAX_COMPOUNDS,'(a11," "),', nb_outputs, '(es12.4e2),i5)'
write(time_format, '(a,i0,a,i0,a)') '(',12*(MAX_COMPOUNDS-1),'(" "),a12,', nb_outputs, '(es12.4e2))'

! We write ASCII output file
open(10, file='rates.out')
! all MAX_COMPOUNDS species involved ('XXX' if no species) ; Each column is the flux for several output times'
! The first line list time for each column of flux
write(10, time_format) 'Time (year)', time(1:nb_outputs)
do reaction=1, nb_reactions
  write(10,rate_format) REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_fluxes(1:nb_outputs, reaction), &
                        REACTION_ID(reaction)
enddo
close(10)

! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Writing rates ASCII files... Done'

contains 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Initialize reactants indexes for all reactions. Usefull to retrieve abundances quickly
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine set_only_reactants()

implicit none

! Locals
integer :: I,J,L

integer :: no_species

no_species = nb_species + 1

! By default, non existing reactants (dummy species) will be assigned (nb_species+1)
REACTION_COMPOUNDS_ID(1:3, 1:nb_reactions) = no_species

do I=1,nb_reactions
  do J=1,nb_species

    do L=1,MAX_REACTANTS
      if (REACTION_COMPOUNDS_NAMES(L,I).EQ.species_name(J)) then
        REACTION_COMPOUNDS_ID(L,I) = J
      endif
    enddo

  enddo
enddo   

return
end subroutine set_only_reactants

end program nautilus_rates