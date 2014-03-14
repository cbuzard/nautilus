program nautilus_rates

use numerical_types
use iso_fortran_env
use utilities

implicit none

! Locals
character(len=80) :: filename_output
integer :: species, output, reaction ! index for loops
integer :: reactant1, reactant2, reactant3 ! to store indexes of the 3 possible reactants of a given reaction
integer :: error ! to store the state of a read instruction
logical :: isDefined

integer :: nptmax = 1

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
real(double_precision), dimension(:,:), allocatable :: temperature
real(double_precision), dimension(:,:), allocatable :: density
real(double_precision), dimension(:,:), allocatable :: visual_extinction !< visual extinction
real(double_precision), dimension(:,:), allocatable :: zeta
real(double_precision) :: zspace

! For rates
real(double_precision), allocatable, dimension(:,:) :: XK ! (nb_outputs, nb_reactions)
character(len=11), allocatable, dimension(:,:) :: SYMBOL
integer, allocatable, dimension(:) :: NUM !< index of the reactions (one of the columns of the concerned file, 
!! declaring a given number for each reaction, like a hashtag.
integer, allocatable, dimension(:,:) :: reaction_substances

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

write(*,'(a, i5)') 'Spatial resolution: ', nptmax
write(*,'(a, i5)') 'Number of time outputs: ', nb_outputs
write(*,'(a, i5)') 'Number of species: ', nb_species
write(*,'(a, i6)') 'Number of reactions: ', nb_reactions

! We allocate the output arrays
allocate(species_name(nb_species))

allocate(time(nb_outputs))
allocate(temperature(nptmax, nb_outputs))
allocate(density(nptmax, nb_outputs))
allocate(visual_extinction(nptmax, nb_outputs))
allocate(zeta(nptmax, nb_outputs))

allocate(abundances(nb_outputs, nb_species+1, nptmax)) ! We create an extra species that will always have an abundance of 1

allocate(symbol(7,nb_reactions))
allocate(xk(nb_outputs, nb_reactions))
allocate(num(nb_reactions))
allocate(reaction_substances(3, nb_reactions))

! Outputs
allocate(reaction_fluxes(nb_outputs, nb_reactions))

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Reading unformatted outputs...'
! We read output files
do output=1,nb_outputs
  write(filename_output, '(a,i0.6,a)') 'abundances.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) time(output), zspace, species_name(1:nb_species)
  read(10) temperature(1:nptmax, output), density(1:nptmax, output), visual_extinction(1:nptmax, output), zeta(1:nptmax, output)
  read(10) abundances(output,1:nb_species, 1:nptmax)
  close(10)
  
  write(filename_output, '(a,i0.6,a)') 'rates.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) species_name(1:nb_species)
  read(10) symbol(1:7, 1:nb_reactions)
  read(10) xk(output,1:nb_reactions)
  read(10) num(1:nb_reactions)
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'

! For non existing reactants in reactions, we create a new species whose abundance is always 1, so that we can calculate the fluxes 
!! more easily.
abundances(1:nb_outputs, nb_species+1, 1:nptmax) = 1.d0

! We retrieve indexes of each reactants for each reaction
call set_only_reactants()

! We replace blanck species by 'XXX' for the outputs constrains
do reaction=1,nb_reactions
  do species=1,7
    if (symbol(species, reaction).eq.'   ') then
      symbol(species, reaction) = 'XXX'
    endif
  enddo
enddo

! We write fluxes for all reactions and all output times
do output=1, nb_outputs
  do reaction=1, nb_reactions
    reactant1 = reaction_substances(1, reaction)
    reactant2 = reaction_substances(2, reaction)
    reactant3 = reaction_substances(3, reaction)

    reaction_fluxes(output, reaction) = xk(output, reaction) * abundances(output, reactant1, 1) * &
                                        abundances(output, reactant2, 1) * abundances(output, reactant3, 1)
  enddo
enddo

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Writing rates ASCII files...'

write(rate_format, '(a,i5,a)') '(7(a11," "),', nb_outputs, '(es12.4e2),i5)'
write(time_format, '(a,i5,a)') '(73(" "),a11,', nb_outputs, '(es12.4e2))'

! We write ASCII output file
open(10, file='rates.out')
! all 7 species involved ('XXX' if no species) ; Each column is the flux for several output times'
! The first line list time for each column of flux
write(10, time_format) ' Time (yr)', time(1:nb_outputs)
do reaction=1, nb_reactions
  write(10,rate_format) symbol(1:7, reaction), reaction_fluxes(1:nb_outputs, reaction), num(reaction)
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
reaction_substances(1:3, 1:nb_reactions) = no_species

do I=1,nb_reactions
  do J=1,nb_species

    do L=1,3
      if (SYMBOL(L,I).EQ.species_name(J)) then
        reaction_substances(L,I) = J
      endif
    enddo

  enddo
enddo   

return
end subroutine set_only_reactants

end program nautilus_rates