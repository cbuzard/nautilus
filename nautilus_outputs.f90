program nautilus_outputs

use numerical_types
use iso_fortran_env
use utilities

implicit none

! Locals
character(len=80) :: filename_output
integer :: species, output ! index for loops
integer :: error ! to store the state of a read instruction
logical :: isDefined

integer :: nptmax = 1
integer :: nb_outputs !< The total number of outputs
integer :: nb_species !< the total number of species

character(len=11), dimension(:), allocatable :: species_name !< species name as a string
real(double_precision), dimension(:,:,:), allocatable :: abundances !< abundances over time for each species. (nb_outputs, nb_species)

integer :: nb_species_for_gas, nb_species_for_grain !< temporary values to store the total number of species

real(double_precision), dimension(:), allocatable :: time
real(double_precision), dimension(:,:), allocatable :: temperature
real(double_precision), dimension(:,:), allocatable :: density
real(double_precision), dimension(:,:), allocatable :: visual_extinction !< visual extinction
real(double_precision), dimension(:,:), allocatable :: zeta
real(double_precision) :: zspace

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
call get_linenumber(filename='grain_species.in', nb_lines=nb_species_for_grain)

nb_species = nb_species_for_gas + nb_species_for_grain ! The total number of species, sum of species in gas and grain

write(*,'(a, i5)') 'Spatial resolution: ', nptmax
write(*,'(a, i5)') 'Number of time outputs: ', nb_outputs
write(*,'(a, i5)') 'Number of species: ', nb_species

! We allocate the output arrays
allocate(species_name(nb_species))

allocate(time(nb_outputs))
allocate(temperature(nptmax, nb_outputs))
allocate(density(nptmax, nb_outputs))
allocate(visual_extinction(nptmax, nb_outputs))
allocate(zeta(nptmax, nb_outputs))

allocate(abundances(nb_outputs, nb_species, nptmax))


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
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Writing *.ab ASCII files...'
! We write ASCII output file
do species=1, nb_species
  write(filename_output, '(a,a)') trim(species_name(species)), '.ab'
  open(10, file=filename_output)
  write(10,'(a)') '! time ; Each column is the abundance for several spatial positions'
  do output=1, nb_outputs
    write(10,*) time(output), abundances(output, species, 1:nptmax)
  enddo
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Writing *.ab ASCII files... Done'


end program nautilus_outputs