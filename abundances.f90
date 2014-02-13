program abundances

use numerical_types
use iso_fortran_env

implicit none

! Locals
character(len=80) :: filename_output
integer :: i ! index for loops
integer :: error ! to store the state of a read instruction
logical :: isDefined

integer :: nb_outputs !< The total number of outputs
integer :: nb_species !< the total number of species

character(len=11), dimension(:), allocatable :: species_name !< species name as a string
real(double_precision), dimension(:,:), allocatable :: abundances !< abundances over time for each species. (nb_outputs, nb_species)

integer :: nb_species_for_gas, nb_species_for_grain !< temporary values to store the total number of species

real(double_precision), dimension(:), allocatable :: time
real(double_precision), dimension(:), allocatable :: temperature
real(double_precision), dimension(:), allocatable :: density
real(double_precision), dimension(:), allocatable :: visual_extinction !< visual extinction
real(double_precision), dimension(:), allocatable :: zeta
real(double_precision) :: zspace

!~ INTEGER,PARAMETER :: ntime = 9
!~ INTEGER,PARAMETER :: nsmax = 684
!~ CHARACTER (len=30) :: filename
!~ CHARACTER (len=30) :: file_output
!~ CHARACTER (len=11), DIMENSION(ntime,nsmax) :: spec
!~ REAL(KIND=8), DIMENSION(ntime,nsmax) :: ab
!~ REAL(KIND=8), DIMENSION(ntime) :: time, temp,dens, tau, zeta
!~ REAL(KIND=8) :: zspace

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

write(*,*) nb_outputs, nb_species

! We allocate the output arrays
allocate(species_name(nb_species))

allocate(time(nb_outputs))
allocate(temperature(nb_outputs))
allocate(density(nb_outputs))
allocate(visual_extinction(nb_outputs))
allocate(zeta(nb_outputs))

allocate(abundances(nb_outputs, nb_species))

! We read output files
do i=1,nb_outputs
  write(filename_output, '(a,i0.6,a)') 'abundances.',i,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) time(i), zspace, species_name(1:nb_species)
  read(10) temperature(i), density(i), visual_extinction(i), zeta(i)
  read(10) abundances(i,1:nb_species)
  close(10)
enddo

! We write ASCII output file
open(10, file='abundances.00001.ascii')
write(10,'(a)') '! Species name ; Each column is the abundance for several times. The first line is the list of times'
write(10,*) 'Output times ', time(1:nb_outputs)
do i=1, nb_species
  write(10,*) species_name(i), abundances(1:nb_outputs, i)
enddo
close(10)

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
  character(len=80) :: line
  character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
  integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string

  !------------------------------------------------------------------------------
  nb_lines = 0
  
  ! Read in filenames and check for duplicate filenames
  inquire (file=filename, exist=test)
  if (.not.test) then
    write(Error_Unit,'(a,a,a)') 'Error: the file "',trim(filename),'" does not exist.'
    call exit(1)
  end if
  open(15, file=filename, status='old')
  
  error = 0
  do 
    read(15,'(a)',iostat=error) line
    
    if(error.ne.0) then
      exit
    endif
    
    ! We get only what is on the left of an eventual comment parameter
    comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    
    
    if (trim(line).ne.'') then
      nb_lines = nb_lines + 1
    endif
  enddo
  close(15)
  
end subroutine get_linenumber

end program abundances