program nautilus_outputs

use numerical_types
use iso_fortran_env
use utilities

implicit none

! Locals
character(len=80) :: filename_output
integer :: species, output, i, idx_1D ! index for loops
integer :: error ! to store the state of a read instruction
logical :: isDefined

character(len=80) :: filename !< name of the file to be read
character(len=80) :: output_format !< format used to output data
character(len=200) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string

real(double_precision), parameter :: YEAR = 3.15576d7 !< one year in seconds

integer :: nb_sample_1D = 1 !< 1D resolution of the structure
integer :: nb_outputs !< The total number of outputs
integer :: nb_species !< the total number of species

character(len=11), dimension(:), allocatable :: species_name !< species name as a string
real(double_precision), dimension(:,:,:), allocatable :: abundances !< abundances over time for each species. (nb_outputs, nb_species)

integer :: nb_species_for_gas, nb_species_for_grain !< temporary values to store the total number of species

character(len=11), dimension(:), allocatable :: gas_species_label 
character(len=11), dimension(:), allocatable :: surface_species_label 

real(double_precision), dimension(:), allocatable :: time !< Simulation time [s]
real(double_precision), dimension(:,:), allocatable :: gas_temperature !< [K]
real(double_precision), dimension(:,:), allocatable :: dust_temperature !< [K]
real(double_precision), dimension(:,:), allocatable :: density !< [part/cm^3] 
real(double_precision), dimension(:,:), allocatable :: visual_extinction !< visual extinction [mag]
real(double_precision), dimension(:,:), allocatable :: x_rate !< X ionisation rate [s-1]

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
allocate(gas_species_label(nb_species_for_gas))
allocate(surface_species_label(nb_species_for_grain))

nb_species = nb_species_for_gas + nb_species_for_grain ! The total number of species, sum of species in gas and grain

write(*,'(a, i5)') 'Spatial resolution: ', nb_sample_1D
write(*,'(a, i5)') 'Number of time outputs: ', nb_outputs
write(*,'(a, i5)') 'Number of species: ', nb_species

! We allocate the output arrays
allocate(species_name(nb_species))

allocate(time(nb_outputs))
allocate(gas_temperature(nb_sample_1D, nb_outputs))
allocate(dust_temperature(nb_sample_1D, nb_outputs))
allocate(density(nb_sample_1D, nb_outputs))
allocate(visual_extinction(nb_sample_1D, nb_outputs))
allocate(x_rate(nb_sample_1D, nb_outputs))

allocate(abundances(nb_outputs, nb_species, nb_sample_1D))

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Reading species name...'


! Reading list of species for gas phase
filename = 'gas_species.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

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
      read(line, '(a11,i3,13(I3))')  gas_species_label(I)
    
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if


! Reading list of species for grain surface
filename = 'grain_species.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

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
      read(line, '(a11)')  surface_species_label(I)
    
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if


! putting everything back into the big tables
do I=1,nb_species_for_gas 
  species_name(I) = gas_species_label(I)
enddo
do I=1,nb_species_for_grain 
  species_name(nb_species_for_gas+I) = surface_species_label(I)
enddo

! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading species name... Done'

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Reading unformatted outputs...'
! We read output files
do output=1,nb_outputs
  write(filename_output, '(a,i0.6,a)') 'abundances.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) time(output)
  read(10) gas_temperature(1:nb_sample_1D, output), dust_temperature(1:nb_sample_1D, output), density(1:nb_sample_1D, output), &
           visual_extinction(1:nb_sample_1D, output), x_rate(1:nb_sample_1D, output)
  read(10) abundances(output,1:nb_species, 1:nb_sample_1D)
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'

! Test if the folder exists
inquire(file='ab', exist=isDefined)

! We create the folder 'ab' if he doesn't exists.
if (.not.isDefined) then
  call system("mkdir ab")
end if

! Remove all existing *.ab if needed. Will return a warning in standard output if nothing exists
call system("rm ab/*.ab")

!####################################################@@
! This part is to write one file per species, each line being one output time
!####################################################@@

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Writing *.ab ASCII files...'
! We write ASCII output file, one file per species
do species=1, nb_species
  write(filename_output, '(a,a,a)') 'ab/', trim(species_name(species)), '.ab'
  open(10, file=filename_output)
  write(10,'(a)') '! time [year]; Each column is the abundance (relative to H) [number ratio] for several spatial positions'
  
  write(output_format, *) '(es10.3e2,',nb_sample_1D,'(es13.6e2," "))'
  do output=1, nb_outputs
    write(10,output_format) time(output)/YEAR, abundances(output, species, 1:nb_sample_1D)
  enddo
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Writing output files... Done'

!####################################################@@
! This part is to write one file per output time, each line being one species
!####################################################@@

!~ ! The next write will be written in the same line
!~ write(*,'(a)', advance='no') 'Writing output ASCII files...'
!~ ! We write ASCII output file, one file per output
!~ do output=1, nb_outputs
!~   write(filename_output, '(a,i0.5,a)') 'ab/abundances.', output, '.ab'
!~   open(10, file=filename_output)
!~   write(10,'(a,es10.2e2, a)') '! time =', time(output) / YEAR, ' years'
!~   write(10,'(a)') '! species name ; Abundance'
!~   do species=1, nb_species
!~     write(10,*) species_name(species), abundances(output, species, 1:nb_sample_1D)
!~   enddo
!~   close(10)
!~ enddo
!~ ! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
!~ write(*,'(a,a)') achar(13), 'Writing output files... Done'

! Test if the folder exists
inquire(file='struct', exist=isDefined)

! We create the folder 'ab' if he doesn't exists.
if (.not.isDefined) then
  call system("mkdir struct")
end if

! Remove all existing *.ab if needed. Will return a warning in standard output if nothing exists
call system("rm struct/*.struct")

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Writing *.struct ASCII files...'
! We write ASCII output file, one file per species
do idx_1D=1, nb_sample_1D
  write(filename_output, '(a,i0.5,a)') 'struct/output.', idx_1D, '.struct'
  open(10, file=filename_output)
  write(10,'(a)') '! time     ; gas temperature ; dust temperature&
                   & ; log10(H2 density)  ; visual extinction ; x ionization rate'
  write(10,'(a)') '!  [year]  ;       [K]       ;         [K]     &
                   & ; [log10(part/cm^3)] ;           [mag]   ;       [s-1] '
  
  do output=1, nb_outputs
    write(10,'(6(es10.3e2," "))') time(output)/YEAR, gas_temperature(idx_1D, output), dust_temperature(idx_1D, output), &
           density(idx_1D, output), visual_extinction(idx_1D, output), x_rate(idx_1D, output)
  enddo
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Writing structure output files... Done'

end program nautilus_outputs