!******************************************************************************
! MODULE: structure
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines linked to the physical structure
!! (cloud and such) and its properties
!
!******************************************************************************

module structure

use iso_fortran_env
use numerical_types
use global_variables

implicit none

integer :: structure_sample !< the number of points from structure_evolution.dat, about the time evolution of the physical structure properties

real(double_precision), allocatable, dimension(:) :: structure_time !< dim(structure_sample) time [s] read from structure_evolution.dat
real(double_precision) :: structure_sample_step !< time [s] between each sample point for the structure evolution

real(double_precision), allocatable, dimension(:) :: structure_log_Av !< dim(structure_sample) log10(Av) [log10(mag)] read from structure_evolution.dat
real(double_precision), allocatable, dimension(:) :: structure_log_density !< dim(structure_sample) log10(density) [log10(part/cm^3)] read from structure_evolution.dat
real(double_precision), allocatable, dimension(:) :: structure_log_gas_temperature !< dim(structure_sample) log10(gas temperature) [log10(K)] read from structure_evolution.dat

logical :: read_dust !< If dust temperature must be read directly from the data file (or not)
real(double_precision), allocatable, dimension(:) :: structure_log_dust_temperature !< dim(structure_sample) log10(dust temperature) [log10(K)] read from structure_evolution.dat


contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read time evolution for the physical structure in structure_evolution.dat
!! to initialize arrays. Thus, interpolation of physical structure properties
!! throughout the simulation be possible.
!! \n
!! structure_evolution.dat will have a two lines header. Then columns will be as follow :
!! \n Time (Myr) ; Number density (part/cm3) ; Temperature (K) ; Av (mag)
!
!> @warning Time sample MUST be linearly and equally spaced
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine init_structure_evolution()

use utilities

implicit none


character(len=80) :: filename = 'structure_evolution.dat' !< name of the file in which time evolution is stored
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction
integer :: nb_columns

integer :: i !< loop index
logical :: isDefined
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then
  
  call get_linenumber(filename=filename, nb_lines=structure_sample)
  nb_columns = get_nb_columns(filename)
  
  if (nb_columns.ge.5) then
    read_dust = .true.
  else
    read_dust = .false.
  endif
  
  ! If by any chance the routine is run several times, this test is here to avoid bugs
  if (allocated(structure_time)) then
    deallocate(structure_time)
    deallocate(structure_log_density)
    deallocate(structure_log_gas_temperature)
    deallocate(structure_log_Av)
  end if
  allocate(structure_time(structure_sample))
  allocate(structure_log_density(structure_sample))
  allocate(structure_log_gas_temperature(structure_sample))
  allocate(structure_log_Av(structure_sample))
  
  
  structure_time(1:structure_sample) = 0.d0
  structure_log_density(1:structure_sample) = 0.d0
  structure_log_gas_temperature(1:structure_sample) = 0.d0
  structure_log_Av(1:structure_sample) = 0.d0
  
  if (read_dust) then
    if (GRAIN_TEMPERATURE_TYPE.ne.'table') then
      write(Error_unit, *) 'Error: The grain temperature column exist'
      write(Error_unit, *) 'in structure_evolution.dat but '
      write(Error_unit, *) 'GRAIN_TEMPERATURE_TYPE = ',trim(GRAIN_TEMPERATURE_TYPE)
      call exit(10)
    endif
    
    if (allocated(structure_log_dust_temperature)) then
      deallocate(structure_log_dust_temperature)
    endif
    
    allocate(structure_log_dust_temperature(structure_sample))
    structure_log_dust_temperature(1:structure_sample) = 0.d0
  endif
  
  open(10, file=filename, status='old')
  i = 1
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
      if (read_dust) then
        read(line, *) structure_time(i), structure_log_Av(i), structure_log_density(i), structure_log_gas_temperature(i), &
                      structure_log_dust_temperature(i)
      else
        read(line, *) structure_time(i), structure_log_Av(i), structure_log_density(i), structure_log_gas_temperature(i)
      endif
      i = i + 1
    endif
  enddo
  close(10)
  
  ! We convert time in million year to time in seconds
  structure_time(1:structure_sample) = (1.d6 * YEAR) * structure_time(1:structure_sample)
  
  ! We get the space between all times for the structure evolution. This has some sense only if times are linearly and equally spaced.
  structure_sample_step = (structure_time(structure_sample) - structure_time(1)) / dfloat(structure_sample - 1)
  
else
  write (Error_Unit,*) 'Error: The file "structure_evolution.dat" does not exist.'
  call exit(7)
end if

return
end subroutine init_structure_evolution

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get properties of the structure (Av, temperature, density)
!! in function of the time. We must launch 'init_structure_evolution' 
!! beforehand to initialize arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_structure_properties_table(time, Av, density, gas_temperature, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]
  
  ! Local
  integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine. 
  real(double_precision) :: t1, t2, log_y1, log_y2
  real(double_precision) :: interpolation_tmp ! Tmp value that contain a small part of the linear interpolation, to avoid computing it three times
  !------------------------------------------------------------------------------

  if (time .lt. structure_time(structure_sample-1)) then
    
    ! in the range
    closest_low_id = 1 + int((time - structure_time(1)) / structure_sample_step)
    
    ! We get the closest times (left and right) from the current one
    t1 = structure_time(closest_low_id)
    t2 = structure_time(closest_low_id + 1)
    
    interpolation_tmp = (time - t2) / (t1 - t2)
    
    ! For the density
    log_y1 = structure_log_density(closest_low_id)
    log_y2 = structure_log_density(closest_low_id + 1)
    density = 10.0d0**(log_y2 + (log_y1 - log_y2) * interpolation_tmp)
    
    ! For the visual extinction
    log_y1 = structure_log_Av(closest_low_id)
    log_y2 = structure_log_Av(closest_low_id + 1)
    Av = 10.0d0**(log_y2 + (log_y1 - log_y2) * interpolation_tmp)
    
    ! For the gas temperature
    log_y1 = structure_log_gas_temperature(closest_low_id)
    log_y2 = structure_log_gas_temperature(closest_low_id + 1)
    gas_temperature = 10.0d0**(log_y2 + (log_y1 - log_y2) * interpolation_tmp)
    
  else if (time .ge. structure_time(structure_sample-1)) then
    density = 10.0d0**(structure_log_density(structure_sample))
    Av = 10.0d0**(structure_log_Av(structure_sample))
    gas_temperature = 10.0d0**(structure_log_gas_temperature(structure_sample))
  end if
  
  call get_grain_temperature(gas_temperature=gas_temperature, grain_temperature=grain_temperature)

  return
end subroutine get_structure_properties_table

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get properties of the structure (Av, temperature, density)
!! In this routine, everything is constant.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_structure_properties_fixed(time, Av, density, gas_temperature, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]

  !------------------------------------------------------------------------------

  density = initial_gas_density
  Av = INITIAL_VISUAL_EXTINCTION
  gas_temperature = initial_gas_temperature
  
  call get_grain_temperature(gas_temperature=gas_temperature, grain_temperature=grain_temperature)

  return
end subroutine get_structure_properties_fixed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get the grain temperature, here assumed to be equal to 
!! the gas temperature
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature_gas(gas_temperature, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  
  ! Outputs
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  grain_temperature = gas_temperature
  
  return
end subroutine get_grain_temperature_gas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get the grain temperature, here fixed to the initial 
!! dust temperature
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature(gas_temperature, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  
  ! Outputs
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  grain_temperature = initial_dust_temperature
  
  return
end subroutine get_grain_temperature
end module structure