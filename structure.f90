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
subroutine get_structure_properties_table(time, Av, density, gas_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]
  
  ! Local
  integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine. 
  real(double_precision) :: t1, t2, log_y1, log_y2
  real(double_precision) :: interpolation_tmp ! Tmp value that contain a small part of the linear interpolation, to avoid computing it three times
  !------------------------------------------------------------------------------

  if (time .lt. structure_time(structure_sample)) then
    
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
    
  else if (time .ge. structure_time(structure_sample)) then
    density = 10.0d0**(structure_log_density(structure_sample))
    Av = 10.0d0**(structure_log_Av(structure_sample))
    gas_temperature = 10.0d0**(structure_log_gas_temperature(structure_sample))
  end if
  
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
subroutine get_structure_properties_fixed(time, Av, density, gas_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]

  !------------------------------------------------------------------------------

  density = initial_gas_density
  Av = INITIAL_VISUAL_EXTINCTION
  gas_temperature = initial_gas_temperature
  
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
subroutine get_grain_temperature_gas(time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
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
subroutine get_grain_temperature_fixed(time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
  ! Outputs
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  grain_temperature = initial_dust_temperature
  
  return
end subroutine get_grain_temperature_fixed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get the grain temperature, here interpolated from 
!! the data read in structure_evolution.dat
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature_table(time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
  ! Outputs
  real(double_precision), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  ! Local
  integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine. 
  real(double_precision) :: t1, t2, log_y1, log_y2
  real(double_precision) :: interpolation_tmp ! Tmp value that contain a small part of the linear interpolation, to avoid computing it three times
  !------------------------------------------------------------------------------

  if (time .lt. structure_time(structure_sample)) then
    
    ! in the range
    closest_low_id = 1 + int((time - structure_time(1)) / structure_sample_step)
    
    ! We get the closest times (left and right) from the current one
    t1 = structure_time(closest_low_id)
    t2 = structure_time(closest_low_id + 1)
    
    interpolation_tmp = (time - t2) / (t1 - t2)
    
    ! For the grain temperature
    log_y1 = structure_log_dust_temperature(closest_low_id)
    log_y2 = structure_log_dust_temperature(closest_low_id + 1)
    grain_temperature = 10.0d0**(log_y2 + (log_y1 - log_y2) * interpolation_tmp)
    
  else if (time .ge. structure_time(structure_sample)) then
    grain_temperature = 10.0d0**(structure_log_dust_temperature(structure_sample))
  end if  
  return
end subroutine get_grain_temperature_table

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the diffusion routine for 0D structure. Thus, 
!! nothing is done here, just to point toward a "do nothing" procedure
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine structure_diffusion_0D(timestep, temp_abundances)
  
  implicit none
  ! Inputs
  real(double_precision), intent(in) :: timestep !<[in] timestep for the diffusion process [s]
  
  ! Inputs/Outputs
  real(double_precision), dimension(:,:), intent(inout) :: temp_abundances !<[in,out] dim(nb_species, nb_sample_1D) 
  !! The abundances for all species, and 
  !! all 1D mesh points (relative to H) [number ratio]
  
  ! Abundances stay equal here
  
end subroutine structure_diffusion_0D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant & Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief 1D diffusion using a Crank Nicholson Scheme. Routine inspired by
!! Numerical Recipes. In particular the tridiag routine to solve tridiagonal
!! matrices
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine crank_nicholson_1D(f,ny,dt,dy,nu,rho, ibc)
! A Crank-Nicholson scheme
! ibc is a flag for boundary conditions
! ibc = 0 -> no flux boundaries (bc is not used then)
! ibc = 1 -> user supplied boundary conditions (bc)
! ibc > 1 -> stops the code, insulting the user
! Computes the evolution for a single timestep of the equation
! df/dt=a d/dy b d/dy c f + S
! On an homogeneous mesh, for user specified a,b,c
! The discretized equation takes the form
! xf(i+1)+yf(i)+zf(i-1)=W
!
! Tu ne dois pas avoir de a, b et c dans ta version. Si je te dis pas de
! betise, avec les notations de ce header :
! a = 1/rho
! b = nu * rho
! c = 1.
! S = 0.
implicit none

! Inputs
integer, intent(in) :: ny
real(double_precision), intent(in), dimension(0:ny) :: rho
real(double_precision), intent(in), dimension(0:ny) :: nu
real(double_precision), intent(in) :: dt
real(double_precision), intent(in) :: dy
integer, intent(in) :: ibc

! Outputs
real(double_precision), intent(out), dimension(0:ny) :: f

! Locals
integer :: ind
real(double_precision), dimension(0:ny) :: s,Q,W,x,y,z,u,v, dd1d
real(double_precision) :: d !, nu


dd1d(:)=nu(:)*rho(:)

d=dt/(dy**2)
Q(:)=rho(:)

s(:)=0.d0

do ind = 1, ny-1
  W(ind) = s(ind)*dt + d/4*(dd1d(ind+1)+dd1d(ind))*f(ind+1)/rho(ind) + (Q(ind)-d/4*(dd1d(ind+1)+2*dd1d(ind) &
  +dd1d(ind-1)))*f(ind)/rho(ind)+d/4*(dd1d(ind)+dd1d(ind-1))*f(ind-1)/rho(ind)
enddo

do ind = 1,ny-1
  x(ind) = -d/4*(dd1d(ind+1)+dd1d(ind))/rho(ind)
  y(ind) = Q(ind)/rho(ind) + d/4*(dd1d(ind+1)+2*dd1d(ind)+dd1d(ind-1))/rho(ind)
  z(ind) = -d/4*(dd1d(ind)+dd1d(ind-1))/rho(ind)
enddo

! Test
u(ny)=1.d0
v(ny)=0.d0
x(ny) = -d/2*dd1d(ny)/rho(ny)
y(ny) = Q(ny)/rho(ny) + d/4*(3*dd1d(ny)+dd1d(ny-1))/rho(ny)
z(ny) = -d/4*(dd1d(ny)+dd1d(ny-1))/rho(ny)
W(ny) = d/2*dd1d(ny)*f(ny)/rho(ny) + (Q(ny)-d/4*(3*dd1d(ny) &
+dd1d(ny-1)))*f(ny)/rho(ny)+d/4*(dd1d(ny)+dd1d(ny-1))*f(ny-1)/rho(ny)

do ind = ny, 1, -1
  u(ind-1) = -z(ind)/ (x(ind)*u(ind) + y(ind))
  v(ind-1) = ( W(ind) - x(ind)*v(ind) )/( x(ind)*u(ind) + y(ind) )
enddo

if (ibc.eq.0) f(0) =  v(0)/( 1. - u(0) )

do ind = 0, ny-1
  f(ind+1) = u(ind)*f(ind) + v(ind)
enddo

return
end subroutine crank_nicholson_1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the diffusion routine for 1D structure of a disk object. 
!! The diffusion is done vertically, on the 'z' dimension. The first point is on the outside 
!! TODO (check that once the procedure is written)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine structure_diffusion_1D_disk_z(timestep, temp_abundances)
  
  implicit none
  ! Inputs
  real(double_precision), intent(in) :: timestep !<[in] timestep for the diffusion process [s]
  
  ! Inputs/Outputs
  real(double_precision), dimension(:,:), intent(inout) :: temp_abundances !<[in,out] dim(nb_species, nb_sample_1D) 
  !! The abundances for all species, and 
  !! all 1D mesh points (relative to H) [number ratio]
  
  ! Locals
  integer :: reaction !< For loops
  
  do reaction=1,nb_species
    call crank_nicholson_1D(f=temp_abundances(reaction, 1:nb_sample_1D), ny=nb_sample_1D-1, dt=timestep, dy=grid_cell_size, &
    nu=diffusion_coefficient(1:nb_sample_1D), rho=H_number_density(1:nb_sample_1D), ibc=0)
  enddo  
  
end subroutine structure_diffusion_1D_disk_z

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the timestep routine for 0D structure.
!! Nothing is done here, there will be only one timestep, equal to the output one. 
!! Indeed, subtimestep are defined by the diffusion process, and there's none in 0D
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_timestep_0D(current_time, final_time, next_timestep)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: current_time !<[in] current time [s]
  real(double_precision), intent(in) :: final_time !<[in] Final output time of the current 
  !! loop occurence. The last sub-step must lead exactly to this time [s]

  ! Outputs
  real(double_precision), intent(out) :: next_timestep !<[out] The next integration sub timestep withing an output integration step [s]

  next_timestep = final_time - current_time

end subroutine get_timestep_0D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the timestep routine for 1D structure of a disk object 
!! whose 1D diffusion is on the vertical dimension. 
!! This routine provide a sub-timestep that allow accurate computation of the diffusion.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_timestep_1D_disk_z(current_time, final_time, next_timestep)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: current_time !<[in] current time [s]
  real(double_precision), intent(in) :: final_time !<[in] Final output time of the current 
  !! loop occurence. The last sub-step must lead exactly to this time [s]

  ! Outputs
  real(double_precision), intent(out) :: next_timestep !<[out] The next integration sub timestep withing an output integration step [s]

  next_timestep = grid_cell_size**2 / maxval(diffusion_coefficient)

  ! TODO write the calculation of the diffusion timestep before that test
  if (current_time+next_timestep.gt.final_time) then
    next_timestep = final_time - current_time
  endif
end subroutine get_timestep_1D_disk_z

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read 1D evolution in 1D_evolution.dat
!! to initialize arrays. Thus, interpolation of physical structure properties
!! throughout the simulation be possible.
!! \n
!! structure_evolution.dat will have a two lines header. Then columns will be as follow :
!! \n Z position (AU) ; Number density (part/cm3) ; Temperature (K) ; diffusion coeff (cm^2/s) ; Av (mag)
!
!> @warning Time sample MUST be linearly and equally spaced
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine init_1D_evolution()

use utilities

implicit none


character(len=80) :: filename = '1D_evolution.dat' !< name of the file in which time evolution is stored
character(len=200) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: error !< to store the state of a read instruction

real(double_precision), dimension(:), allocatable :: tmp_grid !< Distance [cm] /!\ But read column is in AU
real(double_precision), dimension(:), allocatable :: tmp_density !< Gas density [g/cm^3]
real(double_precision), dimension(:), allocatable :: tmp_gas_temperature !< Gas temperature [K]
real(double_precision), dimension(:), allocatable :: tmp_av !< Visual extinction [mag]
real(double_precision), dimension(:), allocatable :: tmp_kappa !< Diffusion coefficient [cm^2/s]
integer :: closest_low_id, nb_values
real(double_precision) :: x1, x2, y1, y2

integer :: i !< loop index
logical :: isDefined
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then
  
  ! We get the total lines of the file
  open(10, file=filename, status='old')
  i = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
    
    if (line(1:1).ne.comment_character) then
      i = i + 1
    end if
  end do
  close(10)
  
  ! We define the sizes of the arrays
  nb_values = i
  if (allocated(tmp_grid)) then
    deallocate(tmp_grid)
    deallocate(tmp_density)
    deallocate(tmp_gas_temperature)
    deallocate(tmp_av)
    deallocate(tmp_kappa)
  end if
  allocate(tmp_grid(nb_values))
  allocate(tmp_density(nb_values))
  allocate(tmp_gas_temperature(nb_values))
  allocate(tmp_av(nb_values))
  allocate(tmp_kappa(nb_values))
  
  ! We get the values of the torque profile in the file
  open(10, file=filename, status='old')
  i = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
    
    if(line(1:1).ne.comment_character) then
      i = i + 1
      read(line, *, iostat=error) tmp_grid(i), tmp_density(i), tmp_gas_temperature(i), tmp_av(i), tmp_kappa(i)
    end if
  end do
  
  ! Convert distances from AU to cm
  tmp_grid(1:nb_values) = tmp_grid(1:nb_values) * AU
  
  ! We now want to interpolate the values from the input sampling to the desired 1D sampling that is 
  !! defined solely by z_max and nb_sample_1D
  closest_low_id = 1
  do i=1,nb_sample_1D
    
    if ((grid_sample(i) .ge. tmp_grid(1)) .and. (grid_sample(i) .lt. tmp_grid(nb_values))) then
      ! we do not initialize closest_low_id at each step, because the sample is sorted, 
      ! so we know that the id will at least be the one of the previous timestep
      do while (grid_sample(i).gt.tmp_grid(closest_low_id+1))
        closest_low_id = closest_low_id + 1
      end do
      
      x1 = tmp_grid(closest_low_id)
      x2 = tmp_grid(closest_low_id + 1)
      
      ! density
      y1 = tmp_density(closest_low_id)
      y2 = tmp_density(closest_low_id + 1)
      H_number_density(i) = (y2 + (y1 - y2) * (grid_sample(i) - x2) / (x1 - x2))
      
      ! Gas temperature
      y1 = tmp_gas_temperature(closest_low_id)
      y2 = tmp_gas_temperature(closest_low_id + 1)
      gas_temperature(i) = (y2 + (y1 - y2) * (grid_sample(i) - x2) / (x1 - x2))
      
      ! Visual Extinction
      y1 = tmp_av(closest_low_id)
      y2 = tmp_av(closest_low_id + 1)
      visual_extinction(i) = (y2 + (y1 - y2) * (grid_sample(i) - x2) / (x1 - x2))
      
      ! Diffusion coefficient
      y1 = tmp_kappa(closest_low_id)
      y2 = tmp_kappa(closest_low_id + 1)
      diffusion_coefficient(i) = (y2 + (y1 - y2) * (grid_sample(i) - x2) / (x1 - x2))
      
    else if (grid_sample(i) .lt. tmp_grid(1)) then
      H_number_density(i) = tmp_density(1)
      gas_temperature(i) = tmp_gas_temperature(1)
      visual_extinction(i) = tmp_av(1)
      diffusion_coefficient(i) = tmp_kappa(1)
    else if (grid_sample(i) .ge. tmp_grid(nb_values)) then
      H_number_density(i) = tmp_density(nb_values)
      gas_temperature(i) = tmp_gas_temperature(nb_values)
      visual_extinction(i) = tmp_av(nb_values)
      diffusion_coefficient(i) = tmp_kappa(nb_values)
    end if
  end do
  
else
  write (Error_Unit,*) 'Error: The file "',trim(filename),'" does not exist.'
  call exit(25)
end if

return
end subroutine init_1D_evolution

end module structure