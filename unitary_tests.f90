!******************************************************************************
! PROGRAM: unitary_tests
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Program that test different functions implemented in nautilus. 
!! Especially routines that were added afterwards, and are not  present in 
!! the original version, thus can't be compared between the two versions.
!!\n\n
!! This program need a test simulation in the sub-folder 'tests'. Output results
!! will be stored in the same folder.
!!\n\n
!! This binary is conceived so that it is run by the python script unitary_tests.py
!
!******************************************************************************
! 

program unitary_tests

  use numerical_types
  use iso_fortran_env, only : error_unit
  use global_variables
  use nautilus_main

    
  implicit none
  
  call initialisation()
  
  call test_read_time_evolution()
  call test_time_interpolation()
  call test_diffusion()
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 7 may 2014
!
! DESCRIPTION: 
!> @brief Test the interpolation of time evolution properties of the structure
!! then compare them to the raw profile, just to be sure the interpolation is
!! correct.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine test_time_interpolation()
  
    implicit none
    
    integer, parameter :: nb_sample = 100
    real(double_precision), parameter :: t_min = 0.1d0 ! in AU
    real(double_precision), parameter :: t_max = 100.d0! in AU
    real(double_precision), parameter :: step = (t_max / t_min)**(1.d0 / (dfloat(nb_sample)))
    
    real(double_precision) :: time ! Time in Myr
    real(double_precision) :: av, density, gas_temperature, grain_temperature

    logical :: isDefined !< true or false, to test if some files exists or not, for instance
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the time evolution interpolation'
    
    inquire(file='structure_evolution.dat', exist=isDefined)
    
    if (.not.isDefined) then
      write(Error_unit, *) "Warning: The file 'structure_evolution.dat' doesn't exist."
      write(Error_unit, *) "A default one is being generated"
      
      open(10, file='structure_evolution.dat')
      write(10,'(a)') '! time    log(Av)    log(n)    log(T)   '
      write(10,'(a)') '! (Myr)   (mag)      (cm-3)    (K)      '
      write(10,'(a)') '0.000e+00 -1.231e+00 1.813e+00 1.698e+00'
      write(10,'(a)') '5.000e+00 -1.221e+00 1.760e+00 1.715e+00'
      write(10,'(a)') '1.000e+01 -1.511e+00 1.327e+00 1.946e+00'
      write(10,'(a)') '1.500e+01 -1.617e+00 1.166e+00 2.063e+00'
      write(10,'(a)') '2.000e+01 -1.906e+00 7.300e-01 2.463e+00'
      write(10,'(a)') '2.500e+01 -1.781e+00 8.800e-01 2.340e+00'
      write(10,'(a)') '3.000e+01 -1.448e+00 1.367e+00 1.930e+00'
      write(10,'(a)') '3.500e+01 -1.091e+00 1.909e+00 1.656e+00'
      write(10,'(a)') '4.000e+01 -1.202e+00 1.748e+00 1.719e+00'
      write(10,'(a)') '4.500e+01 -1.307e+00 1.583e+00 1.805e+00'
      write(10,'(a)') '5.000e+01 -7.650e-01 2.432e+00 1.480e+00'
      write(10,'(a)') '5.300e+01 -1.024e+00 2.058e+00 1.594e+00'
      close(10)
      
    end if
    
    call init_structure_evolution()
    
    open(10, file='test_time_interpolation.dat')
    write(10,*) '# time, av, density, gas_temperature, grain_temperature'
    do j=1,nb_sample
      time = t_min * step**(j - 1.d0)
      ! We generate cartesian coordinate for the given Semi-major axis
      
      call get_structure_properties_table(time=time*(1e6*YEAR),&!Inputs
        Av=av, density=density, gas_temperature=gas_temperature, grain_temperature=grain_temperature)
      
      
      write(10,*) time, av, density, gas_temperature, grain_temperature
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="av_interpolation.gnuplot")
    open(11, file="density_interpolation.gnuplot")
    open(12, file="gas_temperature_interpolation.gnuplot")
    open(13, file="grain_temperature_interpolation.gnuplot")
    
    do j=10,13
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) '!rm "av_interpolation.pdf"'
    write(10,*) "set output 'av_interpolation.pdf'"
    write(10,*) 'set ylabel "Visual extinction [mag]"'
    
    write(11,*) '!rm "density_interpolation.pdf"'
    write(11,*) "set output 'density_interpolation.pdf'"
    write(11,*) 'set ylabel "Density [part/cm^3]"'
    write(11,*) 'set logscale y'
    
    write(12,*) '!rm "gas_temperature_interpolation.pdf"'
    write(12,*) "set output 'gas_temperature_interpolation.pdf'"
    write(12,*) 'set ylabel "Temperature [K]"'
    write(12,*) 'set logscale y'
    
    write(13,*) '!rm "grain_temperature_interpolation.pdf"'
    write(13,*) "set output 'grain_temperature_interpolation.pdf'"
    write(13,*) 'set ylabel "Temperature [K]"'
    write(13,*) 'set logscale y'
    
    do j=10,13
      write(j,*) 'set xlabel "Time [Myr]"'
      write(j,*) 'set grid'
      write(j,*) 'set logscale x'
    end do

    write(10,*) 'plot "test_time_interpolation.dat" using 1:2 with points linetype 1 pointtype 2 title "Interpolation",\'
    write(10,*) '     "structure_evolution.dat" using 1:(10**$2) with lines title "Profile"'

    write(11,*) 'plot "test_time_interpolation.dat" using 1:3 with points linetype 1 pointtype 2 title "Interpolation",\'
    write(11,*) '     "structure_evolution.dat" using 1:(10**$3) with lines title "Profile"'

    write(12,*) 'plot "test_time_interpolation.dat" using 1:4 with points linetype 1 pointtype 2 title "Interpolation",\'
    write(12,*) '     "structure_evolution.dat" using 1:(10**$4) with lines title "Profile"'

    write(13,*) 'plot "test_time_interpolation.dat" using 1:5 with points linetype 1 pointtype 2 title "Interpolation",\'
    write(13,*) '     "structure_evolution.dat" using 1:(10**$5) with lines title "Profile"'
    
    close(10)
    close(11)
    close(12)
    close(13)
  
  end subroutine test_time_interpolation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 13 may 2014
!
! DESCRIPTION: 
!> @brief Test the reading of time evolution, compare the actual profile 
!! with the read one (in global variable)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine test_read_time_evolution()
  
    implicit none
    
    logical :: isDefined !< true or false, to test if some files exists or not, for instance
    
    integer :: j ! for loops
    
    write(*,*) 'Test of read profile from time evolution data file'
    
    inquire(file='structure_evolution.dat', exist=isDefined)
    
    if (.not.isDefined) then
      write(Error_unit, *) "Warning: The file 'structure_evolution.dat' doesn't exist."
      write(Error_unit, *) "A default one is being generated"
      
      open(10, file='structure_evolution.dat')
      write(10,'(a)') '! time    log(Av)    log(n)    log(T)   '
      write(10,'(a)') '! (Myr)   (mag)      (cm-3)    (K)      '
      write(10,'(a)') '0.000e+00 -1.231e+00 1.813e+00 1.698e+00'
      write(10,'(a)') '5.000e+00 -1.221e+00 1.760e+00 1.715e+00'
      write(10,'(a)') '1.000e+01 -1.511e+00 1.327e+00 1.946e+00'
      write(10,'(a)') '1.500e+01 -1.617e+00 1.166e+00 2.063e+00'
      write(10,'(a)') '2.000e+01 -1.906e+00 7.300e-01 2.463e+00'
      write(10,'(a)') '2.500e+01 -1.781e+00 8.800e-01 2.340e+00'
      write(10,'(a)') '3.000e+01 -1.448e+00 1.367e+00 1.930e+00'
      write(10,'(a)') '3.500e+01 -1.091e+00 1.909e+00 1.656e+00'
      write(10,'(a)') '4.000e+01 -1.202e+00 1.748e+00 1.719e+00'
      write(10,'(a)') '4.500e+01 -1.307e+00 1.583e+00 1.805e+00'
      write(10,'(a)') '5.000e+01 -7.650e-01 2.432e+00 1.480e+00'
      write(10,'(a)') '5.500e+01 -1.024e+00 2.058e+00 1.594e+00'
      close(10)
      
    end if
    
    call init_structure_evolution()
    
    open(10, file='test_time_read.dat')
    write(10,*) '# time (Myr), av [mag], density [part/cm^3], gas_temperature (K), grain_temperature (K)'
    do j=1,structure_sample
      
      if (read_dust) then
        write(10,*) structure_time(j)/(1e6*YEAR), 10.d0**(structure_log_Av(j)), 10.d0**(structure_log_density(j)), &
                    10.d0**(structure_log_gas_temperature(j)), 10.d0**(structure_log_dust_temperature(j))
      else
        write(10,*) structure_time(j)/(1e6*YEAR), 10.d0**(structure_log_Av(j)), 10.d0**(structure_log_density(j)), &
                    10.d0**(structure_log_gas_temperature(j))
      endif
      ! if the grain temperature is not defined in the structure_evolution.dat file, grain temperature will appear to be 1K 
      !! since the log will be set to 0
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="test_av_read.gnuplot")
    open(11, file="test_density_read.gnuplot")
    open(12, file="test_gas_temperature_read.gnuplot")
    
    do j=10,12
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) '!rm "test_av_read.pdf"'
    write(10,*) "set output 'test_av_read.pdf'"
    write(10,*) 'set ylabel "Visual extinction [mag]"'
    
    write(11,*) '!rm "test_density_read.pdf"'
    write(11,*) "set output 'test_density_read.pdf'"
    write(11,*) 'set ylabel "Density [part/cm^3]"'
    write(11,*) 'set logscale y'
    
    write(12,*) '!rm "test_gas_temperature_read.pdf"'
    write(12,*) "set output 'test_gas_temperature_read.pdf'"
    write(12,*) 'set ylabel "Temperature [K]"'
    write(12,*) 'set logscale y'
    
    do j=10,12
      write(j,*) 'set xlabel "Time [Myr]"'
      write(j,*) 'set grid'
      write(j,*) 'set logscale x'
    end do

    write(10,*) 'plot "test_time_read.dat" using 1:2 with points linetype 1 pointtype 2 title "Read profile",\'
    write(10,*) '     "structure_evolution.dat" using 1:(10**$2) with lines title "Profile"'

    write(11,*) 'plot "test_time_read.dat" using 1:3 with points linetype 1 pointtype 2 title "Read profile",\'
    write(11,*) '     "structure_evolution.dat" using 1:(10**$3) with lines title "Profile"'

    write(12,*) 'plot "test_time_read.dat" using 1:4 with points linetype 1 pointtype 2 title "Read profile",\'
    write(12,*) '     "structure_evolution.dat" using 1:(10**$4) with lines title "Profile"'
    
    close(10)
    close(11)
    close(12)
    
    if (read_dust) then
      open(13, file="test_grain_temperature_read.gnuplot")
      
      write(13,*) "set terminal pdfcairo enhanced"
      write(13,*) '!rm "test_grain_temperature_read.pdf"'
      write(13,*) "set output 'test_grain_temperature_read.pdf'"
      write(13,*) 'set ylabel "Temperature [K]"'
      write(13,*) 'set logscale y'
      write(13,*) 'set xlabel "Time [Myr]"'
      write(13,*) 'set grid'
      write(13,*) 'set logscale x'
      
      select case(GRAIN_TEMPERATURE_TYPE)
      case('gas')
        write(13,*) 'plot "test_time_read.dat" using 1:5 with points linetype 1 pointtype 2 title "Read profile",\'
        write(13,*) '     "structure_evolution.dat" using 1:(10**$4) with lines title "Profile"'
      
      case('table')
        write(13,*) 'plot "test_time_read.dat" using 1:5 with points linetype 1 pointtype 2 title "Read profile",\'
        write(13,*) '     "structure_evolution.dat" using 1:(10**$5) with lines title "Profile"'

      case default
        write(13,*) 'plot "test_time_read.dat" using 1:5 with points linetype 1 pointtype 2 title "Read profile"'
      end select
      
      close(13)
    endif
    
  
  end subroutine test_read_time_evolution

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 24 july 2014
!
! DESCRIPTION: 
!> @brief Uses the diffusion type defined in the code (constrained in
!! parameters.in).This routine must be run at the end of your tests, because you override 
!! abundances\n
!! Files are generated in a sub folder "dissipation" of the test folder. 
!! Gnuplot file is also here.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine test_diffusion()
  
  implicit none
    
    
    ! time sample
    real(double_precision), parameter :: t_max = 1.d6 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    integer :: time_size ! the size of the array 'time'. 
    
    character(len=80) :: filename_density
    character(len=80) :: output_density, output_time, time_format, purcent_format
    integer :: time_length ! the length of the displayed time, usefull for a nice display
    real(double_precision) :: timestep
    
    integer :: i,k ! for loops
    integer :: nb_time ! The total number of 't' values. 
    integer :: error ! to retrieve error, especially during allocations
    logical :: isDefined

    real(double_precision), dimension(nb_species, nb_sample_1D) :: fake_abundances ! Fill an array where the first species is 
    !! the one we want to test. All the others are put to equal values. 
    !------------------------------------------------------------------------------
    write(*,*) 'Test of the diffusion processes'
    
    fake_abundances(1:nb_species, 1:nb_sample_1D) = 0.d0
    fake_abundances(1, nb_sample_1D/2) = 1.d0 ! Dirac function for the first species    
    
    if (structure_type.eq.'0D') then
      write (error_unit,*) 'Error: There is currently no dissipation (structure_type=0) which is a problem to test it.'
      write (error_unit,*) 'Please set a dissipation_type in "disk.in" before testing it'
      call exit(3)
    end if
    
    inquire(file='dissipation', exist=isDefined)
    
    ! We create the folder 'dissipation' if he doesn't exists.
    if (.not.isDefined) then
      call system("mkdir dissipation")
    end if
    
    call system("rm dissipation/*")
    
    ! We want to know the max size of the time display in order to have a nice display, with filled spaces in the final plots
    write(output_time, '(i0)') int(t_max / 1e6)
    time_length = len(trim(output_time))
    write(time_format, '(a,i1,a,i1,a)') '(f',time_length+2,'.',1,')'
    write(purcent_format, *) '(', trim(time_format), ',"/",',trim(time_format),'," Myr ; k = ",i5)'
    
    !------------------------------------------------------------------------------
    k = 1
    time_size = 512 ! the size of the array. 
    allocate(time(time_size), stat=error)
    
    current_time = 0.d0 ! s
    
    do while (current_time.lt.(t_max*YEAR))
      time(k) = current_time
    
      ! We open the file where we want to write the outputs
      write(filename_density, '(a,i0.5,a)') 'dissipation/abundance.',k,'.dat'    
      
      ! We take on purpose an extra large final time to ensure we get the pure diffusive timestep, without any problems with output times.
      call get_timestep(current_time=current_time, final_time=t_max*YEAR, next_timestep=timestep)
      
      current_time = time(k) + timestep
      ! We generate cartesian coordinate for the given Semi-major axis
      
      call structure_diffusion(timestep=timestep, temp_abundances=fake_abundances(1:nb_species, 1:nb_sample_1D))
      
      open(10, file=filename_density)
      do i=1, nb_sample_1D
        write(10,*) grid_sample(i)/AU, fake_abundances(1, i)
      enddo
      close(10)
      
      ! we expand the 'time' array if the limit is reached
      if (k.eq.time_size) then
        ! If the limit of the 'time' array is reach, we copy the values in a temporary array, allocate with a double size, et paste the 
        ! old values in the new bigger array
        allocate(time_temp(time_size), stat=error)
        time_temp(1:time_size) = time(1:time_size)
        deallocate(time, stat=error)
        time_size = time_size * 2
        allocate(time(time_size), stat=error)
        time(1:time_size/2) = time_temp(1:time_size/2)
        deallocate(time_temp, stat=error)
      end if
      
      write(*,purcent_format) time(k)/(YEAR * 1e6), t_max / 1e6, k
      
      
      k = k + 1 ! We increment the integer that point the time in the array (since it's a 'while' and not a 'do' loop)
    end do

    nb_time = k - 1 ! since for the last step with incremented 'k' by one step that is beyond the limit.
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the density
    open(13, file="dissipation/dissipation.gnuplot")
    write(13,*) "set terminal pngcairo crop enhanced size 1200, 1000 font ',20'"
    write(13,*) 'set xlabel "Distance (AU)"'
    write(13,*) 'set ylabel "Quantity"'
    write(13,*) 'set grid'
    write(13,*) 'set xrange [', grid_sample(1)/AU, ':', grid_sample(nb_sample_1D)/AU, ']'
    write(13,*) 'set yrange [', 0, ':', 1, ']'
    
    do k=1, nb_time
      write(filename_density, '(a,i0.5,a)') 'abundance.',k,'.dat'
      write(output_density, '(a,i0.5,a)') 'abundance.',k,'.png'
      write(output_time, time_format) (time(k)/(YEAR * 1e6))
      
      write(13,*) "set output '",trim(output_density),"'"
      write(13,*) 'set title "T=', trim(output_time),' Myr"'
      write(13,*) "plot '",trim(filename_density),"' using 1:2 with lines linetype -1 notitle"
      write(13,*) ""
    end do
    close(13)
  
  end subroutine test_diffusion

end program unitary_tests
