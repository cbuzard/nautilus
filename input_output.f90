!******************************************************************************
! MODULE: input_output
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines linked to reading input files
!! or writing output files. \n\n
!!
!! Input files tends to be named *.in
!! Output files tends to be named *.out
!! Temporary files that are overwritten at each timestep are named *.tmp
!
!******************************************************************************

module input_output

use iso_fortran_env


implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Read all the simulation parameters in a file named 
!! 'parameters.in'
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_input_files()

use global_variables

implicit none

! Locals
integer :: i,k,j, jk

character(len=80) :: filename !< name of the file to be read
character(len=200) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isDefined

! Variables for the unordered reaction file
character(len=11), dimension(7,nb_gas_phase_reactions) :: SYMBOLUO1 
real(double_precision), dimension(nb_gas_phase_reactions) :: AUO1,BUO1,CUO1 
integer, dimension(nb_gas_phase_reactions) :: itypeUO1,Tmin1,Tmax1,FORMULA1,NUM1 

character (len=11), dimension(7,nb_surface_reactions) :: SYMBOLUO2 
real(double_precision), dimension(nb_surface_reactions) :: AUO2,BUO2,CUO2 
integer, dimension(nb_surface_reactions) :: itypeUO2,Tmin2,Tmax2,FORMULA2,NUM2 

character (len=11), dimension(7,nb_reactions) :: SYMBOLUO
real(double_precision), dimension(nb_reactions) :: AUO,BUO,CUO
integer, dimension(nb_reactions) :: itypeUO,TminUO,TmaxUO,FORMULAUO,NUMUO

call read_parameters_in()

call read_species()

! Reading list of reaction for gas phase
filename = 'gas_reactions.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  j = 0
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
      j = j + 1
      read(line, '(3A11,1x,4A11,11x,3D11.3,23x,I3,2i7,i3,i6)')  (SYMBOLUO1(I,J),I=1,7),AUO1(J),BUO1(J),CUO1(J), &
ITYPEUO1(J),Tmin1(j),Tmax1(j),FORMULA1(J),NUM1(J)
    
    end if
  end do
  close(10)
  
else
  write (Error_unit,*) 'Error: The file ', filename,' does not exist.'
  call exit(1)
end if

! Reading list of reaction for grains
filename = 'grain_reactions.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  j = 0
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
      j = j + 1
      read(line, '(3A11,1x,4A11,11x,3D11.3,23x,I3,2i7,i3,i6)')  (SYMBOLUO2(I,J),I=1,7),AUO2(J),BUO2(J),CUO2(J), &
ITYPEUO2(J),Tmin2(j),Tmax2(j),FORMULA2(J),NUM2(J)
    
    end if
  end do
  close(10)
  
else
  write (Error_unit,*) 'Error: The file ', filename,' does not exist.'
  call exit(1)
end if

! putting everything back into the big tables
do I=1,nb_gas_phase_reactions 
  do k=1,7
    SYMBOLUO(k,I)=SYMBOLUO1(k,I)
  enddo
  AUO(I)=AUO1(I)
  BUO(I)=BUO1(I)
  CUO(I)=CUO1(I)
  ITYPEUO(I)=ITYPEUO1(I)
  TminUO(I) = Tmin1(I)
  TmaxUO(I) = Tmax1(I)
  FORMULAUO(I) = FORMULA1(I)
  NUMUO(I) = NUM1(I)
enddo

do I=1,nb_surface_reactions 
  do  k=1,7
    SYMBOLUO(k,nb_gas_phase_reactions+I)=SYMBOLUO2(k,I)
  enddo
  AUO(nb_gas_phase_reactions+I)=AUO2(I)
  BUO(nb_gas_phase_reactions+I)=BUO2(I)
  CUO(nb_gas_phase_reactions+I)=CUO2(I)
  ITYPEUO(nb_gas_phase_reactions+I)=ITYPEUO2(I)
  TminUO(nb_gas_phase_reactions+I) = Tmin2(I)
  TmaxUO(nb_gas_phase_reactions+I) = Tmax2(I)
  FORMULAUO(nb_gas_phase_reactions+I) = FORMULA2(I)
  NUMUO(nb_gas_phase_reactions+I) = NUM2(I)
enddo

! Reorder reaction file entries with ITYPE
jk=1
do i=0,nitype
  do j=1,nb_reactions
    if (itypeuo(j).eq.i) then
      SYMBOL(:,jk)=SYMBOLUO(:,j)     
      A(jk)=AUO(j)
      B(jk)=BUO(j)
      C(jk)=CUO(j)
      ITYPE(jk)=itypeuo(j)
      Tmin(jk) = dble(TminUO(j))
      Tmax(jk) = dble(TmaxUO(j))
      FORMULA(jk) = FORMULAUO(J)
      NUM(jk) = NUMUO(j)
      jk=jk+1
    endif
  enddo
enddo


if (jk.ne.nb_reactions+1) then
  write(*,*) 'Some reaction was not found by the reorder process'
  write(*,*) jk,'=/',nb_reactions+1 
  stop
endif

!       replace the species names by blanks for non chemical species                                                                        
do j=1,nb_reactions-1
  do i=1,7
    select case(symbol(i,j))
      case ('CR', 'CRP', 'Photon')
        symbol(i,j) = '           '
    end select
  enddo

enddo 


return
end subroutine read_input_files

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read all species information from gas_species.in and grain_species.in
!! files.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_species()

use global_variables
use iso_fortran_env

implicit none

! Locals
integer :: i,k

character(len=80) :: filename !< name of the file to be read
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isDefined

! Variables for the unordered reaction file
character(len=11), dimension(nb_species_for_gas) :: gas_species_label 
integer, dimension(nb_species_for_gas) :: ICG1 
integer, dimension(nemax, nb_species_for_gas) :: IELM1 

character(len=11), dimension(nb_species_for_grain) :: surface_species_label 
integer, dimension(nb_species_for_grain) :: ICG2 
integer, dimension(nemax, nb_species_for_grain) :: IELM2 


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
      read(line, '(A11,i3,13(I3))')  gas_species_label(I),ICG1(I),(IELM1(K,I),K=1,NEMAX) 
    
    end if
  end do
  close(10)
  
else
  write (Error_unit,*) 'Error: The file ', filename,' does not exist.'
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
      read(line, '(A11,i3,13(I3))')  surface_species_label(I),ICG2(I),(IELM2(K,I),K=1,NEMAX) 
    
    end if
  end do
  close(10)
  
else
  write (Error_unit,*) 'Error: The file ', filename,' does not exist.'
  call exit(1)
end if


! putting everything back into the big tables
do I=1,nb_species_for_gas 
  species_name(I)=gas_species_label(I)
  ICG(I)=ICG1(I)
  do k=1,NEMAX
    IELM(K,I)=IELM1(K,I)
  enddo
enddo
do I=1,nb_species_for_grain 
  species_name(nb_species_for_gas+I)=surface_species_label(I)
  ICG(nb_species_for_gas+I)=ICG2(I)
  do k=1,NEMAX
    IELM(K,nb_species_for_gas+I)=IELM2(K,I)
  enddo
enddo

return
end subroutine read_species

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read simulation parameters from the file parameters.in
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_parameters_in()
use global_variables

implicit none

character(len=80) :: filename = 'parameters.in' !< name of the file in which parameters are stored
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isParameter, isDefined
character(len=80) :: identificator, value
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    call get_parameter_value(line, isParameter, identificator, value)
      
    if (isParameter) then
      select case(identificator)
      ! Solver
      case('relative_tolerance')
        read(value, '(e12.6)') RELATIVE_TOLERANCE
      
      !Switches
      case('is_grain_reactions')
        read(value, '(i2)') IS_GRAIN_REACTIONS
      
      case('is_absorption')
        read(value, '(i2)') IS_ABSORPTION
      
      case('grain_tunneling_diffusion')
        read(value, '(i2)') GRAIN_TUNNELING_DIFFUSION
      
      case('IMODH')
        read(value, '(i2)') IMODH
      
      case('conservation_type')
        read(value, '(i2)') CONSERVATION_TYPE
      
      ! Gas phase
      case('initial_gas_density')
        read(value, '(e12.6)') initial_gas_density
      
      case('initial_gas_temperature')
        read(value, '(e12.6)') initial_gas_temperature
      
      case('initial_visual_extinction')
        read(value, '(e12.6)') INITIAL_VISUAL_EXTINCTION
      
      case('cr_ionisation_rate')
        read(value, '(e12.6)') CR_IONISATION_RATE
      
      case('x_ionisation_rate')
        read(value, '(e12.6)') X_IONISATION_RATE
      
      case('uv_flux')
        read(value, '(e12.6)') UV_FLUX
      
      ! Grain
      case('initial_dust_temperature')
        read(value, '(e12.6)') initial_dust_temperature
      
      case('initial_dtg_mass_ratio')
        read(value, '(e12.6)') initial_dtg_mass_ratio
      
      case('sticking_coeff_neutral')
        read(value, '(e12.6)') sticking_coeff_neutral
      
      case('sticking_coeff_positive')
        read(value, '(e12.6)') sticking_coeff_positive
      
      case('sticking_coeff_negative')
        read(value, '(e12.6)') sticking_coeff_negative
      
      case('grain_density')
        read(value, '(e12.6)') GRAIN_DENSITY
      
      case('grain_radius')
        read(value, '(e12.6)') GRAIN_RADIUS
        
      case('site_spacing')
        read(value, '(e12.6)') SITE_SPACING
      
      case('site_density')
        read(value, '(e12.6)') SITE_DENSITY
      
      case('EBFAC')
        read(value, '(e12.6)') EBFAC
      
      case('ACT')
        read(value, '(e12.6)') ACT
      
      case('TSMAX')
        read(value, '(e12.6)') TSMAX
      
      case('CRT')
        read(value, '(e12.6)') CRT
      
      case('CRFE')
        read(value, '(e12.6)') CRFE
      
      case('LAYERS')
        read(value, '(e12.6)') LAYERS
      
      case('ARRK')
        read(value, '(e12.6)') ARRK
      
      ! Outputs
      case('output_per_decade')
        read(value, '(i4)') OUTPUT_PER_DECADE
      
      case('start_time')
        read(value, '(e12.6)') START_TIME
      
      case('stop_time')
        read(value, '(e12.6)') STOP_TIME
      
      case('WSTEP')
        read(value, '(i5)') WSTEP
      
      case('WSTEPR')
        read(value, '(i5)') WSTEPR
      
      case('IRATEOUT')
        read(value, '(i5)') IRATEOUT
      
      ! Initial abundances
      case('minimum_initial_abundance')
        read(value, '(e12.6)') MINIMUM_INITIAL_ABUNDANCE
      
      ! 1D Parameters
      case('is_diffusivity')
        read(value, '(i5)') IS_DIFFUSIVITY
      
      case('turbulent_diffusivity')
        read(value, '(e12.6)') TURBULENT_DIFFUSIVITY
      
      case('box_size')
        read(value, '(e12.6)') BOX_SIZE
        
      case('central_mass')
        read(value, '(e12.6)') CENTRAL_MASS
      
      case('radial_distance')
        read(value, '(e12.6)') RADIAL_DISTANCE
      
      case('TESTJAC')
        read(value, '(i5)') TESTJAC
        
      case('NJAC')
        read(value, '(i5)') NJAC
         
      case default
        write(*,*) 'Warning: An unknown parameter has been found'
        write(*,*) "identificator='", trim(identificator), "' ; value(s)='", trim(value),"'"
      end select
    end if
  end do
  close(10)
  
else
  write (*,*) 'Warning: The file "parameters.in" does not exist. Default values have been used'
end if

if ((IS_DIFFUSIVITY.lt.0).or.(IS_DIFFUSIVITY.gt.2)) then
  write(*,*) 'This value for IS_DIFFUSIVITY is not implemented: ',IS_DIFFUSIVITY
  stop
endif

START_TIME = START_TIME * TYEAR
STOP_TIME = STOP_TIME * TYEAR

return
end subroutine read_parameters_in

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine that write the simulation parameters into the file 'parameters.out'
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_parameters()
use global_variables

  implicit none
  
  character(len=80) :: filename = 'parameters.in'
  
  open(10, file=filename)
  write(10,'(a)') "!# ------------------------------------------------"
  write(10,'(a)') "!# Parameter file for various properties of the disk."
  write(10,'(a)') "!# ------------------------------------------------"
  write(10,'(a)') "!# blanck line or with spaces will be skipped."
  write(10,'(a)') "!# In fact, the only lines that matter are non commented lines with a"
  write(10,'(a)') "!# '=' character to distinguish the identificator and the value(s)"
  write(10,'(a)') "!# (each value must be separated with at least one space."
  write(10,'(a)') "!# Line must not be longer than 80 character, but comments can be far"
  write(10,'(a)') "!# bigger than that, even on line with a parameter to read."
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*     Solver Parameters     *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2, a)') 'relative_tolerance = ',RELATIVE_TOLERANCE, ' ! Relative tolerance'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*          Switches         *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i2,a)') 'is_grain_reactions = ', IS_GRAIN_REACTIONS, ' ! Accretion, grain surface reactions'
  write(10,'(a,i2,a)') 'is_absorption = ', IS_ABSORPTION, ' ! H2 AND CO SELF-SHIELDING'
  write(10,'(a,i2,a)') 'grain_tunneling_diffusion = ', GRAIN_TUNNELING_DIFFUSION, &
  ' ! 0=thermal; For H,H2: 1=QM1; 2=QM2; 3=choose fastest'
  write(10,'(a,i2,a)') 'IMODH = ', IMODH, ' ! 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only'
  write(10,'(a,i2,a)') 'conservation_type = ', CONSERVATION_TYPE, ' ! 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*    Gas phase parameters   *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'initial_gas_density = ', initial_gas_density, ' ! initial gas density (part/cm-3)'
  write(10,'(a,es10.3e2,a)') 'initial_gas_temperature = ', initial_gas_temperature, ' ! initial gas temperature (K)'
  write(10,'(a,es10.3e2,a)') 'initial_visual_extinction = ', INITIAL_VISUAL_EXTINCTION, ' ! initial visual extinction'
  write(10,'(a,es10.3e2,a)') 'cr_ionisation_rate = ', CR_IONISATION_RATE, ' ! cosmic ray ionisation rate (1.3e-17 standard value)'
  write(10,'(a,es10.3e2,a)') 'x_ionisation_rate = ', X_IONISATION_RATE, ' ! Ionisation rate due to X-rays (s-1)'
  write(10,'(a,es10.3e2,a)') 'uv_flux = ', UV_FLUX, ' ! Scale factor for the UV flux, in unit of the reference flux (1.=nominal)'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*      Grain parameters     *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'initial_dust_temperature = ', initial_dust_temperature, ' ! initial dust temp'
  write(10,'(a,es10.3e2,a)') 'initial_dtg_mass_ratio = ', initial_dtg_mass_ratio, ' ! dust-to-gas ratio by mass'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_neutral = ', sticking_coeff_neutral, ' ! sticking coeff for neutral species'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_positive = ', sticking_coeff_positive, ' ! sticking coeff for positive species'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_negative = ', sticking_coeff_negative, ' ! sticking coeff for negative species'
  write(10,'(a,es10.3e2,a)') 'grain_density = ', GRAIN_DENSITY, ' ! mass density of grain material'
  write(10,'(a,es10.3e2,a)') 'grain_radius = ', grain_radius, ' ! grain radius (cm)'
  write(10,'(a,es10.3e2,a)') 'site_spacing = ', SITE_SPACING, ' ! site spacing (cm)'
  write(10,'(a,es10.3e2,a)') 'site_density = ', SITE_DENSITY, ' ! site density (cm-2)'
  write(10,'(a,es10.3e2,a)') 'EBFAC = ', EBFAC, ' ! ratio Eb(I):Ed(I) (excludes H,H2); -ve means use given values'
  write(10,'(a,es10.3e2,a)') 'ACT = ', ACT, ' ! grain rxn activation energy constant'
  write(10,'(a,es10.3e2,a)') 'TSMAX = ', TSMAX, ' ! peak grain temp (CR heating)'
  write(10,'(a,es10.3e2,a)') 'CRT = ', CRT, ' ! duration (s) of peak grain T'
  write(10,'(a,es10.3e2,a)') 'CRFE = ', CRFE, ' ! Fe-ion--grain encounter s-1 grain-1 (for 0.1 micron grain)'
  write(10,'(a,es10.3e2,a)') 'LAYERS = ', LAYERS, ' ! number of monolayers for ITYPE 17-20 currently not used in the code'
  write(10,'(a,es10.3e2,a)') 'ARRK = ', ARRK, ' ! a-coefficient for RRK-style formation-desorption'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*        Output times       *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i4,a)') 'output_per_decade = ', OUTPUT_PER_DECADE, ' ! outputs per decade (Only without diffusion)'
  write(10,'(a,es10.3e2,a)') 'start_time = ', START_TIME/TYEAR, ' ! first output time, after zero (yrs)'
  write(10,'(a,es10.3e2,a)') 'stop_time = ', STOP_TIME/TYEAR, ' ! last output time (yrs)'
  write(10,'(a,i5,a)') 'WSTEP = ', WSTEP, ' ! Outputs every WSTEP timesteps (/=1 only for 1D outputs)'
  write(10,'(a,i5,a)') 'WSTEPR = ', WSTEPR, ' ! Outputs every WSTEPR timesteps for the rate coefficients'
  write(10,'(a,i5,a)') 'IRATEOUT = ', IRATEOUT, ' ! Spatial point for the rate output'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*     Initial Abundances    *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'minimum_initial_abundance = ', MINIMUM_INITIAL_ABUNDANCE, ' ! default minimum initial frac abun'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*       1D parameters       *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i5,a)') 'is_diffusivity = ', IS_DIFFUSIVITY, ' ! Diffusivity flag (1. constant 2. alpha)'
  write(10,'(a,es10.3e2,a)') 'turbulent_diffusivity = ', TURBULENT_DIFFUSIVITY, &
  ' ! Turbulent diffusivity (cgs if 1) (value or parameter)'
  write(10,'(a,es10.3e2,a)') 'box_size = ', BOX_SIZE, ' ! Computing box"s size (cm)'
  write(10,'(a,es10.3e2,a)') 'central_mass = ', CENTRAL_MASS, ' ! Central mass (g)'
  write(10,'(a,es10.3e2,a)') 'radial_distance = ', RADIAL_DISTANCE, ' ! Radial distance (cm) for the disk model'
  write(10,'(a,i5,a)') 'TESTJAC = ', TESTJAC, ' ! Testing the number of non-zero jacobian elements ?'
  write(10,'(a,i5,a)') 'NJAC = ', NJAC, ' ! Number of non-zero jacobian elements (per line, result of TESTJAC=1)'
  write(10,*) ''
  close(10)
  
end subroutine write_parameters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write the list of species and the corresponding index in an
!! output file 'species.out'.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_species()
use global_variables

implicit none

! Locals
integer :: i

open(4, file='species.out')
! Write 'ggo_spec.d': 5 columns of numbered species=====================
write(4,'(5(I4,")",1X,A11,1X))') (I,species_name(I),I=1,nb_species)
close(4)

return
end subroutine write_species

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine that try to split the line in two part, given a 
!! separator value (set in parameter of the subroutine)
!
!> @warning The first character of the parameter value MUST NOT be a string. All spaces will be truncated
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_parameter_value(line, isParameter, id, value)

  implicit none
  
  ! Input
  character(len=80), intent(in) :: line !< [in] Input line in which we want to retrieve a parameter and its value
  
  ! Output
  logical, intent(out) :: isParameter !< [out] a boolean to say whether or not there is a parameter on this line. 
!! i.e if there is an occurence of the separator in the input line
  character(len=80), intent(out) :: id !< [out] the name of the parameter
  character(len=80), intent(out) :: value !< [out] a string that contains the value(s) associated with the parameter name. 
!!         Note that a special attention is given to the fact that the first character of 'value' must NOT be a 'space'
  
  ! Local
  character(len=1), parameter :: SEP = '=' ! the separator of a parameter line
  character(len=1) :: first_character
  integer :: id_first_char
  integer :: sep_position ! an integer to get the position of the separator

  !------------------------------------------------------------------------------

  sep_position = index(line, SEP)
  
  if (sep_position.ne.0) then
    isParameter = .true.
    id = line(1:sep_position-1)
    
    id_first_char = sep_position +1
    first_character = line(id_first_char:id_first_char)
    do while (first_character.eq.' ')
      id_first_char = id_first_char +1
      first_character = line(id_first_char:id_first_char)
    end do
    value = line(id_first_char:)
  else
    isParameter = .false.
  end if

end subroutine get_parameter_value

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read abundances from abundances.in file. All abundances not defined
!! here will have the default value MINIMUM_INITIAL_ABUNDANCE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_abundances()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80), parameter :: filename='abundances.in'
integer :: i, j

real(double_precision), allocatable, dimension(:) :: temp_abundances
character(len=11), allocatable, dimension(:) :: temp_names
integer :: nb_lines

character(len=80) :: line
character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
integer :: comment_position ! the index of the comment character on the line. if zero, there is none on the current string
integer :: error ! to store the state of a read instruction

logical :: isParameter, isDefined
character(len=80) :: identificator, value

call get_linenumber(filename, nb_lines)

allocate(temp_abundances(nb_lines))
allocate(temp_names(nb_lines))

!~ open(5, file=filename)
!~ read(5,15) (temp_names(I),temp_abundances(I),I=1,nb_lines)
!~ close(5)
!~ 15 format(A11,3X,E12.6)


  !------------------------------------------------------------------------------
  
inquire(file=filename, exist=isDefined)
if (isDefined) then

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
    
    call get_parameter_value(line, isParameter, identificator, value)
      
    if (isParameter) then
      read(value, '(e12.6)') temp_abundances(I)
      read(identificator, *) temp_names(I)
      i = i + 1
    end if
  enddo
  close(10)
endif

! We check if all species in abundances.in exists in the simulation
do j=1,nb_lines
  error = 1
  do i=1,nb_species
    if (temp_names(j).eq.species_name(i)) then
      error = 0 ! The species exist
    endif
  enddo
  
  if (error.eq.1) then
    write(*,*) j
    write(*,*) temp_names
    write(Error_Unit,*) 'Input species "', trim(temp_names(j)), '" in "', trim(filename), '" do not match those in reaction file'
    stop
  endif
enddo

! Set initial abundances================================================
do I=1,nb_species
  abundances(I) = MINIMUM_INITIAL_ABUNDANCE
  do j=1,nb_lines
    if (species_name(I).EQ.temp_names(j)) then
      abundances(I)=temp_abundances(j)
    endif
  enddo
enddo

deallocate(temp_abundances)
deallocate(temp_names)

return
end subroutine read_abundances

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write all abundances for all species in an output file at each
!! output time. The total number of output files related to abundances 
!! will be equal to the number of timestep, not equally spaced in time.\n\n
!! Output filename is of the form : abundances.000001.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_current_output()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'abundances.',timestep,'.out'


open(UNIT=35, file=filename_output, form='unformatted')

write(35) TIME, zspace, species_name
write(35) TEMP1D, DENS1D, TAU1D, X_IONISATION_RATE1D
write(35) ZXN
close(35)

return
end subroutine write_current_output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write rates of all chemical reactions for the current timestep.
!! The total number of files will be equal to the total number of timesteps, the
!! routine being called at the end of each timestep.\n\n
!! Output filename is of the form : rates.000001.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_current_rates()

use global_variables

implicit none

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'rates.',timestep,'.out'

open(45, file=filename_output, form='unformatted')

write(45) species_name
write(45) SYMBOL
write(45) XK
write(45) NUM

close(45)

return 
end subroutine write_current_rates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write the total chemical composition of the actual timestep in 
!! a file whose name is given as an input parameter. This allow to use the
!! same routine to write input, temporary and output files
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_abundances(filename)

use global_variables

implicit none

! Input
character(len=*), intent(in) :: filename !< [in] the name of the output file

! Locals
integer :: i

open(13, file=filename)
write(13,'("!DEPTH POINT=",I2,"/",I2,", TIME =",1PD10.3," s",&
&", XNT=",1PD10.3," cm-3",", TEMP=",1PD10.3," K",&
&", TAU=",0PF8.3,", ZETA=",1PD10.3," s-1")') 00,00,TIME,XNT,TEMP,TAU,CR_IONISATION_RATE

do i=1,nb_species
  write(13,'(a," = ",1PE12.5)') species_name(I),abundances(I)
enddo

write(13,*)
close(13)

return
end subroutine write_abundances




end module input_output