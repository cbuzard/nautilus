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

! Variables for the unordered reaction file

character(len=11), dimension(nb_species_for_gas) :: gas_species_label 
integer, dimension(nb_species_for_gas) :: ICG1 
integer, dimension(nemax, nb_species_for_gas) :: IELM1 

character(len=11), dimension(nb_species_for_grain) :: surface_species_label 
integer, dimension(nb_species_for_grain) :: ICG2 
integer, dimension(nemax, nb_species_for_grain) :: IELM2 

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

! Read CO and H2 shielding factors=====================
open(unit=17, file='gg_H2_Photodiss.d', status='OLD')
do I=1,NL1
  read(17,'(E9.3,3X,E9.3)') N1H2(I),T1H2(I)
enddo
close(17)

open(UNIT=16,FILE='gg_CO_Photodiss.d',STATUS='OLD')
do I=1,NL3
  read(16,'(6(E9.3,3X))') N2CO(I),T2CO(I),N2H2(I),T2H2(I),AV2(I),T2AV(I)              
enddo

do I=NL3+1,NL2
  read(16,'(E9.3,3X,E9.3)') N2CO(I),T2CO(I)
enddo
close(16)

! Read species & reaction info from reactions file======================
! WV fev 2012
! There are now two different files in which the reactions and species are

! Reading the gas phase network
open(unit=9, file='gas_species.in',status='OLD')
do I=1,nb_species_for_gas
  read(9,'(A11,i3,13(I3))') gas_species_label(I),ICG1(I),(IELM1(K,I),K=1,NEMAX) 
enddo
close(9)

open(unit=9, file='gas_reactions.in',status='OLD')
read(9,'(3A11,1x,4A11,11x,3D11.3,23x,I3,2i7,i3,i6)') ((SYMBOLUO1(I,J),I=1,7),AUO1(J),BUO1(J),CUO1(J), &
ITYPEUO1(J),Tmin1(j),Tmax1(j),FORMULA1(J),NUM1(J),J=1,nb_gas_phase_reactions) 
close(9)

! Reading the grain network
open(unit=19, file='grain_species.in', status='OLD')
do I=1,nb_species_for_grain
  read(19,'(A11,i3,13(I3))') surface_species_label(I),ICG2(I),(IELM2(K,I),K=1,NEMAX) 
enddo
close(19)

open(unit=19, file='grain_reactions.in', status='OLD')
read(19,'(3A11,1x,4A11,11x,3D11.3,23x,I3,2i7,i3,i6)') ((SYMBOLUO2(I,J),I=1,7),AUO2(J),BUO2(J),CUO2(J), &
ITYPEUO2(J),Tmin2(j),Tmax2(j),FORMULA2(J),NUM2(J),J=1,nb_surface_reactions) 
close(19)

! putting everything back into the big tables

do I=1,nb_species_for_gas 
  SPEC(I)=gas_species_label(I)
  ICG(I)=ICG1(I)
  do k=1,NEMAX
    IELM(K,I)=IELM1(K,I)
  enddo
enddo
do I=1,nb_species_for_grain 
  SPEC(nb_species_for_gas+I)=surface_species_label(I)
  ICG(nb_species_for_gas+I)=ICG2(I)
  do k=1,NEMAX
    IELM(K,nb_species_for_gas+I)=IELM2(K,I)
  enddo
enddo

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
      Tmin(jk) = real(TminUO(j))
      Tmax(jk) = real(TmaxUO(j))
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
!> @brief Read simulation parameters from the file parameters.in
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_parameters_in()
use global_variables

implicit none

character(len=80) :: filename = 'parameters.in' !< name of the file in which parameters are stored
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position ! the index of the comment character on the line. if zero, there is none on the current string
integer :: error ! to store the state of a read instruction
integer :: boolean ! integer value used to define a logical value (a bit complicated to define directly a boolean)

logical :: isParameter, isDefined
character(len=80) :: identificator, value
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  do
    read(10, '(a80)', iostat=error) line
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
      case('IS_GRAIN_REACTIONS')
        read(value, '(i2)') IS_GRAIN_REACTIONS
      
      case('IS_ABSORPTION')
        read(value, '(i2)') IS_ABSORPTION
      
      case('IGRQM')
        read(value, '(i2)') IGRQM
      
      case('IMODH')
        read(value, '(i2)') IMODH
      
      case('ICONS')
        read(value, '(i2)') ICONS
      
      case('IREAD')
        read(value, '(i2)') IREAD
      
      ! Gas phase
      case('XNT0')
        read(value, '(e12.6)') initial_gas_density
      
      case('TEMP0')
        read(value, '(e12.6)') initial_gas_temperature
      
      case('TAU0')
        read(value, '(e12.6)') TAU0
      
      case('ZETA0')
        read(value, '(e12.6)') ZETA0
      
      case('ZETAX')
        read(value, '(e12.6)') ZETAX
      
      case('UVGAS')
        read(value, '(e12.6)') UVGAS
      
      ! Grain
      case('DTEMP0')
        read(value, '(e12.6)') initial_dust_temperature
      
      case('DTOGM')
        read(value, '(e12.6)') initial_dtg_mass_ratio
      
      case('STICK0')
        read(value, '(e12.6)') sticking_coeff_neutral
      
      case('STICKP')
        read(value, '(e12.6)') sticking_coeff_positive
      
      case('STICKN')
        read(value, '(e12.6)') sticking_coeff_negative
      
      case('RHOD')
        read(value, '(e12.6)') RHOD
      
      case('RD')
        read(value, '(e12.6)') grain_radius
        
      case('ACM')
        read(value, '(e12.6)') ACM
      
      case('SNS')
        read(value, '(e12.6)') SNS
      
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
      case('OTPD')
        read(value, '(i4)') OTPD
      
      case('TSTART')
        read(value, '(e12.6)') TSTART
      
      case('TFINAL')
        read(value, '(e12.6)') TFINAL
      
      case('WSTEP')
        read(value, '(i5)') WSTEP
      
      case('WSTEPR')
        read(value, '(i5)') WSTEPR
      
      case('IRATEOUT')
        read(value, '(i5)') IRATEOUT
      
      ! Initial abundances
      case('XNMIN')
        read(value, '(e12.6)') XNMIN
      
      ! 1D Parameters
      case('IDIFF')
        read(value, '(i5)') IDIFF
      
      case('DIFFTY')
        read(value, '(e12.6)') DIFFTY
      
      case('HSIZE')
        read(value, '(e12.6)') HSIZE
        
      case('MCENTER')
        read(value, '(e12.6)') MCENTER
      
      case('DISTR')
        read(value, '(e12.6)') DISTR
      
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

TSTART = TSTART * TYEAR
TFINAL = TFINAL * TYEAR

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
  
  character(len=80) :: filename = 'parameters.out'
  
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
  write(10,'(a,i2,a)') 'IS_GRAIN_REACTIONS = ', IS_GRAIN_REACTIONS, ' ! Accretion, grain surface reactions'
  write(10,'(a,i2,a)') 'IS_ABSORPTION = ', IS_ABSORPTION, ' ! H2 AND CO SELF-SHIELDING'
  write(10,'(a,i2,a)') 'IGRQM = ', IGRQM, ' ! 0=thermal; For H,H2: 1=QM1; 2=QM2; 3=choose fastest'
  write(10,'(a,i2,a)') 'IMODH = ', IMODH, ' ! 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only'
  write(10,'(a,i2,a)') 'ICONS = ', ICONS, ' ! 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc '
  write(10,'(a,i2,a)') 'IREAD = ', IREAD, ' ! 0 read initial abundances from the nls_control'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*    Gas phase parameters   *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'XNT0 = ', initial_gas_density, ' ! initial gas density'
  write(10,'(a,es10.3e2,a)') 'TEMP0 = ', initial_gas_temperature, ' ! initial gas temp'
  write(10,'(a,es10.3e2,a)') 'TAU0 = ', TAU0, ' ! initial visual extinction'
  write(10,'(a,es10.3e2,a)') 'ZETA0 = ', ZETA0, ' ! cosmic ray ionisation rate (1.3e-17 standard value)'
  write(10,'(a,es10.3e2,a)') 'ZETAX = ', ZETAX, ' ! Ionisation rate due to X-rays (s-1)'
  write(10,'(a,es10.3e2,a)') 'UVGAS = ', UVGAS, ' ! scale fac for UV radiation field'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*      Grain parameters     *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'DTEMP0 = ', initial_dust_temperature, ' ! initial dust temp'
  write(10,'(a,es10.3e2,a)') 'DTOGM = ', initial_dtg_mass_ratio, ' ! dust-to-gas ratio by mass'
  write(10,'(a,es10.3e2,a)') 'STICK0 = ', sticking_coeff_neutral, ' ! sticking coeff for neutral species'
  write(10,'(a,es10.3e2,a)') 'STICKP = ', sticking_coeff_positive, ' ! sticking coeff for positive species'
  write(10,'(a,es10.3e2,a)') 'STICKN = ', sticking_coeff_negative, ' ! sticking coeff for negative species'
  write(10,'(a,es10.3e2,a)') 'RHOD = ', RHOD, ' ! mass density of grain material'
  write(10,'(a,es10.3e2,a)') 'RD = ', grain_radius, ' ! grain radius (cm)'
  write(10,'(a,es10.3e2,a)') 'ACM = ', ACM, ' ! site spacing (cm)'
  write(10,'(a,es10.3e2,a)') 'SNS = ', SNS, ' ! site density (cm-2)'
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
  write(10,'(a,i4,a)') 'OTPD = ', OTPD, ' ! outputs per decade  2  8  64  128 (Only without diffusion)'
  write(10,'(a,es10.3e2,a)') 'TSTART = ', TSTART/TYEAR, ' ! first output time, after zero (yrs)'
  write(10,'(a,es10.3e2,a)') 'TFINAL = ', TFINAL/TYEAR, ' ! last output time (yrs)'
  write(10,'(a,i5,a)') 'WSTEP = ', WSTEP, ' ! Outputs every WSTEP timesteps (/=1 only for 1D outputs)'
  write(10,'(a,i5,a)') 'WSTEPR = ', WSTEPR, ' ! Outputs every WSTEPR timesteps for the rate coefficients'
  write(10,'(a,i5,a)') 'IRATEOUT = ', IRATEOUT, ' ! Spatial point for the rate output'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*     Initial Abundances    *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'XNMIN = ', XNMIN, ' ! default minimum initial frac abun'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*       1D parameters       *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i5,a)') 'IDIFF = ', IDIFF, ' ! Diffusivity flag (1. constant 2. alpha)'
  write(10,'(a,es10.3e2,a)') 'DIFFTY = ', DIFFTY, ' ! Turbulent diffusivity (cgs if 1) (value or parameter)'
  write(10,'(a,es10.3e2,a)') 'HSIZE = ', HSIZE, ' ! Computing box"s size (cm)'
  write(10,'(a,es10.3e2,a)') 'MCENTER = ', MCENTER, ' ! Central mass (g)'
  write(10,'(a,es10.3e2,a)') 'DISTR = ', DISTR, ' ! Radial distance (cm) for the disk model 1.500D+15 cm = 100 AU'
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
write(4,'(5(I4,")",1X,A11,1X))') (I,SPEC(I),I=1,nb_species)
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
!! here will have the default value XNMIN
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
integer :: boolean ! integer value used to define a logical value (a bit complicated to define directly a boolean)

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
    read(10, '(a80)', iostat=error) line
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

! Set initial abundances================================================
do I=1,nb_species
  XN(I)=XNMIN
  do j=1,nb_lines
    if (SPEC(I).EQ.temp_names(j)) then
      XN(I)=temp_abundances(j)
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
subroutine write_abundances()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'abundances.',timestep,'.out'


open(UNIT=35, file=filename_output, form='unformatted')

write(35) TIME, zspace, SPEC
write(35) TEMP1D, DEnb_species_for_gasD, TAU1D, ZETAX1D
write(35) ZXN
close(35)

return
end subroutine write_abundances

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
subroutine write_rates()

use global_variables

implicit none

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'rates.',timestep,'.out'

open(45, file=filename_output, form='unformatted')

write(45) SPEC
write(45) SYMBOL
write(45) XK
write(45) NUM

close(45)

return 
end subroutine write_rates

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
subroutine write_chemical_composition(filename)

use global_variables

implicit none

! Input
character(len=*), intent(in) :: filename !< [in] the name of the output file

! Locals
integer :: i

open(13, file=filename)
write(13,'("DEPTH POINT=",I2,"/",I2,", TIME =",1PD10.3," s",&
&", XNT=",1PD10.3," cm-3",", TEMP=",1PD10.3," K",&
&", TAU=",0PF8.3,", ZETA=",1PD10.3," s-1")') 00,00,TIME,XNT,TEMP,TAU,ZETA0

write(13,'(5(A11,":",1X,1PE12.5,2X)) ') (SPEC(I),XN(I),I=1,nb_species)
write(13,*)
close(13)

return
end subroutine write_chemical_composition




end module input_output