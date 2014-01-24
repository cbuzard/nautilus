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

! H2 Shielding factors
! NH2 (cm-2)
N1H2 = (/0.000d+00, 3.690d+11, 3.715d+12, 3.948d+13, 1.233d+14, 2.536d+14, 4.342d+14, 6.653d+14, 6.689d+14, &
9.075d+14, 1.234d+15, 1.631d+15, 2.105d+15, 2.363d+15, 2.899d+15, 3.207d+15, 3.848d+15, 4.636d+15, &
5.547d+15, 6.604d+15, 7.855d+15, 9.368d+15, 1.122d+16, 1.352d+16, 1.643d+16, 2.017d+16, 2.515d+16, &
3.190d+16, 4.128d+16, 5.439d+16, 7.315d+16, 1.009d+17, 1.432d+17, 2.092d+17, 3.123d+17, 4.738d+17, &
5.388d+17, 8.935d+17, 1.381d+18, 2.164d+18, 3.330d+18, 5.024d+18, 7.404d+18, 9.029d+18, 1.316d+19, &
1.813d+19, 2.453d+19, 3.248d+19, 4.216d+19, 5.370d+19, 6.722d+19, 8.277d+19, 9.894d+19, 1.186d+20, &
1.404d+20, 1.644d+20, 1.908d+20, 2.197d+20, 2.510d+20, 2.849d+20, 3.214d+20, 3.604d+20, 4.019d+20, &
4.456d+20, 4.915d+20, 5.393d+20, 5.886d+20, 6.392d+20, 6.909d+20, 8.505d+20, 7.433d+20, 7.965d+20, &
9.056d+20, 9.627d+20, 1.011d+21, 1.068d+21, 1.125d+21, 1.185d+21, 1.250d+21, 1.327d+21, 1.428d+21, &
1.578d+21, 1.851d+21, 2.128d+21, 2.298d+21, 2.389d+21, 2.459d+21, 2.519d+21, 2.571d+21, 2.618d+21, &
2.707d+21, 2.790d+21, 2.887d+21, 3.001d+21, 3.139d+21, 3.303d+21, 3.497d+21, 3.722d+21, 3.983d+21, &
4.283d+21, 4.644d+21, 5.127d+21, 5.945d+21, 8.205d+21, 1.015d+22/)

! teta H2
T1H2 = (/1.000d+00, 9.983d-01, 9.853d-01, 8.761d-01, 7.199d-01, 5.728d-01, 4.455d-01, 3.431d-01, 3.418d-01, &
2.732d-01, 2.110d-01, 1.619d-01, 1.236d-01, 1.084d-01, 8.447d-02, 7.410d-02, 5.774d-02, 4.416d-02, &
3.390d-02, 2.625d-02, 2.048d-02, 1.606d-02, 1.264d-02, 9.987d-03, 7.937d-03, 6.343d-03, 5.088d-03, &
4.089d-03, 3.283d-03, 2.640d-03, 2.130d-03, 1.725d-03, 1.397d-03, 1.129d-03, 9.097d-04, 7.340d-04, &
6.883d-04, 5.377d-04, 4.352d-04, 3.475d-04, 2.771d-04, 2.205d-04, 1.753d-04, 1.549d-04, 1.210d-04, &
9.666d-05, 7.705d-05, 6.148d-05, 4.904d-05, 3.909d-05, 3.112d-05, 2.473d-05, 1.997d-05, 1.578d-05, &
1.244d-05, 9.769d-06, 7.634d-06, 5.932d-06, 4.581d-06, 3.515d-06, 2.679d-06, 2.029d-06, 1.527d-06, &
1.144d-06, 8.523d-07, 6.332d-07, 4.693d-07, 3.475d-07, 2.574d-07, 1.047d-07, 1.907d-07, 1.413d-07, &
7.739d-08, 5.677d-08, 4.386d-08, 3.227d-08, 2.385d-08, 1.750d-08, 1.248d-08, 8.389d-09, 5.026d-09, &
2.382d-09, 6.259d-10, 1.653d-10, 7.399d-11, 4.824d-11, 3.474d-11, 2.633d-11, 2.069d-11, 1.663d-11, &
1.099d-11, 7.506d-12, 4.825d-12, 2.864d-12, 1.534d-12, 7.324d-13, 3.087d-13, 1.135d-13, 3.591d-14, &
9.689d-15, 2.045d-15, 2.618d-16, 8.918d-18, 3.041d-21, 1.739d-23/)

! CO shidlding factors
! NCO (cm-2)
N2CO = (/ 0.000d+00, 1.000d+12, 1.650d+12, 2.995d+12, 5.979d+12, 1.313d+13, 3.172d+13, 8.429d+13, 2.464d+14, &
7.923d+14, 1.670d+15, 2.595d+15, 4.435d+15, 6.008d+15, 8.952d+15, 1.334d+16, 1.661d+16, 2.274d+16, &
3.115d+16, 4.266d+16, 5.843d+16, 8.002d+16, 1.096d+17, 1.501d+17, 2.055d+17, 2.815d+17, 4.241d+17, &
6.389d+17, 9.625d+17, 1.450d+18, 2.184d+18, 3.291d+18, 4.124d+18, 5.685d+18, 7.838d+18, 1.080d+19, &
1.285d+19, 1.681d+19, 2.199d+19, 2.538d+19, 3.222d+19, 4.091d+19, 5.193d+19, 5.893d+19, 7.356d+19, &
8.269d+19, 9.246d+19, 1.031d+20, 1.148d+20, 1.277d+20, 1.419d+20, 1.578d+20/)

! Teta CO
T2CO = (/ 1.000d+00, 9.990d-01, 9.981d-01, 9.961d-01, 9.912d-01, 9.815d-01, 9.601d-01, 9.113d-01, 8.094d-01, &
6.284d-01, 4.808d-01, 3.889d-01, 2.827d-01, 2.293d-01, 1.695d-01, 1.224d-01, 1.017d-01, 7.764d-02, &
5.931d-02, 4.546d-02, 3.506d-02, 2.728d-02, 2.143d-02, 1.700d-02, 1.360d-02, 1.094d-02, 8.273d-03, &
6.283d-03, 4.773d-03, 3.611d-03, 2.704d-03, 1.986d-03, 1.657d-03, 1.258d-03, 9.332d-04, 6.745d-04, &
5.596d-04, 4.123d-04, 2.982d-04, 2.490d-04, 1.827d-04, 1.324d-04, 9.473d-05, 7.891d-05, 5.668d-05, &
4.732d-05, 3.967d-05, 3.327d-05, 2.788d-05, 2.331d-05, 1.944d-05, 1.619d-05/)

!NH2 (cm-2)
N2H2 = (/0.000d+00, 2.666d+13, 3.801d+14, 6.634d+15, 8.829d+16, 9.268d+17, 1.007d+18, 2.021d+18, 3.036d+18, &
4.051d+18, 5.066d+18, 6.082d+18, 7.097d+18, 8.112d+18, 9.341d+18, 1.014d+19, 2.030d+19, 3.045d+19, &
4.061d+19, 5.076d+19, 6.092d+19, 7.107d+19, 8.123d+19, 9.353d+19, 1.015d+20, 2.031d+20, 3.047d+20, &
4.062d+20, 5.078d+20, 6.094d+20, 7.109d+20, 8.125d+20, 9.355d+20, 1.016d+21, 2.031d+21, 3.047d+21, &
4.063d+21, 5.078d+21, 6.094d+21, 7.110d+21, 8.125d+21, 9.355d+21, 1.016d+22/)

! Teta H2
T2H2 = (/1.000d+00, 9.999d-01, 9.893d-01, 9.678d-01, 9.465d-01, 9.137d-01, 9.121d-01, 8.966d-01, 8.862d-01, &
8.781d-01, 8.716d-01, 8.660d-01, 8.612d-01, 8.569d-01, 8.524d-01, 8.497d-01, 8.262d-01, 8.118d-01, &
8.010d-01, 7.921d-01, 7.841d-01, 7.769d-01, 7.702d-01, 7.626d-01, 7.579d-01, 7.094d-01, 6.712d-01, &
6.378d-01, 6.074d-01, 5.791d-01, 5.524d-01, 5.271d-01, 4.977d-01, 4.793d-01, 2.837d-01, 1.526d-01, &
7.774d-02, 3.952d-02, 2.093d-02, 1.199d-02, 7.666d-03, 5.333d-03, 4.666d-03/)

! AV
AV2 = (/0.000d+00, 1.000d-07, 1.000d-06, 1.000d-05, 1.000d-04, 1.000d-03, 1.086d-03, 2.171d-03, 3.257d-03, &
4.343d-03, 5.429d-03, 6.514d-03, 7.600d-03, 8.686d-03, 1.000d-02, 1.086d-02, 2.171d-02, 3.257d-02, &
4.343d-02, 5.429d-02, 6.514d-02, 7.600d-02, 8.686d-02, 1.000d-01, 1.086d-01, 2.171d-01, 3.257d-01, &
4.343d-01, 5.429d-01, 6.514d-01, 7.600d-01, 8.686d-01, 1.000d+00, 1.086d+00, 2.171d+00, 3.257d+00, &
4.343d+00, 5.429d+00, 6.514d+00, 7.600d+00, 8.686d+00, 1.000d+01, 1.086d+01/)

! teta AV
T2AV = (/1.000d+00, 1.000d+00, 1.000d+00, 1.000d+00, 9.991d-01, 9.887d-01, 9.877d-01, 9.755d-01, 9.638d-01, &
9.524d-01, 9.413d-01, 9.306d-01, 9.202d-01, 9.101d-01, 8.983d-01, 8.908d-01, 8.081d-01, 7.422d-01, &
6.875d-01, 6.404d-01, 5.990d-01, 5.620d-01, 5.284d-01, 4.916d-01, 4.696d-01, 2.757d-01, 1.705d-01, &
1.084d-01, 7.013d-02, 4.593d-02, 3.037d-02, 2.023d-02, 1.247d-02, 9.130d-03, 2.133d-04, 6.269d-06, &
2.108d-07, 7.633d-09, 2.873d-10, 1.102d-11, 4.273d-13, 8.412d-15, 6.510d-16/)

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
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction
integer :: boolean !< integer value used to define a logical value (a bit complicated to define directly a boolean)

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
      case('is_grain_reactions')
        read(value, '(i2)') IS_GRAIN_REACTIONS
      
      case('is_absorption')
        read(value, '(i2)') IS_ABSORPTION
      
      case('IGRQM')
        read(value, '(i2)') IGRQM
      
      case('IMODH')
        read(value, '(i2)') IMODH
      
      case('ICONS')
        read(value, '(i2)') ICONS
      
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
      
      case('UVGAS')
        read(value, '(e12.6)') UVGAS
      
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
  write(10,'(a,i2,a)') 'IGRQM = ', IGRQM, ' ! 0=thermal; For H,H2: 1=QM1; 2=QM2; 3=choose fastest'
  write(10,'(a,i2,a)') 'IMODH = ', IMODH, ' ! 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only'
  write(10,'(a,i2,a)') 'ICONS = ', ICONS, ' ! 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc '
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*    Gas phase parameters   *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'initial_gas_density = ', initial_gas_density, ' ! initial gas density'
  write(10,'(a,es10.3e2,a)') 'initial_gas_temperature = ', initial_gas_temperature, ' ! initial gas temp'
  write(10,'(a,es10.3e2,a)') 'initial_visual_extinction = ', INITIAL_VISUAL_EXTINCTION, ' ! initial visual extinction'
  write(10,'(a,es10.3e2,a)') 'cr_ionisation_rate = ', CR_IONISATION_RATE, ' ! cosmic ray ionisation rate (1.3e-17 standard value)'
  write(10,'(a,es10.3e2,a)') 'x_ionisation_rate = ', X_IONISATION_RATE, ' ! Ionisation rate due to X-rays (s-1)'
  write(10,'(a,es10.3e2,a)') 'UVGAS = ', UVGAS, ' ! scale fac for UV radiation field'
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
  write(10,'(a,i4,a)') 'output_per_decade = ', OUTPUT_PER_DECADE, ' ! outputs per decade  2  8  64  128 (Only without diffusion)'
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
  XN(I)=MINIMUM_INITIAL_ABUNDANCE
  do j=1,nb_lines
    if (species_name(I).EQ.temp_names(j)) then
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
  write(13,'(a," = ",1PE12.5)') species_name(I),XN(I)
enddo

write(13,*)
close(13)

return
end subroutine write_abundances




end module input_output