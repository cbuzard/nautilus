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
subroutine read_parameters()

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



open(5, file='parameters.in')

read(5,10)
read(5,14) RELATIVE_TOLERANCE
read(5,11)
read(5,12) IS_GRAIN_REACTIONS
read(5,12) ISABS
read(5,12) IGRQM
read(5,12) IMODH
read(5,12) ICONS
read(5,12) IREAD
read(5,11)
read(5,14) initial_gas_density
read(5,14) initial_gas_temperature
read(5,14) TAU0
read(5,14) ZETA0
read(5,14) ZETAX
read(5,14) UVGAS
read(5,11)
read(5,14) initial_dust_temperature
read(5,14) initial_dtg_mass_ratio
read(5,14) sticking_coeff_neutral
read(5,14) sticking_coeff_positive
read(5,14) sticking_coeff_negative
read(5,14) RHOD
read(5,14) grain_radius
read(5,14) ACM
read(5,14) SNS
read(5,14) EBFAC
read(5,14) ACT
read(5,14) TSMAX
read(5,14) CRT
read(5,14) CRFE
read(5,14) LAYERS
read(5,14) ARRK
read(5,11)
read(5,13) OTPD
read(5,14) TSTART
read(5,14) TFINAL
read(5,'(18X,I5)') WSTEP
read(5,'(18X,I5)') WSTEPR
read(5,'(18X,I5)') IRATEOUT
read(5,11) 
read(5,14) XNMIN
read(5,13) NS0
read(5,*) 
10 format(///////)
11 format(//)
12 format(21X,I2)
13 format(19X,I4)
14 format(11X,E12.6)

TSTART=TSTART*TYEAR
TFINAL=TFINAL*TYEAR

! read 1D and jacobian parameters

read(5,'(///)')
read(5,'(18X,I5)') IDIFF
read(5,'(11X,D12.6)') DIFFTY
read(5,'(11X,D12.6)') HSIZE
read(5,'(11X,D12.6)') MCENTER
read(5,'(11X,D12.6)') DISTR
read(5,'(18X,I5)') TESTJAC
read(5,'(18X,I5)') NJAC

close(5)

call read_abundances()

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
end subroutine read_parameters

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
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write all abundances for all species in an output file at each
!! output time. The total number of output files related to abundances 
!! will be equal to the number of timestep, not equally spaced in time.\n\n
!! Output filename is of the form : abundances.000001.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_abundances()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80) :: filename='abundances.in'
integer :: i, j

real(double_precision), allocatable, dimension(:) :: temp_abundances
character(len=11), allocatable, dimension(:) :: temp_names

call get_linenumber(filename, NS0)

allocate(temp_abundances(NS0))
allocate(temp_names(NS0))

open(5, file=filename)
read(5,15) (temp_names(I),temp_abundances(I),I=1,NS0)
close(5)
15 format(A11,3X,E12.6)

XS0(1:NS0) = temp_names(1:NS0)
XN0(1:NS0) = temp_abundances(1:NS0)

!~ ! Set initial abundances================================================
!~ do I=1,nb_species
!~   XN(I)=XNMIN
!~   do j=1,NS0
!~     if (SPEC(I).EQ.XS0(j)) then
!~       XN(I)=XN0(j)
!~     endif
!~   enddo
!~ enddo

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