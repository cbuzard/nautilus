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
SUBROUTINE read_parameters()

use global_variables
use constants

implicit none

! Locals
integer :: i,k,j, jk

! Variables for the unordered reaction file

CHARACTER (len=11), dimension(ns1) :: SPECUO1
CHARACTER (len=11), dimension(ns2) :: SPECUO2
character (len=11), dimension(7,nkmax) :: SYMBOLUO
real(double_precision), dimension(nkmax) :: AUO,BUO,CUO
integer, dimension(nkmax) :: itypeUO,TminUO,TmaxUO,FORMULAUO,NUMUO
character (len=11), dimension(7,nk1) :: SYMBOLUO1
real(double_precision), dimension(nk1) :: AUO1,BUO1,CUO1
integer, dimension(nk1) :: itypeUO1,Tmin1,Tmax1,FORMULA1,NUM1
character (len=11), dimension(7,nk2) :: SYMBOLUO2
real(double_precision), dimension(nk2) :: AUO2,BUO2,CUO2
integer, dimension(nk2) :: itypeUO2,Tmin2,Tmax2,FORMULA2,NUM2
integer, dimension(ns1) :: ICG1
integer, dimension(ns2) :: ICG2
integer, dimension(nemax,ns1) :: IELM1
integer, dimension(nemax,ns2) :: IELM2

open(5, file='parameters.in')

READ(5,10)
READ(5,14) RTOL
READ(5,11)
READ(5,12) IDUST
READ(5,12) ISABS
READ(5,12) IGRQM
READ(5,12) IMODH
READ(5,12) ICONS
READ(5,12) IREAD
READ(5,11)
READ(5,14) XNT0
READ(5,14) TEMP0
READ(5,14) TAU0
READ(5,14) ZETA0
READ(5,14) ZETAX
READ(5,14) UVGAS
READ(5,11)
READ(5,14) DTEMP0
READ(5,14) DTOGM
READ(5,14) STICK0
READ(5,14) STICKP
READ(5,14) STICKN
READ(5,14) RHOD
READ(5,14) RD
READ(5,14) ACM
READ(5,14) SNS
READ(5,14) EBFAC
READ(5,14) ACT
READ(5,14) TSMAX
READ(5,14) CRT
READ(5,14) CRFE
READ(5,14) LAYERS
READ(5,14) ARRK
READ(5,11)
READ(5,13) OTPD
READ(5,14) TSTART
READ(5,14) TFINAL
read(5,'(18X,I5)') WSTEP
read(5,'(18X,I5)') WSTEPR
read(5,'(18X,I5)') IRATEOUT
READ(5,11) 
READ(5,14) XNMIN
READ(5,13) NS0
READ(5,*) 
READ(5,15) (XS0(I),XN0(I),I=1,NS0)
10 FORMAT (///////)
11 FORMAT (//)
12 FORMAT (21X,I2)
13 FORMAT (19X,I4)
14 FORMAT (11X,E12.6)
15 FORMAT (A11,3X,E12.6)

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

! Read CO and H2 shielding factors=====================

DO I=1,NL1
  READ (H2DIS,81) N1H2(I),T1H2(I)
ENDDO
81        FORMAT (E9.3,3X,E9.3)

DO I=1,NL3
  READ (CODIS,82) N2CO(I),T2CO(I),N2H2(I),T2H2(I),&
  AV2(I),T2AV(I)              
ENDDO
82        FORMAT (6(E9.3,3X))

DO I=NL3+1,NL2
  READ (CODIS,81) N2CO(I),T2CO(I)
ENDDO


! Read species & reaction info from reactions file======================
! WV fev 2012
! There are now two different files in which the reactions and species are

! Reading the gas phase network

DO I=1,NS1
  READ (NJR,80) SPECUO1(I),ICG1(I),(IELM1(K,I),K=1,NEMAX) 
  80    FORMAT (A11,i3,13(I3)) 
ENDDO

READ (NJR,90) ((SYMBOLUO1(I,J),I=1,7),AUO1(J),BUO1(J),CUO1(J), &
ITYPEUO1(J),Tmin1(j),Tmax1(j),FORMULA1(J),NUM1(J),J=1,NK1) 
90 FORMAT (3A11,1x,4A11,11x,3D11.3,23x,I3,2i7,i3,i6)

! Reading the grain network

DO I=1,NS2
  READ (NJR2,80) SPECUO2(I),ICG2(I),(IELM2(K,I),K=1,NEMAX) 
ENDDO

READ (NJR2,90) ((SYMBOLUO2(I,J),I=1,7),AUO2(J),BUO2(J),CUO2(J), &
ITYPEUO2(J),Tmin2(j),Tmax2(j),FORMULA2(J),NUM2(J),J=1,NK2) 
! putting everything back into the big tables

DO I=1,NS1 
  SPEC(I)=SPECUO1(I)
  ICG(I)=ICG1(I)
  do k=1,NEMAX
    IELM(K,I)=IELM1(K,I)
  enddo
ENDDO
DO I=1,NS2 
  SPEC(NS1+I)=SPECUO2(I)
  ICG(NS1+I)=ICG2(I)
  do k=1,NEMAX
    IELM(K,NS1+I)=IELM2(K,I)
  enddo
ENDDO

DO I=1,NK1 
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
ENDDO


DO I=1,NK2 
  do  k=1,7
    SYMBOLUO(k,NK1+I)=SYMBOLUO2(k,I)
  enddo
  AUO(NK1+I)=AUO2(I)
  BUO(NK1+I)=BUO2(I)
  CUO(NK1+I)=CUO2(I)
  ITYPEUO(NK1+I)=ITYPEUO2(I)
  TminUO(NK1+I) = Tmin2(I)
  TmaxUO(NK1+I) = Tmax2(I)
  FORMULAUO(NK1+I) = FORMULA2(I)
  NUMUO(NK1+I) = NUM2(I)
ENDDO

! Reorder reaction file entries with ITYPE
jk=1
do i=0,nitype
  do j=1,nkmax
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


if (jk.ne.nkmax+1) then
  write(*,*) 'Some reaction was not found by the reorder process'
  write(*,*) jk,'=/',nkmax+1 
  stop
endif

!       replace the species names by blanks for non chemical species                                                                        
do j=1,nkmax-1
  do i=1,7
    select case(symbol(i,j))
      case ('CR', 'CRP', 'Photon')
        symbol(i,j) = '           '
    end select
  enddo

ENDDO 


RETURN
END SUBROUTINE read_parameters

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
SUBROUTINE write_species()
use global_variables
use constants

implicit none

! Locals
integer :: i

open(4, file='species.out')
! Write 'ggo_spec.d': 5 columns of numbered species=====================
WRITE (4,270) (I,SPEC(I),I=1,NSMAX)
270 FORMAT (5(I4,')',1X,A11,1X))
close(4)

RETURN
END SUBROUTINE write_species

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
use constants
implicit none
character(len=80) :: filename_output
integer :: i

write(filename_output, '(a,i0.6,a)') 'abundances.',IT,'.out'


open(UNIT=35, file=filename_output, form='unformatted')

write(35) TIME, zspace, SPEC
write(35) TEMP1D, DENS1D, TAU1D, ZETAX1D
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
use constants

implicit none

! Locals
character(len=80) :: filename_output
integer :: i

write(filename_output, '(a,i0.6,a)') 'rates.',IT,'.out'

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
use constants

implicit none

! Input
character(len=*), intent(in) :: filename !< [in] the name of the output file

! Locals
integer :: i

open(13, file=filename)
WRITE (13,690) 00,00,TIME,XNT,TEMP,TAU,ZETA0
690 FORMAT ('DEPTH POINT=',I2,'/',I2,', TIME =',1PD10.3,' s',&
', XNT=',1PD10.3,' cm-3',', TEMP=',1PD10.3,' K',&
', TAU=',0PF8.3,', ZETA=',1PD10.3,' s-1')

WRITE (13,700) (SPEC(I),XN(I),I=1,NSMAX)
700 FORMAT (5(A11,':',1X,1PE12.5,2X)) 
WRITE (13,*)
close(13)

RETURN
END subroutine write_chemical_composition

end module input_output