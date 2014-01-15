module input_output

implicit none

contains

SUBROUTINE READINPUT
! Reads the main chemical control file nls_control.d
use global_variables
use constants
implicit none
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


! Read parameters from 'gg_control.d'===================================
READ (NCON,10)
READ (NCON,14) RTOL
READ (NCON,11)
READ (NCON,12) IDUST
READ (NCON,12) ISABS
READ (NCON,12) IGRQM
READ (NCON,12) IMODH
READ (NCON,12) ICONS
READ (NCON,12) IREAD
READ (NCON,11)
READ (NCON,14) XNT0
READ (NCON,14) TEMP0
READ (NCON,14) TAU0
READ (NCON,14) ZETA0
READ (NCON,14) ZETAX
READ (NCON,14) UVGAS
READ (NCON,11)
READ (NCON,14) DTEMP0
READ (NCON,14) DTOGM
READ (NCON,14) STICK0
READ (NCON,14) STICKP
READ (NCON,14) STICKN
READ (NCON,14) RHOD
READ (NCON,14) RD
READ (NCON,14) ACM
READ (NCON,14) SNS
READ (NCON,14) EBFAC
READ (NCON,14) ACT
READ (NCON,14) TSMAX
READ (NCON,14) CRT
READ (NCON,14) CRFE
READ (NCON,14) LAYERS
READ (NCON,14) ARRK
READ (NCON,11)
READ (NCON,13) OTPD
READ (NCON,14) TSTART
READ (NCON,14) TFINAL
read(NCON,'(18X,I5)') WSTEP
read(NCON,'(18X,I5)') WSTEPR
read(NCON,'(18X,I5)') IRATEOUT
READ (NCON,11) 
READ (NCON,14) XNMIN
READ (NCON,13) NS0
READ (NCON,*) 
READ (NCON,15) (XS0(I),XN0(I),I=1,NS0)
10 FORMAT (///////)
11 FORMAT (//)
12 FORMAT (21X,I2)
13 FORMAT (19X,I4)
14 FORMAT (11X,E12.6)
15 FORMAT (A11,3X,E12.6)

TSTART=TSTART*TYEAR
TFINAL=TFINAL*TYEAR

! read 1D and jacobian parameters

read(NCON,'(///)')
read(NCON,'(18X,I5)') IDIFF
read(NCON,'(11X,D12.6)') DIFFTY
read(NCON,'(11X,D12.6)') HSIZE
read(NCON,'(11X,D12.6)') MCENTER
read(NCON,'(11X,D12.6)') DISTR
read(NCON,'(18X,I5)') TESTJAC
read(NCON,'(18X,I5)') NJAC

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
  print *,'Some reaction was not found by the reorder process'
  print *,jk,'=/',nkmax+1 
  stop
endif

!       replace the species names by blanks for non chemical species                                                                        
do j=1,nkmax-1
  if ((SYMBOL(1,J).eq.'CR         ').or.(SYMBOL(1,J).eq.'CRP        ').or.&
  & (SYMBOL(1,J).eq.'Photon     ')) SYMBOL(1,J)='            '
  if ((SYMBOL(2,J).eq.'CR         ').or.(SYMBOL(2,J).eq.'CRP        ').or.&
  & (SYMBOL(2,J).eq.'Photon     ')) SYMBOL(2,J)='           '
  if ((SYMBOL(3,J).eq.'CR         ').or.(SYMBOL(3,J).eq.'CRP        ').or.&
  & (SYMBOL(3,J).eq.'Photon     ')) SYMBOL(3,J)='           '
  if ((SYMBOL(4,J).eq.'CR         ').or.(SYMBOL(4,J).eq.'CRP        ').or.&
  & (SYMBOL(4,J).eq.'Photon     ')) SYMBOL(4,J)='           '
  if ((SYMBOL(5,J).eq.'CR         ').or.(SYMBOL(5,J).eq.'CRP        ').or.&
  & (SYMBOL(5,J).eq.'Photon     ')) SYMBOL(5,J)='           '
  if ((SYMBOL(6,J).eq.'CR         ').or.(SYMBOL(6,J).eq.'CRP        ').or.&
  & (SYMBOL(6,J).eq.'Photon     ')) SYMBOL(6,J)='           '
  if ((SYMBOL(7,J).eq.'CR         ').or.(SYMBOL(7,J).eq.'CRP        ').or.&
  & (SYMBOL(7,J).eq.'Photon     ')) SYMBOL(7,J)='           '

ENDDO 


RETURN
END SUBROUTINE READINPUT

! ======================================================================
! ======================================================================
SUBROUTINE WRITESPEC
use global_variables
use constants
implicit none

integer :: i

! Write 'ggo_spec.d': 5 columns of numbered species=====================
WRITE (NSP,270) (I,SPEC(I),I=1,NSMAX)
270 FORMAT (5(I4,')',1X,A11,1X))

RETURN
END SUBROUTINE WRITESPEC

! ======================================================================
! ======================================================================
subroutine write1D
! Writes 1D outputs
use global_variables
use constants
implicit none
character (len=6) :: charit
integer :: i

write(CHARIT,'(I6)') IT

do i=1,6
  if (CHARIT(i:i).eq." ") CHARIT(i:i)="0"
enddo

open(UNIT=NOUT2,file='output_1D.'//CHARIT//'', form='unformatted')

write(NOUT2) TIME, zspace, SPEC
write(NOUT2) TEMP1D, DENS1D, TAU1D, ZETAX1D
write(NOUT2) ZXN

close(NOUT2)

return
end subroutine write1D

! ======================================================================
! ======================================================================
subroutine rates1D
! Writes rate coefficient for a particular mesh point
use global_variables
use constants
implicit none
character (len=6) :: charit
integer :: i

write(CHARIT,'(I6)') IT

do i=1,6
  if (CHARIT(i:i).eq." ") CHARIT(i:i)="0"
enddo

open(NORD2,file='rates1D.'//CHARIT//'',form='unformatted')

write(NORD2) SPEC
write(NORD2) SYMBOL
write(NORD2) XK
write(NORD2) NUM

close(NORD2)

return 
end subroutine rates1D

! ======================================================================
! ======================================================================
subroutine writetail
! Writes the final output for 0D runs in readable format
use global_variables
use constants
implicit none
integer :: i

WRITE (NTAI,690) 00,00,TIME,XNT,TEMP,TAU,ZETA0
690 FORMAT ('DEPTH POINT=',I2,'/',I2,', TIME =',1PD10.3,' s',&
', XNT=',1PD10.3,' cm-3',', TEMP=',1PD10.3,' K',&
', TAU=',0PF8.3,', ZETA=',1PD10.3,' s-1')

WRITE (NTAI,700) (SPEC(I),XN(I),I=1,NSMAX)
700 FORMAT (5(A11,':',1X,1PE12.5,2X)) 
WRITE (NTAI,*)

RETURN
END subroutine writetail

! ======================================================================
! ======================================================================

end module input_output