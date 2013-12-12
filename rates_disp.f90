program taux
implicit none
integer, parameter :: NKMAX=8335
integer, parameter :: NSMAX=684
integer, parameter :: nptmax=1
integer, parameter :: IPT=1
integer, parameter :: NT=29
integer, parameter :: WSTEP=1
real(kind=8), dimension(NKMAX)::XK
character (len=11), dimension(7,nkmax) :: SYMBOL
CHARACTER (len=11), dimension(nsmax) :: SPEC
real(kind=8), dimension(nptmax):: TEMP1D, DENS1D, TAU1D 
real(kind=8), dimension(nsmax,nptmax) ::  ZXN
real(kind=8), dimension(0:nsmax,nptmax) ::  ZXN2
real(kind=8) :: time
integer, dimension(NKMAX,7) :: REACT
real(kind=8) :: destr,prod, XNT
integer :: i,j,k,l
integer :: IT, KIT
character (len=6) :: charit
INTEGER, dimension(nkmax):: NUM

do KIT=WSTEP,NT,WSTEP

if (KIT.eq.0) then
IT=1
else
IT=KIT
endif

write(CHARIT,'(I6)') IT

do i=1,6
if (CHARIT(i:i).eq." ") CHARIT(i:i)="0"
enddo

open(10,file='rates1D.'//CHARIT//'',form='unformatted')
read(10) SPEC
read(10) SYMBOL
read(10) XK
read(10) NUM

open(21,file='rate_coeff_modif.dat',form='formatted')
do i=1,NKMAX
	write(21,*) (SYMBOL(j,i), j=1,7),XK(i),NUM(i)
enddo
close(21)	

open(30,file='rates.d.'//CHARIT//'',form='formatted')

REACT(:,:)=0

        DO I=1,NKMAX
                DO J=1,NSMAX
                DO L=1,3
                        IF (SYMBOL(L,I).EQ.SPEC(J)) REACT(I,L)=J
                ENDDO
                DO L=1,4
                        IF (SYMBOL(L+3,I).EQ.SPEC(J)) REACT(I,L+3)=J
                ENDDO
                ENDDO
        ENDDO

open(20,file='output_1D.'//CHARIT//'',form='unformatted')
read(20) TIME
read(20) TEMP1D, DENS1D, TAU1D
read(20) ZXN

XNT=2.*DENS1D(IPT)

ZXN2(1:nsmax,1:nptmax)=ZXN(1:nsmax,1:nptmax)
ZXN2(0,1:nptmax)=1.d0/XNT

do j=1,nkmax
do i=1,3
!if (symbol(i,j).eq.'CO      ') then
if (symbol(i,j).eq.'JN2        ') then
destr=XK(j)*ZXN2(react(j,1),IPT)*ZXN2(react(j,2),IPT)*ZXN2(react(j,3),IPT)*XNT**3
write(30,('(7A9,1pd10.3,1pd10.3)')) (symbol(l,j),l=1,7),XK(j), destr
endif
enddo
do i=4,7
!if (symbol(i,j).eq.'CO      ') then
if (symbol(i,j).eq.'JN2        ') then
prod=XK(j)*ZXN2(react(j,1),IPT)*ZXN2(react(j,2),IPT)*ZXN2(react(j,3),IPT)*XNT**3
write(30,('(8A9,1pd10.3,1pd10.3)')) (symbol(l,j),l=1,7), "       ",XK(j),prod
endif
enddo
enddo

close(10)
close(20)
close(30)

enddo

end
