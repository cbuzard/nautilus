subroutine CHEMSETUP
! Initialize the reactants and products
use header
implicit none
integer :: I,J,L,NSP1

! NSP1 is a dummy species (with a blank name)

NSP1=NSMAX+1

SPEC2(1:NSMAX)=SPEC(1:NSMAX)
SPEC2(NSP1)='           '

DO I=1,NKMAX
   DO J=1,NSP1

      DO L=1,3
      IF (SYMBOL(L,I).EQ.SPEC2(J)) REACT(I,L)=J
      ENDDO

      DO L=1,4
      IF (SYMBOL(L+3,I).EQ.SPEC2(J)) REACT(I,L+3)=J
      ENDDO           

   ENDDO
ENDDO   

return
end

subroutine FCHEMVW(N,Y,YDOT)
! Computes the chemical evolution
use header
implicit none
integer :: N, NSP1
REAL(KIND=8), dimension(NSMAX) :: Y, YDOT
REAL(KIND=8), dimension(NSMAX+1) :: YD2
!REAL(KIND=16), dimension(NSMAX+1) :: YD2
integer :: i
integer :: IR1, IR2, IR3, IPROD1, IPROD2, IPROD3, IPROD4
real(kind=8) :: rate

NSP1=NSMAX+1

ydot(:)=0.d0
yd2(:)=0.d0

! The differential equations are calculated in a loop here
DO I=1,NKMAX

IR1=REACT(I,1)
IR2=REACT(I,2)
IR3=REACT(I,3)

IPROD1=REACT(I,4)
IPROD2=REACT(I,5)
IPROD3=REACT(I,6)
IPROD4=REACT(I,7)
                
if (IR3.ne.NSP1) then
RATE=XK(I)*Y(IR1)*Y(IR2)*Y(IR3)*XNT*XNT
endif
                
if ((IR3.eq.NSP1).and.(IR2.ne.NSP1)) then
RATE=XK(I)*Y(IR1)*Y(IR2)*XNT
endif
        
if (IR2.eq.NSP1) then
RATE=XK(I)*Y(IR1)  
endif

YD2(IPROD1)=YD2(IPROD1)+RATE
YD2(IPROD2)=YD2(IPROD2)+RATE
YD2(IPROD3)=YD2(IPROD3)+RATE
YD2(IPROD4)=YD2(IPROD4)+RATE
YD2(IR1)=YD2(IR1)-RATE
YD2(IR2)=YD2(IR2)-RATE
YD2(IR3)=YD2(IR3)-RATE
ENDDO   

YDOT(1:NSMAX)=YD2(1:NSMAX)

return
end

subroutine JACVW(Y,J,PDJ)
! Computes columns of the chemical jacobian
use header
implicit none
integer :: N, NSP1
integer J
REAL(KIND=8), dimension(NSMAX) :: Y, PDJ
REAL(KIND=8), dimension(NSMAX+1) :: PDJ2
!REAL(KIND=16), dimension(NSMAX+1) :: PDJ2
integer :: i
integer :: IR1, IR2, IR3, IPROD1, IPROD2, IPROD3, IPROD4
integer :: NUMBERJAC

NSP1=NSMAX+1

PDJ2(:)=0.d0

DO I=1,NKMAX
!write(*,*) I

  IR1=REACT(I,1)
  IR2=REACT(I,2)
  IR3=REACT(I,3)
  IPROD1=REACT(I,4)
  IPROD2=REACT(I,5)
  IPROD3=REACT(I,6)
  IPROD4=REACT(I,7)

if (IR3.ne.NSP1) then

 if (IR1.eq.J) then 
 PDJ2(IPROD1)=PDJ2(IPROD1)+XK(I)*Y(IR2)*Y(IR3)*XNT*XNT
 PDJ2(IPROD2)=PDJ2(IPROD2)+XK(I)*Y(IR2)*Y(IR3)*XNT*XNT
 PDJ2(IPROD3)=PDJ2(IPROD3)+XK(I)*Y(IR2)*Y(IR3)*XNT*XNT
 PDJ2(IPROD4)=PDJ2(IPROD4)+XK(I)*Y(IR2)*Y(IR3)*XNT*XNT
 PDJ2(IR1)=PDJ2(IR1)-XK(I)*Y(IR2)*Y(IR3)*XNT*XNT
 PDJ2(IR2)=PDJ2(IR2)-XK(I)*Y(IR2)*Y(IR3)*XNT*XNT
 PDJ2(IR3)=PDJ2(IR3)-XK(I)*Y(IR2)*Y(IR3)*XNT*XNT
 endif

 if (IR2.eq.J) then 
 PDJ2(IPROD1)=PDJ2(IPROD1)+XK(I)*Y(IR1)*Y(IR3)*XNT*XNT
 PDJ2(IPROD2)=PDJ2(IPROD2)+XK(I)*Y(IR1)*Y(IR3)*XNT*XNT
 PDJ2(IPROD3)=PDJ2(IPROD3)+XK(I)*Y(IR1)*Y(IR3)*XNT*XNT
 PDJ2(IPROD4)=PDJ2(IPROD4)+XK(I)*Y(IR1)*Y(IR3)*XNT*XNT
 PDJ2(IR1)=PDJ2(IR1)-XK(I)*Y(IR1)*Y(IR3)*XNT*XNT
 PDJ2(IR2)=PDJ2(IR2)-XK(I)*Y(IR1)*Y(IR3)*XNT*XNT
 PDJ2(IR3)=PDJ2(IR3)-XK(I)*Y(IR1)*Y(IR3)*XNT*XNT
 endif

 if (IR3.eq.J) then 
 PDJ2(IPROD1)=PDJ2(IPROD1)+XK(I)*Y(IR1)*Y(IR2)*XNT*XNT
 PDJ2(IPROD2)=PDJ2(IPROD2)+XK(I)*Y(IR1)*Y(IR2)*XNT*XNT
 PDJ2(IPROD3)=PDJ2(IPROD3)+XK(I)*Y(IR1)*Y(IR2)*XNT*XNT
 PDJ2(IPROD4)=PDJ2(IPROD4)+XK(I)*Y(IR1)*Y(IR2)*XNT*XNT
 PDJ2(IR1)=PDJ2(IR1)-XK(I)*Y(IR1)*Y(IR2)*XNT*XNT
 PDJ2(IR2)=PDJ2(IR2)-XK(I)*Y(IR1)*Y(IR2)*XNT*XNT
 PDJ2(IR3)=PDJ2(IR3)-XK(I)*Y(IR1)*Y(IR2)*XNT*XNT
 endif

endif

if ((IR3.eq.NSP1).and.(IR2.ne.NSP1)) then

 if (IR1.eq.J) then 
 PDJ2(IPROD1)=PDJ2(IPROD1)+XK(I)*Y(IR2)*XNT
 PDJ2(IPROD2)=PDJ2(IPROD2)+XK(I)*Y(IR2)*XNT
 PDJ2(IPROD3)=PDJ2(IPROD3)+XK(I)*Y(IR2)*XNT
 PDJ2(IPROD4)=PDJ2(IPROD4)+XK(I)*Y(IR2)*XNT
 PDJ2(IR1)=PDJ2(IR1)-XK(I)*Y(IR2)*XNT
 PDJ2(IR2)=PDJ2(IR2)-XK(I)*Y(IR2)*XNT
 endif

 if (IR2.eq.J) then 
 PDJ2(IPROD1)=PDJ2(IPROD1)+XK(I)*Y(IR1)*XNT
 PDJ2(IPROD2)=PDJ2(IPROD2)+XK(I)*Y(IR1)*XNT
 PDJ2(IPROD3)=PDJ2(IPROD3)+XK(I)*Y(IR1)*XNT
 PDJ2(IPROD4)=PDJ2(IPROD4)+XK(I)*Y(IR1)*XNT
 PDJ2(IR1)=PDJ2(IR1)-XK(I)*Y(IR1)*XNT
 PDJ2(IR2)=PDJ2(IR2)-XK(I)*Y(IR1)*XNT
 endif

endif

if (IR2.eq.NSP1) then

 if (IR1.eq.J) then 

 PDJ2(IPROD1)=PDJ2(IPROD1)+XK(I)
 PDJ2(IPROD2)=PDJ2(IPROD2)+XK(I)
 PDJ2(IPROD3)=PDJ2(IPROD3)+XK(I)
 PDJ2(IPROD4)=PDJ2(IPROD4)+XK(I)
 PDJ2(IR1)=PDJ2(IR1)-XK(I)
 endif

endif

ENDDO

PDJ(1:NSMAX)=PDJ2(1:NSMAX)

IF (TESTJAC.EQ.1) THEN
NUMBERJAC=0
do i=1,NSMAX
if (PDJ(i).ne.0.d0) NUMBERJAC=NUMBERJAC+1
enddo
PRINT *, 'Number of non-zero values in JAC: ', NUMBERJAC
RETURN
ENDIF

return
end

subroutine computeIAJA(Y)
use header
implicit none
integer :: i,j,k
real(kind=8), dimension(nsmax) :: Y, PDJ

call ratcon(Y)
call ratcon2(Y)

k=1

do j=1,NSMAX
call JACVW(Y,j,PDJ)

IA(j)=k

do i=1,NSMAX
if (abs(PDJ(i)).gt.1.d-99) then
JA(k)=i
k=k+1
endif
enddo

enddo

IA(NSMAX+1)=k

return
end
