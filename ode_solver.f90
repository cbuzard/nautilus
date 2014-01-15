module ode_solver

implicit none

contains

subroutine CHEMSETUP
! Initialize the reactants and products
use global_variables
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
end subroutine CHEMSETUP

subroutine FCHEMVW(N,Y,YDOT)
! Computes the chemical evolution
use global_variables
implicit none
integer :: N, NSP1
real(double_precision), dimension(NSMAX) :: Y, YDOT
real(double_precision), dimension(NSMAX+1) :: YD2
!REAL(KIND=16), dimension(NSMAX+1) :: YD2
integer :: i
integer :: IR1, IR2, IR3, IPROD1, IPROD2, IPROD3, IPROD4
real(double_precision) :: rate

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
end subroutine FCHEMVW

subroutine JACVW(Y,J,PDJ)
! Computes columns of the chemical jacobian
use global_variables
implicit none
integer :: N, NSP1
integer J
real(double_precision), dimension(NSMAX) :: Y, PDJ
real(double_precision), dimension(NSMAX+1) :: PDJ2
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
end subroutine JACVW

subroutine computeIAJA(Y)
use global_variables
implicit none
integer :: i,j,k
real(double_precision), dimension(nsmax) :: Y, PDJ

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
end subroutine computeIAJA

  ! ======================================================================
  ! ======================================================================
  subroutine JAC(N, T, Y, J, IAN, JAN, PDJ)
  use global_variables
  implicit none
  integer N,J
  real(double_precision) :: T
  real(double_precision), dimension(N) :: IAN, JAN
  real(double_precision), dimension(NSMAX) :: Y,PDJ

  !      CALL RATCON2(Y)

  call JACVW(Y,J,PDJ) 

  return
  end subroutine JAC
  
  ! ======================================================================
  ! ======================================================================
  SUBROUTINE FCHEM (N,T,Y,YP)
  use global_variables
  use constants
  implicit none

  real(double_precision), dimension(NSMAX) :: Y,YP
  real(double_precision) :: T
  integer :: N

  CALL RATCON2(Y)
  CALL FCHEMVW(N,Y,YP) 

  RETURN
  END SUBROUTINE FCHEM

  ! ======================================================================
  ! ======================================================================
  SUBROUTINE RATCON(Y)

  ! Reactions coefficient formally dependent on the abundances Y are 
  ! computed in a companion subroutine: RATCON2
  ! Grain surface reactions, self shielding, etc...
  ! VW modification of everything for the new gas-phase network
  ! Fev 2012

  use global_variables 
  use constants
  implicit none

  real(double_precision) :: T300, TI, TSQ
  real(double_precision), dimension(NSMAX) :: Y
  integer :: nsta, nfin
  integer :: k, j, w, m, n
  integer, dimension(10) :: indice
  real(double_precision), dimension(10) :: distmin, distmax

  T300=TEMP/300.D+0
  TI=1.0D+00/TEMP
  TSQ=SQRT(TEMP)

  ! ====== Rxn ITYPE 0
  ! ITYPE 0: Gas phase reactions with GRAINS =) 
  DO J=IRXSTA(0),IRXFIN(0)
    XK(J)=A(J)*(T300**B(J))
  ENDDO

  ! In the case idust eq 0, we still need the formation of H on the grains
  ! this is done with the XH species, reaction types 10 and 11

  ! ====== Rxn ITYPE 10 and 11
  ! ITYPE 10 and 11: H2 formation on the grains when idust eq 0
  if (idust.eq.0) then
    DO J=IRXSTA(10),IRXFIN(10)
      XK(J)=A(J)*1.186D7*exp(225.D0/TEMP)**(-1)*GTODN/XNT
      !                     print *,'GTODN', GTODN
    enddo
    DO J=IRXSTA(11),IRXFIN(11)       
      XK(J)=A(J)*(T300**B(J))*XNT/GTODN
    enddo       
  endif


  ! ====== Rxn ITYPE 1
  ! ITYPE 1: Photodissoc/ionisation with cosmic rays
  ! Add X-rays in case ZETAX is not 0
  DO J=IRXSTA(1),IRXFIN(1)
    XK(J)=A(J)*(ZETA0+ZETAX)
  ENDDO
  DO J=IRXSTA(2),IRXFIN(2)
    XK(J)=A(J)*(ZETA0+ZETAX)
  ENDDO

  ! ====== Rxns ITYPE 4 - 8
  ! Bimolecular gas phase reactions - several possible formula 
  W=1
  NSTA=0
  NFIN=0
  distmin(:)=9999.d0
  distmax(:)=9999.d0
  DO J=4,8
    IF ((IRXSTA(J).NE.0).AND.(NSTA.EQ.0)) NSTA=J
    IF (IRXFIN(J).NE.0) NFIN=J
  ENDDO
  DO J=IRXSTA(NSTA),IRXFIN(NFIN)

    !---------------  KOOIJ FORMULA
    IF (FORMULA(J).eq.3) THEN

      XK(J)=A(J)*(T300**B(J))*EXP(-C(J)/TEMP)

      !              Check for temperature bounderies
      if (TEMP.LT.Tmin(J)) XK(J)=A(J)*((Tmin(J)/300.D0)**B(J))*EXP(-C(J)/Tmin(J))
      if (TEMP.GT.Tmax(J)) XK(J)=A(J)*((Tmax(J)/300.D0)**B(J))*EXP(-C(J)/Tmax(J))

      !              Check for the presence of several rate coefficients present in the network for the
      !              the same reaction
      IF (NUM(J+1).EQ.NUM(J)) THEN
        INDICE(W)=J
        distmin(w)=tmin(j)-temp
        distmax(w)=temp-tmax(j)
        W = W + 1
      ENDIF

      IF ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) THEN

        INDICE(W)=J
        distmin(w)=tmin(j)-temp
        distmax(w)=temp-tmax(j)

        DO M=1,W
          N=INDICE(M)
          !IF(IT==1) PRINT*,N,M, SYMBOL(:,N), tmin(N), tmax(N),distmin(M),distmax(M)
          IF (TEMP.LT.Tmin(N)) XK(N)=0.d0
          IF (TEMP.GT.Tmax(N)) XK(N)=0.d0
        ENDDO

        IF (maxval(XK(indice(1:w))).lt.1.d-99) THEN

          IF (minval(abs(distmin)).lt.minval(abs(distmax))) THEN
            N=indice(minloc(abs(distmin),dim=1))
            XK(N)=A(N)*((Tmin(N)/300.D0)**B(N))*EXP(-C(N)/Tmin(N))
          ELSE
            N=indice(minloc(abs(distmax),dim=1))
            XK(N)=A(N)*((Tmax(N)/300.D0)**B(N))*EXP(-C(N)/Tmax(N))
          ENDIF
        ENDIF

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      ENDIF
    ENDIF

    !---------------  IONPOL1 FORMULA
    IF (FORMULA(J).EQ.4) THEN
      XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/TEMP)**0.5))

      !              Check for temperature bounderies
      IF (TEMP.LT.Tmin(J)) XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/Tmin(J))**0.5))
      IF (TEMP.GT.Tmax(J)) XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/TMAX(J))**0.5))

      !              Check for the presence of several rate coefficients present in the network for the
      !              the same reaction
      IF (NUM(J+1).EQ.NUM(J)) THEN
        INDICE(W)=J
        distmin(w)=tmin(j)-temp
        distmax(w)=temp-tmax(j)
        W = W + 1
      ENDIF

      IF ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) THEN

        INDICE(W)=J
        distmin(w)=tmin(j)-temp
        distmax(w)=temp-tmax(j)

        DO M=1,W
          N=INDICE(M)
          IF (TEMP.LT.Tmin(N)) XK(N)= 0.d0
          IF (TEMP.GT.Tmax(N)) XK(N)= 0.d0
        ENDDO

        IF (maxval(XK(indice(1:w))).lt.1.d-99) THEN
          IF (minval(abs(distmin)).lt.minval(abs(distmax))) THEN
            N=indice(minloc(abs(distmin),dim=1))
            XK(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/Tmin(N))**0.5))
          ELSE
            N=indice(minloc(abs(distmax),dim=1))
            XK(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/TMAX(N))**0.5))
          ENDIF
        ENDIF

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      ENDIF
    ENDIF

    !---------------  IONPOL2 FORMULA
    IF (FORMULA(J).EQ.5) then
      XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TEMP)**0.5)+(C(J)**2*300.d0/(10.526*TEMP)))

      !               Check for temperature bounderies
      IF (TEMP.LT.Tmin(J)) XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMIN(J))**0.5)+(C(J)**2*300.d0/(10.526*TMIN(J))))
      IF (TEMP.GT.Tmax(J)) XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMAX(J))**0.5)+(C(J)**2*300.d0/(10.526*TMAX(J))))

      !               Check for the presence of several rate coefficients present in the network for the
      !               the same reaction
      IF (NUM(J+1).EQ.NUM(J)) THEN
        INDICE(W)=J
        distmin(w)=tmin(j)-temp
        distmax(w)=temp-tmax(j)
        W = W + 1
      ENDIF

      IF ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) THEN

        INDICE(W)=J
        distmin(w)=tmin(j)-temp
        distmax(w)=temp-tmax(j)

        DO M=1,W
          N=INDICE(M)
          IF (TEMP.LT.Tmin(N)) XK(N)=       0.d0
          IF (TEMP.GT.Tmax(N)) XK(N)=       0.d0
        ENDDO

        IF (maxval(XK(indice(1:w))).lt.1.d-99) THEN
          IF (minval(abs(distmin)).lt.minval(abs(distmax))) THEN
            N=indice(minloc(abs(distmin),dim=1))
            XK(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMIN(N))**0.5)+(C(N)**2*300.d0/(10.526*TMIN(N))))
          ELSE
            N=indice(minloc(abs(distmax),dim=1))
            XK(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMAX(N))**0.5)+(C(N)**2*300.d0/(10.526*TMAX(N))))
          ENDIF
        ENDIF

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      ENDIF
    ENDIF
  ENDDO

  ! === Grain surfaces
  IF (IDUST.NE.0) THEN

    ! ========= Set diffusion and evaporation rates (s-1)
    DO K=1,NSMAX
      TINDIF(K)=CHF(K)*EXP(-EB(K)/DTEMP)/TNS
      TINEVA(K)=CHF(K)*EXP(-ED(K)/DTEMP)
    ENDDO

    ! ========= Rxn ITYPE 15 - thermal evaporation
    ! ITYPE 15: Thermal evaporation
    DO J=IRXSTA(15),IRXFIN(15)
      XK(J)=A(J)*XJ(J)*TINEVA(JSP1(J))
    ENDDO

    ! ========= Rxn ITYPE 16
    ! ITYPE 16: Cosmic-ray evaporation
    DO J=IRXSTA(16),IRXFIN(16)
      XK(J)=A(J)*XJ(J)*((ZETA0+ZETAX)/1.3D-17)&
      *CHF(JSP1(J))*CRFE*CRT*EXP(-ED(JSP1(J))/TSMAX)
    ENDDO


    ! Photodesorption, when used appears through ITYPES 66 and 67
    if (irxsta(66).ne.0) then
      ! ========= Rxn ITYPE 66
      ! ITYPE 66: CO photodesorption by external UV
      ! 1.d8 is I_ISRF-FUV from Oberg et al. 2007, ApJ, 662, 23
      do j=irxsta(66),irxfin(66)
        XK(J)=A(J)/SNS*UVGAS*1.d8*exp(-2.*TAU) 
      enddo

      ! ========= Rxn ITYPE 67
      ! ITYPE 67: CO photodesorption by CR generated UV
      do j=irxsta(67),irxfin(67)
        XK(J)=A(J)/SNS*1.d4
      enddo

    endif

    if (irxsta(98).ne.0) then
      ! ========= Rxn ITYPE 98 test the storage of H2S under a refractory form
      ! ITYPE 98: storage of H2S under a refractory form
      do j=irxsta(98),irxfin(98)
        XK(J)=A(J)*(T300**B(J))*EXP(-C(J)*TI)
      enddo
    endif


    ! ====== Rxn ITYPE 17
    ! ITYPE 17: Photodissociations by Cosmic rays on grain surfaces
    ! Add X-rays in case ZETAX is not 0
    DO J=IRXSTA(17),IRXFIN(17)
      XK(J)=A(J)*(ZETA0 + ZETAX)
      !            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
    ENDDO

    ! ====== Rxn ITYPE 18
    ! ITYPE 18: Photodissociations by Cosmic rays on grain surfaces
    DO J=IRXSTA(18),IRXFIN(18)
      XK(J)=A(J)*(ZETA0+ZETAX)
      !            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
    ENDDO

  ENDIF

  ! When dust is turned off, zero all dust rates==========================
  IF ((IDUST.EQ.0).AND.(IT.EQ.1)) THEN
    DO J=IRXSTA(14),IRXFIN(99)
      XK(J)=0.0D+0
      XJ(J)=0.0D+0
    ENDDO
  ENDIF


  RETURN
  END SUBROUTINE RATCON
  
    ! ======================================================================
  ! ======================================================================
  SUBROUTINE RATCON2(Y)

  use global_variables 
  use constants
  implicit none

  real(double_precision) :: ACTIV,BARR,MONLAY,DIFF
  real(double_precision) :: XNH2,XNCO
  real(double_precision) :: TETABIS,TETABIS1,TETABIS2,TETABIS3
  real(double_precision) :: T300, TI, TSQ
  real(double_precision), dimension(NSMAX) :: Y
  real(double_precision) :: YMOD1, YMOD2
  integer IMOD1,IMOD2
  integer :: j, l

  T300=TEMP/300.D+0
  TI=1.0D+00/TEMP
  TSQ=SQRT(TEMP)
  MONLAY=LAYERS*TNS/GTODN

  XNH2=Y(indH2)
  XNCO=Y(indCO)

  ! Density/Av-dependent==================================================
  !
  ! VW 02/07 the treatment for the H2 and CO photodissociation was corrected
  ! for larger and smaller values of H2, CO column densities and Av not in the 
  ! data tables. I also added the choice of using this approximation in gg_control.d
  !
  ! ====== Rxn ITYPE 13
  ! ITYPE 2: Gas phase photodissociations/ionisations by UV
  DO J=IRXSTA(3),IRXFIN(3)
    XK(J)=A(J)*EXP(-C(J)*TAU)*UVGAS

    ! MODIFY THE H2 AND CO PHOTODISSOCIATION IF ISABS EQ 1
    IF (ISABS.EQ.1) THEN 

      ! ====== Compute the H2 self-shielding
      IF (SYMBOL(1,J).EQ.YH2) THEN
        TETABIS=1.D0
        if (iptstore.eq.1) then
          ZNH2(iptstore) = TAU/RAVNH * XNH2
        endif
        NH2=ZNH2(iptstore)

        ! ======= Linear extrapolation of the shielding factors
        DO L=1,NL1-1
          IF ((N1H2(L).LE.NH2).AND.(N1H2(L+1).GE.NH2)) THEN
            TETABIS=T1H2(L)+(NH2-N1H2(L))*(T1H2(L+1)-T1H2(L))/(N1H2(L+1)-N1H2(L))
          ENDIF
        ENDDO
        IF (NH2.GT.N1H2(NL1)) TETABIS = T1H2(NL1)

        XK(J)=2.54D-11*TETABIS

        XK(J)=XK(J)*UVGAS
      ENDIF

      ! ====== Compute the CO self-shielding
      IF (SYMBOL(1,J).EQ.YCO) THEN

        TETABIS1=1.D0
        TETABIS2=1.D0
        TETABIS3=1.D0

        if (iptstore.eq.1) then
          ZNH2(iptstore) = TAU/RAVNH * XNH2
          ZNCO(iptstore) = TAU/RAVNH * XNCO
        endif
        NH2=ZNH2(iptstore)
        NCO=ZNCO(iptstore)

        ! ======= Linear extrapolation of the three shileding factors
        DO L=1,NL2-1
          IF ((N2CO(L).LE.NCO).AND.(N2CO(L+1).GE.NCO))  &
          TETABIS2=T2CO(L)+(NCO-N2CO(L))*(T2CO(L+1)-T2CO(L))&
          /(N2CO(L+1)-N2CO(L))
        ENDDO

        DO L=1,NL3-1
          IF ((N2H2(L).LE.NH2).AND.(N2H2(L+1).GE.NH2)) &
          TETABIS1=T2H2(L)+(NH2-N2H2(L))*(T2H2(L+1)-T2H2(L))&
          /(N2H2(L+1)-N2H2(L))
          IF ((AV2(L).LE.TAU).AND.(AV2(L+1).GE.TAU)) &
          TETABIS3=T2AV(L)+(TAU-AV2(L))*(T2AV(L+1)-T2AV(L))&
          /(AV2(L+1)-AV2(L))
        ENDDO

        ! Saturate the rate coefficient if necessary (when density or Av are out of 
        ! the shielding array, from the photodiss files)

        IF (NCO.GT.N2CO(NL2)) TETABIS2 = T2CO(NL2)
        IF (NH2.GT.N2H2(NL3)) TETABIS1 = T2H2(NL3)
        IF (TAU.GT.AV2(NL3))  TETABIS3 = T2AV(NL3)

        XK(J)=1.03D-10*TETABIS1*TETABIS2*TETABIS3
        XK(J)=XK(J)*UVGAS

      ENDIF

    ENDIF

  ENDDO

  ! Continually time-dependent grain rates================================
  IF (IDUST.NE.0) THEN

    ! ====== Rxn ITYPE 99
    ! ITYPE 99: Adsorption on grains
    DO J=IRXSTA(99),IRXFIN(99)
      ! ========= Set accretion rates
      TINACC(JSP1(J))=CONDSP(JSP1(J))*TSQ*Y(JSP1(J))*XNT
      TINACC(JSP2(J))=TINACC(JSP1(J))
      XK(J)=A(J)*XJ(J)*TINACC(JSP1(J))/Y(JSP1(J))/GTODN
    ENDDO

    ! ====== Rxn ITYPE 14
    ! ITYPE 14: Grain surface reactions
    DO J=IRXSTA(14),IRXFIN(14)
      IMOD1=0
      IMOD2=0
      BARR=1.0D+0
      ! --------- Calculate activation energy barrier multiplier
      IF (EA(J).GE.1.0D-40) THEN
        ACTIV=EA(J)/DTEMP
        ! ------------ Choose fastest of classical or tunnelling
        IF (ACTIV.GT.ACT1(J)) ACTIV=ACT1(J)
        BARR=EXP(-ACTIV)
      ENDIF

      ! --------- Thermal hopping diffusion method
      RDIF1(J)=TINDIF(JSP1(J))
      RDIF2(J)=TINDIF(JSP2(J))

      ! --------- Check for JH,JH2
      IF (SYMBOL(1,J).EQ.YJH)  IMOD1=1
      IF (SYMBOL(1,J).EQ.YJH2) IMOD1=2
      IF (SYMBOL(2,J).EQ.YJH)  IMOD2=1
      IF (SYMBOL(2,J).EQ.YJH2) IMOD2=2

      ! --------- QM for JH,JH2 only - others are too heavy
      IF (IMOD1+IMOD2.NE.0) THEN
        ! ------------ QM1 - Tunnelling (if it's faster than thermal)
        IF (IGRQM.EQ.1) THEN
          IF ((IMOD1.NE.0).AND.&
          (RQ1(JSP1(J)).GT.RDIF1(J))) RDIF1(J)=RQ1(JSP1(J))
          IF ((IMOD2.NE.0).AND.&
          (RQ1(JSP2(J)).GT.RDIF2(J))) RDIF2(J)=RQ1(JSP2(J))
        ENDIF
        ! ------------ QM2 - Tunnelling: use estimated width of lowest energy band (if it's faster than thermal)
        IF (IGRQM.EQ.2) THEN
          IF ((IMOD1.NE.0).AND.&
          (RQ2(JSP1(J)).GT.RDIF1(J))) RDIF1(J)=RQ2(JSP1(J))
          IF ((IMOD2.NE.0).AND.&
          (RQ2(JSP2(J)).GT.RDIF2(J))) RDIF2(J)=RQ2(JSP2(J))
        ENDIF
        ! ------------ QM3 - Fastest out of thermal, QM1, QM2 rates
        IF (IGRQM.EQ.3) THEN
          IF (IMOD1.NE.0) THEN
            IF (RQ1(JSP1(J)).GT.RDIF1(J)) RDIF1(J)=RQ1(JSP1(J))
            IF (RQ2(JSP1(J)).GT.RDIF1(J)) RDIF1(J)=RQ2(JSP1(J))
          ENDIF
          IF (IMOD2.NE.0) THEN
            IF (RQ1(JSP2(J)).GT.RDIF2(J)) RDIF2(J)=RQ1(JSP2(J))
            IF (RQ2(JSP2(J)).GT.RDIF2(J)) RDIF2(J)=RQ2(JSP2(J))
          ENDIF
        ENDIF
      ENDIF

      ! --------- Modify according to IMODH switch:
      IF (IMODH.NE.0) THEN
        ! ------------ If H+H->H2 is only modified rxn:
        IF ((IMODH.EQ.-1).AND.(IMOD1.NE.1.OR.IMOD2.NE.1)) THEN
          IMOD1=0
          IMOD2=0
        ENDIF

        ! ------------ If only H is modified:
        IF ((IMODH.EQ.1).AND.(IMOD1.NE.1)) IMOD1=0
        IF ((IMODH.EQ.1).AND.(IMOD2.NE.1)) IMOD2=0

        ! ------------ Set to modify all rates, if selected (just atoms)
        IF (IMODH.EQ.3) THEN
          IF ((SYMBOL(1,J).EQ.YJH).OR.&
          (SYMBOL(1,J).EQ.'JHe        ').OR.&
          (SYMBOL(1,J).EQ.'JC         ').OR.&
          (SYMBOL(1,J).EQ.'JN         ').OR.&
          (SYMBOL(1,J).EQ.'JO         ').OR.&
          (SYMBOL(1,J).EQ.'JS         ').OR.&
          (SYMBOL(1,J).EQ.'JSi        ').OR.&
          (SYMBOL(1,J).EQ.'JFe        ').OR.&
          (SYMBOL(1,J).EQ.'JNa        ').OR.&
          (SYMBOL(1,J).EQ.'JMg        ').OR.&
          (SYMBOL(1,J).EQ.'JP         ').OR.&
          (SYMBOL(1,J).EQ.'JCl        ')) IMOD1=3
          IF ((SYMBOL(2,J).EQ.YJH).OR.&
          (SYMBOL(2,J).EQ.'JHe        ').OR.&
          (SYMBOL(2,J).EQ.'JC         ').OR.&
          (SYMBOL(2,J).EQ.'JN         ').OR.&
          (SYMBOL(2,J).EQ.'JO         ').OR.&
          (SYMBOL(2,J).EQ.'JS         ').OR.&
          (SYMBOL(2,J).EQ.'JSi        ').OR.&
          (SYMBOL(2,J).EQ.'JFe        ').OR.&
          (SYMBOL(2,J).EQ.'JNa        ').OR.&
          (SYMBOL(2,J).EQ.'JMg        ').OR.&
          (SYMBOL(2,J).EQ.'JP         ').OR.&
          (SYMBOL(2,J).EQ.'JCl        ')) IMOD2=3
        ENDIF

        ! ------------ Modify rates (RDIF1 & RDIF2) according to their own evap/acc rates
        YMOD1=Y(JSP1(J))
        YMOD2=Y(JSP2(J))
        CALL MODIF(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
      ENDIF

      DIFF=RDIF1(J)+RDIF2(J)

      XK(J)=A(J)*XJ(J)*BARR*DIFF*GTODN/XNT
      !              XK(J)=0.D0
      ! --------- Allow only 1 monolayer of each to react
      ! Not used for the time being
      !            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*TNS/GTODN/Y(JSP1(J))
      !            IF (Y(JSP2(J)).GT.MONLAY) XK(J)=XK(J)*TNS/GTODN/Y(JSP2(J))

    ENDDO

    ! ====== Rxn ITYPE 19 - 20
    ! ITYPE 19: Photodissociations by UV photons on grain surfaces
    ! ITYPE 20: Photodissociations by UV photons on grain surfaces
    DO J=IRXSTA(19),IRXFIN(20)
      XK(J)=A(J)*EXP(-C(J)*TAU)*UVGAS
      !            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
    ENDDO

    ! Useful for testing
    ! To disable some reaction types
    !do j=irxsta(14),irxfin(15)
    !xk(j)=0.
    !enddo

    !do j=irxsta(17),irxfin(20)
    !xk(j)=0.
    !enddo
  ENDIF

  ! Continually time-dependent gas phase rates============================
  ! H2 formation
  ! XJ(1) and XJ(2) are zero if IDUST=1
  ! cf GRAINRATE
  ! VW Fev 2012 - this process has been removed
  !      do j=irxsta(0),irxfin(0)
  !      if ((SYMBOL(1,J).eq.YH).and.(SYMBOL(2,j).eq.YH)) then
  !      XK(j)=XJ(j)*A(j)*(T300**B(j))*GTODN/XNT/Y(JSP1(j))
  !      endif
  !      enddo

  ! If rate acoefficients are too small, put them to 0
  where (XK.lt.RXNMIN) XK=0.d0

    RETURN
    END subroutine ratcon2
    
  ! ======================================================================
  ! ======================================================================
  SUBROUTINE MODIF(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
  use global_variables
  use constants
  implicit none

  integer :: J,IMOD1,IMOD2
  real(double_precision) :: BARR,YMOD1,YMOD2, PICK, TESTREF1, TESTREF2, TESTNUM

  EX1(J)=0.0D+0
  EX2(J)=0.0D+0

  ! --- Check value of x = t_acc/t_evap
  ! TINEVA = 1/t_evap
  ! TINACC = 1/t_acc
  IF (TINACC(JSP1(J)).GT.0.0D+0) THEN
    EX1(J)=TINEVA(JSP1(J))/TINACC(JSP1(J))
  ENDIF
  IF (TINACC(JSP2(J)).GT.0.0D+0) THEN
    EX2(J)=TINEVA(JSP2(J))/TINACC(JSP2(J))
  ENDIF
  ! Hence x = 0 if t_evap or t_acc = 0

  ! --- Assign max rates

  IF (BARR.EQ.1.0D+0) THEN
    IF (IMOD1.NE.0) THEN
      IF (EX1(J).LT.1.0D+0) THEN
        IF (RDIF1(J).GT.TINACC(JSP1(J))) RDIF1(J)=TINACC(JSP1(J))
      ENDIF
      IF (EX1(J).GE.1.0D+0) THEN
        IF (RDIF1(J).GT.TINEVA(JSP1(J))) RDIF1(J)=TINEVA(JSP1(J))
      ENDIF
    ENDIF

    IF (IMOD2.NE.0) THEN
      IF (EX2(J).LT.1.0D+0) THEN
        IF (RDIF2(J).GT.TINACC(JSP2(J))) RDIF2(J)=TINACC(JSP2(J))
      ENDIF
      IF (EX2(J).GE.1.0D+0) THEN
        IF (RDIF2(J).GT.TINEVA(JSP2(J))) RDIF2(J)=TINEVA(JSP2(J))
      ENDIF
    ENDIF
  ENDIF

  ! --- Species rate to compare chosen by fastest diffusion rate
  IF (BARR.NE.1.0D+0) THEN
    PICK=0.d0

    TESTREF1=TINACC(JSP1(J))
    IF (EX1(J).GE.1.0D+0) TESTREF1=TINEVA(JSP1(J))
    TESTREF2=TINACC(JSP2(J))
    IF (EX2(J).GE.1.0D+0) TESTREF2=TINEVA(JSP2(J))

    IF (RDIF1(J).GE.RDIF2(J)) THEN
      TESTNUM=(RDIF1(J)+RDIF2(J))*BARR*YMOD2*GTODN
      IF (YMOD2*GTODN.LT.1.0D+0) TESTNUM=(RDIF1(J)+RDIF2(J))*BARR
      IF (TESTNUM.GT.TESTREF1) PICK=1.d0
    ENDIF
    IF (RDIF2(J).GT.RDIF1(J)) THEN
      TESTNUM=(RDIF1(J)+RDIF2(J))*BARR*YMOD1*GTODN
      IF (YMOD1*GTODN.LT.1.0D+0) TESTNUM=(RDIF1(J)+RDIF2(J))*BARR
      IF (TESTNUM.GT.TESTREF2) PICK=2.d0
    ENDIF

    IF (PICK.EQ.1) THEN
      RDIF1(J)=TESTREF1/BARR/YMOD2/GTODN
      IF (YMOD2*GTODN.LT.1.0D+0) RDIF1(J)=TESTREF1/BARR
      RDIF2(J)=0.0D+0
    ENDIF

    IF (PICK.EQ.2) THEN
      RDIF2(J)=TESTREF2/BARR/YMOD1/GTODN
      IF (YMOD1*GTODN.LT.1.0D+0) RDIF2(J)=TESTREF2/BARR
      RDIF1(J)=0.0D+0
    ENDIF

  ENDIF

  RETURN
  END SUBROUTINE MODIF

end module ode_solver