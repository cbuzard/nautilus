!******************************************************************************
! MODULE: dust_temperature
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines needed to compute consistently
!! the dust temperature. Only dust_temp is public.
!!\n\n The grain temperature is only a function of UV flux and visual extinction
!!. If those parameters are constant, you only need to calculate the grain 
!! temperature once, during initialisation.
!
!******************************************************************************

module dust_temperature

use iso_fortran_env, only : error_unit
use global_variables

implicit none

integer, parameter :: nwlen=5000
real(double_precision), dimension(nwlen) :: wlen
real(double_precision), dimension(nwlen) :: Iinc_dust !< incident dust flux (rayonnement ; TODO)
real(double_precision), dimension(nwlen) :: Qabs
real(double_precision), dimension(nwlen) :: Iinc
real(double_precision) :: rate_abs

private

public :: get_grain_temperature_computed !< Only dust_temp will be available outside the module.

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE get_grain_temperature_computed(time, gas_temperature, grain_temperature)

IMPLICIT NONE

! Inputs
real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]

! Outputs
real(double_precision), intent(out) :: grain_temperature !<[out] Steady state dust temperature [K]

REAL(double_precision)                  :: Tleft, Tright ! Temperature boundary
REAL(double_precision), PARAMETER       :: eps = 1.0D-08 ! Tolerance
REAL(double_precision)                  :: fss           ! Rq: fss should be equal to zero
INTEGER                       :: i
INTEGER                       :: etat          ! flag


!----
! Initialise:
!    - the wavelength in Angstrom
!    - Qabs from Draine calculations
!    - the dust emiision intensity from Draine calculations
CALL initial_tdust(wlen,Qabs,Iinc_dust)

!----
! Compute the incident radiation filed in term of specific
! intensity (erg cm-2 s-1 Angstrom-1 sr-1) after screening 
! by dust
CALL compute_Iinc(wlen,Iinc)

!----
! Compute the absorption rate by dust
CALL compute_abs(wlen,Iinc,Qabs,rate_abs)

!----
! Find the steady state dust temperature using a 
! dichotomy algorithm.

! Temperature boundary
Tleft = 2.73D+00
Tright = 1.0D+03

CALL dicho(Tleft,Tright,eps,grain_temperature,fss,etat)

!----
! Check the dust emission
CALL initial_tdust(wlen,Qabs,Iinc_dust)
!~ OPEN(10,FILE='dust_em.dat',STATUS='unknown')
!~ DO i = 1, nwlen
!~    WRITE(10,*), wlen(i),Qabs(i)*bbody(wlen(i),grain_temperature)
!~ ENDDO
!~ CLOSE(10)

END SUBROUTINE get_grain_temperature_computed
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE initial_tdust(x,Q,Iindust)

!----
! This subroutine initialise some parameters:
!    - the wavelength in Angstrom
!    - Qabs from Draine calculations
!    - the dust emission intensity from Draine calculations
!
! Draine's calculations are interpolated to the computed wavelength




IMPLICIT NONE

INTEGER, PARAMETER                         :: nQabs=241
INTEGER, PARAMETER                         :: ndustem=1001
INTEGER                                    :: i, j
REAL(KIND=8), DIMENSION(nwlen),INTENT(OUT) :: x
REAL(KIND=8), DIMENSION(nwlen),INTENT(OUT) :: Iindust
REAL(KIND=8), DIMENSION(nQabs)             :: x1
REAL(KIND=8), DIMENSION(ndustem)           :: x2
REAL(KIND=8), DIMENSION(nQabs)             :: Qdr
REAL(KIND=8), DIMENSION(nwlen),INTENT(OUT) :: Q
REAL(KIND=8), DIMENSION(ndustem)           :: Idustdr


CALL compute_wlen(wlen)

!--------------------------------------------------------------------------
! --- Read the absorption efficiency Qabs for a 0.1 micrometer grain radius
OPEN(20,FILE='../data/si_a1m1.dat', STATUS='unknown')
READ(20,*)
READ(20,*)
DO i=1,nQabs
   READ(20,40) x1(i),Qdr(i)
40 FORMAT(2(ES9.3,1X))
ENDDO
CLOSE(20)

! --- Wavelength: micrometer to Angtroms
x1 = 1e4*x1

! --- Interpolate Qabs
Q=0.0D+00
DO i=1,nQabs-1
   DO j=1,nwlen
      IF(x1(i).ge.x(j) .AND. x1(i+1).lt. x(j)) THEN
         Q(j) = Qdr(i) + (x(j) - x1(i)) * &
            & ( Qdr(i+1) - Qdr(i))/( x1(i+1) - x(i))
      ENDIF
   ENDDO
ENDDO

!--------------------------------------------------------------------------
! --- Read the dust emission from Draine & Li 2007 models
OPEN(20,FILE='../data/U1.00_1.00_MW3.1_60.txt', STATUS='unknown')
DO i = 1, 2
   READ(20,*)
ENDDO

DO i = 1, ndustem
   READ(20,50) x2(i), Idustdr(i)
50 FORMAT(2(ES9.3,2X))
ENDDO
CLOSE(20)

! --- Interpolate Idustdr
Iindust=0.0D+00
x2 = x2*1e4
DO i=ndustem-1,1,-1
   DO j=1,nwlen
      IF(x2(i).ge.x(j) .AND. x2(i+1).lt. x(j)) THEN
         Iindust(j) = Idustdr(i) + (x(j) - x2(i)) * &
                  & ( Idustdr(i+1) - Idustdr(i))/( x2(i+1) - x2(i))
      ENDIF
   ENDDO
ENDDO

!--- From erg.s-1.H-1 to erg.cm-2.s-1.A-1.sr-1
Iindust = Iindust * 5.8e21/(x*4.0*pi)

END SUBROUTINE initial_tdust
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE compute_wlen(x)

!----
! This subroutine compute the wavelength in Angstrom



IMPLICIT NONE
INTEGER                        :: i
REAL(KIND=8), DIMENSION(nwlen) :: x

x = 911.76 *  exp(log(2.2e10/911.76)/nwlen)**(/ (i, i=0,nwlen-1) /)

END SUBROUTINE compute_wlen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE compute_Iinc(x,Iout)



IMPLICIT NONE

REAL(KIND=8), DIMENSION(nwlen), INTENT(IN)    :: x
REAL(KIND=8)                                  :: Iin
REAL(KIND=8), DIMENSION(nwlen), INTENT(OUT)   :: Iout
INTEGER                                       :: i

Iout = 0.0
!OPEN(10,FILE="extinction_curve.dat",STATUS="UNKNOWN")
DO i=1,nwlen
   Iin = max(Iinc1(x(i)),Iinc2(x(i)),Iinc3(x(i)),Iinc_dust(i))
   Iout(i) = Iin * 10**(-ext(x(i))*visual_extinction/2.5)
!~    WRITE(10,*) x(i), ext(x(i))
ENDDO
!CLOSE(10)


END SUBROUTINE compute_Iinc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE compute_abs(x,Iin,Q,dustabs)



IMPLICIT NONE
REAL(KIND=8), DIMENSION(nwlen), INTENT(IN)    :: x
REAL(KIND=8), DIMENSION(nwlen), INTENT(IN)    :: Iin
REAL(KIND=8), DIMENSION(nwlen), INTENT(IN)    :: Q
REAL(KIND=8), INTENT(OUT)                     :: dustabs
REAL(KIND=8), DIMENSION(nwlen)                :: Uin
REAL(KIND=8)                                  :: sigmabs
REAL(KIND=8)                                  :: test
REAL(KIND=8), DIMENSION(nwlen)                :: Utest
INTEGER :: i

! --- Convert the specific intensity into specific energy density
Uin = 4.0 * pi * Iin / SPEED_OF_LIGHT
!OPEN(10,FILE="ISRF.dat",STATUS="UNKNOWN")
!DO i =1,nwlen
!   WRITE(10,*), wlen(i), Uin(i)
!ENDDO
!CLOSE(10)

! --- Test: Compute the integrated energy density
Utest=0.0
DO i = 1,nwlen
   IF(x(i).ge.911.0 .AND. x(i).le. 2460.0) Utest(i) = Uin(i)
ENDDO

CALL trap_int(nwlen,x,Utest,test)
!~ WRITE(*,70), "ISRF from 912 to 2460 Angstroms = ", test,"erg.cm-3"
!~ 70   FORMAT(A35,2X,ES14.8,2X,A8)

CALL trap_int(nwlen,x,Uin,test)
!~ WRITE(*,70), "ISRF from 912 to 2E10 Angstroms = ", test,"erg.cm-3"
!---

! --- Compute the heating rate by absorption of radiation
CALL trap_int(nwlen,x,Uin*Q,sigmabs)
dustabs = sigmabs * pi * grain_radius**2 * SPEED_OF_LIGHT

END SUBROUTINE compute_abs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(KIND=8)  FUNCTION get_ss_grain(Tss)



IMPLICIT NONE

REAL(KIND=8), INTENT(IN)                   :: Tss
REAL(KIND=8)                               :: rate_emi
REAL(KIND=8), DIMENSION(nwlen)             :: dustem
REAL(KIND=8)                               :: sigmaem
INTEGER :: i

! --- Compute the cooling rate by infrared emission
DO i = 1, nwlen
   dustem(i) = bbody(wlen(i),Tss)
ENDDO

CALL trap_int(nwlen,wlen,dustem*Qabs, sigmaem)
rate_emi = 4.0D+00 * pi**2 * grain_radius**2 * sigmaem

! --- cooling balance heating
get_ss_grain = rate_emi - rate_abs

END FUNCTION get_ss_grain
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(KIND=8) FUNCTION  Iinc1(x)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------
! 1st component
!   Local UV background between 912 and 2460 Angstroms (From Mathis et al. 1983)
!   The fit parameters come from "Physics of the ISM and IGM", Draine, 2011
!      - x in Angstroms
!      - starlight(x) in erg cm-2 s-1 Angstrom-1 sr-1
!------------------------------------------------------------------------------



IMPLICIT NONE

REAL (KIND = 8), INTENT (IN)      :: x          ! Wavelength (Angstrom)
REAL (KIND = 8)                   :: aux

IF (x .GE. 1340.0 .AND. x .LT. 2460.0) THEN
   Iinc1 = 2.373D-14 * (x/1e4)**(-0.6678)
ELSEIF (x .GE. 1100.0 .AND. x .LT. 1340.0) THEN
   Iinc1 = 6.825D-13*(x/1e4)
ELSEIF (x .GE. 912.0 .AND. x .LT. 1100.0) THEN
   Iinc1 = 1.287D-9 * (x/1e4)**(4.4172)
ELSE
   Iinc1 = 0.0D+00
ENDIF

Iinc1 = Iinc1*SPEED_OF_LIGHT/4.0D+00/pi/x*UV_FLUX

END FUNCTION Iinc1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(KIND=8) FUNCTION  Iinc2(x)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------
! 2nd component
!   Starlight emmission from near-UV, visible and near-IR (From Mathis et al. 1983)
!   3 diluted Black Body (From Mathis et al. 1983)
!      - x in Angstroms
!      - starlight(x) in erg cm-2 s-1 Angstrom-1 sr-1
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL (KIND = 8), INTENT (IN)      :: x          ! Wavelength (Angstrom)

Iinc2 =  1.0e-14*bbody(x,7500D+00) &
       + 1.0e-13*bbody(x,4000D+00) &
       + 4.0e-13*bbody(x,3000D+00)

END FUNCTION Iinc2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(KIND=8) FUNCTION  Iinc3(x)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------
! 3rd component
!   CMB component at T=2.725 K
!      - x in Angstroms
!      - starlight(x) in erg cm-2 s-1 Angstrom-1 sr-1
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL (KIND = 8), INTENT (IN)      :: x          ! Wavelength (Angstrom)


Iinc3 = bbody(x,2.725D+00)

END FUNCTION Iinc3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(KIND=8) FUNCTION  bbody(x,T)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------
! Black body. x in Angstrom, T in K, bbody in erg cm-2 s-1 Angstrom-1 sr-1
!------------------------------------------------------------------------------



IMPLICIT NONE

REAL (KIND = 8), INTENT (IN)      :: x          ! Wavelength (Angstrom)
REAL (KIND = 8), INTENT (IN)      :: T          ! Temperature (K)

bbody = 2.0 * PLANCK_CONSTANT * SPEED_OF_LIGHT**2 * 1.0e32 / &
       (x**5 * (exp((PLANCK_CONSTANT*SPEED_OF_LIGHT/K_B)*1.0e8 / (x*T)) - 1.0))

END FUNCTION bbody
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(KIND=8) FUNCTION  ext(wl)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fit parameter from Cardelli et al. 1989 (1989ApJ...345..245C)
! -- The x range is optimized to have a smooth transition between each
!    components of the extinction curve

IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: wl
REAL(KIND=8) :: x
REAL(KIND=8) :: a, b, fa, fb,y

! convert x which is in Angstrom to micrometer-1
x = 1.0e4/wl

a = 0.0
b = 0.0
! --- Infrared: 0.0 to 1.4 micrometer-1
IF(x .lt. 1.4) THEN
   a =  0.574 * x**1.61
   b = -0.527 * x**1.61
! Optical/NIR: 1.4 to 2.7 micrometer-1
ELSEIF(x .ge. 1.4 .and. x .le. 2.7) THEN
  y = x - 1.82
  a = 1.00 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + &
      0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
  b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - &
      0.62251*y**5 + 5.30260*y**5 - 2.09002*y**7
! UV: 2.7 to 8.0 micrometer-1
ELSEIF(x .ge. 2.7 .and. x .le. 8.0) THEN
  IF(x .ge. 5.9 .and. x .le. 8.0) THEN
     fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
     fb =  0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3
  ELSEIF(x .le. 5.9) THEN
     fa = 0.0
     fb = 0.0
  ENDIF
  a = 1.752 - 0.316*x - 0.104/((x-4.47)**2 + 0.341) + fa
  b = -3.090 + 1.825*x + 1.206/((x-4.62)**2 + 0.263) + fb
! Far-UV: 8.0 to ...  micrometer-1
ELSEIF(x .ge. 8.0) THEN
  a = -1.073 - 0.628*(x-8.0) + 0.137*(x-8.0)**2 - 0.070*(x-8.0)**3
  b = 13.670 + 4.257*(x-8.0) - 0.420*(x-8.0)**2 + 0.374*(x-8.0)**3
ENDIF

ext = a + b/3.1

END FUNCTION ext
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE trap_int(n,x,f,res)

IMPLICIT NONE
INTEGER :: i
REAL(KIND=8) :: h, int
REAL(KIND=8), INTENT(OUT) :: res
REAL(KIND=8),DIMENSION(n), INTENT(IN) :: x, f
INTEGER, INTENT(IN) :: n

res = 0.0
int = 0.0

DO i=1,n-1
h = abs((x(i+1) - x(i))/2.0)
int = h * (f(i) + f(i+1))
res = res + int
ENDDO

END SUBROUTINE trap_int
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE dicho(xl_in,xr_in,eps_in,xtry,ftry,etat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IMPLICIT NONE

REAL(KIND=8), INTENT(IN)  :: xl_in, xr_in       ! Borne gauche et droite initiale
REAL(KIND=8), INTENT(IN)  :: eps_in             ! Precision
REAL(KIND=8)              :: xl, xr             ! Borne gauche et droite
REAL(KIND=8)              :: fl, fr             ! f(xl) et f(xr)
REAL(KIND=8), INTENT(OUT) :: xtry, ftry         ! Point milieu
INTEGER,       INTENT(OUT):: etat               ! Flag
REAL(KIND=8)              :: dxrel              ! Variation
INTEGER                   :: i                  ! Compteur
INTEGER ,      PARAMETER  :: imax = 10000
REAL(KIND=8) :: ss_grain

!--- Initialisation

xl = xl_in
xr = xr_in

fl = get_ss_grain(xl)
fr = get_ss_grain(xr)

!--- consistency test
IF (xl >= xr) THEN
   PRINT *, " Il faut xl < xr", xl, xr
   STOP
ENDIF

IF (fl*fr < 0.0D+00) THEN                       ! The zero is between xl and xr
   etat = 0                                     ! No solution found yet
   i    = 0
   DO WHILE (etat == 0)
      i = i + 1
      xtry = 0.5D+00 * (xl + xr)
      ftry = get_ss_grain(xtry)
      IF (fl*ftry > 0.0D+00) THEN
         xl = xtry
         fl = ftry
      ELSE
         xr = xtry
         fr = ftry
      ENDIF


      IF (xtry /= 0.0D+00) THEN
         dxrel = ABS((xr-xl)/xtry)
      ELSE
         dxrel = ABS(xr-xl)
      ENDIF
      IF (dxrel < eps_in .AND. ABS(ftry) < eps_in) THEN
         etat = 1
      ELSE IF (i > imax) THEN
         etat = 2
      ENDIF
   ENDDO

   !PRINT*, "Number of iter =",i

ELSE
   etat = 3
ENDIF

IF (etat == 3)THEN
   write(Error_unit,*) " "
   write(Error_unit,*) "--------------------------------------------------------  "
   write(Error_unit,*) " Error in subroutine DICHO of dust_temperature.f90:"
   write(Error_unit,*) " the sign of the function doesn't change in the interval:"
   write(Error_unit,*) " x(left)  :", xl," f(left)   :", fl
   write(Error_unit,*) " x(right) :", xr," f(right)  :", fr
   write(Error_unit,*) "--------------------------------------------------------  "
   write(Error_unit,*) " "
   call exit(11)
ENDIF

END SUBROUTINE dicho
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module dust_temperature