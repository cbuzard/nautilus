! Contains most of the common variable declarations and a few parameters
module header
implicit none
      integer, parameter :: nptmax = 1 ! number of space points (TODO : allocatable)
      integer, PARAMETER :: NKMAX=8335,NSMAX=684,NSGAS=485,NSGRAIN=199,NEMAX=13
      integer, PARAMETER :: NS1=468,NS2=216,NK1=6467,NK2=1868
      integer, PARAMETER :: NITYPE=100, NOPMAX=1000
      double precision, parameter :: RXNMIN=1.0D-99
      integer, PARAMETER :: NL1=105,NL2=52,NL3=43
      double precision, PARAMETER :: BOLTZ=1.38054D-16,GRAV=6.668D-8,PI=3.1415926535898d0
      double precision, PARAMETER :: PLANCK=6.62565D-27,HBAR=1.05459D-27,AMU=1.66043D-24
      double precision, PARAMETER :: ECHARGE=1.60219D-19,AVOGADRO=6.0221415D+23
      double precision, PARAMETER :: TYEAR=3.1536D+07
      double precision, PARAMETER :: AU=1.5D+13
 
      real(kind=8), save :: meanw, Omega2

      real(kind=8) :: RTOL

      CHARACTER (len=11), dimension(nsmax), save :: SPEC
     
      character (len=11), dimension(7,nkmax), save :: SYMBOL
      CHARACTER (len=11), save :: YJH,YJH2,YH,YH2,YHE,YHEP,YE,YGRAIN,YCO
      integer, save :: INDCO, INDH2, INDHE
      CHARACTER (len=11), dimension(nsmax), save :: XS0

      REAL(kind=8), dimension(nkmax), save :: XJ,A,B,C,R
      REAL(kind=8), dimension(nkmax) :: XK
      REAL(kind=8), dimension(nsmax), save :: XN,XNI,XN0,CTS
      REAL(kind=8), dimension(nsmax), save :: DXDT,DXDTP,DXDTN,AWT
      REAL(kind=8), dimension(nemax), save :: ELEMS,ESUM,EMERR
      REAL(kind=8), dimension(nopmax), save :: DENS,AV,TIMS
      REAL(kind=8), dimension(nsmax, nopmax), save :: XNOP,CTSOP
      REAL(kind=8), dimension(nsmax, nopmax), save :: DXDTOP
      real(kind=8), dimension(6, nopmax), save :: EQP
      REAL(kind=8), dimension(nkmax), save :: RDIF1,RDIF2,EX1,EX2,EA, Tmin,Tmax
      REAL(kind=8), dimension(nkmax), save :: ACT1
      real(kind=8), dimension(nsmax), save :: TINDIF,TINACC,TINEVA
      REAL(kind=8), dimension(nsmax), save :: ED,EB,DEB,DHF
      REAL(kind=8), dimension(nsmax), save :: CHF,CONDSP,RQ1,RQ2
      REAL(kind=8), save :: CHARGE,TATOM,CSUM,TSUM,CERR,TERR,DXDTS,DTOGM,GTODN, RAVNH
      REAL(kind=8), save :: RD,RHOD,STICK0,STICKP,STICKN,XNMIN
      REAL(kind=8), save :: XNT,XNT0,TEMP
      REAL(kind=8), save :: TEMP0,DTEMP,DTEMP0,TAU,TAU0,ZETA,ZETA0,XNTI,ZETAX
      REAL(kind=8), save :: UVGAS,UVGRA,LAYERS,ACM,SNS,TNS,ACT,TSMAX,CRT,CRFE,EBFAC
      REAL(kind=8), save :: TSTART,TFINAL,TIME,ALPHA,BFAC,NF,A1,B1,C1,A2,B2,C2,ARRK
      REAL(kind=8), dimension(nl1), save :: N1H2,T1H2
      real(kind=8), dimension(nl2), save :: N2CO,T2CO
      real(kind=8), dimension(nl3), save :: N2H2
      REAL(kind=8), dimension(nl3), save :: T2H2,AV2,T2AV

      INTEGER, dimension(4,nopmax), save :: IEQP
      integer, dimension(nopmax), save :: IORDTM
      INTEGER, dimension(nkmax), save :: INUM, itype, jsp1, jsp2, FORMULA,NUM
      integer, dimension(nsmax,nkmax), save :: ICRTBL,ICROCC
      INTEGER, dimension(nsmax), save :: ICRNUM,ORDSP, icg
      INTEGER, dimension(nemax,nsmax), save :: IELM
      integer, dimension(nemax), save :: ISPELM
      INTEGER, save :: ISPE
      INTEGER, dimension(0:nitype-1), save :: IRXSTA,IRXFIN
      INTEGER, save :: IT,NT,ITFLAG,OTPD,DP,NDP,NS0,IODR,IREFSP, ISORD

      integer, save :: IDENS,ITEMP,IDUST,IGRQM,ICONS,IMODH,IREAD
      integer, save :: IPOUT,IPMON,IPLOT,IPDET,IPRXN,IPORD,ISABS


! Diffusion and 1D variables
real(kind=8), dimension(nptmax), save :: zspace, zaspace ! space variables
real(kind=8), dimension(nsmax,nptmax), save ::  ZXN
real(kind=8), save :: zdt ! diffusive timestep
real(kind=8), save :: zstepsize ! Spatial resolution
real(kind=8), save :: Hsize ! Size of the computing box
real(kind=8), save :: diffty ! Turbulent diffusivity
real(kind=8), save :: Mcenter ! Central mass
real(kind=8), save :: Distr ! Radial distance
real(kind=8), save :: Densmax ! Maximum density of the profile
real(kind=8), save :: TAUBC ! Av at the edge of the computing box
integer, save :: idiff ! Diffusivity flag
real(kind=8), dimension(nptmax), save :: TEMP1D, DTEMP1D, DENS1D, TAU1D, ZETAX1D ! 1D physical structure
real(kind=8), dimension(nptmax), save :: DIFF1D ! 1D diffusivity profile
real(kind=8), dimension(nptmax), save :: ZNCO,ZNH2 ! 1D column density (for the self shielding)
real(kind=8), save :: NCO,NH2 ! column density (for the self shielding)
integer, save :: istep, wstep, wstepr, irateout
integer, save :: testjac, njac

integer, save :: iptstore
integer :: ipts

! For FCHEMVW

CHARACTER(len=11), dimension(nsmax+1) :: SPEC2
integer, dimension(NKMAX,7) :: REACT

! For IA and JA

integer, dimension(NSMAX+1), save :: IA
integer, allocatable, save :: JA(:)

end module header
module unitvar
implicit none
integer, parameter ::      NSP=4
integer, parameter ::      NCON=5
integer, parameter ::      NERR=6
integer, parameter ::      NMOD=8
integer, parameter ::      NJR=9
integer, parameter ::      NJR2=19
integer, parameter ::      NGR=10
integer, parameter ::      NINF=12
integer, parameter ::      NTAI=13
integer, parameter ::      NINI=14
integer, parameter ::      NORD=15
integer, parameter ::      CODIS=16
integer, parameter ::      H2DIS=17
integer, parameter ::      NCON2=25
integer, parameter ::      NOUT2=35
integer, parameter ::      NORD2=45
end module unitvar
