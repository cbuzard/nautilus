! Contains most of the common variable declarations and a few parameters
module global_variables
use numerical_types

implicit none
integer, parameter :: nptmax = 1 ! number of space points (TODO : allocatable)
integer, parameter :: NKMAX=8335
integer, parameter :: NSMAX=684
integer, parameter :: NSGAS=485
integer, parameter :: NSGRAIN=199
integer, parameter :: NEMAX=13
integer, parameter :: NS2=216 !< number of species on the grain surface
integer, parameter :: NK2=1868 !< number of reactions on the grain surface
integer, parameter :: NS1=468 !< number of species in gas phase
integer, parameter :: NK1=6467 !< number of reactions in gas phase
integer, parameter :: NITYPE=100, NOPMAX=1000
integer, parameter :: NL1=105,NL2=52,NL3=43
real(double_precision), parameter :: RXNMIN=1.0D-99
real(double_precision), parameter :: BOLTZ=1.38054D-16,GRAV=6.668D-8,PI=3.1415926535898d0
real(double_precision), parameter :: PLANCK=6.62565D-27,HBAR=1.05459D-27,AMU=1.66043D-24
real(double_precision), parameter :: ECHARGE=1.60219D-19,AVOGADRO=6.0221415D+23
real(double_precision), parameter :: TYEAR=3.1536D+07
real(double_precision), parameter :: AU=1.5D+13

real(double_precision), save :: meanw, Omega2

real(double_precision) :: RTOL

CHARACTER (len=11), dimension(nsmax), save :: SPEC

character (len=11), dimension(7,nkmax), save :: SYMBOL
CHARACTER (len=11), save :: YJH,YJH2,YH,YH2,YHE,YHEP,YE,YGRAIN,YCO
integer, save :: INDCO, INDH2, INDHE
CHARACTER (len=11), dimension(nsmax), save :: XS0

real(double_precision), dimension(nkmax), save :: XJ,A,B,C,R
real(double_precision), dimension(nkmax) :: XK
real(double_precision), dimension(nsmax), save :: XN,XNI,XN0,CTS
real(double_precision), dimension(nsmax), save :: DXDT,DXDTP,DXDTN,AWT
real(double_precision), dimension(nemax), save :: ELEMS,ESUM,EMERR
real(double_precision), dimension(nopmax), save :: DENS,AV,TIMS
real(double_precision), dimension(nsmax, nopmax), save :: XNOP,CTSOP
real(double_precision), dimension(nsmax, nopmax), save :: DXDTOP
real(double_precision), dimension(6, nopmax), save :: EQP
real(double_precision), dimension(nkmax), save :: RDIF1,RDIF2,EX1,EX2,EA, Tmin,Tmax
real(double_precision), dimension(nkmax), save :: ACT1
real(double_precision), dimension(nsmax), save :: TINDIF,TINACC,TINEVA
real(double_precision), dimension(nsmax), save :: ED,EB,DEB,DHF
real(double_precision), dimension(nsmax), save :: CHF,CONDSP,RQ1,RQ2
real(double_precision), save :: CHARGE,TATOM,CSUM,TSUM,CERR,TERR,DXDTS,DTOGM,GTODN, RAVNH
real(double_precision), save :: RD,RHOD,STICK0,STICKP,STICKN,XNMIN
real(double_precision), save :: XNT,XNT0,TEMP
real(double_precision), save :: TEMP0,DTEMP,DTEMP0,TAU,TAU0,ZETA,ZETA0,XNTI,ZETAX
real(double_precision), save :: UVGAS,UVGRA,LAYERS,ACM,SNS,TNS,ACT,TSMAX,CRT,CRFE,EBFAC
real(double_precision), save :: TSTART,TFINAL,TIME,ALPHA,BFAC,NF,A1,B1,C1,A2,B2,C2,ARRK
real(double_precision), dimension(nl1), save :: N1H2,T1H2
real(double_precision), dimension(nl2), save :: N2CO,T2CO
real(double_precision), dimension(nl3), save :: N2H2
real(double_precision), dimension(nl3), save :: T2H2,AV2,T2AV

integer, dimension(4,nopmax), save :: IEQP
integer, dimension(nopmax), save :: IORDTM
integer, dimension(nkmax), save :: INUM, itype, jsp1, jsp2, FORMULA,NUM
integer, dimension(nsmax,nkmax), save :: ICRTBL,ICROCC
integer, dimension(nsmax), save :: ICRNUM,ORDSP, icg
integer, dimension(nemax,nsmax), save :: IELM
integer, dimension(nemax), save :: ISPELM
integer, save :: ISPE
integer, dimension(0:nitype-1), save :: IRXSTA,IRXFIN
integer, save :: IT,NT,ITFLAG,OTPD,DP,NDP,NS0,IODR,IREFSP, ISORD

integer, save :: IDENS,ITEMP,IDUST,IGRQM,ICONS,IMODH,IREAD
integer, save :: IPOUT,IPMON,IPLOT,IPDET,IPRXN,IPORD,ISABS


! Diffusion and 1D variables
real(double_precision), dimension(nptmax), save :: zspace, zaspace ! space variables
real(double_precision), dimension(nsmax,nptmax), save ::  ZXN
real(double_precision), save :: zdt ! diffusive timestep
real(double_precision), save :: zstepsize ! Spatial resolution
real(double_precision), save :: Hsize ! Size of the computing box
real(double_precision), save :: diffty ! Turbulent diffusivity
real(double_precision), save :: Mcenter ! Central mass
real(double_precision), save :: Distr ! Radial distance
real(double_precision), save :: Densmax ! Maximum density of the profile
real(double_precision), save :: TAUBC ! Av at the edge of the computing box
integer, save :: idiff ! Diffusivity flag
real(double_precision), dimension(nptmax), save :: TEMP1D, DTEMP1D, DENS1D, TAU1D, ZETAX1D ! 1D physical structure
real(double_precision), dimension(nptmax), save :: DIFF1D ! 1D diffusivity profile
real(double_precision), dimension(nptmax), save :: ZNCO,ZNH2 ! 1D column density (for the self shielding)
real(double_precision), save :: NCO,NH2 ! column density (for the self shielding)
integer, save :: istep, wstep, wstepr, irateout
integer, save :: testjac, njac

integer, save :: iptstore
integer :: ipts

! For FCHEMVW

character(len=11), dimension(nsmax+1) :: SPEC2
integer, dimension(NKMAX,7) :: REACT

! For IA and JA

integer, dimension(NSMAX+1), save :: IA
integer, allocatable, save :: JA(:)

end module global_variables
