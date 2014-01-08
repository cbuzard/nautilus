nls_control.d
by FH/VW
gas-grain code
nautilus

------------------------------------------------------------------------
SOLVER PARAMETERS
---------+--------------+-----------------------------------------------
RTOL     |    1.000D-04 | Relative tolerance
------------------------------------------------------------------------
SWITCHES
---------+--------------+-----------------------------------------------
IDUST    |            1 | accretion, grain surface reactions
ISABS    |            1 | H2 AND CO SELF-SHILEDING
IGRQM    |            0 | 0=thermal; For H,H2: 1=QM1; 2=QM2; 3=choose fastest
IMODH    |            1 | 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only
ICONS    |            0 | 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc 
IREAD    |            0 | 0 read initial abundances from the nls_control
---------+--------------+-----------------------------------------------
GAS PHASE PARAMETERS
---------+--------------+-----------------------------------------------
XNT0     |    2.000D+04 | initial gas density  3.000D+03  1.000D+07
TEMP0    |    1.000D+01 | initial gas temp
TAU0     |    1.500D+01 | initial visual extinction  1.000D+00  2.231D+02
ZETA0    |    1.300D-17 | cosmic ray ionisation rate (1.3e-17 standard value)
ZETAX    |    0.000D+00 | Ionisation rate due to X-rays (s-1)
UVGAS    |    1.000D+00 | scale fac for UV radiation field
---------+--------------+-----------------------------------------------
GRAIN PARAMETERS
---------+--------------+-----------------------------------------------
DTEMP0   |    1.000D+01 | initial dust temp
DTOGM    |    1.000D-02 | dust-to-gas ratio by mass
STICK0   |    1.000D-00 | sticking coeff for neutral species
STICKP   |    0.000D+00 | sticking coeff for positive species
STICKN   |    0.000D+00 | sticking coeff for negative species
RHOD     |    3.000D+00 | mass density of grain material
RD       |    1.000D-05 | grain radius (cm)
ACM      |    1.000D-08 | site spacing (cm)
SNS      |    1.500D+15 | site density (cm-2)
EBFAC    |    0.500D-00 | ratio Eb(I):Ed(I) (excludes H,H2); -ve means use given values
ACT      |    1.000D-08 | grain rxn activation energy constant
TSMAX    |    7.000D+01 | peak grain temp (CR heating)
CRT      |    1.000D-05 | duration (s) of peak grain T
CRFE     |    3.000D-14 | Fe-ion--grain encounter s-1 grain-1 (for 0.1 micron grain)
LAYERS   |    1.000D+02 | number of monolayers for ITYPE 17-20 currently not used in the code
ARRK     |    1.000D-02 | a-coefficient for RRK-style formation-desorption
---------+--------------+-----------------------------------------------
OUTPUT TIMES
---------+--------------+-----------------------------------------------
OTPD     |            2 | outputs per decade  2  8  64  128 (Only without diffusion)
TSTART   |    1.000D+00 | first output time, after zero (yrs)
TFINAL   |    1.000D+01 | last output time (yrs)  9.380D+05    5.000D+04     2.000D+05    1.000D+06
WSTEP    |            1 | Outputs every WSTEP timesteps (/=1 only for 1D outputs)
WSTEPR   |           10 | Outputs every WSTEPR timesteps for the rate coefficients
IRATEOUT |            1 | Spatial point for the rate output
---------+--------------+-----------------------------------------------
INITIAL ABUNDANCES, ETC
---------+--------------+-----------------------------------------------
XNMIN    |    1.000D-40 | default minimum initial frac abun
NS0      |           14 | number of initial abundances given below
---------+--------------+-----------------------------------------------
He          | 9.000000D-02 | 1.400000D-01  9.000000D-02
N           | 6.200000D-05 | 2.140000D-05  7.500000D-05
O           | 1.400000D-04 | 1.760000D-04  1.400000D-04
H           | 1.000000D-03 | 5.000000D-05
H2          | 0.499500D-00 | 0.499975D-00
C+          | 1.700000D-04 | 7.300000D-05  1.400000D-04
S+          | 8.000000D-08 | 8.000000D-08  1.500000D-05
Si+         | 8.000000D-09 | 8.000000D-09  1.950000D-05
Fe+         | 3.000000D-09 | 3.000000D-09  7.400000D-06
Na+         | 2.000000D-09 |
Mg+         | 7.000000D-09 | 7.000000D-09  2.550000D-05
P+          | 2.000000D-10 | 3.000000D-09  2.300000D-07
Cl+         | 1.000000D-09 | 4.000000D-09  1.400000D-07
F           | 6.680000D-09 |

------------------------------------------------------------------------
1D PARAMETERS
---------+--------------+-----------------------------------------------
IDIFF    |            0 | Diffusivity flag (1. constant 2. alpha)
DIFFTY   |    1.000D-02 | Turbulent diffusivity (cgs if 1) (value or parameter)
HSIZE    |    4.500D+14 | Computing box's size (cm)
MCENTER  |    1.054D+33 | Central mass (g)
DISTR    |    1.500D+15 | Radial distance (cm) for the disk model 1.500D+15 cm = 100 AU
TESTJAC  |            0 | Testing the number of non-zero jacobian elements ?
NJAC     |          314 | Number of non-zero jacobian elements (per line, result of TESTJAC=1)
-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+------------
