module model_1D

implicit none

contains

subroutine mesh()
! Gives the spatial mesh of 1D integrations
use global_variables
implicit none

do ipts=1,nptmax
  zaspace(ipts) = 1.d0 - 2.d0 * dble(ipts-1) / (2.d0 * nptmax - 1)
enddo

zspace(:)=zaspace(:) * BOX_SIZE

if (nptmax.ne.1) then
  !zstepsize = abs(zspace(2)-zspace(1))
  zstepsize=2./(2*nptmax-1) * BOX_SIZE
else
  zstepsize = 0.
endif

return
end subroutine mesh

! Warning !!!
! This parametric disk model overwrites some of the parameters defined in gg_control_1D and gg_control
! BOX_SIZE, DEnb_species, UVGAS, UVGRA, TAUBC
subroutine phys_1D()
use global_variables
use diffusion

implicit none

real(double_precision) :: KFACTOR  
real(double_precision) :: HTAU, TAUEST, NHEST, Hcold
real(double_precision), dimension(nptmax) :: ld1d
real(double_precision) :: TCOLD , TWARM, T100, UV100
real(double_precision) :: ZETAZERO, Nsec, LXR, RXR, COLDENS 

TEMP1D(:)=initial_gas_temperature
DTEMP1D(:)=initial_dust_temperature
TAU1D(:)=INITIAL_VISUAL_EXTINCTION
DENS1D(:)=initial_gas_density/2.
X_IONISATION_RATE1D(:)=X_IONISATION_RATE

if (nptmax.ne.1) then
  ! Here: provide TEMP1D, DTEMP1D, DENS1D, TAU1D for your physical model

  ! Custom model for the structure of the object
  ! Here the DM Tau disk model
  ! Disabled if nptmax = 1

  ! Omega**2 for the vertical GRAVITATIONAL_CONSTANTity

  Omega2 = GRAVITATIONAL_CONSTANT * CENTRAL_MASS / RADIAL_DISTANCE**3

  ! Mean molecular weight
  ! NB: if NH > NH2 then this is wrong  

  mean_molecular_weight = 2.4d0 ! cgs

  ! User provided temperature in z

  temp1d=0.d0

  ! DM Tau parametric model

  TCOLD = 10.D0
  T100 = 30.D0 !Temperature at 100 AU
  TWARM = T100*(RADIAL_DISTANCE/(100.*AU))**(-0.5)

  Hcold = sqrt(K_B*TCOLD/(mean_molecular_weight*amu*Omega2))

  ! Change the box size
  BOX_SIZE = 4.*HCOLD
  call mesh()

  ! Estimated DEnb_species using the inner gaussian (considered to be dominant in the sigma) and 
  ! and the observed sigma (0.8 g cm-2 at 100 AU with a 1.5 power law)

  write(*,*) 'Cold height scale (AU) = ',Hcold/AU
  write(*,*) 'Computing half box size (AU) = ',BOX_SIZE/AU
  write(*,*) 'Estimated DEnb_species = ',0.8*(RADIAL_DISTANCE/(100.*AU))**(-1.5)/(mean_molecular_weight*amu)/Hcold/sqrt(2.*pi)

  ! Overwrite DEnb_species with this estimate
  DEnb_species = 0.8*(RADIAL_DISTANCE/(100.*AU))**(-1.5)/(mean_molecular_weight*amu)/Hcold/sqrt(2.*pi)

  do ipts=1, nptmax
    !TEMP1D(ipts) = 8. + (20. - 8.) * 0.5*(1.+tanh((abs(zspace(ipts))-2.*Hcold)/(Hcold/3.d0)))
    TEMP1D(ipts) = TCOLD + (TWARM-TCOLD) * 0.5d0*(1.d0+tanh((abs(zspace(ipts))-2.d0*Hcold)/(Hcold/3.d0)))
    !write(*,*) zspace(ipts)/AU,Temp1d(ipts)
  enddo

  ! Temperature for dust is the same as gas (FH)

  DTEMP1D(:) = TEMP1D(:)

  ! ---------------------------------------------------------------
  ! Computations using 1st order Euler (to change if necessary) 
  ! ---------------------------------------------------------------

  ! Computation of the density through hydrostatic equilibrium

  ld1d(1)=0.d0
  do ipts=2,nptmax
    ld1d(ipts) = ld1d(ipts-1) - (log(TEMP1D(ipts))-log(TEMP1D(ipts-1))) - Omega2 * mean_molecular_weight &
    * amu / (K_B * TEMP1D(ipts)) * zspace(ipts)*(zspace(ipts)-zspace(ipts-1))
  enddo

  DENS1D(1:nptmax) = exp(ld1d(:))

  ! Rescaling using the fact that whatever A constant, rho*A is still solution
  ! Rescaling so that maxval(rho) = DEnb_species

  DENS1D(1:nptmax) = DENS1D(1:nptmax)/maxval(DENS1D(1:nptmax)) * DEnb_species
  !DENS1D(1:nptmax/2) = DEnb_species

  ! Computation of the opacity for constant absorption
  ! TAU = NH * 5.34E-22
  ! TAU = 2 * NH2 * 5.34E-22

  KFACTOR = AV_NH_ratio

  ! Computation of an estimated TAUBC using and erf function
  ! if T is constant outside, rho has a gaussian tail
  ! BOX_SIZE is usually larger than the true scale height
  ! We use an approximated expression for the erf function
  ! cf "Handbook of mathematical functions" inequality 7.1.13

  HTAU = sqrt(K_B*temp1D(1)/(mean_molecular_weight*amu*Omega2))
  NHEST = 2. * DENS1D(1) * HTAU * sqrt(2.d0)*exp(-(BOX_SIZE/HTAU/sqrt(2.d0))**2) &
  /((BOX_SIZE/HTAU/sqrt(2.d0))+sqrt((BOX_SIZE/HTAU/sqrt(2.d0))**2+2.d0))
  TAUEST = NHEST * KFACTOR

  ! Overwrite TAUBC
  TAUBC=TAUEST

  write(*,*) 'Estimated TAUBC = ',TAUEST
  write(*,*) 'Used TAUBC = ', TAUBC

  TAU1D(1)=TAUBC
  do ipts=2,nptmax
    TAU1D(ipts) = TAU1D(ipts-1) + DENS1D(ipts) * 2. * KFACTOR * (zspace(ipts-1) - zspace(ipts))
  enddo

  ! UV at 4 Hcold
  ! UV100 is divided by a factor of two to take into account
  ! the fraction of the flux entering the disk that is ejected

  UV100=410.D0/2.D0

  UVGAS = UV100/((RADIAL_DISTANCE/(100.*AU))**2+(4.*HCOLD/(100.*AU))**2)
  !UVGAS = 1.d10/((RADIAL_DISTANCE/(100.*AU))**2+(4.*HCOLD/(100.*AU))**2)

  ! X ray approximation from Maloney et al. (1996) and 

  Nsec= 26.D-3*1.6D19 ! number of second ionization of hydrogen per unit of energy eV-1 / converted in Joules
  RXR = 10.D0*6.96D10   ! radius of the X-ray source in cm
  LXR = 1.D31*1.D-7  ! X luminosity in erg.s-1 / converted in Joules
  !kTXR = 5.D3*1.6D-19  ! X energy in eV  / converted in Joules

  ZETAZERO = Nsec/2.D0 * RXR/RADIAL_DISTANCE * LXR /(4.D0*PI*RADIAL_DISTANCE*RADIAL_DISTANCE) 
  write(*,*) ZETAZERO


  do ipts=1,nptmax
    ! column density of H above the point ipts
    COLDENS=TAU1D(ipts)/KFACTOR
    X_IONISATION_RATE1D(ipts) = ZETAZERO * &
    (1.d0/(3.d0*COLDENS) * (EXP(-COLDENS*6.d-23*0.5D0**(-3))-EXP(-COLDENS*6.d-23*0.014D0**(-3))) + &
    3.d0/(8.d0*COLDENS) * (EXP(-COLDENS*2.6d-22*7.D0**(-8/3))-EXP(-COLDENS*2.6d-22*0.5D0**(-8/3))) + &
    3.d0/(8.d0*COLDENS) * (1.D0-EXP(-COLDENS*4.4d-22*7.D0**(-8/3))))
  enddo


  ! Write physical outputs in a separate file
  open(unit=666,file='dmtau_model.dat',form='formatted')
  write(666,*) 'HCOLD = ',HCOLD/AU
  write(666,*) 'DEnb_species = ',DEnb_species
  write(666,*) 'TAUBC = ',TAUBC
  write(666,*) 'UVGAS = ',UVGAS
  close(666)

  ! End of the if test for the 1D structure
endif

call diffusion_setup()

return
end subroutine phys_1D

end module model_1D