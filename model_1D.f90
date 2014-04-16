module model_1D

implicit none

contains

subroutine mesh()
! Gives the spatial mesh of 1D integrations
use global_variables
implicit none

integer :: spatial_index ! index for spatial loops

do spatial_index=1,nptmax
  zaspace(spatial_index) = 1.d0 - 2.d0 * dble(spatial_index-1) / (2.d0 * nptmax - 1)
enddo

zspace(:)=zaspace(:) * BOX_SIZE

zstepsize = 0.d0


return
end subroutine mesh

! Warning !!!
! This parametric disk model overwrites some of the parameters defined in gg_control_1D and gg_control
! BOX_SIZE, DEnb_species, UV_FLUX, UVGRA, TAUBC
subroutine phys_1D()
use global_variables
use diffusion

implicit none

integer :: spatial_index ! index for spatial loops


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

call diffusion_setup()

return
end subroutine phys_1D

end module model_1D