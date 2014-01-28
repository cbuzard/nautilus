module diffusion

use numerical_types

implicit none

contains

subroutine diffusion_setup()
use global_variables
! Computes the value of the diffusion coefficient (and its variations in space)

implicit none

if (nptmax.eq.1) then

  diff1D(:)= TURBULENT_DIFFUSIVITY

  ! No diffusion in 0D
  if (IS_DIFFUSIVITY.ne.0) then
    write(*,*) 'No diffusion allowed in 0D. Please change the IS_DIFFUSIVITY value in parameters.in'
    stop
  endif

else

  ! No diffusion
  if (IS_DIFFUSIVITY.eq.0) diff1D(:) = TURBULENT_DIFFUSIVITY ! diffty then only defines the timescale

  ! Constant diffusivity
  if (IS_DIFFUSIVITY.eq.1) diff1D(:) = TURBULENT_DIFFUSIVITY

  ! Alpha diffusivity (for disks)
  if (IS_DIFFUSIVITY.eq.2) then
    Omega2 = GRAVITATIONAL_CONSTANT * CENTRAL_MASS / RADIAL_DISTANCE**3
    do ipts=1,nptmax
      diff1D(ipts)=TURBULENT_DIFFUSIVITY*K_B*temp1D(ipts)/(mean_molecular_weight*amu)/sqrt(Omega2)
    enddo
  endif

endif

return
end subroutine diffusion_setup

subroutine zdiffusion()
use global_variables
implicit none

! Locals
real(double_precision), dimension(nptmax) :: temp_abundances
integer :: j

if (IS_DIFFUSIVITY.eq.0) then
  call nothinghappens ! the name is explicit enough I guess
  return
endif

!zstepsize = BOX_SIZE/(nptmax-1)

do j = 1, nb_species 
  temp_abundances(:)=zxn(j,:)
  call crank(temp_abundances, nptmax-1, DIFFUSIVE_TIMESTEP, zstepsize, diff1D, DENS1D, 0)
  zxn(j,:) = temp_abundances(:)
enddo

return
end subroutine zdiffusion

subroutine ztimestep()
use global_variables

implicit none

if (IS_DIFFUSIVITY.eq.0) then
  ! Log spacing of time outputs when there is no diffusion
  if (TIME.gt.1.d-2) then
    DIFFUSIVE_TIMESTEP = TIME*(10.**(1.d0/OUTPUT_PER_DECADE)-1.)
  else
    DIFFUSIVE_TIMESTEP = 1.d0 * TYEAR
  endif
else
  ! Diffusion controlled timestep
  ! Required for the Operator splitting procedure
  DIFFUSIVE_TIMESTEP = (BOX_SIZE/nptmax)**2 / maxval(diff1D) ! smallest timescale you can think of
  DIFFUSIVE_TIMESTEP = 0.25d0 * DIFFUSIVE_TIMESTEP ! To ensure a good numerical precision
endif

return
end subroutine ztimestep

subroutine nothinghappens
! dummy subroutine to disable diffusion in a readable way
implicit none

return 
end subroutine nothinghappens

subroutine crank(f,ny,dt,dy,nu,rho, ibc)
! A Crank-Nicholson scheme
! ibc is a flag for boundary conditions
! ibc = 0 -> no flux boundaries (bc is not used then)
! ibc = 1 -> user supplied boundary conditions (bc)
! ibc > 1 -> stops the code, insulting the user
implicit none

! Inputs
integer, intent(in) :: ny
real(double_precision), intent(in), dimension(0:ny) :: rho
real(double_precision), intent(in), dimension(0:ny) :: nu
real(double_precision), intent(in) :: dt
real(double_precision), intent(in) :: dy
integer, intent(in) :: ibc

! Outputs
real(double_precision), intent(out), dimension(0:ny) :: f

! Locals
integer :: ind
real(double_precision), dimension(0:ny) :: s,Q,W,x,y,z,u,v, dd1d
real(double_precision) :: d !, nu


dd1d(:)=nu(:)*rho(:)

d=dt/(dy**2)
Q(:)=rho(:)

s(:)=0.d0

do ind = 1, ny-1
  W(ind) = s(ind)*dt + d/4*(dd1d(ind+1)+dd1d(ind))*f(ind+1)/rho(ind) + (Q(ind)-d/4*(dd1d(ind+1)+2*dd1d(ind) &
  +dd1d(ind-1)))*f(ind)/rho(ind)+d/4*(dd1d(ind)+dd1d(ind-1))*f(ind-1)/rho(ind)
enddo

do ind = 1,ny-1
  x(ind) = -d/4*(dd1d(ind+1)+dd1d(ind))/rho(ind)
  y(ind) = Q(ind)/rho(ind) + d/4*(dd1d(ind+1)+2*dd1d(ind)+dd1d(ind-1))/rho(ind)
  z(ind) = -d/4*(dd1d(ind)+dd1d(ind-1))/rho(ind)
enddo

! Test
u(ny)=1.d0
v(ny)=0.d0
x(ny) = -d/2*dd1d(ny)/rho(ny)
y(ny) = Q(ny)/rho(ny) + d/4*(3*dd1d(ny)+dd1d(ny-1))/rho(ny)
z(ny) = -d/4*(dd1d(ny)+dd1d(ny-1))/rho(ny)
W(ny) = d/2*dd1d(ny)*f(ny)/rho(ny) + (Q(ny)-d/4*(3*dd1d(ny) &
+dd1d(ny-1)))*f(ny)/rho(ny)+d/4*(dd1d(ny)+dd1d(ny-1))*f(ny-1)/rho(ny)

do ind = ny, 1, -1
  u(ind-1) = -z(ind)/ (x(ind)*u(ind) + y(ind))
  v(ind-1) = ( W(ind) - x(ind)*v(ind) )/( x(ind)*u(ind) + y(ind) )
enddo

if (ibc.eq.0) f(0) =  v(0)/( 1. - u(0) )

do ind = 0, ny-1
  f(ind+1) = u(ind)*f(ind) + v(ind)
enddo

return
end subroutine crank


subroutine euler(f,ny,dt,dy,nu,rho)
implicit none

! Inputs
integer, intent(in) :: ny
real(double_precision), intent(in), dimension(0:ny) :: nu
real(double_precision), intent(in), dimension(0:ny) :: rho
real(double_precision), intent(in) :: dt
real(double_precision), intent(in) :: dy

! Outputs
real(double_precision), intent(out), dimension(0:ny) :: f

! Locals
integer :: ind
real(double_precision), dimension(0:ny) :: Q,W, dd1d
real(double_precision) :: d !, nu

dd1d(:)=nu(:)*rho(:)

d=dt/(dy**2)
Q(:)=rho(:)

do ind = 1, ny-1
  W(ind) = d/2*(dd1d(ind+1)+dd1d(ind))*f(ind+1)/rho(ind) + (Q(ind)-d/2*(dd1d(ind+1)+2*dd1d(ind) &
  +dd1d(ind-1)))*f(ind)/rho(ind)+d/2*(dd1d(ind)+dd1d(ind-1))*f(ind-1)/rho(ind)
enddo

f(1:ny-1)=W(1:ny-1)

return
end subroutine euler

end module diffusion