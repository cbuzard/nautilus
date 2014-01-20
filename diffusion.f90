module diffusion

use numerical_types

implicit none

contains

subroutine diffusion_setup()
use global_variables
! Computes the value of the diffusion coefficient (and its variations in space)

implicit none

if (nptmax.eq.1) then

  diff1D(:)=diffty

  ! No diffusion in 0D
  if (idiff.ne.0) then
    write(*,*) 'No diffusion allowed in 0D. Please change the IDIFF value in nls_control.d'
    stop
  endif

else

  ! No diffusion
  if (idiff.eq.0) diff1D(:)=diffty ! diffty then only defines the timescale

  ! Constant diffusivity
  if (idiff.eq.1) diff1D(:)=diffty

  ! Alpha diffusivity (for disks)
  if (idiff.eq.2) then
    Omega2 = GRAVITATIONAL_CONSTANT * MCENTER / DISTR**3
    do ipts=1,nptmax
      diff1D(ipts)=diffty*K_B*temp1D(ipts)/(meanw*amu)/sqrt(Omega2)
    enddo
  endif

endif

if ((idiff.lt.0).or.(idiff.gt.2)) then
  write(*,*) 'This value for idiff is not implemented: ',IDIFF
  stop
endif

return
end subroutine diffusion_setup

subroutine zdiffusion()
use global_variables
implicit none

! Locals
real(double_precision), dimension(nptmax) :: Y
integer :: j

if (idiff.eq.0) then
  call nothinghappens ! the name is explicit enough I guess
  return
endif

!zstepsize = Hsize/(nptmax-1)

do j = 1, nb_species 
  Y(:)=zxn(j,:)
  call crank(Y, nptmax-1, zdt, zstepsize, diff1D, denb_species_for_gasd, 0)
  zxn(j,:) = Y(:)
enddo

return
end subroutine zdiffusion

subroutine ztimestep()
use global_variables

implicit none

if (idiff.eq.0) then
  ! Log spacing of time outputs when there is no diffusion
  if (TIME.gt.1.d-2) then
    zdt = TIME*(10.**(1.d0/OTPD)-1.)
  else
    zdt=1.d0*TYEAR
  endif
else
  ! Diffusion controlled timestep
  ! Required for the Operator splitting procedure
  zdt = (Hsize/nptmax)**2/maxval(diff1D) ! smallest timescale you can think of
  zdt = zdt/4 ! To ensure a good numerical precision
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

! Input
integer, intent(in) :: ny
real(double_precision), intent(in), dimension(0:ny) :: rho
real(double_precision), intent(in), dimension(0:ny) :: nu
real(double_precision), intent(in) :: dt
real(double_precision), intent(in) :: dy
integer, intent(in) :: ibc

! Output
real(double_precision), intent(out), dimension(0:ny) :: f

! Locals
integer :: ind
real(double_precision), dimension(0:ny) :: s,Q,W,x,y,z,u,v, dd1d
real(double_precision) :: d !, nu


dd1d(:)=nu(:)*rho(:)

d=dt/(dy**2)
Q(:)=rho(:)

s(:)=0.

do ind = 1, ny-1
  W(ind) = s(ind)*dt + d/4*(dd1d(ind+1)+dd1d(ind))*f(ind+1)/rho(ind) + (Q(ind)-d/4*(dd1d(ind+1)+2*dd1d(ind) &
  +dd1d(ind-1)))*f(ind)/rho(ind)+d/4*(dd1d(ind)+dd1d(ind-1))*f(ind-1)/rho(ind)
enddo

do ind = 1,ny-1
  x(ind) = -d/4*(dd1d(ind+1)+dd1d(ind))/rho(ind)
  y(ind) = Q(ind)/rho(ind) + d/4*(dd1d(ind+1)+2*dd1d(ind)+dd1d(ind-1))/rho(ind)
  z(ind) = -d/4*(dd1d(ind)+dd1d(ind-1))/rho(ind)
enddo

! Boundary conditions

!if (ibc.eq.0) then
!u(ny-1) = 1.
!v(ny-1) = 0.
!u(ny-1)=0.
!v(ny-1)=f(0)
!else
!if (ibc.eq.1) then
!f(0) = bc
!u(ny-1) = 0.
!v(ny-1) = bc
!else
!write(*,*)  '**** Error: ibc must be either 0 or 1 *****'
!stop
!endif
!endif

! Test
u(ny)=1.
v(ny)=0.
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

! Output
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