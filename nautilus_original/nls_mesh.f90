subroutine mesh
! Gives the spatial mesh of 1D integrations
use header
implicit none

do ipts=1,nptmax
zaspace(ipts) = 1.d0 - 2.*real(ipts-1)/(2*nptmax-1)
enddo

zspace(:)=zaspace(:)*Hsize

if (nptmax.ne.1) then
!zstepsize = abs(zspace(2)-zspace(1))
zstepsize=2./(2*nptmax-1)*Hsize
else
zstepsize = 0.
endif

return
end
