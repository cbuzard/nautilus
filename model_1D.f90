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

!~ zspace(:)=zaspace(:) * BOX_SIZE

zstepsize = 0.d0


return
end subroutine mesh

end module model_1D