subroutine uvw
	
! update the velocity field 

use vars
use params
use grid, only: rank
implicit none

integer :: i,j,k
	
u(1:nx,1:ny,1:nzm) = dudt(1:nx,1:ny,1:nzm,nc)
v(1:nx,1:ny,1:nzm) = dvdt(1:nx,1:ny,1:nzm,nc)
w(1:nx,1:ny,1:nzm) = dwdt(1:nx,1:ny,1:nzm,nc)

! openbl, seperately treat the northmost bc
if (doopeny) then
  if (rank.gt.nsubdomains_x*(nsubdomains_y-1)-1) then
    v(1:nx,ny+1,1:nzm) = dvdt(1:nx,ny+1,1:nzm,nc)
  end if
end if

end subroutine uvw
