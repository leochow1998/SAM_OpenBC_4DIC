subroutine advect_mom

use vars
use params, only: docolumn, doopeny ! doopeny is for openybc
use domain, only: nsubdomains_x,nsubdomains_y ! below variables are for openybc
use grid, only: rank
use open_boundary_y, only: advection_momentum

implicit none
integer i,j,k
real du(nx,ny,nz,3)

! open y-boundary condition variables
real :: dudtbackup(nx+1,1,nzm) ! for southmost MPI domain, y is y=1: southmost grid
real :: dwdtbackup(nx,1,nz)    ! for northmost MPI domain, y is y=ny: northmost grid

if(docolumn) return

call t_startf ('advect_mom')

if(dostatis) then
	
  do k=1,nzm
   do j=1,ny
    do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na) 
     du(i,j,k,2)=dvdt(i,j,k,na) 
     du(i,j,k,3)=dwdt(i,j,k,na) 
    end do
   end do
  end do

endif

! openybc: backup the dudt and dwdt before the call of advection momentum
! initially written as 2nd-order-center-differencing in flux form
! by SAM (kharioutdinov)
if (doopeny) then
  ! southmost y-boudnary
  if (rank.lt.nsubdomains_x) then
    dudtbackup(:,1,:) = dudt(:,1,:,na)
    dwdtbackup(:,1,:) = dwdt(:,1,:,na)
  ! northmost y-boundary
  else if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
    dudtbackup(:,1,:) = dudt(:,ny,:,na)
    dwdtbackup(:,1,:) = dwdt(:,ny,:,na)
  end if
end if

call advect2_mom_xy()
call advect2_mom_z()

! openybc: replace the dudt and dwdt at the northmost and
! the southmost y boundary back to the values
! before the original advective scheme
if (doopeny) then
  ! southmost y-boundary
  if (rank.lt.nsubdomains_x) then
    dudt(:,1,:,na) = dudtbackup(:,1,:)
    dwdt(:,1,:,na) = dwdtbackup(:,1,:)
  ! northmost y-boundary
  else if (rank.gt.(nsubdomains_x*(nsubdomains_y-1))) then
    dudt(:,ny,:,na) = dudtbackup(:,1,:)
    dwdt(:,ny,:,na) = dwdtbackup(:,1,:)
  end if
  !!! call the advection of momentum at the y boundary
  call advection_momentum()
end if

if(dostatis) then
	
  do k=1,nzm
   do j=1,ny
    do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na)-du(i,j,k,1)
     du(i,j,k,2)=dvdt(i,j,k,na)-du(i,j,k,2)
     du(i,j,k,3)=dwdt(i,j,k,na)-du(i,j,k,3)
    end do
   end do
  end do

  call stat_tke(du,tkeleadv)
  call stat_mom(du,momleadv)
  call setvalue(twleadv,nzm,0.)
  call setvalue(qwleadv,nzm,0.)
  call stat_sw1(du,twleadv,qwleadv)

endif

call t_stopf ('advect_mom')

end subroutine advect_mom

