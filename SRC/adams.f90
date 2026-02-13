
subroutine adams

!       Adams-Bashforth scheme

use vars
use domain, only: nsubdomains_x, nsubdomains_y ! leo: ybl
use grid, only: nyp1
use params, only: doopeny ! leo: ybl

implicit none

real dtdx, dtdy, dtdz, rhox, rhoy, rhoz	
integer i,j,k


dtdx = dtn/dx
dtdy = dtn/dy
dtdz = dtn/dz

! back up the v at the current timestep
! for advect_scalar in first order center differencing scheme
if (doopeny) then
  ubackup(:,:,:) = u(:,:,:)
  vbackup(:,:,:) = v(:,:,:)
  wbackup(:,:,:) = w(:,:,:)
end if

! the dvdt at the boundaries already have the openbc and the damping terms
do k=1,nzm 
  rhox = rho(k)*dtdx  
  rhoy = rho(k)*dtdy  
  rhoz = rhow(k)*dtdz  
  do j=1,ny
   do i=1,nx
  
     dudt(i,j,k,nc) = u(i,j,k) + dt3(na) & 
              *(at*dudt(i,j,k,na)+bt*dudt(i,j,k,nb)+ct*dudt(i,j,k,nc))
	   
     dvdt(i,j,k,nc) = v(i,j,k) + dt3(na) &
              *(at*dvdt(i,j,k,na)+bt*dvdt(i,j,k,nb)+ct*dvdt(i,j,k,nc))
	   
     dwdt(i,j,k,nc) = w(i,j,k) + dt3(na) &
              *(at*dwdt(i,j,k,na)+bt*dwdt(i,j,k,nb)+ct*dwdt(i,j,k,nc))
     
    ! second order along t dims (dimensionless)
     u(i,j,k) = 0.5*(u(i,j,k)+dudt(i,j,k,nc)) * rhox
     v(i,j,k) = 0.5*(v(i,j,k)+dvdt(i,j,k,nc)) * rhoy
     misc(i,j,k) = 0.5*(w(i,j,k)+dwdt(i,j,k,nc)) 
     w(i,j,k) = 0.5*(w(i,j,k)+dwdt(i,j,k,nc)) * rhoz
	   
   end do 
  end do 
end do 

! for the yBL
! seperately compute the nyp1 grid for v
if (doopeny) then
  ! for the southmost boundary:
  ! the dvdt already consist of the terms of open BL at y=1
  ! for the northmost boundary:
  if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
    do k=1,nzm
      rhoy = rho(k)*dtdy
      do i=1,nx
        dvdt(i,nyp1,k,nc) = v(i,nyp1,k) + dt3(na) &
                *(at*dvdt(i,nyp1,k,na)+bt*dvdt(i,nyp1,k,nb)+ct*dvdt(i,nyp1,k,nc))
        ! PS. no need to backup as dvdt(:,:,:,nc) will not be updated later
        !     the script update the dimensionless y-velocity v
        
        ! second order along t dims (dimensionless)
        v(i,nyp1,k) = 0.5*(v(i,nyp1,k)+dvdt(i,nyp1,k,nc)) * rhoy
      
      end do
    end do
  end if
end if

end subroutine adams

	
