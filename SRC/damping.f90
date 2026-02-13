
subroutine damping()

!  "Spange"-layer damping at the domain top region

use vars
use microphysics, only: micro_field, index_water_vapor
use params, only: doopeny ! leo: ybl
use grid, only: rank
use domain, only: nsubdomains_x, nsubdomains_y
implicit none

real tau_min	! minimum damping time-scale (at the top)
real tau_max    ! maxim damping time-scale (base of damping layer)
real damp_depth ! damping depth as a fraction of the domain height
parameter(tau_min=60., tau_max=450., damp_depth=0.4) ! leo (follows that of SPCAM)
!parameter(tau_min=60., tau_max=1800., damp_depth=0.3)
!parameter(tau_min=120., tau_max=7200., damp_depth=0.3)
real tau(nzm)   
integer i, j, k, n_damp

if(tau_min.lt.2*dt) then
   print*,'Error: in damping() tau_min is too small!'
   call task_abort()
end if

do k=nzm,1,-1
 if(z(nzm)-z(k).lt.damp_depth*z(nzm)) then 
   n_damp=nzm-k+1
 endif
end do

do k=nzm,nzm-n_damp,-1
 tau(k) = tau_min *(tau_max/tau_min)**((z(nzm)-z(k))/(z(nzm)-z(nzm-n_damp)))
 tau(k)=1./tau(k)
end do

do k = nzm, nzm-n_damp, -1
   do j=1,ny
    do i=1,nx
      dudt(i,j,k,na)= dudt(i,j,k,na)-(u(i,j,k)-u0(k)) * tau(k)
      dvdt(i,j,k,na)= dvdt(i,j,k,na)-(v(i,j,k)-v0(k)) * tau(k)
      dwdt(i,j,k,na)= dwdt(i,j,k,na)-w(i,j,k) * tau(k)
!      t(i,j,k)= t(i,j,k)-dtn*(t(i,j,k)-t0(k)) * tau(k)
!      micro_field(i,j,k,index_water_vapor)= micro_field(i,j,k,index_water_vapor)- &
!                                    dtn*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k)) * tau(k)
    end do! i 
   end do! j
end do ! k

! leo: ybl
! seperately treat the dvdt at y=nyp1
if (doopeny) then
  ! MPI, northmost grid only
  if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
    do k = nzm, nzm-n_damp, -1
      do i=1,nx
        dvdt(i,ny+1,k,na) = dvdt(i,ny+1,k,na) - (v(i,ny+1,k)-v0(k)) * tau(k)
      end do
    end do
  end if
  ! backup the value of dvdt
  if (rank.lt.nsubdomains_x) then
    dvdtsouthmost(:,1,:) = dvdt(:,1,:,na)
  elseif (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
    dvdtnorthmost(:,1,:) = dvdt(:,ny+1,:,na)
  end if
end if

end subroutine damping
