
subroutine kurant

use vars
use sgs, only: kurant_sgs

implicit none

integer i, j, k, ncycle1(1),ncycle2(1)
real wm(nz)  ! maximum vertical wind velocity
real uhm(nz) ! maximum horizontal wind velocity
real cfl, cfl_sgs

! leo: stabillity check
! local
real :: wmax, uhmax, wcfl, uhcfl
! global
real :: gwmax, guhmax, gwcfl, guhcfl
! MPI buffers
real :: buffer(4), buffer1(4)

ncycle = 1
	
wm(nz)=0.
do k = 1,nzm
 wm(k) = maxval(abs(w(1:nx,1:ny,k)))
 uhm(k) = sqrt(maxval(u(1:nx,1:ny,k)**2+YES3D*v(1:nx,1:ny,k)**2))
end do
w_max=max(w_max,maxval(w(1:nx,1:ny,1:nz)))
u_max=max(u_max,maxval(uhm(1:nzm)))

cfl = 0.
do k=1,nzm
  cfl = max(cfl,uhm(k)*dt*sqrt((1./dx)**2+YES3D*(1./dy)**2), &
                   max(wm(k),wm(k+1))*dt/(dz*adzw(k)) )
end do

call kurant_sgs(cfl_sgs)
cfl = max(cfl,cfl_sgs)
	
ncycle = max(1,ceiling(cfl/0.7))

if(dompi) then
  ncycle1(1)=ncycle
  call task_max_integer(ncycle1,ncycle2,1)
  ncycle=ncycle2(1)
end if
if(ncycle.gt.4) then
   if(masterproc) print *,'the number of cycles exceeded 4.'
   call task_abort()
end if

!!! leo: stability check
! compute the local maximum w and uh, and the corresponding cfl
! initialize
wmax = 0.
uhmax = 0.
wcfl = 0.
uhcfl = 0.
do k=1,nzm
  wmax = max(wmax,wm(k))
  uhmax = max(uhmax,uhm(k))
  wcfl = max(wcfl,max(wm(k),wm(k+1))*dt/(dz*adzw(k)))
  uhcfl = max(uhcfl,uhm(k)*dt*sqrt((1./dx)**2+YES3D*(1./dy)**2))
end do
! global maximum w and uh, and the corresponding cfl
buffer(1) = wmax
buffer(2) = uhmax
buffer(3) = wcfl
buffer(4) = uhcfl
call task_max_real(buffer,buffer1,4)
gwmax = buffer1(1)
guhmax = buffer1(2)
gwcfl = buffer1(3)
guhcfl = buffer1(4)
if (masterproc) then
  write(*,*) 'Stability check:'
  write(*,*) 'maximum horizontal wind: ',guhmax
  write(*,*) 'maximum vertical wind: ',gwmax
  write(*,*) 'Courant number, horizontal contribution: ',guhcfl
  write(*,*) 'Courant number, vertical contribution: ',gwcfl
  write(*,*) 'maximum courant number allowed for variable timestep: 0.7'
end if

end subroutine kurant	
