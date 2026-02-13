module open_boundary_y

  !--------------------------------------------------------------------|
  ! open boundary condition at y boundary (i.e. y=1 and y=nyp1)        |
  ! for radiate out gravity wave                                       |
  ! reference:                                                         |
  ! 1. Klemp and Wilhelmson 1978,                                      |
  ! The Simulation of Three-Dimensional Convective Storm Dynamics      |
  ! https://doi.org/10.1175/1520-0469(1978)035<1070:TSOTDC>2.0.CO;2    |
  ! 2. Skamarock et. al. 2008,                                         |
  ! A Description of the Advanced Research WRF Version 3               |
  ! https://www2.mmm.ucar.edu/wrf/users/docs/technote/v3_technote.pdf  |
  ! written by Leo Chow (20230725)                                     |
  ! email: leonardotnchow@cuhk.edu.hk                                  |
  !--------------------------------------------------------------------|

  use vars
  use grid
  use domain, only: nsubdomains_x, nsubdomains_y
  use microphysics, only: nmicro_fields
  use sgs, only: nsgs_fields

  implicit none

  ! gravity wave speed
  real :: cb = 30. ! m/s

  !----------------------------------------------------------------------
  ! other variables
  integer :: nprev, ncurr

  !----------------------------------------------------------------------
  ! variables for previous and current timestep 
  ! (y:1,2 -> 1,2)
  ! (y:3,4 -> ny-1, ny for variables except v)
  ! (y:3,4 -> ny, ny+1 for v)
  ! (nprev = previous timestep, ncurr = current timestep)
  real :: uprev   (dimx1_u:dimx2_u,6,nzm,2) ! x-wind
  real :: vprev   (dimx1_v:dimx2_v,6,nzm,2) ! y-wind
  real :: wprev   (dimx1_w:dimx2_w,6,nz,2 ) ! z-wind
  real :: tprev   (dimx1_s:dimx2_s,6,nzm,2) ! liquid/ice water static energy

  !----------------------------------------------------------------------
  ! microphysics variables
  real :: micro_fieldprev   (dimx1_s:dimx2_s,6,nzm,9,2)

  !----------------------------------------------------------------------
  ! sgs variables
  real :: sgs_fieldprev     (dimx1_s:dimx2_s,6,nzm,nsgs_fields,2)

  !----------------------------------------------------------------------
  ! dtn (dt in variables AB-scheme cycle)
  real :: dtnprev(2)


  CONTAINS

!========================================================================
!!! replace the y-momentum equation at j = 1 and j = nyp1 
! by the open boundary 
! 1st order
  subroutine ymom_openbl()

    implicit none

    real :: advvb, advvt
    integer :: it, jt
    integer :: i,j,k

    ! southern boundary
    if (rank.lt.nsubdomains_x) then
      ! zero out the dvdt at the southmost grid
      dvdt(:,1,:,na) = 0.
      ! from the backup, retain some of the terms (lateral and vertical damping)
      dvdt(:,1,:,na) = dvdtsouthmost(:,1,:)
      ! open boundary
      do i=1,nx
        do k=1,nzm
          ! compute the outflow condition at j = 1
          advvb = min(0.,v(i,1,k)-cb)
          ! compute the source term (advection only now)
          dvdt(i,1,k,na) = dvdt(i,1,k,na) - advvb * &
                                            (v(i,2,k) - v(i,1,k)) / dy
        end do
      end do
    end if

    ! northern boundary
    if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
      ! zero out the dvdt at the northmost grid
      dvdt(:,nyp1,:,na) = 0.
      ! from the backup, retain some of the terms (damping terms)
      dvdt(:,nyp1,:,na) = dvdtnorthmost(:,1,:)
      do i=1,nx
        do k=1,nzm
          ! compute the outflow condition at j = nyp1
          advvt = max(0.,v(i,nyp1,k)+cb)
          ! compute the source term (advection only now)
          dvdt(i,nyp1,k,na) = dvdt(i,nyp1,k,na) - advvt * &
                                                  (v(i,nyp1,k) - v(i,ny,k)) / dy
        end do
      end do
    end if

  end subroutine ymom_openbl

!====================================================================================
!!! change the advection of non-normal momentum 
! at the southmost and the northmost y boundaries 
! to be onesided differencing (upwind)
! e.g. for j= 1, use the value at j = 1 and j = 2
! the idea is not to allow influx/reflection of scalar variables due to advective
! dynamics at the boundaries

  subroutine advection_momentum()

    !--------------------------------------------------|
    ! first order (forward differencing at y = 0) and  |
    ! (backward differencing at y = ny)                |
    ! to estimate the boundary                         |
    ! noted that an upwind scheme is needed            |
    ! for stability (i.e. V-cb < 0 at y=0,             |
    !                     V+cb > 0 at y=ny )           |
    ! notation: Xue and Lin 2001, eq5 and 7            |
    ! https://doi.org/10.1175/1520-0493(2001)          |
    ! 129<0561:NEOAIF>2.0.CO;2                         |
    !--------------------------------------------------|

    implicit none

    integer :: i,j,k, kb, kc, kcu, ic, ib, jc, jb, jtemp
    real :: dx5, dy5, dz5, irho, rhoi
    real :: fuf(nx,2,nzm), fub(nx,2,nzm), fwf(nx,2,nzm), fwb(nx,2,nzm),&
            ff(nx,1,nzm)
    real :: fu(nx,1,nzm), fv(nx,1,nzm), fw(nx,1,nzm)

    dx5 = 0.5 / dx
!    dy5 = 0.5 / dy
    dz5 = 0.5 / dz

    !------------------------------ j = 1 ------------------------------------
    !!! j = 1
    if (rank.lt.nsubdomains_x) then

      !!! second-order center in space differencing
      ! with first-order advection 
      ! advective form, ignore the continuity of momentum (elastic correction)
      ! by applying the anelastic continuity equation 

      do k=1,nzm
        !!! advection in x dir
        do i=1,nx
          ic = i + 1
          ib = i - 1
          ! advection of u wind (second order-center differencing)
          fu(i,1,k) = dx5/2. * ((u(ic,1,k)+u(i,1,k)) * (u(ic,1,k) - u(i,1,k)) + &
                                (u(i,1,k)+u(ib,1,k)) * (u(i,1,k) - u(ib,1,k)))
          ! advection of v wind (first order-center differencing) 
          ! (no need as it is y = 1 is y boundary, no advection)
          !!! add the advection term to the d()dt
          dudt(i,1,k,na) = dudt(i,1,k,na) - fu(i,1,k)
        end do
        !!! advection in y dir (first-order upwind (or zero-advection))
        do i=1,nx
          fu(i,1,k) = min(0.,v(i,1,k)) * (u(i,2,k) - u(i,1,k)) / dy
          !!! add the advection term to the d()dt
          dudt(i,1,k,na) = dudt(i,1,k,na) - fu(i,1,k)
        end do
      end do
      ! advection of z-wind in x and y dir
      do k=2,nzm
        kc = k + 1
        kb = k - 1
        irho = 1./(rhow(k)*adzw(k))
        do i=1,nx
          ic = i + 1
          ib = i - 1
          ! advection of z wind in x-dir (second order-center differencing)
          fw(i,1,k) = dx5/2. * ((u(ic,1,k)*rho(k)*adz(k)+u(ic,1,kb)*rho(kb)*adz(kb)) * &
                                (w(ic,1,k)-w(i,1,k)) + &
                                (u(i,1,k)*rho(k)*adz(k)+u(i,1,kb)*rho(kb)*adz(kb)) * &
                                (w(i,1,k)-w(ib,1,k)))
          ! add the advection term to the d()dt
          dwdt(i,1,k,na) = dwdt(i,1,k,na) - irho*fw(i,1,k)
          ! advection of z wind in y-dir (first order-upwind)
          fw(i,1,k) = 1./dy * min(0.,v(i,1,k))*rho(k)*adz(k) * &
                              (w(i,2,k) - w(i,1,k))
          ! add the advection term to the d()dt
          dwdt(i,1,k,na) = dwdt(i,1,k,na) - irho*fw(i,1,k)
        end do
      end do
      !!!------------------------------------------
      !!! advection in z dir
      ! x-wind and y-wind
      ! at k = 1 and k = nzm
      do i=1,nx
        ic = i + 1
        ! k = 1 (PS. need min? for upwind)
        fu(i,1,1) = 1./(dz*2.) * (min(0.,(w(ic,1,2)+w(i,1,2))*rhow(2)) &
                             * (u(i,1,2) - u(i,1,1)))
        ! k = nzm (PS. need max? for upwind)
        fu(i,1,nzm) = 1./(dz*2.) * (max(0.,(w(ic,1,nzm)+w(i,1,nzm))*rhow(nzm)) * &
                                 (u(i,1,nzm) - u(i,1,nzm-1)))
      end do
      ! other levels
      do k=2,nzm-1
        kc = k + 1
        kb = k - 1
        do i=1,nx
          ic = i + 1
          ! advection of x-wind along z dir (second order-center differencing)
          fu(i,1,k) = dz5/2.*((w(ic,1,kc)+w(i,1,kc))*rhow(kc) * (u(i,1,kc) - u(i,1,k)) + &
                               (w(ic,1,k)+w(i,1,k))*rhow(k) * (u(i,1,k) - u(i,1,kb)))
          ! advection of y-wind along z dir (first order-center differencing)
          ! (no need as y = 1 is y boundary, no need advection)
        end do
      end do
      !!! add the advection term to the d()dt
      do k=1,nzm
        irho = 1./(rho(k)*adz(k))
        do i=1,nx
          dudt(i,1,k,na) = dudt(i,1,k,na) - irho*fu(i,1,k)
        end do
      end do
      !!! advection of z-wind at z dir (second-order center differencing)
!      fw(:,1,1) = 0.  ! these 2 will not appear later at dwdt, so just comment it,
!      fw(:,1,nz) = 0. 
      do k=2,nzm
        kc = k + 1
        kb = k - 1
        irho = 1./(rhow(k)*adzw(k))
        do i=1,nx
          fw(i,1,k) = dz5/2.*((w(i,1,kc)*rhow(kc)+w(i,1,k)*rhow(k)) * &
                              (w(i,1,kc) - w(i,1,k)) + &
                              (w(i,1,k)*rhow(k)+w(i,1,kb)*rhow(kb)) * &
                              (w(i,1,k) - w(i,1,kb)))
                              ! PS problem at k = 2 and k = nzm?
          !!! add the advection term to the d()dt
          dwdt(i,1,k,na) = dwdt(i,1,k,na) - irho*fw(i,1,k)
        end do
      end do
    end if ! if (rank.lt.nsubdomains_x)
        
    !------------------------------ j = ny ------------------------------------
    !!! j = ny
    if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then

      !!! second-order center in space differencing
      ! with first-order / second-order advection 
      ! advective form, ignore the continuity of momentum (elastic correction)
      ! by applying the anelastic continuity equation 

      do k=1,nzm
        !!! advection in x dir
        do i=1,nx
          ic = i + 1
          ib = i - 1
          ! advection of u wind (second order-center differencing)
          fu(i,1,k) = dx5/2. * ((u(ic,ny,k)+u(i,ny,k)) * (u(ic,ny,k) - u(i,ny,k)) + &
                                (u(i,ny,k)+u(ib,ny,k)) * (u(i,ny,k) - u(ib,ny,k)))
          ! advection of v wind (first order-center differencing)
          ! (no need as y = ny is not the y-boundary, can use 2nd-order center differencing)
          ! advection of w wind (first order center differencing)
          !!! add the advection term to the d()dt
          dudt(i,ny,k,na) = dudt(i,ny,k,na) - fu(i,1,k)
        end do
        !!! advection in y dir (first order-upwind differencing, or zero-advection)
        do i=1,nx
          fu(i,1,k) = max(0.,v(i,ny+1,k)) * (u(i,ny,k) - u(i,ny-1,k)) / dy
          !!! add the advection term to the d()dt
          dudt(i,ny,k,na) = dudt(i,ny,k,na) - fu(i,1,k)
        end do
      end do
      ! advection of z-wind in x and y dir
      do k=2,nzm
        kc = k + 1
        kb = k - 1
        irho = 1./(rhow(k)*adzw(k))
        do i=1,nx
          ic = i + 1
          ib = i - 1
          ! advection of z-wind in x dir (second order-center differencing)
          fw(i,1,k) = dx5/2. * ((u(ic,ny,k)*rho(k)*adz(k)+u(ic,ny,kb)*rho(kb)*adz(kb))*&
                                (w(ic,ny,k)-w(i,ny,k)) + &
                                (u(i,ny,k)*rho(k)*adz(k)+u(i,ny,kb)*rho(kb)*adz(kb))*&
                                (w(i,ny,k)-w(ib,ny,k)))
          ! add the advection term to the d()dt
          dwdt(i,ny,k,na) = dwdt(i,ny,k,na) - irho*fw(i,1,k)
          ! advection of z-wind in y dir (first order-upwind)
          fw(i,1,k) = 1./dy * (max(0.,v(i,ny,k))*rho(k)*adz(k)*&
                               (w(i,ny,k)-w(i,ny-1,k)))
          ! add the advection term to the d()dt
          dwdt(i,ny,k,na) = dwdt(i,ny,k,na) - irho*fw(i,1,k)
        end do
      end do
      !!!------------------------------------------
      !!! advection in z dir
      ! x-wind and y-wind
      ! at k = 1 and k = nzm
      do i=1,nx
        ic = i + 1
        ! k = 1 (PS. need min?)
        fu(i,1,1) = 1./(dz*2.) * (min(0.,(w(ic,ny,2)+w(i,ny,2))*rhow(2)) * &
                                (u(i,ny,2) - u(i,ny,1)))
        ! k = nzm (PS. need max?)
        fu(i,1,nzm) = 1./(dz*2.) * (max(0.,(w(ic,ny,nzm)+w(i,ny,nzm))*rhow(nzm)) * &
                                  (u(i,ny,nzm) - u(i,ny,nzm-1)))
      end do
      ! other levels
      do k=2,nzm-1
        kc = k + 1
        kb = k - 1
        do i=1,nx
          ic = i + 1
          ! advection of x-wind in z dir (second order-center differencing)
          fu(i,1,k) = dz5/2.*((w(ic,ny,kc)+w(i,ny,kc))*rhow(kc) * &
                              (u(i,ny,kc) - u(i,ny,k)) + &
                              (w(ic,ny,k)+w(i,ny,k))*rhow(k) * &
                              (u(i,ny,k) - u(i,ny,kb)))
          ! advection of y-wind in z dir (first order-center differencing)
          ! (no need as y = ny is not the y-boundary, can use 2nd-order center differencing)
        end do
      end do
      !!! add the advection term to the d()dt
      do k=1,nzm
        irho = 1./(rho(k)*adz(k))
        do i=1,nx
          dudt(i,ny,k,na) = dudt(i,ny,k,na) - irho*fu(i,1,k)
        end do
      end do
      !!! advection of z-wind at z dir (second order-center differencing)
!      fw(:,1,1) = 0.  ! these 2 will not appear later at dwdt, so just comment it,
!      fw(:,1,nz) = 0. 
      do k=2,nzm
        kc = k + 1
        kb = k - 1
        irho = 1./(rhow(k)*adzw(k))
        do i=1,nx
          fw(i,1,k) = dz5/2.*((w(i,ny,kc)*rhow(kc)+w(i,ny,k)*rhow(k)) * &
                              (w(i,ny,kc) - w(i,ny,k)) + &
                              (w(i,ny,k)*rhow(k)+w(i,ny,kb)*rhow(kb)) * &
                              (w(i,ny,k) - w(i,ny,kb))) 
                              ! PS problem at k = 2 and k = nzm?
          !!! add the advection term to the d()dt
          dwdt(i,ny,k,na) = dwdt(i,ny,k,na) - irho*fw(i,1,k)
        end do
      end do
    end if ! if (rank.lt.nsubdomains_x*(nsubdomains_y-1)-1)

  end subroutine advection_momentum

!====================================================================================
!!! change the advection of non-normal scalar
! at the southmost and the northmost y boundaries 
! to be onesided differencing 
! e.g. for j = 1, use the value at j = 1 and j = 2
! the idea is not to allow influx/reflection of scalar variables due to advective
! dynamics at the boundaries

  subroutine advection_scalar(f,fprev,uinout,vinout,winout,bufferinout)

    !--------------------------------------------------|
    ! first order (forward differencing at y = 0) and  |
    ! (backward differencing at y = nyp1)              |
    ! to estimate the boundary                         |
    ! noted that an upwind scheme is needed            |
    ! for stability (i.e. V-cb < 0 at y=0,             |
    !                     V+cb > 0 at y=ny )           |
    ! notation: Xue and Lin 2001, eq5 and 7            |
    ! https://doi.org/10.1175/1520-0493(2001)          |
    ! 129<0561:NEOAIF>2.0.CO;2                         |
    !                                                  |
    ! 20230816:                                        |
    ! P.s. in MPDATA, y = 1,2,3 and ny-2, ny-1, ny     |
    !      will be affected by the open boundary at    |
    !      y = 1 and y = ny+1, so we need to use       |
    !      second order forward differencing to compute|
    !      the values at these grids                   |
    ! 20230829:                                        |
    ! update the time-differencing scheme to           |
    ! leapfrog scheme                                  | 
    !--------------------------------------------------|

    implicit none

    integer :: i,j,k,jb,jc, jtemp, jtemp1, ib,ic
    integer :: kb, kc, kcu ! z-dir
    real :: ffx(1:nx,3,nzm), ffy(nx,3,nzm), ffz(nx,3,nzm)

    real, intent(inout) :: f(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),&
                           fprev(dimx1_s:dimx2_s,6,nzm,2),&
                           uinout(dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm), &
                           vinout(dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm), &
                           winout(dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz)
    real, intent(out) :: bufferinout(dimx1_s:dimx2_s,6,nzm) ! for output
                                                            ! y=1 is y=1, y=2 is y=ny

    real :: dx5, dy5, dz5
    real :: rhoi, irho

    ! debug
    integer :: it, jt

    dx5 = 0.5 / dx
    dy5 = 0.5 / dy
    dz5 = 0.5 / dz

    !!! y = 1
    !------------------------------ y = 1 ---------------------------------
    if (rank.lt.nsubdomains_x) then
      do k=1,nzm
        do i=1,nx
          ic = i + 1
          ib = i - 1
          ! advection in x dir (second order-center differencing)
          ffx(i,1,k) = dx5 * (uinout(ic,1,k)*(f(ic,1,k) - f(i,1,k)) + &
                              uinout(i,1,k)*(f(i,1,k) - f(ib,1,k)))
          ! advection in y dir (first order-upwind)
          ffy(i,1,k) = 1./dy * min(0.,vinout(i,1,k)) * (f(i,2,k) - f(i,1,k))
        end do
      end do
      !!! advection in z dir
      do i=1,nx
        ! k = 1 and k = nzm (upwind differencing)
        ffz(i,1,1) = 1./dz*rhow(2) * (min(0.,winout(i,1,2))*&
                                      (f(i,1,2) - f(i,1,1)))
        ffz(i,1,nzm) = 1./dz*rhow(nzm) * (max(0.,winout(i,1,nzm))*&
                                          (f(i,1,nzm) - f(i,1,nzm-1)))
        ! other levels
        do k=2,nzm-1
          kc = k + 1
          kb = k - 1
          ffz(i,1,k) = dz5 * (winout(i,1,kc)*rhow(kc) * (f(i,1,kc) - f(i,1,k)) + &
                               winout(i,1,k)*rhow(k) * (f(i,1,k) - f(i,1,kb)))
        end do
      end do
      !!! leapfrog scheme / forward in time scheme
      if (nstep.eq.1) then
        ! nstep = 1, forward in time
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do i=1,nx
            bufferinout(i,1,k) = f(i,1,k) - &
                                 dtn * (ffx(i,1,k) + ffy(i,1,k) + rhoi*ffz(i,1,k))
            ! remove the truncation error, which result in -ve scalar
            bufferinout(i,1,k) = max(0.,bufferinout(i,1,k))
          end do
        end do
      else
        ! nstep.ne.1, leapfrog for u and w direction
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do i=1,nx
            bufferinout(i,1,k) = fprev(i,1,k,nprev) - &
                                 (dtn+dtnprev(nprev)) * (ffx(i,1,k) + ffy(i,1,k) + rhoi*ffz(i,1,k))
            ! remove the truncation error, which result in -ve scalar
            bufferinout(i,1,k) = max(0.,bufferinout(i,1,k))
          end do
        end do
      end if ! if (nstep.eq.1)
    end if ! if (rank.eq.nsubdomains_x) (y=1)
    
    !!! y = ny
    !------------------------------ y = ny --------------------------------
    if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
      do k=1,nzm
        do i=1,nx
          ic = i + 1
          ib = i - 1
          ! advection in x dir (second order-center differencing)
          ffx(i,1,k) = dx5 * (uinout(ic,ny,k)*(f(ic,ny,k) - f(i,ny,k)) + &
                              uinout(i,ny,k)*(f(i,ny,k) - f(ib,ny,k)))
          ! advection in y dir (first order-upwind)
          ffy(i,1,k) = 1./dy * max(0.,vinout(i,ny+1,k)) * (f(i,ny,k) - f(i,ny-1,k))
        end do
      end do
      !!! advection in z dir
      do i=1,nx
        ! k = 1 and k = nzm (upwind differencing)
        ffz(i,1,1) = 1./dz*rhow(2) * (min(0.,winout(i,ny,2))*&
                                      (f(i,ny,2) - f(i,ny,1)))
        ffz(i,1,nzm) = 1./dz*rhow(nzm) * (max(0.,winout(i,ny,nzm))*&
                                          (f(i,ny,nzm) - f(i,ny,nzm-1)))
        ! other levels
        do k=2,nzm-1
          kc = k + 1
          kb = k - 1
          ffz(i,1,k) = dz5 * (winout(i,ny,kc)*rhow(kc) * (f(i,ny,kc) - f(i,ny,k)) + &
                               winout(i,ny,k)*rhow(k) * (f(i,ny,k) - f(i,ny,kb)))
        end do
      end do
      !!! leapfrog scheme / forward in time scheme
      if (nstep.eq.1) then
        ! nstep = 1, forward in time
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do i=1,nx
            bufferinout(i,6,k) = f(i,ny,k) - &
                                 dtn * (ffx(i,1,k) + ffy(i,1,k) + rhoi*ffz(i,1,k))
            ! remove the truncation error, which result in -ve scalar
            bufferinout(i,6,k) = max(0.,bufferinout(i,6,k))
          end do
        end do
      else
        ! nstep.ne.1, leapfrog
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do i=1,nx
            bufferinout(i,6,k) = fprev(i,6,k,nprev) - &
                                 (dtn+dtnprev(nprev)) * (ffx(i,1,k) + ffy(i,1,k) + rhoi*ffz(i,1,k))
            ! remove the truncation error, which result in -ve scalar
            bufferinout(i,6,k) = max(0.,bufferinout(i,6,k))
          end do
        end do
      end if ! if (nstep.eq.1)
    end if ! if (rank.eq.(nsubdomains_x*(nsubdomains_y)-1)-1) y = ny

    !!! y = 2, 3
    !------------------------------ y = 2,3 -------------------------------
    if (rank.lt.nsubdomains_x) then
      do k=1,nzm
        do i=1,nx
          ic = i + 1
          ib = i - 1
          do j=2,3
            jc = j + 1
            jb = j - 1
            ! advection in x dir (second order-center differencing)
            ffx(i,j,k) = dx5 * (uinout(ic,j,k)*(f(ic,j,k) - f(i,j,k)) + &
                                uinout(i,j,k)*(f(i,j,k) - f(ib,j,k)))
            ! advection in y dir (second order-center differencing)
            ffy(i,j,k) = dy5 * (vinout(i,jc,k)*(f(i,jc,k) - f(i,j,k)) + &
                                vinout(i,j,k)*(f(i,j,k) - f(i,jb,k)))
          end do
        end do
      end do
      !!! advection in z dir
      do i=1,nx
        do j=2,3
          ! k = 1 and k = nzm (upwind differencing)
          ffz(i,j,1) = 1./dz*rhow(2) * (min(0.,winout(i,j,2))*&
                                       (f(i,j,2) - f(i,j,1)))
          ffz(i,j,nzm) = 1./dz*rhow(nzm) * (max(0.,winout(i,j,nzm))*&
                                           (f(i,j,nzm) - f(i,j,nzm-1)))
          ! other levels
          do k=2,nzm-1
            kc = k + 1
            kb = k - 1
            ffz(i,j,k) = dz5 * (winout(i,j,kc)*rhow(kc) * (f(i,j,kc) - f(i,j,k)) + &
                                 winout(i,j,k)*rhow(k) * (f(i,j,k) - f(i,j,kb)))
          end do
        end do
      end do
      !!! leapfrog scheme / forward in time scheme
      if (nstep.eq.1) then
        ! nstep = 1, forward in time
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do j=2,3 ! bufferinout(:,2-3,:)
            do i=1,nx
              bufferinout(i,j,k) = f(i,j,k) - &
                                   dtn * (ffx(i,j,k) + ffy(i,j,k) + rhoi*ffz(i,j,k))
              ! remove the truncation error, which result in -ve scalar
              bufferinout(i,j,k) = max(0.,bufferinout(i,j,k))
            end do
          end do
        end do
      else
        ! nstep.ne.1, leapfrog
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do j=2,3 ! bufferinout, fprev(:,2-3,:)
            do i=1,nx
              bufferinout(i,j,k) = fprev(i,j,k,nprev) - &
                                   (dtn+dtnprev(nprev)) * (ffx(i,j,k) + ffy(i,j,k) + rhoi*ffz(i,j,k))
              ! remove the truncation error, which result in -ve scalar
              bufferinout(i,j,k) = max(0.,bufferinout(i,j,k))
            end do
          end do
        end do
      end if ! if (nstep.eq.1)
    end if ! if (rank.eq.nsubdomains_x) (y=2,3)

    !!! y = ny-1,ny-2
    !-------------------------- y = ny-1, ny-2 ----------------------------
    if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
      do k=1,nzm
        do j=ny-2,ny-1
          jc = j + 1
          jb = j - 1
          jtemp = j - ny + 4 ! ffx, ffy(:,2-3,:)
          do i=1,nx
            ic = i + 1
            ib = i - 1
            ! advection in x dir (second order-center differencing)
            ffx(i,jtemp,k) = dx5 * (uinout(ic,j,k)*(f(ic,j,k) - f(i,j,k)) + &
                                    uinout(i,j,k)*(f(i,j,k) - f(ib,j,k)))
            ! advection in y dir (second order-center differencing)
            ffy(i,jtemp,k) = dy5 * (vinout(i,jc,k)*(f(i,jc,k) - f(i,j,k)) + &
                                    vinout(i,j,k)*(f(i,j,k) - f(i,jb,k)))
          end do
        end do
      end do
      !!! advection in z dir
      do i=1,nx
        do j=ny-2,ny-1
          jtemp = j - ny + 4 ! ffz(:,2-3,:)
          ! k = 1 and k = nzm (upwind differencing)
          ffz(i,jtemp,1) = 1./dz*rhow(2) * (min(0.,winout(i,j,2))*&
                                            (f(i,j,2) - f(i,j,1)))
          ffz(i,jtemp,nzm) = 1./dz*rhow(nzm) * (max(0.,winout(i,j,nzm)) * &
                                                (f(i,j,nzm) - f(i,j,nzm-1)))
          ! other levels
          do k=2,nzm-1
            kc = k + 1
            kb = k - 1
            ffz(i,jtemp,k) = dz5 * (winout(i,j,kc)*rhow(kc) * (f(i,j,kc) - f(i,j,k)) + &
                                     winout(i,j,k)*rhow(k) * (f(i,j,k) - f(i,j,kb)))
          end do
        end do
      end do
      !!! leapfrog scheme / forward in time scheme
      if (nstep.eq.1) then
        ! nstep = 1, forward in time
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do j=ny-2,ny-1
            jtemp = j - ny + 4 ! ffx, ffy, ffz(:,2-3,:)
            jtemp1 = jtemp + 2 ! bufferinout(:,4-5,k)
            do i=1,nx
              bufferinout(i,jtemp1,k) = f(i,j,k) - &
                                        dtn * (ffx(i,jtemp,k) + &
                                               ffy(i,jtemp,k) + &
                                               rhoi*ffz(i,jtemp,k))
              ! remove the truncation error, which result in -ve scalar
              bufferinout(i,jtemp1,k) = max(0.,bufferinout(i,jtemp1,k))
            end do
          end do
        end do
      else
        ! nstep.ne.1, leapfrog
        do k=1,nzm
          rhoi = 1./(rho(k)*adz(k))
          do j=ny-2, ny-1
            jtemp = j - ny + 4 ! ffx, ffy, ffz(:,2-3,:)
            jtemp1 = jtemp + 2 ! bufferinout(:,4-5,:)
            do i=1,nx
              bufferinout(i,jtemp1,k) = fprev(i,jtemp1,k,nprev) - &
                                        (dtn+dtnprev(nprev)) * (ffx(i,jtemp,k) + &
                                                                ffy(i,jtemp,k) + &
                                                                rhoi*ffz(i,jtemp,k))
              ! remove the truncation error, which result in -ve scalar
              bufferinout(i,jtemp1,k) = max(0.,bufferinout(i,jtemp1,k))
            end do
          end do
        end do
      end if ! if (nstep.eq.1)
    end if ! if (rank.eq.(nsubdomains_x*(nsubdomains_y)-1)-1) y = ny-1,ny-2

  ! debug
!  do j=1,ny
!    if (rank.eq.28.and.j.gt.1077) then
!      call task_rank_to_index(rank,it,jt)
!      write(*,*) 'before advection_scalar: ',rank,it+11,jt+j,'2',f(11,j,2), vinout(11,j,2)
!    end if
!  end do
!    
!  ! debug
!  if (rank.eq.8) then
!    call task_rank_to_index(rank,it,jt)
!    write(*,*) 'advection_scalar: ', rank, it+11,'1','2',bufferinout(11,1,2)
!  else if (rank.eq.28) then
!    call task_rank_to_index(rank,it,jt)
!    write(*,*) 'advection_scalar: ', rank, it+11,'1080','2',bufferinout(11,6,2)
!  end if

  end subroutine advection_scalar

!===============================================================================
  subroutine prevbackup()

    !-------------------------------------|
    ! back up the prognostic variables    |
    ! at the current timestep             |
    !-------------------------------------|
    
    use microphysics, only: micro_field
    use sgs, only: sgs_field

    implicit none

    integer :: jtemp,j

    if (rank.lt.nsubdomains_x) then
      do j=1,3
        ! prognostic variables
        uprev(:,j,:,ncurr) = u(:,j,:)
        vprev(:,j,:,ncurr) = v(:,j,:)
        wprev(:,j,:,ncurr) = w(:,j,:)
        tprev(:,j,:,ncurr) = t(:,j,:)

        ! microphysics variables in MICRO_M2005
        micro_fieldprev(:,j,:,:,ncurr) = micro_field(:,j,:,:)

        ! sgs variables
        sgs_fieldprev(:,j,:,:,ncurr) = sgs_field(:,j,:,:)
      end do
    else if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
      do j=ny-2,ny
        jtemp = j - ny + 6 ! 4-6

        ! prognostic variables
        uprev(:,jtemp,:,ncurr) = u(:,j,:)
        vprev(:,jtemp,:,ncurr) = v(:,j+1,:) ! ny-1, ny, ny+1
        wprev(:,jtemp,:,ncurr) = w(:,j,:)
        tprev(:,jtemp,:,ncurr) = t(:,j,:)

        ! microphysics variables in MICRO_M2005
        micro_fieldprev(:,jtemp,:,:,ncurr) = micro_field(:,j,:,:)

        ! sgs variables
        sgs_fieldprev(:,jtemp,:,:,ncurr) = sgs_field(:,j,:,:)
      end do
    end if

    ! backup the dtn in current variable Adam-Bashforth scheme
    dtnprev(ncurr) = dtn

  end subroutine prevbackup

!===============================================================================
  subroutine openbl_init()

    !----------------------------------------|
    ! initialize the ncurr and nprev         |
    ! at each timestep,                      |
    ! if nstep = 1, ncurr = 2 and nprev = 1  |
    ! other, shift ncurr and nprev           |
    !----------------------------------------|

    implicit none

    integer :: ntemp

    if (nstep.eq.1) then
      ncurr = 2
      nprev = 1
    else
      ntemp = ncurr
      ncurr = nprev
      nprev = ntemp
    end if
  
  end subroutine openbl_init

!===============================================================================
  subroutine openbc_diagnose()

    !-------------------------------------------|
    ! diagnose the stability of the model       |
    ! related to the open boundary condition    |
    ! inposed                                   |
    ! by finding out the global max and min     |
    ! w and its location                        |
    ! at each vertical level k                  |
    !-------------------------------------------|

    use grid, only: rank, nx,ny,nzm
    use vars, only: w

    implicit none

    include 'mpif.h'

    integer :: i,j,k
    integer :: itt, jtt
    ! results
    real :: wmax(nzm-1), wmin(nzm-1)
    integer :: wmaxloc(2,nzm-1), wminloc(2,nzm-1) ! i and j of the local maximum 
                                                  ! and local minimum value
    integer :: irankmin(nzm-1),irankmax(nzm-1) ! store the rank of the output
    ! mpi buffers
    real :: buffermin(2,nzm-1),buffermax(2,nzm-1),&
            buffer1min(2,nzm-1), buffer1max(2,nzm-1)
            ! 1st is the value of max/min, second is the rankof the max/min

    ! compute the local maximum and minimum of w
    do k=2,nzm
      ! local max/min
      buffermin(1,k-1) = MINVAL(w(1:nx,1:ny,k))
      buffermax(1,k-1) = MAXVAL(w(1:nx,1:ny,k))
      ! store the rank to the buffer
      buffermin(2,k-1) = rank ! convert it to real
      buffermax(2,k-1) = rank 
    end do

    ! get the i and j of the local maximum and minimum
    do k=2,nzm
      wminloc(:,k-1) = MINLOC(w(1:nx,1:ny,k))
      wmaxloc(:,k-1) = MAXLOC(w(1:nx,1:ny,k))
    end do

    ! compute the global maximum and minimum value
    ! and the corresponding rank
    call MPI_ALLREDUCE(buffermin, buffer1min, nzm-1,&
                       MPI_2REAL, MPI_MINLOC, MPI_COMM_WORLD)
    call MPI_ALLREDUCE(buffermax,buffer1max, nzm-1,&
                       MPI_2REAL, MPI_MAXLOC, MPI_COMM_WORLD)

    ! distribute the results to the arrays
    do k=1,nzm-1
      ! max/min value
      wmin(k) = buffer1min(1,k)
      wmax(k) = buffer1max(1,k)
      ! rank of the max/min
      irankmin(k) = buffer1min(2,k) ! convert it back to integer
      irankmax(k) = buffer1max(2,k)
    end do

    !!! write the results (min)
    if (mod(nstep,5).eq.0) then ! write the results by some frequency of NSTEP
      if (masterproc) then
        write(*,*) 'maximum and minimum z-wind diagnostics'
        write(*,*) ' rank     i     j     k        min'
      end if
      do k=nzm-1,1,-1
        call task_barrier() ! synchronize the speed for each rank
        if (rank.eq.irankmin(k)) then
          call task_rank_to_index(rank,itt,jtt)
          write(*,'(4i6,1f11.4)') irankmin(k),itt+wminloc(1,k),jtt+wminloc(2,k),k+1,wmin(k)
        end if
      end do
      call task_barrier() ! synchronize the speed for each rank
      if (masterproc) then
      !!! write the results (max)
        write(*,*) ' rank     i     j     k        max'
      end if
      do k=nzm-1,1,-1
        call task_barrier() ! synchronize the speed for each rank
        if (rank.eq.irankmax(k)) then
          call task_rank_to_index(rank,itt,jtt)
          write(*,'(4i6,1f11.4)') irankmax(k),itt+wmaxloc(1,k),jtt+wmaxloc(2,k),k+1,wmax(k)
        end if
      end do
    end if

  end subroutine openbc_diagnose

!=========================================================================
  subroutine write_restart_openbc()

    !====================================|
    ! store the restart file             |
    ! containing the variables defined   |
    ! at the module in openybc           |
    !====================================|

    implicit none

    character(len=4) :: rankchar
    integer :: irank
    character(len=256) :: filename

    if (masterproc) then
      write(*,*) 'Writing restart file of openbc...'
    end if

    if (restart_sep) then
      write(*,*) 'Currently not supporting restart_sep'
      return
    else
      write(rankchar,'(i4)') nsubdomains
      filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_openybc_restart.bin'
      ! loop to write
      do irank=0,nsubdomains-1
        call task_barrier()
        if (irank.eq.rank) then
          if (masterproc) then
            open(86,file=trim(filename),status='unknown',form='unformatted')
            write(86) nsubdomains,nsubdomains_x,nsubdomains_y
          else
            open(86,file=trim(filename),status='unknown',form='unformatted',&
                 position='append')
          endif

          ! write
          write(86) nprev,ncurr,uprev,vprev,wprev,tprev,micro_fieldprev,&
                    sgs_fieldprev,dtnprev
          close(86)
        end if
      end do
    end if

    ! print the prorgess
    if (masterproc) then
      write(*,*) 'Saved openybc restart file, nstep=',nstep
    end if
    call task_barrier()

  end subroutine write_restart_openbc

!============================================================================

  subroutine read_restart_openbc()

    !=========================================|
    ! reasd the restart file of the openybc   |
    ! variables when using nrestart.ne.0      |
    !=========================================|

    implicit none

    character(len=4) :: rankchar
    integer :: irank,ii
    character(len=256) :: filename

    ! print the progress 
    if (masterproc) then
      write(*,*) 'Reading openybc restart file ...'
    end if
    
    if (restart_sep) then
      write(*,*) 'Currently not supporting restart_sep'
      return
    else
      write(rankchar,'(i4)') nsubdomains
      filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                 rankchar(5-lenstr(rankchar):4)//'_openybc_restart.bin'
      open(86,file=trim(filename),status='unknown',form='unformatted')
      
      ! loop out the rank to read the values
      do irank=0,nsubdomains-1
        call task_barrier()
        if (irank.eq.rank) then
          read(86)
          ! skip records
          do ii=0,irank-1
            read(86)
          end do
          ! read the data
          read(86) nprev,ncurr,uprev,vprev,wprev,tprev,micro_fieldprev,&
                   sgs_fieldprev,dtnprev
          close(86)
        end if
      end do
    end if

  end subroutine read_restart_openbc

end module open_boundary_y