module damping_lateral

  !------------------------------------------------------------------------|
  ! sponge layer on lateral boundary (y boundary)                          |
  ! assuming a rayleigh type damping term added to the prognostic equation |
  ! written by Leo Chow (20230718)                                         |
  ! email: leonardotnchow@cuhk.edu.hk                                      |
  !------------------------------------------------------------------------|

  use grid, only: nx, ny, nzm, dt, nyp1

  implicit none

  ! rayleigh damping coefficient
  real :: Ry(nx,ny,nzm)

  ! reference state
  ! v, w, qv0, qn0, qp0 at the y BL is 0
  ! others: u, t, tends to the initial value
  real :: u0initb(nzm), t0initb(nzm)
  real :: u0initt(nzm), t0initt(nzm)
  ! 2D arrays that is damped to: y, z (symmetric along x)
  real :: u0init(ny, nzm), t0init(ny,nzm)

  ! number of grid at y dimension (bottom and top) to apply the coef
  integer :: yd = 200

  ! number of grid at y dimension (bottom and top) to 
  ! compute the average of the reference state
  integer :: yave = 5

  ! damping constant
  real :: alpha = 1. / 200. ! at 20s, 
!  real :: alpha = 0.1 / dt

  ! global average variables
  ! logical for initializing the smoothstep function
  logical :: smy_init=.false.

  CONTAINS

!============================================================================
!!! initialize the rayleigh damping function
  subroutine damping_lateral_init()

    !-------------------------------------------------------------------|
    ! initialize the rayleigh damping coefficient                       |
    ! for each grid in MPI configuration                                |
    ! reference:                                                        |
    ! https://www2.mmm.ucar.edu/wrf/users/docs/technote/v3_technote.pdf |
    ! Section 4.4.2                                                     |
    !-------------------------------------------------------------------|

    use vars, only: u, t
    use grid, only: rank, dy, ny_gl, masterproc, nsubdomains
    use domain, only: nsubdomains_x,nsubdomains_y

    implicit none
    
    include 'mpif.h'

    real(8) :: pii
    real :: tmp
    real(8) :: coef, coef1, buffer(nzm,4), buffer1(nzm,4) ! MPI stuff
    integer :: ierr ! debug use of MPI

    integer :: i, j, k
    integer :: it, jt

    ! pi:
    pii = 4.d0 * atan(1.d0)

    ! print the progress
    if (masterproc) write(*,*) 'Initializing lateral rayleigh damping at yBL...'

    ! initialize the rayleigh coefficient 
    Ry(:,:,:) = 0.

    ! compute the rayleigh damping coefficient
    ! y < yd, ry = alpha * (sin(pii/2 * (1-y*dy/yd*dy))**2
    ! loop the y to compute the global y index
    do j=1,ny
      call task_rank_to_index(rank,it,jt)
      ! bottom y boundary
      if (jt+j.lt.yd+1) then
        tmp = sin(pii/2. * (1. - float(j+jt-1)/float(yd)))
        Ry(:,j,:) = alpha * tmp * tmp
      elseif (jt+j.gt.(ny_gl-yd+1-1)) then
        tmp = sin(pii/2. * (1. + (float(j+jt-ny_gl)/float(yd))))
        Ry(:,j,:) = alpha * tmp * tmp
      end if
    end do

  !------------------------------------------
    !!! compute the global reference value of u and t at the y-boundaries
    ! (initialize, before reaching to the first point of recomputing zonal average)
    
    call zonal_ave()

  !-------------------------------------------

    ! debug
!    do j=1,ny
!      call task_rank_to_index(rank,it,jt)
!      if (rank.eq.0.and.jt+j.lt.yd+10) then
!        write(*,*) 'Rayleigh lateral damping coef and ref u, t: ', rank, jt+j, Ry(1,j,1), & 
!                                                                   u0init(j,1), t0init(j,1)
!      elseif (rank.eq.20.and.jt+j.gt.ny_gl-yd-10) then
!        write(*,*) 'Rayleigh lateral damping coef and ref u, t: ', rank, jt+j, Ry(11,j,11), &
!                                                                   u0init(j,11), t0init(j,11)
!      end if
!    end do
!    if (masterproc) then
!      write(*,*) '          k         u0_bottom         t0_bottom          u0_top          t0_top'
!      do k=nzm,1,-1
!        write(*,*) k, u0initb(k), t0initb(k), u0initt(k), t0initt(k)
!      end do
!    end if

  end subroutine damping_lateral_init

!============================================================================

  subroutine average_to_yb()

    !---------------------------------------------------------|
    ! do a smooth transition of prognostic variables (u, t)   |
    ! from the values at the some y_min to the y boundary     |
    ! by applying a smoothstep function which                 |
    ! smoothly transform the values from 0 to 1               |
    ! that smoothly add values solved by intergation          |
    ! to an zonal averaged values of (u, t)                   |
    ! by aF(y) + b\bar{F(y)}                                  |
    ! reference:                                              |
    ! 1. https://en.wikipedia.org/wiki/Smoothstep             |
    ! 2. https://stackoverflow.com/questions/45165452/        |
    !    how-to-implement-a-smooth-clamp-function-in-python   |
    ! with some timescale of averaging                        |
    !---------------------------------------------------------|

    use vars, only: u,t
    use grid, only: rank,nstep,masterproc,icycle
    use domain, only: nsubdomains_x, nsubdomains_y

    implicit none

    integer :: it, jt ! compute the global index
    integer :: i,j,k
    integer :: ave_time = 21600 ! unit: second
    
    ! --------------
    ! number of grid to replace by zonal averaged values adjecent to y boundary
    integer :: num_y = 5 !  
    ! number of grid near the num_y to smoothly transit the values to zonal ave. 
    integer :: num_decay_y = 125 
    ! array to store the smoothstep function
    real(8) :: smy(nx,ny,nzm)

    !!! compute the smoothstep function
    if (.not.smy_init) then
      call smoothstep(smy,num_y,num_decay_y)
      smy_init=.true.
    end if

    !!!! if reached multiples of average time scale (and skip subcycle)
    if (mod(int(nstep*dt),ave_time).eq.0.and.icycle.eq.1) then
      
      ! update the zonal average at the y boundary
      call zonal_ave()

!      !!! print the progress
!      if (masterproc) then
!        write(*,*) 'doing smoothstep transition to u and t...'
!      end if
!  
!      !!! replace the values near the y boundary by the zonal average
!      ! southmost boundary
!      if (rank.lt.nsubdomains_x) then
!        do j=1,num_y
!          do i=1,nx
!            do k=1,nzm
!              u(i,j,k) = u0init(j,k)
!              t(i,j,k) = t0init(j,k)
!            end do
!          end do
!        end do
!      end if
!      ! northmost boundary
!      if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
!        do j=ny-num_y+1,ny
!          do i=1,nx
!            do k=1,nzm
!              u(i,j,k) = u0init(j,k)
!              t(i,j,k) = t0init(j,k)
!            end do
!          end do
!        end do
!      end if
!  
!      !!! replace the buffer zone by sum of zonal average and original values
!      ! with rating as the smoothstep function
!      ! aF(y) + b\bar{F(y)}
!      call task_rank_to_index(rank,it,jt)
!      do i=1,nx
!        do j=1,ny
!          do k=1,nzm
!            u(i,j,k) = (1.-smy(i,j,k))*u(i,j,k) + (smy(i,j,k))*u0init(j,k)
!            t(i,j,k) = (1.-smy(i,j,k))*t(i,j,k) + (smy(i,j,k))*t0init(j,k)
!          end do
!        end do
!      end do
    end if ! if (mod(int(nstep*dt),ave_time).eq.0.and.icycle.eq.1)

  end subroutine average_to_yb

!============================================================================
  subroutine smoothstep(smy,num_yi,num_decay_yi)

    ! compute the smoothstep function

    use grid, only: rank, nx, ny, nzm,masterproc
    use domain, only: ny_gl

    implicit none

    integer, intent(in) :: num_yi
    integer, intent(in) :: num_decay_yi
    real(8), intent(inout) :: smy(nx,ny,nzm) ! input: empty
                                             ! output: smoothstep values for each grid points
    real :: ytemp 
    integer :: it, jt
    integer :: i,j,k

    call task_rank_to_index(rank,it,jt)

    ! initialize
    smy(:,:,:) = 0.
    ! compute the smoothstemp function
    do j=1,ny
      do i=1,nx
        do k=1,nzm
          ! south y boundaries (from 1 to 0)
          if ((jt+j).gt.(num_yi-1).and.(jt+j).lt.(num_decay_yi+num_yi+1)) then
            ytemp = float(num_yi + num_decay_yi - (j + jt)) * &
                    1. / float(num_decay_yi + num_yi - num_yi)
            smy(i,j,k) = 35.*ytemp*ytemp*ytemp*ytemp - &
                         84.*ytemp*ytemp*ytemp*ytemp*ytemp + &
                         70.*ytemp*ytemp*ytemp*ytemp*ytemp*ytemp - &
                         20.*ytemp*ytemp*ytemp*ytemp*ytemp*ytemp*ytemp
            ! correct truncation error
            if (smy(i,j,k).gt.1.) then
              smy(i,j,k) = 1.
            else if (smy(i,j,k).lt.0.) then
              smy(i,j,k) = 0.
            end if
          ! northern y boudaries (from 0 to 1)
          else if ((j+jt).gt.(ny_gl-num_yi-num_decay_yi-1).and.&
                   (j+jt).lt.(ny_gl-num_yi+1)) then
            ytemp = float(j+jt-(ny_gl-num_yi-num_decay_yi)) * &
                    1./float((ny_gl-num_yi) - (ny_gl-num_yi-num_decay_yi))
            smy(i,j,k) = 35.*ytemp*ytemp*ytemp*ytemp - &
                         84.*ytemp*ytemp*ytemp*ytemp*ytemp + &
                         70.*ytemp*ytemp*ytemp*ytemp*ytemp*ytemp - &
                         20.*ytemp*ytemp*ytemp*ytemp*ytemp*ytemp*ytemp
            ! correct truncation error
            if (smy(i,j,k).gt.1.) then
              smy(i,j,k) = 1.
            else if (smy(i,j,k).lt.0.) then
              smy(i,j,k) = 0.
            end if
          ! zonal averaged area
          else if ((j+jt).lt.num_yi.or.(j+jt).gt.(ny_gl-num_yi)) then
            smy(i,j,k) = 1.
          end if
        end do
      end do
    end do

  end subroutine smoothstep

!============================================================================

  subroutine zonal_ave()

    !------------------------------------|
    ! compute the zonal average          |
    ! of prognostic variables u, t       |
    ! from yave to the y boundaries      |
    !------------------------------------|

    use vars, only: u, t
    use grid, only: nx,ny,nzm,rank,masterproc
    use domain, only: nsubdomains_x, nsubdomains_y, ny_gl

    implicit none

    include 'mpif.h'

    real(8) :: coef, coef1, buffer(nzm,4), buffer1(nzm,4) ! MPI stuff
    integer :: ierr ! debug use of MPI

    integer :: i, j, k
    integer :: it, jt

    !!! print the progress
    if (masterproc) then
      write(*,*) 'updating nudging values of u and t...'
    end if

    !!! compute the global reference value of u and t at the y-boundaries
    ! local number of grid points
    coef = 1./float(yave*nx)
    ! initialize
    u0initb = 0.
    t0initb = 0.
    u0initt = 0.
    t0initt = 0.
    !!! south y boundary:
    if (rank.lt.nsubdomains_x) then
      ! store the initial u0, t0 (local)
      do k=1,nzm
        do i=1,nx
          do j=1,yave
            u0initb(k) = u0initb(k) + u(i,j,k)
            t0initb(k) = t0initb(k) + t(i,j,k)
          end do
        end do
        ! do the local averave
        u0initb(k) = u0initb(k) * coef
        t0initb(k) = t0initb(k) * coef
      end do
    end if
    !!! north y boundary
    if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
      ! store the initial u0, t0 (local)
      do k=1,nzm
        do i=1,nx
          do j=ny-yave+1,ny
            u0initt(k) = u0initt(k) + u(i,j,k)
            t0initt(k) = t0initt(k) + t(i,j,k)
          end do
        end do
        ! do the local averave
        u0initt(k) = u0initt(k) * coef
        t0initt(k) = t0initt(k) * coef
      end do
    end if

    ! store the data to buffer
    ! becuase the initial value of u0 and t0 is 0, so will not affect the 
    ! buffer value, and can be divided directly by nsubdomains_x
    do k=1,nzm
      buffer(k,1) = u0initb(k)
      buffer(k,2) = t0initb(k)
      buffer(k,3) = u0initt(k)
      buffer(k,4) = t0initt(k)
    end do
    ! global number of processors
    coef1 = 1./float(nsubdomains_x)
    ! compute the global sum
    call task_barrier()
    call MPI_ALLREDUCE(buffer, buffer1, nzm*4, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    ! store it back to local arrays
    do k=1,nzm
      u0initb(k) = buffer1(k,1) * coef1
      t0initb(k) = buffer1(k,2) * coef1
      u0initt(k) = buffer1(k,3) * coef1
      t0initt(k) = buffer1(k,4) * coef1
    end do
    
    !!! store the data to 2D arrays that will be damped to
    ! initialize
    u0init(:,:) = 0.
    t0init(:,:) = 0.
    call task_rank_to_index(rank,it,jt)
    do k=1,nzm
      do j=1,ny
        ! southmost
        if (j+jt.lt.yd+1) then
          u0init(j,k) = u0initb(k)
          t0init(j,k) = t0initb(k)
        ! northmost
        else if (j+jt.gt.(ny_gl-yd+1-1)) then
          u0init(j,k) = u0initt(k)
          t0init(j,k) = t0initt(k)
        end if
      end do
    end do

  end subroutine zonal_ave

!============================================================================

  ! add the rayleigh damping term to the prognostic equations
  ! i.e. dudt, dvdt, dwdt, t
  subroutine damping_lateral_proc()

    use vars, only: dudt, dvdt, dwdt, t, u, v, w, &
                    qv, qcl, qpl, qci, qpi, &
                    dvdtnorthmost, dvdtsouthmost !, u0, v0, t0
    use grid, only: dtn, na, rank
    use domain, only: nsubdomains_x, nsubdomains_y
    use params, only: doopeny ! leo: ybl

    implicit none

    integer :: i,j,k

    ! loop out the index to add the rayleigh damping term
    do k=1,nzm
      do j=1,ny
        do i=1,nx

          dudt(i,j,k,na) = dudt(i,j,k,na) - Ry(i,j,k) * (u(i,j,k) - u0init(j,k))
          dvdt(i,j,k,na) = dvdt(i,j,k,na) - Ry(i,j,k) * v(i,j,k) ! at yBL, v -> 0
          dwdt(i,j,k,na) = dwdt(i,j,k,na) - Ry(i,j,k) * w(i,j,k) ! the reference state 
                                                                 ! of w is zero
!          t(i,j,k) = t(i,j,k) - dtn * Ry(i,j,k) * (t(i,j,k) - t0init(j,k))
!          ! microphysics variables
!          qv(i,j,k) = qv(i,j,k) - dtn * Ry(i,j,k) * qv(i,j,k)
!          qcl(i,j,k) = qcl(i,j,k) - dtn * Ry(i,j,k) * qcl(i,j,k)
!          qpl(i,j,k) = qpl(i,j,k) - dtn * Ry(i,j,k) * qpl(i,j,k)
!          qci(i,j,k) = qci(i,j,k) - dtn * Ry(i,j,k) * qci(i,j,k)
!          qpi(i,j,k) = qpi(i,j,k) - dtn * Ry(i,j,k) * qpi(i,j,k)

        end do
      end do
    end do
  
    ! leo: ybl
    ! seperately treat the nyp1 value for v
    if (doopeny) then
      ! MPI, northmost grid cell only
      if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
        do k=1,nzm
          do i=1,nx
            ! just naively assume the damping coefficient at nyp1 = ny
            dvdt(i,nyp1,k,na) = dvdt(i,nyp1,k,na) - Ry(i,ny,k) * v(i,nyp1,k)
          end do
        end do
      end if
      ! backup the value of dvdt at the northmost and southmost grid cell
      if (rank.lt.nsubdomains_x) then
        dvdtsouthmost(:,1,:) = dvdt(:,1,:,na)
      else if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
        dvdtnorthmost(:,1,:) = dvdt(:,ny+1,:,na)
      end if
    end if
  
  end subroutine damping_lateral_proc

!============================================================================

  !!! damp the scalar variables (t and microphysics variables)
  ! after the adams.f90
  ! by first-order forward-in-time differencing
  subroutine damping_lateral_scalar_proc()

    use vars, only: t
    use microphysics, only: micro_field,nmicro_fields
    use grid, only: dtn, rank

    implicit none

    integer :: i,j,k,nm

    !!! loop out the index to do the rayleigh damping
    ! forward-in-time differencing
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          ! scalar variables
          t(i,j,k) = t(i,j,k) - dtn * Ry(i,j,k) * (t(i,j,k) - t0init(j,k))
          ! microphysics variables in MICRO_M2005
          do nm=1,nmicro_fields
            micro_field(i,j,k,nm) = micro_field(i,j,k,nm) - & 
                                    dtn * Ry(i,j,k) * micro_field(i,j,k,nm)
                                    ! reference at the BL is 0 for 
                                    ! microphysics variables
          end do
        end do
      end do
    end do

    end subroutine damping_lateral_scalar_proc

end module damping_lateral

