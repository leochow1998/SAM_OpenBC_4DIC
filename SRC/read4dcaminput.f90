!--------------------------------------------------------------------
! file: read4dcaminput.f
! Author: Leo Chow (leonardotnchow@cuhk.edu.hk)
! aim:
!     to read in data from CAM history file to initialize
!     the 3D initial condition of SAM
!     with 2 stage:
!     1. initialize IC except the microphysics variables
!     2. initialize IC for the microphysics scheme (M2005)
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! modified for the use of mpi (20230614) (Leo Chow)
! 1. from the masterproc (rank=0), read all the data to &
!    local variables that is only defined in the masterproc
! 2. compute the index from each mpi rank, locate the i and j &
!    read and store the variables to each rank for each cpu
!    by looping out the rank=0 together with each rank
!--------------------------------------------------------------------

subroutine read4dcaminput ()
 
  ! use vars
  use vars
  use params
  use microphysics
  ! variables from vars
!  use vars, only: u, v, w, t, &
!                  tabs, &
!                  fluxbu, fluxbv, fluxbt, fluxbq, &
!                  fluxtu, fluxtv, fluxtt, fluxtq, &
!                  precsfc, rhow
  use grid, only: case, rank, masterproc, nsubdomains, &
                  docam4diopdata, iopfile4d
    
  implicit none

  include 'netcdf.inc'
  include 'mpif.h'

!----------------------------- input --------------------------------
    
  integer error_code        ! return netcdf errors

!------------------------- local variables --------------------------

  ! 3D variables from CAM
  ! prognostic variables
  real, allocatable :: u_in(:,:,:) ! zonal wind
  real, allocatable :: v_in(:,:,:) ! meridional wind
  real, allocatable :: w_in(:,:,:) ! vertical wind in height coordinate
  real, allocatable :: s_in(:,:,:) ! dry static energy (T + g / cp * z)

  ! diagnostic variables
  real, allocatable :: t_in(:,:,:) ! absolute temperature
  ! delected, because it is diagnostic, computed every timestep
  ! by solving the poisson equation on u, v, w, and hydrometers concentration
!-  real, allocatable :: ppert_in ! perturbation pressure
!-  real, allocatable :: rhow_in(:,:,:) ! density at vertical velocity 
  !                                     level surface (to compute the fluxes only)
  
  ! hydrometer and concentration variables
  real, allocatable :: q_in(:,:,:) ! water vapor mixing ratio (kg/kg)
  real, allocatable :: qcl_in(:,:,:) ! cloud liquid water mixing ratio (kg/kg)
  real, allocatable :: ncl_in(:,:,:) ! cloud water number mixing ratio (#/kg)
  real, allocatable :: qr_in(:,:,:) ! rain mixing ratio (kg/kg)
  real, allocatable :: nr_in(:,:,:) ! rain number mixing ratio (#/kg)
  real, allocatable :: qci_in(:,:,:) ! ice mixing ratio (kg/kg)
  real, allocatable :: nci_in(:,:,:) ! ice number mixing ratio (#/kg)
  real, allocatable :: qs_in(:,:,:) ! snow mixing ratio (kg/kg)
  real, allocatable :: ns_in(:,:,:) ! snow number mixing ratio (#/kg)
  ! there is no graupel output in the CAM hist and rest file
  ! yet, the model compute the graupel concentration and mixing ratio

  ! 2d variables from CAM
  ! convective precip. rate, large-scale precip. rate (liq + ice)
  ! delected, because it is computed from the start of every day
!-  real, allocatable :: precc_in(:,:), precl_in(:,:) 
  ! surface sensible heat flux, surface latent heat flux (W/m2)
  ! fluxbt = shf_in / cp / rhows, fluxbq = lhf_in / lcond / rhows
  ! diagnostic, can be ignored (20230627)
!-  real, allocatable :: shf_in(:,:), lhf_in(:,:)
  ! surface zonal and meridional stress
!-  real, allocatable :: taux_in(:,:), tauy_in(:,:)
  ! surface pressre
  

  ! intermediate array
  ! precipitation
  ! delected, because it is computed from the start of each day
!-  real, allocatable :: prec_in(:,:)
  ! fluxes
!-  real, allocatable :: fluxbt_in(:,:), fluxbq_in(:,:), &
!-                       fluxbu_in(:,:), fluxbv_in(:,:)
  
  ! dimensions
  integer :: nx_in, ny_in, nz_in
  real, allocatable :: x_in(:), y_in(:), zlev_in(:)

  ! mpich temperary variables for receive
  real :: u_rec(nx,ny,nzm), v_rec(nx,ny,nzm), s_rec(nx,ny,nzm), &
          t_rec(nx,ny,nzm), ppert_rec(nx,ny,nzm), &
!-          fluxbu_rec(nx,ny), fluxbv_rec(nx,ny), fluxbt_rec(nx,ny), &
!-          fluxbq_rec(nx,ny), &
!-          prec_rec(nx,ny), &
          w_rec(nx,ny,nzm+1), q_rec(nx,ny,nzm), &
          qcl_rec(nx,ny,nzm), ncl_rec(nx,ny,nzm), &
          qr_rec(nx,ny,nzm), nr_rec(nx,ny,nzm), &
          qci_rec(nx,ny,nzm), nci_rec(nx,ny,nzm), &
          qs_rec(nx,ny,nzm), ns_rec(nx,ny,nzm)

  ! variables for satuation adjustment
  real :: tmpqv(nzm), tmpqcl(nzm), tmptabs(nzm), tmppres(nzm) ! this pressure is profile
  real :: qv_tmp(nx,ny,nzm), qcl_tmp(nx,ny,nzm), t_tmp(nx,ny,nzm)

  ! mpich related variables
  integer :: it, jt
  integer :: irank_send, n_send, ierr
  real(8) :: coef, coef1, &
             buffer(nzm,7), buffer1(nzm,7) ! store and compute the global vertical profiles

  ! others
  integer :: i,j,k
  integer :: ncid, status ! for netcdf file opening
  character(len=120) :: iopfilepath
  logical :: use_nf_real = .true. ! use 4 byte real, but not double for the nf variables
  integer :: tot_num_var = 20 ! actually using only 17
  ! fill_value
  real, parameter :: missing_value = -99999.

  ! testing diagnostic variables
!  real :: tabs_diag(nzm), qcl1_diag(nzm), qcl_diag(nzm), qci_diag(nzm)

  ! temporary
!  character(len=120) :: iopfile4d = 'Qobs10_init_new.cam.h0.0001-03-29-00000_apc.nc'
!  character(len=120) :: iopfileloc = '/home/tnchow/Desktop/'

!---------------------- read CAM history file -----------------------

  ! the iopfile4d should be setted in the namelist (leo, 20230331)
!  iopfilepath = trim(iopfileloc)//'/'//trim(iopfile4d)
  iopfilepath = './'//trim(case)//'/'//trim(iopfile4d)
  if (masterproc) then
    write(*,*) 'Opening ', iopfilepath
    status = NF_OPEN(iopfilepath, NF_NOWRITE, ncid)
    if (status.ne.NF_NOERR) then
        if (masterproc) then
            write (6,*) & 
            'ERROR(read4dcaminput.f90):Cant open iop dataset: ', iopfilepath
            call task_abort()
        endif 
    endif 
  endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read and allocate data
  ! masterproc, read and allocate the data to internal arrays
  if (masterproc) then
    write(*,*) 'Reading data from ', iopfilepath, ' ...'
!-------------------- read vertical level data ----------------------

    ! get the dim of the height coordinate
    call get_netcdf_dimlength(ncid, 'zlev', nz_in, status, .true.)

    ! allocate the space to the zlev variables
!    ALLOCATE(zlev_in(nz_in),STAT=status) 
!    if (status.ne.0) then
!      write(6,*) 'Could not allocate zlev in read4dcaminput'
!      call task_abort()
!    endif
  
!    ! read the vertical level in height coordinate
!    call get_netcdf_var1d_real(ncid, 'zlev', zlev_in, use_nf_real, status, .true.)

!--------------- read horizontal and vertical coords ----------------

    ! get the dim of the horizontal and the vertical coords
    call get_netcdf_dimlength(ncid, 'x', nx_in, status, .true.)
    call get_netcdf_dimlength(ncid, 'y', ny_in, status, .true.)
  
    ! allocate the space to the horizontal and the vertical coordinate
!    ALLOCATE(x_in(nx_in), y_in(ny_in), STAT=status)
!    if (status.ne.0) then
!      write(6,*) 'Could not allocate x/y in read4dcaminput'
!      call task_abort()
!    endif
  
    ! get x and y
!    call get_netcdf_var1d_real(ncid, 'x', x_in, use_nf_real, status, .true.)
!    call get_netcdf_var1d_real(ncid, 'y', y_in, use_nf_real, status, .true.)

!========================= 3d variables =============================
!---------------- allocate space to the 3d variables ----------------
  
    ALLOCATE(u_in(nx_in, ny_in, nz_in), v_in(nx_in, ny_in, nz_in), &
             w_in(nx_in, ny_in, nz_in+1), t_in(nx_in, ny_in, nz_in), &
             s_in(nx_in, ny_in, nz_in), &
!-             ppert_in(nx_in, ny_in, nz_in), &
!-             rhow_in(nx_in, ny_in, nz_in+1), &
             q_in(nx_in, ny_in, nz_in), qcl_in(nx_in, ny_in, nz_in), &
             ncl_in(nx_in, ny_in, nz_in), qr_in(nx_in, ny_in, nz_in), &
             nr_in(nx_in, ny_in, nz_in), qci_in(nx_in, ny_in, nz_in), &
             nci_in(nx_in, ny_in, nz_in), qs_in(nx_in, ny_in, nz_in), &
             ns_in(nx_in, ny_in, nz_in), &
             STAT=status)
    if (status.ne.0) then
      write(6,*) 'Could not allocate 3D variables in read4dcaminput'
      call task_abort()
    endif

!------------------------ read 3d variables -------------------------

    ! zonal wind 
    u_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'U', u_in, use_nf_real, status, .true.)
  
    ! meridional wind
    v_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'V', v_in, use_nf_real, status, .true.)
  
    ! vertical wind (z = nz + 1, on vertical wind level)
    w_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'W', w_in, use_nf_real, status, .true.)
  
    ! absolute temperature
    t_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'T', t_in, use_nf_real, status, .true.)
  
    ! dry static energy (Tabs + g / cp * z)
    s_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'S', s_in, use_nf_real, status, .true.)
    
    ! perturbation pressure
    ! computed by solving the poisson equation in the model
    ! for every timestep (with u, v, w, du, dv, dw, and hydrometer variable)
!-    ppert_in(:,:,:) = missing_value
!-    call get_netcdf_var3d_real(ncid, 'p_prime', ppert_in, use_nf_real, status, .true.)
  
    ! density on vertical velocity surface (z = nz + 1)
!-    rhow_in(:,:,:) = missing_value
!-    call get_netcdf_var3d_real(ncid, 'rhow', rhow_in, use_nf_real, status, .true.)
  
    !--- microphysics (double moment microphysics schemes)
    ! water vapor mixing ratio (kg/kg)
    q_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid,'Q', q_in, use_nf_real, status, .true.)
    
    ! cloud liquid water mixing ratio (kg/kg)
    qcl_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'CLDLIQ', qcl_in, use_nf_real, status, .true.)
  
    ! cloud water number mixing ratio (#/kg)
    ncl_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'NUMLIQ', ncl_in, use_nf_real, status, .true.)
  
    ! rain mixing ratio (kg/kg)
    qr_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'AQRAIN', qr_in, use_nf_real, status, .true.)
  
    ! rain number mixing ratio (#/kg)
    nr_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'ANRAIN', nr_in, use_nf_real, status, .true.)
  
    ! ice mixing ratio (kg/kg)
    qci_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'CLDICE', qci_in, use_nf_real, status, .true.)
  
    ! ice number mixing ratio (#/kg)
    nci_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'NUMICE', nci_in, use_nf_real, status, .true.)
  
    ! snow mixing ratio (kg/kg)
    qs_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'AQSNOW', qs_in, use_nf_real, status, .true.)
  
    ! snow number mixing ratio (#/kg)
    ns_in(:,:,:) = missing_value
    call get_netcdf_var3d_real(ncid, 'ANSNOW', ns_in, use_nf_real, status, .true.)

!========================= 2d variables =============================
!-------------- allocate space to the 2d variables ------------------
!-    ALLOCATE(shf_in(nx_in,ny_in),lhf_in(nx_in,ny_in), & 
!-             precc_in(nx_in,ny_in),precl_in(nx_in,ny_in), &
!-             taux_in(nx_in,ny_in), tauy_in(nx_in,ny_in), &
!-             STAT=STATUS)
    
!-    if (status.ne.0) then
!-      write(6,*) 'Could not allocate 2d variables in read4dcaminput'
!-      call task_abort()
!-    end if

!---------------------- read 2d variables ---------------------------

    ! surface heat flux unit: W/m2
!-    shf_in(:,:) = missing_value
!-    call get_netcdf_var2d_real(ncid, 'SHFLX', shf_in, use_nf_real, status, .true.)
  
    ! latent heat flux
!-    lhf_in(:,:) = missing_value
!-    call get_netcdf_var2d_real(ncid, 'LHFLX', lhf_in, use_nf_real, status, .true.)
  
    ! convective precipitation ! now the unit is ms-1
!-    precc_in(:,:) = missing_value
!-    call get_netcdf_var2d_real(ncid, 'PRECC', precc_in, use_nf_real, status, .true.)
  
    ! large-scale precipitation
!-    precl_in(:,:) = missing_value
!-    call get_netcdf_var2d_real(ncid, 'PRECL', precl_in, use_nf_real, status, .true.)
  
    ! surface zonal momentum stress
    ! fluxu0 = taux / rhow_1
!-    taux_in(:,:) = missing_value
!-    call get_netcdf_var2d_real(ncid, 'TAUX', taux_in, use_nf_real, status, .true.)
  
    ! surface meridional momentum stress
!-    tauy_in(:,:) = missing_value
!-    call get_netcdf_var2d_real(ncid, 'TAUY', tauy_in, use_nf_real, status, .true.)

!===================================================================
!            convert some of the read data to SAM's units           
!===================================================================

    ! print the progress
!-    write(*,*) 'Computing intermediate arrays ...'

!-    ! allocate the intermediate arrays
!-    ALLOCATE(fluxbt_in(nx_in, ny_in), &
!-             fluxbq_in(nx_in, ny_in), &
!-             fluxbu_in(nx_in, ny_in), &
!-             fluxbv_in(nx_in, ny_in), &
!-             prec_in(nx_in, ny_in), &
!-             STAT=STATUS)
  
!-    if (status.ne.0) then
!-      write(6,*) 'Could not allocate 2d temp variables in read4dcaminput'
!-      call task_abort()
!-    end if 
  
    ! precipitation
!-    prec_in = precc_in + precl_in
  
    ! sensible heat and latent heat fluxes
!-    fluxbt_in = shf_in / rhow_in(:,:,1) / cp
!-    fluxbq_in = lhf_in / rhow_in(:,:,1) / lcond
  
    ! surface momentum suress
!-    fluxbu_in = taux_in / rhow_in(:,:,1)
!-    fluxbv_in = tauy_in / rhow_in(:,:,1)

  end if ! if masterproc, read and allocate data to local arrays

  call task_barrier() ! let all rank (processors) wait here until 
                      ! the reading completed

  ! diagnostic testing
!  if (masterproc) then
!    do k=1,nzm
!      qcl_diag(k) = sum(dble(qcl_in(1:nx,1:ny,k)))/float(nx*ny)
!      write(6,'(i4,1x,f8.4)') k,qcl_diag(k)*1e3
!    enddo
!  endif

!=============== store the read data to SAM's array ================
  
  ! store the data from X_in to X_rec
  ! for masterproc
  if (masterproc) then

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          u_rec(i,j,k) = u_in(i,j,k)
          v_rec(i,j,k) = v_in(i,j,k)
          s_rec(i,j,k) = s_in(i,j,k)
          t_rec(i,j,k) = t_in(i,j,k)
!-          ppert_rec(i,j,k) = ppert_in(i,j,k)
        end do
      end do
    end do

    ! diagnostic, no need to compute!*
!-    do j=1,ny
!-      do i=1,nx
!-        fluxbu_rec(i,j) = fluxbu_in(i,j)
!-        fluxbv_rec(i,j) = fluxbv_in(i,j)
!-        fluxbt_rec(i,j) = fluxbt_in(i,j)
!-        fluxbq_rec(i,j) = fluxbq_in(i,j)
!-        prec_rec(i,j) = prec_in(i,j)
!-      end do
!-    end do

    do k=1,nzm+1
      do j=1,ny
        do i=1,nx
          w_rec(i,j,k) = w_in(i,j,k)
        end do
      end do
    end do

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          q_rec(i,j,k) = q_in(i,j,k)
          qcl_rec(i,j,k) = qcl_in(i,j,k)
          ncl_rec(i,j,k) = ncl_in(i,j,k)
          qr_rec(i,j,k) = qr_in(i,j,k)
          nr_rec(i,j,k) = nr_in(i,j,k)
          qci_rec(i,j,k) = qci_in(i,j,k)
          nci_rec(i,j,k) = nci_in(i,j,k)
          qs_rec(i,j,k) = qs_in(i,j,k)
          ns_rec(i,j,k) = ns_in(i,j,k)
        end do
      end do
    end do
  
  end if ! if masterproc
  call task_barrier()

  ! sending data to each subprocesses
  do irank_send=1,nsubdomains-1
    if (masterproc) then

      ! compute the tag for the first variables -1 to send to rank = irank_send
      ! so each variable in each subprocesses with have a unique tag
      n_send = (irank_send - 1) * tot_num_var
      ! compute the start i and j for each subprocesses globally
      call task_rank_to_index(irank_send,it,jt)
      write(*,*) 'Sending data to subrank: ',irank_send, it,jt

      ! u_in
!*      write(*,*) 'Start to send u_in to rank: ', irank_send
      call MPI_SEND(u_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send, (n_send+1),MPI_COMM_WORLD,ierr)
      
      ! v_in
!*      write(*,*) 'Start to send v_in to rank: ', irank_send
      call MPI_SEND(v_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+2),MPI_COMM_WORLD,ierr)

      ! s_in
!*      write(*,*) 'Start to send s_in to rank: ', irank_send
      call MPI_SEND(s_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+3),MPI_COMM_WORLD,ierr)

      ! t_in
!*      write(*,*) 'Start to send t_in to rank: ', irank_send
      call MPI_SEND(t_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+4),MPI_COMM_WORLD,ierr)

      ! ppert_in
!-      write(*,*) 'Start to send ppert_in to rank: ', irank_send
!-      call MPI_SEND(ppert_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
!-                    irank_send,(n_send+5),MPI_COMM_WORLD,ierr)

      ! fluxbu_in
!*      write(*,*) 'Start to send fluxbu_in to rank: ', irank_send
!*      call MPI_SEND(fluxbu_in(it+1:it+nx,jt+1:jt+ny),nx*ny,MPI_REAL,&
!*                    irank_send,(n_send+6),MPI_COMM_WORLD,ierr)

      ! fluxbv_in
!*      write(*,*) 'Start to send fluxbv_in to rank: ', irank_send
!*      call MPI_SEND(fluxbv_in(it+1:it+nx,jt+1:jt+ny),nx*ny,MPI_REAL,&
!*                    irank_send,(n_send+7),MPI_COMM_WORLD,ierr)

      ! fluxbt_in
!*      write(*,*) 'Start to send fluxbt_in to rank: ', irank_send
!*      call MPI_SEND(fluxbt_in(it+1:it+nx,jt+1:jt+ny),nx*ny,MPI_REAL,&
!*                    irank_send,(n_send+8),MPI_COMM_WORLD,ierr)

      ! fluxbq_in
!*      write(*,*) 'Start to send fluxbq_in to rank: ', irank_send
!*      call MPI_SEND(fluxbq_in(it+1:it+nx,jt+1:jt+ny),nx*ny,MPI_REAL,&
!*                    irank_send,(n_send+9),MPI_COMM_WORLD,ierr)

      ! prec_in
!-      write(*,*) 'Start to send prec_in to rank: ', irank_send
!-      call MPI_SEND(prec_in(it+1:it+nx,jt+1:jt+ny),nx*ny,MPI_REAL,&
!-                    irank_send,(n_send+10),MPI_COMM_WORLD,ierr)

      ! variables on vertical velocity vertical coordinate
      ! w_in
!*      write(*,*) 'Start to send w_in to rank: ', irank_send
      call MPI_SEND(w_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*(nzm+1),MPI_REAL,&
                    irank_send,(n_send+11),MPI_COMM_WORLD,ierr)

      ! microphysics variables
      ! q_in
!*      write(*,*) 'Start to send q_in to rank: ', irank_send
      call MPI_SEND(q_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+12),MPI_COMM_WORLD,ierr)

      ! qcl_in
!*      write(*,*) 'Start to send qcl_in to rank: ', irank_send
      call MPI_SEND(qcl_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+13),MPI_COMM_WORLD,ierr)

      ! ncl_in
!*      write(*,*) 'Start to send ncl_in to rank: ', irank_send
      call MPI_SEND(ncl_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+14),MPI_COMM_WORLD,ierr)

      ! qr_in
!*      write(*,*) 'Start to send qr_in to rank: ', irank_send
      call MPI_SEND(qr_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+15),MPI_COMM_WORLD,ierr)

      ! nr_in
!*      write(*,*) 'Start to send nr_in to rank: ', irank_send
      call MPI_SEND(nr_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+16),MPI_COMM_WORLD,ierr)

      ! qci_in
!*      write(*,*) 'Start to send qci_in to rank: ', irank_send
      call MPI_SEND(qci_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+17),MPI_COMM_WORLD,ierr)

      ! nci_in
!*      write(*,*) 'Start to send nci_in to rank: ', irank_send
      call MPI_SEND(nci_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+18),MPI_COMM_WORLD,ierr)

      ! qs_in
!*      write(*,*) 'Start to send qs_in to rank: ', irank_send
      call MPI_SEND(qs_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+19),MPI_COMM_WORLD,ierr)

      ! ns_in
!*      write(*,*) 'Start to send ns_in to rank: ', irank_send
      call MPI_SEND(ns_in(it+1:it+nx,jt+1:jt+ny,:),nx*ny*nzm,MPI_REAL,&
                    irank_send,(n_send+20),MPI_COMM_WORLD,ierr)

    elseif(irank_send.eq.rank) then
      
      ! compute the tag corresponding to each variable in each subprocesses
      n_send = (rank - 1) * 20
      
      ! u_in
      call MPI_RECV(u_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+1),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving u_in done, here is ', rank

      ! v_in
      call MPI_RECV(v_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+2),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving v_in done, here is ', rank
      
      ! s_in
      call MPI_RECV(s_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+3),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving s_in done, here is ', rank

      ! t_in
      call MPI_RECV(t_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+4),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving t_in done, here is ', rank

      ! ppert_in
!-      call MPI_RECV(ppert_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
!-                    0,(n_send+5),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!-      write(*,*) 'receiving ppert_in done, here is ', rank

      ! fluxbu_in
!*      call MPI_RECV(fluxbu_rec(:,:),nx*ny,MPI_REAL,&
!*                    0,(n_send+6),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving fluxbu_in done, here is ', rank
      
      ! fluxbv_in
!*      call MPI_RECV(fluxbv_rec(:,:),nx*ny,MPI_REAL,&
!*                    0,(n_send+7),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving fluxbv_in done, here is ', rank

      ! fluxbt_in
!*      call MPI_RECV(fluxbt_rec(:,:),nx*ny,MPI_REAL,&
!*                    0,(n_send+8),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving fluxbt_in done, here is ', rank
      
      ! fluxbq_in
!*      call MPI_RECV(fluxbq_rec(:,:),nx*ny,MPI_REAL,&
!*                    0,(n_send+9),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving fluxbq_in done, here is ', rank

      ! prec_in
!-      call MPI_RECV(prec_rec(:,:),nx*ny,MPI_REAL,&
!-                    0,(n_send+10),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!-      write(*,*) 'receiving prec_in done, here is ', rank

      ! variables on vertical velocity vertical coordinates
      ! w_in
      call MPI_RECV(w_rec(:,:,:),nx*ny*(nzm+1),MPI_REAL,&
                    0,(n_send+11),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving w_in done, here is ', rank

      ! microphysics variables
      ! q_in
      call MPI_RECV(q_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+12),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving q_in done, here is ', rank

      ! qcl_in
      call MPI_RECV(qcl_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+13),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving qcl_in done, here is ', rank

      ! ncl_in
      call MPI_RECV(ncl_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+14),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving ncl_in done, here is ', rank

      ! qr_in
      call MPI_RECV(qr_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+15),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving qr_in done, here is ', rank

      ! nr_in
      call MPI_RECV(nr_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+16),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving nr_in done, here is ', rank

      ! qci_in
      call MPI_RECV(qci_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+17),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving qci_in done, here is ', rank

      ! nci_in
      call MPI_RECV(nci_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+18),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving nci_in done, here is ', rank
      
      ! qs_in
      call MPI_RECV(qs_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+19),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving qs_in done, here is ', rank

      ! ns_in
      call MPI_RECV(ns_rec(:,:,:),nx*ny*nzm,MPI_REAL,&
                    0,(n_send+20),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!*      write(*,*) 'receiving ns_in done, here is ', rank
      
!*      if (irank_send.eq.0) then
!*        do k=1,nz
!*          write(*,*) t_rec(11,11,k)
!*        end do
!*      end if
    
    end if
  end do

  call task_barrier()

  if (masterproc) then 
    write(*,*) 'done for sending and receiving'
  end if
  
  ! print the progress
  if (masterproc) write(*,*) 'Storing the data to SAMs array...'

  ! variables on the pressure grid
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        !! store !!
        u(i,j,k) = u_rec(i,j,k)
        v(i,j,k) = v_rec(i,j,k)
        ! here just leave the w
        t(i,j,k) = s_rec(i,j,k)
        tabs(i,j,k) = t_rec(i,j,k)
        ! diagnostic variabeles
!-        p(i,j,k) = ppert_rec(i,j,k)
!-!-        fluxbu(i,j) = fluxbu_rec(i,j) ! diagnostic, no need to compute
!-        fluxbv(i,j) = fluxbv_rec(i,j)
!-        fluxbt(i,j) = fluxbt_rec(i,j)
!-        fluxbq(i,j) = fluxbq_rec(i,j)
!-        precsfc(i,j) = prec_rec(i,j) ! precsfc is computed from the start of each day
!-                                     ! no unit, if insist to use it
      end do 
    end do
  end do

  ! variables on the vertical wind grid
  do k=1,nzm+1
    do j=1,ny
      do i=1,nx
        w(i,j,k) = w_rec(i,j,k)
      end do
    end do
  end do

  !**********************************************************|
  ! do the satuation adjustment to the qv and qc             |
  ! here the reading of qv and qc is seperate                |
  !**********************************************************|

  ! print the progress
  if (masterproc) write(*,*) 'Doing saturation adjustment...'

  ! loop out the i and j 
  do j=1,ny
    do i=1,nx
      !bloss/qt: saturation adjustment to compute cloud liquid water content.
      !          Note: tmpqv holds qv+qcl on input, qv on output.
      !                tmptabs hold T-(L/Cp)*qcl on input, T on output.
      !                tmpqcl hold qcl on output.
      !                tmppres is unchanged on output, should be in Pa.
      tmpqv = q_rec(i,j,:) + qcl_rec(i,j,:)
      tmptabs = t_rec(i,j,:) - fac_cond * qcl_rec(i,j,:)
      tmppres = 100.*pres(:)
      ! call satuation adjustment
      call task_rank_to_index(rank,it,jt)
      call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres,i+it,j+jt)
      ! store the returned array
      do k=1,nzm
        qv_tmp(i,j,k) = tmpqv(k)
        qcl_tmp(i,j,k) = tmpqcl(k)
        t_tmp(i,j,k) = tmptabs(k)
      end do
    end do 
  end do

  ! microphysics variables
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        ! micro_field
!        micro_field(i,j,k,iqv) = q_rec(i,j,k) + qcl_rec(i,j,k)
        micro_field(i,j,k,iqv) = qv_tmp(i,j,k) + qcl_tmp(i,j,k)
        micro_field(i,j,k,incl) = ncl_rec(i,j,k)
        micro_field(i,j,k,iqr) = qr_rec(i,j,k)
        micro_field(i,j,k,inr) = nr_rec(i,j,k)
        micro_field(i,j,k,iqci) = qci_rec(i,j,k)
        micro_field(i,j,k,inci) = nci_rec(i,j,k)
        micro_field(i,j,k,iqs) = qs_rec(i,j,k)
        micro_field(i,j,k,ins) = ns_rec(i,j,k)
        ! temp array for cloud liquid
        cloudliq(i,j,k) = qcl_tmp(i,j,k)
        ! update the temperature after satuation adjustment
        tabs(i,j,k) = t_tmp(i,j,k)
      end do
    end do
  end do

  ! write the diagnostic absolute temperature, qcl and qci
!  if (rank.eq.1) then
!    do k=1,nzm
!      tabs_diag(k) = sum(dble(tabs(1:nx,1:ny,k)))/float(nx*ny)
!      qcl1_diag(k) = sum(dble(qcl_rec(1:nx,1:ny,k)))/float(nx*ny)
!      qcl_diag(k) = sum(dble(cloudliq(1:nx,1:ny,k)))/float(nx*ny)
!      qci_diag(k) = sum(dble(micro_field(1:nx,1:ny,k,iqci)))/float(nx*ny)
!    write (6,'(i4,1x,f8.2,3f8.4)') k, tabs_diag(k), qcl1_diag(k)*1e3, &
!                                   qcl_diag(k)*1e3, qci_diag(k)*1e3
!    end do
!  endif
      
  if (docloud) call micro_diagnose()

  ! update the dry static energy to liquid/ice static energy
  ! leo (20230914)
  do i=1,nx
    do j=1,ny
      do k=1,nzm
        t(i,j,k) = t(i,j,k) - fac_cond*(qcl(i,j,k)+qpl(i,j,k)) - &
                              fac_sub*(qci(i,j,k)+qpi(i,j,k))
      end do
    end do
  end do

  ! write the diagnostic absolute temperature, qcl and qci
!  if (rank.eq.1) then
!    do k=1,nzm
!      tabs_diag(k) = sum(dble(tabs(1:nx,1:ny,k)))/float(nx*ny)
!      qcl_diag(k) = sum(dble(cloudliq(1:nx,1:ny,k)))/float(nx*ny)
!      qci_diag(k) = sum(dble(micro_field(1:nx,1:ny,k,iqci)))/float(nx*ny)
!      write (6,'(i4,1x,f8.2,3f8.4)') k, tabs_diag(k), qcl1_diag(k)*1e3, &
!                                     qcl_diag(k)*1e3, qci_diag(k)*1e3
!    end do
!  endif

  ! test
!  if (rank.eq.30) then
!    do j=1, ny
!      write(*,*) w(12,j,12)
!    end do
!  end if

  ! print the progress
  if (masterproc) write(*,*) 'done for storing variables'

  call task_barrier()

!================= compute global vertical profile =================

  ! print the progress
  if (masterproc) write(*,*) 'Start computing global averaged vertical profile'

  coef = 1./float(nx*ny)

  do k=1, nzm
    ! initialize
    u0(k) = 0.
    v0(k) = 0.
    t01(k) = tabs0(k) ! for safe, backup the values
    q01(k) = q0(k) ! for safe, backup the values
    t0(k) = 0.
    tabs0(k) = 0.
    q0(k) = 0.
    qn0(k) = 0.
    qp0(k) = 0.
!    p0(k) = 0. ! actually not used

    ! compute the local total
    do j=1,ny
      do i=1,nx
        tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond*(qcl(i,j,k)+qpl(i,j,k)) & 
                                          + fac_sub*(qci(i,j,k)+qpi(i,j,k))
        u0(k) = u0(k) + u(i,j,k)
        v0(k) = v0(k) + v(i,j,k)
!        p0(k) = p0(k) + p(i,j,k)
        t0(k) = t0(k) + t(i,j,k)
        tabs0(k) = tabs0(k) + tabs(i,j,k)
        q0(k) = q0(k) + qv(i,j,k) + qcl(i,j,k) + qci(i,j,k)
        qn0(k) = qn0(k) + qcl(i,j,k) + qci(i,j,k)
        qp0(k) = qp0(k) + qpl(i,j,k) + qpi(i,j,k)
      end do 
    end do
    ! compute the local average
    u0(k) = u0(k) * coef
    v0(k) = v0(k) * coef
!    p0(k) = p0(k) * coef
    t0(k) = t0(k) * coef
    tabs0(k) = tabs0(k) * coef
    q0(k) = q0(k) * coef
    qn0(k) = qn0(k) * coef
    qp0(k) = qp0(k) * coef
  end do

  ! compute the global average
  coef1 = 1./float(nsubdomains)
  do k=1,nzm
    buffer(k,1) = u0(k)
    buffer(k,2) = v0(k)
!    buffer(k,3) = p0(k)
    buffer(k,3) = t0(k)
    buffer(k,4) = tabs0(k)
    buffer(k,5) = q0(k)
    buffer(k,6) = qn0(k)
    buffer(k,7) = qp0(k)
  end do
  call task_barrier()
  call MPI_ALLREDUCE(buffer,buffer1,nzm*7,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  do k=1,nzm
    u0(k) = buffer1(k,1) * coef1
    v0(k) = buffer1(k,2) * coef1
!    p0(k) = buffer1(k,3) * coef1
    t0(k) = buffer1(k,3) * coef1
    tabs0(k) = buffer1(k,4) * coef1
    q0(k) = buffer1(k,5) * coef1
    qn0(k) = buffer1(k,6) * coef1
    qp0(k) = buffer1(k,7) * coef1
  end do
  ! compute global qv0
  qv0(:) = q0(:) - qn0(:)

  if (masterproc) write(*,*) 'Done computing global average vertical profile'

!================== close the data and deallocate ==================
  if (masterproc) then 
    
    ! print the progress
    write(*,*) 'Closing ', iopfilepath

    ! close the dataset
    status = NF_CLOSE(NCID)
    error_code = 0
  
    ! deallocate the arrays
    deallocate(u_in, v_in, w_in, s_in, t_in, &
!-               ppert_in, &
!-               rhow_in, &
               STAT=status)
    if(status.ne.0) then
      write(6,*) 'Processor ', rank, &
                 'Could not de-allocate dimensions in read4dcaminput'
      call task_abort()
    end if
  
!-    deallocate(shf_in, &
!-               lhf_in, taux_in, tauy_in, &
!-               precc_in, precl_in, &
!-               STAT=status)
!-    if(status.ne.0) then
!-      write(6,*) 'Processor ', rank, &
!-                 'Could not de-allocate dimensions in read4dcaminput'
!-      call task_abort()
!-    end if
  
!-    deallocate(fluxbu_in, fluxbv_in, fluxbt_in, fluxbq_in, &
!-               prec_in, &
!-               STAT=status)
!-    if(status.ne.0) then
!-      write(6,*) 'Processor ', rank, &
!-                 'Could not de-allocate dimensions in read4dcaminput'
!-      call task_abort()
!-    end if
  
    ! microphysics arrays
    deallocate(q_in, qcl_in, ncl_in, qr_in, nr_in, qci_in, nci_in, &
               qs_in, ns_in, STAT=status)
    if (status.ne.0) then
      write(6,*) 'Processor ', rank, &
                 'Could not de-allocate dimensions in read4dcaminput'
      call task_abort()
    end if
  
  end if ! if masterproc, close the data

  call task_barrier()

  return
contains
!====================================================================

  subroutine get_netcdf_dimlength( NCID, dimName, dimlength, status, required)
    !====================================================================
    ! subroutine to get the dim length                                  |
    ! input:                                                            |
    ! NCID: id of the dataset (i)                                       |
    ! dimname: name of the dimension (i)                                |
    ! dimlength: variable to store the dim lenght (o)                   |
    ! status: NF error status                                           |
    ! required: check whether the variable is required or optional      |
    !====================================================================
    implicit none
    include 'netcdf.inc'

    ! input/output variables
    integer, intent(in)   :: NCID
    character, intent(in) :: dimName*(*)
    logical, intent(in) :: required

    integer, intent(out) :: status, dimlength

    ! local variables
    integer :: dimID

    ! get variable ID
    STATUS = NF_INQ_DIMID( NCID, dimName, dimID )
    if (STATUS .NE. NF_NOERR ) then
       if(required) then
          if(masterproc) write(6,*) &
               'ERROR(readiopdata.f90):Could not find dimension ID for ', &
               dimName
          STATUS = NF_CLOSE( NCID )
          call task_abort()
       else
          if(masterproc) write(6,*) &
               'Note(readiopdata.f90): No dimension ID for ', dimName
          return
       endif
    endif
    
    ! get the length of the dimension
    STATUS = NF_INQ_DIMLEN( NCID, dimID, dimlength )
    if (STATUS .NE. NF_NOERR ) then
       if(required) then
          if(masterproc) write(6,*) &
               'ERROR(readiopdata.f90):Could not find length of ',dimName
          STATUS = NF_CLOSE( NCID )
          call task_abort()
       else
          if(masterproc) write(6,*) &
               'Note - readiopdata.f90 : Could not find length of ',&
               dimName
       endif
    endif

  end subroutine get_netcdf_dimlength

!====================================================================

  subroutine get_netcdf_var1d_real( NCID, varName, var, use_nf_real, &
                                    status, required)
    !====================================================================
    ! subroutine to get the 1d variables from the dataset               |
    ! input:                                                            |
    ! NCID: id of the dataset (i)                                       |
    ! varName: name of the variable                                     |
    ! var: variables to store the read dataarray                        |
    ! use_nf_real: use real / double to open the data                   |
    ! status: NF error status                                           |
    ! required: check whether the variable is required or optional      |
    !====================================================================
    
    implicit none
    include 'netcdf.inc'
    
    ! input/output variables
    integer, intent(in)   :: NCID
    character, intent(in) :: varName*(*)
    logical, intent(in) :: required, use_nf_real
    
    integer, intent(out) :: status
    real, intent(inout) :: var(:)
    
    ! local variables
    integer :: varID

    ! print the progress
    if (masterproc) write(*,*) 'reading ',varName
    
    ! get variable ID
    STATUS = NF_INQ_VARID( NCID, varName, varID )
    if (STATUS .NE. NF_NOERR ) then
       if(required) then
          if(masterproc) write(6,*) &
               'ERROR(readiopdata.f90):Could not find variable ID for ', &
               varName
          STATUS = NF_CLOSE( NCID )
          call task_abort()
       else
          if(masterproc) write(6,*) &
               'Note(readiopdata.f90): Optional variable ', varName,&
               ' not found in ', TRIM(iopfile)
          return
       endif
    endif
    
    ! read the variables
    if (use_nf_real) then
       STATUS = NF_GET_VAR_REAL( NCID, varID, var )
    else
       STATUS = NF_GET_VAR_DOUBLE( NCID, varID, var )
    endif
    
    if (STATUS .NE. NF_NOERR ) then
       if(required) then
          if(masterproc) write(6,*) &
               'ERROR(readiopdata.f90):Could not find variable ', varName
          STATUS = NF_CLOSE( NCID )
          call task_abort()
       else
          if(masterproc) write(6,*) &
               'Note (readiopdata.f90) : Could not find ', varName
       endif
    endif
    
  end subroutine get_netcdf_var1d_real

!====================================================================

  subroutine get_netcdf_var2d_real( NCID, varName, var, &
                                    use_nf_real, status, required)
    !====================================================================
    ! subroutine to read 2d variables                                   |
    ! ordered in the dims of lat * lon                                  |
    ! to the var array: var(nx, ny)                                     |
    ! input:                                                            |
    ! NCID: id of the dataset (i)                                       |
    ! varName: name of the variable                                     |
    ! var: variables to store the read dataarray                        |
    ! nx, ny, nz: dimension of the x, y, z to read to the var variables |
    ! use_nf_real: use real / double to open the data                   |
    ! status: NF error status                                           |
    ! required: check whether the variable is required or optional      |
    !====================================================================

    implicit none
    include 'netcdf.inc'

    ! input / output variables
    integer, intent(in) :: NCID
    character, intent(in) :: varName*(*)
    logical, intent(in) :: required, use_nf_real

    integer, intent(out) :: status
    real, intent(inout) :: var(:,:)

    ! local variables
    integer :: varID
    character :: dim_name*(NF_MAX_NAME)
    integer :: var_dimIDs(NF_MAX_VAR_DIMS)
    integer :: var_ndims, dim_size, dims_set, i, n, var_type
    logical :: usable_var

    ! print the progress
    if (masterproc) write(*,*) 'reading ',varName

    ! get variable id
    STATUS = NF_INQ_VARID(NCID, varName, varID)
    if (STATUS .NE. NF_NOERR ) then
      if(required) then
        if(masterproc) write(6,*) &
          'ERROR(read4dcaminput.f90):Could not find variable ID for ',&
          varName
        STATUS = NF_CLOSE( NCID )
        call task_abort()
      else
        if(masterproc) write(6,*) &
          'Note(read4dcaminput.f90): Optional variable ', varName,&
          ' not found in ', TRIM(iopfile)
        return
      endif
    endif

    ! check the number of variables dimension
    STATUS = NF_INQ_VARNDIMS(NCID, varID, var_ndims)
    if (var_ndims.ne.2) then
      if (masterproc) then
        write(6,*) &
        'ERROR(read4dcaminput.f90): the number of dimnesion is ', & 
        var_ndims, ' for ', varName, ' not 2 dims'
      endif
      STATUS = -1
      return
    endif

    ! get the dimension ID in the variables
    STATUS = NF_INQ_VARDIMID(NCID, varID, var_dimIDs)
    if (STATUS.ne.NF_NOERR) then
      if (masterproc) then
        write(6,*) & 
        'ERROR(read4dcaminput.f90): cannot get dimension ID for ', &
        varName
      endif
      return
    endif

    ! read the variables
    ! assume the variables are already organized in nx, ny, nz
    ! i.e. in xarray: z, y, x
    if (use_nf_real) then
      STATUS = NF_GET_VAR_REAL(NCID, varID, var)
    else
      STATUS = NF_GET_VAR_DOUBLE(NCID, varID, var)
    endif

    ! error of reading the variables
    if (STATUS.ne.NF_NOERR) then
      if(required) then
        if (masterproc) then
          write(6,*) & 
          'ERROR(read3dcaminput.f90): Could not find variable ', &
          varName
        endif
        STATUS = NF_CLOSE(NCID)
        call task_abort()
      else
        if(masterproc) then
          write(6,*) &
          'Note (read3dcaminput.f90): Could not find ', varName
        endif
      endif
    endif

  end subroutine get_netcdf_var2d_real

!====================================================================

  subroutine get_netcdf_var3d_real(NCID, varName, var, &
                                   use_nf_real, status, required)
    !====================================================================
    ! subroutine to read 3d variables                                   |
    ! ordered in the dims of lev * lat * lon                            |
    ! to the var array: var(nx, ny, nz)                                 |
    ! input:                                                            |
    ! NCID: id of the dataset (i)                                       |
    ! varName: name of the variable                                     |
    ! var: variables to store the read dataarray                        |
    ! nx, ny, nz: dimension of the x, y, z to read to the var variables |
    ! use_nf_real: use real / double to open the data                   |
    ! status: NF error status                                           |
    ! required: check whether the variable is required or optional      |
    !====================================================================

    implicit none
    include 'netcdf.inc'

    ! input / output variables
    integer, intent(in) :: NCID
    character, intent(in) :: varName*(*)
    logical, intent(in) :: required, use_nf_real

    integer, intent(out) :: status
    real, intent(inout) :: var(:,:,:)

    ! local variables
    integer :: varID
    character :: dim_name*(NF_MAX_NAME)
    integer :: var_dimIDs(NF_MAX_VAR_DIMS)
    integer :: start(3), count(3)
    integer :: var_ndims, dim_size, dims_set, i, n, var_type
    logical :: usable_var

    ! print the progress
    if (masterproc) write(*,*) 'reading ',varName

    ! get variable id
    STATUS = NF_INQ_VARID(NCID, varName, varID)
    if (STATUS .NE. NF_NOERR ) then
        if(required) then
           if(masterproc) write(6,*) &
                'ERROR(read4dcaminput.f90):Could not find variable ID for ',&
                varName
           STATUS = NF_CLOSE( NCID )
           call task_abort()
        else
           if(masterproc) write(6,*) &
                'Note(read4dcaminput.f90): Optional variable ', varName,&
                ' not found in ', TRIM(iopfile)
           return
        endif
    endif

    ! check the number of variables dimension
    STATUS = NF_INQ_VARNDIMS(NCID, varID, var_ndims)
    if (var_ndims.ne.3) then
        if (masterproc) then
            write(6,*) &
            'ERROR(read4dcaminput.f90): the number of dimnesion is ', & 
            var_ndims, ' for ', varName, ' not 3 dims'
        endif
        STATUS = -1
        return
    endif

    ! get the dimension ID in the variables
    STATUS = NF_INQ_VARDIMID(NCID, varID, var_dimIDs)
    if (STATUS.ne.NF_NOERR) then
        if (masterproc) then
            write(6,*) & 
            'ERROR(read4dcaminput.f90): cannot get dimension ID for ', &
            varName
        endif
        return
    endif

    ! read the variables
    ! assume the variables are already organized in nx, ny, nz
    ! i.e. in xarray: z, y, x
    if (use_nf_real) then
        STATUS = NF_GET_VAR_REAL(NCID, varID, var)
    else
        STATUS = NF_GET_VAR_DOUBLE(NCID, varID, var)
    endif
    
    ! error of reading the variables
    if (STATUS.ne.NF_NOERR) then
        if(required) then
            if (masterproc) then
                write(6,*) & 
                'ERROR(read3dcaminput.f90): Could not find variable ', &
                varName
            endif
        STATUS = NF_CLOSE(NCID)
        call task_abort()
        else
            if(masterproc) then
                write(6,*) &
                'Note (read3dcaminput.f90): Could not find ', varName
            endif
        endif
    endif
  
  end subroutine get_netcdf_var3d_real

end subroutine read4dcaminput
