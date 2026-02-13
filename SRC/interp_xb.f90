module interp_xb

  use vars
  use grid
  use params

  implicit none

  ! how many grids near to the x-BL in 
  ! each subdomain is computed by interpolation
  integer :: window=3 
  ! only do interpolation at X grids from the y boundary
  integer :: numy=7

  CONTAINS

!==================================================================================

  subroutine interp_xb_set()

    use sgs, only: nsgs_fields,sgs_field,advect_sgs
    use microphysics, only: nmicro_fields,micro_field,&
                            index_water_vapor,flag_precip

    implicit none

    integer :: n ! scalars

    call t_startf('interp_xb')

    ! call the interp_exchange() to do the interpolation
    call interp_xb_exchange(u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,1)
    if (rank.gt.nsubdomains_x*(nsubdomains_y-1)-1) then
      call interp_xb_exchange(v(1:nx,1:ny+1,1:nzm),1,nx,1,(ny+1),nzm,2)
    else
      call interp_xb_exchange(v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,2)
    end if
    call interp_xb_exchange(w(1:nx,1:ny,1:nz),1,nx,1,ny,nz,3)
    call interp_xb_exchange(t(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,4)
    ! sgs fields
    do n=1,nsgs_fields
      if (dosgs.and.advect_sgs) then
        call interp_xb_exchange(sgs_field(1:nx,1:ny,1:nzm,n),1,nx,1,ny,nzm,4+n)
      end if
    end do
    ! microphysics fields
    do n=1,nmicro_fields
      if (n.eq.index_water_vapor &
          .or.docloud.and.flag_precip(n).ne.1 &
          .or.doprecip.and.flag_precip(n).eq.1) then
          call interp_xb_exchange(micro_field(1:nx,1:ny,1:nzm,n),1,nx,1,ny,nzm,4+nsgs_fields+n)
      end if
    end do

    call t_stopf('interp_xb')
  
  end subroutine interp_xb_set

!==================================================================================

  subroutine interp_xb_exchange(finout,dimx1,dimx2,dimy1,dimy2,dimz,nin)

    implicit none

    include 'mpif.h'

    integer :: i,j
    integer,intent(in) :: dimx1,dimx2,dimy1,dimy2,dimz
    real,intent(inout) :: finout(dimx1:dimx2,dimy1:dimy2,1:dimz)
    ! exchange boundary stuff
    integer,intent(in) :: nin
    integer :: iip,iiptemp
    real :: fprepinterp(dimx1-window-1:dimx2+window+1,dimy1:dimy2,1:dimz)
    real :: finterp(dimx1-window-1:dimx2+window+1,dimy1:dimy2,1:dimz)
    real :: buffertoee(1:window+1,dimy1:dimy2,1:dimz), &
            buffertoww(1:window+1,dimy1:dimy2,1:dimz),&
            bufferrec(1:window+1,dimy1:dimy2,1:dimz)
    ! mpi stuff
    integer :: ntagsend,ntagrecv
    integer :: request(2)
    integer :: status(MPI_STATUS_SIZE,2)
    integer :: ierr


    !!! prepare the array for interpolation
    ! initialize
    fprepinterp(:,:,:) = 0.
    ! prepare
    fprepinterp(dimx1:dimx2,dimy1:dimy2,1:dimz) = finout(dimx1:dimx2,dimy1:dimy2,1:dimz)

    !!! to rankee, from ww
    ! prepare the data to send to rankee
    do iip=dimx2-window,dimx2
      iiptemp = iip - (dimx2-window) + 1
      buffertoee(iiptemp,:,:) = finout(iip,dimy1:dimy2,1:dimz)
    end do
    ! compute the tag: nin*1000 + rank of the source
    ! (rank is less than 1000, so to avoid overlaping, nin*1000)
    ntagsend = nin*1000 + rank ! use in send: current rank
    ntagrecv = nin*1000 + rankww ! use in recv: the rank of ww
    ! call the receive first (from ww)
    call MPI_IRECV(bufferrec,(window+1)*(dimy2-(dimy1-1))*dimz,&
                   MPI_REAL,rankww,ntagrecv,MPI_COMM_WORLD,request(1),ierr)
    ! call the send (to rankee)
    call MPI_ISEND(buffertoee,(window+1)*(dimy2-(dimy1-1))*dimz,&
                   MPI_REAL,rankee,ntagsend,MPI_COMM_WORLD,request(2),ierr)
    ! and wait until the send and receive finish
    call MPI_WAITALL(2,request,status,ierr)
    ! save the data to the fprepinterp
    do iip=dimx1-window-1,dimx1-1
      iiptemp = iip - (dimx1-window-1) + 1
      fprepinterp(iip,dimy1:dimy2,1:dimz) = bufferrec(iiptemp,dimy1:dimy2,1:dimz)
    end do
    
    call task_barrier()

    !!! to rankww, from ee
    ! prepare the data to send to rankww
    do iip=dimx1,window+dimx1
      buffertoww(iip,dimy1:dimy2,1:dimz) = finout(iip,dimy1:dimy2,1:dimz)
    end do
    ! compute the tag: nin*1000 + rank of the source
    ! (rank is less than 1000, so to avoid overlaping, nin*1000)
    ntagsend = nin*1000 + rank
    ntagrecv = nin*1000 + rankee
    ! call the receive first (from ee)
    call MPI_IRECV(bufferrec,(window+1)*(dimy2-(dimy1-1))*dimz,&
                   MPI_REAL,rankee,ntagrecv,MPI_COMM_WORLD,request(1),ierr)
    ! call the send (to rankww)
    call MPI_ISEND(buffertoww,(window+1)*(dimy2-(dimy1-1))*dimz,&
                   MPI_REAL,rankww,ntagsend,MPI_COMM_WORLD,request(2),ierr)
    ! and wait until the send and receive finish
    call MPI_WAITALL(2,request,status,ierr)
    ! save the data to the fprepinterp
    do iip=dimx2+1,dimx2+window+1
      iiptemp = iip - (dimx2+1) + 1
      fprepinterp(iip,dimy1:dimy2,1:dimz) = bufferrec(iiptemp,dimy1:dimy2,1:dimz)
    end do

    call task_barrier()

    !!!! call the subroutine to do the interpolation
    call interp_xb_interpx(fprepinterp,finterp,dimx1,dimx2,dimy1,dimy2,dimz)

    !!! save the data back to finout
    ! southmost boundary
    if (rank.lt.nsubdomains_x) then
      do j=dimy1,dimy1+numy-1
        do i=dimx1,dimx1-1+window
          finout(i,j,1:dimz) = finterp(i,j,1:dimz)
        end do
      end do
    end if
    ! northmost boundary
    if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
      do j=dimy2-numy+1,dimy2
        do i=dimx2-window+1,dimx2
          finout(i,j,1:dimz) = finterp(i,j,1:dimz)
        end do
      end do
    end if

  end subroutine interp_xb_exchange

!==================================================================================

  subroutine interp_xb_interpx(fpre,fpost,dimx1,dimx2,dimy1,dimy2,dimz)

    implicit none
    
    integer,intent(in) :: dimx1,dimx2,dimy1,dimy2,dimz
    real, intent(in)  :: fpre(dimx1-window-1:dimx2+window+1,dimy1:dimy2,1:dimz)
    real, intent(out) :: fpost(dimx1-window-1:dimx2+window+1,dimy1:dimy2,1:dimz)
    ! in-subroutine variables
    integer :: i,j,k
    integer :: it,jt
    real :: x0, y0, x1, y1,x


    !!! southmost y boundary
    if (rank.lt.nsubdomains_x) then
      !!! eastern x boundary
      do k=1,dimz
        do j=dimy1,dimy1+numy-1
          do i=dimx1,dimx1-1+window
            call task_rank_to_index(rank,it,jt)
            x0 = (it - window - 1) * dx
            x1 = (it + window + 1 - 1) * dx
            y0 = fpre(dimx1-window-1,j,k)
            y1 = fpre(dimx1-1+window+1,j,k)
            x = (it + i - 1) * dx
            fpost(i,j,k) = (y1 - y0) / (x1 - x0) * (x - x0) + y0
          end do
        end do
      end do
      !!! western x boundary
      do k=1,dimz
        do j=dimy1,dimy1+numy-1
          do i=dimx2-window+1,dimx2
            call task_rank_to_index(rank,it,jt)
            x0 = (it + dimx2 - window -1) * dx
            x1 = (it + dimx2 + window + 1 - 1) * dx
            y0 = fpre(dimx2-window,j,k)
            y1 = fpre(dimx2+window+1,j,k)
            x = (it + i - 1) * dx
            fpost(i,j,k) = (y1 - y0) / (x1 - x0) * (x - x0) + y0
          end do
        end do
      end do
    end if

    !!! northmost y boundary
    if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
      !!! eastern x boundary
      do k=1,dimz
        do j=dimy2-numy+1,dimy2
          do i=dimx1,dimx1-1+window
            call task_rank_to_index(rank,it,jt)
            x0 = (it - window - 1) * dx
            x1 = (it + window + 1 - 1) * dx
            y0 = fpre(dimx1-window-1,j,k)
            y1 = fpre(dimx1-1+window+1,j,k)
            x = (it + i - 1) * dx
            fpost(i,j,k) = (y1 - y0) / (x1 - x0) * (x - x0) + y0
          end do
        end do
      end do
      !!! western x boundary
      do k=1,dimz
        do j=dimy2-numy+1,dimy2
          do i=dimx2-window+1,dimx2
            call task_rank_to_index(rank,it,jt)
            x0 = (it + dimx2 - window -1) * dx
            x1 = (it + dimx2 + window + 1 - 1) * dx
            y0 = fpre(dimx2-window,j,k)
            y1 = fpre(dimx2+window+1,j,k)
            x = (it + i - 1) * dx
            fpost(i,j,k) = (y1 - y0) / (x1 - x0) * (x - x0) + y0
          end do
        end do
      end do
    end if

  end subroutine interp_xb_interpx

end module interp_xb