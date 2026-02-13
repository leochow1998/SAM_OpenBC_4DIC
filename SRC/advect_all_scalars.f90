subroutine advect_all_scalars()

  use vars
  use microphysics
  use sgs
  use tracers
  use params, only: dotracers,doopeny ! leo: openbl (last 2) 
  use domain, only: nsubdomains_x, nsubdomains_y ! leo: openbl
  use grid, only: rank
  use open_boundary_y
  implicit none
  real dummy(nz)
  integer k
  
  ! leo: openbl
  ! define temperary array to store the scalors before
  ! advection
  real :: scalarbackup(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real :: openyout(dimx1_s:dimx2_s,6,nzm) ! 1-3 is y=1 to 3, 4-6 is y=ny-2 to ny
  integer :: j, jtemp
  integer :: it, jt

!---------------------------------------------------------
!      advection of scalars :
     
     ! backup t
     if (doopeny) then
       scalarbackup(:,:,:) = t(:,:,:)
     end if
     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)
     ! advection of scalar at y bl
     if (doopeny) then
          call advection_scalar(scalarbackup,tprev,ubackup,vbackup,wbackup,openyout)
!          call advection_scalar(scalarbackup,ubackup,vbackup,wbackup,openyout)
          ! store it back to the scalar
          if (rank.lt.nsubdomains_x) then
               do j=1,3
                    t(:,j,:) = openyout(:,j,:)
               end do
          else if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
               do j=ny-2,ny
                    jtemp = j - ny + 6
                    t(:,j,:) = openyout(:,jtemp,:)
               end do
          end if
     end if

!
!    Advection of microphysics prognostics:
!

     do k = 1,nmicro_fields
       if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
         if (doopeny) then
             ! backup micro_field
             scalarbackup(:,:,:) = micro_field(:,:,:,k)
         end if
         call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
         ! advection of scalar at y bc
         if (doopeny) then
             call advection_scalar(scalarbackup,micro_fieldprev(:,:,:,k,:),ubackup,vbackup,wbackup,openyout)
!             call advection_scalar(scalarbackup,ubackup,vbackup,wbackup,openyout)
             ! store it back to the scalar
             if (rank.lt.nsubdomains_x) then
               do j=1,3
                    micro_field(:,j,:,k) = openyout(:,j,:)
               end do
             else if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
               do j=ny-2,ny
                    jtemp = j - ny + 6
                    micro_field(:,j,:,k) = openyout(:,jtemp,:)
               end do
             end if
         end if
       end if
     end do

!
!    Advection of sgs prognostics:
!

     if(dosgs.and.advect_sgs) then
       do k = 1,nsgs_fields
          if (doopeny) then
               ! backup sgs_field
               scalarbackup(:,:,:) = sgs_field(:,:,:,k)
          end if
           call advect_scalar(sgs_field(:,:,:,k),sgsadv(:,k),sgswle(:,k),dummy,dummy,dummy,.false.)
          ! advection of scalar at y bc
           if (doopeny) then
               call advection_scalar(scalarbackup,sgs_fieldprev(:,:,:,k,:),ubackup,vbackup,wbackup,openyout)
!               call advection_scalar(scalarbackup,ubackup,vbackup,wbackup,openyout)
               ! store it back to the scalar
               if (rank.lt.nsubdomains_x) then
                    do j=1,3
                         sgs_field(:,j,:,k) = openyout(:,j,:)
                    end do
               else if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
                    do j=ny-2,ny
                         jtemp = j - ny + 6
                         sgs_field(:,j,:,k) = openyout(:,jtemp,:)
                    end do
               end if
          end if
       end do
     end if


!
!   Precipitation fallout:
!
    if(doprecip) then

       total_water_prec = total_water_prec + total_water()

       call micro_precip_fall()

       total_water_prec = total_water_prec - total_water()


    end if

 ! advection of tracers:

     if(dotracers) then

        do k = 1,ntracers
         call advect_scalar(tracer(:,:,:,k),tradv(:,k),trwle(:,k),dummy,dummy,dummy,.false.)
        end do

     end if

end subroutine advect_all_scalars
