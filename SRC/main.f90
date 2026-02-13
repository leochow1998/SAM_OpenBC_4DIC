program crm

!       Main module.

use vars
use hbuffer
use microphysics
use sgs
use tracers
use movies, only: init_movies
use damping_lateral ! leo
use open_boundary_y ! leo: ybc
use grid, only: rank ! leo: ybc
!use interp_xb ! leo: interpolation x-boundary
implicit none

integer k, icyc, nn, nstatsteps

integer :: id,jd,kd
!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init() 
!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()	
!------------------------------------------------------------------
! Initialize timing library.  2nd arg 0 means disable, 1 means enable

   call t_setoptionf (1, 0)
   call t_initializef ()

   call t_startf ('total')
   call t_startf ('initialize')
!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants

!------------------------------------------------------------------
! Initialize or restart from the save-dataset:

if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables
elseif(nrestart.eq.1) then
   call read_all()
   call read_restart_openbc() ! leo: read the previous timestep values in openybc
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call sgs_init()
   call micro_init()  !initialize microphysics
elseif(nrestart.eq.2) then  ! branch run
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call setparm() ! overwrite the parameters
   call sgs_init()
   call micro_init()  !initialize microphysics
   nstep = 0
   day0 = day
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif

call init_movies()
call stat_2Dinit()
call tracers_init() ! initialize tracers
call setforcing()
! leo: initialize Rayleigh damping sponge layer at yBL
if (dodampinglateral) call damping_lateral_init()
if(masterproc) call printout()
!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()
	
!------------------------------------------------------------------
nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0

!------------------------------------------------------------------
! initialize the y boundary condition ! leo: ybc (debug)
! doopeny = .true.
 if (doopeny) then
   dowally = .true.
 end if
!------------------------------------------------------------------

call t_stopf ('initialize')
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------

do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  nstep = nstep + 1
  time = time + dt
  day = day0 + nstep*dt/86400.
  nelapse = nelapse - 1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()

  total_water_before = 0.
  total_water_after = 0.
  total_water_evap = 0.
  total_water_prec = 0.
  total_water_ls = 0.

  do icyc=1,ncycle

     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif

     !bloss:make special statistics flag for radiation,since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.

!--------------------------------------------
!    see whether to do the zonal smoothing
!    near the y boundary for nudging
     call average_to_yb()

!--------------------------------------------
!    initialize the openbl and backup the 
!    previous timestep for the prognostic 
!    variables
!    for leapfrog scheme of BL scalar
     call openbl_init()
     call prevbackup()

!    also call the y-boundary spurious vertical motion
!    diagnostics for testing the stability
     call openbc_diagnose()

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

!-----------------------------------------------------------

     total_water_before = total_water_before + total_water()

!-----------------------------------------------------------
!       Buoyancy term:

     call buoyancy()

!------------------------------------------------------------

     total_water_ls =  total_water_ls - total_water()

!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()

!----------------------------------------------------------
!       Nadging to sounding:

     call nudging()

!----------------------------------------------------------
!    leo: sponge-layer at the yBL:

     if (dodampinglateral) call damping_lateral_proc()

!----------------------------------------------------------

!   	spange-layer damping near the upper boundary:

     if(dodamping) call damping()

!----------------------------------------------------------

     total_water_ls =  total_water_ls + total_water()

!---------------------------------------------------------
!   Ice fall-out
   
      if(docloud) then
          call ice_fall()
      end if

!----------------------------------------------------------
!     Update scalar boundaries after large-scale processes:

     call boundaries(3)

!---------------------------------------------------------
!     Update boundaries for velocities:

     ! yBL: backup the northmost v before updating the BL
     if (doopeny) then
          if (rank.gt.nsubdomains_x*(nsubdomains_y-1)-1) then
               vnorthmost(:,1,:) = v(:,ny+1,:)
          end if
     end if
     
      call boundaries(0)

      ! yBL: from the backup, store back the v at northmost grids
      ! after exchange of boundaries
      if (doopeny) then
          if (rank.gt.(nsubdomains_x*(nsubdomains_y-1)-1)) then
               v(:,ny+1,:) = vnorthmost(:,1,:)
          end if
      end if

!-----------------------------------------------
!     surface fluxes:

     if(dosurface) call surface()

!-----------------------------------------------------------
!  SGS physics:

     if (dosgs) call sgs_proc()

!----------------------------------------------------------
!     Fill boundaries for SGS diagnostic fields:

     call boundaries(4)

!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!----------------------------------------------------------
!	SGS effects on momentum:

     if(dosgs) call sgs_mom()

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()

!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
!  Note that at the end of the call, the velocities are in nondimensional form.
	 
     call adams()

!----------------------------------------------------------
!     Update boundaries for all prognostic scalar fields for advection:

     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :

     call advect_all_scalars()

!---------------------------------------------------------
! leo: rayleigh lateral damping for scalars :

     if (dodampinglateral) call damping_lateral_scalar_proc()

!-----------------------------------------------------------
!    Convert velocity back from nondimensional form:

      call uvw()

!----------------------------------------------------------
!     Update boundaries for scalars to prepare for SGS effects:

     call boundaries(3)
   
!---------------------------------------------------------
!      SGS effects on scalars :

     if (dosgs) call sgs_scalars()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:

      if(docloud.or.dosmoke) call micro_proc()

!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()

!-----------------------------------------------------------
!	Radiation

      if(dolongwave.or.doshortwave) then 
	call radiation()     
      end if

!-----------------------------------------------------------
!    Compute diagnostic fields:

      call diagnose()

!----------------------------------------------------------

      total_water_after = total_water_after + total_water()

!----------------------------------------------------------
! clear all backup storage due to openbl
      if (doopeny) then
           dvdtnorthmost(:,:,:) = 0.
           dvdtsouthmost(:,:,:) = 0.
           vnorthmost(:,:,:) = 0.
           ubackup(:,:,:) = 0.
           vbackup(:,:,:) = 0.
           wbackup(:,:,:) = 0.
      end if

!--------------------------------------------
!    interpolate to replace the prognostic
!    value at the x boundary of each subdomain
     
!      call interp_xb_set()

!----------------------------------------------------------
! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

   end do ! icycle	
          
!----------------------------------------------------------
!  collect statistics, write save-file, etc.

   call stepout(nstatsteps)

!----------------------------------------------------------
!----------------------------------------------------------

end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------

   call t_stopf('total')
   if(masterproc) call t_prf(rank)

call task_stop()

end program crm
