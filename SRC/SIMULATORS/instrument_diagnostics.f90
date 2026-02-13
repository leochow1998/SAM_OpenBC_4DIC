module instrument_diagnostics 
  !
  ! Compute results from ISCCP and MODIS simulators at pixel and domain-mean levels
  !    Robert Pincus, Feb 2010
  !
  use vars,   only: nx, ny, nz, nzm, sstxy, t00, qv, tabs, qcl, qci
  use grid,   only: pres, presi, doisccp, domodis, domisr, dosimfilesout 
  use params, only: coszrs, ggr
  use microphysics, only: Get_reffc, Get_reffi
  use modis_simulator, only: modis_L2_simulator, modis_L3_simulator, &
                    modis_num_taus => numTauHistogramBins, & 
                    modis_num_ctps => numPressureHistogramBins  
  implicit none
  private
  
  !
  ! Parameters we might want to revisit 
  !
  real, parameter :: emsfc_lw = 0.99
  real, parameter :: solarZenithCosLimit = 0.2 
  real, parameter :: minMixingRatio = 1.e-9
  
  integer, parameter :: numIsccpPressureIntervals     = 7, &
                        numIsccpOpticalDepthIntervals = 7, &
  ! Dimensions of MISR joint histograms
                        misr_num_cths = 16
  real,    parameter :: minIsccpOpticalDepth = 0.3
  real, dimension(numIsccpPressureIntervals + 1),      parameter ::       &
           isccpPressureBinEdges = (/ tiny(minIsccpOpticalDepth),         &
                                      180., 310., 440., 560., 680., 800., &
                                      huge(minIsccpOpticalDepth) /)
  real, dimension(numIsccpOpticalDepthIntervals + 1),  parameter ::                    &
          isccpOpticalDepthBinEdges = (/ tiny(minIsccpOpticalDepth),                   &
                                        minIsccpOpticalDepth, 1.3, 3.6, 9.4, 23., 60., &
                                        huge(minIsccpOpticalDepth) /)

  integer, parameter :: isccp_top_height_direction = 2, isccp_top_height = 1, &
                        isccp_num_taus = numIsccpOpticalDepthIntervals, &
                        isccp_num_ctps = numIsccpPressureIntervals

  real, parameter :: fillvalue = -1.     ! fill value for reduced grid and others

  real fq_isccp(isccp_num_taus,isccp_num_ctps)  
  real isccp_totalcldarea 
  real isccp_lowcldarea ! low (p>700 mb) cloud fractions
  real isccp_midcldarea ! mid (400<p<700 mb) cloud fractions
  real isccp_hghcldarea ! high (p<400 mb) cloud fractions
   ! The following three means are averages over the cloudy areas only.  
  real isccp_meantaucld !  mean optical thickness (dimensionless)
                                        !  linear averaging in albedo performed.
  real isccp_meanptop !  mean cloud top pressure (mb) - linear averaging
  real isccp_meanalbedocld ! mean cloud albedo - linear averaging
  real isccp_meantb  ! mean all-sky 10.5 micron brightness temperature 
  real isccp_meantbclr  ! mean clear-sky 10.5 micron brightness temperature 
                                        !  in cloud top pressure.
  real fq_modis(modis_num_taus, modis_num_ctps)  
  real modis_totalcldarea 
  real modis_totalcldarea_l 
  real modis_totalcldarea_i 
  real modis_lowcldarea 
  real modis_midcldarea 
  real modis_hghcldarea 
   ! The following three means are averages over the cloudy areas only.  
  real modis_meantaucld 
  real modis_meantaucld_l 
  real modis_meantaucld_i 
  real modis_meanrel 
  real modis_meanrei 
  real modis_meanptop
  real modis_meanlwp 
  real modis_meaniwp 

  real fq_misr(isccp_num_taus, misr_num_cths) 
  real misr_meanztop, misr_totalcldarea
  real, parameter ::  MISR_CTH_boundaries(misr_num_cths) = (/  0., 0.5, 1., 1.5, 2., 2.5, 3., &
                                     4., 5., 7., 9., 11., 13., 15., 17., 99. /)

  integer nisccp, nmodis, nmisr ! number of collected samples. Greater than 0 only if at least one
                     ! day-time sample was obtained
  integer nmodis_l, nmodis_i
  integer nsun ! number of day-time samples.

  public :: zero_instr_diag, compute_instr_diags, isccp_write, modis_write, misr_write

contains 

  ! -------------------------------------------------------------------------------
subroutine zero_instr_diag()
    fq_isccp = 0.
    isccp_totalcldarea = 0.
    isccp_lowcldarea = 0.
    isccp_midcldarea = 0.
    isccp_hghcldarea = 0.
    isccp_meantaucld = 0.
    isccp_meanptop = 0.
    isccp_meanalbedocld = 0.
    isccp_meantb = 0.
    isccp_meantbclr = 0.
    fq_modis = 0.
    modis_totalcldarea = 0.
    modis_totalcldarea_l = 0.
    modis_totalcldarea_i = 0.
    modis_lowcldarea = 0.
    modis_midcldarea = 0.
    modis_hghcldarea = 0.
    modis_meantaucld = 0.
    modis_meantaucld_l = 0.
    modis_meantaucld_i = 0.
    modis_meanptop = 0.
    modis_meanlwp = 0.
    modis_meaniwp = 0.
    modis_meanrel = 0.
    modis_meanrei = 0.
    fq_misr = 0.
    misr_totalcldarea = 0.
    misr_meanztop = 0.
    nisccp = 0
    nmodis = 0
    nmodis_l = 0
    nmodis_i = 0
    nmisr = 0
    nsun = 0
end subroutine zero_instr_diag
  ! -------------------------------------------------------------------------------
  subroutine init_instr_diags
    !
    ! Stub for the moment 
    ! 
  end subroutine init_instr_diags
  ! -------------------------------------------------------------------------------
  subroutine end_instr_diags
    !
    ! Stub for the moment 
    ! 
  end subroutine end_instr_diags
  ! -------------------------------------------------------------------------------
  subroutine compute_instr_diags
    !
    ! Call instrument simulators on local SAM domain
    !   Simulators expect input arrays with dimension npoints, nsubcols, nlevels
    !   Our entire domain would be considered a single point 
    !
    
    use vars, only: z, pres, pres0, tabs0, compute_reffc, compute_reffi 
    use params, only: ocean
    ! Model state 
    real, dimension(nx * ny, nzm) :: re_l, re_i, q_l, q_i
    real, dimension(nx * ny, nzm) :: lwp, iwp
    ! Domain means ("gridpoint" as opposed to "subgrid")
    real, dimension(1, nzm)   :: p_mid, z_mid, q_v, temp
    real, dimension(1, nzm+1) :: p_int
    real, dimension(nzm)      :: deltap
    
    ! Inputs for simulators
    real,    dimension(nx * ny, nzm) :: dtau_l, dtau_i, emiss 
    integer, dimension(1) :: sunlit = 1
    real,    dimension(1) :: skt

    ! Pixel-level outputs
    real,    dimension(1, nx * ny) :: isccp_tau, isccp_ctp
    real,    dimension(1, nx * ny) :: modis_ctp, modis_tau, modis_re
    integer, dimension(1, nx * ny) :: modis_phase
    
    ! Domain mean output
    real, dimension(1) :: isccp_cf, isccp_cf_low, isccp_cf_mid, isccp_cf_high, &
                          isccp_mean_ctp, isccp_mean_tau, &
                          isccp_mean_albedo, isccp_mean_tb, isccp_mean_tbclr
    real, dimension(1, isccp_num_taus, isccp_num_ctps) :: isccp_jpdf
    
    real, dimension(1) :: modis_cf,            modis_cf_l,            modis_cf_i,            &
                          modis_cf_high,       modis_cf_mid,          modis_cf_low,          &
                          modis_mean_tau,      modis_mean_tau_l,      modis_mean_tau_i,      &
                          modis_mean_log10tau, modis_mean_log10tau_l, modis_mean_log10tau_i, &
                          modis_mean_re_l,     modis_mean_re_i,                              &
                          modis_mean_ctp,      modis_lwp,             modis_iwp
    real, dimension(1, modis_num_taus, modis_num_ctps) :: modis_jpdf
    real, dimension(1, isccp_num_taus, misr_num_cths) :: misr_jpdf
    real, dimension(1,                 misr_num_cths) :: misr_dist_layertops
    real, dimension(1) :: misr_mean_ztop, misr_cldarea

    real rel(nzm), rei(nzm)
    integer k

    if(.not.(doisccp.or.domodis.or.domisr)) return

    ! ----------------------------------------------------------------
    if(.true.) then  ! just do the diagnostics all the time
!    if(coszrs >= solarZenithCosLimit) then 

      nsun = nsun + 1
      !
      ! The physical description - mixing ratio and particle size
      !
      q_l  = reshape(qcl, shape = (/ nx * ny, nzm /) ) 
      q_i  = reshape(qci, shape = (/ nx * ny, nzm /) ) 
        call cldefrint(1,1,0.,tabs0(nzm:1:-1),rel(nzm:1:-1),rei(nzm:1:-1), &
                       pres0*100., pres(nzm:1:-1)*100., merge(0.,1.,OCEAN),0.,0.)
      if(compute_reffc) then
        re_l = reshape(Get_reffc(), shape = (/ nx * ny, nzm /) ) 
      else
        re_l = spread(rel(:), dim = 1, ncopies = nx*ny)
      end if
      if(compute_reffi) then
        re_i = reshape(Get_reffi(), shape = (/ nx * ny, nzm /) ) 
      else
        re_i = spread(rei(:), dim = 1, ncopies = nx*ny)
      end if
      !
      ! Set floors on thinnest clouds to include
      !
      where(q_l < minMixingRatio) re_l = 0. 
      where(q_i < minMixingRatio) re_i = 0. 
      where(q_l < minMixingRatio) q_l  = 0. 
      where(q_i < minMixingRatio) q_i  = 0. 
      
      !
      ! Domain means - stand-ins for "grid point averages" 
      !
      z_mid(1, :) = z(:nzm)
      p_mid(1, :) = pres(:)  * 100. ! Convert from mb to Pa
      p_int(1, :) = presi(:) * 100. 
      deltap(1:nzm) = (presi(1:nzm) - presi(2:nz)) *100.
      q_v  (1, :) = sum(reshape(qv,   shape = (/ nx * ny, nzm /)), dim = 1) / real(nx * ny) 
      temp (1, :) = sum(reshape(tabs, shape = (/ nx * ny, nzm /)), dim = 1) / real(nx * ny) ! Temperature
      skt  (1)    = sum(sstxy(1:nx,1:ny))/real(nx * ny) + t00
      
      !
      ! Instrument simulators are ordered top to bottom; SAM is reverse
      !
      q_l(:, :)  = q_l (:, nzm:1:-1)
      q_i(:, :)  = q_i (:, nzm:1:-1)
      re_l(:, :) = re_l(:, nzm:1:-1) * 1.e-6
      re_i(:, :) = re_i(:, nzm:1:-1) * 1.e-6
  
      z_mid(:, :) = z_mid(:, nzm:1:-1) 
      p_mid(:, :) = p_mid(:, nzm:1:-1) 
      p_int(:, :) = p_int(:, nz :1:-1)
      deltap(  :) = deltap(nzm:1:-1)
      q_v(1, :)   = q_v  (1, nzm:1:-1)
      temp(1, :)  = temp (1, nzm:1:-1)
      
      ! ----------------------------------------------------------------
      !
      ! Compute optical depth and emissivity from size and mass
      !   Liquid and ice water paths are needed first
      !
      lwp(:, :) = spread(deltap(:), dim = 1, ncopies = nx * ny)/ggr * q_l(:, :) 
      iwp(:, :) = spread(deltap(:), dim = 1, ncopies = nx * ny)/ggr * q_i(:, :) 
      dtau_l(:, :)  = compute_tau_l(lwp, re_l)
      dtau_i(:, :)  = compute_tau_i(iwp, re_i)
      emiss(:, :) = compute_emiss(lwp, iwp, re_l, re_i)
!      print*,'liq: z p T ql lwp tau_l  rel'
!      write(6,'(7f12.4)') (z_mid(1,k),p_mid(1,k),temp(1,k),q_l(1,k)*1.e3,lwp(1,k)*1.e3,dtau_l(1,k),re_l(1,k)*1.e6,k=1,nzm)
!      print*,'ice:  z p T  qi lwp tau_i  rei'
!      write(6,'(7f12.4)') (z_mid(1,k),p_mid(1,k),temp(1,k),q_i(1,k)*1.e3,iwp(1,k)*1.e3,dtau_i(1,k),re_i(1,k)*1.e6,k=1,nzm)
      
      ! ----------------------------------------------------------------
      !
      ! Invoke ISCCP simulator
      !
      if(doisccp) then

          call icarus(0, 0, 1, sunlit, nzm, nx*ny, &
                  p_mid, p_int, q_v,             & 
                  reshape(dtau_l + dtau_i, shape = (/ 1, nx*ny, nzm /) ), & 
                  isccp_top_height, isccp_top_height_direction, &
                  skt, emsfc_lw, temp,           &
                  reshape(emiss,           shape = (/ 1, nx*ny, nzm /) ), &
                  ! Domain-mean outputs
                  isccp_jpdf,                    &
                  isccp_cf, isccp_cf_low, isccp_cf_mid, isccp_cf_high, &
                  isccp_mean_ctp, isccp_mean_tau, isccp_mean_albedo, isccp_mean_tb, isccp_mean_tbclr, &
                  ! Pixel-level outputs
                  isccp_tau, isccp_ctp)

          fq_isccp(:,:) = fq_isccp(:,:) + isccp_jpdf(1,:,:)
          isccp_totalcldarea = isccp_totalcldarea + isccp_cf(1)
          isccp_lowcldarea = isccp_lowcldarea + isccp_cf_low(1)
          isccp_midcldarea = isccp_midcldarea + isccp_cf_mid(1)
          isccp_hghcldarea = isccp_hghcldarea + isccp_cf_high(1)
          if(isccp_cf(1).gt.0.) then
             isccp_meanptop = isccp_meanptop + isccp_mean_ctp(1)*isccp_cf(1)
             isccp_meanalbedocld = isccp_meanalbedocld + isccp_mean_albedo(1)*isccp_cf(1)
             isccp_meantb = isccp_meantb + isccp_mean_tb(1)*isccp_cf(1)
             isccp_meantbclr = isccp_meantbclr + isccp_mean_tbclr(1)*isccp_cf(1)
             nisccp = nisccp + 1
          endif

      
      end if ! doisccp
      ! ----------------------------------------------------------------
      !
      ! Invoke MODIS simulator
      !
      if(domodis) then
          call modis_L2_simulator(temp(1, :), p_mid(1, :), p_int(1, :), &
                              dtau_l, dtau_i, re_l, re_i,           & 
                              isccp_tau(1, :), isccp_ctp(1, :),     &
                              modis_phase(1, :), modis_ctp(1, :), modis_tau(1, :), modis_re(1, :))
   
         !
         ! MODIS averaging is done separately
         !
         call modis_L3_simulator(modis_phase, modis_ctp, modis_tau, modis_re,                       &
                              modis_cf,            modis_cf_l,            modis_cf_i,            &
                              modis_cf_high,       modis_cf_mid,          modis_cf_low,          &
                              modis_mean_tau,      modis_mean_tau_l,      modis_mean_tau_i,      &
                              modis_mean_log10tau, modis_mean_log10tau_l, modis_mean_log10tau_i, &
                                                   modis_mean_re_l,       modis_mean_re_i,       &
                              modis_mean_ctp,                                                    &
                                                   modis_lwp,             modis_iwp,             &    
                              modis_jpdf)
         fq_modis(:,:) = fq_modis(:,:) + modis_jpdf(1,:,:)
         modis_totalcldarea = modis_totalcldarea + modis_cf(1)
         modis_lowcldarea = modis_lowcldarea + modis_cf_low(1)
         modis_midcldarea = modis_midcldarea + modis_cf_mid(1)
         modis_hghcldarea = modis_hghcldarea + modis_cf_high(1)
         modis_totalcldarea_l = modis_totalcldarea_l + modis_cf_l(1)
         modis_totalcldarea_i = modis_totalcldarea_i + modis_cf_i(1)
!         print*,'modis:',modis_cf_low,modis_cf_mid,modis_cf_high,modis_mean_tau,modis_mean_re_l*1.e6,modis_mean_re_i*1.e6
         if(modis_cf(1).gt.0.) then
           modis_meanptop = modis_meanptop + modis_mean_ctp(1)*0.01*modis_cf(1)
           modis_meantaucld = modis_meantaucld + modis_mean_tau(1)*modis_cf(1)
           nmodis = nmodis + 1
         end if
         if(modis_cf_l(1).gt.0.) then
           modis_meantaucld_l = modis_meantaucld_l + modis_mean_tau_l(1)*modis_cf_l(1)
           modis_meanlwp = modis_meanlwp + modis_lwp(1)*1000.*modis_cf_l(1)
           modis_meanrel = modis_meanrel + modis_mean_re_l(1)*1.e6*modis_cf_l(1)
           nmodis_l = nmodis_l + 1
         end if
         if(modis_cf_i(1).gt.0.) then
           modis_meantaucld_i = modis_meantaucld_i + modis_mean_tau_i(1)*modis_cf_i(1)
           modis_meaniwp = modis_meaniwp + modis_iwp(1)*1000.*modis_cf_i(1)
           modis_meanrei = modis_meanrei + modis_mean_re_i(1)*1.e6*modis_cf_i(1)
           nmodis_i = nmodis_i + 1
         end if

      end if ! domodis
      ! ----------------------------------------------------------------
      !
      ! MISR simulator - aggregate output only
      !
      if(domisr) then
         call MISR_simulator(1, nzm, nx*ny, sunlit,                                  &
                          z_mid, temp,                                            &
                          reshape(dtau_l + dtau_i, shape = (/ 1, nx*ny, nzm /) ), &
                          misr_jpdf,  misr_dist_layertops, misr_mean_ztop, misr_cldarea)
!         print*,'misr_cldarea:',misr_cldarea
          fq_misr(:,:) = fq_misr(:,:) + misr_jpdf(1,:,:)
          misr_totalcldarea = misr_totalcldarea + misr_cldarea(1)
         if(misr_cldarea(1).gt.0.) then
           misr_meanztop = misr_meanztop + misr_mean_ztop(1)*misr_cldarea(1)
           nmisr = nmisr + 1
!           print*,'misr_layertops:',misr_dist_layertops
!           print*,'misr z top:',misr_mean_ztop
         end if

      !

      end if
    end if ! coszrs > ....
    
  end subroutine compute_instr_diags


  ! -------------------------------------------------------------------------------
  elemental function compute_tau_l(lwp, re_l)
    real, intent(in) :: lwp, re_l
    real             :: compute_tau_l
    !
    ! Diagnose optical thickness (nominally at 0.67 microns) from water clouds
    !
    ! This version comes from radcswmx.f90 in the CAM3 radiation package, which 
    !   expects L/IWP in g/m2 and particle sizes in microns
    !   Our inputs are in MKS
    
    real, parameter :: abarl = 2.817e-02, bbarl = 1.305, & 
                       abari = 3.448e-03, bbari = 2.431

    compute_tau_l = 0.
    if(re_l > 0.) compute_tau_l= (abarl + bbarl/(re_l*1.e6)) * lwp*1000.  
    
  end function compute_tau_l

  ! -------------------------------------------------------------------------------
  elemental function compute_tau_i(iwp, re_i)
    real, intent(in) :: iwp, re_i
    real             :: compute_tau_i
    !
    ! Diagnose optical thickness (nominally at 0.67 microns) from ice clouds
    !
    ! This version comes from radcswmx.f90 in the CAM3 radiation package, which 
    !   expects L/IWP in g/m2 and particle sizes in microns
    !   Our inputs are in MKS
    
    real, parameter :: abarl = 2.817e-02, bbarl = 1.305, & 
                       abari = 3.448e-03, bbari = 2.431

    compute_tau_i = 0.
    if(re_i > 0.) compute_tau_i = (abari + bbari/(re_i*1.e6)) * iwp*1000. 
    
  end function compute_tau_i

  ! -------------------------------------------------------------------------------
  elemental function compute_emiss(lwp, iwp, re_l, re_i)
    real, intent(in) :: lwp, iwp, re_l, re_i
    real :: compute_emiss
    !
    ! Diagnose longwave emissivity (nominally at 10.5 microns) 
    !
    ! This version comes from pkg_cldoptics.f90 in the CAM3 radiation package, which 
    !   expects L/IWP in g/m2 and particle sizes in microns
    !   Our inputs are in MKS
    
    compute_emiss = 1. - exp(-(0.15 * (lwp*1000.) + 1.66 * (0.005 + 1./max(25.,(re_i*1.e6)))*(iwp*1000.)))

  end function compute_emiss

  ! -------------------------------------------------------------------------------

subroutine isccp_write()

   use grid, only: nstep, nstat, day, caseid, nsubdomains, dompi, case, masterproc
   use vars, only: s_acldisccp, s_acldlisccp, s_acldmisccp, s_acldhisccp, &
                   s_tbisccp, s_tbclrisccp, s_cldtauisccp, s_cldalbisccp, s_ptopisccp
   implicit none
   integer, parameter :: ntau = numIsccpOpticalDepthIntervals , &
                         npres = numIsccpPressureIntervals
   integer, parameter :: buf_len = ntau*npres + 8
   real coef, buf(buf_len), buf1(buf_len)
   integer, parameter :: ntape = 88
   integer i,j,nisccp_max,nsun_max, itmp1(2), itmp2(2)
   real dummy
   integer nsteplast

   if(.not.doisccp) return

   nisccp_max = nisccp
   nsun_max = nsun
   if(dompi) then
!bloss: use arrays to match type of dummy argument in task_max_integer
!bloss     call task_max_integer (nisccp,nisccp_max,1)
!bloss     call task_max_integer (nsun,nsun_max,1)
     itmp1(1) = nisccp
     itmp1(2) = nsun
     call task_max_integer (itmp1,itmp2,2)
     nisccp_max = itmp2(1)
     nsun_max = itmp2(2)
   end if

   if(nisccp_max.gt.0.and.nsun.eq.nsun_max) then
     if(dompi) then
       if(nisccp.gt.0) then
         buf(1:ntau*npres) = reshape(fq_isccp, shape = (/ ntau * npres /) )
         buf(ntau*npres+1) = isccp_totalcldarea
         buf(ntau*npres+2) = isccp_meanalbedocld
         buf(ntau*npres+3) = isccp_meanptop
         buf(ntau*npres+4) = isccp_lowcldarea
         buf(ntau*npres+5) = isccp_midcldarea
         buf(ntau*npres+6) = isccp_hghcldarea
         buf(ntau*npres+7) = isccp_meantb
         buf(ntau*npres+8) = isccp_meantbclr
         coef = 1./float(nsubdomains*nisccp)
         buf = buf * coef
       else
         buf(1:buf_len) = 0.
       endif
       call task_sum_real(buf,buf1,buf_len)
       fq_isccp(:,:) = reshape(buf1(1:ntau*npres),shape = (/ntau,npres/))
       isccp_totalcldarea = buf1(ntau*npres+1)
       isccp_meanalbedocld = buf1(ntau*npres+2)
       isccp_meanptop = buf1(ntau*npres+3)
       isccp_lowcldarea = buf1(ntau*npres+4)
       isccp_midcldarea = buf1(ntau*npres+5)
       isccp_hghcldarea = buf1(ntau*npres+6)
       isccp_meantb = buf1(ntau*npres+7)
       isccp_meantbclr = buf1(ntau*npres+8)
     else
       coef = 1./float(nisccp)
       fq_isccp(:,:) = fq_isccp(:,:) * coef
       isccp_totalcldarea = isccp_totalcldarea * coef
       isccp_lowcldarea = isccp_lowcldarea * coef
       isccp_midcldarea = isccp_midcldarea * coef
       isccp_hghcldarea = isccp_hghcldarea * coef
       isccp_meanptop = isccp_meanptop * coef
       isccp_meanalbedocld = isccp_meanalbedocld * coef
       isccp_meantb =  isccp_meantb * coef
       isccp_meantbclr =  isccp_meantbclr * coef
     end if
     if(isccp_totalcldarea.gt.0.01) then
       coef = 1./(isccp_totalcldarea)
       isccp_meanalbedocld = isccp_meanalbedocld * coef
       isccp_meanptop =  isccp_meanptop * coef
       isccp_meantb =  isccp_meantb * coef
       isccp_meantbclr =  isccp_meantbclr * coef
       isccp_meantaucld =  (6.82/((1./isccp_meanalbedocld)-1.))**(1./0.895)
     else
       fq_isccp = 0.
       isccp_totalcldarea = 0.
       isccp_lowcldarea = 0.
       isccp_midcldarea = 0.
       isccp_hghcldarea = 0.
       isccp_meantaucld = fillvalue
       isccp_meanptop = fillvalue
       isccp_meanalbedocld = fillvalue
       isccp_meantb = fillvalue
       isccp_meantbclr = fillvalue
     end if
   else if(nsun.gt.0.and.nsun.eq.nsun_max) then
     fq_isccp = 0.
     isccp_totalcldarea = 0.
     isccp_lowcldarea = 0.
     isccp_midcldarea = 0.
     isccp_hghcldarea = 0.
     isccp_meantaucld = fillvalue
     isccp_meanptop = fillvalue
     isccp_meanalbedocld = fillvalue
     isccp_meantb = fillvalue
     isccp_meantbclr = fillvalue
   else
     fq_isccp = fillvalue
     isccp_totalcldarea = fillvalue
     isccp_lowcldarea = fillvalue
     isccp_midcldarea = fillvalue
     isccp_hghcldarea = fillvalue
     isccp_meantaucld = fillvalue
     isccp_meanptop = fillvalue
     isccp_meanalbedocld = fillvalue
     isccp_meantb = fillvalue
     isccp_meantbclr = fillvalue
   end if
   
   s_acldisccp = isccp_totalcldarea
   s_acldlisccp = isccp_lowcldarea
   s_acldmisccp = isccp_midcldarea
   s_acldhisccp = isccp_hghcldarea
   s_tbisccp = isccp_meantb
   s_tbclrisccp = isccp_meantbclr
   s_cldtauisccp = isccp_meantaucld 
   s_cldalbisccp = isccp_meanalbedocld 
   s_ptopisccp = isccp_meanptop

   
   if(masterproc) then
     print*,'--------------- ISCCP Simulator Diagnostics: -----------------' 
     print*,'ISCCP PDF:'
     write(*,'(a,7(3x,f7.1))') '     ',isccpOpticalDepthBinEdges(1:ntau)
     do j=1,npres
       write(*,'(f6.0,7(3x,f7.3))') isccpPressureBinEdges(j),(fq_isccp(i,j),i=1,ntau)
     end do
     print*,'low cld area (tau > 0.3)=',isccp_lowcldarea
     print*,'mid cld area (tau > 0.3)=',isccp_midcldarea
     print*,'hgh cld area (tau > 0.3)=',isccp_hghcldarea
     print*,'total cld area (tau > 0.3)=',isccp_totalcldarea
     print*,'mean albedo cld (tau > 0.3) ',isccp_meanalbedocld
     print*,'mean tau cld (tau > 0.3) ',isccp_meantaucld
     print*,'mean top pressure (tau > 0.3)',isccp_meanptop
     print*,'mean brightness temp (tau > 0.3)',isccp_meantb
     print*,'mean clearsky brightness temp (tau > 0.3)',isccp_meantbclr

     if(dosimfilesout) then
     
       open (ntape,file='./OUT_STAT/'// &
                  trim(case)//'_'// &
                  trim(caseid)//'.isccp', &
                  status='unknown',form='unformatted')
       if(nstep.ne.nstat) then
         do while(.true.)
            read (ntape,end=222)
            read(ntape)  dummy,nsteplast
            if(nsteplast.ge.nstep) then
              backspace(ntape)
              backspace(ntape)  ! these two lines added because of
              read(ntape)       ! a bug in gfrotran compiler
              print*,'isccp file at nstep ',nsteplast
              goto 222
            end if
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
         end do
222      continue
         backspace(ntape)
         backspace(ntape)  ! these two lines added because of
         read(ntape)       ! a bug in gfrotran compiler
       endif

       print *,'Writting isccp file ',caseid
       write(ntape) caseid
       write(ntape) day,nstep
       write(ntape) ntau,npres
       write(ntape) isccpOpticalDepthBinEdges(1:ntau)
       write(ntape) isccpPressureBinEdges(1:npres)
       write(ntape) fq_isccp
       write(ntape) isccp_totalcldarea,isccp_lowcldarea,isccp_midcldarea,isccp_hghcldarea
       write(ntape) isccp_meantaucld,isccp_meanptop,isccp_meanalbedocld,isccp_meantb,isccp_meantbclr
       close(ntape)
     endif

   end if

end subroutine isccp_write

  ! -------------------------------------------------------------------------------

subroutine modis_write()
  
   use grid, only: nstep, nstat, day, caseid, nsubdomains, dompi, case, masterproc
   use vars, only: s_acldmodis, s_acldlmodis, s_acldmmodis, s_acldhmodis, &
                   s_relmodis, s_reimodis, s_lwpmodis, s_iwpmodis, &
                   s_acldliqmodis, s_acldicemodis, s_ptopmodis, &
                   s_cldtaumodis, s_cldtaulmodis, s_cldtauimodis
   implicit none     
   integer, parameter :: ntau = modis_num_taus , npres = modis_num_ctps
   integer, parameter :: buf_len = ntau*npres + 8
   real coef, buf(buf_len), buf1(buf_len), cldarea, cldarea_l, cldarea_i
   integer, parameter :: ntape = 88
   integer i,j,nmodis_max,nmodis_l_max, nmodis_i_max,nsun_max, itmp1(4), itmp2(4)
   real dummy
   integer nsteplast

   if(.not.domodis) return
                       
   nmodis_max = nmodis
   nmodis_l_max = nmodis_l
   nmodis_i_max = nmodis_i
   nsun_max = nsun
   if(dompi) then 
     itmp1(1) = nmodis
     itmp1(2) = nmodis_l
     itmp1(3) = nmodis_i
     itmp1(4) = nsun 
     call task_max_integer (itmp1,itmp2,4)
     nmodis_max = itmp2(1)
     nmodis_l_max = itmp2(2)
     nmodis_i_max = itmp2(3)
     nsun_max = itmp2(4)
   end if

   cldarea_l = modis_totalcldarea_l
   cldarea_i = modis_totalcldarea_i
    
   if(nmodis_max.gt.0.and.nsun.eq.nsun_max) then
     if(dompi) then
       if(nmodis.gt.0) then
         buf(1:ntau*npres) = reshape(fq_modis, shape = (/ ntau * npres /) )
         buf(ntau*npres+1) = modis_totalcldarea 
         buf(ntau*npres+2) = modis_meantaucld
         buf(ntau*npres+3) = modis_meanptop
         buf(ntau*npres+4) = modis_lowcldarea
         buf(ntau*npres+5) = modis_midcldarea
         buf(ntau*npres+6) = modis_hghcldarea
         buf(ntau*npres+7) = modis_totalcldarea_l 
         buf(ntau*npres+8) = modis_totalcldarea_i 
         coef = 1./float(nsubdomains*nmodis)
         buf = buf * coef
       else
         buf(1:buf_len) = 0.
       endif
       call task_sum_real(buf,buf1,buf_len)
       fq_modis(:,:) = reshape(buf1(1:ntau*npres),shape = (/ntau,npres/))
       modis_totalcldarea = buf1(ntau*npres+1)
       modis_meantaucld = buf1(ntau*npres+2)
       modis_meanptop = buf1(ntau*npres+3)
       modis_lowcldarea = buf1(ntau*npres+4)
       modis_midcldarea = buf1(ntau*npres+5)
       modis_hghcldarea = buf1(ntau*npres+6)
       modis_totalcldarea_l = buf1(ntau*npres+7)
       modis_totalcldarea_i = buf1(ntau*npres+8)
     else
       coef = 1./float(nmodis)
       fq_modis(:,:) = fq_modis(:,:) * coef
       modis_totalcldarea = modis_totalcldarea * coef
       modis_lowcldarea = modis_lowcldarea * coef
       modis_midcldarea = modis_midcldarea * coef
       modis_hghcldarea = modis_hghcldarea * coef
       modis_meanptop = modis_meanptop * coef
       modis_meantaucld = modis_meantaucld * coef
       modis_totalcldarea_l = modis_totalcldarea_l * coef
       modis_totalcldarea_i = modis_totalcldarea_i * coef
     end if
     if(modis_totalcldarea.gt.0.01) then
       coef = 1./(modis_totalcldarea)
       modis_meantaucld = modis_meantaucld * coef
       modis_meanptop =  modis_meanptop * coef
     else
       fq_modis = 0.
       modis_totalcldarea = 0.
       modis_lowcldarea = 0.
       modis_midcldarea = 0.
       modis_hghcldarea = 0.
       modis_totalcldarea_l = 0.
       modis_totalcldarea_i = 0.
       modis_meantaucld = fillvalue
       modis_meanptop = fillvalue
     end if
   else if(nsun.gt.0.and.nsun.eq.nsun_max) then
     fq_modis = 0.
     modis_totalcldarea = 0. 
     modis_lowcldarea = 0.
     modis_midcldarea = 0.
     modis_hghcldarea = 0.
     modis_totalcldarea_l = 0. 
     modis_totalcldarea_i = 0. 
     modis_meantaucld = fillvalue 
     modis_meanptop = fillvalue
   else
     fq_modis = fillvalue
     modis_totalcldarea = fillvalue
     modis_lowcldarea = fillvalue
     modis_midcldarea = fillvalue
     modis_hghcldarea = fillvalue
     modis_meantaucld = fillvalue
     modis_meanptop = fillvalue
     modis_totalcldarea_l = fillvalue
     modis_totalcldarea_i = fillvalue
   end if

   s_acldmodis = modis_totalcldarea
   s_acldlmodis = modis_lowcldarea
   s_acldmmodis = modis_midcldarea
   s_acldhmodis = modis_hghcldarea
   s_acldliqmodis = modis_totalcldarea_l
   s_acldicemodis = modis_totalcldarea_i
   s_cldtaumodis = modis_meantaucld  
   s_ptopmodis = modis_meanptop  

   if(nmodis_l_max.gt.0.and.nsun.eq.nsun_max) then
     if(dompi) then
       if(nmodis_l.gt.0) then
         buf(1) = cldarea_l
         buf(2) = modis_meantaucld_l
         buf(3) = modis_meanlwp
         buf(4) = modis_meanrel
         coef = 1./float(nsubdomains*nmodis_l)
         buf = buf * coef
       else
         buf(1:buf_len) = 0.
       endif
       call task_sum_real(buf,buf1,4)
       cldarea = buf1(1)
       modis_meantaucld_l = buf1(2)
       modis_meanlwp = buf1(3)
       modis_meanrel = buf1(4)
     else
       coef = 1./float(nmodis_l)
       cldarea = cldarea_l * coef
       modis_meantaucld_l = modis_meantaucld_l * coef
       modis_meanlwp = modis_meanlwp * coef
       modis_meanrel = modis_meanrel * coef
     end if
     if(cldarea.gt.0.01) then
       coef = 1./(cldarea)
       modis_meantaucld_l = modis_meantaucld_l * coef
       modis_meanlwp = modis_meanlwp * coef
       modis_meanrel = modis_meanrel * coef
     else
       modis_meantaucld_l = 0.
       modis_meanlwp = 0.
       modis_meanrel = 0.
     end if
   else if(nsun.gt.0.and.nsun.eq.nsun_max) then
     modis_meantaucld_l = fillvalue
     modis_meanlwp = fillvalue
     modis_meanrel = fillvalue
   else
     modis_meantaucld_l = fillvalue
     modis_meanlwp = fillvalue
     modis_meanrel = fillvalue
   end if

   if(nmodis_i_max.gt.0.and.nsun.eq.nsun_max) then
     if(dompi) then
       if(nmodis_i.gt.0) then
         buf(1) = cldarea_i
         buf(2) = modis_meantaucld_i
         buf(3) = modis_meaniwp
         buf(4) = modis_meanrei
         coef = 1./float(nsubdomains*nmodis_i)
         buf = buf * coef
       else
         buf(1:buf_len) = 0.
       endif
       call task_sum_real(buf,buf1,4)
       cldarea = buf1(1)
       modis_meantaucld_i = buf1(2)
       modis_meaniwp = buf1(3)
       modis_meanrei = buf1(4)
     else
       coef = 1./float(nmodis_i)
       cldarea = cldarea_i * coef
       modis_meantaucld_i = modis_meantaucld_i * coef
       modis_meaniwp = modis_meaniwp * coef
       modis_meanrei = modis_meanrei * coef
     end if
     if(cldarea.gt.0.01) then
       coef = 1./(cldarea)
       modis_meantaucld_i = modis_meantaucld_i * coef
       modis_meaniwp = modis_meaniwp * coef
       modis_meanrei = modis_meanrei * coef
     else
       modis_meantaucld_i = 0.
       modis_meaniwp = 0.
       modis_meanrei = 0.
     end if
   else if(nsun.gt.0.and.nsun.eq.nsun_max) then
     modis_meantaucld_i = fillvalue
     modis_meaniwp = fillvalue
     modis_meanrei = fillvalue
   else
     modis_meantaucld_i = fillvalue
     modis_meaniwp = fillvalue
     modis_meanrei = fillvalue
   end if

   s_relmodis = modis_meanrel
   s_reimodis = modis_meanrei
   s_lwpmodis = modis_meanlwp  
   s_iwpmodis = modis_meaniwp  
   s_cldtaulmodis = modis_meantaucld_l  
   s_cldtauimodis = modis_meantaucld_i  
   
   if(masterproc) then      
     print*,'--------------- MODIS Simulator Diagnostics: -----------------' 
     print*,'MODIS PDF:' 
     write(*,'(a,10x,6(3x,f7.1))') '     ',isccpOpticalDepthBinEdges(2:ntau+1)
     do j=1,npres
       write(*,'(f6.0,10x,6(3x,f7.3))') isccpPressureBinEdges(j),(fq_modis(i,j),i=1,ntau)
     end do
     print*,'low cld area =',modis_lowcldarea
     print*,'mid cld area =',modis_midcldarea
     print*,'hgh cld area =',modis_hghcldarea
     print*,'total cld area =',modis_totalcldarea
     print*,'mean tau cld ',modis_meantaucld
     print*,'mean top pressure ',modis_meanptop
     print*,'total cld area (liquid) =',modis_totalcldarea_l
     print*,'total cld area (ice) =',modis_totalcldarea_i
     print*,'mean tau (liquid) =',modis_meantaucld_l
     print*,'mean tau (ice) =',modis_meantaucld_i
     print*,'mean lwp (g/m2) =',modis_meanlwp
     print*,'mean iwp (g/m2) =',modis_meaniwp
     print*,'mean reff (liquid)  =',modis_meanrel
     print*,'mean reff (ice)  =',modis_meanrei
       
     if(dosimfilesout) then
     
       open (ntape,file='./OUT_STAT/'// &
                  trim(case)//'_'// &
                  trim(caseid)//'.modis', &
                  status='unknown',form='unformatted')
       if(nstep.ne.nstat) then
         do while(.true.)
            read (ntape,end=222)
            read(ntape)  dummy,nsteplast
            if(nsteplast.ge.nstep) then
              backspace(ntape)
              backspace(ntape)  ! these two lines added because of
              read(ntape)       ! a bug in gfrotran compiler
              print*,'imodis file at nstep ',nsteplast
              goto 222
            end if
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
         end do
222      continue
         backspace(ntape)
         backspace(ntape)  ! these two lines added because of
         read(ntape)       ! a bug in gfrotran compiler
       endif
     
       print *,'Writting modis file ',caseid
       write(ntape) caseid
       write(ntape) day,nstep
       write(ntape)  ntau,npres
       write(ntape) isccpOpticalDepthBinEdges(1:ntau)
       write(ntape) isccpPressureBinEdges(1:npres)
       write(ntape) fq_modis
       write(ntape) modis_totalcldarea,modis_lowcldarea,modis_midcldarea,modis_hghcldarea
       write(ntape) modis_meantaucld,modis_meanptop
       write(ntape) modis_totalcldarea_l,modis_meantaucld_l, modis_meanlwp, modis_meanrel
       write(ntape) modis_totalcldarea_i,modis_meantaucld_i, modis_meaniwp, modis_meanrei
       close(ntape)
     endif 

   end if

end subroutine modis_write
       


subroutine misr_write()

   use grid, only: nstep, nstat, day, caseid, nsubdomains, dompi, case, masterproc
   use vars, only: s_acldmisr, s_ztopmisr
   implicit none
   integer, parameter :: ntau = isccp_num_taus , npres = misr_num_cths
   integer, parameter :: buf_len = ntau*npres + 2
   real coef, buf(buf_len), buf1(buf_len)
   integer, parameter :: ntape = 88
   integer i,j,nmisr_max,nsun_max, itmp1(2), itmp2(2)
   real dummy
   integer nsteplast

   if(.not.domisr) return

   nmisr_max = nmisr
   nsun_max = nsun
   if(dompi) then
     itmp1(1) = nmisr
     itmp1(2) = nsun
     call task_max_integer (itmp1,itmp2,2)
     nmisr_max = itmp2(1)
     nsun_max = itmp2(2)
   end if

   if(nmisr_max.gt.0.and.nsun.eq.nsun_max) then
     if(dompi) then
       if(nmisr.gt.0) then
         buf(1:ntau*npres) = reshape(fq_misr, shape = (/ ntau * npres /) )
         buf(ntau*npres+1) = misr_totalcldarea
         buf(ntau*npres+2) = misr_meanztop
         coef = 1./float(nsubdomains*nmisr)
         buf = buf * coef
       else
         buf(1:buf_len) = 0.
       endif
       call task_sum_real(buf,buf1,buf_len)
       fq_misr(:,:) = reshape(buf1(1:ntau*npres),shape = (/ntau,npres/))
       misr_totalcldarea = buf1(ntau*npres+1)
       misr_meanztop = buf1(ntau*npres+2)
     else
       coef = 1./float(nmisr)
       fq_misr(:,:) = fq_misr(:,:) * coef
       misr_totalcldarea = misr_totalcldarea * coef
       misr_meanztop = misr_meanztop * coef
     end if
     if(misr_totalcldarea.gt.0.01) then
       coef = 1./(misr_totalcldarea)
       misr_meanztop = misr_meanztop * coef
     else
       fq_misr = 0.
       misr_totalcldarea = 0.
       misr_meanztop = 0.
     end if
   else if(nsun.gt.0.and.nsun.eq.nsun_max) then
     fq_misr = 0.
     misr_totalcldarea = 0.
     misr_meanztop = fillvalue
   else
     fq_misr(:,:) = fillvalue
     misr_totalcldarea = fillvalue
     misr_meanztop = fillvalue
   end if

   s_acldmisr = misr_totalcldarea
   s_ztopmisr = misr_meanztop

   if(masterproc) then
     print*,'--------------- MISR Simulator Diagnostics: -----------------'
     print*,'MISR PDF:'
     write(*,'(a,7(3x,f7.1))') '     ',isccpOpticalDepthBinEdges(1:ntau)
     do j=npres,1,-1
       write(*,'(f6.2,7(3x,f7.3))') MISR_CTH_boundaries(j),(fq_misr(i,j),i=1,ntau)
     end do
     print*,'total cld area =',misr_totalcldarea
     print*,'mean top height (m)',misr_meanztop

     if(dosimfilesout) then

       open (ntape,file='./OUT_STAT/'// &
                  trim(case)//'_'// &
                  trim(caseid)//'.misr', &
                  status='unknown',form='unformatted')
       if(nstep.ne.nstat) then
         do while(.true.)
            read (ntape,end=222)
            read(ntape)  dummy,nsteplast
            if(nsteplast.ge.nstep) then
              backspace(ntape)
              backspace(ntape)  ! these two lines added because of
              read(ntape)       ! a bug in gfrotran compiler
              print*,'misr file at nstep ',nsteplast
              goto 222
            end if
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
            read(ntape)
         end do
222      continue
         backspace(ntape)
         backspace(ntape)  ! these two lines added because of
         read(ntape)       ! a bug in gfrotran compiler
       endif


       print *,'Writting misr file ',caseid
       write(ntape) caseid
       write(ntape)  day,nstep
       write(ntape)  ntau,npres
       write(ntape) isccpOpticalDepthBinEdges(1:ntau)
       write(ntape) MISR_CTH_boundaries(1:npres)
       write(ntape) fq_misr
       write(ntape) misr_totalcldarea,misr_meanztop
       close(ntape)
     endif
   end if

end subroutine misr_write



end module instrument_diagnostics 

