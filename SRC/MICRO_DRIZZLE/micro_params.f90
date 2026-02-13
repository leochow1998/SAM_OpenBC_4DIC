module micro_params

use grid, only: nzm

implicit none

!  Microphysics stuff:

! Densities of hydrometeors

real, parameter :: rhor = 1000. ! Density of water, kg/m3

real, parameter :: qp_threshold = 1.e-8 ! minimal rain/snow water content

real, parameter :: rd_min = 25.e-6 ! Minimum drizzle drop size
real, parameter :: rd_max = 200.e-6 ! Maximum drizzle volume radius
real, parameter :: coefrv = 3./(4.*3.1415*1000.) ! coef to compute volume radius: rv = (coefrv * qr * N)^1/3
real, parameter :: coefconpmin = coefrv/rd_max**3
real, parameter :: coefconpmax = coefrv/rd_min**3
real :: Nc0 ! (cm-3) Prescribed cloud drop concentration 

!bloss: Nc0 can now be set in the MICRO_DRIZZLE namelist
!bloss  real, parameter ::  Nc0 = 40. ! (cm-3) Prescribed cloud drop concentration 

real evapr1(nzm),evapr2(nzm)

end module micro_params
