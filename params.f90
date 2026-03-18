module params

use grid, only: nzm

implicit none

!   Constants:

real, parameter :: cp = 1004.             ! Specific heat of air, J/kg/K
real, parameter :: ggr = 9.81             ! Gravity acceleration, m/s2
real, parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
real, parameter :: lfus = 0.3336e+06      ! Latent heat of fusion, J/kg
real, parameter :: lsub = 2.8440e+06      ! Latent heat of sublimation, J/kg
real, parameter :: rv = 461.              ! Gas constant for water vapor, J/kg/K
real, parameter :: rgas = 287.            ! Gas constant for dry air, J/kg/K
real, parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
real, parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
real, parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

real, parameter :: fac_cond = lcond/cp 
real, parameter :: fac_fus = lfus/cp
real, parameter :: fac_sub = lsub/cp

real, parameter ::  pi = 3.141592653589793

!
!----------------------------------------------
! internally set parameters:

real   epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
logical:: dosubsidence = .false.
real fcorz      ! Vertical Coriolis parameter
real :: coszrs = 0.
real salt_factor ! correction factor for water vapor saturation over sea-water

integer:: ncycle_max = 50  ! maximum number of subcycling within dt

!----------------------------------------------
! Parameters set by PARAMETERS namelist:
! Initialized to default values.
!----------------------------------------------

real:: ug = 0.        ! Velocity of the Domain's drift in x direction
real:: vg	= 0.        ! Velocity of the Domain's drift in y direction
real:: fcor = -999.   ! Coriolis parameter	
real:: longitude0 = 0.    ! latitude of the domain's center 
real:: latitude0  = 0.    ! longitude of the domain's center 
real:: nxco2 = 1         ! factor to modify co2 concentration
logical:: doradlat = .false.
logical:: doradlon = .false.

real(8):: tabs_s =0.	! surface temperature,K
real:: delta_sst = 0.   ! amplitude of sin-pattern of sst about tabs_s (ocean_type=1)
real:: depth_slab_ocean = 2. ! thickness of the slab-ocean (m)
real:: Szero = 0.  ! mean ocean transport (W/m2)
real:: deltaS = 0. ! amplitude of linear variation of ocean transport (W/m2)
real:: timesimpleocean = 0. ! time to start simple ocean

real::   fluxt0 =0.  ! surface sensible flux, Km/s
real::   fluxq0 =0.  ! surface latent flux, m/s
real::   tau0   =0.  ! surface stress, m2/s2
real::   z0     =0.035	! roughness length
real::   soil_wetness =1.! wetness coeff for soil (from 0 to 1.)
integer:: ocean_type =0 ! type of SST forcing
logical:: OCEAN =.false.  ! flag indicating that surface is water
logical:: LAND =.false.   ! flag indicating that surface is land
logical:: SFC_FLX_FXD =.false. ! surface sensible flux is fixed
logical:: SFC_TAU_FXD =.false.! surface drag is fixed
logical:: SLM = .false. ! interactive land model

logical:: LES_S = .true. ! .true. sample cloud and core statistics as for PBL clouds 
                         ! .false. - as for deep clouds (ql>0.01qsat, w>1 m/s, etc)
real:: timelargescale =0. ! time to start large-scale forcing

! nudging boundaries (between z1 and z2, where z2 > z1): 
real:: nudging_uv_z1 =-1., nudging_uv_z2 = 1000000.
real:: nudging_t_z1 =-1., nudging_t_z2 = 1000000.
real:: nudging_q_z1 =-1., nudging_q_z2 = 1000000.
real:: tauls = 99999999.    ! nudging-to-large-scaler-profile time-scale
real:: tautqls = 99999999.! nudging-to-large-scaler-profile time-scale for scalars

logical:: dodamping = .false.
logical:: doupperbound = .false. 
logical:: docloud = .false. 
logical:: doprecip = .false.
logical:: dolongwave = .false. 
logical:: doshortwave = .false.
logical:: dosgs = .false.
logical:: docoriolis = .false. 
logical:: docoriolisz = .false. 
logical:: dofplane = .true.
logical:: dosurface = .false. 
logical:: dolargescale = .false. 
logical:: doradforcing = .false.
logical:: dosfcforcing = .false. 
logical:: doradsimple = .false. 
logical:: donudging_uv = .false. 
logical:: donudging_tq = .false.
logical:: donudging_t = .false. 
logical:: donudging_q = .false.
logical:: doensemble = .false. 
logical:: dowallx = .false. 
logical:: dowally = .false. 
logical:: docolumn = .false. 
logical:: docup = .false.
logical:: doperpetual = .false. 
logical:: doseasons = .false. 
logical:: doradhomo = .false. 
logical:: dosfchomo = .false.
logical:: dossthomo = .false. 
logical:: dodynamicocean = .false. 
logical:: dosolarconstant = .false.
logical:: dotracers = .false. 
logical:: dosmoke = .false. 
logical:: notracegases = .false.
logical:: doseawater = .false.

! Specify solar constant and zenith angle for perpetual insolation.
! Based onn Tompkins and Graig (1998)
! Note that if doperpetual=.true. and dosolarconstant=.false.
! the insolation will be set to the daily-averaged value on day0.
real:: solar_constant = 685. ! solar constant (in W/m2)
real:: zenith_angle = 51.7   ! zenith angle (in degrees)

integer:: nensemble =0   ! the number of subensemble set of perturbations
integer:: perturb_type  = 0 ! type of initial noise in setperturb()

! Initial bubble parameters. Activated when perturb_type = 2
real:: bubble_x0 = 0.
real:: bubble_y0 = 0.
real:: bubble_z0 = 0.
real:: bubble_radius_hor = 0.
real:: bubble_radius_ver = 0.
real:: bubble_dtemp = 0.
real:: bubble_dq = 0.

!=====================================================
! Kuang-Lab Additions Begin Here

! Options
logical:: dompiensemble = .false. ! Subdomains defined in domains.f90 are run separately

! linear response perturbation: layer by layer (Song Qiyu, 2022)
logical:: dolayerperturb = .false.
integer:: tperturbi = 0
integer:: qperturbi = 0
real:: tperturbA = 1.     ! Default perturbation 1 time positive
real:: qperturbA = 1.

! Radiative tendencies as per Pauluis & Garner [2006]
! Added by Nathanael Wong on 2023/07/05
logical :: doradtendency = .false. 
real :: troptend = 1.5 ! Convective tendency in Pauluis & Garner [2006]

! Option to fix wind speed used in calculation of bulk surface fluxes
! Taken from Peter Blossey's version of SAM
! Added by Nathanael Wong on 2023/07/05
logical :: dobulksfcflx = .false.
real :: bulksfcflx_u = 5.

! Option to divide output stat file every some steps
! to avoid large stat file
logical :: dosepstat = .false.
integer :: nstep_sepstat = 99999999

! Damped Gravity Wave and Temperature Gradient Relaxation Implementations
! Added by Nathanael Wong on 2023/07/05
logical :: dodgw = .false.
logical :: dotgr = .false.
logical :: dowtg_decomp = .false.
real :: wtgscale_time = 0. ! period over which theta relaxation timescale scales from infinity to ttheta_wtg.  Express as fraction of time over which WTG large-scale forcing is implemented.  So if WTG/Large-scale is turned on for 100 days, twtg_scale = 1/4 means that the scaling up to WTG occurs over 25 days.

logical :: dowtg_timedependence = .false.
! Add time-dependence to WTG schemes, which could be important for small damping values

logical :: dowtg_blossey_etal_JAMES2009  = .false.
logical :: dowtg_raymondzeng_QJRMS2005   = .false. 
logical :: dowtg_hermanraymond_JAMES2014 = .false.
logical :: dowtg_kuang_JAS2008           = .false.
logical :: dowtg_decompdgw  = .false.
logical :: dowtg_decomptgr  = .false.

real :: am_wtg = 1. ! momentum damping rate in 1/d -- note must be non-zero.
real :: am_wtg_exp = 0. ! exponent of p/p0 in momentum damping rate.
real :: lambda_wtg = 650.e3 ! quarter wavelength in m. default = 650.e3 (=650 km).

real :: tau_wtg = 1. ! Relaxation timescale (in hours) for WTG Approximation of Raymond and Zeng [2005]
logical :: dowtgLBL = .false.
logical :: boundstatic = .true. ! Restrict the static stability lower bound to prevent unrealistically large values of w_wtg
real :: dthetadz_min = 1.e-3 ! if boundstatic = .true., what is the minimum bound? Default from Raymond & Zeng [2005] is 1.e-3 K/km
real :: wtgscale_vertmodepwr = 1. ! Spectral decomposition power, default is 1 as per Herman and Raymond [2014]

integer :: wtgscale_vertmodenum = 2! number of vertical modes
real, dimension(2) :: wtgscale_vertmodescl = (/1., 1./) ! strength scaling for vertical modes (number of items = wtgscale_vertmodenum)

! Specify a "island" within which SST is allowed to vary
! If dosstisland = .false. and  dodynamicocean = .true. the entire domain SST varies
logical :: dosstislands = .false. ! specify an island within which SST is allowed to vary
real    :: sstislands_oceanmld = 0. ! "ocean" slab depth, if 0, ocean SST is constant
real    :: sstislands_landmld  = 0. ! "island" slab depth, set to depth_slab_ocean if 0

! Specify round islands using formula. If readlsm = true, then lsm file will override this
real    :: sstislands_radius   = 0. ! "island" radii in meters
integer :: sstislands_nrow = 1.  ! number of island rows
integer :: sstislands_ncol = 1.  ! number of island columns
real    :: sstislands_sep  = 0.  ! spacing between island centers, should be at least 2*sstisland_radius

! Alternatively, specify a file to read the land-sea mask data
! File is a binary file, with variables in this order:
! (1) nx_lsm, which is an integer specifying number of x points
! (2) ny_lsm, which is an integer specifying number of y points
! (3) lsm, which is an array of 1s and 0s, with 0s denoting ocean and 1s denoting land
logical       :: readlsm = .false. ! read land-sea mask from file
character(80) :: lsmfile = ""

! If nrestart = 2 and dodynamicocean = false, if nrestart_resetsst = true, set all sst back to tabs_s
logical :: nrestart_resetsst = .false.

! Try to simulate a Hadley Cell (with varying tropopause heights)
logical :: dohadley = .false.
!logical :: dodrivenequilibrium = .false.
real    :: whadmax  = 0.
real    :: zhadmax  = 15000.
real    :: hadscale_time = 0.

! Option for whether to advect 1d profile for wtg advection
logical :: doadv3d = .true. ! this was true in the vanilla SAM

! Options for which 1d profile to advect if advecting 1d
logical :: doadvinic = .true. ! advect initial profile
logical :: docalcwtgbg = .false. ! even if not doadvbg, may still want to calculate the wtg bg
logical :: doadvbg = .false. ! advect computed background (time-invariant)
logical :: doadvensnoise = .false. ! works with dompiensemble, use the large-scale w from one ensemble member
! if none of the above: advect the mean profile, which will also evolve with time

! temporal parameters for wtg
integer :: nstartwtg = 0 ! nstep to start wtg calculation
integer :: nstepwtgbg = 0 ! nstep after nstartwtg to calculate wtg background profile
integer :: nstepwtg = 999999999 ! keep constant wtg forcing in this many steps
logical :: dowtgtimestep = .false.

! Options for large-scale vertical advection of temperature/moisture, u/v wind
logical :: dotqlsvadv = .true.
logical :: douvlsvadv = .true.

! Kuang-Lab Additions End Here
!=====================================================

end module params
