!> \file
!! Contains CLASSIC globally accessible parameters
!! @author J. Melton
!!
module classicParams

  !> \ingroup classicparams_main
  !! @{

  ! J. Melton
  ! Jun 23 2013

  ! Change History:

  ! April 2021 - GM - Add in new parameters for dry surface layer (DSL) parameterization (DSLinit, DELTAZ, PSIAIR)
  ! Sep 2018 - EC - Relocate variables defined in "readInParams" to the main module and allocate them
  !                 in "allocateParamsCTEM". Cray compiler doesn't currently support namelists using
  !                 dynamic arrays dimensioned by a parameter specified in the main module.
  ! Jun 2017 - JM - Change to using namelist.
  ! Jun 11 2015 - JM - Add in new disturb params (extnmois_duff and extnmois_veg) replacing duff_dry and extnmois.
  ! Mar 12 2014 - JM - Allow for two sets of paramters (competition and prescribed PFT fractions).
  !                    all ctem subroutines except photosynCanopyConduct keep their parameters here.
  !
  ! Jan 17 2014 - JM - Add in more parameters from ctem.f, phenology.f, and allocate.f.

  implicit none

  public :: readInPFTNums
  public :: readInParams
  public :: prepareGlobalParams
  private :: findPFTindexes

  ! Constants

  real, parameter :: zero     = 1.0e-20                  !< Defintion of zero
  real, parameter :: abszero  = 1e-12                    !< Defintion of zero, this one is for allocateCarbon.f90 and applyAllometry.f90
  real, parameter :: epszero  = max(3*epsilon(0.),abszero) !< Precision dependent definition of zero.
                                                         !! FLAG: EC - Only in allocateCarbon.f90 and applyAllometry.f90, but consider 
                                                         !!            using for comparisons of real residuals against zero.

  real, parameter :: earthrad = 6371.22   !< Radius of Earth, km
  real, parameter :: km2tom2  = 1.0e+06   !< Conversion factor from \f$km^2\f$ to \f$m^2\f$
  real, parameter :: gasc = 8.314         !< Gas constant (\f$J mol^{-1} K^{-1}\f$)
  real, parameter :: convertkgC = 1.201e-8 !< Converts from \f$ \mu molCO2/m2/s\f$ to \f$kgC/m2/s\f$
  real, parameter :: convertkgN = 1.15740740741e-08 !< Converts from \f$ g/(m2.day) \f$ to \f$ (kg/m2.s) \f$
  real, parameter :: convertg2kg = 0.001
  ! The next parameters are normally assigned in the GCM but need to be assigned here since this
  ! is an offline run. The names used below are the same as those in the GCM.
  real, parameter :: TFREZ = 273.16       !< Freezing point of water (K) (GCM name: CELZRO)
  real, parameter :: RGAS = 287.04        !< Gas constant (\f$J kg^{-1} K^{-1}\f$) (GCM name: GAS)
  real, parameter :: RGASV = 461.50       !< Gas constant for water vapour (\f$J kg^{-1} K^{-1}\f$) (GCM name: GASV)
  real, parameter :: GRAV = 9.80616       !< Acceleration due to gravity (\f$m s^{-1}\f$) (GCM name: G)
  real, parameter :: SBC = 5.66796E-8     !< Stefan-Boltzmann constant (\f$W m^{-2} K^{-4} \f$) (GCM name: SIGMA)
  real, parameter :: SPHAIR = 1.00464E3   !< Specific heat of air (\f$J kg^{-1} K^{-1}\f$) (GCM name: CPRES)
  real, parameter :: PI = 3.1415926535898 !< pi (-) (GCM name: CPI)
  real, parameter :: STD_PRESS = 101325.0 !< Standard atmospheric pressure (Pa)

  real, parameter :: lambda14C = 8267.  !< radioactive decay rate for 14C, corresponds to a half life of 5730 years. (\f$yr^{-1}\f$)

  ! FLAG below not used anywhere.
  ! AI=2.88053E+6/1004.5
  ! BI=0.167E+3/1004.5
  ! AW=3.15213E+6/1004.5
  ! BW=2.38E+3/1004.5
  ! SLP=1.0/(T1S-T2S)

  ! These month arrays can be overwritten during leap years.
  integer, dimension(12) :: monthdays !< days in each month
  integer, dimension(13) :: monthend  !< calender day at end of each month
  integer, dimension(12) :: mmday     !< mid-month day
  integer, parameter :: nmon = 12 !< Number of months in a year

  ! Additional values for RPN and GCM common blocks:

  real, parameter :: DELTA = 0.608
  real, parameter :: AS = 12.0
  real, parameter :: ASX = 4.7
  real, parameter :: ANGMAX = 0.85
  real, parameter :: CI = 40.0
  real, parameter :: BS = 1.0
  real, parameter :: BETA = 1.0
  real, parameter :: FACTN = 1.2
  real, parameter :: HMIN = 40.0
  real, parameter :: A = 21.656
  real, parameter :: B = 5418.0
  real, parameter :: EPS1 = 0.622
  real, parameter :: EPS2 = 0.378
  real, parameter :: T1S = 273.16
  real, parameter :: T2S = 233.16
  real, parameter :: RW1 = 53.67957
  real, parameter :: RW2 = - 6743.769
  real, parameter :: RW3 = - 4.8451
  real, parameter :: RI1 = 23.33086
  real, parameter :: RI2 = - 6111.72784
  real, parameter :: RI3 = 0.15215
  real, parameter :: X1 = 0.0
  real, parameter :: X2 = 0.0
  real, parameter :: X3 = 0.0
  real, parameter :: X4 = 0.0
  real, parameter :: X5 = 0.0
  real, parameter :: X6 = 0.0
  real, parameter :: X7 = 0.0
  real, parameter :: X8 = 0.0
  real, parameter :: X9 = 0.0
  real, parameter :: X10 = 0.0
  real, parameter :: X11 = 0.0
  real, parameter :: X12 = 0.0
  real, parameter :: X13 = 0.0
  real, parameter :: X14 = 0.0
  real, parameter :: X15 = 0.0
  real, parameter :: X16 = 0.0
  real, parameter :: K1 = 0
  real, parameter :: K2 = 0
  real, parameter :: K3 = 0
  real, parameter :: K4 = 0
  real, parameter :: K5 = 0
  real, parameter :: K6 = 0
  real, parameter :: K7 = 0
  real, parameter :: K8 = 0
  real, parameter :: K9 = 0
  real, parameter :: K10 = 0
  real, parameter :: K11 = 0

  ! Biogeochemical parameters:
  real, parameter :: deltat   = 1.0       !< CTEM's time step in days
  real, parameter :: seed    = 0.001    !< seed pft fraction, same as in competition in mosaic mode, all tiles are given this as a minimum
  real, parameter :: minbare = 1.0e-5   !< minimum bare fraction when running competition on to prevent numerical problems.
  real, parameter :: c2dom   = 450.0    !< gc / kg dry organic matter conversion factor from carbon to dry organic matter value is from Li et al. 2012 biogeosci
  real, parameter :: wtCH4   = 16.044   !< Molar mass of CH4 (\f$g mol^{-1}\f$)

  integer, parameter :: nbs = 4         !< Number of modelled shortwave radiation wavelength bands COMBAK Can be read in from the init file when I have the new snow albedo scheme fully implemented. Leave here for now.

  integer, parameter :: soilcolrinds = 20 !< Number of soil colour index classes used (Affects ALWV, ALDV, ALWN, ALDN)

  ! Vivek increased tolrance  from 0.0001 to 0.00025 so that litter imbalance (now that we have litter in multiple soil layers)
  ! passes, Dec 2020
  !real, parameter  :: tolrance = 0.00025 !< our tolerance for balancing c budget in kg c/m2 in one day (differs when competition on or not)
  real, parameter  :: tolrance = 0.0005

  ! Relative tolerances to work with 32-bit mode. It seems the tolerance of 0.5 Kg C for a whole grid cell
  ! is too much for 32 bit to preserve. In addition, we need to a use a similar tolerance for denstities as wel.
  real, parameter :: rel_tolerance_amount = 1.0E-03 
  real, parameter :: rel_tolerance_density = 1.0E-06 

  ! YW May 12, 2015 in peatland the C balance gap reaches 0.00016.
  real, parameter  :: tolrnce1 = 0.5      !< kg c, tolerance of total c balance (FOR LUC)

  !> Logical switch for using constant allocation factors (default value is false)
  logical :: consallo = .false.
  
  ! --------------------------------------------------------------

  ! Originally read from netcdf initialization now made static
#if defined without_agcm_
  integer, parameter :: nlat = 1        !< Number of cells we are running, set in the code. Offline this always 1.
#else
  integer, parameter :: nlat = 131      !< Size of latitude dimension when coupled to CanESM.
#endif
  integer, parameter :: nmos = 1           !< Number of mosaic tiles, set in the code
  integer, parameter :: ilg = nlat * nmos  !< nlat x nmos
  integer, parameter :: ignd = 20          !< Number of soil layers, set in the code

  ! Originally read in from the run parameters file now made static
#if defined without_agcm_  
  ! Presently the offline benchmarking suite (CBC) is set up with 12 PFTs. If you wish to 
  ! run offline with 9 PFTs (or any other variation), change the numbers below.
  integer, parameter :: ican = 5             !< Number of CLASS (physics) pfts, set in the code
  integer, parameter :: icc = 12             !< Number of CTEM (biogeochemical) pfts, set in the code
#else
  ! CanESM runs with a default of 9 PFTs
  integer, parameter :: ican = 4             !< Number of CLASS (physics) pfts, set in the code
  integer, parameter :: icc = 9              !< Number of CTEM (biogeochemical) pfts, set in the code
#endif
  integer, parameter :: l2max = 3            !< Maximum number of level 2 CTEM PFTs, set in the code

  ! These are used for the model checks for these parameters
  integer :: testnmos        !< Used to test the number of mosaic tiles, read in from the initialization file
  integer :: testignd        !< Used to test the number of soil layers, read in from the initialization file
  integer :: testican        !< Used to test the number of CLASS (physics) pfts, read in from the parameters namelist
  integer :: test2ican       !< Used to test the number of CLASS (physics) pfts a second time, read from the initialization file  
  integer :: testicc         !< Used to test the number of CTEM (biogeochemical) pfts, read in from the parameters namelist
  integer :: test2icc        !< Used to test the number of CTEM (biogeochemical) pfts a second time, read in from the initialization file
  integer :: testl2max       !< Used to test the maximum number of level 2 CTEM PFTs. This is the maximum number of CTEM PFTs, read in from the parameters namelist

  ! ----
  ! Orginally calculated plant related parameters, now made static
  integer, parameter :: icp1 = ican + 1
  integer, parameter :: iccp1 = icc + 1
  integer, parameter :: iccp2 = icc + 2

  integer, parameter :: kk = l2max * ican !< product of class pfts and l2max
  integer :: numcrops          !< number of crop pfts
  integer :: numtreepfts       !< number of tree pfts
  integer :: numgrass          !< number of grass pfts
  integer :: numshrubs         !< number of shrubs pfts
  integer, dimension(ican) :: nol2pfts !< Number of level 2 PFTs calculated in readInParams
  logical, dimension(icc) :: crop     !< simple crop matrix, define number and position of the crops (NOTE: dimension icc)
  logical, dimension(icc) :: grass    !< simple grass matric, define the number and position of grass (NOTE: dimension icc)
  integer, dimension(kk) :: CL4CTEM  !< Indexing of the CTEM-level PFTs into a CLASS PFT-level array.
  integer, dimension(ican,2) :: reindexPFTs !< Reindexing arrays of CLASS variables into the parameters arrays at the CTEM level.

  ! Four-band albedo information
  integer, PARAMETER :: nsmu     = 10
  integer, PARAMETER :: nsalb    = 11
  integer, PARAMETER :: nbc      = 20
  integer, PARAMETER :: nreff    = 10
  integer, PARAMETER :: nswe     = 11
  integer, PARAMETER :: nbnd_lut = 4
  ! The following are look up tables that are filled in modelStateDrivers:getInput.
  real, DIMENSION(NBC,NSWE,NREFF,NSMU,NSALB,NBND_LUT) :: albdif_lut
  real, DIMENSION(NBC,NSWE,NREFF,NSMU,NSALB,NBND_LUT) :: albdir_lut
  real, DIMENSION(NBC,NSWE,NREFF,NSMU,NSALB,NBND_LUT) :: trandif_lut
  real, DIMENSION(NBC,NSWE,NREFF,NSMU,NSALB,NBND_LUT) :: trandir_lut

  ! ============================================================
  ! Read in from the namelist: ---------------------------------

  ! Vegetation related parameters that define what PFTs are being simulated:
  character(5), dimension(ican) :: classpfts !< CLASS (physics) PFTs
  integer, dimension(kk)      :: modelpft  !< Separation of pfts into level 1 (for class) and level 2 (for ctem) pfts.
  character(8), dimension(kk) :: ctempfts  !< List of CTEM-level PFTs

  ! Physics (CLASS) parameters:

  real :: VKC    !< Von Karman Constant (-)
  real :: CT     !< Drag Coefficient for water (-)
  real :: VMIN   !< Minimum wind speed (\f$m s^{-1}\f$)
  real :: TCW    !< Thermal conductivity of water (\f$W m^{-1} K^{-1} \f$)
  real :: TCICE  !< Thermal conductivity of ice  (\f$W m^{-1} K^{-1} \f$)
  real :: TCSAND !< Thermal conductivity of sand particles (\f$W m^{-1} K^{-1} \f$)
  real :: TCCLAY !< Thermal conductivity of fine mineral particles  (\f$W m^{-1} K^{-1} \f$)
  real :: TCOM   !< Thermal conductivity of organic matter (\f$W m^{-1} K^{-1} \f$)
  real :: TCDRYS !< Thermal conductivity of dry mineral soil  (\f$W m^{-1} K^{-1} \f$) FLAG QUESTION ever used?
  real :: RHOSOL !< Particle density of soil mineral matter (\f$kg m^{-3}\f$)
  real :: RHOOM  !< Particle density of soil organic matter (\f$kg m^{-3}\f$)
  real :: HCPW   !< Volumetric heat capacity of water (\f$J m^{-3} K^{-1}\f$)
  real :: HCPICE !< Volumetric heat capacity of ice (\f$J m^{-3} K^{-1}\f$)
  real :: HCPSOL !< Volumetric heat capacity of mineral matter (\f$J m^{-3} K^{-1}\f$) FLAG QUESTION ever used?
  real :: HCPOM  !< Volumetric heat capacity of organic matter (\f$J m^{-3} K^{-1}\f$)
  real :: HCPSND !< Volumetric heat capacity of sand particles (\f$J m^{-3} K^{-1}\f$)
  real :: HCPCLY !< Volumetric heat capacity of fine mineral particles (\f$J m^{-3} K^{-1}\f$)
  real :: SPHW   !< Specific heat of water (\f$J kg^{-1} K^{-1}\f$)
  real :: SPHICE !< Specific heat of ice  (\f$J kg^{-1} K^{-1}\f$)
  real :: SPHVEG !< Specific heat of vegetation matter  (\f$J kg^{-1} K^{-1}\f$)
  real :: RHOW   !< Density of water (\f$kg m^{-3}\f$)
  real :: RHOICE !< Density of ice (\f$kg m^{-3}\f$)
  real :: TCGLAC !< Thermal conductivity of ice sheets (\f$W m^{-1} K^{-1}\f$)
  real :: CLHMLT !< Latent heat of freezing of water (\f$J kg^{-1}\f$)
  real :: CLHVAP !< Latent heat of vaporization of water (\f$J kg^{-1}\f$)
  real :: ZOLNG  !< Natural log of roughness length of soil (-)
  real :: ZOLNS  !< Natural log of roughness length of snow (-)
  real :: ZOLNI  !< Natural log of roughness length of ice (-)
  real :: ZORATG !< Ratio of soil roughness length for momentum to roughness length for heat (-)
  real :: ALVSI  !< Visible albedo of ice (-)
  real :: ALIRI  !< Near-infrared albedo of ice (-)
  real :: ALVSO  !< Visible albedo of organic matter (-)
  real :: ALIRO  !< Near-infrared albedo of organic matter (-)
  real :: ALBRCK !< Albedo of rock (-)

  real, dimension(18,4,2)  :: GROWYR !< The crop growth descriptor array (see calcLandSurfParams.f90 for description)
  real, dimension(ican)  :: ZORAT   !< the ratio of the roughness length for momentum to the roughness length for heat
  real, dimension(ican)  :: CANEXT  !< an attenuation coefficient used in calculating the sky view factor for
                                              !! vegetation canopies (variable c in the documentation for subroutine canopyAlbedoTransmiss)
  real, dimension(ican)  :: XLEAF   !< a leaf dimension factor used in calculating the leaf boundary resistance
                                              !! (variable Cl in the documentation for subroutine calcLandSurfParams.f90)

  !! Six hydraulic parameters associated with the three basic types of organic soils (fibric, hemic and sapric), see the documentation for subroutine soilProperties.

  real, dimension(3)  :: THPORG !< Organic soils (peat) pore volume \f$[m^3 m^{-3} ] ( \theta_p)\f$
  real, dimension(3)  :: THRORG !< Peat liquid water retention capacity for organic soil \f$[m^3 m^{-3} ] (\theta_{ret} )\f$
  real, dimension(3)  :: THMORG !< Peat residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ] (\theta_{min} )\f$
  real, dimension(3)  :: BORG !< Clapp and Hornberger 'b' value for peat soils [ ]
  real, dimension(3)  :: PSISORG !< Peat soil moisture suction at saturation [m] \f$(\Psi_{sat} )\f$
  real, dimension(3)  :: GRKSORG !< Peat hydraulic conductivity of soil at saturation \f$[m s^{-1} ] (K_{sat} )\f$

  ! Model physics timestep (s) (GCM name: DELTIM)
  real :: DELT

  ! canopyAlbedoTransmiss parameters: ----------------------
  real :: ALVSWC     !< Average background visible albedo of snow covered canopy (-)
  real :: ALIRWC     !< Average background NIR albedo of snow covered canopy (-)
  real :: CXTLRG     !< Effective overall exctinction coefficient given if visible transmissivity is very small

  ! soilProperties parameters:  ----------------------

  real, dimension(soilcolrinds)  :: ALWV !< Lookup tables for wet visible soil albedos based on soil colour index
  real, dimension(soilcolrinds)  :: ALWN !< Lookup tables for wet NIR soil albedos based on soil colour index
  real, dimension(soilcolrinds)  :: ALDV !< Lookup tables for dry visible soil albedos based on soil colour index
  real, dimension(soilcolrinds)  :: ALDN !< Lookup tables for dry NIR soil albedos based on soil colour index

  ! DRCOEF parameters: ----------------------
  ! FLAG has some but poorly documented COMBAK

  ! FLXSURFZ parameters: ----------------------
  ! FLAG has some but poorly documented COMBAK

  ! snowInfiltrateRipen parameters: ----------------------
  real :: WSNCAP !< Maximum water retention capacity of the snow pack (weight percentage)

  ! calcLandSurfParams parameters:

  ! NOTE: If IWF > 0 in joboptions file, indicating that lateral flow of soil
  ! water is being modelled, externally derived user-specified values of the
  ! ponding limit for the four subareas are assigned and the parameters below 
  ! will be ignored!
  
  real :: nonSoilPondLim      !< Bare soil ponding limit for rock surface or glacier (m)  
  real :: bareSoilPondLim     !< Bare soil ponding limit for soils (m)    
  real :: buriedGrassPondLim  !< Buried grass under snow ponding limit (m)
  real :: treePondLim         !< Treed area ponding limit (m)  
  real :: grassPondLim        !< Grass area ponding  limit (m)  
  real :: ICAPCWL             !< Interception capacity of canopy over bare ground for liquid water (unitless) 

  ! Four-band albedo information
  real :: REFF0_LAND  ! "climatological" effective radius of fresh snow grain size, SSA = 60 m^2/kg
  real :: ZSNMIN      ! The minimum depth of snowpack, below which the effective radius is kept constant\f$[m]\f$
  real :: ZSNMAX1     ! A shallow homogeneous layer of snow at the top of the snowpack, 
                      ! below which there is no transfer of BC between the layer and 
                      ! underlying snow except for removal of BC by meltwater scavenging (m)
  real :: ZSNMAX2     ! The maximum depth of the surface snow layer in parameterizations of snow grain size, 
                      ! corresponding to a typical depth of solar penetration in snow (m)

  ! Dry surface layer (DSL) parameters: (G. Meyer, April 2021)
  real :: DSLinit    !< Fraction of the porosity that determines at which moisture value the DSL initiates
  real :: DELTAZ     !< Maximum DSL thickness (m) from Swenson and Lawrence (2014)
  real :: PSIAIR     !< Air dry matric potential (m)

  ! ----------------------
  ! ----------------------

  ! Biogeochemical (CTEM) parameters:

  real, dimension(kk) :: kn               !< Canopy light/nitrogen extinction coefficient

  ! Moved definitions in "readInParams" here and allocate via "allocateParamsCTEM".
  
  real, dimension(kk) :: mxmortge_compete

  ! allocate.f90 parameters: ----------

  real, dimension(kk) :: omega            !< omega, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: epsilonl         !< Epsilon leaf, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: epsilons         !< Epsilon stem, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: epsilonr         !< Epsilon root, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: rtsrmin          !< Minimum root:shoot ratio mostly for support and stability
  real, dimension(kk) :: rtsrmin_peat     !< Peatland values for minimum root:shoot ratio mostly for support and stability 
  real, dimension(kk) :: aldrlfon         !< Allocation to leaves during leaf onset
  real, dimension(kk) :: caleaf           !< Constant allocation fractions to leaves if not using dynamic allocation.
                                                      !! (NOT thoroughly tested, and using dynamic allocation is preferable)
  real, dimension(kk) :: castem           !< Constant allocation fractions to stem if not using dynamic allocation.
                                                      !! (NOT thoroughly tested, and using dynamic allocation is preferable)
  real, dimension(kk) :: caroot           !< Constant allocation fractions to roots if not using dynamic allocation.
                                                      !! (NOT thoroughly tested, and using dynamic allocation is preferable)

  ! bio2str.f90 parameters: ---------

  real, dimension(kk) :: abar             !< parameter determining average root profile
  real, dimension(kk) :: abar_peat        !< Peatland values for parameter determining average root profile
  real, dimension(kk) :: avertmas         !< average root biomass (kg c/m2) for ctem's 8 pfts used for estimating rooting profile
  real, dimension(kk) :: alpha            !< parameter determining how the roots grow
  real, dimension(kk) :: prcnslai         !< storage/imaginary lai is this percentage of maximum leaf area index that a given root+stem biomass can support
  real, dimension(kk) :: minslai          !< minimum storage lai. this is what the model uses as lai when growing
                                                      !! vegetation for scratch. consider these as model seeds.
  real, dimension(kk) :: mxrtdpth         !< maximum rooting depth. this is used so that the rooting depths simulated by ctem's variable rooting depth
                                                      !! parameterzation are constrained to realistic values
  real, dimension(kk) :: albvis           !< visible albedos of the ctem pfts
  real, dimension(kk) :: albnir           !< near IR albedos of the 9 ctem pfts

  !! competitionMod.f90 parameters: ------

  !! the model basically uses the temperature of the coldest month as
  !! the major constraint for pft distribution. a range of the coldest
  !! month temperature is prescribed for each pft within which pfts are
  !! allowed to exist. in addition for tropical broadleaf drought
  !! deciduous trees measure(s) of aridity (function of precipitation
  !! and potential evaporation) are used.

  ! existence subroutine: ----------

  real, dimension(kk) :: tcoldmin         !< minimum coldest month temperature
  real, dimension(kk) :: tcoldmax         !< maximum coldest month temperature
  real, dimension(kk) :: twarmmax         !< maximum warmest month temperature
  real, dimension(kk) :: twarmmin         !< minimum warmest month temperature
  real, dimension(kk) :: gdd5lmt          !< minimum gdd above 5 c required to exist
  real, dimension(kk) :: aridmin          !< aridity index minimum
  real, dimension(kk) :: aridmax          !< aridity index maximum
  real, dimension(kk) :: dryseasonmin     !< minimum length of dry season for PFT to exist
  real, dimension(kk) :: dryseasonmax     !< maximum length of dry season for PFT to exist
  real, dimension(kk) :: annsrplsmin      !< minimum annual water surplus
  real, dimension(kk) :: bio2sap          !< multiplying factor for converting biomass density to sapling density
                                                      !! smaller numbers give faster colonization rates.
  real :: bioclimrt                                   !< mortality rate (1/year) for pfts that no longer exist within their pre-defined bioclimatic range

  ! ctem.f90 parameters: ----------

  real, dimension(kk) :: grescoef         !< Growth respiration coefficient
  real, dimension(kk) :: humicfac         !< Humification factor - used for transferring carbon from litter into soil c pool
  real :: humicfac_bg                                 !< Humification factor - used for transferring carbon from litter into soil c pool. This is the
                                                      !! value used for bare ground or LUC pools.
  real, dimension(kk) :: laimin           !< Minimum lai below which a pft doesn't expand
  real, dimension(kk) :: laimax           !< Maximum lai above which a pft always expands and lambdamax fraction of npp is used for expansion
  real :: lambdamax                       !< Max. fraction of npp that is allocated for reproduction/colonization
  real :: repro_fraction                  !< Fraction of NPP that is used to create reproductive tissues

  ! disturbance parameters: ------------

  real, dimension(2) :: bmasthrs_fire !< min. and max. vegetation biomass thresholds to initiate fire, \f$kg c/m^2\f$
  real :: extnmois_veg                !< extinction moisture content for estimating vegetation fire likeliness due to soil moisture
  real :: extnmois_duff               !< extinction moisture content for estimating duff layer fire likeliness due to soil moisture
  real :: lwrlthrs                    !< lower cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
  ! FireMIP value: 0.025
  real :: hgrlthrs                    !< higher cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
  ! FireMIP value: 1.0

  !> parameter m (mean) and b of logistic distribution used for \n
  !! **Parmlght was increased to 0.8 to make it so areas with higher amounts of
  !! lightning have higher lterm. The saturation is still the same, but the
  !! increase is more gradual at low lightning density. JM
  real :: parmlght
  real :: parblght                    !< estimating fire likelihood due to lightning
  real :: reparea                     !< typical area representing ctem's fire parameterization (km2)
  real :: popdthrshld                 !< threshold of population density (people/km2) [Kloster et al., biogeosci. 2010]
  real :: f0                          !< Fire spread rate in the absence of wind
  real, dimension(kk) :: maxsprd      !< max. fire spread rate, km/hr
  real, dimension(kk) :: frco2glf     !< fraction of green leaf biomass converted to gases due to combustion
  real, dimension(kk) :: frco2blf     !< fraction of brown leaf biomass converted to gases due to combustion
  real, dimension(kk) :: frltrglf     !< fraction of green leaf biomass becoming litter after combustion
  real, dimension(kk) :: frltrblf     !< fraction of brown leaf biomass becoming litter after combustion
  real, dimension(kk) :: frco2stm     !< fraction of stem biomass converted to gases due to combustion
  real, dimension(kk) :: frltrstm     !< fraction of stem biomass becoming litter after combustion
  real, dimension(kk) :: frco2rt      !< fraction of root biomass converted to gases due to combustion
  real, dimension(kk) :: frltrrt      !< fraction of root biomass becoming litter after combustion
  real, dimension(kk) :: frltrbrn     !< fraction of litter burned during fire and emitted as gases
  real, dimension(kk) :: standreplace !< pft prevalence for stand replacing fire events (based on resistance to fire damage, ie. cambial kill)(unitless)
  real, dimension(kk) :: emif_co2     !< pft-specific emission factors for CO2, g species / (kg DOM)
  real, dimension(kk) :: emif_co      !< pft-specific emission factors for CO, g species / (kg DOM)
  real, dimension(kk) :: emif_ch4     !< pft-specific emission factors for CH4, g species / (kg DOM)
  real, dimension(kk) :: emif_nmhc    !< pft-specific emission factors for non-methane hydrocarbons, g species / (kg DOM)
  real, dimension(kk) :: emif_h2      !< pft-specific emission factors for H2, g species / (kg DOM)
  real, dimension(kk) :: emif_nox     !< pft-specific emission factors for NOx, g species / (kg DOM)
  real, dimension(kk) :: emif_n2o     !< pft-specific emission factors for N2O, g species / (kg DOM)
  real, dimension(kk) :: emif_nh3     !< pft-specific emission factors for NH3, g species / (kg DOM)
  real, dimension(kk) :: emif_pm25    !< pft-specific emission factors for particles <2.5 micrometers in diameter, g species / (kg DOM)
  real, dimension(kk) :: emif_tpm     !< pft-specific emission factors for total particulate matter, g species / (kg DOM)
  real, dimension(kk) :: emif_tc      !< pft-specific emission factors for total carbon, g species / (kg DOM)
  real, dimension(kk) :: emif_oc      !< pft-specific emission factors for organic carbon, g species / (kg DOM)
  real, dimension(kk) :: emif_bc      !< pft-specific emission factors for black carbon, g species / (kg DOM)

  ! hetres parameters: ----------

  real, dimension(kk) :: bsratelt     !< litter respiration rates at 15 C in in kg C/kg C.year
  real, dimension(kk) :: bsratesc     !< soil carbon respiration rates at 15 C in kg C/kg C.year
  real, dimension(kk) :: bsratesc_peat !< Peatland values for soil carbon respiration rates at 15 C in kg C/kg C.year
  real, dimension(4) :: tanhq10       !< Constants used in tanh formulation of respiration Q10 determination
  real :: alpha_hetres                !< parameter for finding litter temperature as a weighted average of top soil layer temperature and root temperature
  real :: bsratelt_g                  !< bare ground litter respiration rate at 15 c in kg c/kg c.year
  real :: bsratesc_g                  !< bare ground soil c respiration rates at 15 c in kg c/kg c.year
  real :: a_hetr                      !< parameter describing exponential soil carbon profile. used for estimating temperature of the carbon pool
  real :: r_depthredu                 !< Following Lawrence et al. Environ. Res. Lett., 2015. we adopt 10.0 as our value, controls
                                      !! decomposition at depth to account for observed reductions that are independent of temp and moisture.
  real :: tcrit                       !< temperature below which respiration is inhibited. [ C ]
  real :: frozered                    !< factor to reduce respiration by for temps below tcrit

  ! landuseChangeMod.f90 parameters: --------------

  real, dimension(3) :: combust       !< how much deforested/chopped off biomass is combusted (these absolutely must add to 1.0 !)
  real, dimension(3) :: paper         !< how much deforested/chopped off biomass goes into short term storage such as paper
  real, dimension(3) :: furniture     !< how much deforested/chopped off biomass goes into long term storage such as furniture
  real, dimension(2) :: bmasthrs      !< biomass thresholds for determining if deforested area is a forest, a shrubland, or a bush kg c/m2

  ! mainres.f90 parameters: ----------

  !    Base respiration rates for stem and root for ctem pfts in kg c/kg c.year at 15 degrees celcius. note that maintenance
  !    respiration rates of root are higher because they contain both wood (coarse roots) and fine roots.
  !    New parameter values introduced to produce carbon use efficiencies more in
  !    line with literature (zhang et al. 2009, luyssaert et al. gcb 2007)
  !    values changed for bsrtstem and bsrtroot. jm 06.2012

  real, dimension(kk) :: bsrtstem     !< Base respiration rates for stem in kg c/kg c.year at 15 degrees celcius (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: bsrtroot     !< Base respiration rates for root in kg c/kg c.year at 15 degrees celcius (values differ if using prescribed vs. competition run)
  real :: minlvfr                     !< Minimum live wood fraction

  ! peatlandsMod.f90 parameters:

  ! mossPht subroutine parameters
  
  real :: mossMaxSnow                 !< Maximum snow depth under which moss photosynthesis is allowed to occur (m)
  real :: mossMinTemp                 !< Minimum ground surface temperature under which moss photosynthesis is allowed to occur (deg C)
  real :: rmlmoss25                   !< Base dark respiration rate of moss umol CO2/m2/s
  real :: tau25m                      !< tau coefficient (rate at 25 Celcius) iscomputed on the basis of the specificity
                                      !! factor (102.33) times Kco2/Kh2o (28.38) to convert for value in solution to that based in air.
  real :: ektau                       !< \f$J mol^{-1}\f$  ! FLAG this is???
  real :: kc25                        !< kinetic coef for CO2 at 25 C (Pa)
  real :: ko25                        !< kinetic coef for O2 at 25C (rate)  (Pa)
  real :: ec                          !< Activation energy for Kc of CO2 (\f$J mol^{-1}\f$) at 25C
  real :: ej                          !< Activation energy for electron transport, (\f$J mol^{-1}\f$)
  real :: eo                          !< Activation energy for Ko (\f$J mol^{-1}\f$) at 25C
  real :: evc                         !< Activation energy for carboxylation (\f$J mol^{-1}\f$) at 25C
  real :: sj                          !< Constant affecting J at low temperature (\f$J mol^{-1}\f$) at 25C
  real :: hj                          !< Constant affecting J at high temperature (\f$J mol^{-1}\f$) at 25C
  real :: alpha_moss                  !< Efficiency of conversion of incident photons into electrons(mol electron/mol photon),

  ! decp subroutine parameters:

  real ::     dctmin                   !< minimum temperature of soil respiration, K (peatland soils)
  real ::     dcbaset                  !< base temperature for Q10, K (peatland soils)
  real ::     bsrateltms               !< heterotrophic respiration base rate for peatlands (yr-1)

  ! peatland parameters used in other subroutines:

  real :: zolnmoss                    !< natual logarithm of roughness length of moss surface.
  real :: thpmoss                     !< pore volume of moss (\f$m^3/m^3\f$)
  real :: thrmoss                     !< Liquid water retention capacity of moss (\f$m^3/m^3\f$)
  real :: thmmoss                     !< residual liquid water content after freezing or evaporation of moss (\f$m^3/m^3\f$)
  real :: bmoss                       !< Clapp and Hornberger empirical "b" parameter of moss
  real :: psismoss                    !< Soil moisure suction at saturation of moss (m)
  real :: grksmoss                    !< Saturated hydrualic conductivity of moss (m)
  real :: hcpmoss                     !< Volumetric heat capacity of moss (\f$J m^{-3} K^{-1}\f$)

  ! real :: sphms = 2.70E3      !< Specific heat of moss layer (\f$J m^{-2} K^{-1}\f$) ! FLAG NOT USED also check units - written as J/kg/K. JM Sep 26 2016
  !! same as that specific heat for vegetation
  !! (Berlinger et al 2001)

  ! real :: rhoms = 40.0        !< Density of dry moss (\f$kg m^{-3}\f$)  ! FLAG NOT USED JM Oct 2016.
  !! Based upon:
  !! 40.0 Price et al. 2008, Price and
  !! Whittington 2010 value for feather moss
  !! is lower than sphagnmum (20 kg/m3 in
  !! O'Donnell et al. 2009)

  ! real :: slams = 20.0        !< Specific leaf area of moss (\f$m^2 kg^{-1}\f$), ! FLAG NOT USED JM Oct 2016.
  !! Based upon:
  !! S vensson 1995 sphangum value ranges from
  !! 13.5~47.3 m2/kg (Lamberty et al. 2007)

  real ::   grescoefmoss               !< Moss growth respiration coefficient
  real ::   rmortmoss                  !< Moss mortality rate constant (year-1)
  real ::   humicfacmoss               !< Humification ratio of moss litter

  ! phenology.f90 parameters: ----------

  integer, dimension(kk) :: dayschk       !< Number of days over which to check if net photosynthetic rate is positive before initiating leaf onset
  real, dimension(kk) :: drgta            !< Parameter determining how fast soil dryness causes leaves to fall
  real, dimension(kk) :: eta              !< eta, parameter for estimating min. stem+root biomass
  real, dimension(kk) :: eta_peat         !< Peatland values for eta, parameters for estimating min. stem+root biomass
  real, dimension(kk) :: kappa            !< required to support green leaf biomass.
  real, dimension(kk) :: kappa_peat       !< Peatland values for required to support green leaf biomass. 
  real, dimension(2) :: flhrspan                      !< Harvest span (time in days over which crops are harvested, 15 days),
                                                      !! and fall span (time in days over which bdl cold dcd plants shed their leaves, 30 days)
  real :: fracbofg                                    !< Parameter used to estimate lai of brown leaves. We assume that SLA of brown leaves is this fraction of SLA of green leaves
  real, dimension(kk) :: harvthrs         !< LAI threshold for harvesting crops. values are zero for all pftsother than c3 and c4 crops.
  real, dimension(kk) :: specsla          !< CTEM can use user-specified specific leaf areas (SLA) if the following specified values are greater than zero
  real, dimension(kk) :: thrprcnt         !< Percentage of max. LAI that can be supported which is used as a threshold for determining leaf phenology status
  real, dimension(kk) :: lwrthrsh         !< Lower temperature threshold for ctem's 9 pfts. these are used to estimate cold stress related leaf loss rate (degree c)
  real, dimension(kk) :: cdlsrtmx         !< Max. loss rate for cold stress for all 9 pfts, (1/day)
  real :: kmort1                                      !< kmort1, parameter used in growth efficiency mortality formulation

  real, dimension(kk) :: mxmortge !< Maximum mortality when growth efficiency is zero (1/year) (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: maxage   !< Maximum plant age. used to calculate intrinsic mortality rate. (yrs)
                                              !! maximum age for crops is set to zero since they will be harvested
                                              !! anyway. grasses are treated the same way since the turnover time
                                              !! for grass leaves is 1 year and for roots is 2 years. (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: maxage_peat   !< Peatland values for maximum plant age (yrs)

  real, dimension(kk) :: lfespany !< Leaf life span (in years) for biogeochemical pfts
  real, dimension(kk) :: lfespany_peat !< Peatland values for leaf life span (in years) for biogeochemical pfts
  real, dimension(kk) :: drlsrtmx !< Max. loss rate for drought stress for all 9 pfts, (1/day) (values differ if using prescribed vs. competition run)
  real, dimension(kk) :: colda    !< Parameter determining how fast cold temperatures causes leaves to fall
  integer, dimension(2) :: coldlmt !< No. of days for which some temperature has to remain below a given threshold for initiating a process, days
  real, dimension(2) :: coldthrs  !< 1. -5 c threshold for initiating "leaf fall" mode for ndl dcd trees \n
  !! 2.  8 c threshold for initiating "harvest" for crops, the array colddays_leaffall and colddays_harvest tracks days corresponding to these thresholds, respectively.
  real :: roothrsh                !< Root temperature threshold for initiating leaf onset for cold broadleaf deciduous pft, degrees celcius

  ! Turbation (in soilCProcesses.f90) parameters: ---------------------------------

  real :: cryodiffus            !< Diffusivity (or simply the rate of the cryoturbation) \f$(m^2/d)\f$
                                !! value from \cite Koven2011-796 - 5 cm2/yr.
  real :: biodiffus             !< Diffusivity (or simply the rate of the bioturbation) \f$(m^2/d)\f$
                                !! value from \cite Koven2011-796 - 1 cm2/yr.
  real :: kterm                 !< Constant used in determination of depth at which cryoturbation ceases
  real :: eftime                !< e-folding time scale for updating mean active layer depth used in determining the depth
                                !! roots can penetrate. (years)

  real :: efoldfact             !< Calculated factor used in determining the e-folding time of average active layer depth.
                                !! this is calculated based on the eftime parameter.


  ! soil_ch4uptake parameters: -----------------------------

  real :: D_air    !< Diffusivity of CH4 in air (\f$cm^2 s^-1\f$) @ STP
  real :: g_0      !< Scaling factor takes CH4_soills to mg CH4 \f$m^-2 s^-1\f$ (units: \f$mg CH_4 ppmv^{-1} s s^{-1} m^{-2} cm{-1}\f$)
  real :: betaCH4  !< Constant derived in Curry (2007) from comparison against measurements (-)
  real :: k_o      !< Base oxidation rate derived in Curry (2007) from comparison against measurements \f$(s^{-1})\f$

  ! turnover.f90 parameters: ---------------------------------

  real, dimension(kk) :: stemlife !< Stemlife, turnover time scale for stem for different pfts
  real, dimension(kk) :: rootlife !< Rootlife, turnover time scale for root for different pfts
  real :: stmhrspn                !< Stem harvest span. same as crop harvest span. period in days over which crops are harvested.

  ! wetland_methane.f90 parameters: ------------

  real :: ratioch4                !< methane to carbon dioxide flux scaling factor.

  !> ratio of wetland to upland respiration\n
  !! Use the heterotrophic respiration outputs for soil and litter as the ecosystem basis.  These were summed as "hetrores".
  !! This respiration is for upland soils; we multiply by wtdryres as the ratio of wetland to upland respiration
  !! based on literature measurements: Dalva et al. 1997 found 0.5 factor; Segers 1998 found a 0.4 factor. use 0.45 here (unitless)
  real, dimension(3) :: wtdryres
  real :: lat_thrshld1   !< Northern zone for wetland determination (degrees North)
  real :: lat_thrshld2   !< Boundary with southern zone for wetland determination (degrees North)
  real, dimension(2) :: soilw_thrshN   !< Soil wetness threshold in the North zone
  real, dimension(2) :: soilw_thrshE   !< Soil wetness threshold in the Equatorial zone
  real, dimension(2) :: soilw_thrshS   !< Soil wetness threshold in the South zone

  ! Photosynthesis parameters: --------------------------------------------------

  logical, dimension(kk) :: isc4      !< Array telling which vegetation type is c4
  real, dimension(kk) :: tlow         !< lower temperature limits for photosynthesis, (kelvin)
  real, dimension(kk) :: tup          !< upper temperature limits for photosynthesis, (kelvin)
  real, dimension(kk) :: alpha_phtsyn !< quantum efficiencies, values of 0.08 & 0.04 are used for c3 and c4 plants, respectively
  real, dimension(kk) :: omega_phtsyn !< leaf scattering coefficients, values of 0.15 & 0.17 are used for c3 and c4 plants, respectively
  real, dimension(kk) :: mm           !< parameter m used in photosynthesis-stomatal conductance coupling.
  real, dimension(kk) :: bb           !< parameter b used in photosynthesis-stomatal conductance coupling.
  real, dimension(kk) :: vpd0         !< parameter vpd0 used in leuning type photosynthesis - stomatal conductance coupling, (Pa)
  integer, dimension(kk) :: sn        !< exponent for soil moisture stress. for sn equal to 1, photosynthesis decreases
                                                  !! linearly with soil moisture, and of course non-linearly for values higher than 1.
                                                  !! when sn is about 10, photosynthesis does not start decreasing until soil moisture
                                                  !! is about half way between wilting point and field capacity.
  real, dimension(kk) :: smscale      !< additional constrain of soil moisture stress on photosynthesis. this can be used
                                                  !! to simulate the effect of irrigation for crops.
  real, dimension(kk) :: vmax         !< max. photosynthetic rate, (\f$(mol CO_2 m^{-2} s^{-1}\f$) values are mainly derived from
                                                  !! \cite kattge20090c0 which doesn't include c4. also see \cite alton2017-pd
  real, dimension(kk) :: vmax_peat    !< Peatland values for max. photosynthetic rate, (\f$(mol CO_2 m^{-2} s^{-1}\f$)
  integer :: reqiter                              !< no. of iterations for calculating intercellular co2 concentration
  real :: co2imax                                 !< max. intercellular co2 concentration (Pa)
  real :: beta1, beta2                            !< photosynthesis coupling or curvature coefficients
  real, dimension(kk) :: inico2i      !< parameter to initialize intercellular co2 conc.
  real, dimension(kk) :: rmlcoeff     !< leaf maintenance respiration coefficients
  real, dimension(kk) :: chi          !< additional parameters for two-leaf model leaf angle distribution
  real :: gamma_w                                 !< photosynthesis down regulation parameters equivalent co2 fertilization effect that we want model to yield
  real :: gamma_m                                 !< equivalent co2 fertilization effect that model actually gives without any photosynthesis down-regulation

  ! non-structural carbon parameters: ----------
  real, dimension(kk) :: Fmax_leaf    !< Maximum flux from leaf nsc pool to leaf structural pool.
  real, dimension(kk) :: Fmax_stem    !< Maximum flux from stem nsc pool to stem structural pool.
  real, dimension(kk) :: Fmax_root    !< Maximum flux from root nsc pool to root structural pool.
  real, dimension(kk) :: Fmax_leaf_N  !< Maximum flux from leaf non-structural N pool to leaf structural N pool.
  real, dimension(kk) :: Fmax_stem_N  !< Maximum flux from stem non-structural N pool to stem structural N pool.
  real, dimension(kk) :: Fmax_root_N  !< Maximum flux from root non-structural N pool to root structural N pool.
  real, dimension(kk) :: min_ns2t_l   !< Minimum (nsc/total) leaf pool ratio, following table 3 in N. Li et al. (2016).
  real, dimension(kk) :: min_ns2t_s   !< Minimum (nsc/total) stem pool ratio, following table 3 in N. Li et al. (2016).
  real, dimension(kk) :: min_ns2t_r   !< Minimum (nsc/total) root pool ratio, following table 3 in N. Li et al. (2016).
  real, dimension(kk) :: gleafmasmax  !< Leaf pool's max for each PFT estimated from model values during normal growth mode.
  real, dimension(kk) :: stemmassmax  !< Stem pool's max for each PFT estimated from model values during normal growth mode.
  real, dimension(kk) :: rootmassmax  !< Root pool's max for each PFT estimated from model values during normal growth mode.
  real, dimension(kk) :: re_alloc     !< Coeff. used in the reallocation of non-structural C during leaf out period.

  !dynamic tiling parameters: -------------------------------------------------
  real, parameter :: dynTilingTolrance = 0.000001 !< the value used for logical tests in the final mass balance checks where the checks should yield values very near zero
  real, parameter :: dynTilingTolranceEnergy = 0.001 !< the value used for logical tests in the final energy checks (this is a percentage becase these values are much larger)
  real, parameter :: dynTilingMaxDisturb = 0.999
  
  ! nitrogen cycle parameters: -------------------------------------------------
  integer :: fertstd           !< start day of applying fertilizer (or end day if Southern Hemisphere)
  integer :: fertend           !< end day of applying fertilizer (or start day if Souther Hemisphere)
  real :: k_redcoeff_vcmax     !< parameter describing reduction of V_cmax
  real :: bnfdpth              !< maximum depth (m) of occurance for the BNF following Vitousek et al. (2013)
  real :: r_bnf_f              !< rate of free-living BNF
  real :: C_cost_bnf_s         !< carbon cost of symbiotic BNF
  real :: nvol_coeff           !< coeff. for nitrogen volatilization
  real :: nh3atm               !< concentration of NH3 in the free atmosphere (\f$gN/m^3\f$), 0.3 micro.gN/m3 following Riddick et al. (2016)
  real :: mu_nit               !< coeff. for nitrification (\f$day^{-1}\f$), following Lin et al. (2000)
  real :: topt_nit             !< optimum temperature for nitrification in degree C, following Riedo et al. (1998)
  real :: nitdenitdpth         !< maximum depth (m) of occurance of nitrification and denitrification
  real :: mu_no_nit            !< coeff. for nitrification-induced NO loss
  real :: mu_n2o_nit           !< coeff. for nitrification-induced N2O loss
  real :: mu_denit             !< coeff. for denitrification (\f$day^{-1}\f$), following Lin et al. (2000)
  real :: topt_denit           !< optimum temperature for denitrification in degree C, following Li et al. (2000)
  real :: mu_no_denit          !< coeff. for NO loss through denitrification
  real :: mu_n2o_denit         !< coeff. for N2O loss through denitrification
  real :: denitswthrsh         !< soil wetness onset threshold for denitrification.
  real :: ntolerance           !< our tolerance for balancing N budget in g N/m2
  real :: nleach_coeff         !< coeff. for nitrogen leaching
  real, dimension(kk) :: coeff_a              !< coeff. for active nitrogen uptake, following Wania et al. (2012), (\f$g N/g rootC day\f$)
  real :: coeff_p_nh4          !< coeff. to control how much of the NH4+ is accessible for passive nitrogen uptake
  real :: coeff_p_no3          !< coeff. to control how much of the NO3- is accessible for passive nitrogen uptake
  real :: b_nh4                !< sorption/desorption buffer factor for NH4+, following Wania et al. (2012)
  real :: b_no3                !< sorption/desorption buffer factor for NO3-, following Wania et al. (2012)
  real :: kphalf               !< half-saturation constant for nitrogen uptake, following Wania et al. (2012)
  real :: const_c2n_humus      !< constant C:N ratio for the humus pool used in nMineralImm.f90
  real :: initial_c2n          !< arbitrary value for the initial C:N ratios
  real :: ccost_coeff          !< coeff. for carbon cost associated with nitrogen active uptake (\f$ kg C/kg N\f$)
  real, dimension(kk) :: nresorp_coeff !< leaf nitrogen resorption efficiency coefficient, following Yasumura et al. (2005)
  real, dimension(kk) :: c2n_lmax     !< max leaf C:N ratio for individual PFTs, following Wania et al. (2012)
  real, dimension(kk) :: c2n_smax     !< max stem C:N ratio for individual PFTs, following Wania et al. (2012)
  real, dimension(kk) :: c2n_rmax     !< max root C:N ratio for individual PFTs, following Wania et al. (2012)
  real, dimension(kk) :: c2n_lmin     !< min leaf C:N ratio for individual PFTs, following Wania et al. (2012)
  real, dimension(kk) :: c2n_smin     !< min stem C:N ratio for individual PFTs, following Wania et al. (2012)
  real, dimension(kk) :: c2n_rmin     !< min root C:N ratio for individual PFTs, following Wania et al. (2012)
  real, dimension(kk) :: alpha_vcmax  !< PFT-dependent coeff. to relate Vcmax to leaf N content
  real, dimension(kk) :: beta_vcmax   !< PFT-dependent coeff. to relate Vcmax to leaf N content
  real, dimension(kk) :: ngleafmasmax !< L pool's max for each PFT estimated from model values during normal growth mode.
  real, dimension(kk) :: nstemmassmax !< S pool's max for each PFT estimated from model values during normal growth mode.
  real, dimension(kk) :: nrootmassmax !< R pool's max for each PFT estimated from model values during normal growth mode.
  real, dimension(kk) :: r_bnf_s      !< rate of symbiotic BNF for each PFT
  real, dimension(kk) :: max_frac_bnf      !< max. fraction of gleafmas_ns allocated to symbiotic BNF
  real, dimension(kk) :: b_bnf                !< intercept of symbiotic BNF

  ! tiled disturbance parameters:---------------------------------------------------
  integer :: nTilepres                  !< the number of small tiles to preserve for dynamic tiling (nTilepres # tiles)
  real :: percentthreshold              !< height difference for tile to be joined for dynamic tiling as a percentage of max height

  ! Passed variables:

  character(350)    :: runParamsFile
  logical          :: PFTCompetitionSwitch
  real :: zbldJobOpt, zrfhJobOpt, zrfmJobOpt
  real, dimension(ican) :: RSMN, QA50, VPDA, VPDB, PSGA, PSGB ! These are temporary variables storing these parameters until they are
                                                              ! transfered into their ROT structure in read_initialstate
  !! @}

  ! --------------------------------------------------------------------------

contains

  !> \ingroup classicparams_prepareGlobalParams
  !! @{
  !> Initialize and/or read in all global model parameters
  subroutine prepareGlobalParams

    implicit none

    !> Read in the number of PFTs to check against whats set in this file, this reads them in from the namelist file.
    open(10,file = trim(runParamsFile),action = 'read',status = 'old')
    call readInPFTNums(10)
    close(10)

    !> Initialize the CTEM parameters, this reads them in from a namelist file.
    open(10,file = trim(runParamsFile),action = 'read',status = 'old')
    call readInParams(10)
    close(10)

  end subroutine prepareGlobalParams
  !! @}
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  !> \ingroup classicParams_readInPFTNums
  !! @{
  !> Read in the number of PFTs from the parameters namelist file and test whats set in classicParams.f90.
  subroutine readInPFTNums(inunit)

    implicit none
    integer :: inunit

   namelist /classicPFTbasic/ &
      testican, &
      testicc, &
      testl2max

   read(inunit,nml = classicPFTbasic)

   IF (testican /= ican) THEN
    write(6, * ) 'MISMATCH IN NUMBER OF CLASS PFTs'
    write(6, * ) ican, 'set in classicParams.f90'
    write(6, * ) testican, 'set in the parameters file'
    print*,'----- Run aborting. -----'
    stop
   END IF

   IF (testicc /= icc) THEN
    write(6, * ) 'MISMATCH IN NUMBER OF CTEM PFTs'
    write(6, * ) icc, 'set in classicParams.f90'
    write(6, * ) testicc, 'set in the parameters file'
    print*,'----- Run aborting. -----'
    stop
   END IF

   IF (testl2max /= l2max) THEN
    write(6, * ) 'MISMATCH IN MAX NUMBER OF LEVEL 2 CTEM PFTs'
    write(6, * ) l2max, 'set in classicParams.f90'
    write(6, * ) testl2max, 'set in the parameters file'
    print*,'----- Run aborting. -----'
    stop
   END IF

  end subroutine readInPFTNums
  !! @}
  ! --------------------------------------------------------------------------

  !> \ingroup classicParams_readInParams
  !! @{
  !> Read in the CLASSIC parameters from a namelist file. Populate a few parameters
  !! based on what was read in.

  subroutine readInParams(inunit)

    implicit none

    integer :: i, n, m
    character(8) :: pftkind
    integer :: isumc, k1c, k2c, inunit

    namelist /classicparams/ &
        VKC, &
        CT   , &
        VMIN , &
        TCW  , &
        TCICE, &
        TCSAND, &
        TCCLAY, &
        TCOM  , &
        TCDRYS, &
        RHOSOL, &
        RHOOM , &
        HCPW  , &
        HCPICE, &
        HCPSOL, &
        HCPOM , &
        HCPSND, &
        HCPCLY, &
        SPHW  , &
        SPHICE, &
        SPHVEG, &
        RHOW  , &
        RHOICE, &
        TCGLAC, &
        CLHMLT, &
        CLHVAP, &
        ZOLNG  , &
        ZOLNS  , &
        ZOLNI  , &
        ZORATG , &
        ALVSI  , &
        ALIRI  , &
        ALVSO  , &
        ALIRO  , &
        ALBRCK , &
        GROWYR , &
        ZORAT  , &
        CANEXT , &
        XLEAF  , &
        RSMN, &
        QA50, &
        PSGA, &
        PSGB, &
        VPDA, &
        VPDB, &
        THPORG , &
        THRORG , &
        THMORG , &
        BORG , &
        PSISORG, &
        GRKSORG, &
        DELT, &
        ALVSWC, &
        ALIRWC, &
        CXTLRG, &
        ALWV, &
        ALWN, &
        ALDV, &
        ALDN, &
        WSNCAP, &
        nonSoilPondLim, &
        bareSoilPondLim, &
        buriedGrassPondLim, &
        treePondLim, &
        grassPondLim, &
        ICAPCWL, &
        REFF0_LAND, &
        ZSNMIN, &
        ZSNMAX1, &
        ZSNMAX2, &
        DSLinit, &
        DELTAZ, &
        PSIAIR, &
        modelpft, &
        classpfts, &
        ctempfts, &
        kn, &
        omega, &
        epsilonl, &
        epsilons, &
        epsilonr, &
        rtsrmin, &
        rtsrmin_peat, &
        aldrlfon, &
        caleaf, &
        castem, &
        caroot, &
        abar, &
        abar_peat, &
        avertmas, &
        alpha, &
        prcnslai, &
        minslai, &
        mxrtdpth, &
        albvis, &
        albnir, &
        tcoldmin, &
        tcoldmax, &
        twarmmax, &
        twarmmin, &
        gdd5lmt, &
        aridmin, &
        aridmax, &
        dryseasonmin, &
        dryseasonmax, &
        annsrplsmin, &
        bio2sap, &
        bioclimrt, &
        grescoef, &
        humicfac, &
        humicfac_bg, &
        laimin, &
        laimax, &
        lambdamax, &
        repro_fraction, &
        bmasthrs_fire, &
        extnmois_veg, &
        extnmois_duff, &
        lwrlthrs, &
        hgrlthrs, &
        parmlght, &
        parblght, &
        reparea, &
        popdthrshld, &
        f0, &
        maxsprd, &
        frco2glf, &
        frco2blf, &
        frltrglf, &
        frltrblf, &
        frco2stm, &
        frltrstm, &
        frco2rt, &
        frltrrt, &
        frltrbrn, &
        standreplace, &
        emif_co2, &
        emif_co, &
        emif_ch4, &
        emif_nmhc, &
        emif_h2, &
        emif_nox, &
        emif_n2o, &
        emif_nh3, &
        emif_pm25, &
        emif_tpm, &
        emif_tc, &
        emif_oc, &
        emif_bc, &
        bsratelt, &
        bsratesc, &
        bsratesc_peat, &
        tanhq10, &
        alpha_hetres, &
        bsratelt_g, &
        bsratesc_g, &
        a_hetr, &
        r_depthredu, &
        tcrit, &
        frozered, &
        combust, &
        paper, &
        furniture, &
        bmasthrs, &
        bsrtstem, &
        bsrtroot, &
        minlvfr, &
        mossMaxSnow, &
        mossMinTemp, &
        rmlmoss25, &
        tau25m, &
        ektau, &
        kc25, &
        ko25, &
        ec, &
        ej, &
        eo, &
        evc, &
        sj, &
        hj, &
        alpha_moss, &
        dctmin, &
        dcbaset, &
        bsrateltms, &
        zolnmoss, &
        thpmoss, &
        thrmoss, &
        thmmoss, &
        bmoss, &
        psismoss, &
        grksmoss, &
        hcpmoss, &
        grescoefmoss, &
        rmortmoss, &
        humicfacmoss, &
        dayschk, &
        drgta, &
        eta, &
        eta_peat, &
        kappa, &
        kappa_peat, &
        flhrspan, &
        fracbofg, &
        harvthrs, &
        specsla, &
        thrprcnt, &
        lwrthrsh, &
        cdlsrtmx, &
        kmort1, &
        mxmortge, &
        mxmortge_compete, &
        maxage, &
        maxage_peat, &
        lfespany, &
        lfespany_peat, &
        drlsrtmx, &
        colda, &
        coldlmt, &
        coldthrs, &
        roothrsh, &
        cryodiffus, &
        biodiffus, &
        kterm, &
        eftime, &
        D_air, &
        g_0, &
        betaCH4, &
        k_o, &
        stemlife, &
        rootlife, &
        stmhrspn, &
        ratioch4, &
        wtdryres, &
        lat_thrshld1, &
        lat_thrshld2, &
        soilw_thrshN, &
        soilw_thrshE, &
        soilw_thrshS, &
        kn, &
        tup, &
        tlow, &
        alpha_phtsyn, &
        omega_phtsyn, &
        isc4, &
        mm, &
        bb, &
        vpd0, &
        sn, &
        smscale, &
        vmax, &
        vmax_peat, &
        reqiter, &
        co2imax, &
        beta1, &
        beta2, &
        inico2i, &
        chi, &
        rmlcoeff, &
        gamma_w, &
        gamma_m, &
        Fmax_leaf, &
        Fmax_stem, &
        Fmax_root, &
        Fmax_leaf_N, &
        Fmax_stem_N, &
        Fmax_root_N, &
        min_ns2t_l, &
        min_ns2t_s, &
        min_ns2t_r, &
        gleafmasmax, &
        stemmassmax, &
        rootmassmax, &
        ngleafmasmax, &
        nstemmassmax, &
        nrootmassmax, &
        re_alloc, &
        nresorp_coeff, &
        c2n_lmax, &
        c2n_smax, &
        c2n_rmax, &
        c2n_lmin, &
        c2n_smin, &
        c2n_rmin, &
        alpha_vcmax, &
        beta_vcmax, &
        r_bnf_s, &
        max_frac_bnf, &
        b_bnf, &
        fertstd, &
        fertend, &
        k_redcoeff_vcmax, &
        bnfdpth, &
        r_bnf_f, &
        C_cost_bnf_s, &
        nvol_coeff, &
        nh3atm, &
        mu_nit, &      
        topt_nit, &   
        nitdenitdpth, &   
        mu_no_nit, &      
        mu_n2o_nit, &     
        mu_denit, &       
        topt_denit, &     
        mu_no_denit , &   
        mu_n2o_denit  , & 
        denitswthrsh  , & 
        ntolerance    , & 
        nleach_coeff  , & 
        coeff_a       , & 
        coeff_p_nh4   , & 
        coeff_p_no3   , & 
        b_nh4         , & 
        b_no3         , & 
        kphalf        , & 
        const_c2n_humus, &
        initial_c2n   , & 
        ccost_coeff, &
        nTilepres, &
        percentthreshold    


    ! ----------

    read(inunit,nml = classicparams)

    numcrops    = 0        !< number of crop pfts
    numtreepfts = 0        !< number of tree pfts
    numgrass    = 0        !< number of grass pfts
    numshrubs   = 0        !< number of shrubs pfts
    grass = .false.
    crop = .false.

    ! Calculate number of level 2 pfts using modelpft
    do i = 1,ican
      isumc = 0
      k1c = (i - 1) * l2max + 1
      k2c = k1c + (l2max - 1)
      do n = k1c,k2c
        if (modelpft(n) == 1) isumc = isumc + 1
      end do
      nol2pfts(i) = isumc  ! number of level 2 pfts
    end do

    ! Find the reindexing array to go from CLASS PFTs to CTEM PFTs in the parameter
    ! arrays.
    reindexPFTs = findPFTindexes()

    ! Calculate the CL4CTEM which helps index CTEM to CLASS in calcLandSurfParams.
    do i = 1,ican
      do m = reindexPFTs(i,1),reindexPFTs(i,2)
        CL4CTEM(m) = i
      end do
    end do

    do i = 1,icc
      pftkind = ctempfts(i)
      select case (pftkind)
      case ('NdlEvgTr','NdlDcdTr', 'BdlEvgTr','BdlDCoTr', 'BdlDDrTr')
        numtreepfts = numtreepfts + 1
      case ('BdlDCoSh','BdlEvgSh')  ! FLAG NEED FIX.
        numshrubs = numshrubs + 1
      case ('GrassC3 ','GrassC4 ','Sedge')
        numgrass = numgrass + 1
        grass(i) = .true.
      case ('CropC3  ','CropC4  ')
        numcrops = numcrops + 1
        crop(i) = .true.
      case default
        print * ,'Unknown PFT in classicParams ',pftkind
        call errorHandler('classicParams', - 1)
      end select
    end do

    ! Overwrite the prescribed vars with the compete ones if competition is on.
    ! NOTE: maxage has a peatland specific parameter set. This one below is 
    ! overwritten by the peatland one in mortality if a tile is a peatland!
    if (PFTCompetitionSwitch) then
      mxmortge = mxmortge_compete
    end if

    ! Determine the efoldfact parameter which is calculated from eftime.
    efoldfact = exp( - 1.0 / eftime)
    
    ! Assign these values, may be overwritten (in findLeapYears in generalUtils) if using leap years
    monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] ! days in each month
    monthend  = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] ! calender day at end of each month
    mmday     = [ 16,46,75,106,136,167,197,228,259,289,320,350 ] ! mid-month day

  end subroutine readInParams
  !! @}

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup classicparams_reindexPFTs
  !! @{
  !> Provides indexes to relate the CLASS PFTs to the CTEM ones in the
  !! parameter arrays.
  function findPFTindexes ()

    implicit none

    integer :: j, k1c, k2c
    integer :: findPFTindexes(ican,2)

    k1c = 0
    do j = 1,ican
      if (j == 1) then
        k1c = k1c + 1
      else
        k1c = k1c + nol2pfts(j - 1)
      end if
      k2c = k1c + nol2pfts(j) - 1
      findPFTindexes(j,1) = k1c
      findPFTindexes(j,2) = k2c
    end do

  end function findPFTindexes
  !! @}
  !> \namespace classicparams
  !! This module holds CLASSIC globally accessible parameters
  !!
  !! These parameters are used in all CLASSIC subroutines
  !! via use statements pointing to this module.
  !! The structure of this subroutine is variables that are common to competition/prescribe PFT fractions
  !! first, then the remaining variables are assigned different variables if competition is on, or not.
  !!
  !! Please see documentation/overviewCTEM.md for more info.
  !!
  !! PFT parameters
  !!
  !! Note the structure of vectors which clearly shows the CLASS
  !! PFTs (along rows) and CTEM sub-PFTs (along columns)
  !!
  !! The basic version of CLASSIC includes the following PFTs:
  !! \f[
  !! \begin{array} { | l | c | c | c | c | c | }
  !! \hline
  !! \text{CLASS PFTs} &  \text{CTEM PFTs} &     \text{---} &      \text{---}  \\ \hline
  !! \hline
  !! \text{needle leaf} &  \text{evg} &      \text{dcd} &      \text{---}  \\ \hline
  !! \text{broad leaf}  &  \text{evg} &  \text{dcd-cld} &  \text{dcd-dry} \\ \hline
  !! \text{crops}       &   \text{c3} &       \text{c4} &      \text{---}  \\ \hline
  !! \text{grasses}     &   \text{c3} &       \text{c4} &    \text{---} \\ \hline
  !! \end{array}
  !! \f]
  !!
  !! The peatland module expands upon these PFTs with the addition of two
  !! shrub PFTs (only biogeochemistry) and sedges (see Wu et al. 2016) \cite Wu2016-zt Evergreen
  !! shrubs, for example the ericaceous shrubs, are the common dominant vascular
  !! plants in bogs and poor fens while deciduous shrubs, such as the betulaceous
  !! shrubs, often dominate rich fens. Both shrubs are categorized as broadleaf
  !! trees in CLASS morphologically, but their phenological and physiological
  !! characteristics are more similar to those of needleleaf trees. The shrub
  !! tundra ecosystem is situated adjacent to needleleaf forest in the Northern
  !! Hemisphere (Kaplan et al., 2003) and they share similar responses to climate
  !! in ESMs (e.g. Bonan et al., 2002).  (The photosynthesis and autotrophic
  !! respiration of vascular PFTs are modelled the same as in the original CTEM.)
  !!
  !! \f[
  !! \begin{array} { | l | c | c | c | c | c | }
  !! \hline
  !! \text{CLASS PFTs} &  \text{CTEM PFTs} &     \text{---} &      \text{---}  \\ \hline
  !! \hline
  !! \text{needle leaf} &  \text{evg} &      \text{dcd} &      \text{---} & \text{---} & \text{---}  \\ \hline
  !! \text{broad leaf}  &  \text{evg} &  \text{dcd-cld} &  \text{dcd-dry} &  \text{evg-shrubs} & \text{dcd-shrubs}\\ \hline
  !! \text{crops}       &   \text{C3} &       \text{C4} &      \text{---}  & \text{---} & \text{---}\\ \hline
  !! \text{grasses}     &   \text{C3} &       \text{C4} &      \text{sedges} & \text{---} & \text{---}\\ \hline
  !! \end{array}
  !! \f]
  !! Shrubs have now been implimented in the physic subroutines now as a fifth CLASS
  !! PFT as below:
  !!
  !! \f[
  !! \begin{array} { | l | c | c | c | c | c | }
  !! \hline
  !! \text{CLASS PFTs} &  \text{CTEM PFTs} &     \text{---} &      \text{---}  \\ \hline
  !! \hline
  !! \text{needle leaf} &  \text{evg} &      \text{dcd} &      \text{---}  \\ \hline
  !! \text{broad leaf}  &  \text{evg} &  \text{dcd-cld} &  \text{dcd-dry} \\ \hline
  !! \text{crops}       &   \text{C3} &       \text{C4} &      \text{---}  \\ \hline
  !! \text{grasses}     &   \text{C3} &       \text{C4} &      \text{sedges} \\ \hline
  !! \text{broad leaf shrubs} & \text{evg-shrubs} & \text{dcd-shrubs} &      \text{---}  \\ \hline
  !! \end{array}
  !! \f]
  !!
end module classicParams
