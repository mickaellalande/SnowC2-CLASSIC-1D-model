!> \file
!> Contains the biogeochemistry-related variable type structures.
!! @author J. Melton
!! Variable types herein:
!!
!! 1. c_switch (ctem_switches) - Switches for running the model, read from the joboptions file
!! 2. vrot (veg_rot) - CTEM's 'rot' vars
!! 3. vgat (veg_gat) - CTEM's 'gat' vars
!! 4. ctem_tile (ctem_tile_level) - CTEM's variables per tile
!! 5. ctem_mo (ctem_monthly) - CTEM's variables monthly averaged (per pft)
!! 6. ctem_grd_mo (ctem_gridavg_monthly)  - CTEM's grid average monthly values
!! 7. ctem_tile_mo (ctem_tileavg_monthly) - CTEM's variables per tile monthly values
!! 8. ctem_yr (ctem_annual) - CTEM's average annual values (per PFT)
!! 9. ctem_grd_yr (ctem_gridavg_annual) - CTEM's grid average annual values
!! 10. ctem_tile_yr (ctem_tileavg_annual) - CTEM's variables per tile annual values

module ctemStateVars

  ! S.R.C.    April 2022 - This module was modified to statically allocate these variables thus removing allocCtemVars

  ! J. Melton Apr 2015

  use classicParams,  only : nlat, nmos, ilg, ican, ignd, icp1, icc, iccp2, iccp1, monthend, &
                             mmday, modelpft, l2max, deltat, monthdays, seed, crop, NBS

  implicit none

  public :: initRowVarsBioGeoChem     ! Initializes 'row' variables
  public :: resetMonthEnd   ! Resets monthly variables at month end in preparation for next month
  public :: resetYearEnd    ! Resets annual variables in preparation for next year
  public :: resetMosaicAccum ! Resets physics accumulator variables (used as input to CTEM) after CTEM has been called
  public :: ctemdump ! Dumps all the ctemdump statevars for diagnostic purposes.

  !=================================================================================
  !> Switches for running the model, read from the joboptions file
  type ctem_switches

    logical :: agcm_classic_on  !< True if running CLASSIC within CanESM.
    logical :: projectedGrid    !< True if you have a projected lon lat grid, false if not. Projected grids can only have
    !! regions referenced by the indexes, not coordinates, when running a sub-region
    logical :: ctem_on          !< True if this run includes the biogeochemistry parameterizations (CTEM)
    logical :: Ncycle_on        !< True if this run includes Nitrogen Cycle processes
    integer :: metLoop          !< no. of times the meteorological data is to be looped over. this
    !< option is useful to equilibrate CTEM's C pools
    logical :: leap             !< set to true if the meteorological file includes leap years

    integer :: spinfast         !< set this to a number >1 (up to ~10) to spin up soil carbon pool faster, set to 1
    !< for final production simulations
    integer :: useTracer        !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    character(350) :: tracerCO2file !< Tracer CO2 file, this file needs to be correctly chosen for the useTracer option. It uses the transientCO2
    !! and fixedYearCO2 switches to determine how it is read in.
    integer :: readMetStartYear !< First year of meteorological forcing to read in from the met file
    integer :: readMetEndYear   !< Last year of meteorological forcing to read in from the met file

    logical :: transientCO2     !< use \f$CO_2\f$ time series, if false, fixedYearCO2 is used
    character(350) :: CO2File   !< Location of the netcdf file containing atmospheric \f$CO_2\f$ values
    integer :: fixedYearCO2     !< set the year to use for atmospheric \f$CO_2\f$ if transientCO2 is false. (ppmv)

    logical :: useStaticPeatDep !< True will keep the peat depth at the sdep value. This is used for spinup to ensure
                                !! that the disequilibrium between the peatland and the driving climate, time period 
                                !! since the peatland initiation, or sub-grid heterogeneity impacts don't adversely impact the ability 
                                !! of the model to spinup.
    logical :: doMethane        !< Setting to true will enable wetland CH4 and soil uptake of CH4 subroutines and also require the
                                !! [CH4] input file.
    logical :: transientCH4     !< use \f$CH_4\f$ time series, if false, fixedYearCH4 is used
    character(350) :: CH4File   !< Location of the netcdf file containing atmospheric \f$CH_4\f$ values
    integer :: fixedYearCH4     !< set the year to use for atmospheric \f$CH_4\f$ if transientCH4 is false. (ppmv)

    logical :: dofire           !< boolean, if true allow fire disturbance, if false no fire occurs.
    logical :: transientPOPD    !< if set true use time series of population density data to calculate
    !< fire extinguishing probability and probability of fire due to human causes, false use
    !< value from single year (fixedYearPOPD)
    character(350) :: POPDFile  !< Location of the netcdf file containing population density values
    integer :: fixedYearPOPD    !< set the year to use for population density values if transientPOPD is false. (\f$people/km^2\f$)

    logical :: transientLGHT    !< use lightning strike time series, otherwise use fixedYearLGHT
    character(350) :: LGHTFile  !< Location of the netcdf file containing lightning strike values
    integer :: fixedYearLGHT    !< set the year to use for lightning strikes if transientLGHT is false.

    logical :: PFTCompetition   !< logical boolean telling if competition between pfts is on or not
    logical :: start_bare       !< set this to true if competition is on, and if you wish to start from bare ground.
    !< if this is set to false, the init file info will be used to set up the run.
    !< NOTE: This still keeps the crop fractions (while setting all pools to zero)
    logical :: inibioclim       !< switch telling if bioclimatic parameters are being initialized
    !< from scratch (false) or being initialized from some spun up
    !< values(true).

    logical :: lnduseon         !< If true then the land cover is read in from LUCFile and changes annually
    character(350) :: LUCFile   !< Location of the netcdf file containing land use change information
    integer :: fixedYearLUC     !< Set the year to use for land cover if lnduseon is false. If set to -9999,
    !< we use the PFT distribution found in the initialization file. Any other year
    !< we search for that year in the LUCFile

    logical :: fertilizeron     !< If true then fertilizer is read in from FERFile and changes annually
    character(350) :: FERFile   !< Location of the netcdf file containing fertilizer information
    integer :: fixedYearFER     !< Set the year to use for fertilizer if fertilizeron is false.
    logical :: transientFER     !< use fertilizertime series, if false, fixedYearFER is used

    logical :: depositionon     !< If true then deposition is read in from DEPFile and changes annually
    character(350) :: DEPFile   !< Location of the netcdf file containing deposition information
    integer :: fixedYearDEP     !< Set the year to use for deposition if depositionon is false.
    logical :: transientDEP     !< use depositiontime series, if false, fixedYearDEP is used

    logical :: transientOBSWETF    !< use observed wetland fraction time series, otherwise use fixedYearOBSWETF
    character(350) :: OBSWETFFile  !< Location of the netcdf file containing observed wetland fraction
    integer :: fixedYearOBSWETF    !< set the year to use for observed wetland fraction if transientOBSWETF is false.

    logical :: allLocalTime     !< If your gridded meteorology is relative to Greenwich then set this 
                                !! to false. If the meteorology is all in local time set it to true. To
                                !! determine which, look at your shortwave radiation. As you move through time 
                                !! do you see the sun move across longitudes (allLocalTime = .false.) or across 
                                !! latitudes (allLocalTime = .true.). It seems reanalysis will generally be false 
                                !! while climate model outputs are generally true.
    character(350) :: metFileFss        !< location of the incoming total shortwave radiation meteorology file
    character(350) :: metFilefracFsf        !< location of the incoming diffuse shortwave radiation meteorology file
    character(350) :: metFileFdl        !< location of the incoming longwave radiation meteorology file
    character(350) :: metFilePre        !< location of the precipitation meteorology file
    character(350) :: metFileSnow       !< location of the snow meteorology file
    character(350) :: metFileTa         !< location of the air temperature meteorology file
    character(350) :: metFileQa         !< location of the specific humidity meteorology file
    character(350) :: metFileUv         !< location of the wind speed meteorology file
    character(350) :: metFilePres       !< location of the atmospheric pressure meteorology file
    character(350) :: init_file         !< location of the netcdf initialization file
    character(350) :: rs_file_to_overwrite !< location of the netcdf file that will be written for the restart file
    character(350) :: runparams_file    !< location of the namelist file containing the model parameters
    character(350) :: Comment           !< Comment about the run that will be written to the output netcdfs

    character(350) :: output_directory  !< Directory where the output netcdfs will be placed
    character(350) :: xmlFile           !< location of the xml file that outlines the possible netcdf output files

    logical :: doperpftoutput           !< Switch for making extra output files that are at the per PFT level
    logical :: dopertileoutput          !< Switch for making extra output files that are at the per tile level
    logical :: doChecksums              !< Switch for doing checksum calculations for the run

    logical :: doAnnualOutput           !< Switch for making annual output files
    logical :: doMonthOutput            !< Switch for making monthly output files
    integer :: jmosty                   !< Year to start writing out the monthly output files. If you want to write monthly outputs right
    !< from the start then put in a negative number (like -9999)

    logical :: doDayOutput              !< Switch for making daily output files
    integer :: jdstd                    !< day of the year to start writing the daily output
    integer :: jdendd                   !< day of the year to stop writing the daily output
    integer :: jdsty                    !< simulation year (iyear) to start writing the daily output
    integer :: jdendy                   !< simulation year (iyear) to stop writing the daily output

    logical :: doHhOutput               !< Switch for making half hourly output files
    integer :: jhhstd                   !< day of the year to start writing the half-hourly output
    integer :: jhhendd                  !< day of the year to stop writing the half-hourly output
    integer :: jhhsty                   !< simulation year (iyear) to start writing the half-hourly output
    integer :: jhhendy                  !< simulation year (iyear) to stop writing the half-hourly output

    !< Switch for running the model with dynamic tiling
    !< any routines which use this feature should have a switch nested within this one
    logical :: dynamicTilingOn           !< switch to turn on dynamic tiling
    logical :: tileAgeReset             !< switch that resets tile age to zero at the start of the run
    logical :: trackTileAge             !< switch that tells the model to track tile age durring the run
    logical :: timberHarvest            !< switch to simulate timber harvest with dynamic tiling
    character(350) :: timberHarvestFile !< location of the annual gridcell fraction harvested file
    logical :: prescribedFire            !< switch to simulate prescribed fire with dynamic tiling
    character(350) :: prescribedFireFile !< location of the annual gridcell fraction burned file

    logical :: doPeat                    !< logical switch that tells the model if any tiles contain peatlands

    ! Physics switches:

    integer :: idisp    !< if idisp=0, vegetation displacement heights are ignored,
    !< because the atmospheric model considers these to be part
    !< of the "terrain".
    !< if idisp=1, vegetation displacement heights are calculated.

    integer :: izref    !< if izref=1, the bottom of the atmospheric model is taken
    !< to lie at the ground surface.
    !< if izref=2, the bottom of the atmospheric model is taken
    !< to lie at the local roughness height.

    integer :: islfd    !< if islfd=0, drcoef is called for surface stability corrections
    !< and the original gcm set of screen-level diagnostic calculations
    !< is done.
    !< if islfd=1, drcoef is called for surface stability corrections
    !< and sldiag is called for screen-level diagnostic calculations.
    !< if islfd=2, flxsurfz is called for surface stability corrections
    !< and diasurf is called for screen-level diagnostic calculations.

    integer :: ipcp     !< If a metFileSnow is supplied, this option is *IGNORED*. 
    !< if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 c.
    !< if ipcp=2, a linear partitioning of precipitation betweeen
    !< rainfall and snowfall is done between 0 c and 2 c.
    !< if ipcp=3, rainfall and snowfall are partitioned according to
    !< a polynomial curve between 0 c and 6 c.

    integer :: iwf     !< if iwf=0, only overland flow and baseflow are modelled, and
    !< the ground surface slope is not modelled.
    !< if iwf=n (0<n<4), the watflood calculations of overland flow
    !< and interflow are performed; interflow is drawn from the top
    !< n soil layers.

    integer :: ITC !< itc, itcg and itg are switches to choose the iteration scheme to
    !< be used in calculating the canopy or ground surface temperature
    !< respectively.  if the switch is set to 1, a bisection method is
    !< used; if to 2, the newton-raphson method is used.
    integer :: ITCG !< itc, itcg and itg are switches to choose the iteration scheme to
    !< be used in calculating the canopy or ground surface temperature
    !< respectively.  if the switch is set to 1, a bisection method is
    !< used; if to 2, the newton-raphson method is used.
    integer :: ITG !< itc, itcg and itg are switches to choose the iteration scheme to
    !< be used in calculating the canopy or ground surface temperature
    !< respectively.  if the switch is set to 1, a bisection method is
    !< used; if to 2, the newton-raphson method is used.

    integer :: IPAI !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IHGT !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IALC !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IALS !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IALG !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: isnoalb !< if isnoalb is set to 0, the original two-band snow albedo algorithms are used.
    !< if it is set to 1, the new four-band routines are used.
    character(350) :: alb4BandParamsFile ! Look up table file for the 4-band albedo parameterization
    character(350) :: blackCarbonFile    ! Black carbon concentration file for the 4-band albedo scheme BC parameterization.
    logical :: blackCdepon          ! If true, include consideration of black carbon deposition in snow processes (4-band scheme)
    logical :: blackCtransientDep   ! If true, allow time-varying black carbon deposition flux in snow processes (4-band scheme)
    integer :: fixedYearBCDep       ! set the year to use for black carbon deposition flux if blackCtransientDep is false.
    logical :: KsatScalaron         ! If true, include exponential distribution scalar for saturated hydraulic conductivity (GRKSAT)

  end type ctem_switches

  type (ctem_switches), save, target :: c_switch

  !=================================================================================
  !> CTEM's 'rot' vars
  type veg_rot

    logical, dimension(nlat,nmos,icc) :: pftexist  !< logical array indicating pfts exist(t) or not(f)
    real, dimension(nlat,nmos,icc) :: cc        !< colonization rate
    real, dimension(nlat,nmos,icc) :: mm        !< mortality rate
    integer, dimension(nlat,nmos,icc) :: lfstatus  !< leaf phenology status
    integer, dimension(nlat,nmos,icc) :: pandays   !< days with positive net photosynthesis(an) for use in
    !< the phenology subroutine
    real, dimension(nlat,nmos,icc) :: gleafmas     !< green leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: gleafmas_ns  !< non-structural green leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: gleafmas_s   !< structural green leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: leafns2s     !< carbon flux from non-structural to structural leaf pool, \f$kg c/m^2.day\f$
    real, dimension(nlat,nmos,icc) :: stemns2s     !< carbon flux from non-structural to structural stem pool, \f$kg c/m^2.day\f$
    real, dimension(nlat,nmos,icc) :: rootns2s     !< carbon flux from non-structural to structural root pool, \f$kg c/m^2.day\f$
    real, dimension(nlat,nmos,icc) :: re_alloc_s2l !< amount of nsc reallocated from stem to leaves during leafout
    real, dimension(nlat,nmos,icc) :: re_alloc_r2l !< amount of nsc reallocated from root to leaves during leafout
    real, dimension(nlat,nmos,icc) :: re_alloc_sr2l!< amount of nsc reallocated from stem and root to leaves during leafout
    real, dimension(nlat,nmos,icc) :: bleafmas     !< brown leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemmass     !< stem mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemmass_ns  !< non-structural stem mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemmass_s   !< structural stem mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootmass     !< root mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootmass_ns  !< non-structural root mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootmass_s   !< structural root mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: pstemmass    !< stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(nlat,nmos,icc) :: pgleafmass   !< root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(nlat,nmos,icc) :: fcancmx      !< max. fractional coverage of ctem's pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, dimension(nlat,nmos,icc) :: ngleafmas    !< green leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: ngleafmas_ns !< non-structural green leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: ngleafmas_s  !< structural green leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nbleafmas    !< brown leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nstemmass    !< stem nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nstemmass_ns !<non-structural stem nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nstemmass_s  !<structural stem nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nrootmass    !<root nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nrootmass_ns !<non-structural root nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nrootmass_s  !<structural root nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$

    real, dimension(nlat,nmos,icc) :: ailcg        !< Green LAI for CTEM's pfts
    real, dimension(nlat,nmos,icc) :: ailcgs       !< Green LAI for canopy over snow sub-area
    real, dimension(nlat,nmos,icc) :: fcancs       !< Fraction of canopy over snow for ctem's pfts
    real, dimension(nlat,nmos,icc) :: fcanc        !< Fractional coverage of carbon pfts, canopy over snow
    real, dimension(nlat,nmos,icc) :: co2i1cg      !< Intercellular CO2 conc for pfts for canopy over ground subarea(Pa) - for single/sunlit leaf
    real, dimension(nlat,nmos,icc) :: co2i1cs      !< Same as above but for shaded leaf(above being co2i1cg)
    real, dimension(nlat,nmos,icc) :: co2i2cg      !< Intercellular CO2 conc for pfts for canopy over snowsubarea(pa) - for single/sunlit leaf
    real, dimension(nlat,nmos,icc) :: co2i2cs      !< Same as above but for shaded leaf(above being co2i2cg)
    real, dimension(nlat,nmos,icc) :: ancsveg      !< Net photosynthetic rate for CTEM's pfts for canopy over snow subarea
    real, dimension(nlat,nmos,icc) :: ancgveg      !< Net photosynthetic rate for CTEM's pfts for canopy over ground subarea
    real, dimension(nlat,nmos,icc) :: rmlcsveg     !< Leaf respiration rate for CTEM' pfts forcanopy over snow subarea
    real, dimension(nlat,nmos,icc) :: rmlcgveg     !< Leaf respiration rate for CTEM' pfts forcanopy over ground subarea
    real, dimension(nlat,nmos,icc) :: slai         !< storage/imaginary lai for phenology purposes
    real, dimension(nlat,nmos,icc) :: ailcb        !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, dimension(nlat,nmos,icc) :: flhrloss     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(nlat,nmos,icc) :: flhrloss_ns     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(nlat,nmos,icc) :: flhrloss_s     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(nlat,nmos,icc) :: grwtheff     !< growth efficiency. change in biomass per year per unit max.
    !< lai(\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine
    real, dimension(nlat,nmos,icc) :: lystmmas     !< stem mass at the end of last year
    real, dimension(nlat,nmos,icc) :: lyrotmas     !< root mass at the end of last year
    real, dimension(nlat,nmos,icc) :: tymaxlai     !< this year's maximum lai
    real, dimension(nlat,nmos,icc) :: stmhrlos     !< stem harvest loss for crops, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: vgbiomas_veg !< vegetation biomass for each pft
    real, dimension(nlat,nmos,icc) :: emit_co2     !< carbon dioxide
    real, dimension(nlat,nmos,icc) :: emit_co      !< carbon monoxide
    real, dimension(nlat,nmos,icc) :: emit_ch4     !< methane
    real, dimension(nlat,nmos,icc) :: emit_nmhc    !< non-methane hydrocarbons
    real, dimension(nlat,nmos,icc) :: emit_h2      !< hydrogen gas
    real, dimension(nlat,nmos,icc) :: emit_nox     !< nitrogen oxides
    real, dimension(nlat,nmos,icc) :: emit_n2o     !< nitrous oxide
    real, dimension(nlat,nmos,icc) :: emit_nh3     !< ammonia (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(nlat,nmos,icc) :: emit_pm25    !< particulate matter less than 2.5 um in diameter
    real, dimension(nlat,nmos,icc) :: emit_tpm     !< total particulate matter
    real, dimension(nlat,nmos,icc) :: emit_tc      !< total carbon
    real, dimension(nlat,nmos,icc) :: emit_oc      !< organic carbon
    real, dimension(nlat,nmos,icc) :: emit_bc      !< black carbon
    real, dimension(nlat,nmos,icc) :: burnvegf     !< per PFT fraction burned of that PFT's area
    real, dimension(nlat,nmos,icc) :: smfuncveg    !<
    real, dimension(nlat,nmos,icc) :: bterm        !< biomass term for fire probabilty calc
    real, dimension(nlat,nmos,icc) :: mterm        !< moisture term for fire probabilty calc
    real, dimension(nlat,nmos,icc) :: bmasveg      !< total(gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: veghght      !< vegetation height(meters)
    real, dimension(nlat,nmos,icc) :: rootdpth     !< 99% soil rooting depth(meters)
    !< both veghght & rootdpth can be used as diagnostics to see
    !< how vegetation grows above and below ground, respectively
    real, dimension(nlat,nmos,icc) :: tltrleaf     !< total leaf litter fall rate(u-mol co2/m2.sec)
    real, dimension(nlat,nmos,icc) :: tltrstem     !< total stem litter fall rate(u-mol co2/m2.sec)
    real, dimension(nlat,nmos,icc) :: tltrroot     !< total root litter fall rate(u-mol co2/m2.sec)
    real, dimension(nlat,nmos,icc) :: leaflitr     !< leaf litter fall rate(u-mol co2/m2.sec). this leaf litter
    !< does not include litter generated due to mortality/fire
    real, dimension(nlat,nmos,icc) :: roottemp     !< root temperature, k
    real, dimension(nlat,nmos,icc) :: afrleaf      !< allocation fraction for leaves
    real, dimension(nlat,nmos,icc) :: afrstem      !< allocation fraction for stem
    real, dimension(nlat,nmos,icc) :: afrroot      !< allocation fraction for root
    real, dimension(nlat,nmos,icc) :: wtstatus     !< soil water status used for calculating allocation fractions
    real, dimension(nlat,nmos,icc) :: ltstatus     !< light status used for calculating allocation fractions
    real, dimension(nlat,nmos,icc) :: gppveg       !< ! gross primary productity for each pft
    real, dimension(nlat,nmos,icc) :: vcmax0       !< max. photosynthetic rate at the top of canopy(\f$(mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos,icc) :: nppveg       !< npp for individual pfts, u-mol co2/m2.sec
    real, dimension(nlat,nmos,icc) :: autoresveg   !<
    real, dimension(nlat,nmos,icc) :: rmlvegacc    !<
    real, dimension(nlat,nmos,icc) :: rmsveg       !< stem maintenance resp. rate for each pft
    real, dimension(nlat,nmos,icc) :: rmrveg       !< root maintenance resp. rate for each pft
    real, dimension(nlat,nmos,icc) :: rgveg        !< growth resp. rate for each pft
    real, dimension(nlat,nmos,icc) :: litrfallveg  !< litter fall in for each pft(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos,icc) :: rothrlos     !< root death as crops are harvested, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: pfcancmx     !< previous year's fractional coverages of pfts
    real, dimension(nlat,nmos,icc) :: nfcancmx     !< next year's fractional coverages of pfts
    real, dimension(nlat,nmos,icc) :: anveg        !< net photosynthesis rate for each pft
    real, dimension(nlat,nmos,icc) :: rmlveg       !< leaf maintenance resp. rate for each pft

    real, dimension(nlat,nmos,iccp1) :: bnf_free         !< free-living biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,icc) :: bnf_ant          !< anthropogenic biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,icc) :: bnf_nat          !< natural biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,iccp1) :: bnf_tot          !< total biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,icc) :: nstress          !< N stress
    real, dimension(nlat,nmos,iccp1) :: nitrifveg              !< nitrification flux from nh4_mass to no3_mass pool(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: no_nitveg              !< NO loss through nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: no_denitveg            !< NO loss through denitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: no_nitdenitveg         !< total NO loss from denitrification and nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2o_nitveg             !< N2O loss through nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2o_denitveg           !< N2O loss through denitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2o_nitdenitveg        !< total N2O loss from denitrification and nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2_denitveg            !< N2 loss through denitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: nvolveg                !< nitrogen volatilization(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: nleachveg              !< nitrogen leaching(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: appl_fert              !< applied nitrogen fertilizer \f$(g N m^{-2} cropland day^{-1}\f$
    real, dimension(nlat,nmos,iccp1) :: ndep_nh4               !< deposition influx into the Ammonium pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(nlat,nmos,iccp1) :: ndep_no3               !< deposition influx into the Nitrate pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(nlat,nmos,icc) :: ndemandveg_wp_npp      !< whole plant npp-based nitrogen demand for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_p_nh4       !< passive nh4+ uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_p_no3       !< passive no3- uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_a_actl_nh4  !< actual active nh4+ uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_a_actl_no3  !< actual active no3- uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg             !< total N uptake (active+passive, NH4+NO3) for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nleafns2sveg           !< nitrogen flux from non-structural to structural leaf pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nstemns2sveg           !< nitrogen flux from non-structural to structural stem pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nrootns2sveg           !< nitrogen flux from non-structural to structural root pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nallocveg_l            !< nitrogen allocation to leaves for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nallocveg_s            !< nitrogen allocation to stem for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nallocveg_r            !< nitrogen allocation to root for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nresorpedveg_s         !< resorped N from leaves to be allocated to stem(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nresorpedveg_r         !< resorped N from leaves to be allocated to root(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nre_allocveg_s2l       !< reallocated N from S to L during leaf out period(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nre_allocveg_r2l       !< reallocated N from R to L during leaf out period(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg_l             !< leaf N litterfall(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg_s             !< stem N litterfall(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg_r             !< root N litterfall(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg               !< total N litterfall (\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: gl2bl_grass_nflux      !< N flux from ngleafmas to nbleafmas(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_l               !< simulated C:N ratio for leaves(\f$g C/g N\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_s               !< simulated C:N ratio for stem(\f$g C/g N\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_r               !< simulated C:N ratio for roots(\f$g C/g N\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_wp              !< simulated C:N ratio for the whole plant(\f$g C/g N\f$)
    real, dimension(nlat,nmos,iccp1) :: c2nveg_litr            !< simulated C:N ratio for litter mass(\f$g C/g N\f$)
    real, dimension(nlat,nmos,iccp1) :: c2nveg_humus           !< simulated C:N ratio for humus mass(\f$g C/g N\f$)
    real, dimension(nlat,nmos,iccp1) :: nhumtrsveg             !< N humification for individual PFTs + bare(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nmineralveg_litr       !< N mineralization from litter pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nmineralveg_humus      !< N mineralization from organic soil pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nimmobilveg_nh4        !< N immobilization from NH4+ pool to soilnmas(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nimmobilveg_no3        !< N immobilization from NO3- pool to soilnmas(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: netnmineralveg         !< N mineralization (\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nvgbiomas_veg          !< N vegetation biomass for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: fNnetlandveg           !< net terrestrial N flux(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: redcoeff_vcmax         !< Reduction coeff. passed to the Photosynthesis subroutine

    ! allocated with nlat,nmos:
    real, dimension(nlat,nmos) :: gavglai               !< grid averaged green leaf area index
    real, dimension(nlat,nmos) :: co2conc               !< ATMOS. CO2 CONC. IN PPM
    real, dimension(nlat,nmos) :: ch4conc               !<
    real, dimension(nlat,nmos) :: canres                !<
    real, dimension(nlat,nmos) :: vgbiomas              !< grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, dimension(nlat,nmos) :: gavgltms              !< grid averaged litter mass, \f$kg c/m^2\f$
    real, dimension(nlat,nmos) :: gavgscms              !< grid averaged soil c mass, \f$kg c/m^2\f$
    real, dimension(nlat,nmos) :: burnfrac              !< areal :: fraction burned due to fire for every grid cell(%)
    real, dimension(nlat,nmos) :: popdin                !< population density \f$(people / km^2)\f$
    real, dimension(nlat,nmos) :: soilpH                !< soil PH
    real, dimension(nlat,nmos) :: nfertil               !< nitrogen fertilizer \f$(g N m^{-2} cropland yr^{-1}\f$
    real, dimension(nlat,nmos) :: ndeposit              !< nitrogen deposition \f$(g N m^{-2} yr^{-1}\f$
    real, dimension(nlat) :: timharvrow                 !< the annual fractional area where timber was harvested (0-1,read in from the external forcing file)
    real, dimension(nlat,nmos) :: timharvarearow        !< the fractional area of the cell where timber will be harvested (0-1,this is the actual area harvested per tile by the model, which is determined from timharvrow by the harvest subroutines)
    real, dimension(nlat) :: prsfirerow                 !< the annual fractional area where fire burned (0-1, read in from the external forcing file)
    real, dimension(nlat,nmos) :: prsfirearearow        !< the fractional area of the cell where fire burned (0-1,this is the actual area burned per tile by the model, which is determined from prsfirerow by the prescribed fire subroutines)
    real, dimension(nlat,nmos) :: tileAgerow            !< the age of the tile since the start of the run in months this is reset by harvest and fire and incremented at the CTEM timestep by indexTileAge
    real, dimension(nlat,nmos) :: lterm                 !< lightning term for fire probabilty calc
    real, dimension(nlat,nmos) :: extnprob              !< fire extingusinging probability
    real, dimension(nlat,nmos) :: prbfrhuc              !< probability of fire due to human causes
    real, dimension(nlat,nmos) :: rml                   !< leaf maintenance respiration(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: rms                   !< stem maintenance respiration(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: rmr                   !< root maintenance respiration(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: ch4WetSpec            !< methane flux from wetlands calculated using hetrores in umol ch4/m2.s
    real, dimension(nlat,nmos) :: wetfdyn               !< dynamic wetland fraction
    real, dimension(nlat,nmos) :: wetfrac_pres          !< Prescribed wetland fraction read in from OBSWETFFile
    real, dimension(nlat,nmos) :: ch4WetDyn                !< methane flux from wetlands calculated using hetrores and wetfdyn, in umol ch4/m2.s
    real, dimension(nlat,nmos) :: ch4_soills            !< Methane uptake into the soil column(\f$mg CH_4 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: lucemcom              !< land use change(luc) related combustion emission losses(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: lucltrin              !< luc related inputs to litter pool(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: lucsocin              !< luc related inputs to soil c pool(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: lucemcomn   !< land use change (luc) related combustion emission losses of N (g N m^-2 day^-1)
    real, dimension(nlat,nmos) :: lucltrinn   !< luc related inputs to litter pool of N (g N m^-2 day^-1)
    real, dimension(nlat,nmos) :: lucsocinn   !< luc related inputs to soil N pool of N (g N m^-2 day^-1)
    real, dimension(nlat,nmos) :: npp                   !< net primary productivity
    real, dimension(nlat,nmos) :: nep                   !< net ecosystem productivity
    real, dimension(nlat,nmos) :: nepCMIP               !< net ecosystem productivity 
    real, dimension(nlat,nmos) :: nbp                   !< net biome productivity
    real, dimension(nlat,nmos) :: gpp                   !< gross primary productivity
    real, dimension(nlat,nmos) :: hetrores              !< heterotrophic respiration
    real, dimension(nlat,nmos) :: autores               !< autotrophic respiration
    real, dimension(nlat,nmos) :: soilcresp             !<
    real, dimension(nlat,nmos) :: rm                    !< maintenance respiration
    real, dimension(nlat,nmos) :: rg                    !< growth respiration
    real, dimension(nlat,nmos) :: litres                !< litter respiration
    real, dimension(nlat,nmos) :: socres                !< soil carbon respiration
    real, dimension(nlat,nmos) :: dstcemls              !< carbon emission losses due to disturbance, mainly fire
    real, dimension(nlat,nmos) :: litrfall              !< total litter fall(from leaves, stem, and root) due to
    !< all causes(mortality,turnover,and disturbance)(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: humiftrs              !< transfer of humidified litter from litter to soil c pool
    real, dimension(nlat,nmos) :: cfluxcg               !<
    real, dimension(nlat,nmos) :: cfluxcs               !<
    real, dimension(nlat,nmos) :: ROFB                  !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: dstcemls3             !< carbon emission losses due to disturbance(fire at present) from litter pool
    real, dimension(nlat,nmos) :: uvaccrow_m            !<
    real, dimension(nlat,nmos) :: vvaccrow_m            !<
    real, dimension(nlat,nmos) :: qevpacc_m_save        !<
    real, dimension(nlat) :: twarmm                !< temperature of the warmest month(c)
    real, dimension(nlat) :: tcoldm                !< temperature of the coldest month(c)
    real, dimension(nlat) :: gdd5                  !< growing degree days above 5 c
    real, dimension(nlat) :: aridity               !< aridity index, ratio of potential evaporation to precipitation
    real, dimension(nlat) :: srplsmon              !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, dimension(nlat) :: defctmon              !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, dimension(nlat) :: anndefct              !< annual water deficit(mm)
    real, dimension(nlat) :: annsrpls              !< annual water surplus(mm)
    real, dimension(nlat) :: annpcp                !< annual precipitation(mm)
    real, dimension(nlat) :: dry_season_length     !< length of dry season(months)

    integer, dimension(nlat,nmos) :: ipeatland !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    real, dimension(nlat,nmos) :: litrmsmoss
    real, dimension(nlat,nmos) :: Cmossmas
    real, dimension(nlat,nmos) :: dmoss
    real, dimension(nlat,nmos) :: peatSoilC         !< peat soil C mass, \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: nppmoss
    real, dimension(nlat,nmos) :: rmlmoss
    real, dimension(nlat,nmos) :: gppmoss
    real, dimension(nlat,nmos) :: anmoss
    real, dimension(nlat,nmos) :: armoss
    real, dimension(nlat,nmos) :: peatdep
    real, dimension(nlat,nmos) :: pdd
    integer, dimension(nlat,nmos) :: colddays_leaffall          !< cold days counter for tracking days below a certain
    !< temperature threshold for ndl dcd tree.
    integer, dimension(nlat,nmos) :: colddays_harvest          !< cold days counter for tracking days below a certain
    !< temperature threshold for crops.

    ! allocated with nlat,nmos,ican:
    real, dimension(nlat,nmos,ican) :: zolnc            !< lumped log of roughness length for class' 4 pfts
    real, dimension(nlat,nmos,ican) :: ailc             !< lumped lai for class' 4 pfts
    real, dimension(nlat,nmos,ican) :: cmasvegc         !< total canopy mass for each of the 4 class pfts. recall that
    !< class requires canopy mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, dimension(nlat,nmos,ican) :: alvsctm          !<
    real, dimension(nlat,nmos,ican) :: paic             !< plant area index for class' 4 pfts. this is the sum of leaf
    !< area index and stem area index.
    real, dimension(nlat,nmos,ican) :: slaic            !< storage lai. this will be used as min. lai that class sees
    !< so that it doesn't blow up in its stomatal conductance calculations.
    real, dimension(nlat,nmos,ican) :: alirctm          !<

    ! allocated with nlat,nmos,ican,ignd:
    real, dimension(nlat,nmos,ican,ignd) :: rmatc       !< fraction of roots for each of class' 4 pfts in each soil layer

    ! allocated with nlat,nmos,icc,ignd:
    real, dimension(nlat,nmos,icc,ignd) :: rmatctem     !< fraction of roots for each of ctem's 9 pfts in each soil layer

    ! allocated with nlat,nmos,iccp1:
    real, dimension(nlat,nmos,iccp1) :: nepveg      !< net ecosystem productity for each pft
    real, dimension(nlat,nmos,iccp1) :: nbpveg      !< net biome productity for bare fraction OR net biome productity for each pft
    real, dimension(nlat,nmos,iccp1) :: hetroresveg !<

    ! allocated with nlat,nmos,iccp2,ignd:
    real, dimension(nlat,nmos,iccp2,ignd) :: litrmass    !< litter mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcmas    !< soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,iccp2,ignd) :: litresveg   !<
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcresveg !<
    real, dimension(nlat,nmos,iccp2,ignd) :: humiftrsveg !< Transfer of humidified litter from litter to soil C pool 

    real, dimension(nlat,nmos,iccp1) :: nh4_mass    !< ammonium mass for individual PFTs + bare, \f$g N/m^2\f$
    real, dimension(nlat,nmos,iccp1) :: no3_mass    !< nitrate mass for individual PFTs + bare, \f$g N/m^2\f$
    real, dimension(nlat,nmos,iccp2) :: nlitrmass   !< litter nitrogen mass for individual PFTs + bare, \f$g N/m^2\f$
    real, dimension(nlat,nmos,iccp2) :: soilnmas    !< soil organic nitrogen mass for individual PFTs + bare, \f$g N/m^2\f$

    ! allocated with nlat,nmos,{some number}:
    real, dimension(nlat,nmos,8)  :: slopefrac          !< prescribed fraction of wetlands based on slope
    !< only(0.025,0.05,0.1,0.15,0.20,0.25,0.3 and 0.35 percent slope thresholds)

    ! allocated with nlat:
    real, dimension(nlat)    :: dayl_max        !< maximum daylength for that location(hours)
    real, dimension(nlat)    :: dayl            !< daylength for that location(hours)
    real, dimension(nlat)    :: grclarearow     !< area of the grid cell, \f$km^2\f$

    ! allocated with nlat,nmos,icc:
    real, dimension(nlat,nmos,icc) :: lygleafmasmax        !< last year maximum of the gleafmas
    real, dimension(nlat,nmos,icc) :: lystemmassmax        !< last year maximum of the stemmass
    real, dimension(nlat,nmos,icc) :: lyrootmassmax        !< last year maximum of the rootmass
    real, dimension(nlat,nmos,icc) :: lmaxt, smaxt, rmaxt  !< temp vars to find previous year C pool max

    !input variables that control the dynamic tiling subtoutine
    integer, dimension(nlat,nmos) :: controlVector(nlat,nmos) !< this is the control vector which is specified in the initializtion file
    !integer, dimension(nlat,nmos) :: inputVector !< this is the input vector which is generated by the script (1 = operate on tile, 0 = do not operate on tile) 
    !character(len=5), dimension(nlat) :: mode !< mode is a character string that specifies the operation (i.e. move, split or copy) 
    !real, dimension(nlat) :: FAREAROUT !< a real value greater than zero and less than one which specifies the desired size of the output tile
    logical :: DynTilInitializeFlag !< this flag is used to re-initialize the model after dynamic tiling runs

  end type veg_rot

  type (veg_rot),save,target :: vrot
  !$omp threadprivate(vrot)

  !=================================================================================
  !> CTEM's 'gat' vars
  type veg_gat

    ! This is the basic data structure that contains the state variables
    ! for the Plant Functional Type (PFT). The dimensions are ilg,{icc,iccp1,iccp2}

    real, dimension(ilg,icc) :: gleafmas   !< Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: gleafmas_ns  !< non-structural green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$

    
    real, dimension(ilg,icc) :: gleafmas_s   !< structural green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: leafns2s     !< carbon flux from non-structural to structural leaf pool, \f$kg c/m^2.day\f$
    real, dimension(ilg,icc) :: stemns2s     !< carbon flux from non-structural to structural stem pool, \f$kg c/m^2.day\f$
    real, dimension(ilg,icc) :: rootns2s     !< carbon flux from non-structural to structural root pool, \f$kg c/m^2.day\f$
    real, dimension(ilg,icc) :: re_alloc_s2l !< amount of nsc reallocated from stem to leaves during leaf out, \f$g c/m^2\f$.
    real, dimension(ilg,icc) :: re_alloc_r2l !< amount of nsc reallocated from root to leaves during leaf out, \f$g c/m^2\f$.
    real, dimension(ilg,icc) :: re_alloc_sr2l!< amount of nsc reallocated from stem and root to leaves during leaf out, \f$g c/m^2\f$.
    real, dimension(ilg,icc) :: bleafmas   !< Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stemmass   !< Stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stemmass_ns  !< non-structural stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stemmass_s   !< structural stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rootmass   !< Root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rootmass_ns  !< non-structural root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rootmass_s   !< structural root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: pstemmass  !< Stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc) :: pgleafmass !< Root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc) :: fcancmx    !< max. fractional coverage of ctem's 9 pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, dimension(ilg,icc) :: ngleafmas    !< green leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: ngleafmas_ns !< non-structural green leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: ngleafmas_s  !< structural green leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nbleafmas    !< brown leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nstemmass    !< stem nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nstemmass_ns !< non-structural stem nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nstemmass_s  !< structural stem nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nrootmass    !< root nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nrootmass_ns !< non-structural nitrogen root mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nrootmass_s  !< structural root nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nvgbiomas_veg!< N vegetation biomass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg) :: gavglai        !< grid averaged green leaf area index

    real, dimension(ilg) :: lightng        !< total lightning frequency, flashes/km2.year

    real, dimension(ilg,ican) :: zolnc     !< lumped log of roughness length for class' 4 pfts
    real, dimension(ilg,ican) :: ailc      !< lumped lai for class' 4 pfts

    real, dimension(ilg,icc) :: ailcg      !< green lai for ctem's 9 pfts
    real, dimension(ilg,icc) :: ailcgs     !< GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
    real, dimension(ilg,icc) :: fcancs     !< FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
    real, dimension(ilg,icc) :: fcanc      !< FRACTIONAL COVERAGE OF 8 CARBON PFTs, CANOPY OVER GROUND

    real, dimension(ilg)   :: co2conc    !< ATMOS. CO2 CONC. IN PPM
    real, dimension(ilg)   :: ch4conc    !<

    real, dimension(ilg,icc) :: co2i1cg    !< INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER GROUND SUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, dimension(ilg,icc) :: co2i1cs    !< SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i1cg)
    real, dimension(ilg,icc) :: co2i2cg    !< INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER SNOWSUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, dimension(ilg,icc) :: co2i2cs    !< SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i2cg)
    real, dimension(ilg,icc) :: ancsveg    !< net photosynthetic rate for ctems 9 pfts for canopy over snow subarea
    real, dimension(ilg,icc) :: ancgveg    !< net photosynthetic rate for ctems 9 pfts for canopy over ground subarea
    real, dimension(ilg,icc) :: rmlcsveg   !< leaf respiration rate for ctems 9 pfts forcanopy over snow subarea
    real, dimension(ilg,icc) :: rmlcgveg   !< leaf respiration rate for ctems 9 pfts forcanopy over ground subarea
    real, dimension(ilg,icc) :: slai       !< storage/imaginary lai for phenology purposes
    real, dimension(ilg,icc) :: ailcb      !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, dimension(ilg)   :: canres     !<
    real, dimension(ilg,icc) :: flhrloss   !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(ilg,icc) :: flhrloss_ns     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(ilg,icc) :: flhrloss_s     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(ilg,icc) :: grwtheff   !< growth efficiency. change in biomass per year per unit max.
    !< lai (\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine
    real, dimension(ilg,icc) :: lystmmas   !< stem mass at the end of last year
    real, dimension(ilg,icc) :: lyrotmas   !< root mass at the end of last year
    real, dimension(ilg,icc) :: tymaxlai   !< this year's maximum lai
    real, dimension(ilg)   :: vgbiomas   !< grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, dimension(ilg)   :: gavgltms   !< grid averaged litter mass, \f$kg c/m^2\f$
    real, dimension(ilg)   :: gavgscms   !< grid averaged soil c mass, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stmhrlos   !< stem harvest loss for crops, \f$kg c/m^2\f$
    real, dimension(ilg,ican,ignd) :: rmatc    !< fraction of roots for each of class' 4 pfts in each soil layer
    real, dimension(ilg,icc,ignd) :: rmatctem !< fraction of roots for each of ctem's 9 pfts in each soil layer
    real, dimension(ilg,iccp2,ignd) :: litrmass   !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    real, dimension(ilg,iccp2,ignd) :: soilcmas   !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: vgbiomas_veg !< vegetation biomass for each pft
    real, dimension(ilg,iccp1) :: nh4_mass   !< ammonium mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$
    real, dimension(ilg,iccp1) :: no3_mass   !< nitrate mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$
    real, dimension(ilg,iccp2) :: nlitrmass  !< litter nitrogen mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$
    real, dimension(ilg,iccp2) :: soilnmas   !< soil nitrogen mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$

    real, dimension(ilg,icc) :: emit_co2   !< carbon dioxide (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_co    !< carbon monoxide (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_ch4   !< methane (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_nmhc  !< non-methane hydrocarbons (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_h2    !< hydrogen gas (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_nox   !< nitrogen oxides (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_n2o   !< nitrous oxide (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_nh3   !< ammonia (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_pm25  !< particulate matter less than 2.5 um in diameter (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_tpm   !< total particulate matter (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_tc    !< total carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_oc    !< organic carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_bc    !< black carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg)   :: burnfrac   !< areal :: fraction burned due to fire for every grid cell (%)
    real, dimension(ilg,icc) :: burnvegf   !< per PFT fraction burned of that PFT's area
    real, dimension(ilg,icc) :: smfuncveg  !<
    real, dimension(ilg)   :: popdin     !< population density (people / \f$km^2\f$)
    real, dimension(ilg)   :: soilpH     !< soil PH
    real, dimension(ilg)   :: nfertil    !< nitrogen fertilizer \f$(g N m^{-2} cropland yr^{-1}\f$
    real, dimension(ilg)   :: ndeposit   !< nitrogen deposition \f$(g N m^{-2} yr^{-1}\f$
    real, dimension(ilg) :: timharvareagat   !< the same as timharvarea, but in CTEM's 'gat' format 
    real, dimension(ilg) :: tileAgegat       !< the age of the tile since the start of the run in months this is reset by harvest and fire and incremented at the CTEM timestep by indexTileAge
    real, dimension(ilg) :: prsfireareagat   !< the same as prsfirearea, but in CTEM's 'gat' format
    real, dimension(ilg,icc) :: bterm      !< biomass term for fire probabilty calc
    real, dimension(ilg)   :: lterm      !< lightning term for fire probabilty calc
    real, dimension(ilg,icc) :: mterm      !< moisture term for fire probabilty calc
    real, dimension(ilg,icc) :: glcaemls  !< green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: blcaemls  !< brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rtcaemls  !< root carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stcaemls  !< stem carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: ltrcemls  !< litter carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: ntchlveg  !< fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
    real, dimension(ilg,icc) :: ntchsveg  !< fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
    real, dimension(ilg,icc) :: ntchrveg  !< fluxes for each pft: Net change in root biomass,
    !! the net change is the difference between allocation and
    !! autotrophic respiratory fluxes, u-mol CO2/m2.sec

    real, dimension(ilg)   :: extnprob   !< fire extingusinging probability
    real, dimension(ilg)   :: prbfrhuc   !< probability of fire due to human causes
    real, dimension(ilg)   :: dayl_max   !< maximum daylength for that location (hours)
    real, dimension(ilg)   :: dayl       !< daylength for that location (hours)

    real, dimension(ilg,icc) :: bmasveg    !< total (gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, dimension(ilg,ican) :: cmasvegc   !< total canopy mass for each of the 4 class pfts. recall that
    !< class requires canopy mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, dimension(ilg,icc) :: veghght    !< vegetation height (meters)
    real, dimension(ilg,icc) :: rootdpth   !< 99% soil rooting depth (meters)
    !< both veghght & rootdpth can be used as diagnostics to see
    !< how vegetation grows above and below ground, respectively
    real, dimension(ilg)   :: rml        !< leaf maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: rms        !< stem maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: tltrleaf   !< total leaf litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: blfltrdt !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, dimension(ilg,icc) :: glfltrdt !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, dimension(ilg,icc) :: tltrstem   !< total stem litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: tltrroot   !< total root litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: leaflitr   !< leaf litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$). this leaf litter
    !< does not include litter generated due to mortality/fire
    real, dimension(ilg,icc) :: roottemp   !< root temperature, k
    real, dimension(ilg,icc) :: afrleaf    !< allocation fraction for leaves
    real, dimension(ilg,icc) :: afrstem    !< allocation fraction for stem
    real, dimension(ilg,icc) :: afrroot    !< allocation fraction for root
    real, dimension(ilg,icc) :: wtstatus   !< soil water status used for calculating allocation fractions
    real, dimension(ilg,icc) :: ltstatus   !< light status used for calculating allocation fractions
    real, dimension(ilg)   :: rmr        !< root maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,8) :: slopefrac      !< prescribed fraction of wetlands based on slope
    !< only(0.025,0.05,0.1,0.15,0.20,0.25,0.3 and 0.35 percent slope thresholds)
    real, dimension(ilg)   :: wetfrac_pres  !< Prescribed fraction of wetlands in a grid cell
    real, dimension(ilg)   :: ch4WetSpec       !< methane flux from wetlands calculated using hetrores (\f$\mu mol CH_4 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: wetfdyn       !< dynamic wetland fraction
    real, dimension(ilg)   :: ch4WetDyn       !< methane flux from wetlands calculated using hetrores
    !< and wetfdyn, (\f$\mu mol CH_4 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: ch4_soills    !< Methane uptake into the soil column (\f$mg CH_4 m^{-2} s^{-1}\f$)

    real, dimension(ilg)   :: lucemcom   !< land use change (luc) related combustion emission losses, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: lucltrin   !< luc related inputs to litter pool, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: lucsocin   !< luc related inputs to soil c pool, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: lucemcomn   !< land use change (luc) related combustion emission losses of N (g N m^-2 day^-1)
    real, dimension(ilg)   :: lucltrinn   !< luc related inputs to litter pool of N (g N m^-2 day^-1)
    real, dimension(ilg)   :: lucsocinn   !< luc related inputs to soil N pool of N (g N m^-2 day^-1)

    real, dimension(ilg)   :: npp        !< net primary productivity
    real, dimension(ilg)   :: nep        !< net ecosystem productivity
    real, dimension(ilg)   :: nepCMIP        !< net ecosystem productivity
    real, dimension(ilg)   :: nbp        !< net biome productivity
    real, dimension(ilg)   :: gpp        !< gross primary productivity
    real, dimension(ilg)   :: hetrores   !< heterotrophic respiration
    real, dimension(ilg)   :: autores    !< autotrophic respiration
    real, dimension(ilg)   :: soilcresp  !<
    real, dimension(ilg)   :: rm         !< maintenance respiration
    real, dimension(ilg)   :: rg         !< growth respiration
    real, dimension(ilg)   :: litres     !< litter respiration
    real, dimension(ilg)   :: socres     !< soil carbon respiration
    real, dimension(ilg)   :: dstcemls   !< carbon emission losses due to disturbance, mainly fire
    real, dimension(ilg)   :: litrfall   !< total litter fall (from leaves, stem, and root) due to
    !< all causes (mortality,turnover,and disturbance)(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: humiftrs   !< transfer of humidified litter from litter to soil c pool

    real, dimension(ilg,icc) :: gppveg     !< gross primary productity for each pft
    real, dimension(ilg,icc) :: vcmax0     !< max. photosynthetic rate at the top of canopy (\f$(mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,iccp1) :: nepveg     !< net ecosystem productity for bare fraction expnbaln(i)=0.0 amount
    !< of c related to spatial expansion Not used JM Jun 2014
    !< OR net ecosystem productity for each pft

    real, dimension(ilg,iccp1) :: bnf_free         !< free-living biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,icc) :: bnf_ant          !< anthropogenic biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,icc) :: bnf_nat          !< natural biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,iccp1) :: bnf_tot          !< total biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,icc) :: nstress          !< N stress
    real, dimension(ilg,iccp1) :: nitrifveg             !< nitrification for individual PFTs + bareground (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: no_nitveg             !< NO loss through nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: no_denitveg           !< NO loss through denitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: no_nitdenitveg        !< total NO loss from denitrification and nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2o_nitveg            !< N2O loss through nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2o_denitveg          !< N2O loss through denitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2o_nitdenitveg       !< total N2O loss from denitrification and nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2_denitveg           !< N2 loss through denitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: nvolveg               !< volatilization loss from ammonium pool (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: nleachveg             !< leaching from nitrate pool (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: appl_fert             !< applied nitrogen fertilizer \f$(g N m^{-2} cropland day^{-1}\f$
    real, dimension(ilg,iccp1) :: ndep_nh4              !< deposition influx into the Ammonium pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(ilg,iccp1) :: ndep_no3              !< deposition influx into the Nitrate pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(ilg,icc) :: ndemandveg_wp_npp     !< whole plant npp-based nitrogen demand for each of the 9 ctem pfts (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_p_nh4      !< passive nh4+ uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_p_no3      !< passive no3- uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_a_actl_nh4 !< actual active nh4+ uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_a_actl_no3 !< actual active no3- uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg            !< total N uptake (active+passive, NH4+NO3) for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nleafns2sveg          !< nitrogen flux from non-structural to structural leaf pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nstemns2sveg          !< nitrogen flux from non-structural to structural stem pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nrootns2sveg          !< nitrogen flux from non-structural to structural root pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nallocveg_l           !< nitrogen allocation to leaves for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nallocveg_s           !< nitrogen allocation to stem for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nallocveg_r           !< nitrogen allocation to root for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nresorpedveg_s        !< resorped N from leaves to be allocated to stem (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nresorpedveg_r        !< resorped N from leaves to be allocated to root (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nre_allocveg_s2l      !< reallocated N from S to L during leaf out period (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nre_allocveg_r2l      !< reallocated N from R to L during leaf out period (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg_l            !< leaf N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg_s            !< stem N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg_r            !< root N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg              !< total N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: gl2bl_grass_nflux     !< N flux from ngleafmas to nbleafmas (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: c2nveg_l              !< simulated C:N ratio for leaves (\f$g C/g N\f$)
    real, dimension(ilg,icc) :: c2nveg_s              !< simulated C:N ratio for stem (\f$g C/g N\f$)
    real, dimension(ilg,icc) :: c2nveg_r              !< simulated C:N ratio for roots (\f$g C/g N\f$)
    real, dimension(ilg,icc) :: c2nveg_wp             !< simulated C:N ratio for the whole plant (\f$g C/g N\f$)
    real, dimension(ilg,iccp1) :: c2nveg_litr           !< simulated C:N ratio for litter mass (\f$g C/g N\f$)
    real, dimension(ilg,iccp1) :: c2nveg_humus          !< simulated C:N ratio for humus mass (\f$g C/g N\f$)
    real, dimension(ilg,iccp1) :: nhumtrsveg            !< N humification for individual PFTs + bare (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nmineralveg_litr      !< N mineralization from litter pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nmineralveg_humus     !< N mineralization from organic soil pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: netnmineralveg        !< N mineralization (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nimmobilveg_nh4       !< N immobilization from NH4+ pool to soilnmas (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nimmobilveg_no3       !< N immobilization from NO3- pool to soilnmas (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: fNnetlandveg          !< net terrestrial N flux (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: redcoeff_vcmax        !< Reduction coeff. passed to the Photosynthesis subroutine

    integer, dimension(ilg) :: ipeatland !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    real, dimension(ilg) :: peatdep      !< Depth of peat column (m)
    real, dimension(ilg) :: anmoss     !< net photosynthetic rate of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlmoss    !< maintenance respiration rate of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: gppmoss    !< gross primaray production of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: nppmoss    !< net primary production of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: armoss     !< autotrophic respiration of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: litrmsmoss !< moss litter mass, \f$kg C/m^2\f$
    real, dimension(ilg) :: Cmossmas   !< C in moss biomass, \f$kg C/m^2\f$
    real, dimension(ilg) :: dmoss      !< depth of living moss (m)
    real, dimension(ilg) :: peatSoilC  !< peat soil C mass, \f$kg C/m^2\f$
    real, dimension(ilg) :: pdd        !< peatland degree days above 0 deg C.
    real, dimension(ilg) :: ancsmoss   !< moss net photosynthesis in canopy snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: angsmoss   !< moss net photosynthesis in snow ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: ancmoss    !< moss net photosynthesis in canopy ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: angmoss    !< moss net photosynthesis in bare ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlcsmoss  !< moss maintenance respiration in canopy snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlgsmoss  !< moss maintenance respiration in ground snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlcmoss   !< moss maintenance respiration in canopy ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlgmoss   !< moss maintenance respiration in bare ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,iccp1) :: nbpveg     !< net biome productity for bare fraction OR net biome productity for each pft
    real, dimension(ilg,icc) :: nppveg     !< npp for individual pfts, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,iccp1) :: hetroresveg !<
    real, dimension(ilg,icc) :: autoresveg !<
    real, dimension(ilg,iccp2,ignd) :: litresveg  !<
    real, dimension(ilg,iccp2,ignd) :: soilcresveg !<
    real, dimension(ilg,iccp2,ignd) :: humiftrsveg !<
    real, dimension(ilg,icc) :: rmlvegacc  !<
    real, dimension(ilg,icc) :: rmsveg     !< stem maintenance resp. rate for each pft
    real, dimension(ilg,icc) :: rmrveg     !< root maintenance resp. rate for each pft
    real, dimension(ilg,icc) :: rgveg      !< growth resp. rate for each pft
    real, dimension(ilg,icc) :: litrfallveg !< litter fall for each pft (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: reprocost   !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,icc) :: rothrlos !< root death as crops are harvested, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: pfcancmx !< previous year's fractional coverages of pfts
    real, dimension(ilg,icc) :: nfcancmx !< next year's fractional coverages of pfts
    real, dimension(ilg,ican) :: alvsctm  !<
    real, dimension(ilg,ican) :: paic     !< plant area index for class' 4 pfts. this is the sum of leaf
    !< area index and stem area index.
    real, dimension(ilg,ican) :: slaic    !< storage lai. this will be used as min. lai that class sees
    !< so that it doesn't blow up in its stomatal conductance calculations.
    real, dimension(ilg,ican) :: alirctm  !<
    real, dimension(ilg)   :: cfluxcg  !<
    real, dimension(ilg)   :: cfluxcs  !<
    real, dimension(ilg)   :: CFLUX_GA !< Product of surface drag coefficient and wind speed \f$[m s^{-1} ]\f$
    real, dimension(ilg)   :: USTARBS_GA !< Friction velocity to be used in nitrogen volatilization  \f$[m s^{-1} ]\f$
    real, dimension(ilg)   :: ROFB     !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg)   :: dstcemls3 !< carbon emission losses due to disturbance (fire at present) from litter pool
    real, dimension(ilg,icc) :: anveg    !< net photosynthesis rate for each pft
    real, dimension(ilg,icc) :: rmlveg   !< leaf maintenance resp. rate for each pft

    real, dimension(ilg) :: twarmm            !< temperature of the warmest month (c)
    real, dimension(ilg) :: tcoldm            !< temperature of the coldest month (c)
    real, dimension(ilg) :: gdd5              !< growing degree days above 5 c
    real, dimension(ilg) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
    real, dimension(ilg) :: srplsmon          !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, dimension(ilg) :: defctmon          !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, dimension(ilg) :: anndefct          !< annual water deficit (mm)
    real, dimension(ilg) :: annsrpls          !< annual water surplus (mm)
    real, dimension(ilg) :: annpcp            !< annual precipitation (mm)
    real, dimension(ilg) :: dry_season_length !< length of dry season (months)
    integer, dimension(ilg) :: colddays_leaffall !< cold days counter for tracking days below a certain
    !< temperature threshold for ndl dcd.
    integer, dimension(ilg) :: colddays_harvest !< cold days counter for tracking days below a certain
    !< temperature threshold for crop.

    ! These go into CTEM and are used to keep track of the bioclim limits.
    real, dimension(ilg) :: tcurm     !< temperature of the current month (c)
    real, dimension(ilg) :: srpcuryr  !< water surplus for the current year
    real, dimension(ilg) :: dftcuryr  !< water deficit for the current year
    real, dimension(12,ilg) :: tmonth  !< monthly temperatures
    real, dimension(ilg) :: anpcpcur  !< annual precipitation for current year (mm)
    real, dimension(ilg) :: anpecur   !< annual potential evaporation for current year (mm)
    real, dimension(ilg) :: gdd5cur   !< growing degree days above 5 c for current year
    real, dimension(ilg) :: surmncur  !< number of months with surplus water for current year
    real, dimension(ilg) :: defmncur  !< number of months with water deficit for current year
    real, dimension(ilg) :: srplscur  !< water surplus for the current month
    real, dimension(ilg) :: defctcur  !< water deficit for the current month

    real, dimension(ilg,icc) :: geremort !< growth efficiency related mortality (1/day)
    real, dimension(ilg,icc) :: intrmort !< intrinsic (age related) mortality (1/day)
    real, dimension(ilg,icc) :: cc       !< colonization rate & mortality rate
    real, dimension(ilg,icc) :: mm       !< colonization rate & mortality rate

    logical, dimension(ilg,icc) :: pftexist !< logical array indicating pfts exist (t) or not (f)
    integer, dimension(ilg,icc) :: lfstatus !< leaf phenology status
    integer, dimension(ilg,icc) :: pandays  !< days with positive net photosynthesis (an) for use in
    !< the phenology subroutine
    real, dimension(ilg) :: grclarea      !< area of the grid cell, \f$km^2\f$

    integer, dimension(ilg) :: altotcount_ctm ! nlat  !< Counter used for calculating total albedo
    real, dimension(ilg,icc)  :: todfrac  !(ilg,icc)   !< Max. fractional coverage of ctem's 9 pfts by the end of the day, for use by land use subroutine
    real, dimension(ilg)    :: fsinacc_gat !(ilg)    !<
    real, dimension(ilg)    :: flutacc_gat !(ilg)    !<
    real, dimension(ilg)    :: flinacc_gat !(ilg)    !<
    ! real, dimension(ilg)    :: pregacc_gat !(ilg)    !<
    real, dimension(ilg)    :: altotacc_gat !(ilg)   !<
    real, dimension(ilg)    :: netrad_gat !(ilg)     !<
    real, dimension(ilg)    :: preacc_gat !(ilg)     !<
    real, dimension(ilg)    :: sdepgat !(ilg)        !<
    !     real, dimension(ilg,ignd)  :: rgmgat !(ilg,ignd)    !<
    real, dimension(ilg,ignd)  :: sandgat !(ilg,ignd)   !<
    real, dimension(ilg,ignd)  :: claygat !(ilg,ignd)   !<
    real, dimension(ilg,ignd)  :: orgmgat !(ilg,ignd)   !<
    real, dimension(ilg)    :: xdiffusgat !(ilg)
    real, dimension(ilg)    :: faregat !(ilg)

    real, dimension(ilg,icc) :: lygleafmasmax       !< last year maximum of the gleafmas
    real, dimension(ilg,icc) :: lystemmassmax       !< last year maximum of the stemmass
    real, dimension(ilg,icc) :: lyrootmassmax       !< last year maximum of the rootmass
    real, dimension(ilg,icc) :: lmaxt, smaxt, rmaxt !< temp vars to find previous year C pool

  end type veg_gat

  type (veg_gat), save, target :: vgat
  !$omp threadprivate(vgat)

  !=================================================================================
  type tracersType
    !   Simple tracer variables. Only written to if useTracer > 0.

    ! NOTE: Units may vary depending on the tracer used, see convertTracerUnits in tracer.f90

    ! Pools:
    ! allocated with nlat, nmos, ...:
    real, dimension(nlat,nmos) :: mossCMassrot      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: mossLitrMassrot   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: tracerCO2rot     !< Atmopspheric tracer CO2 concentration (units vary)

    ! allocated with nlat, nmos, icc:
    real, dimension(nlat,nmos,icc) :: gLeafMassrot      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, dimension(nlat,nmos,icc) :: bLeafMassrot      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemMassrot       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootMassrot       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    ! allocated with nlat, nmos, iccp2, ignd:
    real, dimension(nlat,nmos,iccp2,ignd) :: litrMassrot       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,iccp2,ignd) :: soilCMassrot      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    ! allocated with ilg, ...:
    real, dimension(ilg) :: mossCMassgat      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, dimension(ilg) :: mossLitrMassgat   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, dimension(ilg) :: tracerCO2gat     !< Atmopspheric tracer CO2 concentration (units vary)
    ! allocated with nlat, nmos, icc:
    real, dimension(ilg,icc) :: gLeafMassgat      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: bLeafMassgat      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stemMassgat       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rootMassgat       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    ! allocated with nlat, nmos, iccp2, ignd:
    real, dimension(ilg,iccp2,ignd) :: litrMassgat       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, dimension(ilg,iccp2,ignd) :: soilCMassgat      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

  end type tracersType

  type (tracersType), save, target :: tracer
  !$omp threadprivate(tracer)

  !=================================================================================
  !> CTEM's variables per tile
  type ctem_tile_level

    !   Tile-level variables (denoted by an ending of "_t")

    real, dimension(ilg) :: fsnowacc_t       !<
    real, dimension(ilg) :: CFLUX_GAacc_t    !<daily accu. boundary layer aerodynamic conductance \f$[m day^{-1} ]\f$
    real, dimension(ilg) :: USTARBS_GAacc_t  !<daily accu. friction velocity to be used in nitrogen volatilization  \f$[m day^{-1} ]\f$
    real, dimension(ilg) :: ROFBacc_t        !<daily accu. base flow from bottom of soil column \f$[kg m^{-2} day^{-1} ]\f$
    real, dimension(ilg) :: tcansacc_t       !<
    real, dimension(ilg) :: tcanoaccgat_t    !<
    real, dimension(ilg) :: taaccgat_t       !<
    real, dimension(ilg) :: uvaccgat_t       !<
    real, dimension(ilg) :: vvaccgat_t       !<
    real, dimension(ilg) :: anmossac_t       !< daily averaged moss net photosynthesis accumulated (\f$\mu mol /m^2 /s\f$)
    real, dimension(ilg) :: rmlmossac_t      !< daily averaged moss maintainence respiration (\f$\mu mol /m^2 /s\f$)
    real, dimension(ilg) :: gppmossac_t      !< daily averaged gross primary production (\f$\mu mol /m^2 /s\f$)

    ! allocated with ilg, ignd:
    real, dimension(ilg,ignd) :: QFCacc_t    !< daily accu. water removed from soil layers by transpiration \f$[kg m^{-2} day^{-1}]\f$
    real, dimension(ilg,ignd) :: tbaraccgat_t !<
    real, dimension(ilg,ignd) :: thliqacc_t  !<
    real, dimension(ilg,ignd) :: thiceacc_t  !< Added in place of YW's thicaccgat_m. EC Dec 23 2016.

    ! allocated with ilg, icc:
    real, dimension(ilg,icc) :: ancgvgac_t  !<
    real, dimension(ilg,icc) :: rmlcgvga_t  !<

  end type ctem_tile_level

  type (ctem_tile_level), save, target :: ctem_tile
  !$omp threadprivate(ctem_tile)

  !=================================================================================
  !> CTEM's variables monthly averaged (per pft)
  type ctem_monthly

    !     Tile-level monthly variables (denoted by name ending in "_mo_t")

    ! allocated with nlat, nmos, icc/iccp1/iccp2:

    real, dimension(nlat,nmos,icc) :: laimaxg_mo    !<
    real, dimension(nlat,nmos,icc) :: gleafmas_mo   !<
    real, dimension(nlat,nmos,icc) :: gleafmas_NS_mo !<
    real, dimension(nlat,nmos,icc) :: gleafmass_mo  !<
    real, dimension(nlat,nmos,icc) :: stemmass_mo   !<
    real, dimension(nlat,nmos,icc) :: stemmass_NS_mo !<
    real, dimension(nlat,nmos,icc) :: stemmasss_mo  !<
    real, dimension(nlat,nmos,icc) :: rootmass_mo   !<
    real, dimension(nlat,nmos,icc) :: rootmass_NS_mo !<
    real, dimension(nlat,nmos,icc) :: rootmasss_mo  !<    
    real, dimension(nlat,nmos,icc) :: rootdpth_mo   !<
    real, dimension(nlat,nmos,icc) :: litrfallveg_mo !<
    real, dimension(nlat,nmos,iccp2) :: humiftrsveg_mo !<
    real, dimension(nlat,nmos,icc) :: tltrleaf_mo !<
    real, dimension(nlat,nmos,icc) :: tltrstem_mo !<
    real, dimension(nlat,nmos,icc) :: tltrroot_mo !<
    real, dimension(nlat,nmos,icc) :: cc_mo        !<
    real, dimension(nlat,nmos,icc) :: mm_mo        !<
    real, dimension(nlat,nmos,icc) :: npp_mo        !<
    real, dimension(nlat,nmos,icc) :: gpp_mo        !<
    real, dimension(nlat,nmos,icc) :: vcmax0_mo     !<
    real, dimension(nlat,nmos,icc) :: vgbiomas_mo   !<
    real, dimension(nlat,nmos,icc) :: autores_mo    !<
    real, dimension(nlat,nmos,icc) :: soilres_mo    !<
    real, dimension(nlat,nmos,iccp1) :: totcmass_mo   !<
    real, dimension(nlat,nmos,iccp1) :: nep_mo        !<
    real, dimension(nlat,nmos,iccp1) :: nepCMIP_mo        !<
    real, dimension(nlat,nmos,iccp1) :: hetrores_mo   !<
    real, dimension(nlat,nmos,iccp1) :: nbp_mo        !<
    real, dimension(nlat,nmos,icc) :: emit_co2_mo  !<
    real, dimension(nlat,nmos,icc) :: emit_co_mo   !<
    real, dimension(nlat,nmos,icc) :: emit_ch4_mo  !<
    real, dimension(nlat,nmos,icc) :: emit_nmhc_mo !<
    real, dimension(nlat,nmos,icc) :: emit_h2_mo   !<
    real, dimension(nlat,nmos,icc) :: emit_nox_mo  !<
    real, dimension(nlat,nmos,icc) :: emit_n2o_mo  !<
    real, dimension(nlat,nmos,icc) :: emit_nh3_mo  !<
    real, dimension(nlat,nmos,icc) :: emit_pm25_mo !<
    real, dimension(nlat,nmos,icc) :: emit_tpm_mo  !<
    real, dimension(nlat,nmos,icc) :: emit_tc_mo   !<
    real, dimension(nlat,nmos,icc) :: emit_oc_mo   !<
    real, dimension(nlat,nmos,icc) :: emit_bc_mo   !<
    real, dimension(nlat,nmos,icc) :: burnfrac_mo  !<
    real, dimension(nlat,nmos) :: timharvarea_mo   !<
    real, dimension(nlat,nmos) :: tileAge_mo   !<
    real, dimension(nlat,nmos,icc) :: bterm_mo     !<
    real, dimension(nlat,nmos,icc) :: mterm_mo     !<
    real, dimension(nlat,nmos,icc) :: smfuncveg_mo !<
    ! allocated with nlat, nmos, icc/iccp1/iccp2, ignd:

    real, dimension(nlat,nmos,iccp1) :: nh4_mass_mo    !<
    real, dimension(nlat,nmos,iccp1) :: no3_mass_mo    !<
    real, dimension(nlat,nmos,icc) :: ngleafmas_mo   !<
    real, dimension(nlat,nmos,icc) :: ngleafmas_NS_mo !<
    real, dimension(nlat,nmos,icc) :: ngleafmass_mo  !<
    real, dimension(nlat,nmos,icc) :: nbleafmas_mo   !<
    real, dimension(nlat,nmos,icc) :: nstemmass_mo   !<
    real, dimension(nlat,nmos,icc) :: nstemmass_NS_mo !<
    real, dimension(nlat,nmos,icc) :: nstemmasss_mo  !<
    real, dimension(nlat,nmos,icc) :: nrootmass_mo   !<
    real, dimension(nlat,nmos,icc) :: nrootmass_NS_mo !<
    real, dimension(nlat,nmos,icc) :: nrootmasss_mo  !<
    real, dimension(nlat,nmos,iccp2) :: nlitrmass_mo   !<
    real, dimension(nlat,nmos,iccp2) :: soilnmas_mo    !<
    real, dimension(nlat,nmos,icc) :: nvgbiomas_mo   !<
    real, dimension(nlat,nmos,iccp2,ignd) :: litrmass_mo   !<
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcmas_mo   !<
    real, dimension(nlat,nmos,iccp2,ignd) :: litres_mo     !<
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcres_mo   !<
    real, dimension(nlat,nmos,icc) :: leafns2s_mo      !<
    real, dimension(nlat,nmos,icc) :: stemns2s_mo      !<
    real, dimension(nlat,nmos,icc) :: rootns2s_mo      !<
    real, dimension(nlat,nmos,icc) :: re_alloc_s2l_mo  !<
    real, dimension(nlat,nmos,icc) :: re_alloc_r2l_mo  !<
    real, dimension(nlat,nmos,icc) :: re_alloc_sr2l_mo !<
    real, dimension(nlat,nmos,iccp1) :: bnf_tot_mo     !<
    real, dimension(nlat,nmos,iccp1) :: bnf_free_mo     !<
    real, dimension(nlat,nmos,icc) :: bnf_ant_mo     !<
    real, dimension(nlat,nmos,icc) :: bnf_nat_mo     !<
    real, dimension(nlat,nmos,icc) :: nstress_mo
    real, dimension(nlat,nmos,iccp1) :: nitrif_mo        !<
    real, dimension(nlat,nmos,iccp1) :: no_nit_mo        !<
    real, dimension(nlat,nmos,iccp1) :: no_denit_mo      !<
    real, dimension(nlat,nmos,iccp1) :: no_nitdenit_mo   !<
    real, dimension(nlat,nmos,iccp1) :: n2o_nit_mo       !<
    real, dimension(nlat,nmos,iccp1) :: n2o_denit_mo     !<
    real, dimension(nlat,nmos,iccp1) :: n2o_nitdenit_mo  !<
    real, dimension(nlat,nmos,iccp1) :: n2_denit_mo      !<
    real, dimension(nlat,nmos,iccp1) :: nvol_mo          !<
    real, dimension(nlat,nmos,iccp1) :: nleach_mo        !<
    real, dimension(nlat,nmos,iccp1) :: appl_fert_mo     !<
    real, dimension(nlat,nmos,iccp1) :: ndep_nh4_mo      !<
    real, dimension(nlat,nmos,iccp1) :: ndep_no3_mo      !<
    real, dimension(nlat,nmos,icc) :: ndemand_wp_npp_mo!<
    real, dimension(nlat,nmos,icc) :: nuptake_p_nh4_mo !<
    real, dimension(nlat,nmos,icc) :: nuptake_p_no3_mo !<
    real, dimension(nlat,nmos,icc) :: nuptake_a_actl_nh4_mo !<
    real, dimension(nlat,nmos,icc) :: nuptake_a_actl_no3_mo !<
    real, dimension(nlat,nmos,icc) :: nuptake_mo !<
    real, dimension(nlat,nmos,icc) :: nleafns2s_mo     !<
    real, dimension(nlat,nmos,icc) :: nstemns2s_mo     !<
    real, dimension(nlat,nmos,icc) :: nrootns2s_mo     !<
    real, dimension(nlat,nmos,icc) :: nalloc_l_mo      !<
    real, dimension(nlat,nmos,icc) :: nalloc_s_mo      !<
    real, dimension(nlat,nmos,icc) :: nalloc_r_mo      !<
    real, dimension(nlat,nmos,icc) :: nresorped_s_mo   !<
    real, dimension(nlat,nmos,icc) :: nresorped_r_mo   !<
    real, dimension(nlat,nmos,icc) :: nre_alloc_s2l_mo !<
    real, dimension(nlat,nmos,icc) :: nre_alloc_r2l_mo !<
    real, dimension(nlat,nmos,icc) :: nlitr_l_mo       !<
    real, dimension(nlat,nmos,icc) :: nlitr_s_mo       !<
    real, dimension(nlat,nmos,icc) :: nlitr_r_mo       !<
    real, dimension(nlat,nmos,icc) :: nlitr_mo       !<
    real, dimension(nlat,nmos,icc) :: gl2bl_grass_nflux_mo !<
    real, dimension(nlat,nmos,icc) :: c2n_l_mo         !<
    real, dimension(nlat,nmos,icc) :: c2n_s_mo         !<
    real, dimension(nlat,nmos,icc) :: c2n_r_mo         !<
    real, dimension(nlat,nmos,icc) :: c2n_wp_mo        !<
    real, dimension(nlat,nmos,iccp1) :: c2n_litr_mo      !<
    real, dimension(nlat,nmos,iccp1) :: c2n_humus_mo     !<
    real, dimension(nlat,nmos,iccp1) :: nhumtrs_mo       !<
    real, dimension(nlat,nmos,iccp1) :: nmineral_litr_mo !<
    real, dimension(nlat,nmos,iccp1) :: nmineral_humus_mo!<
    real, dimension(nlat,nmos,iccp1) :: nimmobil_nh4_mo  !<
    real, dimension(nlat,nmos,iccp1) :: nimmobil_no3_mo  !<
    real, dimension(nlat,nmos,iccp1) :: netnmineral_mo   !<
    real, dimension(nlat,nmos,iccp1) :: fNnetland_mo     !<

  end type ctem_monthly

  type (ctem_monthly), save, target :: ctem_mo

  !=================================================================================
  !> CTEM's grid average monthly values
  type ctem_gridavg_monthly

    !  Grid averaged monthly variables (denoted by name ending in "_mo_g")

    ! allocated with nlat:
    real, dimension(nlat) :: laimaxg_mo_g  !<
    real, dimension(nlat) :: gleafmas_mo_g !<
    real, dimension(nlat) :: gleafmas_NS_mo_g !<
    real, dimension(nlat) :: gleafmass_mo_g !<
    real, dimension(nlat) :: stemmass_mo_g !<
    real, dimension(nlat) :: stemmass_NS_mo_g !<
    real, dimension(nlat) :: stemmasss_mo_g  !<
    real, dimension(nlat) :: rootmass_mo_g !<
    real, dimension(nlat) :: rootmass_NS_mo_g  !<
    real, dimension(nlat) :: rootmasss_mo_g   !<
    real, dimension(nlat) :: nh4_mass_mo_g    !<
    real, dimension(nlat) :: no3_mass_mo_g    !<
    real, dimension(nlat) :: ngleafmas_mo_g   !<
    real, dimension(nlat) :: ngleafmas_NS_mo_g !<
    real, dimension(nlat) :: ngleafmass_mo_g  !<
    real, dimension(nlat) :: nbleafmas_mo_g   !<
    real, dimension(nlat) :: nstemmass_mo_g   !<
    real, dimension(nlat) :: nstemmass_NS_mo_g !<
    real, dimension(nlat) :: nstemmasss_mo_g  !<
    real, dimension(nlat) :: nrootmass_mo_g   !<
    real, dimension(nlat) :: nrootmass_NS_mo_g !<
    real, dimension(nlat) :: nrootmasss_mo_g  !<
    real, dimension(nlat) :: nlitrmass_mo_g   !<
    real, dimension(nlat) :: soilnmas_mo_g    !<
    real, dimension(nlat) :: nvgbiomas_mo_g   !<

    real, dimension(nlat) :: litrfall_mo_g !<
    real, dimension(nlat) :: humiftrs_mo_g !<
    real, dimension(nlat) :: tltrleaf_mo_g !<
    real, dimension(nlat) :: tltrstem_mo_g !<
    real, dimension(nlat) :: tltrroot_mo_g !<
    real, dimension(nlat) :: npp_mo_g      !<
    real, dimension(nlat) :: gpp_mo_g      !<
    real, dimension(nlat) :: vcmax0_mo_g   !<
    real, dimension(nlat) :: nep_mo_g      !<
    real, dimension(nlat) :: nepCMIP_mo_g      !<
    real, dimension(nlat) :: nbp_mo_g      !<
    real, dimension(nlat) :: hetrores_mo_g !<
    real, dimension(nlat) :: autores_mo_g  !<
    real, dimension(nlat) :: soilres_mo_g  !<
    real, dimension(nlat) :: vgbiomas_mo_g !<
    real, dimension(nlat) :: totcmass_mo_g !<
    real, dimension(nlat) :: emit_co2_mo_g !<
    real, dimension(nlat) :: emit_co_mo_g  !<
    real, dimension(nlat) :: emit_ch4_mo_g !<
    real, dimension(nlat) :: emit_nmhc_mo_g !<
    real, dimension(nlat) :: emit_h2_mo_g  !<
    real, dimension(nlat) :: emit_nox_mo_g !<
    real, dimension(nlat) :: emit_n2o_mo_g !<
    real, dimension(nlat) :: emit_nh3_mo_g !<    
    real, dimension(nlat) :: emit_pm25_mo_g !<
    real, dimension(nlat) :: emit_tpm_mo_g !<
    real, dimension(nlat) :: emit_tc_mo_g  !<
    real, dimension(nlat) :: emit_oc_mo_g  !<
    real, dimension(nlat) :: emit_bc_mo_g  !<
    real, dimension(nlat) :: smfuncveg_mo_g !<
    real, dimension(nlat) :: luc_emc_mo_g  !<
    real, dimension(nlat) :: lucltrin_mo_g !<
    real, dimension(nlat) :: lucsocin_mo_g !<
    real, dimension(nlat) :: luc_emcn_mo_g  !<
    real, dimension(nlat) :: lucltrinn_mo_g !<
    real, dimension(nlat) :: lucsocinn_mo_g !<
    real, dimension(nlat) :: burnfrac_mo_g !<
    real, dimension(nlat) :: timharvarea_mo_g !<
    real, dimension(nlat) :: tileAge_mo_g !<
    real, dimension(nlat) :: bterm_mo_g    !<
    real, dimension(nlat) :: lterm_mo_g    !<
    real, dimension(nlat) :: mterm_mo_g    !<
    real, dimension(nlat) :: ch4WetSpec_mo_g  !<
    real, dimension(nlat) :: wetfdyn_mo_g  !<
    real, dimension(nlat) :: wetfpres_mo_g  !<
    real, dimension(nlat) :: ch4WetDyn_mo_g  !<
    real, dimension(nlat) :: ch4soills_mo_g !<
    real, dimension(nlat) :: cProduct_mo_g          !< Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    real, dimension(nlat) :: nProduct_mo_g          !< Nitrogen in the LUC product pools (litter and soil N iccp2 position) \f$[g N m^{-2}]\f$
    real, dimension(nlat) :: fProductDecomp_mo_g    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, ignd:
    real, dimension(nlat,ignd) :: litrmass_mo_g !<
    real, dimension(nlat,ignd) :: soilcmas_mo_g !<
    real, dimension(nlat,ignd) :: litres_mo_g   !<
    real, dimension(nlat,ignd) :: soilcres_mo_g !<
    real, dimension(nlat) :: leafns2s_mo_g      !<
    real, dimension(nlat) :: stemns2s_mo_g      !<
    real, dimension(nlat) :: rootns2s_mo_g      !<
    real, dimension(nlat) :: re_alloc_s2l_mo_g  !<
    real, dimension(nlat) :: re_alloc_r2l_mo_g  !<
    real, dimension(nlat) :: re_alloc_sr2l_mo_g !<
    real, dimension(nlat) :: bnf_tot_mo_g     !<
    real, dimension(nlat) :: bnf_free_mo_g     !<
    real, dimension(nlat) :: bnf_ant_mo_g     !<
    real, dimension(nlat) :: bnf_nat_mo_g     !<
    real, dimension(nlat) :: nstress_mo_g
    real, dimension(nlat) :: nitrif_mo_g        !<
    real, dimension(nlat) :: no_nit_mo_g        !<
    real, dimension(nlat) :: no_denit_mo_g      !<
    real, dimension(nlat) :: no_nitdenit_mo_g   !<
    real, dimension(nlat) :: n2o_nit_mo_g       !<
    real, dimension(nlat) :: n2o_denit_mo_g     !<
    real, dimension(nlat) :: n2o_nitdenit_mo_g  !<
    real, dimension(nlat) :: n2_denit_mo_g      !<
    real, dimension(nlat) :: nvol_mo_g          !<
    real, dimension(nlat) :: nleach_mo_g        !<
    real, dimension(nlat) :: appl_fert_mo_g     !<
    real, dimension(nlat) :: ndep_nh4_mo_g      !<
    real, dimension(nlat) :: ndep_no3_mo_g      !<
    real, dimension(nlat) :: ndemand_wp_npp_mo_g!<
    real, dimension(nlat) :: nuptake_p_nh4_mo_g !<
    real, dimension(nlat) :: nuptake_p_no3_mo_g !<
    real, dimension(nlat) :: nuptake_a_actl_nh4_mo_g !<
    real, dimension(nlat) :: nuptake_a_actl_no3_mo_g !<
    real, dimension(nlat) :: nuptake_mo_g !<
    real, dimension(nlat) :: nleafns2s_mo_g     !<
    real, dimension(nlat) :: nstemns2s_mo_g     !<
    real, dimension(nlat) :: nrootns2s_mo_g     !<
    real, dimension(nlat) :: nalloc_l_mo_g      !<
    real, dimension(nlat) :: nalloc_s_mo_g      !<
    real, dimension(nlat) :: nalloc_r_mo_g      !<
    real, dimension(nlat) :: nresorped_s_mo_g   !<
    real, dimension(nlat) :: nresorped_r_mo_g   !<
    real, dimension(nlat) :: nre_alloc_s2l_mo_g !<
    real, dimension(nlat) :: nre_alloc_r2l_mo_g !<
    real, dimension(nlat) :: nlitr_l_mo_g       !<
    real, dimension(nlat) :: nlitr_s_mo_g       !<
    real, dimension(nlat) :: nlitr_r_mo_g       !<
    real, dimension(nlat) :: nlitr_mo_g       !<
    real, dimension(nlat) :: gl2bl_grass_nflux_mo_g !<
    real, dimension(nlat) :: c2n_l_mo_g         !<
    real, dimension(nlat) :: c2n_s_mo_g         !<
    real, dimension(nlat) :: c2n_r_mo_g         !<
    real, dimension(nlat) :: c2n_wp_mo_g        !<
    real, dimension(nlat) :: c2n_litr_mo_g      !<
    real, dimension(nlat) :: c2n_humus_mo_g     !<
    real, dimension(nlat) :: nhumtrs_mo_g       !<
    real, dimension(nlat) :: nmineral_litr_mo_g !<
    real, dimension(nlat) :: nmineral_humus_mo_g!<
    real, dimension(nlat) :: nimmobil_nh4_mo_g  !<
    real, dimension(nlat) :: nimmobil_no3_mo_g  !<
    real, dimension(nlat) :: netnmineral_mo_g   !<
    real, dimension(nlat) :: fNnetland_mo_g     !<

  end type ctem_gridavg_monthly

  type (ctem_gridavg_monthly), save, target :: ctem_grd_mo

  !=================================================================================
  !> CTEM's variables per tile monthly values
  type ctem_tileavg_monthly

    !     Tile-level monthly variables (denoted by name ending in "_mo_t")

    ! allocated with nlat, nmos:

    real, dimension(nlat,nmos) :: laimaxg_mo_t  !<
    real, dimension(nlat,nmos) :: gleafmas_mo_t    !<
    real, dimension(nlat,nmos) :: gleafmas_NS_mo_t  !<
    real, dimension(nlat,nmos) :: gleafmass_mo_t   !<
    real, dimension(nlat,nmos) :: stemmass_mo_t !<
    real, dimension(nlat,nmos) :: stemmass_NS_mo_t  !<
    real, dimension(nlat,nmos) :: stemmasss_mo_t   !<
    real, dimension(nlat,nmos) :: rootmass_mo_t !<
    real, dimension(nlat,nmos) :: rootmass_NS_mo_t !<
    real, dimension(nlat,nmos) :: rootmasss_mo_t   !<
    real, dimension(nlat,nmos) :: litrfall_mo_t !<
    real, dimension(nlat,nmos) :: humiftrs_mo_t !<
    real, dimension(nlat,nmos) :: tltrleaf_mo_t !<
    real, dimension(nlat,nmos) :: tltrstem_mo_t !<
    real, dimension(nlat,nmos) :: tltrroot_mo_t !<

    real, dimension(nlat,nmos) :: npp_mo_t      !<
    real, dimension(nlat,nmos) :: gpp_mo_t      !<
    real, dimension(nlat,nmos) :: vcmax0_mo_t      !<
    real, dimension(nlat,nmos) :: vgbiomas_mo_t !<
    real, dimension(nlat,nmos) :: autores_mo_t  !<
    real, dimension(nlat,nmos) :: soilres_mo_t  !<
    real, dimension(nlat,nmos) :: totcmass_mo_t !<
    real, dimension(nlat,nmos) :: nh4_mass_mo_t    !<
    real, dimension(nlat,nmos) :: no3_mass_mo_t    !<
    real, dimension(nlat,nmos) :: ngleafmas_mo_t   !<
    real, dimension(nlat,nmos) :: ngleafmas_NS_mo_t !<
    real, dimension(nlat,nmos) :: ngleafmass_mo_t  !<
    real, dimension(nlat,nmos) :: nbleafmas_mo_t   !<
    real, dimension(nlat,nmos) :: nstemmass_mo_t   !<
    real, dimension(nlat,nmos) :: nstemmass_NS_mo_t !<
    real, dimension(nlat,nmos) :: nstemmasss_mo_t  !<
    real, dimension(nlat,nmos) :: nrootmass_mo_t   !<
    real, dimension(nlat,nmos) :: nrootmass_NS_mo_t !<
    real, dimension(nlat,nmos) :: nrootmasss_mo_t  !<
    real, dimension(nlat,nmos) :: nlitrmass_mo_t   !<
    real, dimension(nlat,nmos) :: soilnmas_mo_t    !<
    real, dimension(nlat,nmos) :: nvgbiomas_mo_t   !<
    real, dimension(nlat,nmos) :: nep_mo_t      !<
    real, dimension(nlat,nmos) :: nepCMIP_mo_t      !<
    real, dimension(nlat,nmos) :: hetrores_mo_t !<
    real, dimension(nlat,nmos) :: nbp_mo_t      !<
    real, dimension(nlat,nmos) :: emit_co2_mo_t !<
    real, dimension(nlat,nmos) :: emit_co_mo_t  !<
    real, dimension(nlat,nmos) :: emit_ch4_mo_t !<
    real, dimension(nlat,nmos) :: emit_nmhc_mo_t !<
    real, dimension(nlat,nmos) :: emit_h2_mo_t  !<
    real, dimension(nlat,nmos) :: emit_nox_mo_t !<
    real, dimension(nlat,nmos) :: emit_n2o_mo_t !<
    real, dimension(nlat,nmos) :: emit_nh3_mo_t !<
    real, dimension(nlat,nmos) :: emit_pm25_mo_t !<
    real, dimension(nlat,nmos) :: emit_tpm_mo_t !<
    real, dimension(nlat,nmos) :: emit_tc_mo_t  !<
    real, dimension(nlat,nmos) :: emit_oc_mo_t  !<
    real, dimension(nlat,nmos) :: emit_bc_mo_t  !<
    real, dimension(nlat,nmos) :: burnfrac_mo_t !<
    real, dimension(nlat,nmos) :: timharvarea_mo_t !<
    real, dimension(nlat,nmos) :: tileAge_mo_t !<
    real, dimension(nlat,nmos) :: smfuncveg_mo_t !<
    real, dimension(nlat,nmos) :: bterm_mo_t    !<
    real, dimension(nlat,nmos) :: luc_emc_mo_t  !<
    real, dimension(nlat,nmos) :: lterm_mo_t    !<
    real, dimension(nlat,nmos) :: lucsocin_mo_t !<
    real, dimension(nlat,nmos) :: mterm_mo_t    !<
    real, dimension(nlat,nmos) :: lucltrin_mo_t !<
    real, dimension(nlat,nmos) :: luc_emcn_mo_t  !<
    real, dimension(nlat,nmos) :: lucsocinn_mo_t !<
    real, dimension(nlat,nmos) :: lucltrinn_mo_t !<
    real, dimension(nlat,nmos) :: ch4WetSpec_mo_t  !<
    real, dimension(nlat,nmos) :: wetfdyn_mo_t  !<
    real, dimension(nlat,nmos) :: wetfpres_mo_t  !<
    real, dimension(nlat,nmos) :: ch4WetDyn_mo_t  !<
    real, dimension(nlat,nmos) :: ch4soills_mo_t !<
    real, dimension(nlat,nmos) :: wind_mo_t     !<
    real, dimension(nlat,nmos) :: fProductDecomp_mo_t    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, nmos, ignd:
    real, dimension(nlat,nmos,ignd) :: litrmass_mo_t !<
    real, dimension(nlat,nmos,ignd) :: soilcmas_mo_t !<
    real, dimension(nlat,nmos,ignd) :: litres_mo_t   !<
    real, dimension(nlat,nmos,ignd) :: soilcres_mo_t !<
    real, dimension(nlat,nmos) :: leafns2s_mo_t      !<
    real, dimension(nlat,nmos) :: stemns2s_mo_t      !<
    real, dimension(nlat,nmos) :: rootns2s_mo_t      !<
    real, dimension(nlat,nmos) :: re_alloc_s2l_mo_t  !<
    real, dimension(nlat,nmos) :: re_alloc_r2l_mo_t  !<
    real, dimension(nlat,nmos) :: re_alloc_sr2l_mo_t !<
    real, dimension(nlat,nmos) :: bnf_tot_mo_t     !<
    real, dimension(nlat,nmos) :: bnf_free_mo_t     !<
    real, dimension(nlat,nmos) :: bnf_ant_mo_t     !<
    real, dimension(nlat,nmos) :: bnf_nat_mo_t     !<
    real, dimension(nlat,nmos) :: nstress_mo_t
    real, dimension(nlat,nmos) :: nitrif_mo_t        !<
    real, dimension(nlat,nmos) :: no_nit_mo_t        !<
    real, dimension(nlat,nmos) :: no_denit_mo_t      !<
    real, dimension(nlat,nmos) :: no_nitdenit_mo_t   !<
    real, dimension(nlat,nmos) :: n2o_nit_mo_t       !<
    real, dimension(nlat,nmos) :: n2o_denit_mo_t     !<
    real, dimension(nlat,nmos) :: n2o_nitdenit_mo_t  !<
    real, dimension(nlat,nmos) :: n2_denit_mo_t      !<
    real, dimension(nlat,nmos) :: nvol_mo_t          !<
    real, dimension(nlat,nmos) :: nleach_mo_t        !<
    real, dimension(nlat,nmos) :: appl_fert_mo_t     !<
    real, dimension(nlat,nmos) :: ndep_nh4_mo_t      !<
    real, dimension(nlat,nmos) :: ndep_no3_mo_t      !<
    real, dimension(nlat,nmos) :: ndemand_wp_npp_mo_t!<
    real, dimension(nlat,nmos) :: nuptake_p_nh4_mo_t !<
    real, dimension(nlat,nmos) :: nuptake_p_no3_mo_t !<
    real, dimension(nlat,nmos) :: nuptake_a_actl_nh4_mo_t !<
    real, dimension(nlat,nmos) :: nuptake_a_actl_no3_mo_t !<
    real, dimension(nlat,nmos) :: nuptake_mo_t !<
    real, dimension(nlat,nmos) :: nleafns2s_mo_t     !<
    real, dimension(nlat,nmos) :: nstemns2s_mo_t     !<
    real, dimension(nlat,nmos) :: nrootns2s_mo_t     !<
    real, dimension(nlat,nmos) :: nalloc_l_mo_t      !<
    real, dimension(nlat,nmos) :: nalloc_s_mo_t      !<
    real, dimension(nlat,nmos) :: nalloc_r_mo_t      !<
    real, dimension(nlat,nmos) :: nresorped_s_mo_t   !<
    real, dimension(nlat,nmos) :: nresorped_r_mo_t   !<
    real, dimension(nlat,nmos) :: nre_alloc_s2l_mo_t !<
    real, dimension(nlat,nmos) :: nre_alloc_r2l_mo_t !<
    real, dimension(nlat,nmos) :: nlitr_l_mo_t       !<
    real, dimension(nlat,nmos) :: nlitr_s_mo_t       !<
    real, dimension(nlat,nmos) :: nlitr_r_mo_t       !<
    real, dimension(nlat,nmos) :: nlitr_mo_t       !<
    real, dimension(nlat,nmos) :: gl2bl_grass_nflux_mo_t !<
    real, dimension(nlat,nmos) :: c2n_l_mo_t         !<
    real, dimension(nlat,nmos) :: c2n_s_mo_t         !<
    real, dimension(nlat,nmos) :: c2n_r_mo_t         !<
    real, dimension(nlat,nmos) :: c2n_wp_mo_t        !<
    real, dimension(nlat,nmos) :: c2n_litr_mo_t      !<
    real, dimension(nlat,nmos) :: c2n_humus_mo_t     !<
    real, dimension(nlat,nmos) :: nhumtrs_mo_t       !<
    real, dimension(nlat,nmos) :: nmineral_litr_mo_t !<
    real, dimension(nlat,nmos) :: nmineral_humus_mo_t!<
    real, dimension(nlat,nmos) :: nimmobil_nh4_mo_t  !<
    real, dimension(nlat,nmos) :: nimmobil_no3_mo_t  !<
    real, dimension(nlat,nmos) :: netnmineral_mo_t!<
    real, dimension(nlat,nmos) :: fNnetland_mo_t     !<

  end type ctem_tileavg_monthly

  type (ctem_tileavg_monthly), save, target :: ctem_tile_mo


  !=================================================================================
  !> CTEM's average annual values (per PFT)
  type ctem_annual

    ! c      Annual output for CTEM mosaic variables:
    ! c      (denoted by name ending in "_yr_m")
    !
    ! allocated with nlat, nmos, icc/iccp1/iccp2:

    real, dimension(nlat,nmos,icc) :: laimaxg_yr   !<
    real, dimension(nlat,nmos,icc) :: gleafmas_yr    !<
    real, dimension(nlat,nmos,icc) :: gleafmas_NS_yr  !<
    real, dimension(nlat,nmos,icc) :: gleafmass_yr   !<
    real, dimension(nlat,nmos,icc) :: bleafmas_yr    !<
    real, dimension(nlat,nmos,icc) :: stemmass_yr  !<
    real, dimension(nlat,nmos,icc) :: stemmass_NS_yr  !<
    real, dimension(nlat,nmos,icc) :: stemmasss_yr   !<
    real, dimension(nlat,nmos,icc) :: rootmass_yr  !<
    real, dimension(nlat,nmos,icc) :: rootmass_NS_yr  !<
    real, dimension(nlat,nmos,icc) :: rootmasss_yr   !<
    real, dimension(nlat,nmos,icc) :: cc_yr       !<
    real, dimension(nlat,nmos,icc) :: mm_yr       !<
    real, dimension(nlat,nmos,icc) :: npp_yr       !<
    real, dimension(nlat,nmos,icc) :: gpp_yr       !<
    real, dimension(nlat,nmos,icc) :: vcmax0_yr      !<
    real, dimension(nlat,nmos,icc) :: leafns2s_yr      !<
    real, dimension(nlat,nmos,icc) :: stemns2s_yr      !<
    real, dimension(nlat,nmos,icc) :: rootns2s_yr      !<
    real, dimension(nlat,nmos,icc) :: re_alloc_s2l_yr  !<
    real, dimension(nlat,nmos,icc) :: re_alloc_r2l_yr  !<
    real, dimension(nlat,nmos,icc) :: re_alloc_sr2l_yr !<
    real, dimension(nlat,nmos,icc) :: vgbiomas_yr  !<
    real, dimension(nlat,nmos,icc) :: autores_yr   !<
    real, dimension(nlat,nmos,icc) :: rmrveg_yr   !<
    real, dimension(nlat,nmos,iccp1) :: totcmass_yr !<
    real, dimension(nlat,nmos,iccp1) :: nh4_mass_yr     !<
    real, dimension(nlat,nmos,iccp1) :: no3_mass_yr     !<
    real, dimension(nlat,nmos,icc) :: ngleafmas_yr    !<
    real, dimension(nlat,nmos,icc) :: ngleafmas_NS_yr  !<
    real, dimension(nlat,nmos,icc) :: ngleafmass_yr   !<
    real, dimension(nlat,nmos,icc) :: nbleafmas_yr    !<
    real, dimension(nlat,nmos,icc) :: nstemmass_yr    !<
    real, dimension(nlat,nmos,icc) :: nstemmass_NS_yr  !<
    real, dimension(nlat,nmos,icc) :: nstemmasss_yr   !<
    real, dimension(nlat,nmos,icc) :: nrootmass_yr    !<
    real, dimension(nlat,nmos,icc) :: nrootmass_NS_yr  !<
    real, dimension(nlat,nmos,icc) :: nrootmasss_yr   !<
    real, dimension(nlat,nmos,iccp2) :: nlitrmass_yr    !<
    real, dimension(nlat,nmos,iccp2) :: soilnmas_yr     !<
    real, dimension(nlat,nmos,icc) :: nvgbiomas_yr    !<
    real, dimension(nlat,nmos,iccp1) :: nep_yr     !<
    real, dimension(nlat,nmos,iccp1) :: nepCMIP_yr     !<
    real, dimension(nlat,nmos,icc) :: litrfall_yr    !<
    real, dimension(nlat,nmos,iccp1) :: hetrores_yr !<
    real, dimension(nlat,nmos,iccp1) :: nbp_yr     !<
    real, dimension(nlat,nmos,icc) :: emit_co2_yr  !<
    real, dimension(nlat,nmos,icc) :: emit_co_yr   !<
    real, dimension(nlat,nmos,icc) :: emit_ch4_yr  !<
    real, dimension(nlat,nmos,icc) :: emit_nmhc_yr !<
    real, dimension(nlat,nmos,icc) :: emit_h2_yr   !<
    real, dimension(nlat,nmos,icc) :: emit_nox_yr  !<
    real, dimension(nlat,nmos,icc) :: emit_n2o_yr  !<
    real, dimension(nlat,nmos,icc) :: emit_nh3_yr  !<
    real, dimension(nlat,nmos,icc) :: emit_pm25_yr !<
    real, dimension(nlat,nmos,icc) :: emit_tpm_yr  !<
    real, dimension(nlat,nmos,icc) :: emit_tc_yr   !<
    real, dimension(nlat,nmos,icc) :: emit_oc_yr   !<
    real, dimension(nlat,nmos,icc) :: emit_bc_yr   !<
    real, dimension(nlat,nmos,icc) :: bterm_yr     !<
    real, dimension(nlat,nmos,icc) :: mterm_yr     !<
    real, dimension(nlat,nmos,icc) :: burnfrac_yr  !<
    real, dimension(nlat,nmos) :: timharvarea_yr   !<
    real, dimension(nlat,nmos) :: tileAge_yr !<
    real, dimension(nlat,nmos,icc) :: smfuncveg_yr !<
    real, dimension(nlat,nmos,icc) :: veghght_yr   !<

    ! allocated with nlat, nmos, iccp2, ignd:
    real, dimension(nlat,nmos,iccp2,ignd) :: litrmass_yr !<
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcmas_yr !<
    real, dimension(nlat,nmos,iccp2,ignd) :: litres_yr  !<
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcres_yr !<
    real, dimension(nlat,nmos,iccp1) :: bnf_tot_yr   !<
    real, dimension(nlat,nmos,iccp1) :: bnf_free_yr   !<
    real, dimension(nlat,nmos,icc) :: bnf_ant_yr   !<
    real, dimension(nlat,nmos,icc) :: bnf_nat_yr   !<
    real, dimension(nlat,nmos,icc) :: nstress_yr
    real, dimension(nlat,nmos,iccp1) :: nitrif_yr      !<
    real, dimension(nlat,nmos,iccp1) :: no_nit_yr      !<
    real, dimension(nlat,nmos,iccp1) :: no_denit_yr    !<
    real, dimension(nlat,nmos,iccp1) :: no_nitdenit_yr !<
    real, dimension(nlat,nmos,iccp1) :: n2o_nit_yr     !<
    real, dimension(nlat,nmos,iccp1) :: n2o_denit_yr   !<
    real, dimension(nlat,nmos,iccp1) :: n2o_nitdenit_yr!<
    real, dimension(nlat,nmos,iccp1) :: n2_denit_yr    !<
    real, dimension(nlat,nmos,iccp1) :: nvol_yr        !<
    real, dimension(nlat,nmos,iccp1) :: nleach_yr      !<
    real, dimension(nlat,nmos,iccp1) :: appl_fert_yr   !<
    real, dimension(nlat,nmos,iccp1) :: ndep_nh4_yr    !<
    real, dimension(nlat,nmos,iccp1) :: ndep_no3_yr    !<
    real, dimension(nlat,nmos,icc) :: ndemand_wp_npp_yr     !<
    real, dimension(nlat,nmos,icc) :: nuptake_p_nh4_yr      !<
    real, dimension(nlat,nmos,icc) :: nuptake_p_no3_yr      !<
    real, dimension(nlat,nmos,icc) :: nuptake_a_actl_nh4_yr !<
    real, dimension(nlat,nmos,icc) :: nuptake_a_actl_no3_yr !<
    real, dimension(nlat,nmos,icc) :: nuptake_yr            !<
    real, dimension(nlat,nmos,icc) :: nleafns2s_yr     !<
    real, dimension(nlat,nmos,icc) :: nstemns2s_yr     !<
    real, dimension(nlat,nmos,icc) :: nrootns2s_yr     !<
    real, dimension(nlat,nmos,icc) :: nalloc_l_yr      !<
    real, dimension(nlat,nmos,icc) :: nalloc_s_yr      !<
    real, dimension(nlat,nmos,icc) :: nalloc_r_yr      !<
    real, dimension(nlat,nmos,icc) :: nresorped_s_yr   !<
    real, dimension(nlat,nmos,icc) :: nresorped_r_yr   !<
    real, dimension(nlat,nmos,icc) :: nre_alloc_s2l_yr !<
    real, dimension(nlat,nmos,icc) :: nre_alloc_r2l_yr !<
    real, dimension(nlat,nmos,icc) :: nlitr_l_yr       !<
    real, dimension(nlat,nmos,icc) :: nlitr_s_yr       !<
    real, dimension(nlat,nmos,icc) :: nlitr_r_yr       !<
    real, dimension(nlat,nmos,icc) :: nlitr_yr         !<
    real, dimension(nlat,nmos,icc) :: gl2bl_grass_nflux_yr !<
    real, dimension(nlat,nmos,icc) :: c2n_l_yr         !<
    real, dimension(nlat,nmos,icc) :: c2n_s_yr         !<
    real, dimension(nlat,nmos,icc) :: c2n_r_yr         !<
    real, dimension(nlat,nmos,icc) :: c2n_wp_yr        !<
    real, dimension(nlat,nmos,iccp1) :: c2n_litr_yr      !<
    real, dimension(nlat,nmos,iccp1) :: c2n_humus_yr     !<
    real, dimension(nlat,nmos,iccp1) :: nhumtrs_yr       !<
    real, dimension(nlat,nmos,iccp1) :: nmineral_litr_yr !<
    real, dimension(nlat,nmos,iccp1) :: nmineral_humus_yr!<
    real, dimension(nlat,nmos,iccp1) :: nimmobil_nh4_yr  !<
    real, dimension(nlat,nmos,iccp1) :: nimmobil_no3_yr  !<
    real, dimension(nlat,nmos,iccp1) :: netnmineral_yr   !<
    real, dimension(nlat,nmos,iccp1) :: fNnetland_yr     !<

  end type ctem_annual

  type (ctem_annual), save, target :: ctem_yr

  !=================================================================================
  !> CTEM's grid average annual values
  type ctem_gridavg_annual

    ! Annual output for CTEM grid-averaged variables:
    ! (denoted by name ending in "_yr_g")

    ! allocated with nlat:
    real, dimension(nlat) :: laimaxg_yr_g  !<
    real, dimension(nlat) :: gleafmas_yr_g    !<
    real, dimension(nlat) :: gleafmas_NS_yr_g  !<
    real, dimension(nlat) :: gleafmass_yr_g   !<
    real, dimension(nlat) :: bleafmas_yr_g    !<
    real, dimension(nlat) :: stemmass_yr_g !<
    real, dimension(nlat) :: stemmass_NS_yr_g  !<
    real, dimension(nlat) :: stemmasss_yr_g   !<
    real, dimension(nlat) :: rootmass_yr_g !<
    real, dimension(nlat) :: rootmass_NS_yr_g  !<
    real, dimension(nlat) :: rootmasss_yr_g   !<
    real, dimension(nlat) :: nh4_mass_yr_g    !<
    real, dimension(nlat) :: no3_mass_yr_g    !<
    real, dimension(nlat) :: ngleafmas_yr_g   !<
    real, dimension(nlat) :: ngleafmas_NS_yr_g !<
    real, dimension(nlat) :: ngleafmass_yr_g  !<
    real, dimension(nlat) :: nbleafmas_yr_g   !<
    real, dimension(nlat) :: nstemmass_yr_g   !<
    real, dimension(nlat) :: nstemmass_NS_yr_g !<
    real, dimension(nlat) :: nstemmasss_yr_g  !<
    real, dimension(nlat) :: nrootmass_yr_g   !<
    real, dimension(nlat) :: nrootmass_NS_yr_g !<
    real, dimension(nlat) :: nrootmasss_yr_g  !<
    real, dimension(nlat) :: nlitrmass_yr_g   !<
    real, dimension(nlat) :: soilnmas_yr_g    !<
    real, dimension(nlat) :: nvgbiomas_yr_g   !<
    real, dimension(nlat) :: npp_yr_g      !<
    real, dimension(nlat) :: gpp_yr_g      !<
    real, dimension(nlat) :: vcmax0_yr_g        !<
    real, dimension(nlat) :: leafns2s_yr_g      !<
    real, dimension(nlat) :: stemns2s_yr_g      !<
    real, dimension(nlat) :: rootns2s_yr_g      !<
    real, dimension(nlat) :: re_alloc_s2l_yr_g  !<
    real, dimension(nlat) :: re_alloc_r2l_yr_g  !<
    real, dimension(nlat) :: re_alloc_sr2l_yr_g !<
    real, dimension(nlat) :: nep_yr_g      !<
    real, dimension(nlat) :: nepCMIP_yr_g      !<
    real, dimension(nlat) :: nbp_yr_g      !<
    real, dimension(nlat) :: hetrores_yr_g !<
    real, dimension(nlat) :: autores_yr_g  !<
    real, dimension(nlat) :: rmrveg_yr_g  !<
    real, dimension(nlat) :: litrfall_yr_g    !<
    real, dimension(nlat) :: vgbiomas_yr_g !<
    real, dimension(nlat) :: totcmass_yr_g !<
    real, dimension(nlat) :: emit_co2_yr_g !<
    real, dimension(nlat) :: emit_co_yr_g  !<
    real, dimension(nlat) :: emit_ch4_yr_g !<
    real, dimension(nlat) :: emit_nmhc_yr_g !<
    real, dimension(nlat) :: emit_h2_yr_g  !<
    real, dimension(nlat) :: emit_nox_yr_g !<
    real, dimension(nlat) :: emit_n2o_yr_g !<
    real, dimension(nlat) :: emit_nh3_yr_g !<
    real, dimension(nlat) :: emit_pm25_yr_g !<
    real, dimension(nlat) :: emit_tpm_yr_g !<
    real, dimension(nlat) :: emit_tc_yr_g  !<
    real, dimension(nlat) :: emit_oc_yr_g  !<
    real, dimension(nlat) :: emit_bc_yr_g  !<
    real, dimension(nlat) :: smfuncveg_yr_g !<
    real, dimension(nlat) :: luc_emc_yr_g  !<
    real, dimension(nlat) :: lucltrin_yr_g !<
    real, dimension(nlat) :: lucsocin_yr_g !<
    real, dimension(nlat) :: luc_emcn_yr_g  !<
    real, dimension(nlat) :: lucltrinn_yr_g !<
    real, dimension(nlat) :: lucsocinn_yr_g !<
    real, dimension(nlat) :: burnfrac_yr_g !<
    real, dimension(nlat) :: timharvarea_yr_g !<
    real, dimension(nlat) :: tileAge_yr_g !<
    real, dimension(nlat) :: bterm_yr_g    !<
    real, dimension(nlat) :: lterm_yr_g    !<
    real, dimension(nlat) :: mterm_yr_g    !<
    real, dimension(nlat) :: ch4WetSpec_yr_g  !<
    real, dimension(nlat) :: wetfdyn_yr_g  !<
    real, dimension(nlat) :: ch4WetDyn_yr_g  !<
    real, dimension(nlat) :: ch4soills_yr_g !<
    real, dimension(nlat) :: veghght_yr_g  !<
    real, dimension(nlat) :: peatdep_yr_g  !<
    real, dimension(nlat) :: peatSoilC_yr_g  !<
    real, dimension(nlat) :: cProduct_yr_g          !< Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    real, dimension(nlat) :: nProduct_yr_g          !< Nitrogen in the LUC product pools (litter and soil N iccp2 position) \f$[g N m^{-2}]\f$
    real, dimension(nlat) :: fProductDecomp_yr_g    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, ignd:
    real, dimension(nlat,ignd) :: litrmass_yr_g !<
    real, dimension(nlat, ignd) :: soilcmas_yr_g !<
    real, dimension(nlat,ignd) :: litres_yr_g   !<
    real, dimension(nlat,ignd) :: soilcres_yr_g !<
    real, dimension(nlat) :: bnf_tot_yr_g  !<
    real, dimension(nlat) :: bnf_s_Cveg_yr_g  !<
    real, dimension(nlat) :: bnf_free_yr_g  !<
    real, dimension(nlat) :: bnf_ant_yr_g  !<
    real, dimension(nlat) :: bnf_nat_yr_g  !<
    real, dimension(nlat) :: nstress_yr_g
    real, dimension(nlat) :: nitrif_yr_g      !<
    real, dimension(nlat) :: no_nit_yr_g      !<
    real, dimension(nlat) :: no_denit_yr_g    !<
    real, dimension(nlat) :: no_nitdenit_yr_g !<
    real, dimension(nlat) :: n2o_nit_yr_g     !<
    real, dimension(nlat) :: n2o_denit_yr_g   !<
    real, dimension(nlat) :: n2o_nitdenit_yr_g!<
    real, dimension(nlat) :: n2_denit_yr_g    !<
    real, dimension(nlat) :: nvol_yr_g        !<
    real, dimension(nlat) :: nleach_yr_g      !<
    real, dimension(nlat) :: appl_fert_yr_g   !<
    real, dimension(nlat) :: ndep_nh4_yr_g    !<
    real, dimension(nlat) :: ndep_no3_yr_g    !<
    real, dimension(nlat) :: ndemand_wp_npp_yr_g     !<
    real, dimension(nlat) :: nuptake_p_nh4_yr_g      !<
    real, dimension(nlat) :: nuptake_p_no3_yr_g      !<
    real, dimension(nlat) :: nuptake_a_actl_nh4_yr_g !<
    real, dimension(nlat) :: nuptake_a_actl_no3_yr_g !<
    real, dimension(nlat) :: nuptake_yr_g !<
    real, dimension(nlat) :: nleafns2s_yr_g     !<
    real, dimension(nlat) :: nstemns2s_yr_g     !<
    real, dimension(nlat) :: nrootns2s_yr_g     !<
    real, dimension(nlat) :: nalloc_l_yr_g      !<
    real, dimension(nlat) :: nalloc_s_yr_g      !<
    real, dimension(nlat) :: nalloc_r_yr_g      !<
    real, dimension(nlat) :: nresorped_s_yr_g   !<
    real, dimension(nlat) :: nresorped_r_yr_g   !<
    real, dimension(nlat) :: nre_alloc_s2l_yr_g !<
    real, dimension(nlat) :: nre_alloc_r2l_yr_g !<
    real, dimension(nlat) :: nlitr_l_yr_g       !<
    real, dimension(nlat) :: nlitr_s_yr_g       !<
    real, dimension(nlat) :: nlitr_r_yr_g       !<
    real, dimension(nlat) :: nlitr_yr_g         !<
    real, dimension(nlat) :: gl2bl_grass_nflux_yr_g !<
    real, dimension(nlat) :: c2n_l_yr_g         !<
    real, dimension(nlat) :: c2n_s_yr_g         !<
    real, dimension(nlat) :: c2n_r_yr_g         !<
    real, dimension(nlat) :: c2n_wp_yr_g        !<
    real, dimension(nlat) :: c2n_litr_yr_g      !<
    real, dimension(nlat) :: c2n_humus_yr_g     !<
    real, dimension(nlat) :: nhumtrs_yr_g       !<
    real, dimension(nlat) :: nmineral_litr_yr_g !<
    real, dimension(nlat) :: nmineral_humus_yr_g!<
    real, dimension(nlat) :: nimmobil_nh4_yr_g  !<
    real, dimension(nlat) :: nimmobil_no3_yr_g  !<
    real, dimension(nlat) :: netnmineral_yr_g!<
    real, dimension(nlat) :: fNnetland_yr_g     !<

  end type ctem_gridavg_annual

  type (ctem_gridavg_annual), save, target :: ctem_grd_yr

  !=================================================================================
  !> CTEM's variables per tile annual values
  type ctem_tileavg_annual

    ! c      Annual output for CTEM mosaic variables:
    ! c      (denoted by name ending in "_yr_m")
    !
    ! allocated with nlat, nmos:
    real, dimension(nlat,nmos) :: laimaxg_yr_t  !<
    real, dimension(nlat,nmos) :: gleafmas_yr_t   !<
    real, dimension(nlat,nmos) :: gleafmas_NS_yr_t !<
    real, dimension(nlat,nmos) :: gleafmass_yr_t  !<
    real, dimension(nlat,nmos) :: bleafmas_yr_t   !<
    real, dimension(nlat,nmos) :: stemmass_yr_t !<
    real, dimension(nlat,nmos) :: stemmass_NS_yr_t !<
    real, dimension(nlat,nmos) :: stemmasss_yr_t  !<
    real, dimension(nlat,nmos) :: rootmass_yr_t !<
    real, dimension(nlat,nmos) :: rootmass_NS_yr_t !<
    real, dimension(nlat,nmos) :: rootmasss_yr_t  !<
    real, dimension(nlat,nmos) :: npp_yr_t      !<
    real, dimension(nlat,nmos) :: gpp_yr_t      !<
    real, dimension(nlat,nmos) :: vcmax0_yr_t     !<
    real, dimension(nlat,nmos) :: leafns2s_yr_t      !<
    real, dimension(nlat,nmos) :: stemns2s_yr_t      !<
    real, dimension(nlat,nmos) :: rootns2s_yr_t      !<
    real, dimension(nlat,nmos) :: re_alloc_s2l_yr_t  !<
    real, dimension(nlat,nmos) :: re_alloc_r2l_yr_t  !<
    real, dimension(nlat,nmos) :: re_alloc_sr2l_yr_t !<
    real, dimension(nlat,nmos) :: vgbiomas_yr_t !<
    real, dimension(nlat,nmos) :: autores_yr_t  !<
    real, dimension(nlat,nmos) :: rmrveg_yr_t  !<
    real, dimension(nlat,nmos) :: totcmass_yr_t !<
    real, dimension(nlat,nmos) :: nh4_mass_yr_t    !<
    real, dimension(nlat,nmos) :: no3_mass_yr_t    !<
    real, dimension(nlat,nmos) :: ngleafmas_yr_t   !<
    real, dimension(nlat,nmos) :: ngleafmas_NS_yr_t !<
    real, dimension(nlat,nmos) :: ngleafmass_yr_t  !<
    real, dimension(nlat,nmos) :: nbleafmas_yr_t   !<
    real, dimension(nlat,nmos) :: nstemmass_yr_t   !<
    real, dimension(nlat,nmos) :: nstemmass_NS_yr_t !<
    real, dimension(nlat,nmos) :: nstemmasss_yr_t  !<
    real, dimension(nlat,nmos) :: nrootmass_yr_t   !<
    real, dimension(nlat,nmos) :: nrootmass_NS_yr_t !<
    real, dimension(nlat,nmos) :: nrootmasss_yr_t  !<
    real, dimension(nlat,nmos) :: nlitrmass_yr_t   !<
    real, dimension(nlat,nmos) :: soilnmas_yr_t    !<
    real, dimension(nlat,nmos) :: nvgbiomas_yr_t   !<
    real, dimension(nlat,nmos) :: nep_yr_t      !<
    real, dimension(nlat,nmos) :: nepCMIP_yr_t      !<
    real, dimension(nlat,nmos) :: litrfall_yr_t   !<
    real, dimension(nlat,nmos) :: hetrores_yr_t !<
    real, dimension(nlat,nmos) :: nbp_yr_t      !<
    real, dimension(nlat,nmos) :: emit_co2_yr_t !<
    real, dimension(nlat,nmos) :: emit_co_yr_t  !<
    real, dimension(nlat,nmos) :: emit_ch4_yr_t !<
    real, dimension(nlat,nmos) :: emit_nmhc_yr_t !<
    real, dimension(nlat,nmos) :: emit_h2_yr_t  !<
    real, dimension(nlat,nmos) :: emit_nox_yr_t !<
    real, dimension(nlat,nmos) :: emit_n2o_yr_t !<
    real, dimension(nlat,nmos) :: emit_nh3_yr_t !<
    real, dimension(nlat,nmos) :: emit_pm25_yr_t !<
    real, dimension(nlat,nmos) :: emit_tpm_yr_t !<
    real, dimension(nlat,nmos) :: emit_tc_yr_t  !<
    real, dimension(nlat,nmos) :: emit_oc_yr_t  !<
    real, dimension(nlat,nmos) :: emit_bc_yr_t  !<
    real, dimension(nlat,nmos) :: burnfrac_yr_t !<
    real, dimension(nlat,nmos) :: timharvarea_yr_t !<
    real, dimension(nlat,nmos) :: tileAge_yr_t !<
    real, dimension(nlat,nmos) :: smfuncveg_yr_t !<
    real, dimension(nlat,nmos) :: bterm_yr_t    !<
    real, dimension(nlat,nmos) :: luc_emc_yr_t  !<
    real, dimension(nlat,nmos) :: lterm_yr_t    !<
    real, dimension(nlat,nmos) :: lucsocin_yr_t !<
    real, dimension(nlat,nmos) :: mterm_yr_t    !<
    real, dimension(nlat,nmos) :: lucltrin_yr_t !<
    real, dimension(nlat,nmos) :: luc_emcn_yr_t  !<
    real, dimension(nlat,nmos) :: lucltrinn_yr_t !<
    real, dimension(nlat,nmos) :: lucsocinn_yr_t !<
    real, dimension(nlat,nmos) :: ch4WetSpec_yr_t  !<
    real, dimension(nlat,nmos) :: wetfdyn_yr_t  !<
    real, dimension(nlat,nmos) :: ch4WetDyn_yr_t  !<
    real, dimension(nlat,nmos) :: ch4soills_yr_t !<
    real, dimension(nlat,nmos) :: veghght_yr_t  !<
    real, dimension(nlat,nmos) :: peatdep_yr_t  !<
    real, dimension(nlat,nmos) :: peatSoilC_yr_t  !<
    real, dimension(nlat,nmos) :: fProductDecomp_yr_t    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, nmos, ignd:
    real, dimension(nlat,nmos,ignd) :: litrmass_yr_t !<
    real, dimension(nlat,nmos,ignd) :: soilcmas_yr_t !<
    real, dimension(nlat,nmos,ignd) :: litres_yr_t   !<
    real, dimension(nlat,nmos,ignd) :: soilcres_yr_t !<
    real, dimension(nlat,nmos) :: bnf_tot_yr_t  !<
    real, dimension(nlat,nmos) :: bnf_free_yr_t  !<
    real, dimension(nlat,nmos) :: bnf_ant_yr_t  !<
    real, dimension(nlat,nmos) :: bnf_nat_yr_t  !<
    real, dimension(nlat,nmos) :: nstress_yr_t
    real, dimension(nlat,nmos) :: nitrif_yr_t      !<
    real, dimension(nlat,nmos) :: no_nit_yr_t      !<
    real, dimension(nlat,nmos) :: no_denit_yr_t    !<
    real, dimension(nlat,nmos) :: no_nitdenit_yr_t !<
    real, dimension(nlat,nmos) :: n2o_nit_yr_t     !<
    real, dimension(nlat,nmos) :: n2o_denit_yr_t   !<
    real, dimension(nlat,nmos) :: n2o_nitdenit_yr_t!<
    real, dimension(nlat,nmos) :: n2_denit_yr_t    !<
    real, dimension(nlat,nmos) :: nvol_yr_t        !<
    real, dimension(nlat,nmos) :: nleach_yr_t      !<
    real, dimension(nlat,nmos) :: appl_fert_yr_t   !<
    real, dimension(nlat,nmos) :: ndep_nh4_yr_t    !<
    real, dimension(nlat,nmos) :: ndep_no3_yr_t    !<
    real, dimension(nlat,nmos) :: ndemand_wp_npp_yr_t     !<
    real, dimension(nlat,nmos) :: nuptake_p_nh4_yr_t      !<
    real, dimension(nlat,nmos) :: nuptake_p_no3_yr_t      !<
    real, dimension(nlat,nmos) :: nuptake_a_actl_nh4_yr_t !<
    real, dimension(nlat,nmos) :: nuptake_a_actl_no3_yr_t !<
    real, dimension(nlat,nmos) :: nuptake_yr_t            !<
    real, dimension(nlat,nmos) :: nleafns2s_yr_t     !<
    real, dimension(nlat,nmos) :: nstemns2s_yr_t     !<
    real, dimension(nlat,nmos) :: nrootns2s_yr_t     !<
    real, dimension(nlat,nmos) :: nalloc_l_yr_t      !<
    real, dimension(nlat,nmos) :: nalloc_s_yr_t      !<
    real, dimension(nlat,nmos) :: nalloc_r_yr_t      !<
    real, dimension(nlat,nmos) :: nresorped_s_yr_t   !<
    real, dimension(nlat,nmos) :: nresorped_r_yr_t   !<
    real, dimension(nlat,nmos) :: nre_alloc_s2l_yr_t !<
    real, dimension(nlat,nmos) :: nre_alloc_r2l_yr_t !<
    real, dimension(nlat,nmos) :: nlitr_l_yr_t       !<
    real, dimension(nlat,nmos) :: nlitr_s_yr_t       !<
    real, dimension(nlat,nmos) :: nlitr_r_yr_t       !<
    real, dimension(nlat,nmos) :: nlitr_yr_t         !<
    real, dimension(nlat,nmos) :: gl2bl_grass_nflux_yr_t !<
    real, dimension(nlat,nmos) :: c2n_l_yr_t         !<
    real, dimension(nlat,nmos) :: c2n_s_yr_t         !<
    real, dimension(nlat,nmos) :: c2n_r_yr_t         !<
    real, dimension(nlat,nmos) :: c2n_wp_yr_t        !<
    real, dimension(nlat,nmos) :: c2n_litr_yr_t      !<
    real, dimension(nlat,nmos) :: c2n_humus_yr_t     !<
    real, dimension(nlat,nmos) :: nhumtrs_yr_t       !<
    real, dimension(nlat,nmos) :: nmineral_litr_yr_t !<
    real, dimension(nlat,nmos) :: nmineral_humus_yr_t!<
    real, dimension(nlat,nmos) :: nimmobil_nh4_yr_t  !<
    real, dimension(nlat,nmos) :: nimmobil_no3_yr_t  !<
    real, dimension(nlat,nmos) :: netnmineral_yr_t   !<
    real, dimension(nlat,nmos) :: fNnetland_yr_t     !<

  end type ctem_tileavg_annual

  type (ctem_tileavg_annual), save, target :: ctem_tile_yr

contains

  ! -----------------------------------------------------

  !> \ingroup ctemstatevars_initRowVarsBioGeoChem
  !! @{
  !> Initializes 'row' variables
  subroutine initRowVarsBioGeoChem

    implicit none

    vrot%co2conc = 0.0
    vrot%npp = 0.0
    vrot%nep = 0.0
    vrot%nepCMIP = 0.0
    vrot%hetrores = 0.0
    vrot%autores = 0.0
    vrot%soilcresp = 0.0
    vrot%rm = 0.0
    vrot%rg = 0.0
    vrot%nbp = 0.0
    vrot%litres = 0.0
    vrot%socres = 0.0
    vrot%gpp = 0.0
    vrot%dstcemls = 0.0
    vrot%dstcemls3 = 0.0
    vrot%litrfall = 0.0
    vrot%humiftrs = 0.0
    vrot%canres = 0.0
    vrot%rml = 0.0
    vrot%rms = 0.0
    vrot%rmr = 0.0
    vrot%lucemcom = 0.0
    vrot%lucltrin = 0.0
    vrot%lucsocin = 0.0
    vrot%lucemcomn = 0.0
    vrot%lucltrinn = 0.0
    vrot%lucsocinn = 0.0
    vrot%burnfrac = 0.0
    vrot%lterm = 0.0
    !vrot%cfluxcg = 0.0
    !vrot%cfluxcs = 0.0
    vrot%ROFB = 0.0
    vrot%ch4WetSpec = 0.0
    vrot%wetfdyn = 0.0
    vrot%ch4WetDyn = 0.0
    vrot%ch4_soills = 0.0
    vrot%nppmoss = 0.0
    vrot%rmlmoss = 0.0
    vrot%gppmoss = 0.0
    vrot%anmoss = 0.0
    vrot%armoss = 0.0
    vrot%peatdep = 0.0
    vrot%peatSoilC = 0.0
    vrot%pdd = 0.0
    vrot%ZOLNC = 0.0
    vrot%AILC = 0.0
    vrot%CMASVEGC = 0.0
    vrot%ALVSCTM = 0.0
    vrot%ALIRCTM = 0.0
    vrot%PAIC = 0.0
    vrot%SLAIC = 0.0
    vrot%RMATC = 0.0
    vrot%smfuncveg = 0.0
    vrot%gleafmas = 0.0
    vrot%gleafmas_ns = 0.0
    vrot%gleafmas_s = 0.0
    vrot%bleafmas = 0.0
    vrot%stemmass = 0.0
    vrot%stemmass_ns = 0.0
    vrot%stemmass_s = 0.0
    vrot%rootmass = 0.0
    vrot%rootmass_ns = 0.0
    vrot%rootmass_s = 0.0
    vrot%pstemmass = 0.0
    vrot%pgleafmass = 0.0
    vrot%litrfallveg = 0.0
    vrot%bterm = 0.0
    vrot%mterm = 0.0
    vrot%ailcg = 0.0
    vrot%ailcgs = 0.0
    vrot%fcancs = 0.0
    vrot%fcanc = 0.0
    vrot%fcancmx = 0.0
    !vrot%co2i1cg = 0.0
    !vrot%co2i1cs = 0.0
    !vrot%co2i2cg = 0.0
    !vrot%co2i2cs = 0.0
    vrot%ancsveg = 0.0
    vrot%ancgveg = 0.0
    vrot%rmlcsveg = 0.0
    vrot%rmlcgveg = 0.0
    vrot%ailcb = 0.0
    vrot%grwtheff = 0.0
    vrot%bmasveg = 0.0
    vrot%tltrleaf = 0.0
    vrot%tltrstem = 0.0
    vrot%tltrroot = 0.0
    vrot%leaflitr = 0.0
    vrot%roottemp = 0.0
    vrot%afrleaf = 0.0
    vrot%afrstem = 0.0
    vrot%afrroot = 0.0
    vrot%wtstatus = 0.0
    vrot%ltstatus = 0.0
    vrot%pfcancmx = 0.0
    vrot%nfcancmx = 0.0
    vrot%nppveg = 0.0
    vrot%veghght = 0.0
    vrot%rootdpth = 0.0
    vrot%anveg = 0.0
    vrot%rmlveg = 0.0
    vrot%rmlvegacc = 0.0
    vrot%rmsveg = 0.0
    vrot%rmrveg = 0.0
    vrot%rgveg = 0.0
    vrot%vgbiomas_veg = 0.0
    vrot%gppveg = 0.0
    vrot%vcmax0 = 0.0
    vrot%autoresveg = 0.0
    vrot%emit_co2 = 0.0
    vrot%emit_co = 0.0
    vrot%emit_ch4 = 0.0
    vrot%emit_nmhc = 0.0
    vrot%emit_h2 = 0.0
    vrot%emit_nox = 0.0
    vrot%emit_n2o = 0.0
    vrot%emit_nh3 = 0.0
    vrot%emit_pm25 = 0.0
    vrot%emit_tpm = 0.0
    vrot%emit_tc = 0.0
    vrot%emit_oc = 0.0
    vrot%emit_bc = 0.0
    vrot%burnvegf = 0.0
    vrot%leafns2s = 0.0
    vrot%stemns2s = 0.0
    vrot%rootns2s = 0.0
    vrot%re_alloc_s2l = 0.0
    vrot%re_alloc_r2l = 0.0
    vrot%re_alloc_sr2l = 0.0
    vrot%ngleafmas = 0.0
    vrot%ngleafmas_ns = 0.0
    vrot%ngleafmas_s = 0.0
    vrot%nbleafmas = 0.0
    vrot%nstemmass = 0.0
    vrot%nstemmass_ns = 0.0
    vrot%nstemmass_s = 0.0
    vrot%nrootmass = 0.0
    vrot%nrootmass_ns = 0.0
    vrot%nrootmass_s = 0.0
    vrot%nvgbiomas_veg = 0.0
    vrot%ndemandveg_wp_npp = 0.0
    vrot%nuptakeveg_p_nh4 = 0.0
    vrot%nuptakeveg_p_no3 = 0.0
    vrot%nuptakeveg_a_actl_nh4 = 0.0
    vrot%nuptakeveg_a_actl_no3 = 0.0
    vrot%nuptakeveg = 0.0
    vrot%nallocveg_l = 0.0
    vrot%nallocveg_s = 0.0
    vrot%nallocveg_r = 0.0
    vrot%nresorpedveg_s = 0.0
    vrot%nresorpedveg_r = 0.0
    vrot%nre_allocveg_s2l = 0.0
    vrot%nre_allocveg_r2l = 0.0
    vrot%nleafns2sveg = 0.0
    vrot%nstemns2sveg = 0.0
    vrot%nrootns2sveg = 0.0
    vrot%gl2bl_grass_nflux = 0.0
    vrot%nlitrveg_l = 0.0
    vrot%nlitrveg_s = 0.0
    vrot%nlitrveg_r = 0.0
    vrot%nlitrveg   = 0.0
    vrot%c2nveg_l = 0.0
    vrot%c2nveg_s = 0.0
    vrot%c2nveg_r = 0.0
    vrot%c2nveg_wp = 0.0
    vrot%rmatctem = 0.0
    vrot%hetroresveg = 0.0
    vrot%nepveg = 0.0
    vrot%nbpveg = 0.0
    vrot%litrmass = 0.0
    vrot%soilcmas = 0.0
    vrot%litresveg = 0.0
    vrot%soilcresveg = 0.0
    vrot%humiftrsveg = 0.0
    vrot%bnf_free = 0.0
    vrot%bnf_ant = 0.0
    vrot%bnf_nat = 0.0
    vrot%bnf_tot = 0.0
    vrot%nstress = 0.0
    vrot%nitrifveg = 0.0
    vrot%no_nitveg = 0.0
    vrot%no_denitveg = 0.0
    vrot%no_nitdenitveg = 0.0
    vrot%n2o_nitveg = 0.0
    vrot%n2o_denitveg = 0.0
    vrot%n2o_nitdenitveg = 0.0
    vrot%n2_denitveg = 0.0
    vrot%nvolveg = 0.0
    vrot%nleachveg = 0.0
    vrot%nh4_mass = 0.0
    vrot%no3_mass = 0.0
    vrot%nlitrmass = 0.0
    vrot%soilnmas = 0.0
    vrot%appl_fert = 0.0
    vrot%ndep_nh4 = 0.0
    vrot%ndep_no3 = 0.0
    vrot%c2nveg_litr = 0.0
    vrot%c2nveg_humus = 0.0
    vrot%nhumtrsveg = 0.0
    vrot%nmineralveg_litr = 0.0
    vrot%nmineralveg_humus = 0.0
    vrot%nimmobilveg_nh4 = 0.0
    vrot%nimmobilveg_no3 = 0.0
    vrot%netnmineralveg  = 0.0
    vrot%fNnetlandveg = 0.0
    vrot%redcoeff_vcmax = 0.0
    vrot%twarmm = 0.0
    vrot%tcoldm = 0.0
    vrot%gdd5 = 0.0
    vrot%aridity = 0.0
    vrot%srplsmon = 0.0
    vrot%defctmon = 0.0
    vrot%anndefct = 0.0
    vrot%annsrpls = 0.0
    vrot%annpcp = 0.0
    vrot%dry_season_length = 0.0
    vrot%tileAgerow = 0.0

  end subroutine initRowVarsBioGeoChem
  !! @}

  !==================================================

  !> \ingroup ctemstatevars_resetMonthEnd
  !! @{
  !> Resets monthly variables at month end in preparation for next month
  subroutine resetMonthEnd (nltest, nmtest)

    use classicParams,    only : iccp2, icc, iccp1, ignd

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i, m, j

    ! These are assigned to mid-month, but are not accumulated so can be
    ! zeroed out at the same time as the other month-end vars.
    do i = 1,nltest
      ctem_grd_mo%gleafmas_mo_g(i)  = 0.0
      ctem_grd_mo%gleafmas_NS_mo_g(i) = 0.0
      ctem_grd_mo%gleafmass_mo_g(i) = 0.0
      ctem_grd_mo%stemmass_mo_g(i) = 0.0
      ctem_grd_mo%stemmass_NS_mo_g(i) = 0.0
      ctem_grd_mo%stemmasss_mo_g(i) = 0.0
      ctem_grd_mo%rootmass_mo_g(i) = 0.0
      ctem_grd_mo%rootmass_NS_mo_g(i) = 0.0
      ctem_grd_mo%rootmasss_mo_g(i) = 0.0
      ctem_grd_mo%litrmass_mo_g(i,1:ignd)=0.0
      ctem_grd_mo%soilcmas_mo_g(i,1:ignd)=0.0
      ctem_grd_mo%vgbiomas_mo_g(i) = 0.0
      ctem_grd_mo%totcmass_mo_g(i) = 0.0
      ctem_grd_mo%nh4_mass_mo_g(i) = 0.0
      ctem_grd_mo%no3_mass_mo_g(i) = 0.0
      ctem_grd_mo%ngleafmas_mo_g(i) = 0.0
      ctem_grd_mo%ngleafmas_NS_mo_g(i) = 0.0
      ctem_grd_mo%ngleafmass_mo_g(i) = 0.0
      ctem_grd_mo%nbleafmas_mo_g(i) = 0.0
      ctem_grd_mo%nstemmass_mo_g(i) = 0.0
      ctem_grd_mo%nstemmass_NS_mo_g(i) = 0.0
      ctem_grd_mo%nstemmasss_mo_g(i) = 0.0
      ctem_grd_mo%nrootmass_mo_g(i) = 0.0
      ctem_grd_mo%nrootmass_NS_mo_g(i) = 0.0
      ctem_grd_mo%nrootmasss_mo_g(i) = 0.0
      ctem_grd_mo%nvgbiomas_mo_g(i) = 0.0
      ctem_grd_mo%ndemand_wp_npp_mo_g(i) = 0.0
      ctem_grd_mo%c2n_l_mo_g(i) = 0.0
      ctem_grd_mo%c2n_s_mo_g(i) = 0.0
      ctem_grd_mo%c2n_r_mo_g(i) = 0.0
      ctem_grd_mo%c2n_wp_mo_g(i) = 0.0
      ctem_grd_mo%c2n_litr_mo_g(i) = 0.0
      ctem_grd_mo%c2n_humus_mo_g(i) = 0.0
      ctem_grd_mo%nlitrmass_mo_g(i) = 0.0
      ctem_grd_mo%soilnmas_mo_g(i) = 0.0
      do m = 1,nmtest
        ctem_tile_mo%gleafmas_mo_t(i,m)  = 0.0
        ctem_tile_mo%gleafmas_NS_mo_t(i,m)= 0.0
        ctem_tile_mo%gleafmass_mo_t(i,m) = 0.0
        ctem_tile_mo%stemmass_mo_t(i,m) = 0.0
        ctem_tile_mo%stemmass_NS_mo_t(i,m) = 0.0
        ctem_tile_mo%stemmasss_mo_t(i,m) = 0.0
        ctem_tile_mo%rootmass_mo_t(i,m) = 0.0
        ctem_tile_mo%rootmass_NS_mo_t(i,m) = 0.0
        ctem_tile_mo%rootmasss_mo_t(i,m) = 0.0
        ctem_tile_mo%litrmass_mo_t(i,m,1:ignd)=0.0
        ctem_tile_mo%soilcmas_mo_t(i,m,1:ignd)=0.0
        ctem_tile_mo%vgbiomas_mo_t(i,m) = 0.0
        ctem_tile_mo%totcmass_mo_t(i,m) = 0.0
        ctem_tile_mo%nh4_mass_mo_t(i,m) = 0.0
        ctem_tile_mo%no3_mass_mo_t(i,m) = 0.0
        ctem_tile_mo%ngleafmas_mo_t(i,m) = 0.0
        ctem_tile_mo%ngleafmas_NS_mo_t(i,m) = 0.0
        ctem_tile_mo%ngleafmass_mo_t(i,m) = 0.0
        ctem_tile_mo%nbleafmas_mo_t(i,m) = 0.0
        ctem_tile_mo%nstemmass_mo_t(i,m) = 0.0
        ctem_tile_mo%nstemmass_NS_mo_t(i,m) = 0.0
        ctem_tile_mo%nstemmasss_mo_t(i,m) = 0.0
        ctem_tile_mo%nrootmass_mo_t(i,m) = 0.0
        ctem_tile_mo%nrootmass_NS_mo_t(i,m) = 0.0
        ctem_tile_mo%nrootmasss_mo_t(i,m) = 0.0
        ctem_tile_mo%ndemand_wp_npp_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_l_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_s_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_r_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_wp_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_litr_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_humus_mo_t(i,m) = 0.0
        ctem_tile_mo%nlitrmass_mo_t(i,m) = 0.0
        ctem_tile_mo%soilnmas_mo_t(i,m) = 0.0
        ctem_tile_mo%nvgbiomas_mo_t(i,m) = 0.0
        do j = 1,icc
          ctem_mo%gleafmas_mo(i,m,j) = 0.0
          ctem_mo%gleafmas_NS_mo(i,m,j) = 0.0
          ctem_mo%gleafmass_mo(i,m,j) = 0.0
          ctem_mo%stemmass_mo(i,m,j) = 0.0
          ctem_mo%stemmass_NS_mo(i,m,j) = 0.0
          ctem_mo%stemmasss_mo(i,m,j) = 0.0
          ctem_mo%rootmass_mo(i,m,j) = 0.0
          ctem_mo%rootmass_NS_mo(i,m,j) = 0.0
          ctem_mo%rootmasss_mo(i,m,j) = 0.0
          ctem_mo%rootdpth_mo(i,m,j) = 0.0
          ctem_mo%litrmass_mo(i,m,j,1:ignd)=0.0
          ctem_mo%soilcmas_mo(i,m,j,1:ignd)=0.0
          ctem_mo%vgbiomas_mo(i,m,j) = 0.0
          ctem_mo%totcmass_mo(i,m,j) = 0.0
          ctem_mo%nh4_mass_mo(i,m,j) = 0.0
          ctem_mo%no3_mass_mo(i,m,j) = 0.0
          ctem_mo%ngleafmas_mo(i,m,j) = 0.0
          ctem_mo%ngleafmas_NS_mo(i,m,j)= 0.0
          ctem_mo%ngleafmass_mo(i,m,j) = 0.0
          ctem_mo%nbleafmas_mo(i,m,j) = 0.0
          ctem_mo%nstemmass_mo(i,m,j) = 0.0
          ctem_mo%nstemmass_NS_mo(i,m,j)= 0.0
          ctem_mo%nstemmasss_mo(i,m,j) = 0.0
          ctem_mo%nrootmass_mo(i,m,j) = 0.0
          ctem_mo%nrootmass_NS_mo(i,m,j)= 0.0
          ctem_mo%nrootmasss_mo(i,m,j) = 0.0
          ctem_mo%ndemand_wp_npp_mo(i,m,j) = 0.0
          ctem_mo%c2n_l_mo(i,m,j) = 0.0
          ctem_mo%c2n_s_mo(i,m,j) = 0.0
          ctem_mo%c2n_r_mo(i,m,j) = 0.0
          ctem_mo%c2n_wp_mo(i,m,j) = 0.0
          ctem_mo%c2n_litr_mo(i,m,j) = 0.0
          ctem_mo%c2n_humus_mo(i,m,j) = 0.0
          ctem_mo%nlitrmass_mo(i,m,j) = 0.0
          ctem_mo%soilnmas_mo(i,m,j) = 0.0
          ctem_mo%nvgbiomas_mo(i,m,j) = 0.0
        end do
        ctem_mo%totcmass_mo(i,m,iccp1) = 0.0
        ctem_mo%litrmass_mo(i,m,iccp1,1:ignd)=0.0
        ctem_mo%soilcmas_mo(i,m,iccp1,1:ignd)=0.0
        ctem_mo%litrmass_mo(i,m,iccp2,1:ignd)=0.0
        ctem_mo%soilcmas_mo(i,m,iccp2,1:ignd)=0.0
        ctem_mo%nh4_mass_mo(i,m,iccp1) = 0.0
        ctem_mo%no3_mass_mo(i,m,iccp1) = 0.0
        ctem_mo%nlitrmass_mo(i,m,iccp1)= 0.0
        ctem_mo%soilnmas_mo(i,m,iccp1) = 0.0
        ctem_mo%nlitrmass_mo(i,m,iccp2)= 0.0
        ctem_mo%soilnmas_mo(i,m,iccp2) = 0.0
        ctem_mo%c2n_litr_mo(i,m,iccp1) = 0.0
        ctem_mo%c2n_humus_mo(i,m,iccp1)= 0.0
      end do
    end do

    ! Now zero out the month end vars.
    do i = 1,nltest
      ! Grid avg
      ctem_grd_mo%laimaxg_mo_g(i) = 0.0
      ctem_grd_mo%npp_mo_g(i) = 0.0
      ctem_grd_mo%gpp_mo_g(i) = 0.0
      ctem_grd_mo%vcmax0_mo_g(i) = 0.0
      ctem_grd_mo%leafns2s_mo_g(i) = 0.0
      ctem_grd_mo%stemns2s_mo_g(i) = 0.0
      ctem_grd_mo%rootns2s_mo_g(i) = 0.0
      ctem_grd_mo%re_alloc_s2l_mo_g(i) = 0.0
      ctem_grd_mo%re_alloc_r2l_mo_g(i) = 0.0
      ctem_grd_mo%re_alloc_sr2l_mo_g(i)= 0.0
      ctem_grd_mo%nep_mo_g(i) = 0.0
      ctem_grd_mo%nepCMIP_mo_g(i) = 0.0
      ctem_grd_mo%nbp_mo_g(i) = 0.0
      ctem_grd_mo%hetrores_mo_g(i) = 0.0
      ctem_grd_mo%autores_mo_g(i) = 0.0
      ctem_grd_mo%soilres_mo_g(i) = 0.0
      ctem_grd_mo%litres_mo_g(i,1:ignd)=0.0
      ctem_grd_mo%soilcres_mo_g(i,1:ignd)=0.0
      ctem_grd_mo%litrfall_mo_g(i) = 0.0
      ctem_grd_mo%humiftrs_mo_g(i) = 0.0
      ctem_grd_mo%tltrleaf_mo_g(i) = 0.0
      ctem_grd_mo%tltrstem_mo_g(i) = 0.0
      ctem_grd_mo%tltrroot_mo_g(i) = 0.0
      ctem_grd_mo%emit_co2_mo_g(i) = 0.0
      ctem_grd_mo%emit_co_mo_g(i) = 0.0
      ctem_grd_mo%emit_ch4_mo_g(i) = 0.0
      ctem_grd_mo%emit_nmhc_mo_g(i) = 0.0
      ctem_grd_mo%emit_h2_mo_g(i) = 0.0
      ctem_grd_mo%emit_nox_mo_g(i) = 0.0
      ctem_grd_mo%emit_n2o_mo_g(i) = 0.0
      ctem_grd_mo%emit_nh3_mo_g(i) = 0.0
      ctem_grd_mo%emit_pm25_mo_g(i) = 0.0
      ctem_grd_mo%emit_tpm_mo_g(i) = 0.0
      ctem_grd_mo%emit_tc_mo_g(i) = 0.0
      ctem_grd_mo%emit_oc_mo_g(i) = 0.0
      ctem_grd_mo%emit_bc_mo_g(i) = 0.0
      ctem_grd_mo%smfuncveg_mo_g(i) = 0.0
      ctem_grd_mo%luc_emc_mo_g(i) = 0.0
      ctem_grd_mo%lucsocin_mo_g(i) = 0.0
      ctem_grd_mo%lucltrin_mo_g(i) = 0.0
      ctem_grd_mo%luc_emcn_mo_g(i) = 0.0
      ctem_grd_mo%lucsocinn_mo_g(i) = 0.0
      ctem_grd_mo%lucltrinn_mo_g(i) = 0.0
      ctem_grd_mo%burnfrac_mo_g(i) = 0.0
      ctem_grd_mo%timharvarea_mo_g(i) = 0.0  
      ctem_grd_mo%tileAge_mo_g(i) = 0.0   
      ctem_grd_mo%bterm_mo_g(i)    = 0.0
      ctem_grd_mo%lterm_mo_g(i)    = 0.0
      ctem_grd_mo%mterm_mo_g(i)    = 0.0
      ctem_grd_mo%ch4WetSpec_mo_g(i)  = 0.0
      ctem_grd_mo%wetfdyn_mo_g(i)  = 0.0
      ctem_grd_mo%wetfpres_mo_g(i)  = 0.0
      ctem_grd_mo%ch4WetDyn_mo_g(i)  = 0.0
      ctem_grd_mo%ch4soills_mo_g(i)  = 0.0
      ctem_grd_mo%cProduct_mo_g(i)  = 0.0
      ctem_grd_mo%nProduct_mo_g(i)  = 0.0
      ctem_grd_mo%fProductDecomp_mo_g(i)  = 0.0
      ctem_grd_mo%bnf_tot_mo_g(i) = 0.0
      ctem_grd_mo%bnf_free_mo_g(i) = 0.0
      ctem_grd_mo%bnf_ant_mo_g(i) = 0.0
      ctem_grd_mo%bnf_nat_mo_g(i) = 0.0
      ctem_grd_mo%nstress_mo_g(i) = 0.0
      ctem_grd_mo%nitrif_mo_g(i) = 0.0
      ctem_grd_mo%no_nit_mo_g(i) = 0.0
      ctem_grd_mo%no_denit_mo_g(i) = 0.0
      ctem_grd_mo%no_nitdenit_mo_g(i) = 0.0
      ctem_grd_mo%n2o_nit_mo_g(i) = 0.0
      ctem_grd_mo%n2o_denit_mo_g(i) = 0.0
      ctem_grd_mo%n2o_nitdenit_mo_g(i) = 0.0
      ctem_grd_mo%n2_denit_mo_g(i) = 0.0
      ctem_grd_mo%nvol_mo_g(i) = 0.0
      ctem_grd_mo%nleach_mo_g(i) = 0.0
      ctem_grd_mo%appl_fert_mo_g(i) = 0.0
      ctem_grd_mo%ndep_nh4_mo_g(i) = 0.0
      ctem_grd_mo%ndep_no3_mo_g(i) = 0.0
      ctem_grd_mo%nuptake_p_nh4_mo_g(i) = 0.0
      ctem_grd_mo%nuptake_p_no3_mo_g(i) = 0.0
      ctem_grd_mo%nuptake_a_actl_nh4_mo_g(i) = 0.0
      ctem_grd_mo%nuptake_a_actl_no3_mo_g(i) = 0.0
      ctem_grd_mo%nuptake_mo_g(i) = 0.0
      ctem_grd_mo%nalloc_l_mo_g(i) = 0.0
      ctem_grd_mo%nalloc_s_mo_g(i) = 0.0
      ctem_grd_mo%nalloc_r_mo_g(i) = 0.0
      ctem_grd_mo%nresorped_s_mo_g(i) = 0.0
      ctem_grd_mo%nresorped_r_mo_g(i) = 0.0
      ctem_grd_mo%nre_alloc_s2l_mo_g(i)= 0.0
      ctem_grd_mo%nre_alloc_r2l_mo_g(i)= 0.0
      ctem_grd_mo%nleafns2s_mo_g(i)  = 0.0
      ctem_grd_mo%nstemns2s_mo_g(i) = 0.0
      ctem_grd_mo%nrootns2s_mo_g(i) = 0.0
      ctem_grd_mo%nlitr_l_mo_g(i) = 0.0
      ctem_grd_mo%nlitr_s_mo_g(i) = 0.0
      ctem_grd_mo%nlitr_r_mo_g(i) = 0.0
      ctem_grd_mo%nlitr_mo_g(i) = 0.0
      ctem_grd_mo%gl2bl_grass_nflux_mo_g(i) = 0.0
      ! ctem_grd_mo%c2n_l_mo_g(i) = 0.0
      ! ctem_grd_mo%c2n_s_mo_g(i) = 0.0
      ! ctem_grd_mo%c2n_r_mo_g(i) = 0.0
      ! ctem_grd_mo%c2n_wp_mo_g(i) = 0.0
      ! ctem_grd_mo%c2n_litr_mo_g(i) = 0.0
      ! ctem_grd_mo%c2n_humus_mo_g(i) = 0.0
      ctem_grd_mo%nhumtrs_mo_g(i) = 0.0
      ctem_grd_mo%nmineral_litr_mo_g(i) = 0.0
      ctem_grd_mo%nmineral_humus_mo_g(i) = 0.0
      ctem_grd_mo%netnmineral_mo_g(i) = 0.0
      ctem_grd_mo%nimmobil_nh4_mo_g(i) = 0.0
      ctem_grd_mo%nimmobil_no3_mo_g(i) = 0.0
      ctem_grd_mo%fNnetland_mo_g(i) = 0.0

      do m = 1,nmtest
        ! Tile avg
        ctem_tile_mo%laimaxg_mo_t(i,m) = 0.0
        ctem_tile_mo%npp_mo_t(i,m) = 0.0
        ctem_tile_mo%gpp_mo_t(i,m) = 0.0
        ctem_tile_mo%vcmax0_mo_t(i,m) = 0.0
        ctem_tile_mo%leafns2s_mo_t(i,m) = 0.0
        ctem_tile_mo%stemns2s_mo_t(i,m) = 0.0
        ctem_tile_mo%rootns2s_mo_t(i,m) = 0.0
        ctem_tile_mo%re_alloc_s2l_mo_t(i,m) = 0.0
        ctem_tile_mo%re_alloc_r2l_mo_t(i,m) = 0.0
        ctem_tile_mo%re_alloc_sr2l_mo_t(i,m) = 0.0
        ctem_tile_mo%nep_mo_t(i,m) = 0.0
        ctem_tile_mo%nepCMIP_mo_t(i,m) = 0.0
        ctem_tile_mo%nbp_mo_t(i,m) = 0.0
        ctem_tile_mo%hetrores_mo_t(i,m) = 0.0
        ctem_tile_mo%autores_mo_t(i,m) = 0.0
        ctem_tile_mo%soilres_mo_t(i,m) = 0.0
        ctem_tile_mo%litres_mo_t(i,m,1:ignd)=0.0
        ctem_tile_mo%soilcres_mo_t(i,m,1:ignd)=0.0
        ctem_tile_mo%litrfall_mo_t(i,m) = 0.0
        ctem_tile_mo%humiftrs_mo_t(i,m) = 0.0
        ctem_tile_mo%tltrleaf_mo_t(i,m) = 0.0
        ctem_tile_mo%tltrstem_mo_t(i,m) = 0.0
        ctem_tile_mo%tltrroot_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_co2_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_co_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_ch4_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_nmhc_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_h2_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_nox_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_n2o_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_nh3_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_pm25_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_tpm_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_tc_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_oc_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_bc_mo_t(i,m) = 0.0
        ctem_tile_mo%smfuncveg_mo_t(i,m) = 0.0
        ctem_tile_mo%luc_emc_mo_t(i,m) = 0.0
        ctem_tile_mo%lucsocin_mo_t(i,m) = 0.0
        ctem_tile_mo%lucltrin_mo_t(i,m) = 0.0
        ctem_tile_mo%luc_emcn_mo_t(i,m) = 0.0
        ctem_tile_mo%lucsocinn_mo_t(i,m) = 0.0
        ctem_tile_mo%lucltrinn_mo_t(i,m) = 0.0
        ctem_tile_mo%burnfrac_mo_t(i,m) = 0.0
        ctem_tile_mo%timharvarea_mo_t(i,m) = 0.0  
        ctem_tile_mo%tileAge_mo_t(i,m) = 0.0   
        ctem_tile_mo%bterm_mo_t(i,m)    = 0.0
        ctem_tile_mo%lterm_mo_t(i,m)    = 0.0
        ctem_tile_mo%mterm_mo_t(i,m)    = 0.0
        ctem_tile_mo%ch4WetSpec_mo_t(i,m)  = 0.0
        ctem_tile_mo%wetfdyn_mo_t(i,m)  = 0.0
        ctem_tile_mo%wetfpres_mo_t(i,m)  = 0.0
        ctem_tile_mo%ch4WetDyn_mo_t(i,m)  = 0.0
        ctem_tile_mo%ch4soills_mo_t(i,m)  = 0.0
        ctem_tile_mo%wind_mo_t(i,m) = 0.0
        ctem_tile_mo%fProductDecomp_mo_t(i,m) = 0.0
        ctem_tile_mo%bnf_tot_mo_t(i,m) = 0.0
        ctem_tile_mo%bnf_free_mo_t(i,m) = 0.0
        ctem_tile_mo%bnf_ant_mo_t(i,m) = 0.0
        ctem_tile_mo%bnf_nat_mo_t(i,m) = 0.0
        ctem_tile_mo%nstress_mo_t(i,m) = 0.0
        ctem_tile_mo%nitrif_mo_t(i,m) = 0.0
        ctem_tile_mo%no_nit_mo_t(i,m) = 0.0
        ctem_tile_mo%no_denit_mo_t(i,m) = 0.0
        ctem_tile_mo%no_nitdenit_mo_t(i,m) = 0.0
        ctem_tile_mo%n2o_nit_mo_t(i,m) = 0.0
        ctem_tile_mo%n2o_denit_mo_t(i,m) = 0.0
        ctem_tile_mo%n2o_nitdenit_mo_t(i,m) = 0.0
        ctem_tile_mo%n2_denit_mo_t(i,m) = 0.0
        ctem_tile_mo%nvol_mo_t(i,m) = 0.0
        ctem_tile_mo%nleach_mo_t(i,m) = 0.0
        ctem_tile_mo%appl_fert_mo_t(i,m) = 0.0
        ctem_tile_mo%ndep_nh4_mo_t(i,m) = 0.0
        ctem_tile_mo%ndep_no3_mo_t(i,m) = 0.0
        ctem_tile_mo%nuptake_p_nh4_mo_t(i,m) = 0.0
        ctem_tile_mo%nuptake_p_no3_mo_t(i,m) = 0.0
        ctem_tile_mo%nuptake_a_actl_nh4_mo_t(i,m) = 0.0
        ctem_tile_mo%nuptake_a_actl_no3_mo_t(i,m) = 0.0
        ctem_tile_mo%nuptake_mo_t(i,m) = 0.0
        ctem_tile_mo%nalloc_l_mo_t(i,m) = 0.0
        ctem_tile_mo%nalloc_s_mo_t(i,m) = 0.0
        ctem_tile_mo%nalloc_r_mo_t(i,m) = 0.0
        ctem_tile_mo%nresorped_s_mo_t(i,m) = 0.0
        ctem_tile_mo%nresorped_r_mo_t(i,m) = 0.0
        ctem_tile_mo%nre_alloc_s2l_mo_t(i,m) = 0.0
        ctem_tile_mo%nre_alloc_r2l_mo_t(i,m) = 0.0
        ctem_tile_mo%nleafns2s_mo_t(i,m) = 0.0
        ctem_tile_mo%nstemns2s_mo_t(i,m) = 0.0
        ctem_tile_mo%nrootns2s_mo_t(i,m) = 0.0
        ctem_tile_mo%nlitr_l_mo_t(i,m) = 0.0
        ctem_tile_mo%nlitr_s_mo_t(i,m) = 0.0
        ctem_tile_mo%nlitr_r_mo_t(i,m) = 0.0
        ctem_tile_mo%nlitr_mo_t(i,m) = 0.0
        ctem_tile_mo%gl2bl_grass_nflux_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_l_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_s_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_r_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_wp_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_litr_mo_t(i,m) = 0.0
        ctem_tile_mo%c2n_humus_mo_t(i,m) = 0.0
        ctem_tile_mo%nhumtrs_mo_t(i,m) = 0.0
        ctem_tile_mo%nmineral_litr_mo_t(i,m) = 0.0
        ctem_tile_mo%nmineral_humus_mo_t(i,m) = 0.0
        ctem_tile_mo%netnmineral_mo_t(i,m) = 0.0
        ctem_tile_mo%nimmobil_nh4_mo_t(i,m) = 0.0
        ctem_tile_mo%nimmobil_no3_mo_t(i,m) = 0.0
        ctem_tile_mo%fNnetland_mo_t(i,m) = 0.0
        ctem_mo%timharvarea_mo(i,m) = 0.0  
        ctem_mo%tileAge_mo(i,m) = 0.0   

        do j = 1,icc
          ! per pft
          ctem_mo%laimaxg_mo(i,m,j) = 0.0
          ctem_mo%cc_mo(i,m,j) = 0.0
          ctem_mo%mm_mo(i,m,j) = 0.0
          ctem_mo%npp_mo(i,m,j) = 0.0
          ctem_mo%gpp_mo(i,m,j) = 0.0
          ctem_mo%vcmax0_mo(i,m,j) = 0.0
          ctem_mo%leafns2s_mo(i,m,j) = 0.0
          ctem_mo%stemns2s_mo(i,m,j) = 0.0
          ctem_mo%rootns2s_mo(i,m,j) = 0.0
          ctem_mo%re_alloc_s2l_mo(i,m,j) = 0.0
          ctem_mo%re_alloc_r2l_mo(i,m,j) = 0.0
          ctem_mo%re_alloc_sr2l_mo(i,m,j) = 0.0
          ctem_mo%nep_mo(i,m,j) = 0.0
          ctem_mo%nepCMIP_mo(i,m,j) = 0.0
          ctem_mo%nbp_mo(i,m,j) = 0.0
          ctem_mo%hetrores_mo(i,m,j) = 0.0
          ctem_mo%autores_mo(i,m,j) = 0.0
          ctem_mo%soilres_mo(i,m,j) = 0.0
          ctem_mo%litres_mo(i,m,j,1:ignd)=0.0
          ctem_mo%soilcres_mo(i,m,j,1:ignd)=0.0
          ctem_mo%litrfallveg_mo(i,m,j) = 0.0
          ctem_mo%humiftrsveg_mo(i,m,j) = 0.0
          ctem_mo%tltrleaf_mo(i,m,j) = 0.0
          ctem_mo%tltrstem_mo(i,m,j) = 0.0
          ctem_mo%tltrroot_mo(i,m,j) = 0.0
          ctem_mo%emit_co2_mo(i,m,j) = 0.0
          ctem_mo%emit_co_mo(i,m,j) = 0.0
          ctem_mo%emit_ch4_mo(i,m,j) = 0.0
          ctem_mo%emit_nmhc_mo(i,m,j) = 0.0
          ctem_mo%emit_h2_mo(i,m,j) = 0.0
          ctem_mo%emit_nox_mo(i,m,j) = 0.0
          ctem_mo%emit_n2o_mo(i,m,j) = 0.0
          ctem_mo%emit_nh3_mo(i,m,j) = 0.0
          ctem_mo%emit_pm25_mo(i,m,j) = 0.0
          ctem_mo%emit_tpm_mo(i,m,j) = 0.0
          ctem_mo%emit_tc_mo(i,m,j) = 0.0
          ctem_mo%emit_oc_mo(i,m,j) = 0.0
          ctem_mo%emit_bc_mo(i,m,j) = 0.0
          ctem_mo%burnfrac_mo(i,m,j) = 0.0
          ctem_mo%bterm_mo(i,m,j) = 0.0
          ctem_mo%mterm_mo(i,m,j) = 0.0
          ctem_mo%smfuncveg_mo(i,m,j) = 0.0
          ctem_mo%bnf_tot_mo(i,m,j) = 0.0
          ctem_mo%bnf_free_mo(i,m,j) = 0.0
          ctem_mo%bnf_ant_mo(i,m,j) = 0.0
          ctem_mo%bnf_nat_mo(i,m,j) = 0.0
          ctem_mo%nstress_mo(i,m,j) = 0.0
          ctem_mo%nitrif_mo(i,m,j) = 0.0
          ctem_mo%no_nit_mo(i,m,j) = 0.0
          ctem_mo%no_denit_mo(i,m,j) = 0.0
          ctem_mo%no_nitdenit_mo(i,m,j) = 0.0
          ctem_mo%n2o_nit_mo(i,m,j) = 0.0
          ctem_mo%n2o_denit_mo(i,m,j) = 0.0
          ctem_mo%n2o_nitdenit_mo(i,m,j) = 0.0
          ctem_mo%n2_denit_mo(i,m,j) = 0.0
          ctem_mo%nvol_mo(i,m,j) = 0.0
          ctem_mo%nleach_mo(i,m,j) = 0.0
          ctem_mo%appl_fert_mo(i,m,j) = 0.0
          ctem_mo%ndep_nh4_mo(i,m,j) = 0.0
          ctem_mo%ndep_no3_mo(i,m,j) = 0.0
          ctem_mo%nuptake_p_nh4_mo(i,m,j) = 0.0
          ctem_mo%nuptake_p_no3_mo(i,m,j) = 0.0
          ctem_mo%nuptake_a_actl_nh4_mo(i,m,j) = 0.0
          ctem_mo%nuptake_a_actl_no3_mo(i,m,j) = 0.0
          ctem_mo%nuptake_mo(i,m,j) = 0.0
          ctem_mo%nalloc_l_mo(i,m,j) = 0.0
          ctem_mo%nalloc_s_mo(i,m,j) = 0.0
          ctem_mo%nalloc_r_mo(i,m,j) = 0.0
          ctem_mo%nresorped_s_mo(i,m,j) = 0.0
          ctem_mo%nresorped_r_mo(i,m,j) = 0.0
          ctem_mo%nre_alloc_s2l_mo(i,m,j) = 0.0
          ctem_mo%nre_alloc_r2l_mo(i,m,j) = 0.0
          ctem_mo%nleafns2s_mo(i,m,j) = 0.0
          ctem_mo%nstemns2s_mo(i,m,j) = 0.0
          ctem_mo%nrootns2s_mo(i,m,j) = 0.0
          ctem_mo%nlitr_l_mo(i,m,j) = 0.0
          ctem_mo%nlitr_s_mo(i,m,j) = 0.0
          ctem_mo%nlitr_r_mo(i,m,j) = 0.0
          ctem_mo%nlitr_mo(i,m,j) = 0.0
          ctem_mo%gl2bl_grass_nflux_mo(i,m,j) = 0.0
          ctem_mo%c2n_l_mo(i,m,j) = 0.0
          ctem_mo%c2n_s_mo(i,m,j) = 0.0
          ctem_mo%c2n_r_mo(i,m,j) = 0.0
          ctem_mo%c2n_wp_mo(i,m,j) = 0.0
          ctem_mo%c2n_litr_mo(i,m,j) = 0.0
          ctem_mo%c2n_humus_mo(i,m,j) = 0.0
          ctem_mo%nhumtrs_mo(i,m,j) = 0.0
          ctem_mo%nmineral_litr_mo(i,m,j) = 0.0
          ctem_mo%nmineral_humus_mo(i,m,j) = 0.0
          ctem_mo%netnmineral_mo(i,m,j) = 0.0
          ctem_mo%nimmobil_nh4_mo(i,m,j) = 0.0
          ctem_mo%nimmobil_no3_mo(i,m,j) = 0.0
          ctem_mo%fNnetland_mo(i,m,j) = 0.0
        end do

        ctem_mo%nep_mo(i,m,iccp1) = 0.0
        ctem_mo%nepCMIP_mo(i,m,iccp1) = 0.0
        ctem_mo%nbp_mo(i,m,iccp1) = 0.0
        ctem_mo%hetrores_mo(i,m,iccp1) = 0.0
        ctem_mo%humiftrsveg_mo(i,m,iccp1) = 0.0
        ctem_mo%humiftrsveg_mo(i,m,iccp2) = 0.0


        ctem_mo%bnf_tot_mo(i,m,iccp1) = 0.0
        ctem_mo%bnf_free_mo(i,m,iccp1) = 0.0
        ctem_mo%nitrif_mo(i,m,iccp1) = 0.0
        ctem_mo%no_nit_mo(i,m,iccp1) = 0.0
        ctem_mo%no_denit_mo(i,m,iccp1) = 0.0
        ctem_mo%no_nitdenit_mo(i,m,iccp1) = 0.0
        ctem_mo%n2o_nit_mo(i,m,iccp1) = 0.0
        ctem_mo%n2o_denit_mo(i,m,iccp1) = 0.0
        ctem_mo%n2o_nitdenit_mo(i,m,iccp1) = 0.0
        ctem_mo%n2_denit_mo(i,m,iccp1) = 0.0
        ctem_mo%nvol_mo(i,m,iccp1) = 0.0
        ctem_mo%nleach_mo(i,m,iccp1) = 0.0
        ctem_mo%appl_fert_mo(i,m,iccp1) = 0.0
        ctem_mo%ndep_nh4_mo(i,m,iccp1) = 0.0
        ctem_mo%ndep_no3_mo(i,m,iccp1) = 0.0
        ctem_mo%c2n_litr_mo(i,m,iccp1) = 0.0
        ctem_mo%c2n_humus_mo(i,m,iccp1) = 0.0
        ctem_mo%nhumtrs_mo(i,m,iccp1) = 0.0
        ctem_mo%nmineral_litr_mo(i,m,iccp1) = 0.0
        ctem_mo%nmineral_humus_mo(i,m,iccp1) = 0.0
        ctem_mo%netnmineral_mo(i,m,iccp1) = 0.0
        ctem_mo%nimmobil_nh4_mo(i,m,iccp1) = 0.0
        ctem_mo%nimmobil_no3_mo(i,m,iccp1) = 0.0
        ctem_mo%fNnetland_mo(i,m,iccp1) = 0.0

        ctem_mo%litres_mo(i,m,iccp1,1:ignd)=0.0
        ctem_mo%soilcres_mo(i,m,iccp1,1:ignd)=0.0
        ctem_mo%humiftrsveg_mo(i,m,iccp1)=0.0
        
        ctem_mo%litres_mo(i,m,iccp2,1:ignd)=0.0
        ctem_mo%soilcres_mo(i,m,iccp2,1:ignd)=0.0
        ctem_mo%humiftrsveg_mo(i,m,iccp2)=0.0

      end do ! nmtest
    end do ! nltest

  end subroutine resetMonthEnd
  !! @}
  !==================================================

  !> \ingroup ctemstatevars_resetYearEnd
  !! @{
  !> Resets annual variables in preparation for next year
  subroutine resetYearEnd (nltest, nmtest)

    use classicParams,   only : iccp2, icc, iccp1, ignd

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i, m, j

    do i = 1,nltest
      ! Grid avg
      ctem_grd_yr%laimaxg_yr_g(i) = 0.0
      ctem_grd_yr%gleafmas_yr_g(i) = 0.0
      ctem_grd_yr%gleafmas_NS_yr_g(i) = 0.0
      ctem_grd_yr%gleafmass_yr_g(i) = 0.0
      ctem_grd_yr%bleafmas_yr_g(i) = 0.0
      ctem_grd_yr%stemmass_yr_g(i) = 0.0
      ctem_grd_yr%stemmass_NS_yr_g(i) = 0.0
      ctem_grd_yr%stemmasss_yr_g(i) = 0.0
      ctem_grd_yr%rootmass_yr_g(i) = 0.0
      ctem_grd_yr%rootmass_NS_yr_g(i) = 0.0
      ctem_grd_yr%rootmasss_yr_g(i) = 0.0
      ctem_grd_yr%litrmass_yr_g(i,1:ignd)=0.0
      ctem_grd_yr%soilcmas_yr_g(i,1:ignd)=0.0
      ctem_grd_yr%litres_yr_g(i,1:ignd)=0.0
      ctem_grd_yr%soilcres_yr_g(i,1:ignd)=0.0
      ctem_grd_yr%nh4_mass_yr_g(i) = 0.0
      ctem_grd_yr%no3_mass_yr_g(i) = 0.0
      ctem_grd_yr%ngleafmas_yr_g(i) = 0.0
      ctem_grd_yr%ngleafmas_NS_yr_g(i) = 0.0
      ctem_grd_yr%ngleafmass_yr_g(i) = 0.0
      ctem_grd_yr%nbleafmas_yr_g(i) = 0.0
      ctem_grd_yr%nstemmass_yr_g(i) = 0.0
      ctem_grd_yr%nstemmass_NS_yr_g(i) = 0.0
      ctem_grd_yr%nstemmasss_yr_g(i) = 0.0
      ctem_grd_yr%nrootmass_yr_g(i) = 0.0
      ctem_grd_yr%nrootmass_NS_yr_g(i) = 0.0
      ctem_grd_yr%nrootmasss_yr_g(i) = 0.0
      ctem_grd_yr%nvgbiomas_yr_g(i) = 0.0
      ctem_grd_yr%ndemand_wp_npp_yr_g(i) = 0.0
      ctem_grd_yr%c2n_l_yr_g(i) = 0.0
      ctem_grd_yr%c2n_s_yr_g(i) = 0.0
      ctem_grd_yr%c2n_r_yr_g(i) = 0.0
      ctem_grd_yr%c2n_wp_yr_g(i) = 0.0
      ctem_grd_yr%c2n_litr_yr_g(i) = 0.0
      ctem_grd_yr%c2n_humus_yr_g(i) = 0.0
      ctem_grd_yr%nlitrmass_yr_g(i) = 0.0
      ctem_grd_yr%soilnmas_yr_g(i) = 0.0
      ctem_grd_yr%vgbiomas_yr_g(i) = 0.0
      ctem_grd_yr%totcmass_yr_g(i) = 0.0
      ctem_grd_yr%veghght_yr_g(i) = 0.0
      ctem_grd_yr%npp_yr_g(i) = 0.0
      ctem_grd_yr%gpp_yr_g(i) = 0.0
      ctem_grd_yr%vcmax0_yr_g(i) = 0.0
      ctem_grd_yr%leafns2s_yr_g(i) = 0.0
      ctem_grd_yr%stemns2s_yr_g(i) = 0.0
      ctem_grd_yr%rootns2s_yr_g(i) = 0.0
      ctem_grd_yr%re_alloc_s2l_yr_g(i) = 0.0
      ctem_grd_yr%re_alloc_r2l_yr_g(i) = 0.0
      ctem_grd_yr%re_alloc_sr2l_yr_g(i)= 0.0
      ctem_grd_yr%nep_yr_g(i) = 0.0
      ctem_grd_yr%nepCMIP_yr_g(i) = 0.0
      ctem_grd_yr%nbp_yr_g(i) = 0.0
      ctem_grd_yr%hetrores_yr_g(i) = 0.0
      ctem_grd_yr%autores_yr_g(i) = 0.0
      ctem_grd_yr%rmrveg_yr_g(i) = 0.0
      ctem_grd_yr%litrfall_yr_g(i) = 0.0
      ctem_grd_yr%emit_co2_yr_g(i) = 0.0
      ctem_grd_yr%emit_co_yr_g(i) = 0.0
      ctem_grd_yr%emit_ch4_yr_g(i) = 0.0
      ctem_grd_yr%emit_nmhc_yr_g(i) = 0.0
      ctem_grd_yr%emit_h2_yr_g(i) = 0.0
      ctem_grd_yr%emit_nox_yr_g(i) = 0.0
      ctem_grd_yr%emit_n2o_yr_g(i) = 0.0
      ctem_grd_yr%emit_nh3_yr_g(i) = 0.0
      ctem_grd_yr%emit_pm25_yr_g(i) = 0.0
      ctem_grd_yr%emit_tpm_yr_g(i) = 0.0
      ctem_grd_yr%emit_tc_yr_g(i) = 0.0
      ctem_grd_yr%emit_oc_yr_g(i) = 0.0
      ctem_grd_yr%emit_bc_yr_g(i) = 0.0
      ctem_grd_yr%smfuncveg_yr_g(i) = 0.0
      ctem_grd_yr%luc_emc_yr_g(i) = 0.0
      ctem_grd_yr%lucsocin_yr_g(i) = 0.0
      ctem_grd_yr%lucltrin_yr_g(i) = 0.0
      ctem_grd_yr%luc_emcn_yr_g(i) = 0.0
      ctem_grd_yr%lucsocinn_yr_g(i) = 0.0
      ctem_grd_yr%lucltrinn_yr_g(i) = 0.0
      ctem_grd_yr%burnfrac_yr_g(i) = 0.0
      ctem_grd_yr%timharvarea_yr_g(i) = 0.0  
      ctem_grd_yr%tileAge_yr_g(i) = 0.0   
      ctem_grd_yr%bterm_yr_g(i) = 0.0
      ctem_grd_yr%lterm_yr_g(i) = 0.0
      ctem_grd_yr%mterm_yr_g(i) = 0.0
      ctem_grd_yr%ch4WetSpec_yr_g(i)  = 0.0
      ctem_grd_yr%wetfdyn_yr_g(i)  = 0.0
      ctem_grd_yr%ch4WetDyn_yr_g(i)  = 0.0
      ctem_grd_yr%ch4soills_yr_g(i)  = 0.0
      ctem_grd_yr%peatdep_yr_g(i)  = 0.0
      ctem_grd_yr%peatSoilC_yr_g(i)  = 0.0
      ctem_grd_yr%cProduct_yr_g(i)  = 0.0
      ctem_grd_yr%nProduct_yr_g(i)  = 0.0
      ctem_grd_yr%fProductDecomp_yr_g(i)  = 0.0
      ctem_grd_yr%bnf_tot_yr_g(i) = 0.0
      ctem_grd_yr%bnf_free_yr_g(i) = 0.0
      ctem_grd_yr%bnf_ant_yr_g(i) = 0.0
      ctem_grd_yr%bnf_nat_yr_g(i) = 0.0
      ctem_grd_yr%nstress_yr_g(i) = 0.0
      ctem_grd_yr%nitrif_yr_g(i) = 0.0
      ctem_grd_yr%no_nit_yr_g(i) = 0.0
      ctem_grd_yr%no_denit_yr_g(i) = 0.0
      ctem_grd_yr%no_nitdenit_yr_g(i) = 0.0
      ctem_grd_yr%n2o_nit_yr_g(i) = 0.0
      ctem_grd_yr%n2o_denit_yr_g(i) = 0.0
      ctem_grd_yr%n2o_nitdenit_yr_g(i) = 0.0
      ctem_grd_yr%n2_denit_yr_g(i) = 0.0
      ctem_grd_yr%nvol_yr_g(i) = 0.0
      ctem_grd_yr%nleach_yr_g(i) = 0.0
      ctem_grd_yr%appl_fert_yr_g(i) = 0.0
      ctem_grd_yr%ndep_nh4_yr_g(i) = 0.0
      ctem_grd_yr%ndep_no3_yr_g(i) = 0.0
      ctem_grd_yr%nuptake_p_nh4_yr_g(i) = 0.0
      ctem_grd_yr%nuptake_p_no3_yr_g(i) = 0.0
      ctem_grd_yr%nuptake_a_actl_nh4_yr_g(i) = 0.0
      ctem_grd_yr%nuptake_a_actl_no3_yr_g(i) = 0.0
      ctem_grd_yr%nuptake_yr_g(i) = 0.0
      ctem_grd_yr%nalloc_l_yr_g(i) = 0.0
      ctem_grd_yr%nalloc_s_yr_g(i) = 0.0
      ctem_grd_yr%nalloc_r_yr_g(i) = 0.0
      ctem_grd_yr%nresorped_s_yr_g(i) = 0.0
      ctem_grd_yr%nresorped_r_yr_g(i) = 0.0
      ctem_grd_yr%nre_alloc_s2l_yr_g(i) = 0.0
      ctem_grd_yr%nre_alloc_r2l_yr_g(i) = 0.0
      ctem_grd_yr%nleafns2s_yr_g(i) = 0.0
      ctem_grd_yr%nstemns2s_yr_g(i) = 0.0
      ctem_grd_yr%nrootns2s_yr_g(i) = 0.0
      ctem_grd_yr%nlitr_l_yr_g(i) = 0.0
      ctem_grd_yr%nlitr_s_yr_g(i) = 0.0
      ctem_grd_yr%nlitr_r_yr_g(i) = 0.0
      ctem_grd_yr%nlitr_yr_g(i) = 0.0
      ctem_grd_yr%gl2bl_grass_nflux_yr_g(i) = 0.0
      ctem_grd_yr%nhumtrs_yr_g(i) = 0.0
      ctem_grd_yr%nmineral_litr_yr_g(i) = 0.0
      ctem_grd_yr%nmineral_humus_yr_g(i) = 0.0
      ctem_grd_yr%netnmineral_yr_g(i) = 0.0
      ctem_grd_yr%nimmobil_nh4_yr_g(i) = 0.0
      ctem_grd_yr%nimmobil_no3_yr_g(i) = 0.0
      ctem_grd_yr%fNnetland_yr_g(i) = 0.0

      do m = 1,nmtest
        ! Tile avg
        ctem_tile_yr%laimaxg_yr_t(i,m) = 0.0
        ctem_tile_yr%gleafmas_yr_t(i,m) = 0.0
        ctem_tile_yr%gleafmas_NS_yr_t(i,m) = 0.0
        ctem_tile_yr%gleafmass_yr_t(i,m) = 0.0
        ctem_tile_yr%bleafmas_yr_t(i,m) = 0.0
        ctem_tile_yr%stemmass_yr_t(i,m) = 0.0
        ctem_tile_yr%stemmass_NS_yr_t(i,m) = 0.0
        ctem_tile_yr%stemmasss_yr_t(i,m) = 0.0
        ctem_tile_yr%rootmass_yr_t(i,m) = 0.0
        ctem_tile_yr%rootmass_NS_yr_t(i,m) = 0.0
        ctem_tile_yr%rootmasss_yr_t(i,m) = 0.0
        ctem_tile_yr%litrmass_yr_t(i,m,1:ignd)=0.0
        ctem_tile_yr%soilcmas_yr_t(i,m,1:ignd)=0.0
        ctem_tile_yr%litres_yr_t(i,m,1:ignd)=0.0
        ctem_tile_yr%soilcres_yr_t(i,m,1:ignd)=0.0
        ctem_tile_yr%nh4_mass_yr_t(i,m) = 0.0
        ctem_tile_yr%no3_mass_yr_t(i,m) = 0.0
        ctem_tile_yr%ngleafmas_yr_t(i,m) = 0.0
        ctem_tile_yr%ngleafmas_NS_yr_t(i,m) = 0.0
        ctem_tile_yr%ngleafmass_yr_t(i,m) = 0.0
        ctem_tile_yr%nbleafmas_yr_t(i,m)  = 0.0
        ctem_tile_yr%nstemmass_yr_t(i,m) = 0.0
        ctem_tile_yr%nstemmass_NS_yr_t(i,m)= 0.0
        ctem_tile_yr%nstemmasss_yr_t(i,m) = 0.0
        ctem_tile_yr%nrootmass_yr_t(i,m) = 0.0
        ctem_tile_yr%nrootmass_NS_yr_t(i,m)= 0.0
        ctem_tile_yr%nrootmasss_yr_t(i,m) = 0.0
        ctem_tile_yr%nvgbiomas_yr_t(i,m) = 0.0
        ctem_tile_yr%ndemand_wp_npp_yr_t(i,m) = 0.0
        ctem_tile_yr%c2n_l_yr_t(i,m) = 0.0
        ctem_tile_yr%c2n_s_yr_t(i,m) = 0.0
        ctem_tile_yr%c2n_r_yr_t(i,m) = 0.0
        ctem_tile_yr%c2n_wp_yr_t(i,m) = 0.0
        ctem_tile_yr%c2n_litr_yr_t(i,m) = 0.0
        ctem_tile_yr%c2n_humus_yr_t(i,m) = 0.0
        ctem_tile_yr%nlitrmass_yr_t(i,m) = 0.0
        ctem_tile_yr%soilnmas_yr_t(i,m) = 0.0
        ctem_tile_yr%vgbiomas_yr_t(i,m) = 0.0
        ctem_tile_yr%totcmass_yr_t(i,m) = 0.0
        ctem_tile_yr%veghght_yr_t(i,m) = 0.0
        ctem_tile_yr%npp_yr_t(i,m) = 0.0
        ctem_tile_yr%gpp_yr_t(i,m) = 0.0
        ctem_tile_yr%vcmax0_yr_t(i,m) = 0.0
        ctem_tile_yr%leafns2s_yr_t(i,m) = 0.0
        ctem_tile_yr%stemns2s_yr_t(i,m) = 0.0
        ctem_tile_yr%rootns2s_yr_t(i,m) = 0.0
        ctem_tile_yr%re_alloc_s2l_yr_t(i,m) = 0.0
        ctem_tile_yr%re_alloc_r2l_yr_t(i,m) = 0.0
        ctem_tile_yr%re_alloc_sr2l_yr_t(i,m) = 0.0
        ctem_tile_yr%nep_yr_t(i,m) = 0.0
        ctem_tile_yr%nepCMIP_yr_t(i,m) = 0.0
        ctem_tile_yr%nbp_yr_t(i,m) = 0.0
        ctem_tile_yr%hetrores_yr_t(i,m) = 0.0
        ctem_tile_yr%autores_yr_t(i,m) = 0.0
        ctem_tile_yr%rmrveg_yr_t(i,m) = 0.0
        ctem_tile_yr%litrfall_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_co2_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_co_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_ch4_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_nmhc_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_h2_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_nox_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_n2o_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_nh3_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_pm25_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_tpm_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_tc_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_oc_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_bc_yr_t(i,m) = 0.0
        ctem_tile_yr%smfuncveg_yr_t(i,m) = 0.0
        ctem_tile_yr%luc_emc_yr_t(i,m) = 0.0
        ctem_tile_yr%lucsocin_yr_t(i,m) = 0.0
        ctem_tile_yr%lucltrin_yr_t(i,m) = 0.0
        ctem_tile_yr%luc_emcn_yr_t(i,m) = 0.0
        ctem_tile_yr%lucsocinn_yr_t(i,m) = 0.0
        ctem_tile_yr%lucltrinn_yr_t(i,m) = 0.0
        ctem_tile_yr%burnfrac_yr_t(i,m) = 0.0
        ctem_tile_yr%timharvarea_yr_t(i,m) = 0.0  
        ctem_tile_yr%tileAge_yr_t(i,m) = 0.0   
        ctem_tile_yr%bterm_yr_t(i,m) = 0.0
        ctem_tile_yr%lterm_yr_t(i,m) = 0.0
        ctem_tile_yr%mterm_yr_t(i,m) = 0.0
        ctem_tile_yr%ch4WetSpec_yr_t(i,m)  = 0.0
        ctem_tile_yr%wetfdyn_yr_t(i,m)  = 0.0
        ctem_tile_yr%ch4WetDyn_yr_t(i,m)  = 0.0
        ctem_tile_yr%ch4soills_yr_t(i,m)  = 0.0
        ctem_tile_yr%peatdep_yr_t(i,m)  = 0.0
        ctem_tile_yr%peatSoilC_yr_t(i,m)  = 0.0
        ctem_tile_yr%fProductDecomp_yr_t(i,m) = 0.0
        ctem_tile_yr%bnf_tot_yr_t(i,m) = 0.0
        ctem_tile_yr%bnf_free_yr_t(i,m) = 0.0
        ctem_tile_yr%bnf_ant_yr_t(i,m) = 0.0
        ctem_tile_yr%bnf_nat_yr_t(i,m) = 0.0
        ctem_tile_yr%nstress_yr_t(i,m) = 0.0
        ctem_tile_yr%nitrif_yr_t(i,m) = 0.0
        ctem_tile_yr%no_nit_yr_t(i,m) = 0.0
        ctem_tile_yr%no_denit_yr_t(i,m) = 0.0
        ctem_tile_yr%no_nitdenit_yr_t(i,m) = 0.0
        ctem_tile_yr%n2o_nit_yr_t(i,m) = 0.0
        ctem_tile_yr%n2o_denit_yr_t(i,m) = 0.0
        ctem_tile_yr%n2o_nitdenit_yr_t(i,m) = 0.0
        ctem_tile_yr%n2_denit_yr_t(i,m) = 0.0
        ctem_tile_yr%nvol_yr_t(i,m) = 0.0
        ctem_tile_yr%nleach_yr_t(i,m) = 0.0
        ctem_tile_yr%appl_fert_yr_t(i,m) = 0.0
        ctem_tile_yr%ndep_nh4_yr_t(i,m) = 0.0
        ctem_tile_yr%ndep_no3_yr_t(i,m) = 0.0
        ctem_tile_yr%nuptake_p_nh4_yr_t(i,m) = 0.0
        ctem_tile_yr%nuptake_p_no3_yr_t(i,m) = 0.0
        ctem_tile_yr%nuptake_a_actl_nh4_yr_t(i,m) = 0.0
        ctem_tile_yr%nuptake_a_actl_no3_yr_t(i,m) = 0.0
        ctem_tile_yr%nuptake_yr_t(i,m) = 0.0
        ctem_tile_yr%nalloc_l_yr_t(i,m) = 0.0
        ctem_tile_yr%nalloc_s_yr_t(i,m) = 0.0
        ctem_tile_yr%nalloc_r_yr_t(i,m) = 0.0
        ctem_tile_yr%nresorped_s_yr_t(i,m) = 0.0
        ctem_tile_yr%nresorped_r_yr_t(i,m) = 0.0
        ctem_tile_yr%nre_alloc_s2l_yr_t(i,m) = 0.0
        ctem_tile_yr%nre_alloc_r2l_yr_t(i,m) = 0.0
        ctem_tile_yr%nleafns2s_yr_t(i,m) = 0.0
        ctem_tile_yr%nstemns2s_yr_t(i,m) = 0.0
        ctem_tile_yr%nrootns2s_yr_t(i,m) = 0.0
        ctem_tile_yr%nlitr_l_yr_t(i,m) = 0.0
        ctem_tile_yr%nlitr_s_yr_t(i,m) = 0.0
        ctem_tile_yr%nlitr_r_yr_t(i,m) = 0.0
        ctem_tile_yr%nlitr_yr_t(i,m) = 0.0
        ctem_tile_yr%gl2bl_grass_nflux_yr_t(i,m) = 0.0
        ctem_tile_yr%nhumtrs_yr_t(i,m) = 0.0
        ctem_tile_yr%nmineral_litr_yr_t(i,m) = 0.0
        ctem_tile_yr%nmineral_humus_yr_t(i,m) = 0.0
        ctem_tile_yr%netnmineral_yr_t(i,m) = 0.0
        ctem_tile_yr%nimmobil_nh4_yr_t(i,m) = 0.0
        ctem_tile_yr%nimmobil_no3_yr_t(i,m) = 0.0
        ctem_tile_yr%fNnetland_yr_t(i,m) = 0.0
        ctem_yr%timharvarea_yr(i,m) = 0.0  
        ctem_yr%tileAge_yr(i,m) = 0.0   

        do j = 1,icc
          ! per pft
          ctem_yr%laimaxg_yr(i,m,j) = 0.0
          ctem_yr%gleafmas_yr(i,m,j) = 0.0
          ctem_yr%gleafmas_NS_yr(i,m,j) = 0.0
          ctem_yr%gleafmass_yr(i,m,j) = 0.0
          ctem_yr%bleafmas_yr(i,m,j) = 0.0
          ctem_yr%stemmass_yr(i,m,j) = 0.0
          ctem_yr%stemmass_NS_yr(i,m,j) = 0.0
          ctem_yr%stemmasss_yr(i,m,j) = 0.0
          ctem_yr%rootmass_yr(i,m,j) = 0.0
          ctem_yr%rootmass_NS_yr(i,m,j) = 0.0
          ctem_yr%rootmasss_yr(i,m,j) = 0.0
          ctem_yr%litrmass_yr(i,m,j,1:ignd)=0.0
          ctem_yr%soilcmas_yr(i,m,j,1:ignd)=0.0
          ctem_yr%litres_yr(i,m,j,1:ignd)=0.0
          ctem_yr%soilcres_yr(i,m,j,1:ignd)=0.0
          ctem_yr%nh4_mass_yr(i,m,j) = 0.0
          ctem_yr%no3_mass_yr(i,m,j) = 0.0
          ctem_yr%ngleafmas_yr(i,m,j) = 0.0
          ctem_yr%ngleafmas_NS_yr(i,m,j) = 0.0
          ctem_yr%ngleafmass_yr(i,m,j) = 0.0
          ctem_yr%nbleafmas_yr(i,m,j) = 0.0
          ctem_yr%nstemmass_yr(i,m,j) = 0.0
          ctem_yr%nstemmass_NS_yr(i,m,j) = 0.0
          ctem_yr%nstemmasss_yr(i,m,j) = 0.0
          ctem_yr%nrootmass_yr(i,m,j) = 0.0
          ctem_yr%nrootmass_NS_yr(i,m,j) = 0.0
          ctem_yr%nrootmasss_yr(i,m,j) = 0.0
          ctem_yr%nvgbiomas_yr(i,m,j) = 0.0
          ctem_yr%ndemand_wp_npp_yr(i,m,j) = 0.0
          ctem_yr%c2n_l_yr(i,m,j) = 0.0
          ctem_yr%c2n_s_yr(i,m,j) = 0.0
          ctem_yr%c2n_r_yr(i,m,j) = 0.0
          ctem_yr%c2n_wp_yr(i,m,j) = 0.0
          ctem_yr%c2n_litr_yr(i,m,j) = 0.0
          ctem_yr%c2n_humus_yr(i,m,j) = 0.0
          ctem_yr%nlitrmass_yr(i,m,j) = 0.0
          ctem_yr%soilnmas_yr(i,m,j) = 0.0
          ctem_yr%vgbiomas_yr(i,m,j) = 0.0
          ctem_yr%totcmass_yr(i,m,j) = 0.0
          ctem_yr%veghght_yr(i,m,j) = 0.0
          ctem_yr%cc_yr(i,m,j) = 0.0
          ctem_yr%mm_yr(i,m,j) = 0.0
          ctem_yr%npp_yr(i,m,j) = 0.0
          ctem_yr%gpp_yr(i,m,j) = 0.0
          ctem_yr%vcmax0_yr(i,m,j) = 0.0
          ctem_yr%leafns2s_yr(i,m,j) = 0.0
          ctem_yr%stemns2s_yr(i,m,j) = 0.0
          ctem_yr%rootns2s_yr(i,m,j) = 0.0
          ctem_yr%re_alloc_s2l_yr(i,m,j) = 0.0
          ctem_yr%re_alloc_r2l_yr(i,m,j) = 0.0
          ctem_yr%re_alloc_sr2l_yr(i,m,j) = 0.0
          ctem_yr%nep_yr(i,m,j) = 0.0
          ctem_yr%nepCMIP_yr(i,m,j) = 0.0
          ctem_yr%nbp_yr(i,m,j) = 0.0
          ctem_yr%hetrores_yr(i,m,j) = 0.0
          ctem_yr%autores_yr(i,m,j) = 0.0
          ctem_yr%rmrveg_yr(i,m,j) = 0.0
          ctem_yr%litrfall_yr(i,m,j) = 0.0
          ctem_yr%emit_co2_yr(i,m,j) = 0.0
          ctem_yr%emit_co_yr(i,m,j) = 0.0
          ctem_yr%emit_ch4_yr(i,m,j) = 0.0
          ctem_yr%emit_nmhc_yr(i,m,j) = 0.0
          ctem_yr%emit_h2_yr(i,m,j) = 0.0
          ctem_yr%emit_nox_yr(i,m,j) = 0.0
          ctem_yr%emit_n2o_yr(i,m,j) = 0.0
          ctem_yr%emit_nh3_yr(i,m,j) = 0.0
          ctem_yr%emit_pm25_yr(i,m,j) = 0.0
          ctem_yr%emit_tpm_yr(i,m,j) = 0.0
          ctem_yr%emit_tc_yr(i,m,j) = 0.0
          ctem_yr%emit_oc_yr(i,m,j) = 0.0
          ctem_yr%emit_bc_yr(i,m,j) = 0.0
          ctem_yr%bterm_yr(i,m,j) = 0.0
          ctem_yr%mterm_yr(i,m,j) = 0.0
          ctem_yr%burnfrac_yr(i,m,j) = 0.0
          ctem_yr%smfuncveg_yr(i,m,j) = 0.0
          ctem_yr%bnf_tot_yr(i,m,j) = 0.0
          ctem_yr%bnf_free_yr(i,m,j) = 0.0
          ctem_yr%bnf_ant_yr(i,m,j) = 0.0
          ctem_yr%bnf_nat_yr(i,m,j) = 0.0
          ctem_yr%nstress_yr(i,m,j) = 0.0
          ctem_yr%nitrif_yr(i,m,j) = 0.0
          ctem_yr%no_nit_yr(i,m,j) = 0.0
          ctem_yr%no_denit_yr(i,m,j) = 0.0
          ctem_yr%no_nitdenit_yr(i,m,j) = 0.0
          ctem_yr%n2o_nit_yr(i,m,j) = 0.0
          ctem_yr%n2o_denit_yr(i,m,j) = 0.0
          ctem_yr%n2o_nitdenit_yr(i,m,j) = 0.0
          ctem_yr%n2_denit_yr(i,m,j) = 0.0
          ctem_yr%nvol_yr(i,m,j) = 0.0
          ctem_yr%nleach_yr(i,m,j) = 0.0
          ctem_yr%appl_fert_yr(i,m,j) = 0.0
          ctem_yr%ndep_nh4_yr(i,m,j) = 0.0
          ctem_yr%ndep_no3_yr(i,m,j) = 0.0
          ctem_yr%nuptake_p_nh4_yr(i,m,j) = 0.0
          ctem_yr%nuptake_p_no3_yr(i,m,j) = 0.0
          ctem_yr%nuptake_a_actl_nh4_yr(i,m,j) = 0.0
          ctem_yr%nuptake_a_actl_no3_yr(i,m,j) = 0.0
          ctem_yr%nuptake_yr(i,m,j) = 0.0
          ctem_yr%nalloc_l_yr(i,m,j) = 0.0
          ctem_yr%nalloc_s_yr(i,m,j) = 0.0
          ctem_yr%nalloc_r_yr(i,m,j) = 0.0
          ctem_yr%nresorped_s_yr(i,m,j) = 0.0
          ctem_yr%nresorped_r_yr(i,m,j) = 0.0
          ctem_yr%nre_alloc_s2l_yr(i,m,j) = 0.0
          ctem_yr%nre_alloc_r2l_yr(i,m,j) = 0.0
          ctem_yr%nleafns2s_yr(i,m,j) = 0.0
          ctem_yr%nstemns2s_yr(i,m,j) = 0.0
          ctem_yr%nrootns2s_yr(i,m,j) = 0.0
          ctem_yr%nlitr_l_yr(i,m,j) = 0.0
          ctem_yr%nlitr_s_yr(i,m,j) = 0.0
          ctem_yr%nlitr_r_yr(i,m,j) = 0.0
          ctem_yr%nlitr_yr(i,m,j) = 0.0
          ctem_yr%gl2bl_grass_nflux_yr(i,m,j) = 0.0
          ctem_yr%nhumtrs_yr(i,m,j) = 0.0
          ctem_yr%nmineral_litr_yr(i,m,j) = 0.0
          ctem_yr%nmineral_humus_yr(i,m,j) = 0.0
          ctem_yr%netnmineral_yr(i,m,j) = 0.0
          ctem_yr%nimmobil_nh4_yr(i,m,j) = 0.0
          ctem_yr%nimmobil_no3_yr(i,m,j) = 0.0
          ctem_yr%fNnetland_yr(i,m,j) = 0.0
        end do

        ctem_yr%hetrores_yr(i,m,iccp1) = 0.0
        ctem_yr%nep_yr(i,m,iccp1) = 0.0
        ctem_yr%nepCMIP_yr(i,m,iccp1) = 0.0
        ctem_yr%nbp_yr(i,m,iccp1) = 0.0
        ctem_yr%totcmass_yr(i,m,iccp1) = 0.0

        ctem_yr%litres_yr(i,m,iccp1,1:ignd)=0.0
        ctem_yr%soilcres_yr(i,m,iccp1,1:ignd)=0.0
        ctem_yr%litrmass_yr(i,m,iccp1,1:ignd)=0.0
        ctem_yr%soilcmas_yr(i,m,iccp1,1:ignd)=0.0

        ctem_yr%litrmass_yr(i,m,iccp2,1:ignd)=0.0
        ctem_yr%soilcmas_yr(i,m,iccp2,1:ignd)=0.0
        ctem_yr%litres_yr(i,m,iccp2,1:ignd)=0.0
        ctem_yr%soilcres_yr(i,m,iccp2,1:ignd)=0.0
        ctem_yr%bnf_tot_yr(i,m,iccp1) = 0.0
        ctem_yr%bnf_free_yr(i,m,iccp1) = 0.0
        ctem_yr%nitrif_yr(i,m,iccp1) = 0.0
        ctem_yr%no_nit_yr(i,m,iccp1) = 0.0
        ctem_yr%no_denit_yr(i,m,iccp1) = 0.0
        ctem_yr%no_nitdenit_yr(i,m,iccp1) = 0.0
        ctem_yr%n2o_nit_yr(i,m,iccp1) = 0.0
        ctem_yr%n2o_denit_yr(i,m,iccp1) = 0.0
        ctem_yr%n2o_nitdenit_yr(i,m,iccp1) = 0.0
        ctem_yr%n2_denit_yr(i,m,iccp1) = 0.0
        ctem_yr%nvol_yr(i,m,iccp1) = 0.0
        ctem_yr%nleach_yr(i,m,iccp1) = 0.0
        ctem_yr%nh4_mass_yr(i,m,iccp1) = 0.0
        ctem_yr%no3_mass_yr(i,m,iccp1) = 0.0
        ctem_yr%nlitrmass_yr(i,m,iccp1) = 0.0
        ctem_yr%soilnmas_yr(i,m,iccp1) = 0.0
        ctem_yr%nlitrmass_yr(i,m,iccp2) = 0.0
        ctem_yr%soilnmas_yr(i,m,iccp2) = 0.0
        ctem_yr%appl_fert_yr(i,m,iccp1) = 0.0
        ctem_yr%ndep_nh4_yr(i,m,iccp1) = 0.0
        ctem_yr%ndep_no3_yr(i,m,iccp1) = 0.0
        ctem_yr%c2n_litr_yr(i,m,iccp1) = 0.0
        ctem_yr%c2n_humus_yr(i,m,iccp1) = 0.0
        ctem_yr%nhumtrs_yr(i,m,iccp1) = 0.0
        ctem_yr%nmineral_litr_yr(i,m,iccp1) = 0.0
        ctem_yr%nmineral_humus_yr(i,m,iccp1) = 0.0
        ctem_yr%netnmineral_yr(i,m,iccp1) = 0.0
        ctem_yr%nimmobil_nh4_yr(i,m,iccp1) = 0.0
        ctem_yr%nimmobil_no3_yr(i,m,iccp1) = 0.0
        ctem_yr%fNnetland_yr(i,m,iccp1) = 0.0

      end do ! nmtest
    end do ! nltest

  end subroutine resetYearEnd
  !! @}
  !==================================================
  !> \ingroup ctemstatevars_resetMosaicAccum
  !! @{
  !> Resets physics accumulator variables (used as input to CTEM) after CTEM has been called
  subroutine resetMosaicAccum

    implicit none

    vgat%fsinacc_gat(:) = 0.
    vgat%flinacc_gat(:) = 0.
    vgat%flutacc_gat(:) = 0.
    vgat%preacc_gat(:) = 0.
    ctem_tile%fsnowacc_t(:) = 0.0
    ctem_tile%taaccgat_t(:) = 0.0
    ctem_tile%vvaccgat_t(:) = 0.0
    ctem_tile%uvaccgat_t(:) = 0.0
    ctem_tile%CFLUX_GAacc_t(:) = 0.0
    ctem_tile%USTARBS_GAacc_t(:) = 0.0
    ctem_tile%ROFBacc_t(:) = 0.0
    vgat%altotacc_gat(:) = 0.0
    vgat%altotcount_ctm(:) = 0
    ctem_tile%QFCacc_t(:,:) = 0.0
    ctem_tile%tbaraccgat_t(:,:) = 0.0
    ctem_tile%thliqacc_t(:,:) = 0.0
    ctem_tile%thiceacc_t(:,:) = 0.0

    ctem_tile%ancgvgac_t(:,:) = 0.0
    ctem_tile%rmlcgvga_t(:,:) = 0.0

    !-reset peatland accumulators-------------------------------
    ctem_tile%anmossac_t(:)  = 0.0
    ctem_tile%rmlmossac_t(:) = 0.0
    ctem_tile%gppmossac_t(:) = 0.0

  end subroutine resetMosaicAccum
  !! @}
  !=================================================================================

  !==================================================
  !> \ingroup ctemstatevars_ctemdump
  !! @{
  !> dumps all the ctemdump statevars for diagnostic purposes.
  subroutine ctemdump

    implicit none

    print *, 'begin ctem dump'
    print *, 'agcm_classic_on ', c_switch%agcm_classic_on
    print *, 'projectedGrid ', c_switch%projectedGrid
    print *, 'ctem_on ', c_switch%ctem_on
    print *, 'Ncycle_on ', c_switch%Ncycle_on
    print *, 'metLoop ', c_switch%metLoop
    print *, 'leap ', c_switch%leap
    print *, 'spinfast ', c_switch%spinfast
    print *, 'useTracer ', c_switch%useTracer
    print *, 'tracerCO2file ', c_switch%tracerCO2file
    print *, 'readMetStartYear ', c_switch%readMetStartYear
    print *, 'readMetEndYear ', c_switch%readMetEndYear
    print *, 'transientCO2 ', c_switch%transientCO2
    print *, 'CO2File ', c_switch%CO2File
    print *, 'fixedYearCO2 ', c_switch%fixedYearCO2
    print *, 'doMethane ', c_switch%doMethane
    print *, 'transientCH4 ', c_switch%transientCH4
    print *, 'CH4File ', c_switch%CH4File
    print *, 'fixedYearCH4 ', c_switch%fixedYearCH4
    print *, 'dofire ', c_switch%dofire
    print *, 'transientPOPD ', c_switch%transientPOPD
    print *, 'POPDFile ', c_switch%POPDFile
    print *, 'fixedYearPOPD ', c_switch%fixedYearPOPD
    print *, 'transientLGHT ', c_switch%transientLGHT
    print *, 'LGHTFile ', c_switch%LGHTFile
    print *, 'fixedYearLGHT ', c_switch%fixedYearLGHT
    print *, 'PFTCompetition ', c_switch%PFTCompetition
    print *, 'start_bare ', c_switch%start_bare
    print *, 'inibioclim ', c_switch%inibioclim
    print *, 'lnduseon ', c_switch%lnduseon
    print *, 'LUCFile ', c_switch%LUCFile
    print *, 'fixedYearLUC ', c_switch%fixedYearLUC
    print *, 'fertilizeron ', c_switch%fertilizeron
    print *, 'FERFile ', c_switch%FERFile
    print *, 'fixedYearFER ', c_switch%fixedYearFER
    print *, 'transientFER ', c_switch%transientFER
    print *, 'depositionon ', c_switch%depositionon
    print *, 'DEPFile ', c_switch%DEPFile
    print *, 'fixedYearDEP ', c_switch%fixedYearDEP
    print *, 'transientDEP ', c_switch%transientDEP
    print *, 'transientOBSWETF ', c_switch%transientOBSWETF
    print *, 'OBSWETFFile ', c_switch%OBSWETFFile
    print *, 'fixedYearOBSWETF ', c_switch%fixedYearOBSWETF
    print *, 'allLocalTime ', c_switch%allLocalTime
    print *, 'metFileFss ', c_switch%metFileFss
    print *, 'metFilefracFsf ', c_switch%metFilefracFsf
    print *, 'metFileFdl ', c_switch%metFileFdl
    print *, 'metFileSnow ', c_switch%metFileSnow
    print *, 'metFilePre ', c_switch%metFilePre
    print *, 'metFileTa ', c_switch%metFileTa
    print *, 'metFileQa ', c_switch%metFileQa
    print *, 'metFileUv ', c_switch%metFileUv
    print *, 'metFilePres ', c_switch%metFilePres
    print *, 'init_file ', c_switch%init_file
    print *, 'rs_file_to_overwrite ', c_switch%rs_file_to_overwrite
    print *, 'runparams_file ', c_switch%runparams_file
    print *, 'Comment ', c_switch%Comment
    print *, 'output_directory ', c_switch%output_directory
    print *, 'xmlFile ', c_switch%xmlFile
    print *, 'doperpftoutput ', c_switch%doperpftoutput
    print *, 'dopertileoutput ', c_switch%dopertileoutput
    print *, 'doChecksums ', c_switch%doChecksums
    print *, 'doAnnualOutput ', c_switch%doAnnualOutput
    print *, 'doMonthOutput ', c_switch%doMonthOutput
    print *, 'jmosty ', c_switch%jmosty
    print *, 'doDayOutput ', c_switch%doDayOutput
    print *, 'jdstd ', c_switch%jdstd
    print *, 'jdendd ', c_switch%jdendd
    print *, 'jdsty ', c_switch%jdsty
    print *, 'jdendy ', c_switch%jdendy
    print *, 'doHhOutput ', c_switch%doHhOutput
    print *, 'jhhstd ', c_switch%jhhstd
    print *, 'jhhendd ', c_switch%jhhendd
    print *, 'jhhsty ', c_switch%jhhsty
    print *, 'jhhendy ', c_switch%jhhendy
    print *, 'idisp ', c_switch%idisp
    print *, 'izref ', c_switch%izref
    print *, 'islfd ', c_switch%islfd
    print *, 'ipcp ', c_switch%ipcp
    print *, 'iwf ', c_switch%iwf
    print *, 'ITC ', c_switch%ITC
    print *, 'ITCG ', c_switch%ITCG
    print *, 'ITG ', c_switch%ITG
    print *, 'IPAI ', c_switch%IPAI
    print *, 'IHGT ', c_switch%IHGT
    print *, 'IALC ', c_switch%IALC
    print *, 'IALS ', c_switch%IALS
    print *, 'IALG ', c_switch%IALG
    print *, 'isnoalb ', c_switch%isnoalb
    print *, 'alb4BandParamsFile ', c_switch%alb4BandParamsFile
    print *, 'KsatScalaron ', c_switch%KsatScalaron

    print *, 'pftexist ', vrot%pftexist
    print *, 'lfstatus ', vrot%lfstatus
    print *, 'pandays ', vrot%pandays
    print *, 'gleafmas ', vrot%gleafmas
    print *, 'gleafmas_ns ', vrot%gleafmas_ns
    print *, 'gleafmas_s ', vrot%gleafmas_s
    print *, 'leafns2s ', vrot%leafns2s
    print *, 'stemns2s ', vrot%stemns2s
    print *, 'rootns2s ', vrot%rootns2s
    print *, 're_alloc_s2l ', vrot%re_alloc_s2l
    print *, 're_alloc_r2l ', vrot%re_alloc_r2l
    print *, 're_alloc_sr2l ', vrot%re_alloc_sr2l
    print *, 'bleafmas ', vrot%bleafmas
    print *, 'stemmass ', vrot%stemmass
    print *, 'stemmass_ns ', vrot%stemmass_ns
    print *, 'stemmass_s ', vrot%stemmass_s
    print *, 'rootmass ', vrot%rootmass
    print *, 'rootmass_ns ', vrot%rootmass_ns
    print *, 'rootmass_s ', vrot%rootmass_s
    print *, 'pstemmass ', vrot%pstemmass
    print *, 'pgleafmass ', vrot%pgleafmass
    print *, 'fcancmx ', vrot%fcancmx
    print *, 'ngleafmas ', vrot%ngleafmas
    print *, 'ngleafmas_ns ', vrot%ngleafmas_ns
    print *, 'ngleafmas_s ', vrot%ngleafmas_s
    print *, 'nbleafmas ', vrot%nbleafmas
    print *, 'nstemmass ', vrot%nstemmass
    print *, 'nstemmass_ns ', vrot%nstemmass_ns
    print *, 'nstemmass_s ', vrot%nstemmass_s
    print *, 'nrootmass ', vrot%nrootmass
    print *, 'nrootmass_ns ', vrot%nrootmass_ns
    print *, 'nrootmass_s ', vrot%nrootmass_s
    print *, 'ailcg ', vrot%ailcg
    print *, 'ailcgs ', vrot%ailcgs
    print *, 'fcancs ', vrot%fcancs
    print *, 'fcanc ', vrot%fcanc
    print *, 'co2i1cg ', vrot%co2i1cg
    print *, 'co2i1cs ', vrot%co2i1cs
    print *, 'co2i2cg ', vrot%co2i2cg
    print *, 'co2i2cs ', vrot%co2i2cs
    print *, 'ancsveg ', vrot%ancsveg
    print *, 'ancgveg ', vrot%ancgveg
    print *, 'rmlcsveg ', vrot%rmlcsveg
    print *, 'rmlcgveg ', vrot%rmlcgveg
    print *, 'slai ', vrot%slai
    print *, 'ailcb ', vrot%ailcb
    print *, 'flhrloss ', vrot%flhrloss
    print *, 'flhrloss_ns ', vrot%flhrloss_ns
    print *, 'flhrloss_s ', vrot%flhrloss_s
    print *, 'grwtheff ', vrot%grwtheff
    print *, 'lystmmas ', vrot%lystmmas
    print *, 'lyrotmas ', vrot%lyrotmas
    print *, 'tymaxlai ', vrot%tymaxlai
    print *, 'stmhrlos ', vrot%stmhrlos
    print *, 'vgbiomas_veg ', vrot%vgbiomas_veg
    print *, 'emit_co2 ', vrot%emit_co2
    print *, 'emit_co ', vrot%emit_co
    print *, 'emit_ch4 ', vrot%emit_ch4
    print *, 'emit_nmhc ', vrot%emit_nmhc
    print *, 'emit_h2 ', vrot%emit_h2
    print *, 'emit_nox ', vrot%emit_nox
    print *, 'emit_n2o ', vrot%emit_n2o
    print *, 'emit_pm25 ', vrot%emit_pm25
    print *, 'emit_tpm ', vrot%emit_tpm
    print *, 'emit_tc ', vrot%emit_tc
    print *, 'emit_oc ', vrot%emit_oc
    print *, 'emit_bc ', vrot%emit_bc
    print *, 'burnvegf ', vrot%burnvegf
    print *, 'smfuncveg ', vrot%smfuncveg
    print *, 'bterm ', vrot%bterm
    print *, 'mterm ', vrot%mterm
    print *, 'bmasveg ', vrot%bmasveg
    print *, 'veghght ', vrot%veghght
    print *, 'rootdpth ', vrot%rootdpth
    print *, 'tltrleaf ', vrot%tltrleaf
    print *, 'tltrstem ', vrot%tltrstem
    print *, 'tltrroot ', vrot%tltrroot
    print *, 'leaflitr ', vrot%leaflitr
    print *, 'roottemp ', vrot%roottemp
    print *, 'afrleaf ', vrot%afrleaf
    print *, 'afrstem ', vrot%afrstem
    print *, 'afrroot ', vrot%afrroot
    print *, 'wtstatus ', vrot%wtstatus
    print *, 'ltstatus ', vrot%ltstatus
    print *, 'gppveg ', vrot%gppveg
    print *, 'vcmax0 ', vrot%vcmax0
    print *, 'nppveg ', vrot%nppveg
    print *, 'autoresveg ', vrot%autoresveg
    print *, 'rmlvegacc ', vrot%rmlvegacc
    print *, 'rmsveg ', vrot%rmsveg
    print *, 'rmrveg ', vrot%rmrveg
    print *, 'rgveg ', vrot%rgveg
    print *, 'litrfallveg ', vrot%litrfallveg
    print *, 'rothrlos ', vrot%rothrlos
    print *, 'pfcancmx ', vrot%pfcancmx
    print *, 'nfcancmx ', vrot%nfcancmx
    print *, 'anveg ', vrot%anveg
    print *, 'rmlveg ', vrot%rmlveg
    print *, 'bnf_free ', vrot%bnf_free
    print *, 'bnf_ant ', vrot%bnf_ant
    print *, 'bnf_nat ', vrot%bnf_nat
    print *, 'bnf_tot ', vrot%bnf_tot
    print *, 'nstress ', vrot%nstress
    print *, 'nitrifveg ', vrot%nitrifveg
    print *, 'no_nitveg ', vrot%no_nitveg
    print *, 'no_denitveg ', vrot%no_denitveg
    print *, 'no_nitdenitveg ', vrot%no_nitdenitveg
    print *, 'n2o_nitveg ', vrot%n2o_nitveg
    print *, 'n2o_denitveg ', vrot%n2o_denitveg
    print *, 'n2o_nitdenitveg ', vrot%n2o_nitdenitveg
    print *, 'n2_denitveg ', vrot%n2_denitveg
    print *, 'nvolveg ', vrot%nvolveg
    print *, 'nleachveg ', vrot%nleachveg
    print *, 'appl_fert ', vrot%appl_fert
    print *, 'ndep_nh4 ', vrot%ndep_nh4
    print *, 'ndep_no3 ', vrot%ndep_no3
    print *, 'ndemandveg_wp_npp ', vrot%ndemandveg_wp_npp
    print *, 'nuptakeveg_p_nh4 ', vrot%nuptakeveg_p_nh4
    print *, 'nuptakeveg_p_no3 ', vrot%nuptakeveg_p_no3
    print *, 'nuptakeveg_a_actl_nh4 ', vrot%nuptakeveg_a_actl_nh4
    print *, 'nuptakeveg_a_actl_no3 ', vrot%nuptakeveg_a_actl_no3
    print *, 'nuptakeveg ', vrot%nuptakeveg
    print *, 'nleafns2sveg ', vrot%nleafns2sveg
    print *, 'nstemns2sveg ', vrot%nstemns2sveg
    print *, 'nrootns2sveg ', vrot%nrootns2sveg
    print *, 'nallocveg_l ', vrot%nallocveg_l
    print *, 'nallocveg_s ', vrot%nallocveg_s
    print *, 'nallocveg_r ', vrot%nallocveg_r
    print *, 'nresorpedveg_s ', vrot%nresorpedveg_s
    print *, 'nresorpedveg_r ', vrot%nresorpedveg_r
    print *, 'nre_allocveg_s2l ', vrot%nre_allocveg_s2l
    print *, 'nre_allocveg_r2l ', vrot%nre_allocveg_r2l
    print *, 'nlitrveg_l ', vrot%nlitrveg_l
    print *, 'nlitrveg_s ', vrot%nlitrveg_s
    print *, 'nlitrveg_r ', vrot%nlitrveg_r
    print *, 'nlitrveg ', vrot%nlitrveg
    print *, 'gl2bl_grass_nflux ', vrot%gl2bl_grass_nflux
    print *, 'c2nveg_l ', vrot%c2nveg_l
    print *, 'c2nveg_s ', vrot%c2nveg_s
    print *, 'c2nveg_r ', vrot%c2nveg_r
    print *, 'c2nveg_wp ', vrot%c2nveg_wp
    print *, 'c2nveg_litr ', vrot%c2nveg_litr
    print *, 'c2nveg_humus ', vrot%c2nveg_humus
    print *, 'nhumtrsveg ', vrot%nhumtrsveg
    print *, 'nmineralveg_litr ', vrot%nmineralveg_litr
    print *, 'nmineralveg_humus ', vrot%nmineralveg_humus
    print *, 'netnmineralveg ', vrot%netnmineralveg
    print *, 'nimmobilveg_nh4 ', vrot%nimmobilveg_nh4
    print *, 'nimmobilveg_no3 ', vrot%nimmobilveg_no3
    print *, 'nvgbiomas_veg ', vrot%nvgbiomas_veg
    print *, 'fNnetlandveg ', vrot%fNnetlandveg
    print *, 'redcoeff_vcmax ', vrot%redcoeff_vcmax
    print *, 'gavglai ', vrot%gavglai
    print *, 'co2conc ', vrot%co2conc
    print *, 'ch4conc ', vrot%ch4conc
    print *, 'canres ', vrot%canres
    print *, 'vgbiomas ', vrot%vgbiomas
    print *, 'gavgltms ', vrot%gavgltms
    print *, 'gavgscms ', vrot%gavgscms
    print *, 'burnfrac ', vrot%burnfrac
    print *, 'timharvarearot ', vrot%timharvarearow
    print *, 'tileAgerow ', vrot%tileAgerow
    print *, 'popdin ', vrot%popdin
    print *, 'soilpH ', vrot%soilpH
    print *, 'nfertil ', vrot%nfertil
    print *, 'ndeposit ', vrot%ndeposit
    print *, 'lterm ', vrot%lterm
    print *, 'extnprob ', vrot%extnprob
    print *, 'prbfrhuc ', vrot%prbfrhuc
    print *, 'rml ', vrot%rml
    print *, 'rms ', vrot%rms
    print *, 'rmr ', vrot%rmr
    print *, 'ch4WetSpec ', vrot%ch4WetSpec
    print *, 'wetfdyn ', vrot%wetfdyn
    print *, 'wetfrac_pres ', vrot%wetfrac_pres
    print *, 'ch4WetDyn ', vrot%ch4WetDyn
    print *, 'ch4_soills ', vrot%ch4_soills
    print *, 'lucemcom ', vrot%lucemcom
    print *, 'lucltrin ', vrot%lucltrin
    print *, 'lucsocin ', vrot%lucsocin
    print *, 'npp ', vrot%npp
    print *, 'nep ', vrot%nep
    print *, 'nbp ', vrot%nbp
    print *, 'gpp ', vrot%gpp
    print *, 'hetrores ', vrot%hetrores
    print *, 'autores ', vrot%autores
    print *, 'soilcresp ', vrot%soilcresp
    print *, 'rm ', vrot%rm
    print *, 'rg ', vrot%rg
    print *, 'litres ', vrot%litres
    print *, 'socres ', vrot%socres
    print *, 'dstcemls ', vrot%dstcemls
    print *, 'litrfall ', vrot%litrfall
    print *, 'humiftrs ', vrot%humiftrs
    print *, 'cfluxcg ', vrot%cfluxcg
    print *, 'cfluxcs ', vrot%cfluxcs
    print *, 'ROFB ', vrot%ROFB
    print *, 'dstcemls3 ', vrot%dstcemls3
    print *, 'uvaccrow_m ', vrot%uvaccrow_m
    print *, 'vvaccrow_m ', vrot%vvaccrow_m
    print *, 'qevpacc_m_save ', vrot%qevpacc_m_save
    print *, 'twarmm ', vrot%twarmm
    print *, 'tcoldm ', vrot%tcoldm
    print *, 'gdd5 ', vrot%gdd5
    print *, 'aridity ', vrot%aridity
    print *, 'srplsmon ', vrot%srplsmon
    print *, 'defctmon ', vrot%defctmon
    print *, 'anndefct ', vrot%anndefct
    print *, 'annsrpls ', vrot%annsrpls
    print *, 'annpcp ', vrot%annpcp
    print *, 'dry_season_length ', vrot%dry_season_length
    print *, 'ipeatland ', vrot%ipeatland
    print *, 'litrmsmoss ', vrot%litrmsmoss
    print *, 'Cmossmas ', vrot%Cmossmas
    print *, 'dmoss ', vrot%dmoss
    print *, 'peatSoilC ', vrot%peatSoilC
    print *, 'nppmoss ', vrot%nppmoss
    print *, 'rmlmoss ', vrot%rmlmoss
    print *, 'gppmoss ', vrot%gppmoss
    print *, 'anmoss ', vrot%anmoss
    print *, 'armoss ', vrot%armoss
    print *, 'peatdep ', vrot%peatdep
    print *, 'pdd ', vrot%pdd
    print *, 'colddays_leaffall ', vrot%colddays_leaffall
    print *, 'colddays_harvest ', vrot%colddays_harvest
    print *, 'zolnc ', vrot%zolnc
    print *, 'ailc ', vrot%ailc
    print *, 'cmasvegc ', vrot%cmasvegc
    print *, 'alvsctm ', vrot%alvsctm
    print *, 'paic ', vrot%paic
    print *, 'slaic ', vrot%slaic
    print *, 'alirctm ', vrot%alirctm
    print *, 'rmatc ', vrot%rmatc
    print *, 'rmatctem ', vrot%rmatctem
    print *, 'nepveg ', vrot%nepveg
    print *, 'nbpveg ', vrot%nbpveg
    print *, 'hetroresveg ', vrot%hetroresveg
    print *, 'litrmass ', vrot%litrmass
    print *, 'soilcmas ', vrot%soilcmas
    print *, 'litresveg ', vrot%litresveg
    print *, 'soilcresveg ', vrot%soilcresveg
    print *, 'humiftrsveg ', vrot%humiftrsveg
    print *, 'nh4_mass ', vrot%nh4_mass
    print *, 'no3_mass ', vrot%no3_mass
    print *, 'nlitrmass ', vrot%nlitrmass
    print *, 'soilnmas ', vrot%soilnmas
    print *, 'slopefrac ', vrot%slopefrac
    print *, 'dayl_max ', vrot%dayl_max
    print *, 'dayl ', vrot%dayl
    print *, 'grclarearow ', vrot%grclarearow
    print *, 'lygleafmasmax ', vrot%lygleafmasmax
    print *, 'lystemmassmax ', vrot%lystemmassmax
    print *, 'lyrootmassmax ', vrot%lyrootmassmax
    print *, 'lmaxt ', vrot%lmaxt
    print *, 'smaxt ', vrot%smaxt
    print *, 'rmaxt ', vrot%rmaxt

    print *, 'gleafmas ', vgat%gleafmas
    print *, 'gleafmas_ns ', vgat%gleafmas_ns
    print *, 'gleafmas_s ', vgat%gleafmas_s
    print *, 'leafns2s ', vgat%leafns2s
    print *, 'stemns2s ', vgat%stemns2s
    print *, 'rootns2s ', vgat%rootns2s
    print *, 're_alloc_s2l ', vgat%re_alloc_s2l
    print *, 're_alloc_r2l ', vgat%re_alloc_r2l
    print *, 're_alloc_sr2l ', vgat%re_alloc_sr2l
    print *, 'bleafmas ', vgat%bleafmas
    print *, 'stemmass ', vgat%stemmass
    print *, 'stemmass_ns ', vgat%stemmass_ns
    print *, 'stemmass_s ', vgat%stemmass_s
    print *, 'rootmass ', vgat%rootmass
    print *, 'rootmass_ns ', vgat%rootmass_ns
    print *, 'rootmass_s ', vgat%rootmass_s
    print *, 'pstemmass ', vgat%pstemmass
    print *, 'pgleafmass ', vgat%pgleafmass
    print *, 'fcancmx ', vgat%fcancmx
    print *, 'ngleafmas ', vgat%ngleafmas
    print *, 'ngleafmas_ns ', vgat%ngleafmas_ns
    print *, 'ngleafmas_s ', vgat%ngleafmas_s
    print *, 'nbleafmas ', vgat%nbleafmas
    print *, 'nstemmass ', vgat%nstemmass
    print *, 'nstemmass_ns ', vgat%nstemmass_ns
    print *, 'nstemmass_s ', vgat%nstemmass_s
    print *, 'nrootmass ', vgat%nrootmass
    print *, 'nrootmass_ns ', vgat%nrootmass_ns
    print *, 'nrootmass_s ', vgat%nrootmass_s
    print *, 'nvgbiomas_veg ', vgat%nvgbiomas_veg
    print *, 'gavglai ', vgat%gavglai
    print *, 'lightng ', vgat%lightng
    print *, 'zolnc ', vgat%zolnc
    print *, 'ailc ', vgat%ailc
    print *, 'ailcg ', vgat%ailcg
    print *, 'ailcgs ', vgat%ailcgs
    print *, 'fcancs ', vgat%fcancs
    print *, 'fcanc ', vgat%fcanc
    print *, 'co2conc ', vgat%co2conc
    print *, 'ch4conc ', vgat%ch4conc
    print *, 'co2i1cg ', vgat%co2i1cg
    print *, 'co2i1cs ', vgat%co2i1cs
    print *, 'co2i2cg ', vgat%co2i2cg
    print *, 'co2i2cs ', vgat%co2i2cs
    print *, 'ancsveg ', vgat%ancsveg
    print *, 'ancgveg ', vgat%ancgveg
    print *, 'rmlcsveg ', vgat%rmlcsveg
    print *, 'rmlcgveg ', vgat%rmlcgveg
    print *, 'slai ', vgat%slai
    print *, 'ailcb ', vgat%ailcb
    print *, 'canres ', vgat%canres
    print *, 'flhrloss ', vgat%flhrloss
    print *, 'flhrloss_ns ', vgat%flhrloss_ns
    print *, 'flhrloss_s ', vgat%flhrloss_s
    print *, 'grwtheff ', vgat%grwtheff
    print *, 'lystmmas ', vgat%lystmmas
    print *, 'lyrotmas ', vgat%lyrotmas
    print *, 'tymaxlai ', vgat%tymaxlai
    print *, 'vgbiomas ', vgat%vgbiomas
    print *, 'gavgltms ', vgat%gavgltms
    print *, 'gavgscms ', vgat%gavgscms
    print *, 'stmhrlos ', vgat%stmhrlos
    print *, 'rmatc ', vgat%rmatc
    print *, 'rmatctem ', vgat%rmatctem
    print *, 'litrmass ', vgat%litrmass
    print *, 'soilcmas ', vgat%soilcmas
    print *, 'vgbiomas_veg ', vgat%vgbiomas_veg
    print *, 'nh4_mass ', vgat%nh4_mass
    print *, 'no3_mass ', vgat%no3_mass
    print *, 'nlitrmass ', vgat%nlitrmass
    print *, 'soilnmas ', vgat%soilnmas
    print *, 'emit_co2 ', vgat%emit_co2
    print *, 'emit_co ', vgat%emit_co
    print *, 'emit_ch4 ', vgat%emit_ch4
    print *, 'emit_nmhc ', vgat%emit_nmhc
    print *, 'emit_h2 ', vgat%emit_h2
    print *, 'emit_nox ', vgat%emit_nox
    print *, 'emit_n2o ', vgat%emit_n2o
    print *, 'emit_pm25 ', vgat%emit_pm25
    print *, 'emit_tpm ', vgat%emit_tpm
    print *, 'emit_tc ', vgat%emit_tc
    print *, 'emit_oc ', vgat%emit_oc
    print *, 'emit_bc ', vgat%emit_bc
    print *, 'burnfrac ', vgat%burnfrac
    print *, 'timharvareagat ', vgat%timharvareagat
    print *, 'tileAgegat ', vgat%tileAgegat
    print *, 'burnvegf ', vgat%burnvegf
    print *, 'smfuncveg ', vgat%smfuncveg
    print *, 'popdin ', vgat%popdin
    print *, 'soilpH ', vgat%soilpH
    print *, 'nfertil ', vgat%nfertil
    print *, 'ndeposit ', vgat%ndeposit
    print *, 'bterm ', vgat%bterm
    print *, 'lterm ', vgat%lterm
    print *, 'mterm ', vgat%mterm
    print *, 'glcaemls ', vgat%glcaemls
    print *, 'blcaemls ', vgat%blcaemls
    print *, 'rtcaemls ', vgat%rtcaemls
    print *, 'stcaemls ', vgat%stcaemls
    print *, 'ltrcemls ', vgat%ltrcemls
    print *, 'ntchlveg ', vgat%ntchlveg
    print *, 'ntchsveg ', vgat%ntchsveg
    print *, 'extnprob ', vgat%extnprob
    print *, 'prbfrhuc ', vgat%prbfrhuc
    print *, 'dayl_max ', vgat%dayl_max
    print *, 'dayl ', vgat%dayl
    print *, 'bmasveg ', vgat%bmasveg
    print *, 'cmasvegc ', vgat%cmasvegc
    print *, 'veghght ', vgat%veghght
    print *, 'rootdpth ', vgat%rootdpth
    print *, 'rml ', vgat%rml
    print *, 'rms ', vgat%rms
    print *, 'tltrleaf ', vgat%tltrleaf
    print *, 'blfltrdt ', vgat%blfltrdt
    print *, 'glfltrdt ', vgat%glfltrdt
    print *, 'tltrstem ', vgat%tltrstem
    print *, 'tltrroot ', vgat%tltrroot
    print *, 'leaflitr ', vgat%leaflitr
    print *, 'roottemp ', vgat%roottemp
    print *, 'afrleaf ', vgat%afrleaf
    print *, 'afrstem ', vgat%afrstem
    print *, 'afrroot ', vgat%afrroot
    print *, 'wtstatus ', vgat%wtstatus
    print *, 'ltstatus ', vgat%ltstatus
    print *, 'rmr ', vgat%rmr
    print *, 'slopefrac ', vgat%slopefrac
    print *, 'wetfrac_pres ', vgat%wetfrac_pres
    print *, 'ch4WetSpec ', vgat%ch4WetSpec
    print *, 'wetfdyn ', vgat%wetfdyn
    print *, 'ch4WetDyn ', vgat%ch4WetDyn
    print *, 'ch4_soills ', vgat%ch4_soills
    print *, 'lucemcom ', vgat%lucemcom
    print *, 'lucltrin ', vgat%lucltrin
    print *, 'lucsocin ', vgat%lucsocin
    print *, 'npp ', vgat%npp
    print *, 'nep ', vgat%nep
    print *, 'nbp ', vgat%nbp
    print *, 'gpp ', vgat%gpp
    print *, 'hetrores ', vgat%hetrores
    print *, 'autores ', vgat%autores
    print *, 'soilcresp ', vgat%soilcresp
    print *, 'rm ', vgat%rm
    print *, 'rg ', vgat%rg
    print *, 'litres ', vgat%litres
    print *, 'socres ', vgat%socres
    print *, 'dstcemls ', vgat%dstcemls
    print *, 'litrfall ', vgat%litrfall
    print *, 'humiftrs ', vgat%humiftrs
    print *, 'gppveg ', vgat%gppveg
    print *, 'vcmax0 ', vgat%vcmax0
    print *, 'nepveg ', vgat%nepveg
    print *, 'bnf_free ', vgat%bnf_free
    print *, 'bnf_ant ', vgat%bnf_ant
    print *, 'bnf_nat ', vgat%bnf_nat
    print *, 'bnf_tot ', vgat%bnf_tot
    print *, 'nstress ', vgat%nstress
    print *, 'nitrifveg ', vgat%nitrifveg
    print *, 'no_nitveg ', vgat%no_nitveg
    print *, 'no_denitveg ', vgat%no_denitveg
    print *, 'no_nitdenitveg ', vgat%no_nitdenitveg
    print *, 'n2o_nitveg ', vgat%n2o_nitveg
    print *, 'n2o_denitveg ', vgat%n2o_denitveg
    print *, 'n2o_nitdenitveg ', vgat%n2o_nitdenitveg
    print *, 'n2_denitveg ', vgat%n2_denitveg
    print *, 'nvolveg ', vgat%nvolveg
    print *, 'nleachveg ', vgat%nleachveg
    print *, 'appl_fert ', vgat%appl_fert
    print *, 'ndep_nh4 ', vgat%ndep_nh4
    print *, 'ndep_no3 ', vgat%ndep_no3
    print *, 'ndemandveg_wp_npp ', vgat%ndemandveg_wp_npp
    print *, 'nuptakeveg_p_nh4 ', vgat%nuptakeveg_p_nh4
    print *, 'nuptakeveg_p_no3 ', vgat%nuptakeveg_p_no3
    print *, 'nuptakeveg_a_actl_nh4 ', vgat%nuptakeveg_a_actl_nh4
    print *, 'nuptakeveg_a_actl_no3 ', vgat%nuptakeveg_a_actl_no3
    print *, 'nuptakeveg ', vgat%nuptakeveg
    print *, 'nleafns2sveg ', vgat%nleafns2sveg
    print *, 'nstemns2sveg ', vgat%nstemns2sveg
    print *, 'nrootns2sveg ', vgat%nrootns2sveg
    print *, 'nallocveg_l ', vgat%nallocveg_l
    print *, 'nallocveg_s ', vgat%nallocveg_s
    print *, 'nallocveg_r ', vgat%nallocveg_r
    print *, 'nresorpedveg_s ', vgat%nresorpedveg_s
    print *, 'nresorpedveg_r ', vgat%nresorpedveg_r
    print *, 'nre_allocveg_s2l ', vgat%nre_allocveg_s2l
    print *, 'nre_allocveg_r2l ', vgat%nre_allocveg_r2l
    print *, 'nlitrveg_l ', vgat%nlitrveg_l
    print *, 'nlitrveg_s ', vgat%nlitrveg_s
    print *, 'nlitrveg_r ', vgat%nlitrveg_r
    print *, 'nlitrveg ', vgat%nlitrveg
    print *, 'gl2bl_grass_nflux ', vgat%gl2bl_grass_nflux
    print *, 'c2nveg_l ', vgat%c2nveg_l
    print *, 'c2nveg_s ', vgat%c2nveg_s
    print *, 'c2nveg_r ', vgat%c2nveg_r
    print *, 'c2nveg_wp ', vgat%c2nveg_wp
    print *, 'c2nveg_litr ', vgat%c2nveg_litr
    print *, 'c2nveg_humus ', vgat%c2nveg_humus
    print *, 'nhumtrsveg ', vgat%nhumtrsveg
    print *, 'nmineralveg_litr ', vgat%nmineralveg_litr
    print *, 'nmineralveg_humus ', vgat%nmineralveg_humus
    print *, 'netnmineralveg ', vgat%netnmineralveg
    print *, 'nimmobilveg_nh4 ', vgat%nimmobilveg_nh4
    print *, 'nimmobilveg_no3 ', vgat%nimmobilveg_no3
    print *, 'fNnetlandveg ', vgat%fNnetlandveg
    print *, 'redcoeff_vcmax ', vgat%redcoeff_vcmax
    print *, 'ipeatland ', vgat%ipeatland
    print *, 'peatdep ', vgat%peatdep
    print *, 'anmoss ', vgat%anmoss
    print *, 'rmlmoss ', vgat%rmlmoss
    print *, 'gppmoss ', vgat%gppmoss
    print *, 'nppmoss ', vgat%nppmoss
    print *, 'armoss ', vgat%armoss
    print *, 'litrmsmoss ', vgat%litrmsmoss
    print *, 'Cmossmas ', vgat%Cmossmas
    print *, 'dmoss ', vgat%dmoss
    print *, 'peatSoilC ', vgat%peatSoilC
    print *, 'pdd ', vgat%pdd
    print *, 'ancsmoss ', vgat%ancsmoss
    print *, 'angsmoss ', vgat%angsmoss
    print *, 'ancmoss ', vgat%ancmoss
    print *, 'angmoss ', vgat%angmoss
    print *, 'rmlcsmoss ', vgat%rmlcsmoss
    print *, 'rmlgsmoss ', vgat%rmlgsmoss
    print *, 'rmlcmoss ', vgat%rmlcmoss
    print *, 'rmlgmoss ', vgat%rmlgmoss
    print *, 'nbpveg ', vgat%nbpveg
    print *, 'nppveg ', vgat%nppveg
    print *, 'hetroresveg ', vgat%hetroresveg
    print *, 'autoresveg ', vgat%autoresveg
    print *, 'litresveg ', vgat%litresveg
    print *, 'soilcresveg ', vgat%soilcresveg
    print *, 'humiftrsveg ', vgat%humiftrsveg
    print *, 'rmlvegacc ', vgat%rmlvegacc
    print *, 'rmsveg ', vgat%rmsveg
    print *, 'rmrveg ', vgat%rmrveg
    print *, 'rgveg ', vgat%rgveg
    print *, 'litrfallveg ', vgat%litrfallveg
    print *, 'reprocost ', vgat%reprocost
    print *, 'rothrlos ', vgat%rothrlos
    print *, 'pfcancmx ', vgat%pfcancmx
    print *, 'nfcancmx ', vgat%nfcancmx
    print *, 'alvsctm ', vgat%alvsctm
    print *, 'paic ', vgat%paic
    print *, 'slaic ', vgat%slaic
    print *, 'alirctm ', vgat%alirctm
    print *, 'cfluxcg ', vgat%cfluxcg
    print *, 'cfluxcs ', vgat%cfluxcs
    print *, 'CFLUX_GA ', vgat%CFLUX_GA
    print *, 'USTARBS_GA ', vgat%USTARBS_GA
    print *, 'ROFB ', vgat%ROFB
    print *, 'dstcemls3 ', vgat%dstcemls3
    print *, 'anveg ', vgat%anveg
    print *, 'rmlveg ', vgat%rmlveg
    print *, 'twarmm ', vgat%twarmm
    print *, 'tcoldm ', vgat%tcoldm
    print *, 'gdd5 ', vgat%gdd5
    print *, 'aridity ', vgat%aridity
    print *, 'srplsmon ', vgat%srplsmon
    print *, 'defctmon ', vgat%defctmon
    print *, 'anndefct ', vgat%anndefct
    print *, 'annsrpls ', vgat%annsrpls
    print *, 'annpcp ', vgat%annpcp
    print *, 'dry_season_length ', vgat%dry_season_length
    print *, 'colddays_leaffall ', vgat%colddays_leaffall
    print *, 'colddays_harvest ', vgat%colddays_harvest
    print *, 'tcurm ', vgat%tcurm
    print *, 'srpcuryr ', vgat%srpcuryr
    print *, 'dftcuryr ', vgat%dftcuryr
    print *, 'tmonth ', vgat%tmonth
    print *, 'anpcpcur ', vgat%anpcpcur
    print *, 'anpecur ', vgat%anpecur
    print *, 'gdd5cur ', vgat%gdd5cur
    print *, 'surmncur ', vgat%surmncur
    print *, 'defmncur ', vgat%defmncur
    print *, 'srplscur ', vgat%srplscur
    print *, 'defctcur ', vgat%defctcur
    print *, 'geremort ', vgat%geremort
    print *, 'intrmort ', vgat%intrmort
    print *, 'cc ', vgat%cc
    print *, 'mm ', vgat%mm
    print *, 'pftexist ', vgat%pftexist
    print *, 'lfstatus ', vgat%lfstatus
    print *, 'pandays ', vgat%pandays
    print *, 'grclarea ', vgat%grclarea
    print *, 'altotcount_ctm ', vgat%altotcount_ctm
    print *, 'todfrac ', vgat%todfrac
    print *, 'fsinacc_gat ', vgat%fsinacc_gat
    print *, 'flutacc_gat ', vgat%flutacc_gat
    print *, 'flinacc_gat ', vgat%flinacc_gat
    print *, 'altotacc_gat ', vgat%altotacc_gat
    print *, 'netrad_gat ', vgat%netrad_gat
    print *, 'preacc_gat ', vgat%preacc_gat
    print *, 'sdepgat ', vgat%sdepgat
    print *, 'sandgat ', vgat%sandgat
    print *, 'claygat ', vgat%claygat
    print *, 'orgmgat ', vgat%orgmgat
    print *, 'xdiffusgat ', vgat%xdiffusgat
    print *, 'faregat ', vgat%faregat
    print *, 'lygleafmasmax ', vgat%lygleafmasmax
    print *, 'lystemmassmax ', vgat%lystemmassmax
    print *, 'lyrootmassmax ', vgat%lyrootmassmax
    print *, 'lmaxt ', vgat%lmaxt
    print *, 'smaxt ', vgat%smaxt
    print *, 'rmaxt ', vgat%rmaxt

    print *, 'mossCMassrot ', tracer%mossCMassrot
    print *, 'mossLitrMassrot ', tracer%mossLitrMassrot
    print *, 'tracerCO2rot ', tracer%tracerCO2rot
    print *, 'gLeafMassrot ', tracer%gLeafMassrot
    print *, 'bLeafMassrot ', tracer%bLeafMassrot
    print *, 'stemMassrot ', tracer%stemMassrot
    print *, 'rootMassrot ', tracer%rootMassrot
    print *, 'litrMassrot ', tracer%litrMassrot
    print *, 'soilCMassrot ', tracer%soilCMassrot
    print *, 'mossCMassgat ', tracer%mossCMassgat
    print *, 'mossLitrMassgat ', tracer%mossLitrMassgat
    print *, 'tracerCO2gat ', tracer%tracerCO2gat
    print *, 'gLeafMassgat ', tracer%gLeafMassgat
    print *, 'bLeafMassgat ', tracer%bLeafMassgat
    print *, 'stemMassgat ', tracer%stemMassgat
    print *, 'rootMassgat ', tracer%rootMassgat
    print *, 'litrMassgat ', tracer%litrMassgat
    print *, 'soilCMassgat ', tracer%soilCMassgat

    print *, 'taaccgat_t ', ctem_tile%taaccgat_t
    print *, 'uvaccgat_t ', ctem_tile%uvaccgat_t
    print *, 'vvaccgat_t ', ctem_tile%vvaccgat_t
    print *, 'anmossac_t ', ctem_tile%anmossac_t
    print *, 'rmlmossac_t ', ctem_tile%rmlmossac_t
    print *, 'gppmossac_t ', ctem_tile%gppmossac_t
    print *, 'QFCacc_t ', ctem_tile%QFCacc_t
    print *, 'tbaraccgat_t ', ctem_tile%tbaraccgat_t
    print *, 'thliqacc_t ', ctem_tile%thliqacc_t
    print *, 'thiceacc_t ', ctem_tile%thiceacc_t
    print *, 'ancgvgac_t ', ctem_tile%ancgvgac_t
    print *, 'rmlcgvga_t ', ctem_tile%rmlcgvga_t

    print *, 'laimaxg_mo ', ctem_mo%laimaxg_mo
    print *, 'gleafmas_mo ', ctem_mo%gleafmas_mo
    print *, 'gleafmas_NS_mo ', ctem_mo%gleafmas_NS_mo
    print *, 'gleafmass_mo ', ctem_mo%gleafmass_mo
    print *, 'stemmass_mo ', ctem_mo%stemmass_mo
    print *, 'stemmass_NS_mo ', ctem_mo%stemmass_NS_mo
    print *, 'stemmasss_mo ', ctem_mo%stemmasss_mo
    print *, 'rootmass_mo ', ctem_mo%rootmass_mo
    print *, 'rootmass_NS_mo ', ctem_mo%rootmass_NS_mo
    print *, 'rootmasss_mo ', ctem_mo%rootmasss_mo
    print *, 'rootdpth_mo ', ctem_mo%rootdpth_mo
    print *, 'litrfallveg_mo ', ctem_mo%litrfallveg_mo
    print *, 'humiftrsveg_mo ', ctem_mo%humiftrsveg_mo
    print *, 'tltrleaf_mo ', ctem_mo%tltrleaf_mo
    print *, 'tltrstem_mo ', ctem_mo%tltrstem_mo
    print *, 'tltrroot_mo ', ctem_mo%tltrroot_mo
    print *, 'npp_mo ', ctem_mo%npp_mo
    print *, 'gpp_mo ', ctem_mo%gpp_mo
    print *, 'vcmax0_mo ', ctem_mo%vcmax0_mo
    print *, 'vgbiomas_mo ', ctem_mo%vgbiomas_mo
    print *, 'autores_mo ', ctem_mo%autores_mo
    print *, 'soilres_mo ', ctem_mo%soilres_mo
    print *, 'totcmass_mo ', ctem_mo%totcmass_mo
    print *, 'nep_mo ', ctem_mo%nep_mo
    print *, 'hetrores_mo ', ctem_mo%hetrores_mo
    print *, 'nbp_mo ', ctem_mo%nbp_mo
    print *, 'emit_co2_mo ', ctem_mo%emit_co2_mo
    print *, 'emit_co_mo ', ctem_mo%emit_co_mo
    print *, 'emit_ch4_mo ', ctem_mo%emit_ch4_mo
    print *, 'emit_nmhc_mo ', ctem_mo%emit_nmhc_mo
    print *, 'emit_h2_mo ', ctem_mo%emit_h2_mo
    print *, 'emit_nox_mo ', ctem_mo%emit_nox_mo
    print *, 'emit_n2o_mo ', ctem_mo%emit_n2o_mo
    print *, 'emit_pm25_mo ', ctem_mo%emit_pm25_mo
    print *, 'emit_tpm_mo ', ctem_mo%emit_tpm_mo
    print *, 'emit_tc_mo ', ctem_mo%emit_tc_mo
    print *, 'emit_oc_mo ', ctem_mo%emit_oc_mo
    print *, 'emit_bc_mo ', ctem_mo%emit_bc_mo
    print *, 'burnfrac_mo ', ctem_mo%burnfrac_mo
    print *, 'bterm_mo ', ctem_mo%bterm_mo
    print *, 'mterm_mo ', ctem_mo%mterm_mo
    print *, 'smfuncveg_mo ', ctem_mo%smfuncveg_mo
    print *, 'nh4_mass_mo ', ctem_mo%nh4_mass_mo
    print *, 'no3_mass_mo ', ctem_mo%no3_mass_mo
    print *, 'ngleafmas_mo ', ctem_mo%ngleafmas_mo
    print *, 'ngleafmas_NS_mo ', ctem_mo%ngleafmas_NS_mo
    print *, 'ngleafmass_mo ', ctem_mo%ngleafmass_mo
    print *, 'nbleafmas_mo ', ctem_mo%nbleafmas_mo
    print *, 'nstemmass_mo ', ctem_mo%nstemmass_mo
    print *, 'nstemmass_NS_mo ', ctem_mo%nstemmass_NS_mo
    print *, 'nstemmasss_mo ', ctem_mo%nstemmasss_mo
    print *, 'nrootmass_mo ', ctem_mo%nrootmass_mo
    print *, 'nrootmass_NS_mo ', ctem_mo%nrootmass_NS_mo
    print *, 'nrootmasss_mo ', ctem_mo%nrootmasss_mo
    print *, 'nlitrmass_mo ', ctem_mo%nlitrmass_mo
    print *, 'soilnmas_mo ', ctem_mo%soilnmas_mo
    print *, 'nvgbiomas_mo ', ctem_mo%nvgbiomas_mo
    print *, 'litrmass_mo ', ctem_mo%litrmass_mo
    print *, 'soilcmas_mo ', ctem_mo%soilcmas_mo
    print *, 'litres_mo ', ctem_mo%litres_mo
    print *, 'soilcres_mo ', ctem_mo%soilcres_mo
    print *, 'leafns2s_mo ', ctem_mo%leafns2s_mo
    print *, 'stemns2s_mo ', ctem_mo%stemns2s_mo
    print *, 'rootns2s_mo ', ctem_mo%rootns2s_mo
    print *, 're_alloc_s2l_mo ', ctem_mo%re_alloc_s2l_mo
    print *, 're_alloc_r2l_mo ', ctem_mo%re_alloc_r2l_mo
    print *, 're_alloc_sr2l_mo ', ctem_mo%re_alloc_sr2l_mo
    print *, 'bnf_tot_mo ', ctem_mo%bnf_tot_mo
    print *, 'bnf_free_mo ', ctem_mo%bnf_free_mo
    print *, 'bnf_ant_mo ', ctem_mo%bnf_ant_mo
    print *, 'bnf_nat_mo ', ctem_mo%bnf_nat_mo
    print *, 'nstress_mo ', ctem_mo%nstress_mo
    print *, 'nitrif_mo ', ctem_mo%nitrif_mo
    print *, 'no_nit_mo ', ctem_mo%no_nit_mo
    print *, 'no_denit_mo ', ctem_mo%no_denit_mo
    print *, 'no_nitdenit_mo ', ctem_mo%no_nitdenit_mo
    print *, 'n2o_nit_mo ', ctem_mo%n2o_nit_mo
    print *, 'n2o_denit_mo ', ctem_mo%n2o_denit_mo
    print *, 'n2o_nitdenit_mo ', ctem_mo%n2o_nitdenit_mo
    print *, 'n2_denit_mo ', ctem_mo%n2_denit_mo
    print *, 'nvol_mo ', ctem_mo%nvol_mo
    print *, 'nleach_mo ', ctem_mo%nleach_mo
    print *, 'appl_fert_mo ', ctem_mo%appl_fert_mo
    print *, 'ndep_nh4_mo ', ctem_mo%ndep_nh4_mo
    print *, 'ndep_no3_mo ', ctem_mo%ndep_no3_mo
    print *, 'ndemand_wp_npp_mo ', ctem_mo%ndemand_wp_npp_mo
    print *, 'nuptake_p_nh4_mo ', ctem_mo%nuptake_p_nh4_mo
    print *, 'nuptake_p_no3_mo ', ctem_mo%nuptake_p_no3_mo
    print *, 'nuptake_a_actl_nh4_mo ', ctem_mo%nuptake_a_actl_nh4_mo
    print *, 'nuptake_a_actl_no3_mo ', ctem_mo%nuptake_a_actl_no3_mo
    print *, 'nuptake_mo ', ctem_mo%nuptake_mo
    print *, 'nleafns2s_mo ', ctem_mo%nleafns2s_mo
    print *, 'nstemns2s_mo ', ctem_mo%nstemns2s_mo
    print *, 'nrootns2s_mo ', ctem_mo%nrootns2s_mo
    print *, 'nalloc_l_mo ', ctem_mo%nalloc_l_mo
    print *, 'nalloc_s_mo ', ctem_mo%nalloc_s_mo
    print *, 'nalloc_r_mo ', ctem_mo%nalloc_r_mo
    print *, 'nresorped_s_mo ', ctem_mo%nresorped_s_mo
    print *, 'nresorped_r_mo ', ctem_mo%nresorped_r_mo
    print *, 'nre_alloc_s2l_mo ', ctem_mo%nre_alloc_s2l_mo
    print *, 'nre_alloc_r2l_mo ', ctem_mo%nre_alloc_r2l_mo
    print *, 'nlitr_l_mo ', ctem_mo%nlitr_l_mo
    print *, 'nlitr_s_mo ', ctem_mo%nlitr_s_mo
    print *, 'nlitr_r_mo ', ctem_mo%nlitr_r_mo
    print *, 'nlitr_mo ', ctem_mo%nlitr_mo
    print *, 'gl2bl_grass_nflux_mo ', ctem_mo%gl2bl_grass_nflux_mo
    print *, 'c2n_l_mo ', ctem_mo%c2n_l_mo
    print *, 'c2n_s_mo ', ctem_mo%c2n_s_mo
    print *, 'c2n_r_mo ', ctem_mo%c2n_r_mo
    print *, 'c2n_wp_mo ', ctem_mo%c2n_wp_mo
    print *, 'c2n_litr_mo ', ctem_mo%c2n_litr_mo
    print *, 'c2n_humus_mo ', ctem_mo%c2n_humus_mo
    print *, 'nhumtrs_mo ', ctem_mo%nhumtrs_mo
    print *, 'nmineral_litr_mo ', ctem_mo%nmineral_litr_mo
    print *, 'nmineral_humus_mo ', ctem_mo%nmineral_humus_mo
    print *, 'netnmineral_mo ', ctem_mo%netnmineral_mo
    print *, 'nimmobil_nh4_mo ', ctem_mo%nimmobil_nh4_mo
    print *, 'nimmobil_no3_mo ', ctem_mo%nimmobil_no3_mo
    print *, 'fNnetland_mo ', ctem_mo%fNnetland_mo

    print *, 'laimaxg_mo_g ', ctem_grd_mo%laimaxg_mo_g
    print *, 'gleafmas_mo_g ', ctem_grd_mo%gleafmas_mo_g
    print *, 'gleafmas_NS_mo_g ', ctem_grd_mo%gleafmas_NS_mo_g
    print *, 'gleafmass_mo_g ', ctem_grd_mo%gleafmass_mo_g
    print *, 'stemmass_mo_g ', ctem_grd_mo%stemmass_mo_g
    print *, 'stemmass_NS_mo_g ', ctem_grd_mo%stemmass_NS_mo_g
    print *, 'stemmasss_mo_g ', ctem_grd_mo%stemmasss_mo_g
    print *, 'rootmass_mo_g ', ctem_grd_mo%rootmass_mo_g
    print *, 'rootmass_NS_mo_g ', ctem_grd_mo%rootmass_NS_mo_g
    print *, 'rootmasss_mo_g ', ctem_grd_mo%rootmasss_mo_g
    print *, 'nh4_mass_mo_g ', ctem_grd_mo%nh4_mass_mo_g
    print *, 'no3_mass_mo_g ', ctem_grd_mo%no3_mass_mo_g
    print *, 'ngleafmas_mo_g ', ctem_grd_mo%ngleafmas_mo_g
    print *, 'ngleafmas_NS_mo_g ', ctem_grd_mo%ngleafmas_NS_mo_g
    print *, 'ngleafmass_mo_g ', ctem_grd_mo%ngleafmass_mo_g
    print *, 'nbleafmas_mo_g ', ctem_grd_mo%nbleafmas_mo_g
    print *, 'nstemmass_mo_g ', ctem_grd_mo%nstemmass_mo_g
    print *, 'nstemmass_NS_mo_g ', ctem_grd_mo%nstemmass_NS_mo_g
    print *, 'nstemmasss_mo_g ', ctem_grd_mo%nstemmasss_mo_g
    print *, 'nrootmass_mo_g ', ctem_grd_mo%nrootmass_mo_g
    print *, 'nrootmass_NS_mo_g ', ctem_grd_mo%nrootmass_NS_mo_g
    print *, 'nrootmasss_mo_g ', ctem_grd_mo%nrootmasss_mo_g
    print *, 'nlitrmass_mo_g ', ctem_grd_mo%nlitrmass_mo_g
    print *, 'soilnmas_mo_g ', ctem_grd_mo%soilnmas_mo_g
    print *, 'nvgbiomas_mo_g ', ctem_grd_mo%nvgbiomas_mo_g
    print *, 'litrfall_mo_g ', ctem_grd_mo%litrfall_mo_g
    print *, 'humiftrs_mo_g ', ctem_grd_mo%humiftrs_mo_g
    print *, 'tltrleaf_mo_g ', ctem_grd_mo%tltrleaf_mo_g
    print *, 'tltrstem_mo_g ', ctem_grd_mo%tltrstem_mo_g
    print *, 'tltrroot_mo_g ', ctem_grd_mo%tltrroot_mo_g
    print *, 'npp_mo_g ', ctem_grd_mo%npp_mo_g
    print *, 'gpp_mo_g ', ctem_grd_mo%gpp_mo_g
    print *, 'vcmax0_mo_g ', ctem_grd_mo%vcmax0_mo_g
    print *, 'nep_mo_g ', ctem_grd_mo%nep_mo_g
    print *, 'nbp_mo_g ', ctem_grd_mo%nbp_mo_g
    print *, 'hetrores_mo_g ', ctem_grd_mo%hetrores_mo_g
    print *, 'autores_mo_g ', ctem_grd_mo%autores_mo_g
    print *, 'soilres_mo_g ', ctem_grd_mo%soilres_mo_g
    print *, 'vgbiomas_mo_g ', ctem_grd_mo%vgbiomas_mo_g
    print *, 'totcmass_mo_g ', ctem_grd_mo%totcmass_mo_g
    print *, 'emit_co2_mo_g ', ctem_grd_mo%emit_co2_mo_g
    print *, 'emit_co_mo_g ', ctem_grd_mo%emit_co_mo_g
    print *, 'emit_ch4_mo_g ', ctem_grd_mo%emit_ch4_mo_g
    print *, 'emit_nmhc_mo_g ', ctem_grd_mo%emit_nmhc_mo_g
    print *, 'emit_h2_mo_g ', ctem_grd_mo%emit_h2_mo_g
    print *, 'emit_nox_mo_g ', ctem_grd_mo%emit_nox_mo_g
    print *, 'emit_n2o_mo_g ', ctem_grd_mo%emit_n2o_mo_g
    print *, 'emit_pm25_mo_g ', ctem_grd_mo%emit_pm25_mo_g
    print *, 'emit_tpm_mo_g ', ctem_grd_mo%emit_tpm_mo_g
    print *, 'emit_tc_mo_g ', ctem_grd_mo%emit_tc_mo_g
    print *, 'emit_oc_mo_g ', ctem_grd_mo%emit_oc_mo_g
    print *, 'emit_bc_mo_g ', ctem_grd_mo%emit_bc_mo_g
    print *, 'smfuncveg_mo_g ', ctem_grd_mo%smfuncveg_mo_g
    print *, 'luc_emc_mo_g ', ctem_grd_mo%luc_emc_mo_g
    print *, 'lucltrin_mo_g ', ctem_grd_mo%lucltrin_mo_g
    print *, 'lucsocin_mo_g ', ctem_grd_mo%lucsocin_mo_g
    print *, 'burnfrac_mo_g ', ctem_grd_mo%burnfrac_mo_g
    print *, 'bterm_mo_g ', ctem_grd_mo%bterm_mo_g
    print *, 'lterm_mo_g ', ctem_grd_mo%lterm_mo_g
    print *, 'mterm_mo_g ', ctem_grd_mo%mterm_mo_g
    print *, 'ch4WetSpec_mo_g ', ctem_grd_mo%ch4WetSpec_mo_g
    print *, 'wetfdyn_mo_g ', ctem_grd_mo%wetfdyn_mo_g
    print *, 'wetfpres_mo_g ', ctem_grd_mo%wetfpres_mo_g
    print *, 'ch4WetDyn_mo_g ', ctem_grd_mo%ch4WetDyn_mo_g
    print *, 'ch4soills_mo_g ', ctem_grd_mo%ch4soills_mo_g
    print *, 'cProduct_mo_g ', ctem_grd_mo%cProduct_mo_g
    print *, 'nProduct_mo_g ', ctem_grd_mo%nProduct_mo_g
    print *, 'fProductDecomp_mo_g ', ctem_grd_mo%fProductDecomp_mo_g
    print *, 'litrmass_mo_g ', ctem_grd_mo%litrmass_mo_g
    print *, 'soilcmas_mo_g ', ctem_grd_mo%soilcmas_mo_g
    print *, 'litres_mo_g ', ctem_grd_mo%litres_mo_g
    print *, 'soilcres_mo_g ', ctem_grd_mo%soilcres_mo_g
    print *, 'leafns2s_mo_g ', ctem_grd_mo%leafns2s_mo_g
    print *, 'stemns2s_mo_g ', ctem_grd_mo%stemns2s_mo_g
    print *, 'rootns2s_mo_g ', ctem_grd_mo%rootns2s_mo_g
    print *, 're_alloc_s2l_mo_g ', ctem_grd_mo%re_alloc_s2l_mo_g
    print *, 're_alloc_r2l_mo_g ', ctem_grd_mo%re_alloc_r2l_mo_g
    print *, 're_alloc_sr2l_mo_g ', ctem_grd_mo%re_alloc_sr2l_mo_g
    print *, 'bnf_tot_mo_g ', ctem_grd_mo%bnf_tot_mo_g
    print *, 'bnf_free_mo_g ', ctem_grd_mo%bnf_free_mo_g
    print *, 'bnf_ant_mo_g ', ctem_grd_mo%bnf_ant_mo_g
    print *, 'bnf_nat_mo_g ', ctem_grd_mo%bnf_nat_mo_g
    print *, 'nstress_mo_g ', ctem_grd_mo%nstress_mo_g
    print *, 'nitrif_mo_g ', ctem_grd_mo%nitrif_mo_g
    print *, 'no_nit_mo_g ', ctem_grd_mo%no_nit_mo_g
    print *, 'no_denit_mo_g ', ctem_grd_mo%no_denit_mo_g
    print *, 'no_nitdenit_mo_g ', ctem_grd_mo%no_nitdenit_mo_g
    print *, 'n2o_nit_mo_g ', ctem_grd_mo%n2o_nit_mo_g
    print *, 'n2o_denit_mo_g ', ctem_grd_mo%n2o_denit_mo_g
    print *, 'n2o_nitdenit_mo_g ', ctem_grd_mo%n2o_nitdenit_mo_g
    print *, 'n2_denit_mo_g ', ctem_grd_mo%n2_denit_mo_g
    print *, 'nvol_mo_g ', ctem_grd_mo%nvol_mo_g
    print *, 'nleach_mo_g ', ctem_grd_mo%nleach_mo_g
    print *, 'appl_fert_mo_g ', ctem_grd_mo%appl_fert_mo_g
    print *, 'ndep_nh4_mo_g ', ctem_grd_mo%ndep_nh4_mo_g
    print *, 'ndep_no3_mo_g ', ctem_grd_mo%ndep_no3_mo_g
    print *, 'ndemand_wp_npp_mo_g ', ctem_grd_mo%ndemand_wp_npp_mo_g
    print *, 'nuptake_p_nh4_mo_g ', ctem_grd_mo%nuptake_p_nh4_mo_g
    print *, 'nuptake_p_no3_mo_g ', ctem_grd_mo%nuptake_p_no3_mo_g
    print *, 'nuptake_a_actl_nh4_mo_g ', ctem_grd_mo%nuptake_a_actl_nh4_mo_g
    print *, 'nuptake_a_actl_no3_mo_g ', ctem_grd_mo%nuptake_a_actl_no3_mo_g
    print *, 'nuptake_mo_g ', ctem_grd_mo%nuptake_mo_g
    print *, 'nleafns2s_mo_g ', ctem_grd_mo%nleafns2s_mo_g
    print *, 'nstemns2s_mo_g ', ctem_grd_mo%nstemns2s_mo_g
    print *, 'nrootns2s_mo_g ', ctem_grd_mo%nrootns2s_mo_g
    print *, 'nalloc_l_mo_g ', ctem_grd_mo%nalloc_l_mo_g
    print *, 'nalloc_s_mo_g ', ctem_grd_mo%nalloc_s_mo_g
    print *, 'nalloc_r_mo_g ', ctem_grd_mo%nalloc_r_mo_g
    print *, 'nresorped_s_mo_g ', ctem_grd_mo%nresorped_s_mo_g
    print *, 'nresorped_r_mo_g ', ctem_grd_mo%nresorped_r_mo_g
    print *, 'nre_alloc_s2l_mo_g ', ctem_grd_mo%nre_alloc_s2l_mo_g
    print *, 'nre_alloc_r2l_mo_g ', ctem_grd_mo%nre_alloc_r2l_mo_g
    print *, 'nlitr_l_mo_g ', ctem_grd_mo%nlitr_l_mo_g
    print *, 'nlitr_s_mo_g ', ctem_grd_mo%nlitr_s_mo_g
    print *, 'nlitr_r_mo_g ', ctem_grd_mo%nlitr_r_mo_g
    print *, 'nlitr_mo_g ', ctem_grd_mo%nlitr_mo_g
    print *, 'gl2bl_grass_nflux_mo_g ', ctem_grd_mo%gl2bl_grass_nflux_mo_g
    print *, 'c2n_l_mo_g ', ctem_grd_mo%c2n_l_mo_g
    print *, 'c2n_s_mo_g ', ctem_grd_mo%c2n_s_mo_g
    print *, 'c2n_r_mo_g ', ctem_grd_mo%c2n_r_mo_g
    print *, 'c2n_wp_mo_g ', ctem_grd_mo%c2n_wp_mo_g
    print *, 'c2n_litr_mo_g ', ctem_grd_mo%c2n_litr_mo_g
    print *, 'c2n_humus_mo_g ', ctem_grd_mo%c2n_humus_mo_g
    print *, 'nhumtrs_mo_g ', ctem_grd_mo%nhumtrs_mo_g
    print *, 'nmineral_litr_mo_g ', ctem_grd_mo%nmineral_litr_mo_g
    print *, 'nmineral_humus_mo_g ', ctem_grd_mo%nmineral_humus_mo_g
    print *, 'netnmineral_mo_g ', ctem_grd_mo%netnmineral_mo_g
    print *, 'nimmobil_nh4_mo_g ', ctem_grd_mo%nimmobil_nh4_mo_g
    print *, 'nimmobil_no3_mo_g ', ctem_grd_mo%nimmobil_no3_mo_g
    print *, 'fNnetland_mo_g ', ctem_grd_mo%fNnetland_mo_g

    print *, 'laimaxg_mo_t ', ctem_tile_mo%laimaxg_mo_t
    print *, 'gleafmas_mo_t ', ctem_tile_mo%gleafmas_mo_t
    print *, 'gleafmas_NS_mo_t ', ctem_tile_mo%gleafmas_NS_mo_t
    print *, 'gleafmass_mo_t ', ctem_tile_mo%gleafmass_mo_t
    print *, 'stemmass_mo_t ', ctem_tile_mo%stemmass_mo_t
    print *, 'stemmass_NS_mo_t ', ctem_tile_mo%stemmass_NS_mo_t
    print *, 'stemmasss_mo_t ', ctem_tile_mo%stemmasss_mo_t
    print *, 'rootmass_mo_t ', ctem_tile_mo%rootmass_mo_t
    print *, 'rootmass_NS_mo_t ', ctem_tile_mo%rootmass_NS_mo_t
    print *, 'rootmasss_mo_t ', ctem_tile_mo%rootmasss_mo_t
    print *, 'litrfall_mo_t ', ctem_tile_mo%litrfall_mo_t
    print *, 'humiftrs_mo_t ', ctem_tile_mo%humiftrs_mo_t
    print *, 'tltrleaf_mo_t ', ctem_tile_mo%tltrleaf_mo_t
    print *, 'tltrstem_mo_t ', ctem_tile_mo%tltrstem_mo_t
    print *, 'tltrroot_mo_t ', ctem_tile_mo%tltrroot_mo_t
    print *, 'npp_mo_t ', ctem_tile_mo%npp_mo_t
    print *, 'gpp_mo_t ', ctem_tile_mo%gpp_mo_t
    print *, 'vcmax0_mo_t ', ctem_tile_mo%vcmax0_mo_t
    print *, 'vgbiomas_mo_t ', ctem_tile_mo%vgbiomas_mo_t
    print *, 'autores_mo_t ', ctem_tile_mo%autores_mo_t
    print *, 'soilres_mo_t ', ctem_tile_mo%soilres_mo_t
    print *, 'totcmass_mo_t ', ctem_tile_mo%totcmass_mo_t
    print *, 'nh4_mass_mo_t ', ctem_tile_mo%nh4_mass_mo_t
    print *, 'no3_mass_mo_t ', ctem_tile_mo%no3_mass_mo_t
    print *, 'ngleafmas_mo_t ', ctem_tile_mo%ngleafmas_mo_t
    print *, 'ngleafmas_NS_mo_t ', ctem_tile_mo%ngleafmas_NS_mo_t
    print *, 'ngleafmass_mo_t ', ctem_tile_mo%ngleafmass_mo_t
    print *, 'nbleafmas_mo_t ', ctem_tile_mo%nbleafmas_mo_t
    print *, 'nstemmass_mo_t ', ctem_tile_mo%nstemmass_mo_t
    print *, 'nstemmass_NS_mo_t ', ctem_tile_mo%nstemmass_NS_mo_t
    print *, 'nstemmasss_mo_t ', ctem_tile_mo%nstemmasss_mo_t
    print *, 'nrootmass_mo_t ', ctem_tile_mo%nrootmass_mo_t
    print *, 'nrootmass_NS_mo_t ', ctem_tile_mo%nrootmass_NS_mo_t
    print *, 'nrootmasss_mo_t ', ctem_tile_mo%nrootmasss_mo_t
    print *, 'nlitrmass_mo_t ', ctem_tile_mo%nlitrmass_mo_t
    print *, 'soilnmas_mo_t ', ctem_tile_mo%soilnmas_mo_t
    print *, 'nvgbiomas_mo_t ', ctem_tile_mo%nvgbiomas_mo_t
    print *, 'nep_mo_t ', ctem_tile_mo%nep_mo_t
    print *, 'hetrores_mo_t ', ctem_tile_mo%hetrores_mo_t
    print *, 'nbp_mo_t ', ctem_tile_mo%nbp_mo_t
    print *, 'emit_co2_mo_t ', ctem_tile_mo%emit_co2_mo_t
    print *, 'emit_co_mo_t ', ctem_tile_mo%emit_co_mo_t
    print *, 'emit_ch4_mo_t ', ctem_tile_mo%emit_ch4_mo_t
    print *, 'emit_nmhc_mo_t ', ctem_tile_mo%emit_nmhc_mo_t
    print *, 'emit_h2_mo_t ', ctem_tile_mo%emit_h2_mo_t
    print *, 'emit_nox_mo_t ', ctem_tile_mo%emit_nox_mo_t
    print *, 'emit_n2o_mo_t ', ctem_tile_mo%emit_n2o_mo_t
    print *, 'emit_pm25_mo_t ', ctem_tile_mo%emit_pm25_mo_t
    print *, 'emit_tpm_mo_t ', ctem_tile_mo%emit_tpm_mo_t
    print *, 'emit_tc_mo_t ', ctem_tile_mo%emit_tc_mo_t
    print *, 'emit_oc_mo_t ', ctem_tile_mo%emit_oc_mo_t
    print *, 'emit_bc_mo_t ', ctem_tile_mo%emit_bc_mo_t
    print *, 'burnfrac_mo_t ', ctem_tile_mo%burnfrac_mo_t
    print *, 'smfuncveg_mo_t ', ctem_tile_mo%smfuncveg_mo_t
    print *, 'bterm_mo_t ', ctem_tile_mo%bterm_mo_t
    print *, 'luc_emc_mo_t ', ctem_tile_mo%luc_emc_mo_t
    print *, 'lterm_mo_t ', ctem_tile_mo%lterm_mo_t
    print *, 'lucsocin_mo_t ', ctem_tile_mo%lucsocin_mo_t
    print *, 'mterm_mo_t ', ctem_tile_mo%mterm_mo_t
    print *, 'lucltrin_mo_t ', ctem_tile_mo%lucltrin_mo_t
    print *, 'ch4WetSpec_mo_t ', ctem_tile_mo%ch4WetSpec_mo_t
    print *, 'wetfdyn_mo_t ', ctem_tile_mo%wetfdyn_mo_t
    print *, 'wetfpres_mo_t ', ctem_tile_mo%wetfpres_mo_t
    print *, 'ch4WetDyn_mo_t ', ctem_tile_mo%ch4WetDyn_mo_t
    print *, 'ch4soills_mo_t ', ctem_tile_mo%ch4soills_mo_t
    print *, 'wind_mo_t ', ctem_tile_mo%wind_mo_t
    print *, 'fProductDecomp_mo_t ', ctem_tile_mo%fProductDecomp_mo_t
    print *, 'litrmass_mo_t ', ctem_tile_mo%litrmass_mo_t
    print *, 'soilcmas_mo_t ', ctem_tile_mo%soilcmas_mo_t
    print *, 'litres_mo_t ', ctem_tile_mo%litres_mo_t
    print *, 'soilcres_mo_t ', ctem_tile_mo%soilcres_mo_t
    print *, 'leafns2s_mo_t ', ctem_tile_mo%leafns2s_mo_t
    print *, 'stemns2s_mo_t ', ctem_tile_mo%stemns2s_mo_t
    print *, 'rootns2s_mo_t ', ctem_tile_mo%rootns2s_mo_t
    print *, 're_alloc_s2l_mo_t ', ctem_tile_mo%re_alloc_s2l_mo_t
    print *, 're_alloc_r2l_mo_t ', ctem_tile_mo%re_alloc_r2l_mo_t
    print *, 're_alloc_sr2l_mo_t ', ctem_tile_mo%re_alloc_sr2l_mo_t
    print *, 'bnf_tot_mo_t ', ctem_tile_mo%bnf_tot_mo_t
    print *, 'bnf_free_mo_t ', ctem_tile_mo%bnf_free_mo_t
    print *, 'bnf_ant_mo_t ', ctem_tile_mo%bnf_ant_mo_t
    print *, 'bnf_nat_mo_t ', ctem_tile_mo%bnf_nat_mo_t
    print *, 'nstress_mo_t ', ctem_tile_mo%nstress_mo_t
    print *, 'nitrif_mo_t ', ctem_tile_mo%nitrif_mo_t
    print *, 'no_nit_mo_t ', ctem_tile_mo%no_nit_mo_t
    print *, 'no_denit_mo_t ', ctem_tile_mo%no_denit_mo_t
    print *, 'no_nitdenit_mo_t ', ctem_tile_mo%no_nitdenit_mo_t
    print *, 'n2o_nit_mo_t ', ctem_tile_mo%n2o_nit_mo_t
    print *, 'n2o_denit_mo_t ', ctem_tile_mo%n2o_denit_mo_t
    print *, 'n2o_nitdenit_mo_t ', ctem_tile_mo%n2o_nitdenit_mo_t
    print *, 'n2_denit_mo_t ', ctem_tile_mo%n2_denit_mo_t
    print *, 'nvol_mo_t ', ctem_tile_mo%nvol_mo_t
    print *, 'nleach_mo_t ', ctem_tile_mo%nleach_mo_t
    print *, 'appl_fert_mo_t ', ctem_tile_mo%appl_fert_mo_t
    print *, 'ndep_nh4_mo_t ', ctem_tile_mo%ndep_nh4_mo_t
    print *, 'ndep_no3_mo_t ', ctem_tile_mo%ndep_no3_mo_t
    print *, 'ndemand_wp_npp_mo_t ', ctem_tile_mo%ndemand_wp_npp_mo_t
    print *, 'nuptake_p_nh4_mo_t ', ctem_tile_mo%nuptake_p_nh4_mo_t
    print *, 'nuptake_p_no3_mo_t ', ctem_tile_mo%nuptake_p_no3_mo_t
    print *, 'nuptake_a_actl_nh4_mo_t ', ctem_tile_mo%nuptake_a_actl_nh4_mo_t
    print *, 'nuptake_a_actl_no3_mo_t ', ctem_tile_mo%nuptake_a_actl_no3_mo_t
    print *, 'nuptake_mo_t ', ctem_tile_mo%nuptake_mo_t
    print *, 'nleafns2s_mo_t ', ctem_tile_mo%nleafns2s_mo_t
    print *, 'nstemns2s_mo_t ', ctem_tile_mo%nstemns2s_mo_t
    print *, 'nrootns2s_mo_t ', ctem_tile_mo%nrootns2s_mo_t
    print *, 'nalloc_l_mo_t ', ctem_tile_mo%nalloc_l_mo_t
    print *, 'nalloc_s_mo_t ', ctem_tile_mo%nalloc_s_mo_t
    print *, 'nalloc_r_mo_t ', ctem_tile_mo%nalloc_r_mo_t
    print *, 'nresorped_s_mo_t ', ctem_tile_mo%nresorped_s_mo_t
    print *, 'nresorped_r_mo_t ', ctem_tile_mo%nresorped_r_mo_t
    print *, 'nre_alloc_s2l_mo_t ', ctem_tile_mo%nre_alloc_s2l_mo_t
    print *, 'nre_alloc_r2l_mo_t ', ctem_tile_mo%nre_alloc_r2l_mo_t
    print *, 'nlitr_l_mo_t ', ctem_tile_mo%nlitr_l_mo_t
    print *, 'nlitr_s_mo_t ', ctem_tile_mo%nlitr_s_mo_t
    print *, 'nlitr_r_mo_t ', ctem_tile_mo%nlitr_r_mo_t
    print *, 'nlitr_mo_t ', ctem_tile_mo%nlitr_mo_t
    print *, 'gl2bl_grass_nflux_mo_t ', ctem_tile_mo%gl2bl_grass_nflux_mo_t
    print *, 'c2n_l_mo_t ', ctem_tile_mo%c2n_l_mo_t
    print *, 'c2n_s_mo_t ', ctem_tile_mo%c2n_s_mo_t
    print *, 'c2n_r_mo_t ', ctem_tile_mo%c2n_r_mo_t
    print *, 'c2n_wp_mo_t ', ctem_tile_mo%c2n_wp_mo_t
    print *, 'c2n_litr_mo_t ', ctem_tile_mo%c2n_litr_mo_t
    print *, 'c2n_humus_mo_t ', ctem_tile_mo%c2n_humus_mo_t
    print *, 'nhumtrs_mo_t ', ctem_tile_mo%nhumtrs_mo_t
    print *, 'nmineral_litr_mo_t ', ctem_tile_mo%nmineral_litr_mo_t
    print *, 'nmineral_humus_mo_t ', ctem_tile_mo%nmineral_humus_mo_t
    print *, 'netnmineral_mo_t ', ctem_tile_mo%netnmineral_mo_t
    print *, 'nimmobil_nh4_mo_t ', ctem_tile_mo%nimmobil_nh4_mo_t
    print *, 'nimmobil_no3_mo_t ', ctem_tile_mo%nimmobil_no3_mo_t
    print *, 'fNnetland_mo_t ', ctem_tile_mo%fNnetland_mo_t

    print *, 'laimaxg_yr ', ctem_yr%laimaxg_yr
    print *, 'gleafmas_yr ', ctem_yr%gleafmas_yr
    print *, 'gleafmas_NS_yr ', ctem_yr%gleafmas_NS_yr
    print *, 'gleafmass_yr ', ctem_yr%gleafmass_yr
    print *, 'bleafmas_yr ', ctem_yr%bleafmas_yr
    print *, 'stemmass_yr ', ctem_yr%stemmass_yr
    print *, 'stemmass_NS_yr ', ctem_yr%stemmass_NS_yr
    print *, 'stemmasss_yr ', ctem_yr%stemmasss_yr
    print *, 'rootmass_yr ', ctem_yr%rootmass_yr
    print *, 'rootmass_NS_yr ', ctem_yr%rootmass_NS_yr
    print *, 'rootmasss_yr ', ctem_yr%rootmasss_yr
    print *, 'npp_yr ', ctem_yr%npp_yr
    print *, 'gpp_yr ', ctem_yr%gpp_yr
    print *, 'vcmax0_yr ', ctem_yr%vcmax0_yr
    print *, 'leafns2s_yr ', ctem_yr%leafns2s_yr
    print *, 'stemns2s_yr ', ctem_yr%stemns2s_yr
    print *, 'rootns2s_yr ', ctem_yr%rootns2s_yr
    print *, 're_alloc_s2l_yr ', ctem_yr%re_alloc_s2l_yr
    print *, 're_alloc_r2l_yr ', ctem_yr%re_alloc_r2l_yr
    print *, 're_alloc_sr2l_yr ', ctem_yr%re_alloc_sr2l_yr
    print *, 'vgbiomas_yr ', ctem_yr%vgbiomas_yr
    print *, 'autores_yr ', ctem_yr%autores_yr
    print *, 'rmrveg_yr ', ctem_yr%rmrveg_yr
    print *, 'totcmass_yr ', ctem_yr%totcmass_yr
    print *, 'nh4_mass_yr ', ctem_yr%nh4_mass_yr
    print *, 'no3_mass_yr ', ctem_yr%no3_mass_yr
    print *, 'ngleafmas_yr ', ctem_yr%ngleafmas_yr
    print *, 'ngleafmas_NS_yr ', ctem_yr%ngleafmas_NS_yr
    print *, 'ngleafmass_yr ', ctem_yr%ngleafmass_yr
    print *, 'nbleafmas_yr ', ctem_yr%nbleafmas_yr
    print *, 'nstemmass_yr ', ctem_yr%nstemmass_yr
    print *, 'nstemmass_NS_yr ', ctem_yr%nstemmass_NS_yr
    print *, 'nstemmasss_yr ', ctem_yr%nstemmasss_yr
    print *, 'nrootmass_yr ', ctem_yr%nrootmass_yr
    print *, 'nrootmass_NS_yr ', ctem_yr%nrootmass_NS_yr
    print *, 'nrootmasss_yr ', ctem_yr%nrootmasss_yr
    print *, 'nlitrmass_yr ', ctem_yr%nlitrmass_yr
    print *, 'soilnmas_yr ', ctem_yr%soilnmas_yr
    print *, 'nvgbiomas_yr ', ctem_yr%nvgbiomas_yr
    print *, 'nep_yr ', ctem_yr%nep_yr
    print *, 'litrfall_yr ', ctem_yr%litrfall_yr
    print *, 'hetrores_yr ', ctem_yr%hetrores_yr
    print *, 'nbp_yr ', ctem_yr%nbp_yr
    print *, 'emit_co2_yr ', ctem_yr%emit_co2_yr
    print *, 'emit_co_yr ', ctem_yr%emit_co_yr
    print *, 'emit_ch4_yr ', ctem_yr%emit_ch4_yr
    print *, 'emit_nmhc_yr ', ctem_yr%emit_nmhc_yr
    print *, 'emit_h2_yr ', ctem_yr%emit_h2_yr
    print *, 'emit_nox_yr ', ctem_yr%emit_nox_yr
    print *, 'emit_n2o_yr ', ctem_yr%emit_n2o_yr
    print *, 'emit_pm25_yr ', ctem_yr%emit_pm25_yr
    print *, 'emit_tpm_yr ', ctem_yr%emit_tpm_yr
    print *, 'emit_tc_yr ', ctem_yr%emit_tc_yr
    print *, 'emit_oc_yr ', ctem_yr%emit_oc_yr
    print *, 'emit_bc_yr ', ctem_yr%emit_bc_yr
    print *, 'bterm_yr ', ctem_yr%bterm_yr
    print *, 'mterm_yr ', ctem_yr%mterm_yr
    print *, 'burnfrac_yr ', ctem_yr%burnfrac_yr
    print *, 'smfuncveg_yr ', ctem_yr%smfuncveg_yr
    print *, 'veghght_yr ', ctem_yr%veghght_yr
    print *, 'litrmass_yr ', ctem_yr%litrmass_yr
    print *, 'soilcmas_yr ', ctem_yr%soilcmas_yr
    print *, 'litres_yr ', ctem_yr%litres_yr
    print *, 'soilcres_yr ', ctem_yr%soilcres_yr
    print *, 'bnf_tot_yr ', ctem_yr%bnf_tot_yr
    print *, 'bnf_free_yr ', ctem_yr%bnf_free_yr
    print *, 'bnf_ant_yr ', ctem_yr%bnf_ant_yr
    print *, 'bnf_nat_yr ', ctem_yr%bnf_nat_yr
    print *, 'nstress_yr ', ctem_yr%nstress_yr
    print *, 'nitrif_yr ', ctem_yr%nitrif_yr
    print *, 'no_nit_yr ', ctem_yr%no_nit_yr
    print *, 'no_denit_yr ', ctem_yr%no_denit_yr
    print *, 'no_nitdenit_yr ', ctem_yr%no_nitdenit_yr
    print *, 'n2o_nit_yr ', ctem_yr%n2o_nit_yr
    print *, 'n2o_denit_yr ', ctem_yr%n2o_denit_yr
    print *, 'n2o_nitdenit_yr ', ctem_yr%n2o_nitdenit_yr
    print *, 'n2_denit_yr ', ctem_yr%n2_denit_yr
    print *, 'nvol_yr ', ctem_yr%nvol_yr
    print *, 'nleach_yr ', ctem_yr%nleach_yr
    print *, 'appl_fert_yr ', ctem_yr%appl_fert_yr
    print *, 'ndep_nh4_yr ', ctem_yr%ndep_nh4_yr
    print *, 'ndep_no3_yr ', ctem_yr%ndep_no3_yr
    print *, 'ndemand_wp_npp_yr ', ctem_yr%ndemand_wp_npp_yr
    print *, 'nuptake_p_nh4_yr ', ctem_yr%nuptake_p_nh4_yr
    print *, 'nuptake_p_no3_yr ', ctem_yr%nuptake_p_no3_yr
    print *, 'nuptake_a_actl_nh4_yr ', ctem_yr%nuptake_a_actl_nh4_yr
    print *, 'nuptake_a_actl_no3_yr ', ctem_yr%nuptake_a_actl_no3_yr
    print *, 'nuptake_yr ', ctem_yr%nuptake_yr
    print *, 'nleafns2s_yr ', ctem_yr%nleafns2s_yr
    print *, 'nstemns2s_yr ', ctem_yr%nstemns2s_yr
    print *, 'nrootns2s_yr ', ctem_yr%nrootns2s_yr
    print *, 'nalloc_l_yr ', ctem_yr%nalloc_l_yr
    print *, 'nalloc_s_yr ', ctem_yr%nalloc_s_yr
    print *, 'nalloc_r_yr ', ctem_yr%nalloc_r_yr
    print *, 'nresorped_s_yr ', ctem_yr%nresorped_s_yr
    print *, 'nresorped_r_yr ', ctem_yr%nresorped_r_yr
    print *, 'nre_alloc_s2l_yr ', ctem_yr%nre_alloc_s2l_yr
    print *, 'nre_alloc_r2l_yr ', ctem_yr%nre_alloc_r2l_yr
    print *, 'nlitr_l_yr ', ctem_yr%nlitr_l_yr
    print *, 'nlitr_s_yr ', ctem_yr%nlitr_s_yr
    print *, 'nlitr_r_yr ', ctem_yr%nlitr_r_yr
    print *, 'nlitr_yr ', ctem_yr%nlitr_yr
    print *, 'gl2bl_grass_nflux_yr ', ctem_yr%gl2bl_grass_nflux_yr
    print *, 'c2n_l_yr ', ctem_yr%c2n_l_yr
    print *, 'c2n_s_yr ', ctem_yr%c2n_s_yr
    print *, 'c2n_r_yr ', ctem_yr%c2n_r_yr
    print *, 'c2n_wp_yr ', ctem_yr%c2n_wp_yr
    print *, 'c2n_litr_yr ', ctem_yr%c2n_litr_yr
    print *, 'c2n_humus_yr ', ctem_yr%c2n_humus_yr
    print *, 'nhumtrs_yr ', ctem_yr%nhumtrs_yr
    print *, 'nmineral_litr_yr ', ctem_yr%nmineral_litr_yr
    print *, 'nmineral_humus_yr ', ctem_yr%nmineral_humus_yr
    print *, 'netnmineral_yr ', ctem_yr%netnmineral_yr
    print *, 'nimmobil_nh4_yr ', ctem_yr%nimmobil_nh4_yr
    print *, 'nimmobil_no3_yr ', ctem_yr%nimmobil_no3_yr
    print *, 'fNnetland_yr ', ctem_yr%fNnetland_yr

    print *, 'n2o_nitdenit_yr_g ', ctem_grd_yr%n2o_nitdenit_yr_g
    print *, 'n2_denit_yr_g ', ctem_grd_yr%n2_denit_yr_g
    print *, 'nvol_yr_g ', ctem_grd_yr%nvol_yr_g
    print *, 'nleach_yr_g ', ctem_grd_yr%nleach_yr_g
    print *, 'appl_fert_yr_g ', ctem_grd_yr%appl_fert_yr_g
    print *, 'ndep_nh4_yr_g ', ctem_grd_yr%ndep_nh4_yr_g
    print *, 'ndep_no3_yr_g ', ctem_grd_yr%ndep_no3_yr_g
    print *, 'ndemand_wp_npp_yr_g ', ctem_grd_yr%ndemand_wp_npp_yr_g
    print *, 'nuptake_p_nh4_yr_g ', ctem_grd_yr%nuptake_p_nh4_yr_g
    print *, 'nuptake_p_no3_yr_g ', ctem_grd_yr%nuptake_p_no3_yr_g
    print *, 'nuptake_a_actl_nh4_yr_g ', ctem_grd_yr%nuptake_a_actl_nh4_yr_g
    print *, 'nuptake_a_actl_no3_yr_g ', ctem_grd_yr%nuptake_a_actl_no3_yr_g
    print *, 'nuptake_yr_g ', ctem_grd_yr%nuptake_yr_g
    print *, 'nleafns2s_yr_g ', ctem_grd_yr%nleafns2s_yr_g
    print *, 'nstemns2s_yr_g ', ctem_grd_yr%nstemns2s_yr_g
    print *, 'nrootns2s_yr_g ', ctem_grd_yr%nrootns2s_yr_g
    print *, 'nalloc_l_yr_g ', ctem_grd_yr%nalloc_l_yr_g
    print *, 'nalloc_s_yr_g ', ctem_grd_yr%nalloc_s_yr_g
    print *, 'nalloc_r_yr_g ', ctem_grd_yr%nalloc_r_yr_g
    print *, 'nresorped_s_yr_g ', ctem_grd_yr%nresorped_s_yr_g
    print *, 'nresorped_r_yr_g ', ctem_grd_yr%nresorped_r_yr_g
    print *, 'nre_alloc_s2l_yr_g ', ctem_grd_yr%nre_alloc_s2l_yr_g
    print *, 'nre_alloc_r2l_yr_g ', ctem_grd_yr%nre_alloc_r2l_yr_g
    print *, 'nlitr_l_yr_g ', ctem_grd_yr%nlitr_l_yr_g
    print *, 'nlitr_s_yr_g ', ctem_grd_yr%nlitr_s_yr_g
    print *, 'nlitr_r_yr_g ', ctem_grd_yr%nlitr_r_yr_g
    print *, 'nlitr_yr_g ', ctem_grd_yr%nlitr_yr_g
    print *, 'gl2bl_grass_nflux_yr_g ', ctem_grd_yr%gl2bl_grass_nflux_yr_g
    print *, 'c2n_l_yr_g ', ctem_grd_yr%c2n_l_yr_g
    print *, 'c2n_s_yr_g ', ctem_grd_yr%c2n_s_yr_g
    print *, 'c2n_r_yr_g ', ctem_grd_yr%c2n_r_yr_g
    print *, 'c2n_wp_yr_g ', ctem_grd_yr%c2n_wp_yr_g
    print *, 'c2n_litr_yr_g ', ctem_grd_yr%c2n_litr_yr_g
    print *, 'c2n_humus_yr_g ', ctem_grd_yr%c2n_humus_yr_g
    print *, 'nhumtrs_yr_g ', ctem_grd_yr%nhumtrs_yr_g
    print *, 'nmineral_litr_yr_g ', ctem_grd_yr%nmineral_litr_yr_g
    print *, 'nmineral_humus_yr_g ', ctem_grd_yr%nmineral_humus_yr_g
    print *, 'netnmineral_yr_g ', ctem_grd_yr%netnmineral_yr_g
    print *, 'nimmobil_nh4_yr_g ', ctem_grd_yr%nimmobil_nh4_yr_g
    print *, 'nimmobil_no3_yr_g ', ctem_grd_yr%nimmobil_no3_yr_g
    print *, 'fNnetland_yr_g ', ctem_grd_yr%fNnetland_yr_g

    print *, 'laimaxg_yr_t ', ctem_tile_yr%laimaxg_yr_t
    print *, 'gleafmas_yr_t ', ctem_tile_yr%gleafmas_yr_t
    print *, 'gleafmas_NS_yr_t ', ctem_tile_yr%gleafmas_NS_yr_t
    print *, 'gleafmass_yr_t ', ctem_tile_yr%gleafmass_yr_t
    print *, 'bleafmas_yr_t ', ctem_tile_yr%bleafmas_yr_t
    print *, 'stemmass_yr_t ', ctem_tile_yr%stemmass_yr_t
    print *, 'stemmass_NS_yr_t ', ctem_tile_yr%stemmass_NS_yr_t
    print *, 'stemmasss_yr_t ', ctem_tile_yr%stemmasss_yr_t
    print *, 'rootmass_yr_t ', ctem_tile_yr%rootmass_yr_t
    print *, 'rootmass_NS_yr_t ', ctem_tile_yr%rootmass_NS_yr_t
    print *, 'rootmasss_yr_t ', ctem_tile_yr%rootmasss_yr_t
    print *, 'npp_yr_t ', ctem_tile_yr%npp_yr_t
    print *, 'gpp_yr_t ', ctem_tile_yr%gpp_yr_t
    print *, 'vcmax0_yr_t ', ctem_tile_yr%vcmax0_yr_t
    print *, 'leafns2s_yr_t ', ctem_tile_yr%leafns2s_yr_t
    print *, 'stemns2s_yr_t ', ctem_tile_yr%stemns2s_yr_t
    print *, 'rootns2s_yr_t ', ctem_tile_yr%rootns2s_yr_t
    print *, 're_alloc_s2l_yr_t ', ctem_tile_yr%re_alloc_s2l_yr_t
    print *, 're_alloc_r2l_yr_t ', ctem_tile_yr%re_alloc_r2l_yr_t
    print *, 're_alloc_sr2l_yr_t ', ctem_tile_yr%re_alloc_sr2l_yr_t
    print *, 'vgbiomas_yr_t ', ctem_tile_yr%vgbiomas_yr_t
    print *, 'autores_yr_t ', ctem_tile_yr%autores_yr_t
    print *, 'rmrveg_yr_t ', ctem_tile_yr%rmrveg_yr_t
    print *, 'totcmass_yr_t ', ctem_tile_yr%totcmass_yr_t
    print *, 'nh4_mass_yr_t ', ctem_tile_yr%nh4_mass_yr_t
    print *, 'no3_mass_yr_t ', ctem_tile_yr%no3_mass_yr_t
    print *, 'ngleafmas_yr_t ', ctem_tile_yr%ngleafmas_yr_t
    print *, 'ngleafmas_NS_yr_t ', ctem_tile_yr%ngleafmas_NS_yr_t
    print *, 'ngleafmass_yr_t ', ctem_tile_yr%ngleafmass_yr_t
    print *, 'nbleafmas_yr_t ', ctem_tile_yr%nbleafmas_yr_t
    print *, 'nstemmass_yr_t ', ctem_tile_yr%nstemmass_yr_t
    print *, 'nstemmass_NS_yr_t ', ctem_tile_yr%nstemmass_NS_yr_t
    print *, 'nstemmasss_yr_t ', ctem_tile_yr%nstemmasss_yr_t
    print *, 'nrootmass_yr_t ', ctem_tile_yr%nrootmass_yr_t
    print *, 'nrootmass_NS_yr_t ', ctem_tile_yr%nrootmass_NS_yr_t
    print *, 'nrootmasss_yr_t ', ctem_tile_yr%nrootmasss_yr_t
    print *, 'nlitrmass_yr_t ', ctem_tile_yr%nlitrmass_yr_t
    print *, 'soilnmas_yr_t ', ctem_tile_yr%soilnmas_yr_t
    print *, 'nvgbiomas_yr_t ', ctem_tile_yr%nvgbiomas_yr_t
    print *, 'nep_yr_t ', ctem_tile_yr%nep_yr_t
    print *, 'litrfall_yr_t ', ctem_tile_yr%litrfall_yr_t
    print *, 'hetrores_yr_t ', ctem_tile_yr%hetrores_yr_t
    print *, 'nbp_yr_t ', ctem_tile_yr%nbp_yr_t
    print *, 'emit_co2_yr_t ', ctem_tile_yr%emit_co2_yr_t
    print *, 'emit_co_yr_t ', ctem_tile_yr%emit_co_yr_t
    print *, 'emit_ch4_yr_t ', ctem_tile_yr%emit_ch4_yr_t
    print *, 'emit_nmhc_yr_t ', ctem_tile_yr%emit_nmhc_yr_t
    print *, 'emit_h2_yr_t ', ctem_tile_yr%emit_h2_yr_t
    print *, 'emit_nox_yr_t ', ctem_tile_yr%emit_nox_yr_t
    print *, 'emit_n2o_yr_t ', ctem_tile_yr%emit_n2o_yr_t
    print *, 'emit_pm25_yr_t ', ctem_tile_yr%emit_pm25_yr_t
    print *, 'emit_tpm_yr_t ', ctem_tile_yr%emit_tpm_yr_t
    print *, 'emit_tc_yr_t ', ctem_tile_yr%emit_tc_yr_t
    print *, 'emit_oc_yr_t ', ctem_tile_yr%emit_oc_yr_t
    print *, 'emit_bc_yr_t ', ctem_tile_yr%emit_bc_yr_t
    print *, 'burnfrac_yr_t ', ctem_tile_yr%burnfrac_yr_t
    print *, 'smfuncveg_yr_t ', ctem_tile_yr%smfuncveg_yr_t
    print *, 'bterm_yr_t ', ctem_tile_yr%bterm_yr_t
    print *, 'luc_emc_yr_t ', ctem_tile_yr%luc_emc_yr_t
    print *, 'lterm_yr_t ', ctem_tile_yr%lterm_yr_t
    print *, 'lucsocin_yr_t ', ctem_tile_yr%lucsocin_yr_t
    print *, 'mterm_yr_t ', ctem_tile_yr%mterm_yr_t
    print *, 'lucltrin_yr_t ', ctem_tile_yr%lucltrin_yr_t
    print *, 'ch4WetSpec_yr_t ', ctem_tile_yr%ch4WetSpec_yr_t
    print *, 'wetfdyn_yr_t ', ctem_tile_yr%wetfdyn_yr_t
    print *, 'ch4WetDyn_yr_t ', ctem_tile_yr%ch4WetDyn_yr_t
    print *, 'ch4soills_yr_t ', ctem_tile_yr%ch4soills_yr_t
    print *, 'veghght_yr_t ', ctem_tile_yr%veghght_yr_t
    print *, 'peatdep_yr_t ', ctem_tile_yr%peatdep_yr_t
    print *, 'peatSoilC_yr_t ', ctem_tile_yr%peatSoilC_yr_t
    print *, 'fProductDecomp_yr_t ', ctem_tile_yr%fProductDecomp_yr_t
    print *, 'litrmass_yr_t ', ctem_tile_yr%litrmass_yr_t
    print *, 'soilcmas_yr_t ', ctem_tile_yr%soilcmas_yr_t
    print *, 'litres_yr_t ', ctem_tile_yr%litres_yr_t
    print *, 'soilcres_yr_t ', ctem_tile_yr%soilcres_yr_t
    print *, 'bnf_tot_yr_t ', ctem_tile_yr%bnf_tot_yr_t
    print *, 'bnf_free_yr_t ', ctem_tile_yr%bnf_free_yr_t
    print *, 'bnf_ant_yr_t ', ctem_tile_yr%bnf_ant_yr_t
    print *, 'bnf_nat_yr_t ', ctem_tile_yr%bnf_nat_yr_t
    print *, 'nstress_yr_t ', ctem_tile_yr%nstress_yr_t
    print *, 'nitrif_yr_t ', ctem_tile_yr%nitrif_yr_t
    print *, 'no_nit_yr_t ', ctem_tile_yr%no_nit_yr_t
    print *, 'no_denit_yr_t ', ctem_tile_yr%no_denit_yr_t
    print *, 'no_nitdenit_yr_t ', ctem_tile_yr%no_nitdenit_yr_t
    print *, 'n2o_nit_yr_t ', ctem_tile_yr%n2o_nit_yr_t
    print *, 'n2o_denit_yr_t ', ctem_tile_yr%n2o_denit_yr_t
    print *, 'n2o_nitdenit_yr_t ', ctem_tile_yr%n2o_nitdenit_yr_t
    print *, 'n2_denit_yr_t ', ctem_tile_yr%n2_denit_yr_t
    print *, 'nvol_yr_t ', ctem_tile_yr%nvol_yr_t
    print *, 'nleach_yr_t ', ctem_tile_yr%nleach_yr_t
    print *, 'appl_fert_yr_t ', ctem_tile_yr%appl_fert_yr_t
    print *, 'ndep_nh4_yr_t ', ctem_tile_yr%ndep_nh4_yr_t
    print *, 'ndep_no3_yr_t ', ctem_tile_yr%ndep_no3_yr_t
    print *, 'ndemand_wp_npp_yr_t ', ctem_tile_yr%ndemand_wp_npp_yr_t
    print *, 'nuptake_p_nh4_yr_t ', ctem_tile_yr%nuptake_p_nh4_yr_t
    print *, 'nuptake_p_no3_yr_t ', ctem_tile_yr%nuptake_p_no3_yr_t
    print *, 'nuptake_a_actl_nh4_yr_t ', ctem_tile_yr%nuptake_a_actl_nh4_yr_t
    print *, 'nuptake_a_actl_no3_yr_t ', ctem_tile_yr%nuptake_a_actl_no3_yr_t
    print *, 'nuptake_yr_t ', ctem_tile_yr%nuptake_yr_t
    print *, 'nleafns2s_yr_t ', ctem_tile_yr%nleafns2s_yr_t
    print *, 'nstemns2s_yr_t ', ctem_tile_yr%nstemns2s_yr_t
    print *, 'nrootns2s_yr_t ', ctem_tile_yr%nrootns2s_yr_t
    print *, 'nalloc_l_yr_t ', ctem_tile_yr%nalloc_l_yr_t
    print *, 'nalloc_s_yr_t ', ctem_tile_yr%nalloc_s_yr_t
    print *, 'nalloc_r_yr_t ', ctem_tile_yr%nalloc_r_yr_t
    print *, 'nresorped_s_yr_t ', ctem_tile_yr%nresorped_s_yr_t
    print *, 'nresorped_r_yr_t ', ctem_tile_yr%nresorped_r_yr_t
    print *, 'nre_alloc_s2l_yr_t ', ctem_tile_yr%nre_alloc_s2l_yr_t
    print *, 'nre_alloc_r2l_yr_t ', ctem_tile_yr%nre_alloc_r2l_yr_t
    print *, 'nlitr_l_yr_t ', ctem_tile_yr%nlitr_l_yr_t
    print *, 'nlitr_s_yr_t ', ctem_tile_yr%nlitr_s_yr_t
    print *, 'nlitr_r_yr_t ', ctem_tile_yr%nlitr_r_yr_t
    print *, 'nlitr_yr_t ', ctem_tile_yr%nlitr_yr_t
    print *, 'gl2bl_grass_nflux_yr_t ', ctem_tile_yr%gl2bl_grass_nflux_yr_t
    print *, 'c2n_l_yr_t ', ctem_tile_yr%c2n_l_yr_t
    print *, 'c2n_s_yr_t ', ctem_tile_yr%c2n_s_yr_t
    print *, 'c2n_r_yr_t ', ctem_tile_yr%c2n_r_yr_t
    print *, 'c2n_wp_yr_t ', ctem_tile_yr%c2n_wp_yr_t
    print *, 'c2n_litr_yr_t ', ctem_tile_yr%c2n_litr_yr_t
    print *, 'c2n_humus_yr_t ', ctem_tile_yr%c2n_humus_yr_t
    print *, 'nhumtrs_yr_t ', ctem_tile_yr%nhumtrs_yr_t
    print *, 'nmineral_litr_yr_t ', ctem_tile_yr%nmineral_litr_yr_t
    print *, 'nmineral_humus_yr_t ', ctem_tile_yr%nmineral_humus_yr_t
    print *, 'netnmineral_yr_t ', ctem_tile_yr%netnmineral_yr_t
    print *, 'nimmobil_nh4_yr_t ', ctem_tile_yr%nimmobil_nh4_yr_t
    print *, 'nimmobil_no3_yr_t ', ctem_tile_yr%nimmobil_no3_yr_t
    print *, 'fNnetland_yr_t ', ctem_tile_yr%fNnetland_yr_t
    print *, 'end ctem dump'

  end subroutine ctemdump
  !! @}
  !=================================================================================

  !> \namespace ctemstatevars
  !> Contains the biogeochemistry-related variable type structures.
  !! 1. c_switch - switches for running CTEM, read from the joboptions file
  !! 2. vrot - CTEM's 'rot' vars
  !! 3. vgat - CTEM's 'gat' vars
  !! 4. ctem_grd - CTEM's grid average variables
  !! 5. ctem_tile - CTEM's variables per tile
  !! 6. ctem_mo - CTEM's variables monthly averaged (per pft)
  !! 7. ctem_grd_mo - CTEM's grid average monthly values
  !! 8. ctem_tile_mo - CTEM's variables per tile monthly values
  !! 9. ctem_yr - CTEM's average annual values (per PFT)
  !! 10. ctem_grd_yr - CTEM's grid average annual values
  !! 11. ctem_tile_yr - CTEM's variables per tile annual values
  !
end module ctemStateVars
