!> \file
!> Main model driver for CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing.
!! @author D. Verseghy, V. Arora, J. Melton, E. Chan
!! This driver program initializes the run, reads in CLASSIC input files,
!! calls the time-stepping subroutine, calls subroutines
!! that aggregate and write outputs, and closes the run for this grid cell.

module main

  implicit none

  public :: main_driver

contains

  !> \ingroup main_main_driver
  !> @{
  !>
  subroutine main_driver (longitude, latitude, lonIndex, latIndex, lonLocalIndex, latLocalIndex)

    use mainCore,           only : main_core_driver

    use classicParams,      only : nlat, nmos, ilg, nmon, monthend, DELT, zbldJobOpt, zrfhJobOpt, zrfmJobOpt
    use landuseChange,      only : initializeLandCover
    use ctemStateVars,      only : vrot, vgat, c_switch, initRowVarsBioGeoChem, resetMonthEnd, &
                                   resetYearEnd, resetMosaicAccum
    use classStateVars,     only : class_gat, class_rot, resetAccVars, resetClassMon, resetClassYr,initRowVarsPhysics
    use classGatherScatter, only : classGatherPrep
    use prepareOutputs,     only : class_monthly_aw, ctem_annual_aw, ctem_monthly_aw, &
                                   ctem_daily_aw, class_annual_aw, class_hh_w, class_daily_aw, convertUnitsCTEM
    use modelStateDrivers,  only : read_initialstate, write_restart, deallocInput, getMet, updateMet, getInput, updateInput
    use generalUtils,       only : run_model, findPermafrostVars, initRandomSeed, checksumCalc, makeDOYArray, &
                                   findDaylength, findLeapYears, findCloudiness
    use metDisaggModule,    only : disaggMet
    use outputManager,      only : consecDays

    implicit none

    ! Arguments
    real, intent(in) :: longitude, latitude                 !< Longitude/latitude of grid cell (degrees)
    integer, intent(in) :: lonIndex, latIndex               !< Index of grid cell being run on the input files grid
    integer, intent(in) :: lonLocalIndex, latLocalIndex     !< Index of grid cell being run on the output files grid

    ! Local variables
    integer :: lastDOY             !< Initialized to 365 days, can be overwritten later is leap = true and it is a leap year.
    integer :: metTimeIndex        !< Counter used to move through the meteorological input arrays
    logical :: metDone             !< Logical switch when the end of the stored meteorological array is reached.
    logical :: leapnow             !< Leap year flag (if the switch 'leap' is true, this will be used, otherwise it remains false)
    integer :: runyr               !< Year of the model run (counts up starting with readMetStartYear continously, even if metLoop > 1)

    integer :: N       !< Timestep counter
    integer :: NCOUNT  !< Counter for daily averaging
    integer :: NDAY    !< Number of short (physics) timesteps in one day. e.g., if physics timestep is 15 min this is 48.
    integer :: IYEAR   !< Year of run
    integer :: IMONTH  !< Month of the year simulation is in.
    integer :: IDAY    !< Julian day of the year
    integer :: IHOUR   !< Hour of day
    integer :: IMIN    !< Minutes elapsed in current hour
    integer :: DOM     !< Day of month counter
    integer :: NT      !<

    integer :: NML     !< Counter representing number of mosaic tiles on modelled domain that are land

    integer   :: lopcount, m, i ! xday, month1, month2

    integer, allocatable :: arrayDOYs(:) !< Array of the number of days in each year from met start to met end year.

    !    ================================================================================= 

    associate( &

    ! This driver is set up to handle one grid cell with any number of mosaic tiles.

    nltest            => nlat,                          & !< integer: Number of grid cells being modelled for this run
    nmtest            => nmos,                          & !< integer: Number of mosaic tiles per grid cell being modelled for this run

    ! Associate names with variables defined in derived types.

    readMetStartYear  => c_switch%readMetStartYear,     & !< integer: First year of meteorological forcing to read in from the met file 
    readMetEndYear    => c_switch%readMetEndYear,       & !< integer: Last year of meteorological forcing to read in from the met file 

    useTracer         => c_switch%useTracer,            & !< integer: Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
                                                          !!          useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple 
                                                          !!                        tracer then requires that the tracer values in the init_file and the
                                                          !!                        tracerCO2file are set to meaningful values for the experiment being run.
                                                          !!          useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                                          !!          useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme. 
    metLoop           => c_switch%metLoop,              & !< integer: no. of times the .met file is to be read. this
                                                          !!          option is useful to see how ctem's c pools
                                                          !!          equilibrate when driven with same climate data
                                                          !!          over and over again. 
    jhhstd            => c_switch%jhhstd,               & !< integer: day of the year to start writing the half-hourly output 
    jhhendd           => c_switch%jhhendd,              & !< integer: day of the year to stop writing the half-hourly output 
    jdstd             => c_switch%jdstd,                & !< integer: day of the year to start writing the daily output 
    jdendd            => c_switch%jdendd,               & !< integer: day of the year to stop writing the daily output 
    jhhsty            => c_switch%jhhsty,               & !< integer: simulation year (runyr) to start writing the half-hourly output 
    jhhendy           => c_switch%jhhendy,              & !< integer: simulation year (runyr) to stop writing the half-hourly output 
    jdsty             => c_switch%jdsty,                & !< integer: simulation year (runyr) to start writing the daily output 
    jdendy            => c_switch%jdendy,               & !< integer: simulation year (runyr) to stop writing the daily output 
    jmosty            => c_switch%jmosty,               & !< integer: Year to start writing out the monthly output files. If you want to write monthly outputs right 
    fixedYearLUC      => c_switch%fixedYearLUC,         & !< integer: Set the year to use for land cover if lnduseon is false. If set to -9999,
                                                          !!          we use the PFT distribution found in the initialization file. Any other year
                                                          !!          we search for that year in the LUCFile 
    fixedYearFER      => c_switch%fixedYearFER,         & !< integer: Set the year to use for fertilizer if fertilizeron is false.
    fixedYearDEP      => c_switch%fixedYearDEP,         & !< integer, Set the year to use for deposition if depositionon is false.
    fixedYearOBSWETF  => c_switch%fixedYearOBSWETF,     & !< integer: set the year to use for observed wetland fraction if transientOBSWETF is false. 

    ! Model switches:
    isnoalb           => c_switch%isnoalb,              & !< integer: 
    ctem_on           => c_switch%ctem_on,              & !< logical:
    leap              => c_switch%leap,                 & !< logical:
    Ncycle_on         => c_switch%Ncycle_on,            & !< logical:
    dofire            => c_switch%dofire,               & !< logical: 
    lnduseon          => c_switch%lnduseon,             & !< logical: 
    fertilizeron      => c_switch%fertilizeron,         & !< logical:
    depositionon      => c_switch%depositionon,         & !< logical: 
    transientCO2      => c_switch%transientCO2,         & !< logical: 
    doMethane         => c_switch%doMethane,            & !< logical: 
    transientCH4      => c_switch%transientCH4,         & !< logical: 
    transientPOPD     => c_switch%transientPOPD,        & !< logical: 
    blackCdepon       => c_switch%blackCdepon,          & !< logical: If true, include consideration of black carbon deposition in snow processes (4-band scheme) 
    transientOBSWETF  => c_switch%transientOBSWETF,     & !< logical: 
    transientFER      => c_switch%transientFER,         & !< logical:
    transientDEP      => c_switch%transientDEP,         & !< logical:
    doAnnualOutput    => c_switch%doAnnualOutput,       & !< logical: 
    doMonthOutput     => c_switch%doMonthOutput,        & !< logical: 
    doDayOutput       => c_switch%doDayOutput,          & !< logical: 
    doHhOutput        => c_switch%doHhOutput,           & !< logical: 
    projectedGrid     => c_switch%projectedGrid,        & !< logical: True if you have a projected lon lat grid, false if not. 
                                                          !!          Projected grids can only have regions referenced by the indexes,  
                                                          !!          not coordinates, when running a sub-region 
    dynamicTilingOn     => c_switch%dynamicTilingOn,        & !< logical:
    timberHarvest     => c_switch%timberHarvest,        & !< logical:
    prescribedFire      => c_switch%prescribedFire,      & !< logical:

    ! ROW:
    ipeatlandrow     => vrot%ipeatland,                 & !< integer, dimension(:,:) : This is first set in read_from_ctm. 
    daylrow          => vrot%dayl,                      & !< real, dimension(:) : 
    pfcancmxrow      => vrot%pfcancmx,                  & !< real, dimension(:,:,:) : 
    nfcancmxrow      => vrot%nfcancmx,                  & !< real, dimension(:,:,:) : 
    pddrow           => vrot%pdd,                       & !< real, dimension(:,:) : 

    latIndexROW      => class_rot%latIndexROW,          & !< real, dimension(:) : Index of grid cell being run on the input files grid (latitude)
    lonIndexROW      => class_rot%lonIndexROW,          & !< real, dimension(:) : Index of grid cell being run on the input files grid (longitude)
    DLONROW          => class_rot%DLONROW,              & !< real, dimension(:) : 
    DLATROW          => class_rot%DLATROW,              & !< real, dimension(:) : 
    RADJROW          => class_rot%RADJROW,              & !< real, dimension(:) : Latitude of grid cell (positive north of equator) [rad] 
    ZBLDROW          => class_rot%ZBLDROW,              & !< real, dimension(:) : 
    ZDHROW           => class_rot%ZDHROW,               & !< real, dimension(:) : 
    ZDMROW           => class_rot%ZDMROW,               & !< real, dimension(:) : 
    ZRFHROW          => class_rot%ZRFHROW,              & !< real, dimension(:) : 
    ZRFMROW          => class_rot%ZRFMROW,              & !< real, dimension(:) : 
    FAREROT          => class_rot%FAREROT,              & !< real, dimension(:,:) : 
    TCANROT          => class_rot%TCANROT,              & !< real, dimension(:,:) : Vegetation canopy temperature [K]
    TAROW            => class_rot%TAROW,                & !< real, dimension(:) : Air temperature at reference height [K] 
    BCSNROW          => class_rot%BCSNROW,              & !< real, dimension(:) : Black carbon mixing ratio [kg m-3] 

    ! GAT:
    ILMOS             => class_gat%ILMOS,               & !< integer, dimension(:) : Index of gridcell corresponding to current element of gathered vector of land surface variables [ ] 
    JLMOS             => class_gat%JLMOS,               & !< integer, dimension(:) : Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ] 
    grclarea          => vgat%grclarea,                 & !< real, dimension(:) : 
    wetfrac_presgat   => vgat%wetfrac_pres,             & !< real, dimension(:) : 

    !Disturbance forcing variables
    timharvrow        => vrot%timharvrow,               & !< real, dimension(:) :
    prsfirerow        => vrot%prsfirerow                & !< real(nlat,nmos)
    )

    !    =================================================================================

    ! Variables set in this section are stored via derived types defined in modules, Although they are not used
    ! explicitly in this routine, they need to initialized before being used in other routines.

    DLATROW = latitude
    DLONROW = longitude

    !! Store the index of the cell so that it can be used for a more informative error message
    !! during projected grid simulations. It is much easier to find the problematic cells this way.
    latIndexROW = latIndex
    lonIndexROW = lonIndex

    !> The grid-average height for the momentum diagnostic variables, ZDMROW, and for the
    !! energy diagnostic variables, ZDHROW, are hard-coded to the standard anemometer
    !! height of 10 m and to the screen height of 2 m respectively.
    ZDMROW(:) = 10.0
    ZDHROW(:) = 2.0

    !> ZRFMROW and ZRFHROW, the reference heights at which the momentum variables (wind speed) and energy variables
    !! (temperature and specific humidity) are provided.  In a run using atmospheric model forcing data, these heights
    !! would vary by time step, but since this version of the driver is set up to use field data, ZRFMROW and ZRFHROW
    !! refer to the measurement height of these variables, which is fixed. The value is read in from the job options file.
    ZRFMROW = zrfmJobOpt
    ZRFHROW = zrfhJobOpt

    !> ZBLDROW, the atmospheric blending height.  Technically this variable depends on the length scale of the
    !! patches of roughness elements on the land surface, but this is difficult to ascertain.  Usually it is assigned a value of 50 m.
    !!  The value is read in from the job options file.
    ZBLDROW = zbldJobOpt

    !> The timestep counter N for the run is initialized to 0, the daily
    !! averaging counter NCOUNT is set to 1, and the total number of
    !! timesteps in the day NDAY is calculated as the number of seconds
    !! in a day (86400) divided by the timestep length DELT.
    N = 0
    NCOUNT = 1
    NDAY = 86400/NINT(DELT)
    metTimeIndex = 1    !< Counter used to move through the meteorological input arrays
    metDone = .false.   !< Logical switch when the end of the stored meteorological array is reached.
    run_model = .true.  !< Simple logical switch to either keep run going or finish
    IMONTH = 0          !< Month of the year simulation is in.
    DOM = 1             !< Day of month counter
    iyear = 0           !< Year of simulation (based upon incoming meteorology), set to 0 to start.
    lopcount = 1
    lastDOY = 365
    wetfrac_presgat = -9999. !< If transientOBSWETF or fixedYearOBSWETF != -9999 this variable will be overwritten with
                              !! real wetland fractions. Otherwise the negative is used as a switch so the dynamic
                              !! wetland extent is used instead of the prescribed.


    !> Reset variables in preparation for the run
    call initRowVarsBioGeoChem
    call initRowVarsPhysics
    call resetAccVars(nlat, nmos)
    call resetMosaicAccum
    BCSNROW(:)=0.0
   
    !> Read in the model initial state. 
    call read_initialstate(lonIndex, latIndex)

    !> If using the four band albedo parameterization, read in the look up table for the 
    !! parameters. If black carbon deposition is desired, read in the BC deposition flux
    if (isnoalb == 1) then
      call getInput('FourBandAlbedo')
      if (blackCdepon) call getInput('BlackCarbon',longitude,latitude)
    end if
    
    if (ctem_on) then
      !> Read in the inputs for a run with biogeochemical component turned on
      call getInput('CO2') ! CO2 atmospheric concentration
      if (doMethane) call getInput('CH4') ! CH4 atmospheric concentration
      if (useTracer > 0) call getInput('tracerCO2',longitude,latitude) ! tracer atmospheric values
      if (.not. projectedGrid) then
        ! regular lon/lat grid
        if (dofire) call getInput('POPD',longitude,latitude) ! Population density
        if (dofire) call getInput('LGHT',longitude,latitude) ! Cloud-to-ground lightning frequency
        if (doMethane .and. transientOBSWETF .or. fixedYearOBSWETF /= - 9999) call getInput('OBSWETF',longitude,latitude) ! Observed wetland distribution
        if (lnduseon .or. (fixedYearLUC /= - 9999)) call getInput('LUC',longitude,latitude) ! Land use change
        if (Ncycle_on) then
          if (fertilizeron .or. (fixedYearFER .ne. -9999)) call getInput('FER',longitude,latitude) ! Fertilizer
          if (depositionon .or. (fixedYearDEP .ne. -9999)) call getInput('DEP',longitude,latitude) ! Deposition
        end if
      else
        ! Projected grids use the lon and lat indexes, not the actual coordinates
        if (dofire) call getInput('POPD',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Population density
        if (dofire) call getInput('LGHT',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Cloud-to-ground lightning frequency
        if (doMethane .and. transientOBSWETF .or. fixedYearOBSWETF /= - 9999) &
            call getInput('OBSWETF',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Observed wetland distribution
        if (lnduseon .or. (fixedYearLUC /= - 9999)) &
            call getInput('LUC',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Land use change
        if (Ncycle_on) then
          if (fertilizeron .or. (fixedYearFER .ne. -9999)) call getInput('FER',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Fertilizer
          if (depositionon .or. (fixedYearDEP .ne. -9999)) call getInput('DEP',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Deposition
        end if
      end if
      !> Regardless of whether lnduseon or not, we need to check the land cover that was read in
      !! and assign the CLASS PFTs as they are not read in when ctem_on.
      call initializeLandCover

    end if

    if(timberHarvest) call getInput('TIMHAR',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! annual fractional area harvested
    if(prescribedFire) call getInput('PRSFIRE',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! annual fractional area burned

    !> Read in the meteorological forcing data to a suite of arrays
    if (.not. projectedGrid) then
      ! regular lon lat grid
      call getMet(longitude, latitude, nday)
    else
      ! Projected grids use the lon and lat indexes, not the actual coordinates
      call getMet(longitude, latitude, nday, projLonInd = lonIndex, projLatInd = latIndex)
    end if

    !> In preparation for the use of the random number generator by disaggMet,
    !! we need to provide a seed to allow repeatable results.
    call initRandomSeed()

    !> Create an array of the year length in days for each year from 
    !! readMetStartYear to readMetEndYear. First, allocate the size of the array
    !! that will be filled with the last DOY information. If not leap, then just fill the 
    !! array with 365s.
    allocate(arrayDOYs(readMetEndYear - readMetStartYear + 1))
    if (leap) then
      call makeDOYArray(readMetStartYear, readMetEndYear,arrayDOYs)
    else
      arrayDOYs = 365
    end if 

    !> Now disaggregate the meteorological forcing to the right timestep
    !! for this model run (if needed; this is checked for in the subroutine)
    call disaggMet(longitude, latitude, arrayDOYs)

    ! Initialize accumulated array for monthly & yearly outputs
    call resetClassMon(nltest)
    call resetClassYr(nltest)
    if (ctem_on) then
      call resetMonthEnd(nltest, nmtest)
      call resetYearEnd(nltest, nmtest)
    end if

    !     **** LAUNCH RUN. **** ! EC - moved here from main_core_driver 

    !> The do while loop marks the beginning of the time stepping loop
    !! for the actual run.  N is incremented by 1, and the atmospheric forcing
    !! data for the current time step are updated for each grid cell or modelled
    !! area (see the manual section on “Data Requirements”).

    runyr = readMetStartYear  ! Initialize the runyr as the first year of met forcing.

    ! start up the main model loop

    mainModelLoop: do while (run_model)

      leapnow = .false.
      ! Check if this year is a leap year, and if so adjust the monthdays, monthend and mmday values.
      if (leap) call findLeapYears(iyear,leapnow,lastDOY)

      !
      !> Update the meteorological forcing data for current time step
      !
      call updateMet(metTimeIndex, lastDOY, iyear, iday, ihour, imin, metDone)

      ! Check if we are on the first timestep of the day
      if (ihour == 0 .and. imin == 0) then
  
        ! Find the daylength of this day
        daylrow = findDaylength(real(iday),radjrow(1)) ! following rest of code, radjrow is always given index of 1 offline.
        
        ! Update the lightning if fire is on and transientLGHT is true
        if (dofire .and. ctem_on) call updateInput('LGHT',runyr,imonth = imonth,iday = iday,dom = DOM)
  
        ! Update the wetland fractions if we are using read-in wetland fractions
        if (ctem_on .and. doMethane .and. (transientOBSWETF .or. fixedYearOBSWETF /= - 9999)) then
          call updateInput('OBSWETF',runyr,imonth = imonth,iday = iday,dom = DOM)
        end if
  
        ! Check if this is the first day of the year
        if (iday == 1) then

          ! Check if TCAN is coming into the model with an unreasonable value. If so, correct it
          ! with the air temp for the upcoming timestep.
          do i = 1,nlat ! loop 100
            do m = 1,nmos
              if ((TCANROT(i,m) > 373. .or. TCANROT(i,m) < 173.) .and. FAREROT(i,m) > 0.0) then
                print*,'Bad incoming TCAN value (',TCANROT(i,m),') automatically overwritten'
                print*,' with TA value:', TAROW(i)
                print*, 'lonIndex = ',lonIndex,' latIndex = ',latIndex
                TCANROT(i,m) = TAROW(i)
              end if 
            end do
          end do          
  
          ! If needed, update values that were read in from the accessory input files (popd, wetlands, lightning...)
          if (ctem_on) then
  
            if (transientCO2) call updateInput('CO2',runyr)
            if (useTracer > 0 .and. transientCO2) call updateInput('tracerCO2',runyr)
            if (doMethane .and. transientCH4) call updateInput('CH4',runyr)
            if (dofire .and. transientPOPD) call updateInput('POPD',runyr)
            if (lnduseon) then
              call updateInput('LUC',runyr)
            else ! If landuse change is not on, then set the next years landcover to be
              ! the same as this years.
              nfcancmxrow = pfcancmxrow
            end if
  
            if (Ncycle_on) then
                if (fertilizeron .and. transientFER) call updateInput('FER',runyr)
                if (depositionon .and. transientDEP) call updateInput('DEP',runyr)
            end if
  
            ! Reset the peatland degree days counter (JM:FLAG is this ok for S. Hemi to be Jan 1?)
            pddrow = 0

            !read in fire and harvest forcing if needed
            if (timberHarvest) then
              call updateInput('TIMHAR', runyr, imonth = imonth, iday = iday, dom = DOM)
            else
              timharvrow = 0.0
            end if

            if (prescribedFire) then
              call updateInput('PRSFIRE', runyr, imonth = imonth, iday = iday, dom = DOM)
            else
              prsfirerow = 0.0
            end if
  
          end if
        end if ! first day
      end if ! first timestep
  
      !> The cosine of the solar zenith angle COSZ is calculated from the day of
      !> the year, the hour, the minute and the latitude using basic radiation geometry,
      !> and (avoiding vanishingly small numbers) is assigned to CSZROW.  The fractional
      !> cloud cover FCLOROW is commonly not available so a rough estimate is
      !> obtained by setting it to 1 when precipitation is occurring, and to the fraction
      !> of incoming diffuse radiation XDIFFUS otherwise (assumed to be 1 when the sun
      !> is at the horizon, and 0.10 when it is at the zenith). These calculations are
      !> done in findCloudiness
  
      call findCloudiness(nltest, imin, ihour, iday, lastDOY)

      !> classGatherPrep assigns values to vectors governing the gather-scatter operations
      call classGatherPrep (ILMOS, JLMOS, FAREROT, & ! Formerly GATPREP
                            NML, 1, NMTEST, NMTEST, NLAT, ILG, 1, NLTEST)

      !print*, 'year=', iyear, 'day=', iday, ' hour=', ihour, ' min=', imin

      !> The call to main_core_driver is the same in both CLASSIC and CanAM.
      !> Changes to it and all common routines must be applicable to both models.
      call main_core_driver(nltest, NML, runyr, NDAY, ncount, IDAY, lastDOY, leapnow, N)

      if (ncount == nday) then

        DOM = DOM + 1 ! increment the day of month counter

        ! Reset mosaic accumulator arrays.
        if (ctem_on) call resetMosaicAccum

        ! Determine the active layer depth and depth to the frozen water table.
        ! This only occurs once per day since they don't change rapidly.
        call findPermafrostVars(nmtest, nltest, iday)

      end if ! ncount eq nday

      !=======================================================================
              
      ! Half-hourly physics outputs
      if  (doHhOutput .and. &
          (runyr >= jhhsty) .and. &
          (runyr <= jhhendy) .and. &
          (iday >= jhhstd) .and. &
          (iday <= jhhendd) ) call class_hh_w(lonLocalIndex, latLocalIndex, nltest, &
                                              nmtest, ncount, nday, iday, runyr)

      ! Daily physics outputs
      if (doDayOutput .and. &
          (runyr >= jdsty) .and. &
          (runyr <= jdendy) .and. &
          (iday  >= jdstd) .and. &
          (iday  <= jdendd))  call class_daily_aw(lonLocalIndex, latLocalIndex, iday, nltest, nmtest, &
                                                  ncount, nday, lastDOY, runyr)

      do NT = 1,NMON
        if ((IDAY == monthend(NT + 1)) .and. (NCOUNT == NDAY)) then
          IMONTH = NT
          DOM = 1 ! reset the day of month counter
        end if
      end do

      ! Monthly physics outputs
      if (doMonthOutput .and. (runyr >= jmosty)) call class_monthly_aw(lonLocalIndex, &
                                                        latLocalIndex, IDAY, runyr, NCOUNT, &
                                                        NDAY, nltest, nmtest, lastDOY)

      ! Annual physics outputs
      if (doAnnualOutput) call class_annual_aw(lonLocalIndex, latLocalIndex, IDAY, runyr, NCOUNT, NDAY, &
                                              nltest, nmtest, lastDOY)

      if (ctem_on .and. (ncount == nday)) then

        ! Convert units in preparation for output:
        call convertUnitsCTEM(nltest,nmtest)

        ! Daily outputs from biogeochem (CTEM)
        if (doDayOutput .and. &
            (runyr >= jdsty) .and. &
            (runyr <= jdendy) .and. &
            (iday   >= jdstd) .and. &
            (iday   <= jdendd)) call ctem_daily_aw(lonLocalIndex, latLocalIndex, nltest, &
                                                  nmtest, iday, ncount, nday, &
                                                  runyr, grclarea, ipeatlandrow)

        ! Monthly biogeochem outputs
        if (doMonthOutput .and. &
            (runyr >= jmosty)) call ctem_monthly_aw(lonLocalIndex, latLocalIndex, nltest, &
                                                    nmtest, iday, runyr, nday, lastDOY)

        ! Annual biogeochem outputs
        if (doAnnualOutput) call ctem_annual_aw(lonLocalIndex, latLocalIndex, iday, runyr, nltest, nmtest, lastDOY)
        
      end if

      ! Increment the consecDays if it is the end of the day. This  
      ! is used to create timestamps for the output files.
      if (NCOUNT == NDAY) consecDays = consecDays + 1.

      ! Check if it is the last timestep of the last day of the year
      if ((IDAY == lastDOY) .and. (NCOUNT == NDAY)) then

        write( * , * )'IYEAR = ',IYEAR,'runyr = ',runyr,'Loop count = ',lopcount,'/',metLoop

        if (dynamicTilingOn .or. (nmos > 3)) then
          ! If dynamic tiling is on or more than 3 tiles are used 
          ! only write to the restart file if it's the last timestep 
          ! as the file is large and computationally expensive to write
          if (lopcount == metLoop) then
            call write_restart(lonIndex,latIndex)
          end if 
        else
          ! If dynamic tiling is off the restart file is written more frequently. Note
          ! if you are finding your run is too computationally expensive, manually editing this
          ! so it is only written to at the end, like above, could help speed things up at the 
          ! risk that your run restart could be lost if the run fails midway (i.e. you would have to
          ! restart from the beginning, rather than some midpoint).
          call write_restart(lonIndex,latIndex)
        end if

        ! Increment the runyr
        runyr = runyr + 1

        ! Reset the month
        imonth = 0

      end if ! last day of year check

      ! Increment the counter for timestep of the day, or reset it to 1.
      NCOUNT = NCOUNT + 1
      if (NCOUNT > NDAY) then
        NCOUNT = 1
      end if

      !> Now check if the met file is done, needs to loop more, or just continues to the next timestep
      if (.not. metDone) then
        !> Increment the metTimeIndex to advance to the next timestep on the next time around
        metTimeIndex = metTimeIndex + 1
      else if (metDone) then
        !> End of met array read-in reached, decide what to do
        if (lopcount == metLoop) then
          !> The lopcount is reached so the run must be over
          run_model = .false.
        else
          !> Loop again so reset the metTimeIndex
          lopcount = lopcount + 1
          metTimeIndex = 1
        end if
      end if

    end do mainModelLoop ! MAIN MODEL LOOP

    ! If we've enabled checksums, now is the time to calculate them.
    if (c_switch%doChecksums) call checksumCalc(lonIndex, latIndex)

    ! deallocate arrays used for input files
    call deallocInput

    ! This is stored in outputManager so that means it retains its values between
    ! grid cells run. It must be reset here to ensure the value doesn't carry over to
    ! the next grid cell !
    consecDays = 1.

    deallocate(arrayDOYs)

    end associate
    return

  end subroutine main_driver
  !! @}
  !> \namespace main
  !> Main model driver for CLASSIC in stand-alone mode using specified boundary
  !! conditions and atmospheric forcing.
  !!
  !! This driver program initializes the run, reads in CLASSIC input files,
  !! calls the time-stepping subroutine (mainCore), calls subroutines
  !! that aggregate and write outputs, and closes the run for this grid cell.

end module main
