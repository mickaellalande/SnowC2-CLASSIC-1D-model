!> \file
!> Parses command line arguments to program and reads in joboptions file.
module readJobOpts

  implicit none

  public :: readFromJobOptions
  public :: readOptions
  public :: parsecoords

  !> Stores geographic information about the total domain and the region to be simulated
  type simulationDomain
    real, dimension(:), allocatable     :: lonLandCell, latLandCell     !< Long/Lat values of only the land cells in our model domain
    integer, dimension(:), allocatable  :: lonLandIndex, latLandIndex   !< Indexes of only the land cells in our model domain for our resolution
    real, dimension(:), allocatable     :: allLonValues, allLatValues   !< All long/Lat values in our model domain (including ocean/non-land)
    integer, dimension(:), allocatable  :: lonLocalIndex, latLocalIndex !< The index for only the region that is being simulated
    real, dimension(:), allocatable     :: latUnique, lonUnique         !< The index for only the region that is being simulated with each value only once
    integer                             :: LandCellCount                !< number of land cells that the model will run over
    real, dimension(4) :: domainBounds                                  !< Corners of the domain to be simulated (netcdfs)
    integer :: srtx                                                     !< starting index for this simulation for longitudes
    integer :: srty                                                     !< starting index for this simulation for latitudes
    integer :: cntx                                                     !< number of grid cells for this simulation in the longitude direction
    integer :: cnty                                                     !< number of grid cells for this simulation in the latitude direction
  end type

  type(simulationDomain) :: myDomain

contains

  ! ---------------------------------------------------
  !> \ingroup readjobopts_readFromJobOptions
  !! @{
  !> Opens the joboptions file, calls readOptions, and determines the geographic domain
  !! of the simulation. 
  !> @author Joe Melton, Ed Chan

  subroutine readFromJobOptions

    implicit none

    character(350) :: jobfile
    character(350) :: argbuff
    integer :: argcount, iargc

    !-------------------------
    ! read the joboptions

    argcount = iargc()

    if (argcount /= 2) then
      write( * , * )'Usage is as follows'
      write( * , * )' '
      write( * , * )'bin/CLASSIC joboptions_file longitude/{longitude}/latitude/{latitude}'
      write( * , * )' '
      write( * , * )' - joboptions_file - an example is '
      write( * , * )'  configurationFiles/template_job_options_file.txt.'
      write( * , * )' '
      write( * , * )' - longitude/latitude '
      write( * , * )'  e.g. 105.23/40.91 '
      write( * , * )' '
      write( * , * )' * OR * '
      write( * , * )' if you wish to run a region then you give '
      write( * , * )' the corners of the box you wish to run '
      write( * , * )' '
      write( * , * )' - longitude/longitude/latitude/latitude '
      write( * , * )'  e.g. 90/105/30/45'
      write( * , * )' '
      write( * , * )' ** If you are running a projected grid you must'
      write( * , * )' use the grid cell indices, not coordinates !** '
      write( * , * )' '
      stop
    end if

    !> Argument 1 is the jobfile, which is openned and the namelist is read
    call getarg(1,jobfile)

    open(10,file = jobfile,action = 'read',status = 'old')
    call readOptions(10)
    close(10)

    !> Parse the 2nd argument to get the domain that the simulation should be run over
    call getarg(2, argbuff)
    call parsecoords(argbuff, myDomain%domainBounds)

  end subroutine readFromJobOptions
  !! @}
  ! ----------------------------------------------------------------------------------

  !> \ingroup readOptions_readFromJobOptions
  !! @{
  !> Reads from the joboptions file and assigns the model switches. 
  !! All switches are described in the configurationFiles/template_job_options_file.txt file
  !! and also the user manual (housed in type ctem_switches in ctemStateVars.f90).
  !> @author Joe Melton, Ed Chan

  subroutine readOptions(input_unit)

    use ctemStateVars, only : c_switch
    use classicParams, only : runParamsFile, PFTCompetitionSwitch, &
                              zbldJobOpt, zrfhJobOpt, zrfmJobOpt
    implicit none

    integer :: input_unit 
    integer :: ioerr = 0
    real :: ZBLD, ZRFH, ZRFM

    ! -------------

    logical, pointer :: agcm_classic_on   !< Switch for running CLASSIC in CanESM
    logical, pointer :: projectedGrid  

    ! ctem model switches

    integer, pointer :: metLoop
    logical, pointer :: ctem_on
    logical, pointer :: Ncycle_on
    integer, pointer :: readMetStartYear
    integer, pointer :: readMetEndYear
    logical, pointer :: lnduseon
    logical, pointer :: fertilizeron
    logical, pointer :: depositionon
    integer, pointer :: spinfast
    integer, pointer :: useTracer
    character(:), pointer :: tracerCO2file
    logical, pointer :: transientCO2
    character(:), pointer :: CO2File
    integer, pointer :: fixedYearCO2
    logical, pointer :: useStaticPeatDep
    logical, pointer :: doMethane
    logical, pointer :: transientCH4
    character(:), pointer :: CH4File
    integer, pointer :: fixedYearCH4
    logical, pointer :: transientPOPD
    character(:), pointer :: POPDFile
    integer, pointer :: fixedYearPOPD
    logical, pointer :: dofire
    integer, pointer :: fixedYearLGHT
    logical, pointer :: transientLGHT
    character(:), pointer :: LGHTFile
    character(:), pointer :: LUCFile
    character(:), pointer :: FERFile
    character(:), pointer :: DEPFile
    integer, pointer :: fixedYearLUC
    integer, pointer :: fixedYearFER
    logical, pointer :: transientFER
    integer, pointer :: fixedYearDEP
    logical, pointer :: transientDEP
    logical, pointer :: transientOBSWETF
    character(:), pointer :: OBSWETFFile
    integer, pointer  :: fixedYearOBSWETF
    logical, pointer :: PFTCompetition
    logical, pointer :: inibioclim
    logical, pointer :: start_bare
    logical, pointer :: allLocalTime            !< If your gridded meteorology is relative to Greenwich then set this 
                                                !! to false. If the meteorology is all in local time set it to true. To
                                                !! determine which, look at your shortwave radiation. As you move through time 
                                                !! do you see the sun move across longitudes (allLocalTime = .false.) or across 
                                                !! latitudes (allLocalTime = .true.). It seems reanalysis will generally be false 
                                                !! while climate model outputs are generally true.
    character(:), pointer :: metFileFss
    character(:), pointer :: metFilefracFsf     !< character(350): File name for the diffuse fraction of total shortwave radiation input file
    character(:), pointer :: metFileFdl
    character(:), pointer :: metFilePre
    character(:), pointer :: metFileSnow   !< character(350): location of the snow meteorology file
    character(:), pointer :: metFileTa
    character(:), pointer :: metFileQa
    character(:), pointer :: metFileUv
    character(:), pointer :: metFilePres
    character(:), pointer :: init_file
    character(:), pointer :: rs_file_to_overwrite
    character(:), pointer :: runparams_file
    character(:), pointer :: output_directory
    character(:), pointer :: xmlFile
    logical, pointer :: leap
    logical, pointer :: dynamicTilingOn
    logical, pointer :: tileAgeReset
    logical, pointer :: trackTileAge 
    logical, pointer :: timberHarvest
    character(:), pointer :: timberHarvestFile
    logical, pointer :: prescribedFire
    character(:), pointer :: prescribedFireFile
    ! -------------
    ! class model switches

    integer, pointer :: idisp
    integer, pointer :: izref
    integer, pointer :: islfd
    integer, pointer :: ipcp
    integer, pointer :: iwf
    integer, pointer :: ITC
    integer, pointer :: ITCG
    integer, pointer :: ITG
    integer, pointer :: IPAI
    integer, pointer :: IHGT
    integer, pointer :: IALC
    integer, pointer :: IALS
    integer, pointer :: IALG
    integer, pointer :: isnoalb
    character(:), pointer :: alb4BandParamsFile
    character(:), pointer :: blackCarbonFile
    logical, pointer :: blackCdepon
    logical, pointer :: blackCtransientDep
    integer, pointer :: fixedYearBCDep
    logical, pointer :: KsatScalaron

    ! -------------
    ! Output switches

    integer, pointer :: jhhstd          !< day of the year to start writing the half-hourly output
    integer, pointer :: jhhendd         !< day of the year to stop writing the half-hourly output
    integer, pointer :: jdstd           !< day of the year to start writing the daily output
    integer, pointer :: jdendd          !< day of the year to stop writing the daily output
    integer, pointer :: jhhsty          !< simulation year (iyear) to start writing the half-hourly output
    integer, pointer :: jhhendy         !< simulation year (iyear) to stop writing the half-hourly output
    integer, pointer :: jdsty           !< simulation year (iyear) to start writing the daily output
    integer, pointer :: jdendy          !< simulation year (iyear) to stop writing the daily output
    integer, pointer :: jmosty          !< Year to start writing out the monthly output files. If you want to write monthly outputs right
    logical, pointer :: doperpftoutput  !< Switch for making extra output files that are at the per PFT level
    logical, pointer :: dopertileoutput !< Switch for making extra output files that are at the per tile level
    logical, pointer :: doAnnualOutput  !< Switch for making annual output files
    logical, pointer :: doMonthOutput   !< Switch for making monthly output files
    logical, pointer :: doDayOutput     !< Switch for making daily output files
    logical, pointer :: doHhOutput      !< Switch for making half hourly output files
    logical, pointer :: doChecksums
    character(:), pointer :: Comment    !< Comment about the run that will be written to the output netcdfs

    ! Order of the namelist and order in the file don't have to match.

    namelist /joboptions/ &
        projectedGrid, &
        metLoop, &
        readMetStartYear, &
        readMetEndYear, &
        leap, &
        ctem_on, &
        Ncycle_on, &
        dynamicTilingOn, &
        tileAgeReset, &
        trackTileAge, &
        timberHarvest, &
        timberHarvestFile, &
        prescribedFire, &
        prescribedFireFile, &
        spinfast, &
        useTracer, &
        tracerCO2file, &
        transientCO2, &
        CO2File, &
        fixedYearCO2, &
        useStaticPeatDep, &
        doMethane, &
        transientCH4, &
        CH4File, &
        fixedYearCH4, &
        transientPOPD, &
        POPDFile, &
        fixedYearPOPD, &
        lnduseon, &
        LUCFile, &
        fixedYearLUC, &
        fertilizeron, &
        FERFile,      &
        fixedYearFER, &
        transientFER, &
        depositionon, &
        fixedYearDEP, &
        DEPFile,      &
        transientDEP, &
        PFTCompetition, &
        inibioclim, &
        start_bare, &
        dofire, &
        transientLGHT, &
        fixedYearLGHT, &
        LGHTFile, &
        transientOBSWETF, &
        OBSWETFFile, &
        fixedYearOBSWETF, &
        allLocalTime, &
        metFileFss, &
        metFilefracFsf, &
        metFileFdl, &
        metFilePre, &
        metFileSnow, &
        metFileTa, &
        metFileQa, &
        metFileUv, &
        metFilePres, &
        init_file, &
        rs_file_to_overwrite, &
        runparams_file, &
        IDISP, &
        IZREF, &
        ZBLD, &
        ZRFH, &
        ZRFM, &
        ISLFD, &
        IPCP, &
        ITC, &
        ITCG, &
        ITG, &
        IWF, &
        IPAI, &
        IHGT, &
        IALC, &
        IALS, &
        IALG, &
        isnoalb, &
        alb4BandParamsFile, &
        blackCdepon, &
        blackCtransientDep, &
        blackCarbonFile, &
        fixedYearBCDep, &
        KsatScalaron, &
        output_directory, &
        xmlFile, &
        doperpftoutput, &
        dopertileoutput, &
        doHhOutput, &
        JHHSTD, &
        JHHENDD, &
        JHHSTY, &
        JHHENDY, &
        doDayOutput, &
        JDSTD, &
        JDENDD, &
        JDSTY, &
        JDENDY, &
        doMonthOutput, &
        JMOSTY, &
        doAnnualOutput, &
        doChecksums, &
        Comment

    namelist /canesm_options/ agcm_classic_on

    ! Point pointers:
    agcm_classic_on  => c_switch%agcm_classic_on
    projectedGrid   => c_switch%projectedGrid
    metLoop         => c_switch%metLoop
    readMetStartYear=> c_switch%readMetStartYear
    readMetEndYear  => c_switch%readMetEndYear
    ctem_on         => c_switch%ctem_on
    Ncycle_on       => c_switch%Ncycle_on
    dynamicTilingOn   => c_switch%dynamicTilingOn
    tileAgeReset    => c_switch%tileAgeReset
    trackTileAge    => c_switch%trackTileAge
    timberHarvest   => c_switch%timberHarvest
    timberHarvestFile => c_switch%timberHarvestFile
    prescribedFire  => c_switch%prescribedFire
    prescribedFireFile  => c_switch%prescribedFireFile
    fertilizeron    => c_switch%fertilizeron
    depositionon    => c_switch%depositionon
    lnduseon        => c_switch%lnduseon
    LUCFile         => c_switch%LUCFile
    FERFile         => c_switch%FERFile
    DEPFile         => c_switch%DEPFile
    fixedYearLUC    => c_switch%fixedYearLUC
    fixedYearFER    => c_switch%fixedYearFER
    fixedYearDEP    => c_switch%fixedYearDEP
    spinfast        => c_switch%spinfast
    useTracer       => c_switch%useTracer
    tracerCO2file   => c_switch%tracerCO2file
    transientCO2    => c_switch%transientCO2
    CO2File         => c_switch%CO2File
    fixedYearCO2    => c_switch%fixedYearCO2
    useStaticPeatDep=> c_switch%useStaticPeatDep
    doMethane       => c_switch%doMethane
    transientCH4    => c_switch%transientCH4
    CH4File         => c_switch%CH4File
    fixedYearCH4    => c_switch%fixedYearCH4
    transientPOPD   => c_switch%transientPOPD
    POPDFile        => c_switch%POPDFile
    fixedYearPOPD   => c_switch%fixedYearPOPD
    dofire          => c_switch%dofire
    fixedYearLGHT   => c_switch%fixedYearLGHT
    transientLGHT   => c_switch%transientLGHT
    LGHTFile        => c_switch%LGHTFile
    transientOBSWETF=> c_switch%transientOBSWETF
    transientFER    => c_switch%transientFER
    transientDEP    => c_switch%transientDEP
    OBSWETFFile     => c_switch%OBSWETFFile
    fixedYearOBSWETF=> c_switch%fixedYearOBSWETF
    PFTCompetition  => c_switch%PFTCompetition
    inibioclim      => c_switch%inibioclim
    start_bare      => c_switch%start_bare
    rs_file_to_overwrite => c_switch%rs_file_to_overwrite
    output_directory => c_switch%output_directory
    xmlFile         => c_switch%xmlFile
    allLocalTime    => c_switch%allLocalTime
    metFileFss      => c_switch%metFileFss
    metFilefracFsf  => c_switch%metFilefracFsf
    metFileFdl      => c_switch%metFileFdl
    metFilePre      => c_switch%metFilePre
    metFileSnow     => c_switch%metFileSnow
    metFileTa       => c_switch%metFileTa
    metFileQa       => c_switch%metFileQa
    metFileUv       => c_switch%metFileUv
    metFilePres     => c_switch%metFilePres
    runparams_file  => c_switch%runparams_file
    init_file       => c_switch%init_file
    IDISP           => c_switch%IDISP
    IZREF           => c_switch%IZREF
    ISLFD           => c_switch%ISLFD
    IPCP            => c_switch%IPCP
    ITC             => c_switch%ITC
    ITCG            => c_switch%ITCG
    ITG             => c_switch%ITG
    IWF             => c_switch%IWF
    IPAI            => c_switch%IPAI
    IHGT            => c_switch%IHGT
    IALC            => c_switch%IALC
    IALS            => c_switch%IALS
    IALG            => c_switch%IALG
    isnoalb         => c_switch%isnoalb
    alb4BandParamsFile => c_switch%alb4BandParamsFile
    blackCarbonFile => c_switch% blackCarbonFile
    blackCdepon     => c_switch%blackCdepon
    blackCtransientDep => c_switch%blackCtransientDep
    fixedYearBCDep => c_switch%fixedYearBCDep
    KsatScalaron    => c_switch%KsatScalaron
    leap            => c_switch%leap
    jhhstd          => c_switch%jhhstd
    jhhendd         => c_switch%jhhendd
    jdstd           => c_switch%jdstd
    jdendd          => c_switch%jdendd
    jhhsty          => c_switch%jhhsty
    jhhendy         => c_switch%jhhendy
    jdsty           => c_switch%jdsty
    jdendy          => c_switch%jdendy
    jmosty          => c_switch%jmosty
    doAnnualOutput  => c_switch%doAnnualOutput
    doperpftoutput  => c_switch%doperpftoutput
    dopertileoutput => c_switch%dopertileoutput
    doMonthOutput   => c_switch%doMonthOutput
    doDayOutput     => c_switch%doDayOutput
    doHhOutput      => c_switch%doHhOutput
    doChecksums     => c_switch%doChecksums
    Comment         => c_switch%Comment

    ! Read the job options.
    read(input_unit,nml = joboptions)

    ! Set the default canesm_options.
    agcm_classic_on=.false.

    ! Read the canesm_options namelist, if it is supplied.
    rewind(input_unit)
    read(input_unit,nml = canesm_options, iostat=ioerr)
    !if ( ioerr == 0 ) write(*,canesm_options)

    ! Assign some vars that are passed out
    runParamsFile = runparams_file
    PFTCompetitionSwitch = PFTCompetition
    zbldJobOpt = ZBLD
    zrfhJobOpt = ZRFH
    zrfmJobOpt = ZRFM

  end subroutine readOptions
  !! @}
  ! ----------------------------------------------------------------------------------

  !> \ingroup readjobopts_parsecoords
  !! @{
  !> Parses a coordinate string
  !> @author Joe Melton

  subroutine parsecoords (coordstring, val)

    implicit none

    character(45), intent(in)  :: coordstring
    real, dimension(4), intent(out) :: val

    character(10), dimension(4) :: cval = '0'

    integer :: i
    integer :: lasti = 1
    integer :: part  = 1

    do i = 1,len_trim(coordstring)
      if (coordstring(i:i) == '/') then
        cval(part) = coordstring(lasti:i - 1)
        lasti = i + 1
        part = part + 1
      end if
    end do

    cval(part) = coordstring(lasti:i - 1)

    read(cval, * )val

    if (part < 4) then
      val(3) = val(2)
      val(4) = val(3)
      val(2) = val(1)
    end if

  end subroutine parsecoords
  !! @}

  !> \namespace readjobopts
  !> Parses command line arguments to program and reads in joboptions file.

end module readJobOpts
