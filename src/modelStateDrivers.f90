!> \file
!> Central driver to read in, and write out all model state variables (replacing INI and CTM files)
!! as well as the model inputs such as MET, population density, land use change, CO2 etc.
!> @author J. Melton, A. Asaadi

module modelStateDrivers

  use fileIOModule
  use generalUtils, only : closeEnough, fndloc

  use, intrinsic :: iso_fortran_env, only: r8=>real64

  implicit none

  public  :: read_modelsetup
  public  :: read_initialstate
  public  :: write_restart
  public  :: getInput
  public  :: updateInput
  public  :: getMet
  public  :: updateMet
  public  :: checkTimeUnits
  public  :: deallocInput
  private :: closestCell

  ! Note: the inputs that are daily should be defined r8 for the their time vectors to ensure precision is kept.

  integer, dimension(:), allocatable :: CO2Time           !< The time (years) from the CO2File
  real, dimension(:), allocatable :: CO2FromFile          !< The array of CO2 values (ppm) from the CO2File
  integer, dimension(:), allocatable :: tracerCO2Time     !< The time (years) from the tracerCO2File
  real, dimension(:), allocatable :: tracerCO2FromFile    !< The array of tracerCO2 values (varied units) from the tracerCO2File
  integer, dimension(:), allocatable :: CH4Time           !< The time (years) from the CH4File
  real, dimension(:), allocatable :: CH4FromFile          !< The array of CH4 values (ppm) from the CH4File
  integer, dimension(:), allocatable :: POPDTime          !< The time (years) from the population density file
  real, dimension(:), allocatable :: POPDFromFile         !< The array of CH4 values (ppm) from the POPDFile
  real(r8), dimension(:), allocatable :: LGHTTime         !< The time from the lightning density file (usually months)
  real, dimension(:), allocatable :: LGHTFromFile         !< The array of lightning density from the LGHTFile
  real(r8), dimension(:), allocatable :: BCDepTime         !< The time from the lightning density file (usually months)
  real, dimension(:), allocatable :: BCDepFromFile         !< The array of lightning density from the LGHTFile
  integer, dimension(:), allocatable :: LUCTime           !< The time from the LUC file
  real, dimension(:,:), allocatable :: LUCFromFile        !< The array of LUC from the LUCFile
  integer, dimension(:), allocatable :: FERTime           !< The time from the FER file
  real, dimension(:), allocatable :: FERFromFile          !< The array of FER from the FERFile
  integer, dimension(:), allocatable :: DEPTime           !< The time from the DEP file
  real, dimension(:), allocatable :: DEPFromFile          !< The array of DEP from the DEPFile
  integer, dimension(:), allocatable :: TIMTime         !< The time from the timber harvest file
  real, dimension(:), allocatable :: TIMFromFile          !< The array of timber harvest from the timber harvest file
  integer, dimension(:), allocatable :: FireTime         !< The time from the prescribed fire file
  real, dimension(:), allocatable :: FireFromFile          !< The array of timber harvest from the prescribed fire file
  real(r8), dimension(:), allocatable :: OBSWETFTime      !< The time from the observed wetland distribution file
  real, dimension(:), allocatable :: OBSWETFFromFile      !< The array of observed wetland distribution from the OBSWETFFile

  real(r8), dimension(:), allocatable :: metTime          !< The time from the Met file
  real, dimension(:), allocatable :: metFss               !< Incoming total shortwave radiation from metFile \f$[W m^{-2} ]\f$
  real, dimension(:), allocatable :: metfracFsf           !< Incoming diffuse fraction of total shortwave radiation from metFile \f$[ ]\f$
  real, dimension(:), allocatable :: metFdl               !< Incoming longwave radiation from metFile \f$[W m^{-2} ]\f$
  real, dimension(:), allocatable :: metPre               !< Precipitation from metFile \f$[kg m^{-2} s^{-1} ]\f$
  real, dimension(:), allocatable :: metSnow              !< Precipitation that is snow from metFile \f$[kg m^{-2} s^{-1} ]\f$
  real, dimension(:), allocatable :: metFracSnow          !< Fraction of total precipitation that is snow from metFile \f$[ ]\f$
  real, dimension(:), allocatable :: metTa                !< Air temperature from metFile (Celsius)
  real, dimension(:), allocatable :: metQa                !< Specific humidity from metFile
  real, dimension(:), allocatable :: metUv                !< Wind speed from metFile
  real, dimension(:), allocatable :: metPres              !< Atmospheric pressure from metFile

  integer :: metFssId                             !> netcdf file id for the incoming total shortwave radiation meteorology file
  character(80) :: metFssVarName                  !> Name of variable in file
  integer :: metfracFsfId                         !> netcdf file id for the incoming diffuse fraction of shortwave radiation meteorology file
  character(80) :: metfracFsfVarName              !> Name of variable in file
  integer :: metFdlId                             !> netcdf file id for the incoming longwave radiation meteorology file
  character(80) :: metFdlVarName                  !> Name of variable in file
  integer :: metPreId                             !> netcdf file id for the precipitation meteorology file
  character(80) :: metPreVarName                  !> Name of variable in file
  integer :: metSnowId                            !> netcdf file id for the snow precipitation meteorology file
  character(80) :: metSnowVarName                 !> Name of variable in file
  integer :: metTaId                              !> netcdf file id for the air temperature meteorology file
  character(80) :: metTaVarName                   !> Name of variable in file
  integer :: metQaId                              !> netcdf file id for the specific humidity meteorology file
  character(80) :: metQaVarName                   !> Name of variable in file
  integer :: metUvId                              !> netcdf file id for the wind speed meteorology file
  character(80) :: metUvVarName                   !> Name of variable in file
  integer :: metPresId                            !> netcdf file id for the atmospheric pressure meteorology file
  character(80) :: metPresVarName                 !> Name of variable in file
  integer :: initid                               !> netcdf file id for the model initialization file
  integer :: rsid                                 !> netcdf file id for the model restart file
  integer :: co2id                                !> netcdf file id for the CO2 input file
  character(80) :: co2VarName                     !> Name of variable in file
  integer :: tracerco2id                          !> netcdf file id for the CO2 input file
  character(80) :: tracerco2VarName               !> Name of variable in file
  integer :: ch4id                                !> netcdf file id for the CH4 input file
  character(80) :: ch4VarName                     !> Name of variable in file
  integer :: popid                                !> netcdf file id for the population density input file
  character(80) :: popVarName                     !> Name of variable in file
  integer :: lghtid                               !> netcdf file id for the lightning density input file
  character(80) :: lghtVarName                    !> Name of variable in file
  integer :: lucid                                !> netcdf file id for the land use change input file
  character(80) :: lucVarName                     !> Name of variable in file
  integer :: ferid                                !> netcdf file id for the fertilizer input file
  character(80) :: ferVarName                     !> Name of variable in file
  integer :: depid                                !> netcdf file id for the deposition input file
  character(80) :: depVarName                     !> Name of variable in file
  integer :: timid                                !> netcdf file id for the timber harvest input file
  character(80) :: timVarName                     !> Name of variable in file
  integer :: fireid                               !> netcdf file id for the prescribed fire input file
  character(80) :: fireVarName                    !> Name of variable in file
  integer :: obswetid                             !> netcdf file id for the observed wetland distribution input file
  character(80) :: obswetVarName                  !> Name of variable in file
  integer :: albid                                !> netcdf file id for the 4 band albedo parameter look up table (only if isnoalb = 1)
  character(80) :: bcDepVarName                    !> Name of variable in file
  integer :: bcid                                 !> netcdf file id for the black carbon deposition flux file (only used if isnoabl = 1 & blackCdepon = true)

  real :: metInputTimeStep                        !> The timestep of the read in meteorology (hours)
  integer :: totlon,totlat                        !> Number of lons/lats in the init file.

contains

  !---

  !> \ingroup modelstatedrivers_read_modelsetup
  !! @{
  !> Reads in the model setup from the netcdf initialization file.
  !> The number of latitudes is always 1 offline while the maximum number of
  !> mosaics (nmos), the number of soil layers (ignd), are read from the netcdf.
  !> ilg is then calculated from nlat and nmos.
  !> @author Joe Melton

  subroutine read_modelsetup

    use ctemStateVars, only : c_switch
    use classicParams, only : nmos, testnmos, nlat, ignd, testignd, ilg, zero, ican, test2ican, icc, test2icc  ! These are set in this subroutine !
    use outputManager, only : writeStaticVars
    use readJobOpts,   only : myDomain

    implicit none

    ! Local vars
    real, allocatable, dimension(:,:) :: mask
    integer :: i, j
    integer :: totsize
    integer, dimension(1) :: pos
    integer, dimension(2) :: xpos, ypos
    integer :: lonloc, latloc, flattenedIndex, tempIndex
    character(30) :: row_bounds

    ! Associate names with variables defined in derived types.
    associate( &
    init_file               => c_switch%init_file,            & !< character(350): 
    rs_file_to_overwrite    => c_switch%rs_file_to_overwrite, & !< character(350): 
    metFileFss              => c_switch%metFileFss,           & !< character(350): 
    metFilefracFsf          => c_switch%metFilefracFsf,       & !< character(350): File name for the diffuse fraction of total shortwave radiation input file
    metFileFdl              => c_switch%metFileFdl,           & !< character(350): 
    metFilePre              => c_switch%metFilePre,           & !< character(350):
    metFileSnow             => c_switch%metFileSnow,          & !< character(350): File name for the snow flux file
    metFileTa               => c_switch%metFileTa,            & !< character(350): 
    metFileQa               => c_switch%metFileQa,            & !< character(350): 
    metFileUv               => c_switch%metFileUv,            & !< character(350): 
    metFilePres             => c_switch%metFilePres,          & !< character(350): 
    CO2File                 => c_switch%CO2File,              & !< character(350): 
    tracerCO2File           => c_switch%tracerCO2File,        & !< character(350): 
    CH4File                 => c_switch%CH4File,              & !< character(350): 
    POPDFile                => c_switch%POPDFile,             & !< character(350): 
    LGHTFile                => c_switch%LGHTFile,             & !< character(350): 
    LUCFile                 => c_switch%LUCFile,              & !< character(350): 
    FERFile                 => c_switch%FERFile,              & !< character(350): 
    DEPFile                 => c_switch%DEPFile,              & !< character(350): 
    timberHarvest           => c_switch%timberHarvest,        & !< logical:
    timberHarvestFile       => c_switch%timberHarvestFile,    & !< character(350): 
    prescribedFire           => c_switch%prescribedFire,      & !< logical: 
    prescribedFireFile       => c_switch%prescribedFireFile,  & !< character(350): 
    dynamicTilingOn     => c_switch%dynamicTilingOn,        & !< logical:
    OBSWETFFile             => c_switch%OBSWETFFile,          & !< character(350): 
    ctem_on                 => c_switch%ctem_on,              & !< logical: 
    Ncycle_on               => c_switch%Ncycle_on,            & !< logical: 
    projectedGrid           => c_switch%projectedGrid,        & !< logical: 
    doMethane               => c_switch%doMethane,            & !< logical: 
    dofire                  => c_switch%dofire,               & !< logical: 
    PFTCompetition          => c_switch%PFTCompetition,       & !< logical: 
    lnduseon                => c_switch%lnduseon,             & !< logical: 
    fertilizeron            => c_switch%fertilizeron,         & !< logical: 
    depositionon            => c_switch%depositionon,         & !< logical: 
    transientOBSWETF        => c_switch%transientOBSWETF,     & !< logical: 
    useTracer               => c_switch%useTracer,            & !< integer: 
    fixedYearLUC            => c_switch%fixedYearLUC,         & !< integer: 
    fixedYearFER            => c_switch%fixedYearFER,         & !< integer: 
    fixedYearDEP            => c_switch%fixedYearDEP,         & !< integer: 
    fixedYearOBSWETF        => c_switch%fixedYearOBSWETF,     & !< integer: 
    isnoalb                 => c_switch%isnoalb,              & !< integer: 
    alb4BandParamsFile      => c_switch%alb4BandParamsFile,   & !< character(350):
    blackCarbonFile         => c_switch%blackCarbonFile,      & !< character(350):
    blackCdepon             => c_switch%blackCdepon          & !< integer
    )

    ! ------------

    !> First, open initial conditions file.
    initid = ncOpen(init_file,NF90_NOWRITE)

    if (.not. projectedGrid) then

      !> Next, retrieve dimensions and allocate arrays to hold the lons and lats.
      !! We assume the file has 'lon' and 'lat' for names of longitude and latitude.

      totlon = ncGetDimLen(initid,'lon')
      totlat = ncGetDimLen(initid,'lat')
     
      allocate(myDomain%allLonValues(totlon),myDomain%allLatValues(totlat))

      !> Read in the coordinate variables from the initialization file.      
      myDomain%allLonValues = ncGetDimValues(initid, 'lon', count = (/totlon/))
      myDomain%allLatValues = ncGetDimValues(initid, 'lat', count = (/totlat/))

      !> Try and catch if the user has put in lon values from -180 to 180 or 0 to 360
      !! when the input file expects the opposite.
      if (myDomain%domainBounds(1) < 0. .and. myDomain%allLonValues(1) >= 0.) then
        myDomain%domainBounds(1) = 360. + myDomain%domainBounds(1)
        print * ,'Based on init_file,adjusted your domain (longitude) to',myDomain%domainBounds(1)
      end if
      if (myDomain%domainBounds(2) < 0. .and. myDomain%allLonValues(1) >= 0.) then
        myDomain%domainBounds(2) = 360. + myDomain%domainBounds(2)
        print * ,'Based on init_file,adjusted your domain (longitude) to',myDomain%domainBounds(2)
      end if
      if (myDomain%domainBounds(1) > 180. .and. myDomain%allLonValues(1) < 0.) then
        myDomain%domainBounds(1) = myDomain%domainBounds(1) - 360.
        print * ,'Based on init_file,adjusted your domain (longitude) to',myDomain%domainBounds(1)
      end if
      if (myDomain%domainBounds(2) > 180. .and. myDomain%allLonValues(1) < 0.) then
        myDomain%domainBounds(2) = myDomain%domainBounds(2) - 360.
        print * ,'Based on init_file,adjusted your domain (longitude) to',myDomain%domainBounds(2)
      end if

      ! FLAG - be good to put in a check here but need to do this better.
      !         !> Check that our domain is within the longitude and latitude limits of
      !         !! the input files. Otherwise print a warning. Primarily we are trying to
      !         !! catch instances where the input file runs from 0 to 360 longitude while
      !         !! the user expects -180 to 180.
      !         if (myDomain%domainBounds(1) < myDomain%allLonValues(1)) then ! W most lon
      !             print*,'=>Your domain bound ', myDomain%domainBounds(1),' is outside of', &
      !                 ' the limits of the init_file ',myDomain%allLonValues(1)
      !         else if (myDomain%domainBounds(2) > myDomain%allLonValues(ubound(myDomain%allLonValues,1))) then ! E most lon
      !             print*,'=>Your domain bound ', myDomain%domainBounds(2),' is outside of', &
      !                 ' the limits of the init_file ',myDomain%allLonValues(ubound(myDomain%allLonValues,1))
      !         else if (myDomain%domainBounds(3) < myDomain%allLatValues(1)) then ! S most lat
      !             print*,'=>Your domain bound ', myDomain%domainBounds(3),' is outside of', &
      !                 ' the limits of the init_file ',myDomain%allLatValues(1)
      !         else if (myDomain%domainBounds(4) > myDomain%allLatValues(ubound(myDomain%allLatValues,1))) then ! N most lat
      !             print*,'=>Your domain bound ', myDomain%domainBounds(4),' is outside of', &
      !                 ' the limits of the init_file ',myDomain%allLatValues(ubound(myDomain%allLatValues,1))
      !         end if

      !> Since the domainBounds are coordinates,need to find the indices corresponding to the domain bounds.

      if (myDomain%domainBounds(1) + myDomain%domainBounds(2) + &
          myDomain%domainBounds(3) + myDomain%domainBounds(4) == 0) then
        ! Special case, if the domainBounds are 0/0/0/0 then take whole domain.
        print * , ' domainBounds given = 0/0/0/0 so running whole domain of',totlon,' longitude cells and ',totlat,' latitude cells.'
        xpos(1) = 1
        xpos(2) = totlon
        ypos(1) = 1
        ypos(2) = totlat
      else
        ! Use the domain as given and obtain the indices of the lons and lats closest to the specified domainBounds.
        pos = minloc(abs(myDomain%allLonValues - myDomain%domainBounds(1)))
        xpos(1) = pos(1)

        pos = minloc(abs(myDomain%allLonValues - myDomain%domainBounds(2)))
        xpos(2) = pos(1)

        pos = minloc(abs(myDomain%allLatValues - myDomain%domainBounds(3)))
        ypos(1) = pos(1)

        pos = minloc(abs(myDomain%allLatValues - myDomain%domainBounds(4)))
        ypos(2) = pos(1)
      end if

      !> Obtain the starting indices of the coordinate vectors.

      myDomain%srtx = minval(xpos)
      myDomain%srty = minval(ypos)

      !> Ensure that the starting indices are within the specified bounds of the domain.
      !! Note that when the domain is a single point, the 1st/2nd elements of domainBounds are set to the
      !! the longitude and the 3rd/4th elements are set to the latitude (see readFromJobOptions.f90).
      !! In this case, the closest grid point to the specified coordinates is used and the check should not be done.

      if (myDomain%allLonValues(myDomain%srtx) < myDomain%domainBounds(1) .and. &
          myDomain%domainBounds(2) /= myDomain%domainBounds(1)) myDomain%srtx = myDomain%srtx + 1

      if (myDomain%allLatValues(myDomain%srty) < myDomain%domainBounds(3) .and. &
          myDomain%domainBounds(4) /= myDomain%domainBounds(3)) myDomain%srty = myDomain%srty + 1

      !> Compute the size of the coordinate vectors.

      myDomain%cntx = 1 + abs(maxval(xpos) - myDomain%srtx)
      myDomain%cnty = 1 + abs(maxval(ypos) - myDomain%srty)

      !> Ensure that the last index of each vector is within the domain bounds.

      if (myDomain%allLonValues(maxval(xpos)) > myDomain%domainBounds(2) .and. &
          myDomain%domainBounds(2) /= myDomain%domainBounds(1)) myDomain%cntx = myDomain%cntx - 1

      if (myDomain%allLatValues(maxval(ypos)) > myDomain%domainBounds(4) .and. &
          myDomain%domainBounds(4) /= myDomain%domainBounds(3)) myDomain%cnty = myDomain%cnty - 1

    else ! projected grid

      !> On a projected grid we have to use the grid cell indexes to delineate our domain to run
      !! over. We then use the indexes to determine the values of longitude and latitude for
      !! each grid cell.

      !> Retrieve dimensions. We assume the file has 'lon' and 'lat' for
      !! names of longitude and latitude dimensions.

      totlon = ncGetDimLen(initid,'lon')
      totlat = ncGetDimLen(initid,'lat')

      !> calculate the number and indices of the pixels to be calculated
      allocate(myDomain%allLonValues(totlat * totlon),myDomain%allLatValues(totlat * totlon))

      !> This will get all lon and lat grids as flattened vectors. This actually uses ncGetVar under the hood.
      !> NOTE: here, if you use the short names, it will fail since it will find the dimension which is 1D 
      !> the fail is catastrophic and won't give any useful info so beware.
      myDomain%allLonValues = ncGetDimValues(initid, 'longitude', count2D = (/totlon,totlat/))
      myDomain%allLatValues = ncGetDimValues(initid, 'latitude', count2D = (/totlon,totlat/))

      !> Since the domainBounds are indexes, and not coordinates, we can use them directly.
      xpos(1) = myDomain%domainBounds(1)
      xpos(2) = myDomain%domainBounds(2)
      ypos(1) = myDomain%domainBounds(3)
      ypos(2) = myDomain%domainBounds(4)

      !> Special case, if the domainBounds are 0/0/0/0 then take whole domain
      if (myDomain%domainBounds(1) + myDomain%domainBounds(2) + &
          myDomain%domainBounds(3) + myDomain%domainBounds(4) == 0) then
        print * , ' domainBounds given = 0/0/0/0 so running whole domain of',totlon,' longitude cells and ',totlat,' latitude cells.'
        xpos(1) = 1
        xpos(2) = totlon
        ypos(1) = 1
        ypos(2) = totlat
      end if

      myDomain%srtx = minval(xpos)
      myDomain%srty = minval(ypos)

      myDomain%cntx = 1 + abs(maxval(xpos) - myDomain%srtx)
      myDomain%cnty = 1 + abs(maxval(ypos) - myDomain%srty)

    end if

    !> Save the longitudes and latitudes over the region of interest for making the
    !! output files.
    totsize = myDomain%cntx * myDomain%cnty
    allocate(myDomain%latLandCell(totsize), &
    myDomain%lonLandCell(totsize), &
    myDomain%latLandIndex(totsize), &
    myDomain%lonLandIndex(totsize), &
    myDomain%latLocalIndex(totsize), &
    myDomain%lonLocalIndex(totsize))
    if (.not. projectedGrid) then
      allocate(myDomain%latUnique(myDomain%cnty), &
                   myDomain%lonUnique(myDomain%cntx))
    else
      allocate(myDomain%latUnique(totsize), &
                   myDomain%lonUnique(totsize))
    end if

    !> Retrieve the number of soil layers (set ignd !)
    testignd = ncGetDimLen(initid, 'layer')

    IF (testignd /= ignd) THEN
      write(6, * ) 'MISMATCH IN NUMBER OF SOIL LAYERS'
      write(6, * ) ignd, 'set in classicParams.f90'
      write(6, * ) testignd, 'in the initialization file'
      print*,'----- Run aborting. -----'
      stop
    END IF

    !> Grab the model domain. We use FLND since it is the land cells we want to run the model over.
    !! the 'Mask' variable is all land.
    allocate(mask(myDomain%cntx,myDomain%cnty))
    mask = ncGet2DVar(initid, 'FLND', start = [myDomain%srtx,myDomain%srty], &
           count = [myDomain%cntx,myDomain%cnty],format = [myDomain%cntx,myDomain%cnty])
    myDomain%LandCellCount = 0
    do i = 1,myDomain%cntx
      do j = 1,myDomain%cnty
        if (mask(i,j) > zero) then  ! Can choose a higher threshold but set to model default zero for now.
          ! print*, "(", i, ",", j, ") or (", myDomain%allLonValues(i + myDomain%srtx - 1) &
          ! , ",", myDomain%allLatValues(j + myDomain%srty - 1), ") is land"
          myDomain%LandCellCount = myDomain%LandCellCount + 1
          myDomain%lonLandIndex(myDomain%LandCellCount) = i + myDomain%srtx - 1
          myDomain%lonLocalIndex(myDomain%LandCellCount) = i
          myDomain%latLandIndex(myDomain%LandCellCount) = j + myDomain%srty - 1
          myDomain%latLocalIndex(myDomain%LandCellCount) = j
          if (.not. projectedGrid) then
            myDomain%lonLandCell(myDomain%LandCellCount) = myDomain%allLonValues(i + myDomain%srtx - 1)
            myDomain%latLandCell(myDomain%LandCellCount) = myDomain%allLatValues(j + myDomain%srty - 1)
          else ! projected grid so the lons and lats are flattened vectors representing their 2D grids
            ! print*, "(", i, ",", j, ") or (", myDomain%allLonValues(flattenedIndex) &
            ! , ",", myDomain%allLatValues(flattenedIndex), ") is valid"
            flattenedIndex = (j + myDomain%srty - 2) * totlon + (i + myDomain%srtx - 1)
            myDomain%lonLandCell(myDomain%LandCellCount) = myDomain%allLonValues(flattenedIndex)
            myDomain%latLandCell(myDomain%LandCellCount) = myDomain%allLatValues(flattenedIndex)
          end if
        end if
      end do
    end do

    ! Extract lat/lon for the part of the grid that is being processed. This could be a subgrid or the full domain.
    ! (This has been split out of the previous loop as it is independent of the mask) - EC.

    if (.not. projectedGrid) then
      do i = 1,myDomain%cntx
        myDomain%lonUnique(i) = myDomain%allLonValues(i + myDomain%srtx - 1)
      end do
      do j = 1,myDomain%cnty
        myDomain%latUnique(j) = myDomain%allLatValues(j + myDomain%srty - 1)
      end do
    else
      do j = 1,myDomain%cnty
        do i = 1,myDomain%cntx
          flattenedIndex = (j + myDomain%srty - 2) * totlon + (i + myDomain%srtx - 1)
          tempIndex = (j - 1) * myDomain%cntx + i
          myDomain%lonUnique(tempIndex) = myDomain%allLonValues(flattenedIndex)
          myDomain%latUnique(tempIndex) = myDomain%allLatValues(flattenedIndex)
        end do
      end do
    end if

    if (myDomain%LandCellCount == 0) then
      print * ,'=>Your domain is not land my friend.'
      if (.not. projectedGrid) then
        if (closeEnough(myDomain%domainBounds(1),myDomain%domainBounds(2),1.E-5)) then ! point run
          lonloc = closestCell(initid,'lon',myDomain%domainBounds(1))
          latloc = closestCell(initid,'lat',myDomain%domainBounds(3))
          print * ,'Closest grid cell is ',myDomain%allLonValues(lonloc),'/',myDomain%allLatValues(latloc)
          print * ,'but that may not be land. Check your input files to be sure'
        end if
      end if
    end if

    !> To determine testnmos, based on the tile dimension in init file
    testnmos = ncGetDimLen(initid, 'tile')

    IF (testnmos /= nmos) THEN
      write(6, * ) 'MISMATCH IN NUMBER OF MOSAIC TILES'
      write(6, * ) nmos, 'set in classicParams.f90'
      write(6, * ) testnmos, 'in the initialization file'
      print*,'----- Run aborting. -----'
      stop
    END IF

    !> implement two new checks for ican and icc

    test2ican = ncGetDimLen(initid, 'ic')

    IF (ican /= test2ican) THEN
      write(6, * ) 'MISMATCH IN NUMBER OF CLASS PFTs'
      write(6, * ) ican, 'set in classicParams.f90'
      write(6, * ) test2ican, 'in the initialization file'
      print*,'----- Run aborting. -----'
      stop
    END IF

    test2icc = ncGetDimLen(initid, 'icc')

    IF (icc /= test2icc) THEN
      write(6, * ) 'MISMATCH IN NUMBER OF CTEM PFTs'
      write(6, * ) icc, 'set in classicParams.f90'
      write(6, * ) test2icc, 'in the initialization file'
      print*,'----- Run aborting. -----'
      stop
    END IF

    !> Lastly, open some files so they are ready

    rsid = ncOpen(rs_file_to_overwrite,nf90_write)

    !> Add global attribute to restart file, storing the rows overwritten.
    !> Could be used in the future to simplify stitching of the restart file when the run is split across multiple nodes.
    !write(row_bounds,'(I0,x,I0)') ypos
    write(row_bounds,'(I0,x,I0)') myDomain%srty,myDomain%srty+myDomain%cnty-1
    call ncReDef(rsid)
    call ncPutAtt(rsid,nf90_global,'row_bounds',charvalues = trim(row_bounds))
    call ncEndDef(rsid)
    
    !> If using the four band albedo parameterization, open the file containing the look up table for the 
    !! parameters and possibly the black carbon deposition file.
    if (isnoalb == 1) then 
      albid = ncOpen(alb4BandParamsFile,nf90_nowrite)
      if (blackCdepon) then
        bcid = ncOpen(blackCarbonFile,nf90_nowrite)
        call checkTimeUnits(bcid,blackCarbonFile)  ! Check that the file has the expected time units      
        bcDepVarName = ncGetVarName(bcid)
      end if 
    end if 

    if (ctem_on) then
      co2id = ncOpen(CO2File,nf90_nowrite)
      call checkTimeUnits(co2id,CO2File)  ! Check that the file has the expected time units      
      co2VarName = ncGetVarName(co2id)
      
      if (doMethane) then
        ch4id = ncOpen(CH4File,nf90_nowrite)
        call checkTimeUnits(ch4id,CH4File)  ! Check that the file has the expected time units      
        ch4VarName = ncGetVarName(ch4id)
      end if
      if (useTracer > 0) then
        tracerco2id = ncOpen(tracerCO2File,nf90_nowrite)
        call checkTimeUnits(tracerco2id,tracerCO2File)  ! Check that the file has the expected time units      
        ! 14C has different values depending on latitudional bands. The expected
        ! input file is from CMIP6 and it splits the bands as follows:
        ! Southern Hemisphere (30-90°S), Tropics (30°S-30°N),
        ! and Northern Hemisphere (30-90°N). Since at this stage of
        ! CLASSIC we don't know the latitude of the cell being simulated
        ! we need to get tracerco2VarName later for 14C simulations.
        if (useTracer /= 2) tracerco2VarName = ncGetVarName(tracerco2id)

      end if
      if (dofire) then
        popid = ncOpen(POPDFile,nf90_nowrite)
        call checkTimeUnits(popid,POPDFile)  ! Check that the file has the expected time units      
        popVarName = ncGetVarName(popid)
      
        lghtid = ncOpen(LGHTFile,nf90_nowrite)
        call checkTimeUnits(lghtid,LGHTFile)  ! Check that the file has the expected time units      
        lghtVarName = ncGetVarName(lghtid)
      end if
      if (lnduseon .or. (fixedYearLUC /= - 9999)) then
        lucid = ncOpen(LUCFile,nf90_nowrite)
        call checkTimeUnits(lucid,LUCFile)  ! Check that the file has the expected time units      
        lucVarName = ncGetVarName(lucid)
      end if
      if (transientOBSWETF .or. (fixedYearOBSWETF /= - 9999)) then
        obswetid = ncOpen(OBSWETFFile,nf90_nowrite)
        call checkTimeUnits(obswetid,OBSWETFFile)  ! Check that the file has the expected time units      
        obswetVarName = ncGetVarName(obswetid)
      end if
      if (Ncycle_on .and. (fertilizeron .or. (fixedYearFER /= - 9999))) then
        ferid = ncOpen(FERFile,nf90_nowrite)
        call checkTimeUnits(ferid,FERFile)  ! Check that the file has the expected time units
        ferVarName = ncGetVarName(ferid)
      end if
      if (Ncycle_on .and. (depositionon .or. (fixedYearDEP /= - 9999))) then
        depid = ncOpen(DEPFile,nf90_nowrite)
        call checkTimeUnits(depid,DEPFile)  ! Check that the file has the expected time units
        depVarName = ncGetVarName(depid)
      end if
    end if

    if (timberHarvest) then

        !perform a check and stop the model if CTEM is off, competition is on or agricultural LUC is on
        if ((.not. (ctem_on)) .or. (lnduseon) .or. (PFTCompetition)) then
          write(6, * ) 'dynamic tiling and timber harvest not avaliable with competition, agricultural LUC or when ctem is off'
          print*,'----- Run aborting. -----'
          stop
        end if

        timid = ncOpen(timberHarvestFile,nf90_nowrite)
        call checkTimeUnits(timid,timberHarvestFile)  ! Check that the file has the expected time units      
        timVarName = ncGetVarName(timid)
    end if

    if (prescribedFire) then

        !perform a check and stop the model if CTEM is off, competition is on or agricultural LUC is on
        if ((.not. (ctem_on)) .or. (lnduseon) .or. (PFTCompetition)) then
          write(6, * ) 'dynamic tiling and prescribed fire harvest not avaliable with competition, agricultural LUC or when ctem is off'
          print*,'----- Run aborting. -----'
          stop
        end if

        fireid = ncOpen(prescribedFireFile,nf90_nowrite)
        call checkTimeUnits(fireid,prescribedFireFile)  ! Check that the file has the expected time units      
        fireVarName = ncGetVarName(fireid)
    end if

    !> Open the meteorological forcing files and find the variable name in the file
    metFssId    = ncOpen(metFileFss,nf90_nowrite)
    call checkTimeUnits(metFssId,metFileFss)  ! Check that the file has the expected time units      
    metFssVarName = ncGetVarName(metFssId)
    !> Check if we have a diffuse shortwave radiation file, if so do as the rest of the met vars
    if (trim(metFilefracFsf) /= '') then
      metfracFsfId    = ncOpen(metFilefracFsf,nf90_nowrite)
      call checkTimeUnits(metfracFsfId,metFilefracFsf)  ! Check that the file has the expected time units      
      metfracFsfVarName = ncGetVarName(metfracFsfId)
    end if 

    metFdlId    = ncOpen(metFileFdl,nf90_nowrite)
    call checkTimeUnits(metFdlId,metFileFdl)  ! Check that the file has the expected time units      
    metFdlVarName = ncGetVarName(metFdlId)
    metPreId    = ncOpen(metFilePre,nf90_nowrite)
    call checkTimeUnits(metPreId,metFilePre)  ! Check that the file has the expected time units      
    metPreVarName = ncGetVarName(metPreId)

    !> Check if we have a fraction of precip that is snow file, if so do as the rest of the met vars
    if (trim(metFileSnow) /= '') then
      metSnowId    = ncOpen(metFileSnow,nf90_nowrite)
      call checkTimeUnits(metSnowId,metFileSnow)  ! Check that the file has the expected time units      
      metSnowVarName = ncGetVarName(metSnowId)
    end if 

    metTaId     = ncOpen(metFileTa,nf90_nowrite)
    call checkTimeUnits(metTaId,metFileTa)  ! Check that the file has the expected time units      
    metTaVarName = ncGetVarName(metTaId)
    metQaId     = ncOpen(metFileQa,nf90_nowrite)
    call checkTimeUnits(metQaId,metFileQa)  ! Check that the file has the expected time units      
    metQaVarName = ncGetVarName(metQaId)
    metUvId     = ncOpen(metFileUv,nf90_nowrite)
    call checkTimeUnits(metUvId,metFileUv)  ! Check that the file has the expected time units      
    metUvVarName = ncGetVarName(metUvId)
    metPresId   = ncOpen(metFilePres,nf90_nowrite)
    call checkTimeUnits(metPresId,metFilePres)  ! Check that the file has the expected time units      
    metPresVarName = ncGetVarName(metPresId)

    ! Lastly, create an output file for the model domain (FLND) so that it is readily
    ! accessible for plotting/processing model outputs.
    call writeStaticVars('sftlf', 'Percentage of the grid cell occupied by land (including lakes)', '%', mask * 100.)

    deallocate(mask)

    end associate
  end subroutine read_modelsetup

  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_read_initialstate
  !! @{
  !> Reads in the model initial conditions for both physics and biogeochemistry (if CTEM on)
  !> @author Joe Melton

  subroutine read_initialstate (lonIndex,latIndex)

    ! J. Melton
    ! Nov 2016

    use ctemStateVars,  only : c_switch, vrot, vgat, tracer
    use classStateVars, only : class_rot, class_gat
    use classicParams,  only : icc, iccp1, iccp2, nmos, ignd, icp1, nlat, ican, pi, crop, TFREZ, &
                               RSMN, QA50, VPDA, VPDB, PSGA, PSGB, &
                               albdif_lut, albdir_lut, trandif_lut, trandir_lut, &
                               nsmu, nsalb, nbc, nreff, nswe, nbnd_lut, zero, nTilepres

    implicit none

    ! arguments
    integer, intent(in) :: lonIndex, latIndex

    ! local variables

    integer :: i, m, j, n
    real :: bots
    ! integer :: ipnt !, albdim
    ! integer :: isalb, ismu, isgs, iswe, ibc
    ! real, dimension(:,:,:), allocatable :: tmpalb

    ! Associate names with variables defined in derived types.
    associate( &
    ctem_on           => c_switch%ctem_on,              & !< logical: 
    Ncycle_on         => c_switch%Ncycle_on,            & !< logical: 
    fertilizeron      => c_switch%fertilizeron,         & !< logical: 
    depositionon      => c_switch%depositionon,         & !< logical: 
    dofire            => c_switch%dofire,               & !< logical: 
    PFTCompetition    => c_switch%PFTCompetition,       & !< logical: 
    inibioclim        => c_switch%inibioclim,           & !< logical: 
    start_bare        => c_switch%start_bare,           & !< logical: 
    lnduseon          => c_switch%lnduseon,             & !< logical: 
    useTracer         => c_switch%useTracer,            & !< integer: useTracer = 0, the tracer code is not used.
                                                          !!          useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the
                                                          !!                        tracer values in the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                                          !!          useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                                          !!          useTracer = 3 [Not implemented yet means the tracer is 13C and will then call a 13C fractionation scheme. 
    isnoalb           => c_switch%isnoalb,              & !< integer: 
    dynamicTilingOn     => c_switch%dynamicTilingOn,        & !< logical:
    prescribedFire      => c_switch%prescribedFire,        & !< logical:
    timberHarvest       => c_switch%timberHarvest,        & !< logical:

    tileAgeReset      => c_switch%tileAgeReset,         & !< logical
    trackTileAge      => c_switch%trackTileAge,    & !< logical:
    doPeat            => c_switch%doPeat,          & !< logical
    fcancmxrow        => vrot%fcancmx,                  & !< real, dimension(:,:,:) :! 
    gleafmasrow       => vrot%gleafmas,                 & !< real, dimension(:,:,:) : Green leaf mass for each of the CTEM pfts, \f$kg c/m^2\f$ 
    gleafmas_NSrow     => vrot%gleafmas_ns,             & !< real, dimension(:,:,:) :! 
    gleafmassrow      => vrot%gleafmas_s,               & !< real, dimension(:,:,:) :! 
    bleafmasrow       => vrot%bleafmas,                 & !< real, dimension(:,:,:) : Brown leaf mass for each of the CTEM pfts, \f$kg c/m^2\f$ 
    stemmassrow       => vrot%stemmass,                 & !< real, dimension(:,:,:) : Stem mass for each of the CTEM pfts, \f$kg c/m^2\f$ 
    stemmass_NSrow     => vrot%stemmass_ns,             & !< real, dimension(:,:,:) :! 
    stemmasssrow      => vrot%stemmass_s,               & !< real, dimension(:,:,:) :! 
    rootmassrow       => vrot%rootmass,                 & !< real, dimension(:,:,:) : Root mass for each of the CTEM pfts, \f$kg c/m^2\f$ 
    rootmass_NSrow     => vrot%rootmass_ns,             & !< real, dimension(:,:,:) :! 
    rootmasssrow      => vrot%rootmass_s,               & !< real, dimension(:,:,:) :! 
    lygleafmasmaxrow  => vrot%lygleafmasmax,            & !< real, dimension(:,:,:) : last year maximum of the green leaf mass for each PFT, \f$kg c/m^2\f$ 
    lystemmassmaxrow  => vrot%lystemmassmax,            & !< real, dimension(:,:,:) : last year maximum of the stem mass for each PFT, \f$kg c/m^2\f$ 
    lyrootmassmaxrow  => vrot%lyrootmassmax,            & !< real, dimension(:,:,:) : last year maximum of the root mass for each PFT, \f$kg c/m^2\f$ 
    pstemmassrow      => vrot%pstemmass,                & !< real, dimension(:,:,:) : Stem mass from previous timestep, is value before fire. used by burntobare subroutine 
    pgleafmassrow     => vrot%pgleafmass,               & !< real, dimension(:,:,:) : Green leaf mass from previous timestep, is value before fire. used by burntobare subroutine 
    grwtheffrow       => vrot%grwtheff,                 & !< real, dimension(:,:,:) : growth efficiency. change in biomass per year per unit max.
                                                          !!                          lai (\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine 
    nppvegrow         => vrot%nppveg,                   & !< real, dimension(:,:,:)
    twarmmrow         => vrot%twarmm,                   & !< real, dimension(:,:) : temperature of the warmest month (c) 
    tcoldmrow         => vrot%tcoldm,                   & !< real, dimension(:,:) : temperature of the coldest month (c) 
    gdd5row           => vrot%gdd5,                     & !< real, dimension(:,:) : growing degree days above 5 c 
    aridityrow        => vrot%aridity,                  & !< real, dimension(:,:) : aridity index, ratio of potential evaporation to precipitation 
    srplsmonrow       => vrot%srplsmon,                 & !< real, dimension(:,:) : number of months in a year with surplus water i.e.precipitation more than potential evaporation 
    defctmonrow       => vrot%defctmon,                 & !< real, dimension(:,:) : number of months in a year with water deficit i.e.precipitation less than potential evaporation 
    anndefctrow       => vrot%anndefct,                 & !< real, dimension(:,:) : annual water deficit (mm) 
    annsrplsrow       => vrot%annsrpls,                 & !< real, dimension(:,:) : annual water surplus (mm) 
    annpcprow         => vrot%annpcp,                   & !< real, dimension(:,:) : annual precipitation (mm) 
    dry_season_lengthrow => vrot%dry_season_length,        & !< real, dimension(:,:) : length of dry season (months) 
    litrmassrow       => vrot%litrmass,                 & !< real, dimension(:,:,:,:) : Litter mass for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    soilcmasrow       => vrot%soilcmas,                 & !< real, dimension(:,:,:,:) : Soil C mass for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    slopefrac         => vrot%slopefrac,                & !< real, dimension(:,:,:) : 
    lfstatusrow       => vrot%lfstatus,                 & !< integer, dimension(:,:,:) : 
    pandaysrow        => vrot%pandays,                  & !< integer, dimension(:,:,:) : 
    ipeatlandrow      => vrot%ipeatland,                & !< integer, dimension(:,:) : Peatland switch: 0 = not a peatland, 1 = bog, 2 = fen 
    Cmossmas          => vrot%Cmossmas,                 & !< real, dimension(:,:) : Carbon in moss biomass, \f$kg C/m^2\f$ 
    litrmsmoss        => vrot%litrmsmoss,               & !< real, dimension(:,:) : moss litter mass, \f$kg C/m^2\f$ 
    dmoss             => vrot%dmoss,                    & !< real, dimension(:,:) : depth of living moss (m) 
    peatSoilC         => vrot%peatSoilC,                & !< real, dimension(:,:) : peat soil C mass, \f$kg C/m^2\f$ 
    nh4_massrow       => vrot%nh4_mass,                 & !< real, dimension(:,:,:) :ammonium mass (\f$g N/m^2\f$) 
    no3_massrow       => vrot%no3_mass,                 & !< real, dimension(:,:,:) :nitrate mass (\f$g N/m^2\f$) 
    nlitrmassrow      => vrot%nlitrmass,                & !< real, dimension(:,:,:) :litter nitrogen mass (\f$g N/m^2\f$) 
    soilnmasrow       => vrot%soilnmas,                 & !< real, dimension(:,:,:) :soil organic nitrogen mass (\f$g N/m^2\f$) 
    ngleafmasrow      => vrot%ngleafmas,                & !< real, dimension(:,:,:) :green leaf nitrogen mass (\f$g N/m^2\f$) 
    ngleafmas_NSrow    => vrot%ngleafmas_ns,            & !< real, dimension(:,:,:) :non-structural green leaf nitrogen mass (\f$g N/m^2\f$) 
    ngleafmassrow     => vrot%ngleafmas_s,              & !< real, dimension(:,:,:) :structural green leaf nitrogen mass (\f$g N/m^2\f$) 
    nbleafmasrow      => vrot%nbleafmas,                & !< real, dimension(:,:,:) :brown leaf nitrogen mass (\f$g N/m^2\f$) 
    nstemmassrow      => vrot%nstemmass,                & !< real, dimension(:,:,:) :stem nitrogen mass (\f$g N/m^2\f$) 
    nstemmass_NSrow    => vrot%nstemmass_ns,            & !< real, dimension(:,:,:) :non-structural stem nitrogen mass (\f$g N/m^2\f$) 
    nstemmasssrow     => vrot%nstemmass_s,              & !< real, dimension(:,:,:) :structural stem nitrogen mass (\f$g N/m^2\f$) 
    nrootmassrow      => vrot%nrootmass,                & !< real, dimension(:,:,:) :root nitrogen mass (\f$g N/m^2\f$) 
    nrootmass_NSrow    => vrot%nrootmass_ns,            & !< real, dimension(:,:,:) :non-structural root nitrogen mass (\f$g N/m^2\f$) 
    nrootmasssrow     => vrot%nrootmass_s,              & !< real, dimension(:,:,:) :structural root nitrogen mass (\f$g N/m^2\f$) 
    bnf_natrow        => vrot%bnf_nat,                  & !< real, dimension(:,:,:) :natural bnf (\f$g N/m^2/d\f$) 
    bnf_antrow        => vrot%bnf_ant,                  & !< real, dimension(:,:,:) :anthropogenic bnf (\f$g N/m^2/d\f$) 
    soilpHrow         => vrot%soilpH,                   & !< real, dimension(:,:)   : 
    grclarearow       => vrot%grclarearow,              & !< real, dimension(:)     : area of the grid cell, \f$km^2\f$ 
    flhrlossrow        => vrot%flhrloss,                   & !< real, dimension(:,:) : 
    flhrloss_nsrow  => vrot%flhrloss_ns,                & !< real, dimension(:,:) : 
    flhrloss_srow   => vrot%flhrloss_s,                 & !< real, dimension(:,:) : 
    stmhrlosrow        => vrot%stmhrlos,                 & !< real, dimension(:,:) : 
    rothrlosrow        => vrot%rothrlos,                 & !< real, dimension(:,:) : 
    lystmmasrow       => vrot%lystmmas,                 & !< real, dimension(:,:,:) : 
    lyrotmasrow       => vrot%lyrotmas,                 & !< real, dimension(:,:,:) : 
    tymaxlairow       => vrot%tymaxlai,                 & !< real, dimension(:,:,:) :
    cfluxcgrow        => vrot%cfluxcg,                  & !< real, dimension(:,:) :
    cfluxcsrow        => vrot%cfluxcs,                  & !< real, dimension(:,:) :
    colddays_leaffallrow => vrot%colddays_leaffall,     & !< integer, dimension(:,:) :
    colddays_harvestrow  => vrot%colddays_harvest,      & !< integer, dimension(:,:) : 
    co2i1cgrow        => vrot%co2i1cg,                  & !< real, dimension(:,:,:) :
    co2i1csrow        => vrot%co2i1cs,                  & !< real, dimension(:,:,:) :
    co2i2cgrow        => vrot%co2i2cg,                  & !< real, dimension(:,:,:) :
    co2i2csrow        => vrot%co2i2cs,                  & !< real, dimension(:,:,:) :
    controlVector     => vrot%controlVector,            & !< integer, dimension(nlat,nmos)
    tileAgerow        => vrot%tileAgerow,               & ! !< real, dimension(nlat,nmos)

    tracerGLeafMass   => tracer%gLeafMassrot,           & !< real, dimension(:,:,:) : Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerBLeafMass   => tracer%bLeafMassrot,           & !< real, dimension(:,:,:) : Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerStemMass    => tracer%stemMassrot,            & !< real, dimension(:,:,:) : Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerRootMass    => tracer%rootMassrot,            & !< real, dimension(:,:,:) : Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerLitrMass    => tracer%litrMassrot,            & !< real, dimension(:,:,:,:) : Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    tracerSoilCMass   => tracer%soilCMassrot,           & !< real, dimension(:,:,:,:) : Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    tracerMossCMass   => tracer%mossCMassrot,           & !< real, dimension(:,:) : Tracer mass in moss biomass, \f$kg C/m^2\f$ 
    tracerMossLitrMass => tracer%mossLitrMassrot,       & !< real, dimension(:,:) : Tracer mass in moss litter, \f$kg C/m^2\f$ 
    FCANROT           => class_rot%FCANROT,             & !< real, dimension(:,:,:) : Maximum fractional coverage of modelled 
    FAREROT           => class_rot%FAREROT,             & !< real, dimension(:,:)   :
    RSMNROT           => class_rot%RSMNROT,             & !< real, dimension(:,:,:) : 
    QA50ROT           => class_rot%QA50ROT,             & !< real, dimension(:,:,:) : 
    VPDAROT           => class_rot%VPDAROT,             & !< real, dimension(:,:,:) : 
    VPDBROT           => class_rot%VPDBROT,             & !< real, dimension(:,:,:) : 
    PSGAROT           => class_rot%PSGAROT,             & !< real, dimension(:,:,:) : 
    PSGBROT           => class_rot%PSGBROT,             & !< real, dimension(:,:,:) : 
    DRNROT            => class_rot%DRNROT,              & !< real, dimension(:,:)   : 
    SDEPROT           => class_rot%SDEPROT,             & !< real, dimension(:,:)   : 
    XSLPROT           => class_rot%XSLPROT,             & !< real, dimension(:,:)   : 
    GRKFROT           => class_rot%GRKFROT,             & !< real, dimension(:,:)   : 
    WFSFROT           => class_rot%WFSFROT,             & !< real, dimension(:,:)   : 
    WFCIROT           => class_rot%WFCIROT,             & !< real, dimension(:,:)   : 
    !MIDROT            => class_rot%MIDROT,              & !< integer, dimension(:,:) : 
    DELZ              => class_gat%DELZ,                & !< real, dimension(:)     : 
    ZBOT              => class_gat%ZBOT,                & !< real, dimension(:)     : 
    SANDROT           => class_rot%SANDROT,             & !< real, dimension(:,:,:) : 
    CLAYROT           => class_rot%CLAYROT,             & !< real, dimension(:,:,:) : 
    ORGMROT           => class_rot%ORGMROT,             & !< real, dimension(:,:,:) : 
    TBARROT           => class_rot%TBARROT,             & !< real(r8), dimension(:,:,:) : 
    THLQROT           => class_rot%THLQROT,             & !< real, dimension(:,:,:) : 
    THICROT           => class_rot%THICROT,             & !< real, dimension(:,:,:) : 
    TCANROT           => class_rot%TCANROT,             & !< real, dimension(:,:)   : Vegetation canopy temperature [K]
    TSNOROT           => class_rot%TSNOROT,             & !< real, dimension(:,:)   : 
    TPNDROT           => class_rot%TPNDROT,             & !< real, dimension(:,:)   : 
    ZPNDROT           => class_rot%ZPNDROT,             & !< real, dimension(:,:)   : 
    RCANROT           => class_rot%RCANROT,             & !< real, dimension(:,:)   : 
    SCANROT           => class_rot%SCANROT,             & !< real, dimension(:,:)   : 
    SNOROT            => class_rot%SNOROT,              & !< real, dimension(:,:)   : 
    ALBSROT           => class_rot%ALBSROT,             & !< real, dimension(:,:)   : 
    RHOSROT           => class_rot%RHOSROT,             & !< real, dimension(:,:)   : 
    GROROT            => class_rot%GROROT,              & !< real, dimension(:,:)   : 
    !GCROW             => class_rot%GCROW,               & !< real, dimension(:)     : Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land) 
    ALVCROT           => class_rot%ALVCROT,             & !< real, dimension(:,:,:) : 
    ALICROT           => class_rot%ALICROT,             & !< real, dimension(:,:,:) : 
    PAMNROT           => class_rot%PAMNROT,             & !< real, dimension(:,:,:) : 
    PAMXROT           => class_rot%PAMXROT,             & !< real, dimension(:,:,:) : 
    LNZ0ROT           => class_rot%LNZ0ROT,             & !< real, dimension(:,:,:) : 
    CMASROT           => class_rot%CMASROT,             & !< real, dimension(:,:,:) : 
    ROOTROT           => class_rot%ROOTROT,             & !< real, dimension(:,:,:) : 
    DLATROW           => class_rot%DLATROW,             & !< real, dimension(:)     : 
    DLONROW           => class_rot%DLONROW,             & !< real, dimension(:)     : 
    RADJROW           => class_rot%RADJROW,             & !< real, dimension(:)     : Latitude of grid cell (positive north of equator) [rad] 
    Z0ORROW           => class_rot%Z0ORROW,             & !< real, dimension(:)     : 
    GGEOROW           => class_rot%GGEOROW,             & !< real, dimension(:)     : Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$ 
    SOCIROT           => class_rot%SOCIROT,             & !< real, dimension(:,:)   : 
    TBASROT           => class_rot%TBASROT,             & !< real, dimension(:,:)   : 
    CMAIROT           => class_rot%CMAIROT,             & !< real, dimension(:,:)   : 
    WSNOROT           => class_rot%WSNOROT,             & !< real, dimension(:,:)   : 
    ZSNLROT           => class_rot%ZSNLROT,             & !< real, dimension(:,:)   : Limiting snow depth (m) 
    TSFSROT           => class_rot%TSFSROT,             & !< real, dimension(:,:,:) : Ground surface temperature over subarea [K] 
    TACROT            => class_rot%TACROT,              & !< real, dimension(:,:)   : Temperature of air within vegetation canopy \f$[K] (T_{ac} )\f$ 
    QACROT            => class_rot%QACROT,              & !< real, dimension(:,:)   : Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ] (q_{ac} )\f$ 
    ITCTROT           => class_rot%ITCTROT,             & !< integer, dimension(:,:,:,:) : Counter of number of iterations required to solve surface energy balance for the elements of the four subareas 
    maxAnnualActLyr   => class_rot%maxAnnualActLyrROT   & !< real, dimension(:,:)   : Active layer depth maximum over the e-folding period specified by parameter eftime (m). 
    )

    ! ----------------------------

    do i = 1,nlat
      RADJROW(i) = DLATROW(i) * PI/180.
      Z0ORROW(i) = 0.0
      GGEOROW(i) = 0.0
    end do

    if(dynamicTilingOn) controlVector = ncGet2DVar(initid, 'control_vector', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    if((trackTileAge .or. dynamicTilingOn) .and. .not.(tileAgeReset)) tileAgerow = ncGet2DVar(initid, 'tile_age', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    if((trackTileAge .or. dynamicTilingOn) .and. tileAgeReset) tileAgerow = 0

    !< check the tiling setup
    if(dynamicTilingOn) then
      if(.not. (timberHarvest .or. prescribedFire)) then
        write(6, * ) 'Dynamic tiling is not avaliable'
              write(6, * ) 'without perscribed fire or harvest'
              print*,'----- Run aborting. -----'
                call errorHandler('read_initialstate', - 2)
          stop
      else if(timberHarvest .and. prescribedFire) then
        if((nmos - count(controlVector == 0)) <= (nTilepres + 2)) then
          write(6, * ) 'Not enough tiles to support this scenario'
                print*,'----- Run aborting. -----'
                call errorHandler('read_initialstate', - 2)
          stop
        end if
      else if (timberHarvest .or. prescribedFire) then
        if((nmos - count(controlVector == 0))<= (nTilepres + 1)) then
          write(6, * ) 'Not enough tiles to support this scenario'
                print*,'----- Run aborting. -----'
                call errorHandler('read_initialstate', - 2)
          stop
        end if
      end if
    end if

    DRNROT = ncGet2DVar(initid, 'DRN', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    SDEPROT = ncGet2DVar(initid, 'SDEP', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    SOCIROT = ncGet2DVar(initid, 'SOCI', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    FAREROT = ncGet2DVar(initid, 'FARE', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    ! The following four variables are not presently in use. Comment out read so not needed to be in input file.
    ! XSLPROT = ncGet2DVar(initid, 'XSLP', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    ! GRKFROT = ncGet2DVar(initid, 'GRKF', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    ! WFSFROT = ncGet2DVar(initid, 'WFSF', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    ! WFCIROT = ncGet2DVar(initid, 'WFCI', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    TCANROT = ncGet2DVar(initid, 'TCAN', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    TSNOROT = ncGet2DVar(initid, 'TSNO', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    TPNDROT = ncGet2DVar(initid, 'TPND', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    ZPNDROT = ncGet2DVar(initid, 'ZPND', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    RCANROT = ncGet2DVar(initid, 'RCAN', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    SCANROT = ncGet2DVar(initid, 'SCAN', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    SNOROT = ncGet2DVar(initid, 'SNO', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    ALBSROT = ncGet2DVar(initid, 'ALBS', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    RHOSROT = ncGet2DVar(initid, 'RHOS', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    GROROT = ncGet2DVar(initid, 'GRO', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    maxAnnualActLyr = ncGet2DVar(initid, 'maxAnnualActLyr', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    LNZ0ROT = ncGet3DVar(initid, 'LNZ0', start = [lonIndex,latIndex,1,1], count = [1,1,icp1,nmos], format = [nlat,nmos,icp1])
    ALVCROT = ncGet3DVar(initid, 'ALVC', start = [lonIndex,latIndex,1,1], count = [1,1,icp1,nmos], format = [nlat,nmos,icp1])
    ALICROT = ncGet3DVar(initid, 'ALIC', start = [lonIndex,latIndex,1,1], count = [1,1,icp1,nmos], format = [nlat,nmos,icp1])
    PAMNROT = ncGet3DVar(initid, 'PAMN', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    PAMXROT = ncGet3DVar(initid, 'PAMX', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    CMASROT = ncGet3DVar(initid, 'CMAS', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    ROOTROT = ncGet3DVar(initid, 'ROOT', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])

    ! The following six are parameters that can be made to spatially vary by uncommenting below and including them in the
    ! model init file. However, in practice these parameters are used with spatially invariable values so are read in from
    ! the CLASSIC namelist in classicParams.f90.
    ! RSMNROT = ncGet3DVar(initid, 'RSMN', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    ! QA50ROT = ncGet3DVar(initid, 'QA50', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    ! VPDAROT = ncGet3DVar(initid, 'VPDA', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    ! VPDBROT = ncGet3DVar(initid, 'VPDB', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    ! PSGAROT = ncGet3DVar(initid, 'PSGA', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    ! PSGBROT = ncGet3DVar(initid, 'PSGB', start = [lonIndex,latIndex,1,1], count = [1,1,ican,nmos], format = [nlat,nmos,ican])
    ! Here we apply the values read in from the namelist file:
    do i = 1,nlat
      do m = 1,nmos
        RSMNROT(i,m,:) = RSMN(:)
        QA50ROT(i,m,:) = QA50(:)
        VPDAROT(i,m,:) = VPDA(:)
        VPDBROT(i,m,:) = VPDB(:)
        PSGAROT(i,m,:) = PSGA(:)
        PSGBROT(i,m,:) = PSGB(:)
      end do
    end do

    SANDROT = ncGet3DVar(initid, 'SAND', start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos], format = [nlat,nmos,ignd])
    CLAYROT = ncGet3DVar(initid, 'CLAY', start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos], format = [nlat,nmos,ignd])
    ORGMROT = ncGet3DVar(initid, 'ORGM', start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos], format = [nlat,nmos,ignd])
    TBARROT = ncGet3DVar(initid, 'TBAR', start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos], format = [nlat,nmos,ignd])
    THLQROT = ncGet3DVar(initid, 'THLQ', start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos], format = [nlat,nmos,ignd])
    THICROT = ncGet3DVar(initid, 'THIC', start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos], format = [nlat,nmos,ignd])
    ipeatlandrow = ncGet2DVar(initid, 'ipeatland', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
    
    !set a general flag that tells the model if there are any peatland cells
    if (any(ipeatlandrow == 1 .or. ipeatlandrow == 2)) doPeat = .TRUE.

    DELZ = ncGet1DVar(initid, 'DELZ', start = [1], count = [ignd])

    ! From DELZ we can find ZBOT as:
    bots = 0.
    do n = 1,ignd
      bots = bots + delz(n)
      ZBOT(n) = bots
    end do

    if (.not. ctem_on) then
      FCANROT = ncGet3DVar(initid, 'FCAN', start = [lonIndex,latIndex,1,1], count = [1,1,icp1,nmos], format = [nlat,nmos,icp1])
      ! Error check:
      do i = 1,nlat
        do m = 1,nmos
          if (FAREROT(i,m) > 1.0) then
            print * ,'FAREROT > 1',FAREROT(I,M)
            call errorHandler('read_initialstate', - 1)
          end if
        end do
      end do
      ! else fcancmx is read in instead and fcanrot is derived later.
    end if

    ! Complete some initial set up work. The limiting snow
    ! depth, ZSNL, is assigned its operational value of 0.10 m.
    do I = 1,nlat ! loop 100
      do M = 1,nmos
        do J = 1,IGND
          TBARROT(I,M,J) = TBARROT(I,M,J) + TFREZ
        end do
        TSNOROT(I,M) = TSNOROT(I,M) + TFREZ
        TCANROT(I,M) = TCANROT(I,M) + TFREZ
        TPNDROT(I,M) = TPNDROT(I,M) + TFREZ
        TBASROT(I,M) = TBARROT(I,M,IGND)
        CMAIROT(I,M) = 0.
        WSNOROT(I,M) = 0.
        ZSNLROT(I,M) = 0.10
        TSFSROT(I,M,1) = TFREZ
        TSFSROT(I,M,2) = TFREZ
        TSFSROT(I,M,3) = TBARROT(I,M,1)
        TSFSROT(I,M,4) = TBARROT(I,M,1)
        TACROT (I,M) = TCANROT(I,M)
        QACROT (I,M) = 0.5E-2
      end do
    end do ! loop 100

    ! Set the counter for the number of iterations required to solve surface energy balance for the elements of the four subareas to zero.
    ITCTROT = 0

    ! Check that the THIC and THLQ values are set to zero for soil layers
    ! that are non-permeable (bedrock).
     do i = 1,nlat
      do j = 1,nmos
        if (FAREROT(i,j) > zero) then !  only operate on tiles in use.
          do m = 1,ignd - 1
            if (zbot(m) < SDEPROT(i,j) .and. zbot(m + 1) >= SDEPROT(i,j)) then ! if soil permeable depth lies in layer m+1
              if ( (m+2).le.ignd ) then ! and m+2, i.e. the next layer is not past the last layer
                if (any(SANDROT(i,j,m+2:ignd) > 0.) .or. any(CLAYROT(i,j,m+2:ignd) > 0.) .or. &
                    any(ORGMROT(i,j,m+2:ignd) > 0.) .or. any(THLQROT(i,j,m+2:ignd) > zero) .or. &
                    any(THICROT(i,j,m+2:ignd) > zero)) then
                    write(*,*)' '
                    write(*,*)'For ground layer ',m+2,' in tile ',j,' and onwards ...'
                    write(*,*)'setting SAND, CLAY, & ORGM to -3, and setting THIC and THLQ to zero,'
                    write(*,*)'since soil permeable depth SDEP ',SDEPROT(i,j),' lies in layer ',m+1
                    write(*,*)' '
                    SANDROT(i,j,m+2:ignd) = -3.0
                    CLAYROT(i,j,m+2:ignd) = -3.0
                    ORGMROT(i,j,m+2:ignd) = -3.0
                    THLQROT(i,j,m+2:ignd) = 0.
                    THICROT(i,j,m+2:ignd) = 0.
                end if 
              end if 
              exit
            end if
          end do
        else ! this tile is not in use, set the values to 0 so that the model won't do an NAN operations in soilProperties.
             ! the model will otherwise ignore this tile in other operations. 
          SANDROT(i,j,:) = 0.
          CLAYROT(i,j,:) = 0.
          ORGMROT(i,j,:) = 0.
          THLQROT(i,j,:) = 0.
          THICROT(i,j,:) = 0.
        end if 
      end do
     end do

    if (ctem_on) then

      grclarearow = ncGet1DVar(initid, 'grclarea', start = [lonIndex,latIndex], count = [1,1])

      slopefrac = ncGet3DVar(initid, 'slopefrac', start = [lonIndex,latIndex,1,1], count = [1,1,8,nmos], format = [nlat,nmos,8])
      Cmossmas = ncGet2DVar(initid, 'Cmossmas', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      litrmsmoss = ncGet2DVar(initid, 'litrmsmoss', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      dmoss = ncGet2DVar(initid, 'dmoss', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      peatSoilC = ncGet2DVar(initid, 'peatSoilC', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      fcancmxrow = ncGet3DVar(initid, 'fcancmx', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      gleafmas_NSrow = ncGet3DVar(initid, 'gleafmas_ns', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      gleafmassrow  = ncGet3DVar(initid, 'gleafmas_s', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      lygleafmasmaxrow = ncGet3DVar(initid, 'lygleafmasmax', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      bleafmasrow = ncGet3DVar(initid, 'bleafmas', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      stemmass_NSrow = ncGet3DVar(initid, 'stemmass_ns', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      stemmasssrow = ncGet3DVar(initid, 'stemmass_s', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      lystemmassmaxrow = ncGet3DVar(initid, 'lystemmassmax', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      rootmass_NSrow = ncGet3DVar(initid, 'rootmass_ns', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      rootmasssrow = ncGet3DVar(initid, 'rootmass_s', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      lyrootmassmaxrow = ncGet3DVar(initid, 'lyrootmassmax', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      grwtheffrow = ncGet3DVar(initid, 'grwtheff', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])      
      litrmassrow = ncGet4DVar(initid, 'litrmass', start = [lonIndex,latIndex,1,1,1], count = [1,1,iccp2,ignd,nmos], format = [nlat,nmos,iccp2,ignd])
      soilcmasrow = ncGet4DVar(initid, 'soilcmas', start = [lonIndex,latIndex,1,1,1], count = [1,1,iccp2,ignd,nmos], format = [nlat,nmos,iccp2,ignd])
      flhrloss_nsrow = ncGet3DVar(initid, 'flhrloss_ns', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      flhrloss_srow = ncGet3DVar(initid, 'flhrloss_s', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      stmhrlosrow = ncGet3DVar(initid, 'stmhrlos', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      rothrlosrow = ncGet3DVar(initid, 'rothrlos', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      lystmmasrow = ncGet3DVar(initid, 'lystmmas', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      lyrotmasrow = ncGet3DVar(initid, 'lyrotmas', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      tymaxlairow = ncGet3DVar(initid, 'tymaxlai', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      cfluxcgrow = ncGet2DVar(initid, 'cfluxcg', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      cfluxcsrow = ncGet2DVar(initid, 'cfluxcs', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      colddays_leaffallrow = ncGet2DVar(initid, 'colddays_leaffall', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      colddays_harvestrow = ncGet2DVar(initid, 'colddays_harvest', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      co2i1cgrow = ncGet3DVar(initid, 'co2i1cg', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      co2i1csrow = ncGet3DVar(initid, 'co2i1cs', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      co2i2cgrow = ncGet3DVar(initid, 'co2i2cg', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      co2i2csrow = ncGet3DVar(initid, 'co2i2cs', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      ! Since the green leaf mass, stem mass, and root mass are no longer read in, make
      ! the values here
      gleafmasrow = gleafmas_NSrow + gleafmassrow
      stemmassrow = stemmass_NSrow + stemmasssrow
      rootmassrow = rootmass_NSrow + rootmasssrow
      flhrlossrow = flhrloss_nsrow + flhrloss_srow

      ! If a tracer is being used, read in those values.
      if (useTracer > 0) then
        tracerGLeafMass = ncGet3DVar(initid, 'tracerGLeafMass', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
        tracerBLeafMass = ncGet3DVar(initid, 'tracerBLeafMass', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
        tracerStemMass = ncGet3DVar(initid, 'tracerStemMass', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
        tracerRootMass = ncGet3DVar(initid, 'tracerRootMass', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
        tracerLitrMass = ncGet4DVar(initid, 'tracerLitrMass', start = [lonIndex,latIndex,1,1,1], count = [1,1,iccp2,ignd,nmos], format = [nlat,nmos,iccp2,ignd])
        tracerSoilCMass = ncGet4DVar(initid, 'tracerSoilCMass', start = [lonIndex,latIndex,1,1,1], count = [1,1,iccp2,ignd,nmos], format = [nlat,nmos,iccp2,ignd])
        tracerMossCMass = ncGet2DVar(initid, 'tracerMossCMass', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
        tracerMossLitrMass = ncGet2DVar(initid, 'tracerMossLitrMass', start = [lonIndex,latIndex,1], count = [1,1,nmos], format = [nlat,nmos])
      end if

      lfstatusrow = ncGet3DVar(initid, 'lfstatus', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])
      pandaysrow = ncGet3DVar(initid, 'pandays', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])

      if (Ncycle_on) then
          !ngleafmasrow   = ncGet3DVar(initid, 'ngleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          ngleafmas_NSrow = ncGet3DVar(initid, 'ngleafmas_ns', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          ngleafmassrow  = ncGet3DVar(initid, 'ngleafmas_s', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          nbleafmasrow   = ncGet3DVar(initid, 'nbleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          !nstemmassrow   = ncGet3DVar(initid, 'nstemmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          nstemmass_NSrow = ncGet3DVar(initid, 'nstemmass_ns', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          nstemmasssrow  = ncGet3DVar(initid, 'nstemmass_s', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          !nrootmassrow   = ncGet3DVar(initid, 'nrootmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          nrootmass_NSrow = ncGet3DVar(initid, 'nrootmass_ns', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          nrootmasssrow  = ncGet3DVar(initid, 'nrootmass_s', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          bnf_natrow    = ncGet3DVar(initid, 'bnf_nat', start = [lonIndex, latIndex, 1, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          bnf_antrow    = ncGet3DVar(initid, 'bnf_ant', start = [lonIndex, latIndex, 1, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
          nh4_massrow    = ncGet3DVar(initid, 'nh4_mass', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp1, nmos], format = [nlat, nmos,iccp1])
          no3_massrow    = ncGet3DVar(initid, 'no3_mass', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp1, nmos], format = [nlat, nmos,iccp1])
          nlitrmassrow   = ncGet3DVar(initid, 'nlitrmass', start = [lonIndex, latIndex, 1, 1, 1], count = [1, 1, iccp2, nmos], format = [nlat, nmos, iccp2])
          soilnmasrow    = ncGet3DVar(initid, 'soilnmas', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp2, nmos], format = [nlat, nmos,iccp2])
          soilpHrow      = ncGet2DVar(initid, 'soilpH', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])

        ! Since the green leaf mass, stem mass, and root mass are no longer read in, make
        ! the values here 
        ngleafmasrow = ngleafmas_NSrow + ngleafmassrow
        nstemmassrow = nstemmass_NSrow + nstemmasssrow
        nrootmassrow = nrootmass_NSrow + nrootmasssrow
      end if

      if (PFTCompetition .and. inibioclim) then  ! read in the bioclimatic parameters

        twarmmrow = ncGet1DVar(initid, 'twarmm', start = [lonIndex,latIndex], count = [1,1])
        tcoldmrow = ncGet1DVar(initid, 'tcoldm', start = [lonIndex,latIndex], count = [1,1])
        gdd5row = ncGet1DVar(initid, 'gdd5', start = [lonIndex,latIndex], count = [1,1])
        aridityrow = ncGet1DVar(initid, 'aridity', start = [lonIndex,latIndex], count = [1,1])
        srplsmonrow = ncGet1DVar(initid, 'srplsmon', start = [lonIndex,latIndex], count = [1,1])
        defctmonrow = ncGet1DVar(initid, 'defctmon', start = [lonIndex,latIndex], count = [1,1])
        anndefctrow = ncGet1DVar(initid, 'anndefct', start = [lonIndex,latIndex], count = [1,1])
        annsrplsrow = ncGet1DVar(initid, 'annsrpls', start = [lonIndex,latIndex], count = [1,1])
        annpcprow = ncGet1DVar(initid, 'annpcp', start = [lonIndex,latIndex], count = [1,1])
        dry_season_lengthrow = ncGet1DVar(initid, 'dry_season_length', start = [lonIndex,latIndex], count = [1,1])
        nppvegrow = ncGet3DVar(initid, 'nppveg', start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos], format = [nlat,nmos,icc])

      else if (PFTCompetition .and. .not. inibioclim) then ! set them to zero

        twarmmrow = 0.0
        tcoldmrow = 0.0
        gdd5row = 0.0
        aridityrow = 0.0
        srplsmonrow = 0.0
        defctmonrow = 0.0
        anndefctrow = 0.0
        annsrplsrow = 0.0
        annpcprow = 0.0
        dry_season_lengthrow = 0.0

      end if

      !> if this run uses the competition and starts from bare ground, set up the model state here. this
      !> overwrites what was read in from the initialization file.

      if (PFTCompetition .and. start_bare) then

        ! If useTracer > 0 then the tracer values are left initialized at what they were read in as.
        do i = 1,nlat
          do m = 1,nmos
            do j = 1,icc
              if (.not. crop(j)) fcancmxrow(i,m,j) = 0.0
              gleafmasrow(i,m,j) = 0.0
              gleafmas_NSrow(i,m,j)= 0.0
              gleafmassrow(i,m,j) = 0.0
              bleafmasrow(i,m,j) = 0.0
              stemmassrow(i,m,j) = 0.0
              stemmass_NSrow(i,m,j)= 0.0
              stemmasssrow(i,m,j) = 0.0
              rootmassrow(i,m,j) = 0.0
              rootmass_NSrow(i,m,j)= 0.0
              rootmasssrow(i,m,j) = 0.0
              lfstatusrow(i,m,j) = 4
              pandaysrow(i,m,j) = 0
              if (Ncycle_on) then
                 ngleafmasrow(i,m,j) = 0.0
                 ngleafmas_NSrow(i,m,j) = 0.0
                 ngleafmassrow(i,m,j) = 0.0
                 nbleafmasrow(i,m,j) = 0.0
                 nstemmassrow(i,m,j) = 0.0
                 nstemmass_NSrow(i,m,j) = 0.0
                 nstemmasssrow(i,m,j) = 0.0
                 nrootmassrow(i,m,j) = 0.0
                 nrootmass_NSrow(i,m,j) = 0.0
                 nrootmasssrow(i,m,j) = 0.0
                 bnf_natrow(i,m,j) = 0.0
                 bnf_antrow(i,m,j) = 0.0
              end if
            end do

            lfstatusrow(i,m,1) = 2

            do j = 1,iccp2
              litrmassrow(i,m,j,1:ignd)=0.0
              soilcmasrow(i,m,j,1:ignd)=0.0
               if (Ncycle_on) then
                  nlitrmassrow(i,m,j)= 0.0
                  soilnmasrow(i,m,j) = 0.0
               end if
            end do

            do j = 1,iccp1
              if (Ncycle_on) then
                  nh4_massrow(i,m,j) = 0.0
                  no3_massrow(i,m,j) = 0.0
              end if
            end do

          end do ! nmtest
        end do ! nltest

      end if ! if (PFTCompetition .and. start_bare)

      !> If fire and competition are on, save the stemmass and rootmass for use in burntobare subroutine on the first timestep.
      if (dofire .and. PFTCompetition) then
        do i = 1,nlat
          do m = 1,nmos
            do j = 1,icc
              pstemmassrow(i,m,j) = stemmassrow(i,m,j)
              pgleafmassrow(i,m,j) = rootmassrow(i,m,j)
            end do
          end do
        end do
      end if

    end if ! ctem_on

    end associate
  end subroutine read_initialstate

  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_write_restart
  !! @{
  !> Write out the model restart file to netcdf. We only write out the variables that the model
  !! influences. This overwrites a pre-existing netcdf file.
  !> @author Joe Melton

  subroutine write_restart (lonIndex, latIndex)

    use ctemStateVars,  only : c_switch, vrot, tracer
    use classStateVars, only : class_rot
    use classicParams,  only : icc, nmos, ignd, icp1, modelpft, iccp2, TFREZ, iccp1

    implicit none

    ! arguments
    integer, intent(in) :: lonIndex, latIndex

    ! local
    integer :: k,m

    ! Associate names with variables defined in derived types.
    associate( &
    ctem_on           => c_switch%ctem_on,              & !< logical: 
    Ncycle_on         => c_switch%Ncycle_on,            & !< logical: 
    fertilizeron      => c_switch%fertilizeron,         & !< logical: 
    depositionon      => c_switch%depositionon,         & !< logical: 
    PFTCompetition    => c_switch%PFTCompetition,       & !< logical: 
    lnduseon          => c_switch%lnduseon,             & !< logical: 
    useTracer         => c_switch%useTracer,            & !< integer: 
    dynamicTilingOn     => c_switch%dynamicTilingOn,        & !< logical:
    trackTileAge      => c_switch%trackTileAge,    & !< logical:
    ipeatlandrow      => vrot%ipeatland,                & !< integer, dimension(:,:) : Peatland switch: 0 = not a peatland, 1 = bog, 2 = fen 
    fcancmxrow        => vrot%fcancmx,                  & !< real, dimension(:,:,:) :! 
    gleafmasrow       => vrot%gleafmas,                 & !< real, dimension(:,:,:) :! 
    gleafmas_NSrow    => vrot%gleafmas_ns,             & !< real, dimension(:,:,:) :! 
    gleafmassrow      => vrot%gleafmas_s,               & !< real, dimension(:,:,:) :! 
    bleafmasrow       => vrot%bleafmas,                 & !< real, dimension(:,:,:) :! 
    stemmassrow       => vrot%stemmass,                 & !< real, dimension(:,:,:) :! 
    stemmass_NSrow    => vrot%stemmass_ns,             & !< real, dimension(:,:,:) :! 
    stemmasssrow      => vrot%stemmass_s,               & !< real, dimension(:,:,:) :! 
    rootmassrow       => vrot%rootmass,                 & !< real, dimension(:,:,:) :! 
    rootmass_NSrow    => vrot%rootmass_ns,             & !< real, dimension(:,:,:) :! 
    rootmasssrow      => vrot%rootmass_s,               & !< real, dimension(:,:,:) :!
    nppvegrow         => vrot%nppveg,                 & !< real, dimension(:,:,:) :! 
    flhrlossrow       => vrot%flhrloss,                 & !< fall or harvest loss for deciduous trees and crops, respectively, \f$(kg C/m^2)\f$
    flhrloss_nsrow    => vrot%flhrloss_ns,              & !< fall or harvest loss for deciduous trees and crops, respectively, \f$(kg C/m^2)\f$
    flhrloss_srow     => vrot%flhrloss_s,               & !< fall or harvest loss for deciduous trees and crops, respectively, \f$(kg C/m^2)\f$
    stmhrlosrow       => vrot%stmhrlos,                 & !< real, dimension(:,:) : 
    rothrlosrow       => vrot%rothrlos,                 & !< real, dimension(:,:) : 
    lystmmasrow       => vrot%lystmmas,                 & !< real, dimension(:,:,:) : 
    lyrotmasrow       => vrot%lyrotmas,                 & !< real, dimension(:,:,:) : 
    tymaxlairow       => vrot%tymaxlai,                 & !< real, dimension(:,:,:) :
    cfluxcgrow        => vrot%cfluxcg,                  & !< real, dimension(:,:) :
    cfluxcsrow        => vrot%cfluxcs,                  & !< real, dimension(:,:) :
    colddays_leaffallrow => vrot%colddays_leaffall,     & !< integer, dimension(:,:) : 
    colddays_harvestrow  => vrot%colddays_harvest,      & !< integer, dimension(:,:) : 
    co2i1cgrow        => vrot%co2i1cg,                  & !< real, dimension(:,:,:) :
    co2i1csrow        => vrot%co2i1cs,                  & !< real, dimension(:,:,:) :
    co2i2cgrow        => vrot%co2i2cg,                  & !< real, dimension(:,:,:) :
    co2i2csrow        => vrot%co2i2cs,                  & !< real, dimension(:,:,:) :

    lygleafmasmaxrow  => vrot%lygleafmasmax,            & !< real, dimension(:,:,:) : last year maximum of the green leaf mass for each PFT, \f$kg c/m^2\f$ 
    lystemmassmaxrow  => vrot%lystemmassmax,            & !< real, dimension(:,:,:) : last year maximum of the stem mass for each PFT, \f$kg c/m^2\f$ 
    lyrootmassmaxrow  => vrot%lyrootmassmax,            & !< real, dimension(:,:,:) : last year maximum of the root mass for each PFT, \f$kg c/m^2\f$ 
    twarmm            => vrot%twarmm,                   & !< real, dimension(:,:) : temperature of the warmest month (c) 
    tcoldm            => vrot%tcoldm,                   & !< real, dimension(:,:) : temperature of the coldest month (c) 
    gdd5              => vrot%gdd5,                     & !< real, dimension(:,:) : growing degree days above 5 c 
    aridity           => vrot%aridity,                  & !< real, dimension(:,:) : aridity index, ratio of potential evaporation to precipitation 
    srplsmon          => vrot%srplsmon,                 & !< real, dimension(:,:) : number of months in a year with surplus water i.e.precipitation more than potential evaporation 
    defctmon          => vrot%defctmon,                 & !< real, dimension(:,:) : number of months in a year with water deficit i.e.precipitation less than potential evaporation 
    anndefct          => vrot%anndefct,                 & !< real, dimension(:,:) : annual water deficit (mm) 
    annsrpls          => vrot%annsrpls,                 & !< real, dimension(:,:) : annual water surplus (mm) 
    annpcp            => vrot%annpcp,                   & !< real, dimension(:,:) : annual precipitation (mm) 
    dry_season_length => vrot%dry_season_length,        & !< real, dimension(:,:) : length of dry season (months) 
    litrmassrow       => vrot%litrmass,                 & !< real, dimension(:,:,:,:) : 
    soilcmasrow       => vrot%soilcmas,                 & !< real, dimension(:,:,:,:) : 
    grwtheffrow       => vrot%grwtheff,                 & !< real, dimension(:,:,:) : growth efficiency. change in biomass per year per unit max.
                                                          !!                          lai (\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine 
    lfstatusrow       => vrot%lfstatus,                 & !< integer, dimension(:,:,:) : 
    pandaysrow        => vrot%pandays,                  & !< integer, dimension(:,:,:) : 
    Cmossmas          => vrot%Cmossmas,                 & !< real, dimension(:,:) : C in moss biomass, \f$kg C/m^2\f$ 
    litrmsmoss        => vrot%litrmsmoss,               & !< real, dimension(:,:) : moss litter mass, \f$kg C/m^2\f$ 
    dmoss             => vrot%dmoss,                    & !< real, dimension(:,:) : depth of living moss (m) 
    peatSoilC         => vrot%peatSoilC,                & !< real, dimension(:,:) : peat soil C mass, \f$kg C/m^2\f$ 
    tracerGLeafMass   => tracer%gLeafMassrot,           & !< real, dimension(:,:,:) : Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerBLeafMass   => tracer%bLeafMassrot,           & !< real, dimension(:,:,:) : Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerStemMass    => tracer%stemMassrot,            & !< real, dimension(:,:,:) : Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerRootMass    => tracer%rootMassrot,            & !< real, dimension(:,:,:) : Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerLitrMass    => tracer%litrMassrot,            & !< real, dimension(:,:,:,:) : Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    tracerSoilCMass   => tracer%soilCMassrot,           & !< real, dimension(:,:,:,:) : Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    tracerMossCMass   => tracer%mossCMassrot,           & !< real, dimension(:,:) : Tracer mass in moss biomass, \f$kg C/m^2\f$ 
    tracerMossLitrMass => tracer%mossLitrMassrot,       & !< real, dimension(:,:) : Tracer mass in moss litter, \f$kg C/m^2\f$ 
    nh4_massrow       => vrot%nh4_mass,                 & !< real, dimension(:,:,:) : ammonium mass (\f$g N/m^2\f$) 
    no3_massrow       => vrot%no3_mass,                 & !< real, dimension(:,:,:) : nitrate mass (\f$g N/m^2\f$) 
    ngleafmasrow      => vrot%ngleafmas,                & !< real, dimension(:,:,:) : 
    ngleafmas_NSrow    => vrot%ngleafmas_ns,            & !< real, dimension(:,:,:) : 
    ngleafmassrow     => vrot%ngleafmas_s,              & !< real, dimension(:,:,:) : 
    nbleafmasrow      => vrot%nbleafmas,                & !< real, dimension(:,:,:) : 
    nstemmassrow      => vrot%nstemmass,                & !< real, dimension(:,:,:) : 
    nstemmass_NSrow    => vrot%nstemmass_ns,            & !< real, dimension(:,:,:) : 
    nstemmasssrow     => vrot%nstemmass_s,              & !< real, dimension(:,:,:) : 
    nrootmassrow      => vrot%nrootmass,                & !< real, dimension(:,:,:) : 
    nrootmass_NSrow    => vrot%nrootmass_ns,            & !< real, dimension(:,:,:) : 
    nrootmasssrow     => vrot%nrootmass_s,              & !< real, dimension(:,:,:) : 
    bnf_natrow        => vrot%bnf_nat,                  & !< real, dimension(:,:,:) :
    bnf_antrow        => vrot%bnf_ant,                  & !< real, dimension(:,:,:) :
    nlitrmassrow      => vrot%nlitrmass,                & !< real, dimension(:,:,:) : 
    soilnmasrow       => vrot%soilnmas,                 & !< real, dimension(:,:,:) : 
    controlVector     => vrot%controlVector,            & !< integer, dimension(nlat,nmos)
    tileAgerow        => vrot%tileAgerow,               & !< real, dimension(nlat,nmos)
    FCANROT           => class_rot%FCANROT,             & !< real, dimension(:,:,:) : 
    FAREROT           => class_rot%FAREROT,             & !< real, dimension(:,:)   : 
    TBARROT           => class_rot%TBARROT,             & !< real(r8), dimension(:,:,:) : 
    THLQROT           => class_rot%THLQROT,             & !< real, dimension(:,:,:) : 
    THICROT           => class_rot%THICROT,             & !< real, dimension(:,:,:) : 
    TCANROT           => class_rot%TCANROT,             & !< real, dimension(:,:)   : 
    TSNOROT           => class_rot%TSNOROT,             & !< real, dimension(:,:)   : 
    TPNDROT           => class_rot%TPNDROT,             & !< real, dimension(:,:)   : 
    ZPNDROT           => class_rot%ZPNDROT,             & !< real, dimension(:,:)   : 
    RCANROT           => class_rot%RCANROT,             & !< real, dimension(:,:)   : 
    SCANROT           => class_rot%SCANROT,             & !< real, dimension(:,:)   : 
    SNOROT            => class_rot%SNOROT,              & !< real, dimension(:,:)   : 
    ALBSROT           => class_rot%ALBSROT,             & !< real, dimension(:,:)   : 
    RHOSROT           => class_rot%RHOSROT,             & !< real, dimension(:,:)   : 
    GROROT            => class_rot%GROROT,              & !< real, dimension(:,:)   : 
    maxAnnualActLyr   => class_rot%maxAnnualActLyrROT,  & !< real, dimension(:,:)   : Active layer depth maximum over the e-folding period specified by parameter eftime (m). 
    SDEPROT           => class_rot%SDEPROT              & !< real, dimension(:,:)   : 
    )

    call ncPut2DVar(rsid, 'FARE', FAREROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut3DVar(rsid, 'FCAN', FCANROT,start = [lonIndex,latIndex,1,1], count = [1,1,icp1,nmos])
    call ncPut3DVar(rsid, 'THLQ', THLQROT,start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos])
    call ncPut3DVar(rsid, 'THIC', THICROT,start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos])

    !call ncPut3DVar(rsid, 'TBAR', TBARROT - TFREZ,start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos])
    ! Convert 64-bit precision TBARROT to default real type in writeOutput1D.
    ! FLAG, EC: This means TBARROT will be output as 32-bit reals if compiled in 32-bit mode.
    !           Need to change if want to always output 64-bit reals instead.
    call ncPut3DVar(rsid, 'TBAR', real(TBARROT) - TFREZ,start = [lonIndex,latIndex,1,1], count = [1,1,ignd,nmos])

    call ncPut2DVar(rsid, 'TCAN', TCANROT - TFREZ,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'TSNO', TSNOROT - TFREZ,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'TPND', TPNDROT - TFREZ,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'ZPND', ZPNDROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'RCAN', RCANROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'SCAN', SCANROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'SNO', SNOROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'ALBS', ALBSROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'RHOS', RHOSROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'GRO', GROROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    call ncPut2DVar(rsid, 'maxAnnualActLyr', maxAnnualActLyr,start = [lonIndex,latIndex,1], count = [1,1,nmos])

    ! If the peatland module is on for bogs or fens, the SDEP reflects the peat depth so can change throughout a run. 
    ! as a result, the SDEP needs to be written back to the restart file. If this is the case, for simplicity,
    ! write the SDEP for all tiles in the grid cell. 
    if (any(ipeatlandrow == 1 .or. ipeatlandrow == 2)) call ncPut2DVar(rsid, 'SDEP', SDEPROT,start = [lonIndex,latIndex,1], count = [1,1,nmos])
    
    if (ctem_on) then
      call ncPut3DVar(rsid, 'fcancmx', fcancmxrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      !call ncPut3DVar(rsid, 'gleafmas', gleafmasrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'lygleafmasmax', lygleafmasmaxrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'gleafmas_ns', gleafmas_NSrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'gleafmas_s', gleafmassrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'bleafmas', bleafmasrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      !call ncPut3DVar(rsid, 'stemmass', stemmassrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'lystemmassmax', lystemmassmaxrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'stemmass_ns', stemmass_NSrow, start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'stemmass_s', stemmasssrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      !call ncPut3DVar(rsid, 'rootmass', rootmassrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'lyrootmassmax', lyrootmassmaxrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'rootmass_ns', rootmass_NSrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'rootmass_s', rootmasssrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      !call ncPut3DVar(rsid, 'flhrloss', flhrlossrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'flhrloss_ns', flhrloss_nsrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'flhrloss_s', flhrloss_srow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'stmhrlos', stmhrlosrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'rothrlos', rothrlosrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'lystmmas', lystmmasrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'lyrotmas', lyrotmasrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'tymaxlai', tymaxlairow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'grwtheff', grwtheffrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut2DVar(rsid, 'cfluxcg', cfluxcgrow,start = [lonIndex,latIndex,1], count = [1,1,nmos])
      call ncPut2DVar(rsid, 'cfluxcs', cfluxcsrow,start = [lonIndex,latIndex,1], count = [1,1,nmos])
      call ncPut2DVar(rsid, 'colddays_leaffall', real(colddays_leaffallrow),start = [lonIndex,latIndex,1], count = [1,1,nmos])
      call ncPut2DVar(rsid, 'colddays_harvest', real(colddays_harvestrow),start = [lonIndex,latIndex,1], count = [1,1,nmos])
      call ncPut3DVar(rsid, 'co2i1cg', co2i1cgrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'co2i1cs', co2i1csrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'co2i2cg', co2i2cgrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'co2i2cs', co2i2csrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])

      do m = 1,nmos !s.r.c added outer nmos loop to ensure tiles are written in correct order       
        do k = 1,ignd
          call ncPut2DVar(rsid, 'litrmass', litrmassrow(:,m,:,k),start = [lonIndex,latIndex,1,k,m], count = [1,1,iccp2,1,1])
          call ncPut2DVar(rsid, 'soilcmas', soilcmasrow(:,m,:,k),start = [lonIndex,latIndex,1,k,m], count = [1,1,iccp2,1,1])
        end do
      end do

      call ncPut3DVar(rsid, 'lfstatus', real(lfstatusrow),start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut3DVar(rsid, 'pandays', real(pandaysrow),start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
      call ncPut2DVar(rsid, 'Cmossmas', Cmossmas,start = [lonIndex,latIndex,1], count = [1,1,nmos])
      call ncPut2DVar(rsid, 'litrmsmoss', litrmsmoss,start = [lonIndex,latIndex,1], count = [1,1,nmos])
      
      call ncPut2DVar(rsid, 'peatSoilC', peatSoilC,start = [lonIndex,latIndex,1], count = [1,1,nmos])
      call ncPut2DVar(rsid, 'dmoss', dmoss,start = [lonIndex,latIndex,1], count = [1,1,nmos])

      ! If a tracer is being used,read in those values.
      if (useTracer > 0) then
        call ncPut3DVar(rsid, 'tracerGLeafMass', tracerGLeafMass,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
        call ncPut3DVar(rsid, 'tracerBLeafMass', tracerBLeafMass,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
        call ncPut3DVar(rsid, 'tracerStemMass', tracerStemMass,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
        call ncPut3DVar(rsid, 'tracerRootMass', tracerRootMass,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])
        do m = 1,nmos !s.r.c added outer nmos loop to ensure tiles are written in correct order
          do k = 1,ignd
            call ncPut2DVar(rsid, 'tracerLitrMass', tracerLitrMass(:,m,:,k),start = [lonIndex,latIndex,1,k,m], count = [1,1,iccp2,1,1])
            call ncPut2DVar(rsid, 'tracerSoilCMass', tracerSoilCMass(:,m,:,k),start = [lonIndex,latIndex,1,k,m], count = [1,1,iccp2,1,1])
          end do
        end do
        call ncPut2DVar(rsid, 'tracerMossCMass', tracerMossCMass,start = [lonIndex,latIndex,1], count = [1,1,nmos])
        call ncPut2DVar(rsid, 'tracerMossLitrMass', tracerMossLitrMass,start = [lonIndex,latIndex,1], count = [1,1,nmos])
      end if

      if (Ncycle_on) then
        call ncPut3DVar(rsid, 'nh4_mass', nh4_massrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, iccp1, nmos])
        call ncPut3DVar(rsid, 'no3_mass', no3_massrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, iccp1, nmos])
        !call ncPut3DVar(rsid, 'ngleafmas', ngleafmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'ngleafmas_ns', ngleafmas_NSrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'ngleafmas_s', ngleafmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'nbleafmas', nbleafmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        !call ncPut3DVar(rsid, 'nstemmass', nstemmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'nstemmass_ns', nstemmass_NSrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'nstemmass_s', nstemmasssrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        !call ncPut3DVar(rsid, 'nrootmass', nrootmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'nrootmass_ns', nrootmass_NSrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'nrootmass_s', nrootmasssrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'bnf_nat', bnf_natrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'bnf_ant', bnf_antrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
        call ncPut3DVar(rsid, 'nlitrmass', nlitrmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, iccp2, nmos])
        call ncPut3DVar(rsid, 'soilnmas', soilnmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, iccp2, nmos])
      end if

      if (PFTCompetition) then

        ! Since these climate related variables are only sensible at the gridcell level,we just write out the
        ! value for the first tile (nlat is always 1 offline too).
        call ncPutVar(rsid, 'twarmm', realValues = twarmm, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'tcoldm', realValues = tcoldm, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'gdd5', realValues = gdd5, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'aridity', realValues = aridity, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'srplsmon', realValues = srplsmon, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'defctmon', realValues = defctmon, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'anndefct', realValues = anndefct, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'annsrpls', realValues = annsrpls, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'annpcp', realValues = annpcp, start = [lonIndex,latIndex], count = [1,1])
        call ncPutVar(rsid, 'dry_season_length', realValues = dry_season_length, start = [lonIndex,latIndex], count = [1,1])
        call ncPut3DVar(rsid, 'nppveg', nppvegrow,start = [lonIndex,latIndex,1,1], count = [1,1,icc,nmos])

      end if ! PFTCompetition

      if(dynamicTilingOn) call ncPut2DVar(rsid, 'control_vector', real(controlVector), start = [lonIndex,latIndex,1], count = [1,1,nmos])
      if(trackTileAge .or. dynamicTilingOn) call ncPut2DVar(rsid, 'tile_age', tileAgerow, start = [lonIndex,latIndex,1], count = [1,1,nmos])
    end if ! ctem_on

    end associate
  end subroutine write_restart

  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_getInput
  !! @{
  !>  Read in a model input from a netcdf file and store the file's time array
  !! as well as the input values into memory.
  !> @author Joe Melton, Ed Chan

  subroutine getInput (inputRequested, longitude, latitude, projLonInd, projLatInd)

    use fileIOModule
    use generalUtils,  only : findLeapYears
    use ctemStateVars, only : c_switch, vrot, tracer
    use classicParams, only : icc, nmos, nsmu, nsalb, nbc, nreff, nswe, nbnd_lut, &
                              albdif_lut, albdir_lut, trandif_lut, trandir_lut
    use tracerModule, only : convertTracerUnits

    implicit none

    character( * ), intent(in) :: inputRequested
    real, intent(in), optional :: longitude
    real, intent(in), optional :: latitude
    integer, intent(in), optional :: projLonInd
    integer, intent(in), optional :: projLatInd

    integer :: lengthOfFile
    integer :: lonloc, latloc
    integer :: i, arrindex, arrindex2, ntimes, m, numPFTsinFile, d, j
    real(r8), dimension(:), allocatable :: fileTime
    !real, dimension(5) :: dateTime
    real(r8) :: startTime, endTime
    logical :: dummyVar
    integer :: lastDOY
    integer :: ipnt, albdim
    real, dimension(:,:,:), allocatable :: tmpalb
    integer :: isalb, ismu, isgs, iswe, ibc

    associate( &
    readMetStartYear => c_switch%readMetStartYear,      & !< integer: First year of meteorological forcing to read in from the met file 
    readMetEndYear   => c_switch%readMetEndYear,        & !< integer: Last year of meteorological forcing to read in from the met file 
    projectedGrid   => c_switch%projectedGrid,          & !< logical: 
    transientCO2    => c_switch%transientCO2,           & !< logical: 
    fixedYearCO2    => c_switch%fixedYearCO2,           & !< integer: 
    transientCH4    => c_switch%transientCH4,           & !< logical: 
    fixedYearCH4    => c_switch%fixedYearCH4,           & !< integer: 
    transientPOPD   => c_switch%transientPOPD,          & !< logical: 
    fixedYearPOPD   => c_switch%fixedYearPOPD,          & !< integer: 
    transientLGHT   => c_switch%transientLGHT,          & !< logical: 
    fixedYearLGHT   => c_switch%fixedYearLGHT,          & !< integer: 
    transientOBSWETF=> c_switch%transientOBSWETF,       & !< logical: 
    fixedYearOBSWETF=> c_switch%fixedYearOBSWETF,       & !< integer: 
    lnduseon        => c_switch%lnduseon,               & !< logical: 
    fixedYearLUC    => c_switch%fixedYearLUC,           & !< integer: 
    Ncycle_on       => c_switch%Ncycle_on,              & !< logical: 
    fertilizeron    => c_switch%fertilizeron,           & !< logical: 
    fixedYearFER    => c_switch%fixedYearFER,           & !< integer: 
    transientFER    => c_switch%transientFER,           & !< logical: 
    depositionon    => c_switch%depositionon,           & !< logical: 
    fixedYearDEP    => c_switch%fixedYearDEP,           & !< integer: 
    transientDEP    => c_switch%transientDEP,           & !< logical: 
    blackCtransientDep => c_switch%blackCtransientDep,  & ! logical: If true, allow time-varying black carbon deposition flux in snow processes (4-band scheme)
    fixedYearBCDep  => c_switch%fixedYearBCDep,         & ! integer: set the year to use for black carbon deposition flux if blackCtransientDep is false.
    leap            => c_switch%leap,                   & !< logical: 
    useTracer       => c_switch%useTracer,              & !< integer: useTracer = 0, the tracer code is not used.
                                                          !!          useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the
                                                          !!                        tracer values in the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                                          !!          useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                                          !!          useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme. 
    metLoop         => c_switch%metLoop,                & !< integer: 
    co2concrow      => vrot%co2conc,                    & !< real, dimension(:,:) : 
    tracerco2conc   => tracer%tracerCO2rot,             & !< real, dimension(:,:) : 
    ch4concrow      => vrot%ch4conc,                    & !< real, dimension(:,:) : 
    popdinrow       => vrot%popdin,                     & !< real, dimension(:,:) : 
    nfertilrow      => vrot%nfertil,                    & !< real, dimension(:,:) : 
    ndepositrow     => vrot%ndeposit,                   & !< real, dimension(:,:) : 
    timharvrow      => vrot%timharvrow,                  & !< real, dimension(:) :
    prsfirerow      => vrot%prsfirerow,                 & !< real, prsfirerow(:) :
    fcancmxrow      => vrot%fcancmx                     & !< real, dimension(:,:,:) : 
    )

    select case (trim(inputRequested))

      !> For each of the time varying inputs in this subroutine, we take in the whole dataset
      !! and later determine the year we need (in updateInput). The general approach is that these
      !! files are light enough on memory demands to make this acceptable.

      !! It is important that the files have time as the fastest varying dimension.

    case ('CO2') ! Carbon dioxide concentration

      lengthOfFile = ncGetDimLen(co2id, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(CO2id, 'time', start = 1, count = lengthOfFile)

      ! Parse these into just years (expected format is "day as %Y%m%d.%f")
      fileTime=aint(fileTime/10000)

      if (transientCO2) then

        ! Find the requested years in the file.
        arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The CO2 file does not contain first requested year')

        ! Sometimes it is correct to have transient CO2 but otherwise have constant conditiions (recycling MET),
        ! in this case metloop is >1 but transientCO2 is true. So grab the full length of the CO2 file rather than only
        ! the years requested for the met.
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
        else
          arrindex2 = lengthOfFile
        end if

        if (arrindex2 == 0) stop ('getInput says: The CO2 file does not contain last requested year')
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.

        allocate(CO2FromFile(ntimes))
        CO2FromFile = ncGet1DVar(CO2id,trim(co2VarName),start = [arrindex], count = [ntimes])

        allocate(CO2Time(ntimes))
        CO2Time = int(fileTime(arrindex:arrindex2))
      else
        ! Find the requested year in the file.
        arrindex = fndloc(fileTime,real(fixedYearCO2,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The CO2 file does not contain requested year')

        ! We read in only the suggested year
        ! Presently all cells get the exact same CO2 (we don't do it by cell) so 
        ! just set all the same now.
        do j = 1, nmos
          co2concrow(:,j) = ncGet1DVar(CO2id,trim(co2VarName),start = [arrindex], count = [1])
        end do
      end if

    case ('tracerCO2') ! tracer Carbon dioxide atmospheric values.

      lengthOfFile = ncGetDimLen(tracerco2id, 'time')
      allocate(fileTime(lengthOfFile))
      allocate(tracerCO2Time(lengthOfFile))

      fileTime = ncGetTime(tracerCO2id, 'time', start = 1, count = lengthOfFile)

      if (useTracer == 2) then
        ! 14C has different values depending on latitudional bands. The expected
        ! input file is from CMIP6 and it splits the bands as follows:
        ! Southern Hemisphere (30-90°S), Tropics (30°S-30°N),
        ! and Northern Hemisphere (30-90°N). We now assign the file
        ! variable name here
        if (latitude > 30.) then
          tracerco2VarName = 'NH_D14C'
        else if (latitude <= 30. .and. latitude >= - 30.) then
          tracerco2VarName = 'Tropics_D14C'
        else if (latitude < - 30.) then
          tracerco2VarName = 'SH_D14C'
        end if
      end if
      
      ! Parse these into just years (expected format is "day as %Y%m%d.%f")
      tracerCO2Time = aint(fileTime)/10000

      if (transientCO2) then
        ! We read in the whole CO2 times series and store it.
        allocate(tracerCO2FromFile(lengthOfFile))
        tracerCO2FromFile = ncGet1DVar(tracerCO2id,trim(tracerco2VarName),start = [1], count = [lengthOfFile])

      else
        ! Find the requested year in the file.
        arrindex = fndloc(tracerCO2Time,readMetStartYear,1)
        if (arrindex == 0) stop ('getInput says: The tracer CO2 file does not contain requested year')

        ! We read in only the suggested year
        ! Presently all cells get the exact same tracerCO2 (we don't do it by cell) so 
        ! just set all the same now.
        do j = 1, nmos
          tracerco2conc(:,j) = ncGet1DVar(tracerco2id,trim(tracerco2VarName),start = [arrindex], count = [1])
        end do 
      end if

      ! Convert the units of the tracer depending on the tracer being simulated.
      tracerco2conc = convertTracerUnits(tracerco2conc)

    case ('CH4') ! Methane concentration

      lengthOfFile = ncGetDimLen(ch4id, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(ch4id, 'time', start = 1, count = lengthOfFile)

      ! Parse these into just years (expected format is "day as %Y%m%d.%f")
      fileTime=aint(fileTime/10000)

      if (transientCH4) then
        ! Find the requested years in the file.
        arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The CH4 file does not contain first requested year')
        ! Sometimes it is correct to have transient CH4 but otherwise have constant conditiions (recycling MET),
        ! in this case metloop is >1 but transientCH4 is true. So grab the full length of the CH4 file rather than only
        ! the years requested for the met.
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
        else
          arrindex2 = lengthOfFile
        end if

        if (arrindex2 == 0) stop ('getInput says: The CH4 file does not contain last requested year')
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.

        allocate(CH4FromFile(ntimes))
        CH4FromFile = ncGet1DVar(ch4id,trim(ch4VarName),start = [arrindex], count = [ntimes])

        allocate(CH4Time(ntimes))
        CH4Time = int(fileTime(arrindex:arrindex2))
      else
        ! Find the requested year in the file.
        arrindex = fndloc(fileTime,real(fixedYearCH4,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The CH4 file does not contain requested year')

        ! We read in only the suggested year
        ! Presently all cells get the exact same CH4 (we don't do it by cell) so 
        ! just set all the same now.
        do j = 1, nmos
          ch4concrow(:,j) = ncGet1DVar(ch4id,trim(ch4VarName),start = [arrindex], count = [1])
        end do 
      end if

    case ('POPD') ! Population density

      lengthOfFile = ncGetDimLen(popid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(popid, 'time', start = 1, count = lengthOfFile)

      ! Parse these into just years (expected format is "day as %Y%m%d.%f")
      fileTime=aint(fileTime/10000)

      if (.not. projectedGrid) then
        lonloc = closestCell(popid,'lon',longitude)
        latloc = closestCell(popid,'lat',latitude)
      else
        ! For projected grids,we use the index of the cells,not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if

      if (transientPOPD) then
        ! Find the requested years in the file.
        arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The POPD file does not contain first requested year')
        ! Sometimes it is correct to have transient POPD but otherwise have constant conditiions (recycling MET),
        ! in this case metloop is >1 but transientPOPD is true. So grab the full length of the POPD file rather than only
        ! the years requested for the met.
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
        else
          arrindex2 = lengthOfFile
        end if

        if (arrindex2 == 0) stop ('getInput says: The POPD file does not contain last requested year')
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.

        allocate(POPDFromFile(ntimes))
        POPDFromFile = ncGet1DVar(popid,trim(popVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])

        allocate(POPDTime(ntimes))
        POPDTime = int(fileTime(arrindex:arrindex2))
      else
        ! Find the requested year in the file.
        arrindex = fndloc(fileTime,real(fixedYearPOPD,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The POPD file does not contain requested year')

        ! We read in only the suggested year
        i = 1 ! offline nlat is always 1 so just set
        popdinrow(i,:) = ncGet1DVar(popid,trim(popVarName),start = [lonloc,latloc,arrindex], count = [1,1,1])

      end if

    case ('LGHT') ! Lightning strikes

      lengthOfFile = ncGetDimLen(lghtid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(lghtid, 'time', start = 1, count = lengthOfFile)

      ! The lightning file is daily (expected format is "day as %Y%m%d.%f")
      ! We want to retain all except the partial day.
      fileTime=aint(fileTime)

      if (.not. projectedGrid) then
        lonloc = closestCell(lghtid,'lon',longitude)
        latloc = closestCell(lghtid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if

      ! Units expected are "strikes km-2 yr-1"

      if (transientLGHT) then
        ! Find the beginning and end day in the file.
        ! Assume we are grabbing from first day of start year to last day of last year.

        startTime = real(readMetStartYear,kind=r8) * 10000. + 1. * 100. + 1.
        endTime = real(readMetEndYear,kind=r8) * 10000. + 12. * 100. + 31.

        arrindex = fndloc(fileTime,startTime,1)
        if (arrindex == 0) stop ('getInput says: The LGHT file does not contain first requested day')
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,endTime,1)
          if (arrindex2 == 0) stop ('getInput says: The LGHT file does not contain last requested day')
        else
          arrindex2 = lengthOfFile
        end if
        
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.

        allocate(LGHTFromFile(ntimes))
        LGHTFromFile = ncGet1DVar(lghtid,trim(lghtVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])

        allocate(LGHTTime(ntimes))
        LGHTTime = fileTime(arrindex:arrindex2)
      else
        ! Find the requested day and year in the file.
        ! Assume we are grabbing from day 1
        startTime = real(fixedYearLGHT,kind=r8) * 10000. + 1. * 100. + 1.

        arrindex = fndloc(fileTime,startTime,1)
        if (arrindex == 0) stop ('getInput says: The LGHT file does not contain requested year')

        ! We read in only the suggested year of daily inputs

        ! If we are using leap years, check if that year is a leap year
        if (leap) then
          call findLeapYears(fixedYearLGHT, dummyVar, lastDOY)
        else 
          lastDOY = 365
        end if 

        allocate(LGHTFromFile(lastDOY))
        LGHTFromFile = ncGet1DVar(lghtid,trim(lghtVarName),start = [lonloc,latloc,arrindex], count = [1,1,lastDOY])

        ! Lastly, remake the LGHTTime to be only counting for one year for simplicity
        allocate(LGHTTime(lastDOY))
        do d = 1,lastDOY
          LGHTTime(d) = d
        end do

      end if

    case ('LUC') ! Land use change

      lengthOfFile = ncGetDimLen(lucid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(lucid, 'time', start = 1, count = lengthOfFile)

      ! Parse these into just years (expected format is "day as %Y%m%d.%f")
      fileTime=aint(fileTime/10000)

      if (.not. projectedGrid) then
        lonloc = closestCell(lucid,'lon',longitude)
        latloc = closestCell(lucid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if

      ! Ensure the file has the expected number of PFTs
      numPFTsinFile = ncGetDimLen(lucid, 'lev')
      if (numPFTsinFile /= icc) stop ('getInput says: LUC file does not have expected number of PFTs')

      if (lnduseon) then
        ! Find the requested years in the file.
        arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The LUC file does not contain first requested year')
        ! Sometimes it is correct to have transient LUC but otherwise have constant conditiions (recycling MET),
        ! in this case metloop is >1 but lnduseon is true. So grab the full length of the LUC file rather than only
        ! the years requested for the met.
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
        else
          arrindex2 = lengthOfFile
        end if
        if (arrindex2 == 0) stop ('getInput says: The LUC file does not contain last requested year')
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.

        allocate(LUCFromFile(icc,ntimes))
        LUCFromFile = ncGet2DVar(lucid,trim(lucVarName),start = [lonloc,latloc,1,arrindex], count = [1,1,icc,ntimes])

        allocate(LUCTime(ntimes))
        LUCTime = int(fileTime(arrindex:arrindex2))
      else
        ! Find the requested year in the file.
        arrindex = fndloc(fileTime,real(fixedYearLUC,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The LUC file does not contain requested year')

        ! We read in only the suggested year
        i = 1 ! offline nlat is always 1 so just set
        m = 1 ! FLAG this is set up only for 1 tile at PRESENT ! JM

        if (nmos /= 1) stop ('getInput for LUC is not setup for more than one tile at present !')

        fcancmxrow(i,m,:) = ncGet1DVar(lucid,trim(lucVarName),start = [lonloc,latloc,1,arrindex], count = [1,1,icc,1])

      end if

    case ('OBSWETF') ! Observed wetland fractions

      lengthOfFile = ncGetDimLen(obswetid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(obswetid, 'time', start = 1, count = lengthOfFile)

      ! The obswetf file is daily (expected format is "day as %Y%m%d.%f")
      ! We want to retain all except any partial day info.
      fileTime=aint(fileTime)

      if (.not. projectedGrid) then
        lonloc = closestCell(obswetid,'lon',longitude)
        latloc = closestCell(obswetid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if

      if (transientOBSWETF) then
        ! Find the beginning and end day in the file.
        ! Assume we are grabbing from first day of start year to last day of last year.

        startTime = real(readMetStartYear,kind=r8) * 10000. + 1. * 100. + 1.
        endTime = real(readMetEndYear,kind=r8) * 10000. + 12. * 100. + 31.

        arrindex = fndloc(fileTime,startTime,1)
        if (arrindex == 0) stop ('getInput says: The OBSWETF file does not contain first requested day')
        
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,endTime,1)
          if (arrindex2 == 0) stop ('getInput says: The OBSWETF file does not contain last requested day')
        else
          arrindex2 = lengthOfFile
        end if
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.

        allocate(OBSWETFFromFile(ntimes))
        OBSWETFFromFile = ncGet1DVar(obswetid,trim(obswetVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])

        allocate(OBSWETFTime(ntimes))
        OBSWETFTime = fileTime(arrindex:arrindex2)
      else
        ! Find the requested day and year in the file.
        ! Assume we are grabbing from day 1
        startTime = real(fixedYearOBSWETF,kind=r8) * 10000. + 1. * 100. + 1.

        ! Find the requested year in the file.
        arrindex = fndloc(fileTime,startTime,1)
        if (arrindex == 0) stop ('getInput says: The OBSWETF file does not contain requested year')

        ! We read in only the suggested year's worth of daily data

        ! If we are using leap years, check if that year is a leap year
        if (leap) then
          call findLeapYears(fixedYearOBSWETF, dummyVar, lastDOY)
        else  
          lastDOY = 365
        end if 

        allocate(OBSWETFFromFile(lastDOY))
        OBSWETFFromFile = ncGet1DVar(obswetid,trim(obswetVarName),start = [lonloc,latloc,arrindex], count = [1,1,lastDOY])

        ! Lastly, remake the OBSWETFTime to be only counting for one year for simplicity
        allocate(OBSWETFTime(lastDOY))
        do d = 1,lastDOY
          OBSWETFTime(d) = d
        end do
      end if

    case ('FER') ! Nitrogen Fertilizer

      lengthOfFile = ncGetDimLen(ferid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(ferid, 'time', start = 1, count = lengthOfFile)

      ! Parse these into just years (expected format is "day as %Y%m%d.%f")
      fileTime = aint(fileTime/10000)
   
      if (.not. projectedGrid) then
        lonloc = closestCell(ferid,'lon',longitude)
        latloc = closestCell(ferid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if

      if (fertilizeron .and. transientFER) then

        ! Find the requested years in the file.
        arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The FER file does not contain first requested year')

        ! Sometimes it is correct to have transient FER but otherwise have constant conditiions (recycling MET),
        ! in this case metloop is >1 but transientFER is true. So grab the full length of the FER file rather than only
        ! the years requested for the met.
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
        else
          arrindex2 = lengthOfFile
        end if

        if (arrindex2 == 0) stop ('getInput says: The FER file does not contain last requested year')
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.
        allocate(FERFromFile(ntimes))
        FERFromFile = ncGet1DVar(ferid,trim(ferVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])

        allocate(FERTime(ntimes))
        FERTime = int(fileTime(arrindex:arrindex2))

      else

        ! Find the requested year in the file.
        arrindex = fndloc(fileTime,real(fixedYearFER,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The FER file does not contain first requested year')
        ! We read in only the suggested year
        i = 1 ! offline nlat is always 1 so just set
        nfertilrow(:,i) = ncGet1DVar(ferid, trim(ferVarName), start = [lonloc,latloc,arrindex], count = [1,1,1])
        
      end if

    case ('DEP') ! Nitrogen Deposition

      lengthOfFile = ncGetDimLen(depid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(depid, 'time', start = 1, count = lengthOfFile)

      ! Parse these into just years (expected format is "day as %Y%m%d.%f")
      fileTime = aint(fileTime/10000)

      if (.not. projectedGrid) then
        lonloc = closestCell(depid,'lon',longitude)
        latloc = closestCell(depid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if

      if (depositionon .and. transientDEP) then

        ! Find the requested years in the file.
        arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The DEP file does not contain first requested year')

        ! Sometimes it is correct to have transient DEP but otherwise have constant conditiions (recycling MET),
        ! in this case metloop is >1 but transientDEP is true. So grab the full length of the DEP file rather than only
        ! the years requested for the met.
        if (metLoop == 1) then
          arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
        else
          arrindex2 = lengthOfFile
        end if

        if (arrindex2 == 0) stop ('getInput says: The DEP file does not contain last requested year')
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.
        allocate(DEPFromFile(ntimes))
        DEPFromFile = ncGet1DVar(depid,trim(depVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])

        allocate(DEPTime(ntimes))
        DEPTime = int(fileTime(arrindex:arrindex2))

      else

        ! Find the requested year in the file.
        arrindex = fndloc(fileTime,real(fixedYearDEP,kind=r8),1)
        if (arrindex == 0) stop ('getInput says: The DEP file does not contain first requested year')
        ! We read in only the suggested year
        i = 1 ! offline nlat is always 1 so just set
        ndepositrow(:,i) = ncGet1DVar(depid, trim(depVarName), start = [lonloc,latloc,arrindex], count = [1,1,1])

      end if

    case ('TIMHAR') ! timber harvest with dynamic tiling

      lengthOfFile = ncGetDimLen(timid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(timid, 'time', start = 1, count = lengthOfFile)
      fileTime=aint(fileTime/10000)

      if (.not. projectedGrid) then
        lonloc = closestCell(timid,'lon',longitude)
        latloc = closestCell(timid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if


      ! Find the requested years in the file.
      arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
      if (arrindex == 0) stop ('getInput says: The timber harvest file does not contain first requested year')
      if (metLoop > 1) stop ('getInput says: timber harvest with tiling does not make sense when recycling MET years')
      arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
      if (arrindex2 == 0) stop ('getInput says: The timber harvest file does not contain last requested year')
      ntimes = arrindex2 - arrindex + 1
      ! Read in and keep only the required elements.
      allocate(TIMFromFile(ntimes))
      TIMFromFile = ncGet1DVar(timid,trim(timVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])
      allocate(TIMTime(ntimes))
      TIMTime = int(fileTime(arrindex:arrindex2))

    case ('PRSFIRE') ! prescribed fire with dynamic tiling

      lengthOfFile = ncGetDimLen(fireid, 'time')
      allocate(fileTime(lengthOfFile))

      fileTime = ncGetTime(fireid, 'time', start = 1, count = lengthOfFile)
      fileTime=aint(fileTime/10000)

      if (.not. projectedGrid) then
        lonloc = closestCell(fireid,'lon',longitude)
        latloc = closestCell(fireid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if


      ! Find the requested years in the file.
      arrindex = fndloc(fileTime,real(readMetStartYear,kind=r8),1)
      if (arrindex == 0) stop ('getInput says: The prescribed fire file does not contain first requested year')
      if (metLoop > 1) stop ('getInput says: prescribed fire with tiling does not make sense when recycling MET years')
      arrindex2 = fndloc(fileTime,real(readMetEndYear,kind=r8),1)
      if (arrindex2 == 0) stop ('getInput says: The timber harvest file does not contain last requested year')
      ntimes = arrindex2 - arrindex + 1
      ! Read in and keep only the required elements.
      allocate(FireFromFile(ntimes))
      FireFromFile = ncGet1DVar(fireid,trim(fireVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])
      allocate(FireTime(ntimes))
      FireTime = int(fileTime(arrindex:arrindex2))

    case ('FourBandAlbedo') ! Four-band albedo parameters

      ! Read in the look up table that is used by the four-band albedo parameterization

      ! Determine the size of the arrays to be read in.
      albdim = NSWE * NREFF * NSMU * NSALB * NBC
      allocate(tmpalb(4,4,albdim))

      ! Diffuse albedo
      tmpalb(1,1,:) = ncGet1DVar(albid, 'albedoDiffuse1', start = [1], count = [albdim])
      tmpalb(1,2,:) = ncGet1DVar(albid, 'albedoDiffuse2', start = [1], count = [albdim])
      tmpalb(1,3,:) = ncGet1DVar(albid, 'albedoDiffuse3', start = [1], count = [albdim])
      tmpalb(1,4,:) = ncGet1DVar(albid, 'albedoDiffuse4', start = [1], count = [albdim])

      ! Direct albedo
      tmpalb(2,1,:) = ncGet1DVar(albid, 'albedoDirect1', start = [1], count = [albdim])
      tmpalb(2,2,:) = ncGet1DVar(albid, 'albedoDirect2', start = [1], count = [albdim])
      tmpalb(2,3,:) = ncGet1DVar(albid, 'albedoDirect3', start = [1], count = [albdim])
      tmpalb(2,4,:) = ncGet1DVar(albid, 'albedoDirect4', start = [1], count = [albdim])

      ! Diffuse transmissivity
      tmpalb(3,1,:) = ncGet1DVar(albid, 'transmisDiffuse1', start = [1], count = [albdim])
      tmpalb(3,2,:) = ncGet1DVar(albid, 'transmisDiffuse2', start = [1], count = [albdim])
      tmpalb(3,3,:) = ncGet1DVar(albid, 'transmisDiffuse3', start = [1], count = [albdim])
      tmpalb(3,4,:) = ncGet1DVar(albid, 'transmisDiffuse4', start = [1], count = [albdim])

      ! Direct transmissivity
      tmpalb(4,1,:) = ncGet1DVar(albid, 'transmisDirect1', start = [1], count = [albdim])
      tmpalb(4,2,:) = ncGet1DVar(albid, 'transmisDirect2', start = [1], count = [albdim])
      tmpalb(4,3,:) = ncGet1DVar(albid, 'transmisDirect3', start = [1], count = [albdim])
      tmpalb(4,4,:) = ncGet1DVar(albid, 'transmisDirect4', start = [1], count = [albdim])

      do i = 1,nbnd_lut
        ipnt = 1
        do isalb = 1,nsalb ! nsfa
          do ismu = 1,nsmu
            do isgs = 1,nreff ! nsgs
              do iswe = 1,nswe
                do ibc = 1,nbc
                  albdif_lut(ibc,iswe,isgs,ismu,isalb,i) = tmpalb(1,i,ipnt)
                  albdir_lut(ibc,iswe,isgs,ismu,isalb,i) = tmpalb(2,i,ipnt)
                  trandif_lut(ibc,iswe,isgs,ismu,isalb,i) = tmpalb(3,i,ipnt)
                  trandir_lut(ibc,iswe,isgs,ismu,isalb,i) = tmpalb(4,i,ipnt)
                  ipnt = ipnt + 1
                end do ! ibc
              end do ! iswe
            end do ! isgs
          end do ! ismu
        end do ! isalb
      end do ! nbnd

      deallocate(tmpalb)

    case ('BlackCarbon')

      ! Black carbon is a daily input.
      lengthOfFile = ncGetDimLen(bcid, 'time')
      allocate(fileTime(lengthOfFile))
      
      fileTime = ncGetTime(bcid, 'time', start = 1, count = lengthOfFile)
      
      ! The BC deposition file is daily (expected format is "day as %Y%m%d.%f")
      ! We want to retain all except the partial day.
      fileTime=aint(fileTime)
      
      if (.not. projectedGrid) then
        lonloc = closestCell(bcid,'lon',longitude)
        latloc = closestCell(bcid,'lat',latitude)
      else
        ! For projected grids, we use the index of the cells, not their coordinates.
        lonloc = projLonInd
        latloc = projLatInd
      end if

      ! Units expected are "kg/m2/s"

      if (blackCtransientDep) then  ! Get a multiyear timeseries of BC dep.    
      
        ! Find the beginning and end day in the file.
        ! Assume we are grabbing from first day of start year to last day of last year.
        ! If metLoop > 1 then just grab the whole file.

        startTime = real(readMetStartYear,kind=r8) * 10000. + 1. * 100. + 1.
        arrindex = fndloc(fileTime,startTime,1)
        if (arrindex == 0) stop ('getInput says: The BC dep file does not contain first requested day')
        
        if (metLoop == 1) then
          endTime = real(readMetEndYear,kind=r8) * 10000. + 12. * 100. + 31.
          arrindex2 = fndloc(fileTime,endTime,1)
          if (arrindex2 == 0) stop ('getInput says: The BC dep file does not contain last requested day')
        else 
          arrindex2 = lengthOfFile
        end if
        ntimes = arrindex2 - arrindex + 1

        ! Read in and keep only the required elements.

        allocate(BCDepFromFile(ntimes))
        BCDepFromFile = ncGet1DVar(bcid,trim(bcDepVarName),start = [lonloc,latloc,arrindex], count = [1,1,ntimes])

        allocate(BCDepTime(ntimes))
        BCDepTime = fileTime(arrindex:arrindex2)

      else ! Just get a single year of values

        ! Find the requested day and year in the file.
        ! Assume we are grabbing from day 1
        startTime = real(fixedYearBCDep,kind=r8) * 10000. + 1. * 100. + 1.

        arrindex = fndloc(fileTime,startTime,1)
        if (arrindex == 0) stop ('getInput says: The BC dep. file does not contain requested year')

        ! We read in only the suggested year of daily inputs

        ! If we are using leap years, check if that year is a leap year
        call findLeapYears(fixedYearBCDep, dummyVar, lastDOY)

        allocate(BCDepFromFile(lastDOY))
        BCDepFromFile = ncGet1DVar(bcid,trim(bcDepVarName),start = [lonloc,latloc,arrindex], count = [1,1,lastDOY])

        ! Lastly, remake the BCDepTime to be only counting for one year for simplicity
        allocate(BCDepTime(lastDOY))
        do d = 1,lastDOY
          BCDepTime(d) = d
        end do

      end if

    case default
      stop ('Specify an input kind for getInput')

    end select

    if (allocated(fileTime)) deallocate(fileTime)

    end associate
  end subroutine getInput

  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_updateInput
  !! @{
  !> Update the input field variable based on the present model timestep
  !> @author Joe Melton, Ed Chan

  subroutine updateInput (inputRequested, yearNeeded, imonth, iday, dom)

    use ctemStateVars, only : vrot, c_switch, vgat, tracer
    use classStateVars, only : class_rot
    use classicParams, only : nmos, testnmos
    use generalUtils,  only : abandonCell
    use tracerModule, only : convertTracerUnits

    implicit none

    character( * ), intent(in) :: inputRequested
    integer, intent(in) :: yearNeeded
    integer, intent(in), optional :: imonth
    integer, intent(in), optional :: iday
    integer, intent(in), optional :: dom            ! day of monthncmxrow
    integer :: arrindex, lengthTime, i, m
    real(r8) :: LGHTTimeNow, OBSWTimeNow, BCdepTimeNow
    character(4) :: seqstring

    associate( &
    co2concrow      => vrot%co2conc,                    & !< real, dimension(:,:) : 
    tracerco2conc   => tracer%tracerCO2rot,             & !< real, dimension(:,:) : 
    ch4concrow      => vrot%ch4conc,                    & !< real, dimension(:,:) : 
    popdinrow       => vrot%popdin,                     & !< real, dimension(:,:) : 
    nfertilrow      => vrot%nfertil,                    & !< real, dimension(:,:) : 
    ndepositrow     => vrot%ndeposit,                   & !< real, dimension(:,:) : 
    timharvrow      => vrot%timharvrow,                 & !< real, timharvrow(:) :
    prsfirerow      => vrot%prsfirerow,                 & !< real, prsfirerow(:) :
    nfcancmxrow     => vrot%nfcancmx,                   & !< real, dimension(:,:,:) : 
    transientLGHT   => c_switch%transientLGHT,          & !< logical: 
    fixedYearLGHT   => c_switch%fixedYearLGHT,          & !< integer: 
    transientOBSWETF=> c_switch%transientOBSWETF,       & !< logical: 
    transientFER    => c_switch%transientFER,           & !< logical: 
    fertilizeron    => c_switch%fertilizeron,           & !< logical: 
    transientDEP    => c_switch%transientDEP,           & !< logical: 
    depositionon    => c_switch%depositionon,           & !< logical: 
    Ncycle_on       => c_switch%Ncycle_on,              & !< logical: 
    blackCtransientDep => c_switch%blackCtransientDep,  & ! logical: If true, allow time-varying black carbon deposition flux in snow processes (4-band scheme)
    lightng         => vgat%lightng,                    & !< real, dimension(:) : total \f$lightning, flashes/(km^2 . year)\f$ it is assumed that cloud
                                                          !!                      to ground lightning is some fixed fraction of total lightning. 
    DEPBROW         => class_rot%DEPBROW,               & !< real, dimension(:) : Black carbon deposition flux \f$[kg m^{-2} s^{-1} ]\f$
    wetfrac_presgat => vgat%wetfrac_pres                & !< real, dimension(:) : 
    )

    select case (trim(inputRequested))

    case ('CO2')

      ! Find the requested year in the file.
      arrindex = fndloc(CO2Time,yearNeeded,1)
      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The CO2 file does not contain requested year: '//seqstring)
      else
        i = 1 ! offline nlat is always 1 so just set
        co2concrow(i,:) = CO2FromFile(arrindex)
      end if

    case ('tracerCO2')

      ! Find the requested year in the file.
      arrindex = fndloc(tracerCO2Time,yearNeeded,1)

      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The tracerCO2 file does not contain requested year: '//seqstring)
      else
        i = 1 ! offline nlat is always 1 so just set
        tracerco2conc(i,:) = tracerCO2FromFile(arrindex)
      end if

      ! Convert the units of the tracer depending on the tracer being simulated.
      tracerco2conc = convertTracerUnits(tracerco2conc)

    case ('CH4')

      ! Find the requested year in the file.
      arrindex = fndloc(CH4Time,yearNeeded,1)
      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The CH4 file does not contain requested year: '//seqstring)
      else
        i = 1 ! offline nlat is always 1 so just set
        ch4concrow(i,:) = CH4FromFile(arrindex)
      end if

    case ('POPD')

      ! Find the requested year in the file.
      arrindex = fndloc(POPDTime,yearNeeded,1)
      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The POPD file does not contain requested year: '//seqstring)
      else
        i = 1 ! offline nlat is always 1 so just set
        popdinrow(i,:) = POPDFromFile(arrindex)
      end if

    case ('LUC')

      ! Find the requested year in the file.
      arrindex = fndloc(LUCTime,yearNeeded,1)
      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The LUC file does not contain requested year: '//seqstring)
      else
        i = 1 ! offline nlat is always 1 so just set
        m = 1 ! FLAG this is set up only for 1 tile at PRESENT ! JM
        if (nmos > 1) stop ('updateInput for LUC only set up for 1 tile at present')
        nfcancmxrow(i,m,:) = LUCFromFile(:,arrindex)
      end if

    case ('TIMHAR')

      ! Find the requested year in the file.
      arrindex = fndloc(TIMTime,yearNeeded,1)

      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The timber harvest file does not contain requested year: '//seqstring)
      else
        i = 1 ! offline nlat is always 1 so just set
        timharvrow(i) = TIMFromFile(arrindex)
      end if

    case ('PRSFIRE')

      ! Find the requested year in the file.
      arrindex = fndloc(FireTime,yearNeeded,1)

      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The prescribed fire file does not contain requested year: '//seqstring)
      else
        i = 1 ! offline nlat is always 1 so just set
        prsfirerow(i) = FireFromFile(arrindex)
      end if

    case ('FER') ! Nitrogen Fertilizer
      if (fertilizeron .and. transientFER) then
        ! Find the requested year in the file.
        arrindex = fndloc(FERTime,yearNeeded,1)
        if (arrindex == 0) then
          write (seqstring,'(I0)') yearNeeded
          call abandonCell('updateInput says: The FER file does not contain requested year: '//seqstring)
        else
          i = 1 ! offline nlat is always 1 so just set
          nfertilrow(:,i) = FERFromFile(arrindex)
        end if
      end if
    case ('DEP') ! Nitrogen Deposition
      if (depositionon .and. transientDEP) then ! Deposition
        ! Find the requested year in the file.
        arrindex = fndloc(DEPTime,yearNeeded,1)
        if (arrindex == 0) then
          write (seqstring,'(I0)') yearNeeded
          call abandonCell('updateInput says: The DEP file does not contain requested year: '//seqstring)
        else
          i = 1 ! offline nlat is always 1 so just set
          ndepositrow(:,i) = DEPFromFile(arrindex)
        end if
      end if

    case ('LGHT')

      ! This file is daily so we need to find the day we are looking for.
      ! imonth is starting at 0 so add 1 always.

      if (transientLGHT) then
        LGHTTimeNow = real(yearNeeded,kind=r8) * 10000. + real(imonth + 1) * 100. + real(dom)
      else ! we only need the day
        LGHTTimeNow = iday
      end if

      ! Find the requested year in the file.
      arrindex = fndloc(LGHTTime,LGHTTimeNow,1)
      if (arrindex == 0) then
        write (seqstring,'(I0)') LGHTTimeNow
        call abandonCell('updateInput says: The LGHT file does not contain requested time: '//seqstring)
      else
        lightng(1) = LGHTFromFile(arrindex)
        ! Since lighning is the same for all tiles, and nlat is always 1 offline, then we
        ! can just pass the same values across all ilg.
        do m = 1,size(lightng)
          lightng(m) = lightng(1)
        end do
      end if

    case ('OBSWETF')

      ! This file is daily so we need to find the day we are looking for.
      ! imonth is starting at 0 so add 1 always.

      if (transientOBSWETF) then
        OBSWTimeNow = real(yearNeeded,kind=r8) * 10000. + real(imonth + 1) * 100. + real(dom)
      else ! we only need the day
        OBSWTimeNow = iday
      end if

      ! Find the requested year in the file.
      arrindex = fndloc(OBSWETFTime,OBSWTimeNow,1)
      if (arrindex == 0) then
        write (seqstring,'(I0)') yearNeeded
        call abandonCell('updateInput says: The OBSWETF file does not contain requested year: '//seqstring)
      else
        wetfrac_presgat(1) = OBSWETFFromFile(arrindex)

        ! Since wetland area is presently assumed the same for all tiles, and nlat is
        ! always 1 offline, then we can just pass the same values across all ilg.
        do m = 1,size(wetfrac_presgat)
          wetfrac_presgat(m) = wetfrac_presgat(1)
        end do
      end if

    case ('BlackCarbon')

      ! This file is daily so we need to find the day we are looking for.
      ! imonth is starting at 0 so add 1 always.

      if (blackCtransientDep) then
        BCdepTimeNow = real(yearNeeded,kind=r8) * 10000. + real(imonth + 1) * 100. + real(dom)
      else ! we only need the day
        BCdepTimeNow = iday
      end if

      ! Find the requested year in the file.
      arrindex = fndloc(BCDepTime,BCdepTimeNow,1)
      if (arrindex == 0) then
        write (seqstring,'(I0)') BCdepTimeNow
        call abandonCell('updateInput says: The BC dep file does not contain requested time: '//seqstring)
      else
        DEPBROW(1) = BCDepFromFile(arrindex)
        ! Since BC dep is the same for all tiles, and nlat is always 1 offline, then we
        ! can just pass the same values across all ilg.
        do m = 1,size(DEPBROW)
          DEPBROW(m) = DEPBROW(1)
        end do
      end if

    case default
      stop ('specify an input kind for updateInput')
    end select

    end associate
  end subroutine updateInput

  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_getMet
  !! @{
  !> Read in the meteorological input from a netcdf file
  !! It is **very** important that the files are chunked correctly (for global and regional runs).
  !! There is an orders of magnitude slow-up otherwise !
  !> @author Joe Melton, Ed Chan

  subroutine getMet (longitude, latitude, nday, projLonInd, projLatInd)

    use fileIOModule
    use classicParams, only : delt
    use ctemStateVars, only : c_switch
    use generalUtils,  only : parseTimeStamp, closeEnough
    use readJobOpts,   only : myDomain

    implicit none

    real, intent(in) :: longitude       !< Longitude of grid cell of interest
    real, intent(in) :: latitude        !< Latitude of grid cell of interest
    integer, intent(in) :: nday         !< Maximum number of physics timesteps in one day

    integer, intent(in), optional :: projLonInd !< Longitude index of the cell for projected grid runs
    integer, intent(in), optional :: projLatInd !< Latitude index of the cell for projected grid runs


    real :: moStart, moEnd, domStart, domEnd !< Assumed start and end months and days of month
    real(r8) :: timeStart, timeEnd           !< Calculated start and end in the format:%Y%m%d.%f
    integer :: lengthOfFile
    integer :: metlon,metlat
    integer :: lonloc, latloc, i
    real(r8), dimension(:), allocatable :: fileTime
    integer :: validTimestep
    integer :: firstIndex, lastIndex
    real, dimension(5) :: firstTime, secondTime

    associate( &
    projectedGrid     => c_switch%projectedGrid,        & !< logical: True if you have a projected lon lat grid, false if not. Projected grids can only have
                                                          !!          regions referenced by the indexes, not coordinates, when running a sub-region
    metFilefracFsf    => c_switch%metFilefracFsf,       & !< character(350): File name for the diffuse fraction of total shortwave radiation input file 
    metFileSnow       => c_switch%metFileSnow,          & !< character(350): location of the preciptiation that is snow meteorology file
    readMetStartYear  => c_switch%readMetStartYear,     & !< integer: First year of meteorological forcing to read in from the met file 
    readMetEndYear    => c_switch%readMetEndYear,       & !< integer: Last year of meteorological forcing to read in from the met file
    IPCP              => c_switch%IPCP                  & !< integer: if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 C.
                                                          !!          if ipcp=2, a linear partitioning of precipitation betweeen
                                                          !!          rainfall and snowfall is done between 0 C and 2 C.
                                                          !!          if ipcp=3, rainfall and snowfall are partitioned according to
                                                          !!          a polynomial curve between 0 C and 6 C. 
                                                          !!          if a snow flux file is supplied IPCP is ignored and read-in values are used.
    )

    !! It is very important that the files have time as the fastest varying dimension.
    !! There is a orders of magnitude slow-up if the dimensions are out of order.

    ! Grab the length of time dimension from the total SW met file and write it to an array.
    ! NOTE: We assume the user is careful enough to ensure the time array is the same
    ! across all met files!
    lengthOfFile = ncGetDimLen(metFssId, 'time')
    allocate(fileTime(lengthOfFile))
    fileTime = ncGetTime(metFssId, 'time', start = 1, count = lengthOfFile)

    ! Construct the time bounds that we will look for in the file.
    ! We assume that you will start on the first timestep of the day.
    ! Further the default is to start on (or at least look for) Jan 1
    ! of the yrStart year.
    moStart = 1.
    domStart = 1.
    ! The first time is considered to be the first physics timestep so given a fractional day of 0.
    timeStart = real(readMetStartYear,kind=r8) * 10000. + moStart * 100. + domStart
    moEnd = 12.
    domEnd = 31.
    ! The last time is considered to be the last physics timestep of the day
    timeEnd =  real(readMetEndYear,kind=r8) * 10000. + moEnd * 100. + domEnd + (real(nday - 1) * delt / 86400.)

    if (timeStart < fileTime(1) .or. timeStart > fileTime(lengthOfFile)) then
      print * ,' *** Check readMetStartYear in your joboptions file '
      print * ,'as it appears to be in conflict with your met file'
      call errorHandler('modelStateDrivers:getMet', - 1)
    endif

    ! Find the array indices closest to the given start/end times
    ! and use these to extract the times from the time variable of the met file.

    firstIndex=minloc(abs(fileTime-timeStart),1)
    if (fileTime(firstIndex) < timeStart) firstIndex=firstIndex+1
    lastIndex=minloc(abs(fileTime-timeEnd),1)
    if (fileTime(lastIndex) > timeEnd) lastIndex=lastIndex-1

    validTimestep=lastIndex-firstIndex+1
    allocate(metTime(validTimeStep))
    metTime=fileTime(firstIndex:lastIndex)

    ! Check that the first day is Jan 1, otherwise warn the user
    firstTime =  parseTimeStamp(metTime(1))
    if (.not. closeEnough(firstTime(5),1.,0.001)) then
    !if (mod(int(metTime(1)),100) /= 1) then
      print * ,'Warning,your met file does not start on Jan 1.'
    end if
    
    ! Check that the first day starts at 0 min, otherwise warn the user
    ! and stop the run.
    if (firstTime(4) /= 0.) then
    !if (metTime(1)-aint(metTime(1)) > 0) then
      print * ,'*** Error,your met file does not start at 0 hour and 0 minutes.'
      call errorHandler('modelStateDrivers:getMet', - 2)
    end if

    ! Determine the time step of the met data and
    ! convert from fraction of day to period in seconds
    secondTime =  parseTimeStamp(metTime(2))
    
    metInputTimeStep = (secondTime(4) - firstTime(4)) * 86400.
    
    ! It is extremely unlikely that met data would be in timesteps that are 
    ! fractions of a second. So make sure to round to the nearest whole value 
    metInputTimeStep = float(nint(metInputTimeStep))

    ! Find the closest cell to our lon and lat
    if (.not. projectedGrid) then
      lonloc = closestCell(metFssId,'lon',longitude)
      latloc = closestCell(metFssId,'lat',latitude)
    else
      ! For projected grids, we use the index of the cells, not their coordinates.
      ! So the index has been passed in as a real, convert here to an integer.
      lonloc = projLonInd
      latloc = projLatInd 
      ! If the input met files are subsets of the domain, then the indices need to 
      ! adjusted to match the correct columns/rows in the subsets. 
      ! This is only the case if the lons or lats in the met files are less than in the init file.
      metlon = ncGetDimLen(metFssId,'lon')
      metlat = ncGetDimLen(metFssId,'lat')
      if ( metlon < totlon ) lonloc = lonloc - myDomain%srtx + 1
      if ( metlat < totlat ) latloc = latloc - myDomain%srty + 1
    end if

    ! Now read in the whole MET times series and store it for each variable
    allocate(metFss(validTimestep),metFdl(validTimestep),metPre(validTimestep), &
    metTa(validTimestep),metQa(validTimestep),metUv(validTimestep),metPres(validTimestep))

    ! Only allocate the metFsf if we actually have that as an input file.
    if (trim(metFilefracFsf) /= '' .and. .not. allocated(metfracFsf)) allocate(metfracFsf(validTimestep))

    ! Only allocate the metFsf if we actually have that as an input file.
    if (trim(metFileSnow) /= '') then
       if (.not. allocated(metSnow)) allocate(metSnow(validTimestep))
    end if 

    ! NOTE: Carefully check that your incoming inputs are in the expected units !

    ! WARNING. If you use ncdump on a file it will show the opposite order for the
    ! dimensions of a variable than how fortran reads them in. So var(lat,lon,time) is actually
    ! var(time,lon,lat) from the perspective of fortran. Pay careful attention!

    metFss = ncGet1DVar(metFssId,trim(metFssVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
    if (allocated(metfracFsf)) then
      metfracFsf = ncGet1DVar(metfracFsfId,trim(metfracFsfVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
    end if 
    metFdl = ncGet1DVar(metFdlId,trim(metFdlVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
    metPre = ncGet1DVar(metPreId,trim(metPreVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
    if (allocated(metSnow)) then
      metSnow = ncGet1DVar(metSnowId,trim(metSnowVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
      IPCP = 4 ! Set it to overwrite whatever is in the job options file. IPCP = 4 tells it to use the read in values.
    end if 
    metTa = ncGet1DVar(metTaId,trim(metTaVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
    metQa = ncGet1DVar(metQaId,trim(metQaVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
    metUv = ncGet1DVar(metUvId,trim(metUvVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])
    metPres = ncGet1DVar(metPresId,trim(metPresVarName),start = [lonloc,latloc,firstIndex], count = [1,1,validTimestep])

    end associate
  end subroutine getMet

  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_updateMet
  !! @{
  !> This transfers the met data of this time step from the read-in array to the
  !! instantaneous variables. This also sets iyear to the present year of MET being read in.
  !> @author Joe Melton

  subroutine updateMet (metTimeIndex, lastDOY, iyear, iday, ihour, imin, metDone)

    use classicParams,  only : delt,tfrez
    use classStateVars, only : class_rot
    use generalUtils,   only : parseTimeStamp,abandonCell
    use outputManager,  only : stepct

    implicit none

    integer, intent(in)  :: metTimeIndex        !< Index to read from met file
    integer, intent(in) :: lastDOY             !< Initialized to 365 days, can be overwritten later is leap = true and it is a leap year.
    integer, intent(inout) :: iyear               !< Present year of simulation
    integer, intent(out) :: iday                !< Present day of simulation
    integer, intent(out) :: ihour               !< Present hour of simulation
    integer, intent(out) :: imin                !< Present minute of simulation
    logical, intent(out) :: metDone             !< Switch signalling end of met data
    
    integer :: i, numsteps, stepPerHour
    real, dimension(5) :: theTime
    real :: dayfrac, month, dom, minute
    integer :: laststepyr !< iyear from the last timestep

    associate( &
    FSSROW => class_rot%FSSROW,      & !< real, dimension(:) : Total shortwave radiation \f$[W m^{-2} ]\f$
    fracFSFROW => class_rot%fracFSFROW, & !< real, dimension(:) : Diffuse fraction of total shortwave radiation (if provided) \f$[ ]\f$ 
    FDLROW => class_rot%FDLROW,      & !< real, dimension(:) : Downwelling longwave sky radiation \f$[W m^{-2} ]\f$ 
    FDLROT => class_rot%FDLROT,      & !< real, dimension(:,:) : Downwelling longwave sky radiation (tile version of FDLROW) \f$[W m^{-2} ]\f$ 
    PREROW => class_rot%PREROW,      & !< real, dimension(:) : Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$ 
    RPREROW => class_rot%RPREROW,    & !< real, dimension(:) : Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    SPREROW => class_rot%SPREROW,    & !< real, dimension(:) : Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    TAROW => class_rot%TAROW,        & !< real, dimension(:) : Air temperature at reference height [K] 
    QAROW => class_rot%QAROW,        & !< real, dimension(:) : Specific humidity at reference height \f$[kg kg^{-1}]\f$ 
    UVROW => class_rot%UVROW,        & !< real, dimension(:) : Wind speed at reference height \f$[m s^{-1} ]\f$ 
    PRESROW => class_rot%PRESROW,     & !< real, dimension(:) : Surface air pressure \f$[P_a]\f$ 
    FSVHROW => class_rot%FSVHROW,    & !< real, dimension(:,:) : Visible shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    FSIHROW => class_rot%FSIHROW,    & !< real, dimension(:,:) : Near infrared shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    ULROW => class_rot%ULROW,        & !< real, dimension(:) : Zonal component of wind velocity \f$[m s^{-1} ]\f$
    VLROW => class_rot%VLROW,        & !< real, dimension(:) : Meridional component of wind velocity \f$[m s^{-1} ]\f$
    VMODROW => class_rot%VMODROW,    & !< real, dimension(:) : Wind speed at reference height \f$[m s^{-1} ]\f$
    FSSBROL => class_rot%FSSBROL,     & !< real, dimension(:,:,:) : Total solar radiation in each modelled wavelength band  [W m-2]
    FSDBROL => class_rot%FSDBROL,     & !< real, dimension(:,:,:) : Direct solar radiation in each modelled wavelength band  [W m-2]
    FSFBROL => class_rot%FSFBROL     & !< real, dimension(:,:,:) : Diffuse solar radiation in each modelled wavelength band  [W m-2]
    )

    metDone = .false.
    laststepyr = iyear
    
    ! Find the timestep info from the array already read in.
    theTime =  parseTimeStamp(metTime(metTimeIndex))

    iyear = int(theTime(1))

    ! Check if we have begun a new year and if the year that just finished had the right amount of timesteps.
    if (iyear /= laststepyr) then
      if (laststepyr == 0) then 
        stepct = 0
      else if (nint(24. / (delt / 3600.)) * lastDOY /= stepct) then
        print*,'Year ',laststepyr,'has only ',stepct,'timesteps instead of expected ', nint(24. / (delt / 3600.)) * lastDOY
        call abandonCell('updateMet says: Check your met file')
      else
        stepct = 0
      end if 
    end if 
    
    month = theTime(2)
    dom = theTime(3)
    dayfrac = theTime(4)
    iday = int(theTime(5))
    
    !> The dayfrac can then be parsed to give the hour and minute.
    numsteps = nint(dayfrac * 24. / (delt / 3600.))
    stepPerHour = nint(3600./delt)
    ihour = floor(real(numsteps) / stepPerHour)
    minute = mod(numsteps,stepPerHour)
    imin = nint(minute) * int((delt / 60.))

    !> The meteorological data is then passed to the instantaneous variables
    !! from the larger variables that store the run's met data read in earlier.
    i = 1 ! always 1 offline
    FSSROW(I)   = metFss(metTimeIndex)
    if (allocated(metfracFsf)) then 
      fracFSFROW(I)  = metfracFsf(metTimeIndex) 
    else
      fracFSFROW(I) = -9999. ! Flag used to denote no real value. 
    end if 
    FDLROW(i)   = metFdl(metTimeIndex)

    if (.not. allocated(metSnow)) then ! No snow given so just operate on the total P.
      PREROW(i)   = metPre(metTimeIndex)
    else ! We have the snow so find the rain from the total.
      PREROW(i)   = metPre(metTimeIndex)
      SPREROW(i) = metSnow(metTimeIndex)
      RPREROW(i) = max(0., min(metPre(metTimeIndex), metPre(metTimeIndex) - SPREROW(i)))
    end if 

    TAROW(i)    = metTa(metTimeIndex) ! This is converted from the read-in degree C to K in main_driver !
    QAROW(i)    = metQa(metTimeIndex)
    
    ! To prevent a divide by zero in atmosphericVarsCalc, we set this lower limit on the specific humidity.
    if (QAROW(i) == 0.) then
      QAROW(i) = 1.E-6
      ! print * ,'Warning, specific humidity of 0 in your input file. metTimeindex = ',metTimeIndex
      ! print * ,'setting to 1.E-6 g/kg and moving on (updateMet)'
    end if

    UVROW(i)    = metUv(metTimeIndex)
    PRESROW(i)  = metPres(metTimeIndex)

    ! Sanity check the values read in. This is in C so > 100 or < -100 is unreasonable.
    if (TAROW(i) > 100. .or. TAROW(i) < -100.) then
      print*,'TAROW as read in by updateMet is unreasonable: ',TAROW(i)
      call abandonCell('updateMet says: Check your met file')
    end if 

    ! Sanity check the values read in. This is mm/s so > 100 mm/s or a negative value is unreasonable.
    ! World record is apparently 38mm in 1 min.
    if (PREROW(i) > 100. .or. PREROW(i) < 0.) then
      print*,'PREROW as read in by updateMet is unreasonable: ',PREROW(i)
      call abandonCell('updateMet says: Check your met file')
    end if 

    !> If the end of the timeseries is reached, change the metDone switch to true.
    if (metTimeIndex ==  size(metTime)) metDone = .true.

    !> Generally only the total incoming shortwave radiation FSDOWN
    !! is available; so it is partitioned 50:50 between the incoming visible (FSVHROW)
    !! and near-infrared (FSIHROW) radiation.  The first two elements of the
    !! generalized incoming radiation array, FSSBROL (used for both the ISNOALB=0
    !! and ISNOALB=1 options) are set to FSVHROW and FSIHROW respectively.
    !! The air temperature TAROW is converted from degrees C to K.  The zonal
    !! (ULROW) and meridional (VLROW) components of the wind speed are generally not
    !! used; only the overall wind speed UVROW is
    !! measured.  However, CLASSIC does not require wind direction for its calculations,
    !! so ULROW is arbitrarily assigned the value of UVROW and VLROW is set to zero for
    !! this run.  The input wind speed VMODROW is assigned the value of UVROW.

    i = 1 ! always 1 offline

    FDLROT(I,:) = FDLROW(I)
    FSVHROW(I,:) = 0.5 * FSSROW(I)
    FSIHROW(I,:) = 0.5 * FSSROW(I)
    TAROW(I) = TAROW(I) + TFREZ
    ULROW(I) = UVROW(I)
    VLROW(I) = 0.0
    VMODROW(I) = UVROW(I)
    FSSBROL(I,:,1) = FSVHROW(I,:)
    FSSBROL(I,:,2) = FSIHROW(I,:)

   ! if ISNOALB=1, 4-band albedo requires direct and diffuse components for each band. Both are partitioned 50:50 between the
   ! incoming visible and near-infrared, and the near-infrared is splitted into the 3-band using 0.75, 0.23, and 0.02 provided by
   ! Jason (based on radiative transfer simulations). Use residual for the last band to avoid possible small surplus due to
   ! numerical rounding.
    FSFBROL(I,:,1) = 0.5 * fracFSFROW(I)*FSSROW(I)
    FSFBROL(I,:,2) = FSFBROL(I,:,1)*0.75
    FSFBROL(I,:,3) = FSFBROL(I,:,1)*0.23
    FSFBROL(I,:,4) = FSFBROL(I,:,1)-FSFBROL(I,:,2)-FSFBROL(I,:,3)
    FSDBROL(I,:,1) = 0.5 * (FSSROW(I)-fracFSFROW(I)*FSSROW(I))
    FSDBROL(I,:,2) = FSDBROL(I,:,1)*0.75
    FSDBROL(I,:,3) = FSDBROL(I,:,1)*0.23
    FSDBROL(I,:,4) = FSDBROL(I,:,1)-FSDBROL(I,:,2)-FSDBROL(I,:,3)

    !if (allocated(metfracFsf)) then
      ! If we have the diffuse radiation read in from an external file, we could assign it to the 
      ! FSFBROL variable but I have not done that here as it is only used in the 4-band albedo and
      ! I am not sure what would make the most sense for its various bands... JM July 2021
    !end if
    stepct = stepct + 1

    end associate
    return

  end subroutine updateMet

  !! @}
  !---------------------------------------------------------------------------------------
  
  !> \ingroup modelstatedrivers_checkTimeUnits
  !! @{
  !> Checks the netcdf file attributes to make sure the time units are in the 
  !! expected "day as %Y%m%d.%f" format.
  !! @author Joe Melton
  !!
  subroutine checkTimeUnits(ncid,fileName)

    implicit none

    integer, intent(in)  :: ncid
    character(*), intent(in) :: fileName
    character(80)     :: fileUnits
    character(len=80) :: expected = "day as %Y%m%d.%f"
    character(len=80) :: altexpected = "day as YYYYMMDD.FFFF" ! possible variant.

    ! Get the time units that are in the file
    fileUnits = ncGetAtt(ncid,'time','units')
    ! Check that against what we are expecting
    if (trim(fileUnits) /= trim(expected) .and. trim(fileUnits) /= trim(altexpected)) then
      print*,'Time units in your file',trim(fileName),' are "',trim(fileUnits),'"'
      print*,'but CLASSIC expects: ',expected
      print*,'----- Run aborting. -----'
      stop
    end if     
  end subroutine checkTimeUnits
  !! @}

  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_closestCell
  !! @{
  !> Finds the closest grid cell in the file
  !> @author Joe Melton

  integer function closestCell (ncid, label, gridPoint)

    use fileIOModule

    implicit none

    integer, intent(in) :: ncid
    character( * ), intent(in) :: label
    real, intent(in) :: gridPoint
    integer :: lengthdim
    real, dimension(:), allocatable :: filevals
    integer, dimension(1) :: tempintarr

    lengthdim = ncGetDimLen(ncid,label)
    allocate(filevals(lengthdim))
    filevals = ncGet1DVar(ncid,label,start = [1], count = [lengthdim])
    filevals = filevals - gridPoint
    tempintarr = minloc(abs(filevals))
    closestCell = tempintarr(1)

  end function closestCell
  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup modelstatedrivers_deallocInput
  !! @{
  !> Deallocates the input files arrays
  !> @author Joe Melton

  subroutine deallocInput

    implicit none

    if (allocated(CO2Time))       deallocate(CO2Time)
    if (allocated(CO2FromFile))   deallocate(CO2FromFile)
    if (allocated(CH4Time))       deallocate(CH4Time)
    if (allocated(CH4FromFile))   deallocate(CH4FromFile)
    if (allocated(POPDTime))      deallocate(POPDTime)
    if (allocated(POPDFromFile))  deallocate(POPDFromFile)
    if (allocated(LGHTTime))      deallocate(LGHTTime)
    if (allocated(LGHTFromFile))  deallocate(LGHTFromFile)
    if (allocated(LUCTime))       deallocate(LUCTime)
    if (allocated(LUCFromFile))   deallocate(LUCFromFile)
    if (allocated(FERTime))       deallocate(FERTime)
    if (allocated(FERFromFile))   deallocate(FERFromFile)
    if (allocated(DEPTime))       deallocate(DEPTime)
    if (allocated(DEPFromFile))   deallocate(DEPFromFile)
    if (allocated(OBSWETFTime)) deallocate(OBSWETFTime)
    if (allocated(OBSWETFFromFile)) deallocate(OBSWETFFromFile)
    if (allocated(BCDepTime)) deallocate(BCDepTime)
    if (allocated(BCDepFromFile)) deallocate(BCDepFromFile)
    if (allocated(TIMFromFile)) deallocate(TIMFromFile)
    if (allocated(TIMTime)) deallocate(TIMTime)
    if (allocated(FireFromFile)) deallocate(FireFromFile)
    if (allocated(FireTime)) deallocate(FireTime)

    deallocate(metTime,metFss,metFdl,metPre,metPres,metQa,metTa,metUv)

  end subroutine deallocInput

  !! @}
  !> \namespace modelstatedrivers
  !> Central driver to read in, and write out all model state variables (replacing INI and CTM files)
  !! as well as the model inputs such as MET, population density, land use change, CO2 etc.

end module modelStateDrivers
