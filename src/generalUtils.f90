!> \file
!> Central module for all general utilities
!!
module generalUtils

  use, intrinsic :: iso_fortran_env, only: r8=>real64
  implicit none

  public :: fndloc
  public :: abandonCell
  public :: findDaylength
  public :: findWaterTable
  public :: calcEsat
  public :: findCloudiness
  public :: makeDOYArray
  public :: findLeapYears
  public :: findPermafrostVars
  public :: parseTimeStamp
  public :: closeEnough
  public :: initRandomSeed
  public :: checksumCalc

  logical :: run_model           !< Simple logical switch to either keep run going or finish

contains

  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_fndloc
  !! @{
  !> Mimics the functionality of the intrinsic "findloc" which is included in the 
  !! Fortran 2008 standard, but missing from gfortran below version 9. This is intended 
  !! only for use with 64-bit real and integer inputs. And only for searches of 1D arrays. 
  !! This function should be replaced by "findloc", if it is available (by substituting 
  !! all instances of "fndloc" with "findloc"). The argument list should be the same for both.
  !! @author Ed Chan
  !!
  function fndloc(array,value,axis)

    implicit none

    integer :: fndloc
    class(*), intent(in) :: array(:)
    class(*), intent(in) :: value
    integer, optional, intent(in) :: axis ! not used, only for compatibility with findloc
    integer :: i

    fndloc=-1

    select type (value) 
      type is (integer)
        select type (array)
          type is (integer)   ! both array/value are integer
            fndloc = 0
            do i = 1, size(array)
              if (array(i) == value) then
                fndloc = i
                exit
              end if
            end do
        end select
      type is (real(r8))
        select type (array)
          type is (real(r8)) ! both array/value are r8
            fndloc = 0
            do i = 1, size(array)
              if (array(i) == value) then
                fndloc = i
                exit
              end if
            end do
        end select
    end select

    if ( fndloc == -1 ) stop 'Both input variables must be either real64 or integer'

  end function fndloc

  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_abandonCell
  !! @{
  !> Used to stop running a model grid cell. For errors that need to be caught early in a run,
  !! the fortran intrinsic 'stop' is preferred but for errors later in a run or simple fails on
  !! single grid cells, abandonCell is best since it allows the netcdf files to continue to
  !! written to and won't disrupt the MPI processes (as stop does)
  !! @author Joe Melton and Ed Wisernig
  !!
  subroutine abandonCell (errmsg)

    use classStateVars, only : class_rot
    use ctemStateVars, only : c_switch

    implicit none

    character( * ), intent(in), optional :: errmsg

    associate( &
    DLATROW         => class_rot%DLATROW,               & !< real, dimension(:) : 
    DLONROW         => class_rot%DLONROW,               & !< real, dimension(:) : 
    latIndexROW     => class_rot%latIndexROW,           & !< real, dimension(:) : Index of grid cell being run on the input files grid (latitude)     
    lonIndexROW     => class_rot%lonIndexROW,           & !< real, dimension(:) : Index of grid cell being run on the input files grid (longitude) 
    projectedGrid   => c_switch%projectedGrid           & !< logical: 
    )

    run_model = .false.
    if (present(errmsg)) then
      print*,errmsg
      if (.not. projectedGrid) then 
        print*,'exiting cell: ',DLONROW,DLATROW
      else 
        print*,'exiting index: ',lonIndexROW,latIndexROW
      end if 
    end if
    if (.not. projectedGrid) then 
      print*,'died on',DLONROW,DLATROW
    else 
      print*,'died on index',lonIndexROW,latIndexROW
    end if 
    
    end associate
    return

  end subroutine abandonCell
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findDaylength
  !! @{
  !> Calculate the daylength based on the latitude and day of year
  !! @author Joe Melton
  !!
  real function findDaylength (solday, radl)

    ! Joe Melton Dec 18 2015 (taken from phenlogy.f)

    use classicParams, only : pi

    implicit none

    real, intent(in) :: solday  ! day of year
    real, intent(in) :: radl    ! latitude
    real :: theta               ! temp var
    real :: decli               ! temp var
    real :: term                ! temp var

    theta = 0.2163108 + 2.0 * atan(0.9671396 * tan(0.0086 * (solday - 186.0)))
    decli = asin(0.39795 * cos(theta))      ! declination ! note I see that CLASS does this also but with different formula...
    term = (sin(radl) * sin(decli))  /(cos(radl) * cos(decli))
    term = max( - 1.0,min(term,1.0))
    findDaylength = 24.0 - (24.0/pi) * acos(term)

  end function findDaylength
  !! @}
  
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findWaterTable
  !! @{
  !> Calculate the water table depth based on soil moisture
  !! @author Joe Melton
  !!
  real function findWaterTable(ig,thliq,thice,thpor,zbotw,delzw,thlret)

    ! Joe Melton Apr 2 2020 (adapted from energyBudgetPrep.f90)
    
    implicit none

    integer, intent(in) :: ig  !< number of ground layers
    real, intent(in) :: thliq(:) !< Liquid water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
    real, intent(in) :: thice(:) !< Frozen water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
    real, intent(in) :: thpor(:)    !< Pore volume in soil layer \f$[m^3 m^{-3}] (\theta_p)\f$
    real, intent(in) :: thlret(:)   !< Liquid water retention capacity for organic soil \f$[m^3 m^{-3}] (\theta_{ret})\f$
    real, intent(in) :: delzw(:)    !< Permeable thickness of soil layer \f$[m] (\Delta z_w)\f$
    real, intent(in) :: zbotw(:)    !< Depth to permeable bottom of soil layer \f$[m] (z_{b,w})\f$
    
    integer :: j,indexLastPermLayer
    
    ! Find the index of the last layer that is permeable (bedrock below)
    do j = 2,ig 
      indexLastPermLayer = j - 1
      if (delzw(j) < 0.0001) exit           
    end do
    
    ! Then we find the water table. We start at the surface and move 
    ! down the soil column. The first saturated layer represents the start of 
    ! the water table. We then look at the layer above and consider its level 
    ! of saturation to put the water table partially in that layer depending 
    ! on how close to saturation it is ()   
    
    ! Set the water table to the bottom of the last permeable layer so if the whole 
    ! soil column is unsaturated then the water table has a value. 
    findWaterTable = zbotw(indexLastPermLayer)
    do j = 1, indexLastPermLayer
      if ( (thliq(j) + thice(j)) > (thpor(j) - 0.01) ) then 
        findWaterTable = zbotw(j) - delzw(j)
        if (j /= 1) then 
          findWaterTable = findWaterTable -  delzw(j-1) * &
                        min(1.0, (thliq(j-1) + thice(j-1) - thlret(j-1)) / &
                                  (thpor(j-1) - thlret(j-1)))
        end if
        exit
      end if 
    end do 

  end function findWaterTable
  !! @}
  
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_calcEsat
  !! @{
  !> Calculate the saturated vapour pressure in Pa. Based upon 
  !! the parameterization of Emanuel, 1994 \cite Emanuel1994-dt. 
  !! @author Joe Melton
  !!
  real function calcEsat(ta)

    use classicParams, only : TFREZ

    implicit none

    real, intent(in) :: ta  ! air/canopy temperature (K)
    
    if (ta >= tfrez) then
      calcEsat = exp(53.67957 - 6743.769/ta - 4.8451 * log(TA)) * 100. !100 converts from hPa to Pa.
    else !
      calcEsat = exp(23.33086 - 6111.72784/ta + 0.15215 * log(TA)) * 100. !100 converts from hPa to Pa.
    end if
      
  end function calcEsat
  !! @}

  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_makeDOYArray
  !! @{
  !> Create an array of the number of days in each year for the model run.
  !! @author Joe Melton
  !!
  subroutine makeDOYArray(readMetStartYear, readMetEndYear, arrayDOYs)

    implicit none

    integer, intent(in)    :: readMetStartYear  !< First year of met read in
    integer, intent(in)    :: readMetEndYear    !< Last year of met read in
    integer, intent(inout) :: arrayDOYs(:)      !< Array to fill with last days of year

    integer :: i, k
    logical :: leapnow
    integer :: lastDOY

    k = 1
    do i = readMetStartYear,readMetEndYear
      call findLeapYears(i,leapnow,lastDOY)
      arrayDOYs(k) = lastDOY
      k = k + 1
    end do

  end subroutine makeDOYArray

    !! @}

  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findLeapYears
  !! @{
  !> Check if this year is a leap year
  !! @author Joe Melton
  !!
  subroutine findLeapYears (iyear, leapnow, lastDOY)

    use classicParams, only : monthend, mmday, monthdays

    implicit none

    logical, intent(out)   :: leapnow
    integer, intent(in)    :: iyear
    integer, intent(out) :: lastDOY

    if (mod(iyear,4) /= 0) then ! it is a common year
      leapnow = .false.
    else if (mod(iyear,100) /= 0) then ! it is a leap year
      leapnow = .true.
    else if (mod(iyear,400) /= 0) then ! it is a common year
      leapnow = .false.
    else ! it is a leap year
      leapnow = .true.
    end if

    if (leapnow) then
      lastDOY = 366
    else
      lastDOY = 365
    end if

    ! We do not check the MET files to make sure the incoming MET is in fact
    ! 366 days if leapnow. You must verify this in your own input files. Later
    ! in the code it will fail and print an error message to screen warning you
    ! that your file is not correct.
    if (leapnow) then ! adjust the calendar and set the error check.
      monthdays = (/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
      monthend = (/ 0,31,60,91,121,152,182,213,244,274,305,335,366 /)
      mmday = (/ 16,46,76,107,137,168,198,229,260,290,321,351 /)

    else
      monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] !< days in each month
      monthend  = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] !< calender day at end of each month
      mmday     = [ 16,46,75,106,136,167,197,228,259,289,320,350 ] !< mid-month day
    end if

  end subroutine findLeapYears
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findCloudiness
  !! @{
  !> The cosine of the solar zenith angle COSZ is calculated from the day of
  !> the year, the hour, the minute and the latitude using basic radiation geometry,
  !> and (avoiding vanishingly small numbers) is assigned to CSZROW.  The fractional
  !> cloud cover FCLOROW is commonly not available so a rough estimate is
  !> obtained by setting it to 1 when precipitation is occurring, and to the fraction
  !> of incoming diffuse radiation XDIFFUS otherwise (assumed to be 1 when the sun
  !> is at the horizon, and 0.10 when it is at the zenith). When incoming diffuse SW radiation is 
  !> available it is used, along with incoming total SW radiation to calculate XDIFFUS.
  !> the fractional cloudiness then follows from XDIFFUS as usual. 
  !! @author Diana Verseghy, Joe Melton
  !!
  subroutine findCloudiness (nltest, imin, ihour, iday, lastDOY)

    use classicParams, only : pi
    use classStateVars, only : class_rot

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: imin
    integer, intent(in) :: ihour
    integer, intent(in) :: iday
    integer, intent(in) :: lastDOY

    integer :: i
    real :: day
    real :: decl    !< Declination
    real :: hour
    real :: cosz    !< Cosine of the zenith angle

    associate( &
    RADJROW => class_rot%RADJROW,    & !< real, dimension(:) : Latitude of grid cell (positive north of equator) [rad] 
    CSZROW => class_rot%CSZROW,      & !< real, dimension(:) : Cosine of solar zenith angle [ ] 
    PREROW => class_rot%PREROW,      & !< real, dimension(:) : Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$ 
    XDIFFUS => class_rot%XDIFFUS,    & !< real, dimension(:,:) : Fraction of diffused radiation 
    FCLOROW => class_rot%FCLOROW,    & !< real, dimension(:) : Fractional cloud cover [ ] 
    FSSROW => class_rot%FSSROW,      & !< real, dimension(:) : Total shortwave radiation \f$[W m^{-2} ]\f$
    fracFSFROW => class_rot%fracFSFROW   & !< real, dimension(:) : Diffuse fraction of total shortwave radiation (if provided) \f$[ - ]\f$ 
    )

    day = real(iday) + (real(ihour) + real(imin)/60.)/24.
    decl = sin(2. * pi * (284. + day)/real(lastDOY)) * 23.45 * pi/180.
    hour = (real(ihour) + real(imin)/60.) * pi/12. - pi

    do i = 1,nltest

      cosz = sin(radjrow(i)) * sin(decl) + cos(radjrow(i)) * cos(decl) * cos(hour)

      cszrow(i) = sign(max(abs(cosz),1.0e-3),cosz)

      if (fracfsfrow(i) < -9998.) then
        ! There is no read in diffuse radiation field so we will construct the xdiffus using 
        ! precipitation and zenith angle
        if (prerow(i) > 0.) then
          xdiffus(i,:) = 1.0
        else
          xdiffus(i,:) = max(0.0,min(1.0 - 0.9 * cosz,1.0))
        end if
      else 
        ! there is read-in fraction of incoming diffuse shortwave radiation so use it.
        if (fssrow(i) > 0.) then
          xdiffus(i,:) = fracfsfrow(i)
        else
          xdiffus(i,:) = 0.
        end if 
      end if 

      fclorow(i) = xdiffus(i,1) ! all tiles are the same for offline
    end do

    end associate
  end subroutine findCloudiness
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_parseTimeStamp
  !! @{
  !> Parses a time stamp in the expected form "day as %Y%m%d.%f"
  !! Returns an array with 1) year, 2) month, 3) day, 4) fraction of day
  !! 5) day of year
  !! @author Joe Melton, Ed Wisernig
  !!
  function parseTimeStamp (timeStamp)

    use classicParams, only : monthdays
    use, intrinsic :: iso_fortran_env, only: r8=>real64

    implicit none

    real, dimension(5)   :: parseTimeStamp
    real(r8), intent(in) :: timeStamp
    real(r8)             :: date
    integer              :: intdate, day, month, year, totdays, t

    date = floor(timeStamp) ! remove the part days
    parseTimeStamp(4) = timeStamp - date ! save the part days
    intdate = int(date)
    day = mod(intdate,100);   intdate = intdate / 100
    month = mod(intdate,100); intdate = intdate / 100
    year = intdate
    parseTimeStamp(1) = real(year)
    parseTimeStamp(2) = real(month)
    parseTimeStamp(3) = real(day)
    totdays = 0
    if (month > 1) then
      do t = 1,month - 1
        totdays = totdays + monthdays(t)
      end do
    end if
    parseTimeStamp(5) = real(totdays + day)
    
  end function parseTimeStamp
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findPermafrostVars
  !! @{
  !> Finds the active layer depth and depth to the frozen water table.
  !! @author Joe Melton
  !!
  subroutine findPermafrostVars (nmtest, nltest, iday)

    use classicParams, only : ignd, tfrez, eftime, efoldfact
    use classStateVars, only : class_rot, class_gat

    implicit none

    integer, intent(in) :: nmtest
    integer, intent(in) :: nltest
    integer, intent(in) :: iday

    integer :: i, j, m

    associate( &
    ftable  => class_rot%ftable,                        & !< real, dimension(:,:)  : Depth to frozen water table (m) 
    actlyr  => class_rot%actlyr,                        & !< real, dimension(:,:)  : Active layer depth (m) 
    maxAnnualActLyr => class_rot%maxAnnualActLyrROT,    & !< real, dimension(:,:) : Active layer depth maximum over the e-folding period specified by parameter eftime (m). 
    actLyrThisYr => class_rot%actLyrThisYrROT,          & !< real, dimension(:,:) : Annual active layer depth maximum starting from summer solstice for the present year (m) 
    tbarrot => class_rot%tbarrot,                       & !< real(r8), dimension(:,:,:) : Temperature of soil layers [K] 
    thlqrot => class_rot%thlqrot,                       & !< real, dimension(:,:,:) : Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    thicrot => class_rot%thicrot,                       & !< real, dimension(:,:,:) : Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    isndrot => class_rot%isndrot,                       & !< integer, dimension(:,:,:) : Sand content flag, used to delineate non-soils. 
    dlzwrot => class_rot%dlzwrot,                       & !< real, dimension(:,:,:) : Permeable thickness of soil layer [m] 
    delz    => class_gat%delz,                          & !< real, dimension(:) : Overall thickness of soil layer [m] 
    thmrot => class_rot%thmrot,                         & !< real, dimension(:,:,:) : Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$ 
    dlatrow => class_rot%dlatrow,                       & !< real, dimension(:) : Latitude (degrees) 
    FAREROT => class_rot%FAREROT                        & !< real, dimension(:,:)   :
    )
    !---

    actlyr = 0.0
    ftable = 0.0
    do j = 1,ignd
      do i = 1,nltest
        do m = 1,nmtest
          if(FAREROT(i,m) > 0.0) then
            if (abs(tbarrot(i,m,j) - tfrez) < 0.0001) then
              if (isndrot(i,m,j) > - 3) then
                actlyr(i,m) = actlyr(i,m) + (thlqrot(i,m,j) / (thlqrot(i,m,j) &
                              + thicrot(i,m,j))) * dlzwrot(i,m,j)
                ftable(i,m) = ftable(i,m) + (thicrot(i,m,j) / (thlqrot(i,m,j) &
                              + thicrot(i,m,j) - thmrot(i,m,j))) * dlzwrot(i,m,j)
                ! else if (isndgat(1,j)==-3) then
                !    actlyr=actlyr+delz(j)
                !    ftable=ftable+delz(j)
              end if
            else if (tbarrot(i,m,j) > tfrez) then
              actlyr(i,m) = actlyr(i,m) + delz(j)
              ftable(i,m) = ftable(i,m) + delz(j)
            end if
          end if
        end do
      end do
    end do

    do i = 1,nltest
      do m = 1,nmtest

        ! Once a year we adjust the maximum annual active layer depth
        ! in an e-folding sense with the
        ! present maxmimum active layer depth for the year ending on the summer
        ! solstice. The maximum values are used in bio2str
        ! to ensure roots are not placed into frozen soil layers.

        if ((dlatrow(i) > 0. .and. iday == 355) & ! Boreal :: winter solstice.
            .or. (dlatrow(i) < 0. .and. iday == 172)) then  ! Austral winter solstice.
          maxAnnualActLyr(i,m) = maxAnnualActLyr(i,m) * efoldfact &
                                 + actLyrThisYr(i,m) * (1.0 - efoldfact)
          actLyrThisYr(i,m) = 0.0
        else
          ! Compare the present active layer depth against actLyrThisYr
          actLyrThisYr(i,m) = max(actLyrThisYr(i,m),actlyr(i,m))
        end if
      end do
    end do

    end associate
  end subroutine findPermafrostVars
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findPermafrostVars
  !! @{
  !> As real :: numbers are not precise, this is a simple way to compare two reals
  !! @author Joe Melton
  !!
  logical function closeEnough (num1, num2, error)

    implicit none

    real, intent(in)     :: num1, num2
    real, intent(in)     :: error
    if (abs(num1 - num2) < error) then
      closeEnough = .true.
    else
      closeEnough = .false.
    end if
  end function closeEnough
  !! @}
  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup generalutils_initRandomSeed
  !! @{
  !! This subroutine sets a repeatable seed for the random number generator.
  !! @author J. Melton
  !!
  subroutine initRandomSeed

    implicit none

    ! NOTE: this subroutine will eventually be replaced by "call random_init()" which
    ! is an intrinsic in the  Fortran 2018 standard.
    integer :: n
    integer, allocatable :: seed(:)

    call random_seed(size = n)

    allocate(seed(n))

    seed = 589389089

    call random_seed(put = seed)

  end subroutine initRandomSeed
  !! @}

  !> \ingroup generalutils_checksumCalc
  !! @{
  !! This subroutine takes the lonIndex and latIndex of a cell, accesses many
  !! attributes of the cell after the run, and creates a checksum from those attributes.
  !! This checksum is written to a .csv file in the output directory, which is then
  !! compared against the checksum from a previous run.
  !! @author M. Fortier
  !!
  subroutine checksumCalc (lonIndex, latIndex)

    use ctemStateVars,  only : c_switch, vrot
    use classStateVars, only : class_rot
    use classicParams,  only : icc, nmos, ignd, icp1, modelpft, TFREZ

    implicit none

    ! arguments
    integer, intent(in) :: lonIndex       !< Index of the longitude of this cell
    integer, intent(in) :: latIndex       !< Index of the latitude of this cell

    ! local variables for this subroutine
    integer                             :: checksum, k !, PFTCover_size
    real, allocatable, dimension(:)     :: PFTCover, Soil, Canopy_real, Snow, CPools, Peatlands, CompetClimate !< Checksum values to track
    integer, allocatable, dimension(:)  :: Canopy_int
    character(len = 10)                 :: lonchar, latchar, checksumchar !< string representations
    character(len = 500)                :: buffer, filename

    ! Associate names with variables defined in derived types.
    associate( &
    ctem_on           => c_switch%ctem_on,              & !< logical: 
    PFTCompetition    => c_switch%PFTCompetition,       & !< logical: 
    lnduseon          => c_switch%lnduseon,             & !< logical: 
    fcancmxrow        => vrot%fcancmx,                  & !< real, dimension(:,:,:) : 
    gleafmasrow       => vrot%gleafmas,                 & !< real, dimension(:,:,:) : 
    bleafmasrow       => vrot%bleafmas,                 & !< real, dimension(:,:,:) : 
    stemmassrow       => vrot%stemmass,                 & !< real, dimension(:,:,:) : 
    rootmassrow       => vrot%rootmass,                 & !< real, dimension(:,:,:) : 
    twarmm            => vrot%twarmm,                   & !< real, dimension(:,:) : 
    tcoldm            => vrot%tcoldm,                   & !< real, dimension(:,:) : 
    gdd5              => vrot%gdd5,                     & !< real, dimension(:,:) : 
    aridity           => vrot%aridity,                  & !< real, dimension(:,:) : aridity index, ratio of potential evaporation to precipitation 
    srplsmon          => vrot%srplsmon,                 & !< real, dimension(:,:) : number of months in a year with surplus water i.e.precipitation more than potential evaporation 
    defctmon          => vrot%defctmon,                 & !< real, dimension(:,:) : number of months in a year with water deficit i.e.precipitation less than potential evaporation 
    anndefct          => vrot%anndefct,                 & !< real, dimension(:,:) : annual water deficit (mm) 
    annsrpls          => vrot%annsrpls,                 & !< real, dimension(:,:) : annual water surplus (mm) 
    annpcp            => vrot%annpcp,                   & !< real, dimension(:,:) : annual precipitation (mm) 
    dry_season_length => vrot%dry_season_length,        & !< real, dimension(:,:) : length of dry season (months) 
    litrmassrow       => vrot%litrmass,                 & !< real, dimension(:,:,:,:) : 
    soilcmasrow       => vrot%soilcmas,                 & !< real, dimension(:,:,:,:) : 
    lfstatusrow       => vrot%lfstatus,                 & !< integer, dimension(:,:,:) : 
    pandaysrow        => vrot%pandays,                  & !< integer, dimension(:,:,:) : 
    Cmossmas          => vrot%Cmossmas,                 & !< real, dimension(:,:) : C in moss biomass, \f$kg C/m^2\f$ 
    litrmsmoss        => vrot%litrmsmoss,               & !< real, dimension(:,:) : moss litter mass, \f$kg C/m^2\f$ 
    dmoss             => vrot%dmoss,                    & !< real, dimension(:,:) : depth of living moss (m) 
    peatSoilC         => vrot%peatSoilC,                & !< real, dimension(:,:) : peat soil C mass, \f$kg C/m^2\f$ 
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
    GROROT            => class_rot%GROROT               & !< real, dimension(:,:)   : 
    )


    write(lonchar, '(I4)')lonIndex
    write(latchar, '(I4)')latIndex
    ! generate the proper, formatted filename
    !ignoreLint(3) (it messes with the file path)
    write(filename, "(A,A,A,'_',A,A)") TRIM(adjustl(c_switch%output_directory)), '/checksums/', &
                                       TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), '.csv'
    open(unit = 500, file = TRIM(adjustl(filename)))


    checksum = 0
    allocate(PFTCover(SIZE(FAREROT) + SIZE(FCANROT) + SIZE(GROROT) + SIZE(fcancmxrow))) !< allocate the array
    PFTCover = (/ pack(FAREROT, .true.),pack(FCANROT, .true.),pack(GROROT, .true.),pack(fcancmxrow, .true.) /)
    do k = 1,SIZE(PFTCover)
      checksum = checksum + bitcount(PFTCover(k))
    end do
    write(checksumchar, '(I4)')checksum   !< transfer to a string variable
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "PFTCover",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Soil(SIZE(THLQROT) + SIZE(THICROT) + SIZE(TPNDROT) + SIZE(ZPNDROT)))
    Soil = (/ pack(THLQROT, .true.),pack(THICROT, .true.),pack(TPNDROT, .true.),pack(ZPNDROT, .true.) /)
    do k = 1,SIZE(Soil)
      checksum = checksum + bitcount(Soil(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Soil",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Canopy_real(SIZE(TCANROT) + SIZE(RCANROT) + SIZE(SCANROT)))
    Canopy_real = (/ pack(TCANROT, .true.),pack(RCANROT, .true.),pack(SCANROT, .true.) /)
    allocate(Canopy_int(SIZE(lfstatusrow) + SIZE(pandaysrow)))
    Canopy_int = (/ pack(lfstatusrow, .true.),pack(pandaysrow, .true.) /)
    do k = 1,SIZE(Canopy_real)
      checksum = checksum + bitcount(Canopy_real(k))
    end do
    do k = 1,SIZE(Canopy_int)
      checksum = checksum + bitcount_int(Canopy_int(k))
    end do

    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Canopy",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Snow(SIZE(TSNOROT) + SIZE(SNOROT) +  SIZE(ALBSROT) + SIZE(RHOSROT)))
    Snow = (/ pack(TSNOROT, .true.),pack(SNOROT, .true.),pack(ALBSROT, .true.),pack(RHOSROT, .true.) /)
    do k = 1,SIZE(Snow)
      checksum = checksum + bitcount(Snow(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Snow",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(CPools(SIZE(gleafmasrow) + SIZE(bleafmasrow) + SIZE(stemmassrow) + SIZE(rootmassrow) + SIZE(litrmassrow) + SIZE(soilcmasrow)))
    CPools = (/ pack(gleafmasrow, .true.),pack(bleafmasrow, .true.),pack(stemmassrow, .true.),pack(rootmassrow, .true.),pack(litrmassrow, .true.),pack(soilcmasrow, .true.) /)
    do k = 1,SIZE(Cpools)
      checksum = checksum + bitcount(Cpools(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Cpools",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Peatlands(SIZE(Cmossmas) + SIZE(litrmsmoss) + SIZE (dmoss) + SIZE(peatSoilC)))
    Peatlands = (/ pack(Cmossmas, .true.),pack(litrmsmoss, .true.),pack(dmoss, .true.) ,pack(peatSoilC, .true.) /)
    do k = 1,SIZE(Peatlands)
      checksum = checksum + bitcount(Peatlands(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Peatlands",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(CompetClimate(SIZE(twarmm) + SIZE(tcoldm) + SIZE(gdd5) + SIZE(aridity) + SIZE(srplsmon) + SIZE(defctmon) + SIZE(anndefct) + SIZE(annsrpls) + SIZE(annpcp) + SIZE(dry_season_length)))
    CompetClimate = (/ pack(twarmm, .true.),pack(tcoldm, .true.),pack(gdd5, .true.),pack(aridity, .true.),pack(srplsmon, .true.),pack(defctmon, .true.),pack(anndefct, .true.),pack(annsrpls, .true.),pack(annpcp, .true.),pack(dry_season_length, .true.) /)
    do k = 1,SIZE(CompetClimate)
      checksum = checksum + bitcount(CompetClimate(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "CompetClimate",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    close(500)

    end associate
  end subroutine checksumCalc
  !! @}

  !> \ingroup generalutils_bitcount
  !! @{
  !! This function generates the bitcount of the specified variable
  !! @author M. Fortier
  !!
  integer function bitcount (scalar)

    implicit none
    real    :: scalar     !< Scalar to be bit-counted
    integer(SELECTED_REAL_KIND(15,307)) :: scalar_int !< Integer with memory representation of 'scalar'

    integer :: bit
    bitcount = 0
    scalar_int = TRANSFER(scalar,scalar_int)
    do bit = 0,BIT_SIZE(scalar_int) - 1
      if (BTEST(scalar_int,bit) ) bitcount = bitcount + 1
    end do

  end function bitcount
  !! @}


  !> \ingroup generalutils_bitcount_int
  !! @{
  !! This function generates the bitcount of the specified variable, from an integer
  !! @author M. Fortier
  !!
  integer function bitcount_int (scalar)

    implicit none
    integer    :: scalar     !< Scalar to be bit-counted
    integer    :: bit

    bitcount_int = 0
    do bit = 0,BIT_SIZE(scalar) - 1
      if (BTEST(scalar,bit) ) bitcount_int = bitcount_int + 1
    end do

  end function bitcount_int
  !! @}
  ! -----------------------------------------------------------------------------------------------

  !> \namespace generalutils
  !> Central module for all general utilities
  !!
  !! The checksum subroutine computes content-based checksums of several groups of variables after
  !! a run has completed. These checksums are used when making non-logical changes to
  !! the code-base. This is accomplished by comparing the computed checksums, to the
  !! checksums of an identical run (same parameters, met files, etc.) before the changes
  !! were made. This is a modified form of regression testing specific to non-logical
  !! software changes.
  !!
  !! Checksums are not an infallible means of ensuring no logical changes, as two numbers
  !! may have different binary representations with the same number of flipped bits. However,
  !! with the number of variables we are checking, it is extremely unlikely to render a false-
  !! positive. If we take a data item of bit-length \f$n\f$, having \f$b\f$ flipped bits in
  !! its representation, the number of same-length data items with the same checksum can
  !! be expressed through the binary coefficient
  !!
  !! \f[ \binom{n}{b} \f] 
  !!
  !! If we assume the worst case where \f$b=\frac{n}{2}\f$, then we are left with
  !!
  !! \f[ \binom{n}{n/2} \f]
  !!
  !! Dividing by the total number of possible values for an n-digit binary number,
  !! we find the probability of a false positive checksum to be, in the worst case:
  !!
  !! \f[ \frac{\binom{n}{n/2}}{2^n} \f]
  !!
  !! Due to the cascading effects of any logical changes to the model, all 7 groups
  !! of variables (6 if peatlands are disabled) will be affected. The smallest of these
  !! groups consists of 3 arrays of numbers for which checksums are computed. With the
  !! conservative assumption that each array contains only a single 32-bit number, the
  !! probability of a false positive (in the absolute worst case) becomes:
  !!
  !! \f[ \large \left({\frac{96\choose 48}{2^{96}}}\right)^{6} ~= 2.87105123*10^{-7}       \f]
  !!

end module generalUtils
