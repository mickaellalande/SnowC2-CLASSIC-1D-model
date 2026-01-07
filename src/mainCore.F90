!> \file
!> Time stepping driver for CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing.
!! @author D. Verseghy, V. Arora, J. Melton, E. Chan
!! This driver program calls the model physics and biogeochemical subroutines,
!! manages the run and the coupling between CLASS and CTEM.

module mainCore

  implicit none

  public :: main_core_driver

contains

  !> \ingroup main_core_driver
  !> @{
  !>
  !! ## Dimension statements.
  !!
  !!     ### First set of definitions:
  !!     Background variables, and prognostic and diagnostic
  !!     variables normally provided by and/or used by the GCM.
  !!     The suffix "rot" refers to variables existing on the
  !!     mosaic grid on the current latitude circle.  The suffix
  !!     "gat" refers to the same variables after they have undergone
  !!     a "gather" operation in which the two mosaic dimensions
  !!     are collapsed into one.  The suffix "row" refers both to
  !!     grid-constant input variables. and to grid-averaged
  !!     diagnostic variables.
  !!
  !!     The first dimension element of the "rot" variables
  !!     refers to the number of grid cells on the current
  !!     latitude circle.  In this stand-alone version, this
  !!     number is set to 1, the second dimension
  !!     element of the "rot" variables refers to the maximum
  !!     number of tiles in the mosaic.  The first
  !!     dimension element in the "gat" variables is given by
  !!     the product of the first two dimension elements in the
  !!     "rot" variables.
  !!
  !!     The majority of CTEM parameters are stored in classicParams.f90.
  !!     Also the CLASS and CTEM variables are stored in modules that we point to
  !!     in this driver. We access the variables and parameters
  !!     through use statements for modules:

  subroutine main_core_driver(nltest, NML, runyr, NDAY, ncount, IDAY, lastDOY, leapnow, N)

    use classicParams,      only : nlat, nmos, ilg, ican, ignd, icc, l2max, NBS, nol2pfts, DELT, TFREZ
    use ctemStateVars,      only : vrot, vgat, c_switch, ctem_tile, tracer
    use classStateVars,     only : class_gat, class_rot, initDiagnosticVars, class_out
    use generalUtils,       only : findDaylength, findLeapYears, findCloudiness
    use classGatherScatter, only : classGather, classScatter, classGatherPrep
    use ctemGatherScatter,  only : ctems2, ctemg1, ctemg2
    use ctemUtilities,      only : dayEndCTEMPreparation, accumulateForCTEM, ctemInit
    use ctemDriver,         only : ctem
    use tracerModule,       only : decay14C
    use applyAllometry,     only : allometry
    use snowProcesses,      only : snowProcessesDriver
    use generalUtils,       only : findPermafrostVars
    use dynamicTiling,      only : moveSplitCopyTiles
    use tiledDisturbance,   only : tiledDisturbancePrep

    implicit none

    ! Arguments
    integer, intent(in) :: nltest  !< Number of grid cells being modelled for this run
    integer, intent(inout) :: NML     !< Counter representing number of mosaic tiles on modelled domain that are land
    integer, intent(in) :: runyr   !< Year of the model run (counts up starting with readMetStartYear continously, even if metLoop > 1)
    integer, intent(in) :: NDAY    !< Number of short (physics) timesteps in one day. e.g., if physics timestep is 15 min this is 48.
    integer, intent(in) :: NCOUNT  !< Counter for daily averaging

    integer, intent(in) :: IDAY    !< Julian day of the year (computed in updateMet, which is not done for CanESM, so must be passed in).
    integer, intent(in) :: lastDOY !< Initialized to 365 days, can be overwritten later if leap = true and it is a leap year.
    logical, intent(in) :: leapnow !< Leap year flag
    integer, intent(inout) :: N    !< Time loop counter


    ! Local variables
    integer :: NMW     !< Counter representing number of mosaic tiles on modelled domain that are lakes
    integer :: NLANDCS !< Number of modelled areas that contain subareas of canopy over snow
    integer :: NLANDGS !< Number of modelled areas that contain subareas of snow over bare ground
    integer :: NLANDC  !< Number of modelled areas that contain subareas of canopy over bare ground
    integer :: NLANDG  !< Number of modelled areas that contain subareas of bare ground
    integer :: NLANDI  !< Number of modelled areas that are ice sheets
    integer :: I, J, K, L, M

    !    =================================================================================

    associate( &

    !nltest            => nlat,                          & !< integer: Number of grid cells being modelled for this run
    nmtest            => nmos,                          & !< integer: Number of mosaic tiles per grid cell being modelled for this run
    NTLD              => nmos,                          & !< integer: Number of mosaic tiles, read in from the initialization file

    ! Associate names with CLASS variables defined in derived types.

    JLAT => NINT(class_rot%latIndexROW(1)),             & !< real, dimension(:) : Index of grid cell being run on the input files grid (latitude)

    ! The following are stored in the data structure: class_gat
    ! they are allocatted in allocClassVars in the classStateVars
    ! module and associated with a different name here.

    DELZ    => class_gat%DELZ,                          & !< real, dimension(ignd) : 
    ZBOT    => class_gat%ZBOT,                          & !< real, dimension(ignd) : 

    ! These are allocated the dimension: 'ilg'

    ILMOS   => class_gat%ILMOS,                         & !< integer, dimension(:) : Index of gridcell corresponding to current element of gathered vector of land surface variables [ ] 
    JLMOS   => class_gat%JLMOS,                         & !< integer, dimension(:) : Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ] 
    IGDRGAT => class_gat%IGDRGAT,                       & !< integer, dimension(:) : Index of soil layer in which bedrock is encountered 
    TZSGAT  => class_gat%TZSGAT,                        & !< real, dimension(:) : Vertical temperature gradient in a snow pack 
    PCSNGAT => class_gat%PCSNGAT,                       & !< real, dimension(:) : Snow fall flux \f$[kg m^{-2} s^{-1} ]\f$ 
    ALBSGAT => class_gat%ALBSGAT,                       & !< real, dimension(:) : Snow albedo [ ] 
    CMAIGAT => class_gat%CMAIGAT,                       & !< real, dimension(:) : Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$ 
    GROGAT  => class_gat%GROGAT,                        & !< real, dimension(:) : Vegetation growth index [ ] 
    QACGAT  => class_gat%QACGAT,                        & !< real, dimension(:) : Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$ 
    RCANGAT => class_gat%RCANGAT,                       & !< real, dimension(:) : Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$ 
    RHOSGAT=> class_gat%RHOSGAT,                        & !< real, dimension(:) : Density of snow \f$[kg m^{-3} ]\f$ 
    SCANGAT=> class_gat%SCANGAT,                        & !< real, dimension(:) : Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$ 
    SNOGAT=> class_gat%SNOGAT,                          & !< real, dimension(:) : Mass of snow pack [kg m^{-2} ]\f$ 
    TACGAT=> class_gat%TACGAT,                          & !< real, dimension(:) : Temperature of air within vegetation canopy [K] 
    TBASGAT=> class_gat%TBASGAT,                        & !< real, dimension(:) : Temperature of bedrock in third soil layer [K] 
    TCANGAT=> class_gat%TCANGAT,                        & !< real, dimension(:) : Vegetation canopy temperature [K] 
    TPNDGAT=> class_gat%TPNDGAT,                        & !< real, dimension(:) : Temperature of ponded water [K] 
    TSNOGAT=> class_gat%TSNOGAT,                        & !< real, dimension(:) : Snowpack temperature [K] 
    TSNBGAT=> class_gat%TSNBGAT,                        & !< real, dimension(:) : Bottom snowpack temperature [K] 
    WSNOGAT=> class_gat%WSNOGAT,                        & !< real, dimension(:) : Liquid water content of snow pack \f$[kg m^{-2} ]\f$ 
    ZPNDGAT=> class_gat%ZPNDGAT,                        & !< real, dimension(:) : Depth of ponded water on surface [m] 
    REFGAT=> class_gat%REFGAT,                          & !< real, dimension(:) : Snow grain size (for ISNOALB=1 option)  [m] 
    BCSNGAT => class_gat%BCSNGAT,                       & !< real, dimension(:) : Black carbon mixing ratio in snow \f$[kg m^{-3} ]\f$ 
    AGIDGAT => class_gat%AGIDGAT,                       & !< real, dimension(:) : Optional user-specified value of ground near-infrared albedo to override CLASS-calculated value [ ] 
    AGVDGAT => class_gat%AGVDGAT,                       & !< real, dimension(:) : Optional user-specified value of ground visible albedo to override CLASS-calculated value [ ] 
    ALGDGAT => class_gat%ALGDGAT,                       & !< real, dimension(:) : Reference albedo for dry soil [ ] 
    ALGWGAT => class_gat%ALGWGAT,                       & !< real, dimension(:) : Reference albedo for saturated soil [ ] 
    ASIDGAT => class_gat%ASIDGAT,                       & !< real, dimension(:) : Optional user-specified value of snow near-infrared albedo to override CLASS-calculated value [ ] 
    ASVDGAT => class_gat%ASVDGAT,                       & !< real, dimension(:) : Optional user-specified value of snow visible albedo to override CLASS-calculated value [ ] 
    DRNGAT => class_gat%DRNGAT,                         & !< real, dimension(:) : Drainage index at bottom of soil profile [ ] 
    GRKFGAT => class_gat%GRKFGAT,                       & !< real, dimension(:) : WATROF parameter used when running MESH code [ ] 
    WFCIGAT => class_gat%WFCIGAT,                       & !< real, dimension(:) : WATROF parameter used when running MESH code [ ] 
    WFSFGAT => class_gat%WFSFGAT,                       & !< real, dimension(:) : WATROF parameter used when running MESH code [ ] 
    XSLPGAT => class_gat%XSLPGAT,                       & !< real, dimension(:) : Surface slope (used when running MESH code) [degrees] 
    ZPLGGAT => class_gat%ZPLGGAT,                       & !< real, dimension(:) : Maximum water ponding depth for snow-free subareas (user-specified when running MESH code) [m] 
    ZPLSGAT => class_gat%ZPLSGAT,                       & !< real, dimension(:) : Maximum water ponding depth for snow-covered subareas (user-specified when running MESH code) [m] 
    ZSNLGAT => class_gat%ZSNLGAT,                       & !< real, dimension(:) : Limiting snow depth below which coverage is < 100% [m] 
    ALGWVGAT => class_gat%ALGWVGAT,                     & !< real, dimension(:) : 
    ALGWNGAT => class_gat%ALGWNGAT,                     & !< real, dimension(:) : 
    ALGDVGAT => class_gat%ALGDVGAT,                     & !< real, dimension(:) : 
    ALGDNGAT => class_gat%ALGDNGAT,                     & !< real, dimension(:) : 
    EMISGAT => class_gat%EMISGAT,                       & !< real, dimension(:) : 
    CSZGAT => class_gat%CSZGAT,                         & !< real, dimension(:) : Cosine of solar zenith angle [ ] 
    DLONGAT => class_gat%DLONGAT,                       & !< real, dimension(:) : Longitude of grid cell (east of Greenwich) [degrees] 
    DLATGAT => class_gat%DLATGAT,                       & !< real, dimension(:) : Latitude of grid cell [degrees] 
    FCLOGAT => class_gat%FCLOGAT,                       & !< real, dimension(:) : Fractional cloud cover [ ] 
    FDLGAT => class_gat%FDLGAT,                         & !< real, dimension(:) : Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$ 
    FSIHGAT => class_gat%FSIHGAT,                       & !< real, dimension(:) : Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$ 
    FSVHGAT => class_gat%FSVHGAT,                       & !< real, dimension(:) : Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$ 
    GGEOGAT => class_gat%GGEOGAT,                       & !< real, dimension(:) : Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$ 
    PADRGAT => class_gat%PADRGAT,                       & !< real, dimension(:) : Partial pressure of dry air [Pa] 
    PREGAT => class_gat%PREGAT,                         & !< real, dimension(:) : Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$ 
    PRESGAT => class_gat%PRESGAT,                       & !< real, dimension(:) : Surface air pressure [Pa] 
    QAGAT => class_gat%QAGAT,                           & !< real, dimension(:) : Specific humidity at reference height \f$[kg kg^{-1} ]\f$ 
    RADJGAT => class_gat%RADJGAT,                       & !< real, dimension(:) : Latitude of grid cell (positive north of equator) [rad] 
    RHOAGAT => class_gat%RHOAGAT,                       & !< real, dimension(:) : Density of air \f$[kg m^{-3} ]\f$ 
    RHSIGAT => class_gat%RHSIGAT,                       & !< real, dimension(:) : Density of fresh snow \f$[kg m^{-3} ]\f$ 
    RPCPGAT => class_gat%RPCPGAT,                       & !< real, dimension(:) : Rainfall rate over modelled area \f$[m s^{-1} ]\f$ 
    SPCPGAT => class_gat%SPCPGAT,                       & !< real, dimension(:) : Snowfall rate over modelled area \f$[m s^{-1} ]\f$ 
    TAGAT => class_gat%TAGAT,                           & !< real, dimension(:) : Air temperature at reference height [K] 
    TADPGAT => class_gat%TADPGAT,                       & !< real, dimension(:) : Dew point temperature of air [K] 
    TRPCGAT => class_gat%TRPCGAT,                       & !< real, dimension(:) : Rainfall temperature [K] 
    TSPCGAT => class_gat%TSPCGAT,                       & !< real, dimension(:) : Snowfall temperature [K] 
    ULGAT => class_gat%ULGAT,                           & !< real, dimension(:) : Zonal component of wind velocity \f$[m s^{-1} ]\f$ 
    VLGAT => class_gat%VLGAT,                           & !< real, dimension(:) : Meridional component of wind velocity \f$[m s^{-1} ]\f$ 
    VMODGAT => class_gat%VMODGAT,                       & !< real, dimension(:) : Wind speed at reference height \f$[m s^{-1} ]\f$ 
    VPDGAT => class_gat%VPDGAT,                         & !< real, dimension(:) : Vapour pressure deficit [mb] 
    Z0ORGAT => class_gat%Z0ORGAT,                       & !< real, dimension(:) : Orographic roughness length [m] 
    ZBLDGAT => class_gat%ZBLDGAT,                       & !< real, dimension(:) : Atmospheric blending height for surface roughness length averaging [m] 
    ZDHGAT => class_gat%ZDHGAT,                         & !< real, dimension(:) : User-specified height associated with diagnosed screen-level variables [m] 
    ZDMGAT => class_gat%ZDMGAT,                         & !< real, dimension(:) : User-specified height associated with diagnosed anemometer-level wind speed [m] 
    ZRFHGAT => class_gat%ZRFHGAT,                       & !< real, dimension(:) : Reference height associated with forcing air temperature and humidity [m] 
    ZRFMGAT => class_gat%ZRFMGAT,                       & !< real, dimension(:) : Reference height associated with forcing wind speed [m] 
    FSGGAT => class_gat%FSGGAT,                         & !< real, dimension(:) : 
    FLGGAT => class_gat%FLGGAT,                         & !< real, dimension(:) : 
    GUSTGAT => class_gat%GUSTGAT,                       & !< real, dimension(:) : 
    DEPBGAT => class_gat%DEPBGAT,                       & !< real, dimension(:) : Black carbon deposition flux \f$[kg m^{-2} s^{-1} ]\f$ 
    GTBS => class_gat%GTBS,                             & !< real, dimension(:) : 
    SFCUBS => class_gat%SFCUBS,                         & !< real, dimension(:) : 
    SFCVBS => class_gat%SFCVBS,                         & !< real, dimension(:) : 
    USTARBS => class_gat%USTARBS,                       & !< real, dimension(:) : 
    TCSNOW => class_gat%TCSNOW,                         & !< real, dimension(:) : Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$ 
    GSNOW => class_gat%GSNOW,                           & !< real, dimension(:) : Diagnostic heat flux at snow surface for use in CCCma black carbon deposition scheme \f$[W m^{-2}]\f$ 
    ALIRGAT => class_gat%ALIRGAT,                       & !< real, dimension(:) : Diagnosed total near-infrared albedo of land surface [ ] 
    ALVSGAT => class_gat%ALVSGAT,                       & !< real, dimension(:) : Diagnosed total visible albedo of land surface [ ] 
    CDHGAT => class_gat%CDHGAT,                         & !< real, dimension(:) : Surface drag coefficient for heat [ ] 
    CDMGAT => class_gat%CDMGAT,                         & !< real, dimension(:) : Surface drag coefficient for momentum [ ] 
    DRGAT => class_gat%DRGAT,                           & !< real, dimension(:) : Surface drag coefficient under neutral stability [ ] 
    EFGAT => class_gat%EFGAT,                           & !< real, dimension(:) : Evaporation efficiency at ground surface [ ] 
    FLGGGAT => class_gat%FLGGGAT,                       & !< real, dimension(:) : Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$ 
    FLGSGAT => class_gat%FLGSGAT,                       & !< real, dimension(:) : Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$ 
    FLGVGAT => class_gat%FLGVGAT,                       & !< real, dimension(:) : Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$ 
    FSGGGAT => class_gat%FSGGGAT,                       & !< real, dimension(:) : Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$ 
    FSGSGAT => class_gat%FSGSGAT,                       & !< real, dimension(:) : Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$ 
    FSGVGAT => class_gat%FSGVGAT,                       & !< real, dimension(:) : Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$ 
    FSNOGAT => class_gat%FSNOGAT,                       & !< real, dimension(:) : Diagnosed fractional snow coverage [ ] 
    GAGAT => class_gat%GAGAT,                           & !< real, dimension(:) : Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$ 
    GTGAT => class_gat%GTGAT,                           & !< real, dimension(:) : Diagnosed effective surface black-body temperature [K] 
    HBLGAT => class_gat%HBLGAT,                         & !< real, dimension(:) : Height of the atmospheric boundary layer [m] 
    HEVCGAT => class_gat%HEVCGAT,                       & !< real, dimension(:) : Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$ 
    HEVGGAT => class_gat%HEVGGAT,                       & !< real, dimension(:) : Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$ 
    HEVSGAT => class_gat%HEVSGAT,                       & !< real, dimension(:) : Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$ 
    HFSGAT => class_gat%HFSGAT,                         & !< real, dimension(:) : Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$ 
    HFSCGAT => class_gat%HFSCGAT,                       & !< real, dimension(:) : Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$ 
    HFSGGAT => class_gat%HFSGGAT,                       & !< real, dimension(:) : Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$ 
    HFSSGAT => class_gat%HFSSGAT,                       & !< real, dimension(:) : Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$ 
    HMFCGAT => class_gat%HMFCGAT,                       & !< real, dimension(:) : Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$ 
    HMFNGAT => class_gat%HMFNGAT,                       & !< real, dimension(:) : Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$ 
    HTCCGAT => class_gat%HTCCGAT,                       & !< real, dimension(:) : Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$ 
    HTCSGAT => class_gat%HTCSGAT,                       & !< real, dimension(:) : Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$ 
    ILMOGAT => class_gat%ILMOGAT,                       & !< real, dimension(:) : Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$ 
    PCFCGAT => class_gat%PCFCGAT,                       & !< real, dimension(:) : Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    PCLCGAT => class_gat%PCLCGAT,                       & !< real, dimension(:) : Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    PCPGGAT => class_gat%PCPGGAT,                       & !< real, dimension(:) : Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$ 
    PCPNGAT => class_gat%PCPNGAT,                       & !< real, dimension(:) : Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$ 
    PETGAT => class_gat%PETGAT,                         & !< real, dimension(:) : Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$ 
    QEVPGAT => class_gat%QEVPGAT,                       & !< real, dimension(:) : Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$ 
    QFCFGAT => class_gat%QFCFGAT,                       & !< real, dimension(:) : Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    QFCLGAT => class_gat%QFCLGAT,                       & !< real, dimension(:) : Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    QFGGAT => class_gat%QFGGAT,                         & !< real, dimension(:) : Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$ 
    QFNGAT => class_gat%QFNGAT,                         & !< real, dimension(:) : Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$ 
    QFSGAT => class_gat%QFSGAT,                         & !< real, dimension(:) : Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    QFXGAT => class_gat%QFXGAT,                         & !< real, dimension(:) : Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$ 
    QGGAT => class_gat%QGGAT,                           & !< real, dimension(:) : Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$ 
    ROFGAT => class_gat%ROFGAT,                         & !< real, dimension(:) : Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$ 
    ! FLAG: ROFBGAT is assigned also to vgat%ROFB (see below) -- go with that one for now since it is also used in ctemUtilities.f90.
    !ROFBGAT => class_gat%ROFBGAT,                       & !< real, dimension(:) : Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROFCGAT => class_gat%ROFCGAT,                       & !< real, dimension(:) : Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    ROFNGAT => class_gat%ROFNGAT,                       & !< real, dimension(:) : Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    ROFOGAT => class_gat%ROFOGAT,                       & !< real, dimension(:) : Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROFSGAT => class_gat%ROFSGAT,                       & !< real, dimension(:) : Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROVGGAT => class_gat%ROVGGAT,                       & !< real, dimension(:) : Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    SFCQGAT => class_gat%SFCQGAT,                       & !< real, dimension(:) : Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    SFCTGAT => class_gat%SFCTGAT,                       & !< real, dimension(:) : Diagnosed screen-level air temperature [K]
    SFCUGAT => class_gat%SFCUGAT,                       & !< real, dimension(:) : Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    SFCVGAT => class_gat%SFCVGAT,                       & !< real, dimension(:) : Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    TFXGAT => class_gat%TFXGAT,                         & !< real, dimension(:) : Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    TROBGAT => class_gat%TROBGAT,                       & !< real, dimension(:) : Temperature of base flow from bottom of soil column [K]
    TROFGAT => class_gat%TROFGAT,                       & !< real, dimension(:) : Temperature of total runoff [K]
    TROOGAT => class_gat%TROOGAT,                       & !< real, dimension(:) : Temperature of overland flow from top of soil column [K]
    TROSGAT => class_gat%TROSGAT,                       & !< real, dimension(:) : Temperature of interflow from sides of soil column [K]
    UEGAT => class_gat%UEGAT,                           & !< real, dimension(:) : Friction velocity of air \f$[m s^{-1} ]\f$
    !WTABGAT => class_gat%WTABGAT,                       & !< real, dimension(:) : Depth of water table in soil [m]
    WTRCGAT => class_gat%WTRCGAT,                       & !< real, dimension(:) : Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    WTRGGAT => class_gat%WTRGGAT,                       & !< real, dimension(:) : Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    WTRSGAT => class_gat%WTRSGAT,                       & !< real, dimension(:) : Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    wtableGAT => class_gat%wtableGAT,                   & !< real, dimension(:) : Depth of water table in soil [m]
    maxAnnualActLyrGAT => class_gat%maxAnnualActLyrGAT, & !< real, dimension(:) : Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    QLWOGAT => class_gat%QLWOGAT,                       & !< real, dimension(:) :
    SFRHGAT => class_gat%SFRHGAT,                       & !< real, dimension(:) :
    FTEMP => class_gat%FTEMP,                           & !< real, dimension(:) :
    FVAP => class_gat%FVAP,                             & !< real, dimension(:) :
    RIB => class_gat%RIB,                               & !< real, dimension(:) :
    FC => class_gat%FC,                                 & !< real, dimension(:) :
    FG => class_gat%FG,                                 & !< real, dimension(:) :
    FCS => class_gat%FCS,                               & !< real, dimension(:) :
    FGS => class_gat%FGS,                               & !< real, dimension(:) :
    RBCOEF => class_gat%RBCOEF,                         & !< real, dimension(:) :
    ZSNOW => class_gat%ZSNOW,                           & !< real, dimension(:) : Depth of snow pack \f$[m] (z_s)\f$
    FSVF => class_gat%FSVF,                             & !< real, dimension(:) :
    FSVFS => class_gat%FSVFS,                           & !< real, dimension(:) :
    ALVSCN => class_gat%ALVSCN,                         & !< real, dimension(:) :
    ALIRCN => class_gat%ALIRCN,                         & !< real, dimension(:) :
    ALVSG => class_gat%ALVSG,                           & !< real, dimension(:) :
    ALIRG => class_gat%ALIRG,                           & !< real, dimension(:) :
    ALVSCS => class_gat%ALVSCS,                         & !< real, dimension(:) :
    ALIRCS => class_gat%ALIRCS,                         & !< real, dimension(:) :
    ALVSSN => class_gat%ALVSSN,                         & !< real, dimension(:) :
    ALIRSN => class_gat%ALIRSN,                         & !< real, dimension(:) :
    ALVSGC => class_gat%ALVSGC,                         & !< real, dimension(:) :
    ALIRGC => class_gat%ALIRGC,                         & !< real, dimension(:) :
    ALVSSC => class_gat%ALVSSC,                         & !< real, dimension(:) :
    ALIRSC => class_gat%ALIRSC,                         & !< real, dimension(:) :
    TRVSCN => class_gat%TRVSCN,                         & !< real, dimension(:) :
    TRIRCN => class_gat%TRIRCN,                         & !< real, dimension(:) :
    TRVSCS => class_gat%TRVSCS,                         & !< real, dimension(:) :
    TRIRCS => class_gat%TRIRCS,                         & !< real, dimension(:) :
    RC => class_gat%RC,                                 & !< real, dimension(:) :
    RCS => class_gat%RCS,                               & !< real, dimension(:) :
    FRAINC => class_gat%FRAINC,                         & !< real, dimension(:) :
    FSNOWC => class_gat%FSNOWC,                         & !< real, dimension(:) :
    FRAICS => class_gat%FRAICS,                         & !< real, dimension(:) :
    FSNOCS => class_gat%FSNOCS,                         & !< real, dimension(:) :
    CMASSC => class_gat%CMASSC,                         & !< real, dimension(:) :
    CMASCS => class_gat%CMASCS,                         & !< real, dimension(:) :
    DISP => class_gat%DISP,                             & !< real, dimension(:) :
    DISPS => class_gat%DISPS,                           & !< real, dimension(:) :
    ZOMLNC => class_gat%ZOMLNC,                         & !< real, dimension(:) :
    ZOELNC => class_gat%ZOELNC,                         & !< real, dimension(:) :
    ZOMLNG => class_gat%ZOMLNG,                         & !< real, dimension(:) :
    ZOELNG => class_gat%ZOELNG,                         & !< real, dimension(:) :
    ZOMLCS => class_gat%ZOMLCS,                         & !< real, dimension(:) :
    ZOELCS => class_gat%ZOELCS,                         & !< real, dimension(:) :
    ZOMLNS => class_gat%ZOMLNS,                         & !< real, dimension(:) :
    ZOELNS => class_gat%ZOELNS,                         & !< real, dimension(:) :
    TRSNOWC => class_gat%TRSNOWC,                       & !< real, dimension(:) :
    CHCAP => class_gat%CHCAP,                           & !< real, dimension(:) : Heat capacity of vegetation canopy \f$[J m^{-2} K^{-1} ] (C_c)\f$
    CHCAPS => class_gat%CHCAPS,                         & !< real, dimension(:) :
    GZEROC => class_gat%GZEROC,                         & !< real, dimension(:) : Vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    GZEROG => class_gat%GZEROG,                         & !< real, dimension(:) : Bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    GZROCS => class_gat%GZROCS,                         & !< real, dimension(:) : Snow-covered vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    GZROGS => class_gat%GZROGS,                         & !< real, dimension(:) : Snow-covered bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    groundHeatFlux => class_gat%groundHeatFlux,         & !< real, dimension(:) : Heat flux at soil surface \f$[W m^{-2} ]\f$
    G12C => class_gat%G12C,                             & !< real, dimension(:) :
    G12G => class_gat%G12G,                             & !< real, dimension(:) :
    G12CS => class_gat%G12CS,                           & !< real, dimension(:) :
    G12GS => class_gat%G12GS,                           & !< real, dimension(:) :
    G23C => class_gat%G23C,                             & !< real, dimension(:) :
    G23G => class_gat%G23G,                             & !< real, dimension(:) :
    G23CS => class_gat%G23CS,                           & !< real, dimension(:) :
    G23GS => class_gat%G23GS,                           & !< real, dimension(:) :
    QFREZC => class_gat%QFREZC,                         & !< real, dimension(:) :
    QFREZG => class_gat%QFREZG,                         & !< real, dimension(:) :
    QMELTC => class_gat%QMELTC,                         & !< real, dimension(:) :
    QMELTG => class_gat%QMELTG,                         & !< real, dimension(:) :
    EVAPC => class_gat%EVAPC,                           & !< real, dimension(:) :
    EVAPCG => class_gat%EVAPCG,                         & !< real, dimension(:) :
    EVAPG => class_gat%EVAPG,                           & !< real, dimension(:) :
    EVAPCS => class_gat%EVAPCS,                         & !< real, dimension(:) :
    EVPCSG => class_gat%EVPCSG,                         & !< real, dimension(:) :
    EVAPGS => class_gat%EVAPGS,                         & !< real, dimension(:) :
    TCANO => class_gat%TCANO,                           & !< real, dimension(:) :
    TCANS => class_gat%TCANS,                           & !< real, dimension(:) :
    RAICAN => class_gat%RAICAN,                         & !< real, dimension(:) :
    SNOCAN => class_gat%SNOCAN,                         & !< real, dimension(:) :
    RAICNS => class_gat%RAICNS,                         & !< real, dimension(:) :
    SNOCNS => class_gat%SNOCNS,                         & !< real, dimension(:) :
    CWLCAP => class_gat%CWLCAP,                         & !< real, dimension(:) :
    CWFCAP => class_gat%CWFCAP,                         & !< real, dimension(:) :
    CWLCPS => class_gat%CWLCPS,                         & !< real, dimension(:) :
    CWFCPS => class_gat%CWFCPS,                         & !< real, dimension(:) :
    TSNOCS => class_gat%TSNOCS,                         & !< real, dimension(:) :
    TSNOGS => class_gat%TSNOGS,                         & !< real, dimension(:) :
    RHOSCS => class_gat%RHOSCS,                         & !< real, dimension(:) :
    RHOSGS => class_gat%RHOSGS,                         & !< real, dimension(:) :
    WSNOCS => class_gat%WSNOCS,                         & !< real, dimension(:) :
    WSNOGS => class_gat%WSNOGS,                         & !< real, dimension(:) :
    TPONDC => class_gat%TPONDC,                         & !< real, dimension(:) :
    TPONDG => class_gat%TPONDG,                         & !< real, dimension(:) :
    TPNDCS => class_gat%TPNDCS,                         & !< real, dimension(:) :
    TPNDGS => class_gat%TPNDGS,                         & !< real, dimension(:) :
    ZPLMCS => class_gat%ZPLMCS,                         & !< real, dimension(:) :
    ZPLMGS => class_gat%ZPLMGS,                         & !< real, dimension(:) :
    ZPLIMC => class_gat%ZPLIMC,                         & !< real, dimension(:) :
    ZPLIMG => class_gat%ZPLIMG,                         & !< real, dimension(:) :
    DSL    => class_gat%DSL,                            & !< real, dimension(:) : Thickness of dry surface layer 
    DSLC   => class_gat%DSLC,                           & !< real, dimension(:) : Thickness of dry surface layer under canopy 
    RB     => class_gat%RB,                             & !< real, dimension(:) : Leaf boundary resistance of vegetation 
    LAIPAIRatio => class_gat%LAIPAIRatio,               & !< real, dimension(:) : LAI to PAI ratio for canopy over bare ground [] 
    LAISPAISRatio => class_gat%LAISPAISRatio,           & !< real, dimension(:) : LAI to PAI ratio for canopy over snow [] 
    !
    !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
    !     * BALANCES.
    !
    CTVSTP => class_gat%CTVSTP,                         & !< real, dimension(:) : 
    CTSSTP => class_gat%CTSSTP,                         & !< real, dimension(:) : 
    CT1STP => class_gat%CT1STP,                         & !< real, dimension(:) : 
    CT2STP => class_gat%CT2STP,                         & !< real, dimension(:) : 
    CT3STP => class_gat%CT3STP,                         & !< real, dimension(:) : 
    WTVSTP => class_gat%WTVSTP,                         & !< real, dimension(:) : 
    WTSSTP => class_gat%WTSSTP,                         & !< real, dimension(:) : 
    WTGSTP => class_gat%WTGSTP,                         & !< real, dimension(:) : 

    ! These are allocated the dimension: 'ilg,ignd'

    ISNDGAT => class_gat%ISNDGAT,                       & !< integer, dimension(:,:) : Integer identifier associated with sand content 
    TBARGAT => class_gat%TBARGAT,                       & !< real, dimension(:,:) : Temperature of soil layers [K] 
    TCTOGAT => class_gat%TCTOGAT,                       & !< real, dimension(:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$ 
    TCBOGAT => class_gat%TCBOGAT,                       & !< real, dimension(:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$ 
    THICGAT => class_gat%THICGAT,                       & !< real, dimension(:,:) : Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    THLQGAT => class_gat%THLQGAT,                       & !< real, dimension(:,:) : Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    BIGAT => class_gat%BIGAT,                           & !< real, dimension(:,:) : Clapp and Hornberger empirical “b” parameter [ ] 
    DLZWGAT => class_gat%DLZWGAT,                       & !< real, dimension(:,:) : Permeable thickness of soil layer [m] 
    GRKSGAT => class_gat%GRKSGAT,                       & !< real, dimension(:,:) : Saturated hydraulic conductivity of soil layers \f$[m s^{-1} ]\f$ 
    HCPSGAT => class_gat%HCPSGAT,                       & !< real, dimension(:,:) : Volumetric heat capacity of soil particles \f$[J m^{-3} ]\f$ 
    PSISGAT => class_gat%PSISGAT,                       & !< real, dimension(:,:) : Soil moisture suction at saturation [m] 
    PSIWGAT => class_gat%PSIWGAT,                       & !< real, dimension(:,:) : Soil moisture suction at wilting point [m] 
    TCSGAT => class_gat%TCSGAT,                         & !< real, dimension(:,:) : Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$ 
    THFCGAT => class_gat%THFCGAT,                       & !< real, dimension(:,:) : Field capacity \f$[m^3 m^{-3} ]\f$ 
    THMGAT => class_gat%THMGAT,                         & !< real, dimension(:,:) : Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$ 
    THPGAT => class_gat%THPGAT,                         & !< real, dimension(:,:) : Pore volume in soil layer \f$[m^3 m^{-3} ]\f$ 
    THRGAT => class_gat%THRGAT,                         & !< real, dimension(:,:) : Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$ 
    THRAGAT => class_gat%THRAGAT,                       & !< real, dimension(:,:) : Fractional saturation of soil behind the wetting front [ ] 
    ZBTWGAT => class_gat%ZBTWGAT,                       & !< real, dimension(:,:) : Depth to permeable bottom of soil layer [m] 
    THLWGAT => class_gat%THLWGAT,                       & !< real, dimension(:,:) : Soil water content at wilting point, \f$[m^3 m^{-3} ]\f$ 
    GFLXGAT => class_gat%GFLXGAT,                       & !< real, dimension(:,:) : Heat conduction between soil layers \f$[W m^{-2} ]\f$ 
    HMFGGAT => class_gat%HMFGGAT,                       & !< real, dimension(:,:) : Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$ 
    HTCGAT => class_gat%HTCGAT,                         & !< real, dimension(:,:) : Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$ 
    QFCGAT => class_gat%QFCGAT,                         & !< real, dimension(:,:) : Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$ 
    TBARC => class_gat%TBARC,                           & !< real, dimension(:,:) : Temperature of soil layers for ground under canopy subarea [K] 
    TBARG => class_gat%TBARG,                           & !< real, dimension(:,:) : Temperature of soil layers for bareground subarea [K] 
    TBARCS => class_gat%TBARCS,                         & !< real, dimension(:,:) : Temperature of soil layers for snow-covered ground under canopy subarea [K] 
    TBARGS => class_gat%TBARGS,                         & !< real, dimension(:,:) : Temperature of soil layers for bareground subarea [K] 
    THLIQC => class_gat%THLIQC,                         & !< real, dimension(:,:) : Volumetric liquid water content of soil layers for ground under canopy subarea \f$[m^3 m^{-3} ]\f$ 
    THLIQG => class_gat%THLIQG,                         & !< real, dimension(:,:) : Volumetric liquid water content of soil layers for bareground subarea \f$[m^3 m^{-3} ]\f$ 
    THICEC => class_gat%THICEC,                         & !< real, dimension(:,:) : Volumetric ice content of soil layers for ground under canopy subarea \f$[m^3 m^{-3} ]\f$ 
    THICEG => class_gat%THICEG,                         & !< real, dimension(:,:) : Volumetric ice content of soil layers for bareground subarea \f$[m^3 m^{-3} ]\f$ 
    FROOT => class_gat%FROOT,                           & !< real, dimension(:,:) : 
    HCPC => class_gat%HCPC,                             & !< real, dimension(:,:) : 
    HCPG => class_gat%HCPG,                             & !< real, dimension(:,:) : 
    FROOTS => class_gat%FROOTS,                         & !< real, dimension(:,:) : 
    TCTOPC => class_gat%TCTOPC,                         & !< real, dimension(:,:) : 
    TCBOTC => class_gat%TCBOTC,                         & !< real, dimension(:,:) : 
    TCTOPG => class_gat%TCTOPG,                         & !< real, dimension(:,:) : 
    TCBOTG => class_gat%TCBOTG,                         & !< real, dimension(:,:) : 

    ! These are allocated the dimension: 'ilg,ican'

    ACIDGAT => class_gat%ACIDGAT,                       & !< real, dimension(:,:) : Optional user-specified value of canopy near-infrared albedo to override CLASS-calculated value [ ] 
    ACVDGAT => class_gat%ACVDGAT,                       & !< real, dimension(:,:) : Optional user-specified value of canopy visible albedo to override CLASS-calculated value [ ] 
    CMASGAT => class_gat%CMASGAT,                       & !< real, dimension(:,:) : Maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$ 
    HGTDGAT => class_gat%HGTDGAT,                       & !< real, dimension(:,:) : Optional user-specified values of height of vegetation categories to override CLASS-calculated values [m] 
    PAIDGAT => class_gat%PAIDGAT,                       & !< real, dimension(:,:) : Optional user-specified value of plant area indices of vegetation categories to override CLASS-calculated values [ ] 
    PAMNGAT => class_gat%PAMNGAT,                       & !< real, dimension(:,:) : Minimum plant area index of vegetation category [ ] 
    PAMXGAT => class_gat%PAMXGAT,                       & !< real, dimension(:,:) : Minimum plant area index of vegetation category [ ] 
    PSGAGAT => class_gat%PSGAGAT,                       & !< real, dimension(:,:) : Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ] 
    PSGBGAT => class_gat%PSGBGAT,                       & !< real, dimension(:,:) : Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ] 
    QA50GAT => class_gat%QA50GAT,                       & !< real, dimension(:,:) : Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$ 
    ROOTGAT => class_gat%ROOTGAT,                       & !< real, dimension(:,:) : Maximum rooting depth of vegetation category [m] 
    RSMNGAT => class_gat%RSMNGAT,                       & !< real, dimension(:,:) : Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$ 
    VPDAGAT => class_gat%VPDAGAT,                       & !< real, dimension(:,:) : Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ] 
    VPDBGAT => class_gat%VPDBGAT,                       & !< real, dimension(:,:) : Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ] 

    ! These are allocated the dimension: 'ilg,icp1'

    ALICGAT => class_gat%ALICGAT,                       & !< real, dimension(:,:) : Background average near-infrared albedo of vegetation category [ ] 
    ALVCGAT => class_gat%ALVCGAT,                       & !< real, dimension(:,:) : Background average visible albedo of vegetation category [ ] 
    FCANGAT => class_gat%FCANGAT,                       & !< real, dimension(:,:) : Maximum fractional coverage of modelled area by vegetation category [ ] 
    LNZ0GAT => class_gat%LNZ0GAT,                       & !< real, dimension(:,:) : Natural logarithm of maximum roughness length of vegetation category [ ] 

    ! These are allocated the dimension: 'ilg,nbs'

    FSDBGAT => class_gat%FSDBGAT,                       & !< real, dimension(:,:) : 
    FSFBGAT => class_gat%FSFBGAT,                       & !< real, dimension(:,:) : 
    FSSBGAT => class_gat%FSSBGAT,                       & !< real, dimension(:,:) : 
    SALBGAT => class_gat%SALBGAT,                       & !< real, dimension(:,:) : 
    CSALGAT => class_gat%CSALGAT,                       & !< real, dimension(:,:) : 
    ALTG => class_gat%ALTG,                             & !< real, dimension(:,:) : 
    ALSNO => class_gat%ALSNO,                           & !< real, dimension(:,:) : 
    TRSNOWG => class_gat%TRSNOWG,                       & !< real, dimension(:,:) : 

    ! These are allocated the dimension: 'ilg,4'
    TSFSGAT => class_gat%TSFSGAT,                       & !< real, dimension(:,:) : Ground surface temperature over subarea [K] 

    ! These are allocated the dimension: 'ilg,6,50'
    ITCTGAT => class_gat%ITCTGAT,                       & !< integer, dimension(:,:,:) : Counter of number of iterations required to solve surface energy balance for the elements of the four subareas 

    ! The following are stored in the data structure: class_rot
    ! they are allocatted in allocClassVars in the classStateVars
    ! module and are associated with different names here.

    ! These are allocated the dimension: 'nlat'

    CSZROW => class_rot%CSZROW,                         & !< real, dimension(:) : Cosine of solar zenith angle [ ] 
    DLONROW => class_rot%DLONROW,                       & !< real, dimension(:) : 
    DLATROW => class_rot%DLATROW,                       & !< real, dimension(:) : 
    FCLOROW => class_rot%FCLOROW,                       & !< real, dimension(:) : Fractional cloud cover [ ] 
    FDLROT => class_rot%FDLROT,                         & !< real, dimension(:) : 
    FSIHROW => class_rot%FSIHROW,                       & !< real, dimension(:,:) : Near infrared shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    FSVHROW => class_rot%FSVHROW,                       & !< real, dimension(:,:) : Visible shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    GGEOROW => class_rot%GGEOROW,                       & !< real, dimension(:) : The geothermal heat flux 
    PADRROW => class_rot%PADRROW,                       & !< real, dimension(:) : 
    PREROW => class_rot%PREROW,                         & !< real, dimension(:) : Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$ 
    PRESROW => class_rot%PRESROW,                       & !< real, dimension(:) : 
    QAROW => class_rot%QAROW,                           & !< real, dimension(:) : 
    RADJROW => class_rot%RADJROW,                       & !< real, dimension(:) : Latitude of grid cell (positive north of equator) [rad] 
    RHOAROW => class_rot%RHOAROW,                       & !< real, dimension(:) : 
    RHSIROW => class_rot%RHSIROW,                       & !< real, dimension(:) : 
    RPCPROW => class_rot%RPCPROW,                       & !< real, dimension(:) : 
    RPREROW => class_rot%RPREROW,                       & !< real, dimension(:) : Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    SPCPROW => class_rot%SPCPROW,                       & !< real, dimension(:) : 
    SPREROW => class_rot%SPREROW,                       & !< real, dimension(:) : Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    TAROW => class_rot%TAROW,                           & !< real, dimension(:) : 
    TADPROW => class_rot%TADPROW,                       & !< real, dimension(:) : 
    TRPCROW => class_rot%TRPCROW,                       & !< real, dimension(:) : 
    TSPCROW => class_rot%TSPCROW,                       & !< real, dimension(:) : 
    ULROW => class_rot%ULROW,                           & !< real, dimension(:) : Zonal component of wind velocity \f$[m s^{-1} ]\f$
    VLROW => class_rot%VLROW,                           & !< real, dimension(:) : Meridional component of wind velocity \f$[m s^{-1} ]\f$
    VMODROW => class_rot%VMODROW,                       & !< real, dimension(:) : Wind speed at reference height \f$[m s^{-1} ]\f$
    VPDROW => class_rot%VPDROW,                         & !< real, dimension(:) : 
    ZBLDROW => class_rot%ZBLDROW,                       & !< real, dimension(:) : 
    ZDHROW => class_rot%ZDHROW,                         & !< real, dimension(:) : 
    ZDMROW => class_rot%ZDMROW,                         & !< real, dimension(:) : 
    ZRFHROW => class_rot%ZRFHROW,                       & !< real, dimension(:) : 
    ZRFMROW => class_rot%ZRFMROW,                       & !< real, dimension(:) : 
    UVROW => class_rot%UVROW,                           & !< real, dimension(:) : Wind speed at reference height as read in from file \f$[m s^{-1} ]\f$
    XDIFFUS => class_rot%XDIFFUS,                       & !< real, dimension(:,:) : Fraction of diffused radiation 
    Z0ORROW => class_rot%Z0ORROW,                       & !< real, dimension(:) : The orographic roughness length 
    PRENROW => class_rot%PRENROW,                       & !< real, dimension(:) : 
    FSGROL => class_rot%FSGROL,                         & !< real, dimension(:) : 
    FLGROL => class_rot%FLGROL,                         & !< real, dimension(:) : 
    GUSTROL => class_rot%GUSTROL,                       & !< real, dimension(:) : 
    DEPBROW => class_rot%DEPBROW,                       & !< real, dimension(:) : Black carbon deposition flux \f$[kg m^{-2} s^{-1} ]\f$
    BCSNROW => class_rot%BCSNROW,                       & !< real, dimension(:) : Black carbon mixing ratio in snow \f$[kg m^{-3} ]\f$ 
    ! These are allocated the dimension: 'nlat,nmos'

    IGDRROT => class_rot%IGDRROT,                       & !< integer, dimension(:,:) : 
    ALBSROT => class_rot%ALBSROT,                       & !< real, dimension(:,:) : 
    CMAIROT => class_rot%CMAIROT,                       & !< real, dimension(:,:) : 
    GROROT => class_rot%GROROT,                         & !< real, dimension(:,:) : 
    QACROT => class_rot%QACROT,                         & !< real, dimension(:,:) : 
    RCANROT => class_rot%RCANROT,                       & !< real, dimension(:,:) : 
    RHOSROT => class_rot%RHOSROT,                       & !< real, dimension(:,:) : 
    TCSNROT => class_rot%TCSNROT,                       & !< real, dimension(:,:) : Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
    SCANROT => class_rot%SCANROT,                       & !< real, dimension(:,:) : 
    SNOROT => class_rot%SNOROT,                         & !< real, dimension(:,:) : 
    TACROT => class_rot%TACROT,                         & !< real, dimension(:,:) : 
    TBASROT => class_rot%TBASROT,                       & !< real, dimension(:,:) : 
    TCANROT => class_rot%TCANROT,                       & !< real, dimension(:,:) : 
    TPNDROT => class_rot%TPNDROT,                       & !< real, dimension(:,:) : 
    TSNOROT => class_rot%TSNOROT,                       & !< real, dimension(:,:) : Snowpack temperature [K]
    TSNBROT => class_rot%TSNBROT,                       & !< real, dimension(:,:) : Bottom snowpack temperature [K]
    WSNOROT => class_rot%WSNOROT,                       & !< real, dimension(:,:) : 
    ZPNDROT => class_rot%ZPNDROT,                       & !< real, dimension(:,:) : 
    REFROT => class_rot%REFROT,                         & !< real, dimension(:,:) : Snow grain size (for ISNOALB=1 option)  [m] 
    AGIDROT => class_rot%AGIDROT,                       & !< real, dimension(:,:) : 
    AGVDROT => class_rot%AGVDROT,                       & !< real, dimension(:,:) : 
    ALGDROT => class_rot%ALGDROT,                       & !< real, dimension(:,:) : 
    ALGWROT => class_rot%ALGWROT,                       & !< real, dimension(:,:) : 
    ASIDROT => class_rot%ASIDROT,                       & !< real, dimension(:,:) : 
    ASVDROT => class_rot%ASVDROT,                       & !< real, dimension(:,:) : 
    DRNROT => class_rot%DRNROT,                         & !< real, dimension(:,:) : 
    FAREROT => class_rot%FAREROT,                       & !< real, dimension(:,:) : Fractional coverage of mosaic tile on modelled area 
    GRKFROT => class_rot%GRKFROT,                       & !< real, dimension(:,:) : 
    WFCIROT => class_rot%WFCIROT,                       & !< real, dimension(:,:) : 
    WFSFROT => class_rot%WFSFROT,                       & !< real, dimension(:,:) : 
    XSLPROT => class_rot%XSLPROT,                       & !< real, dimension(:,:) : 
    ZPLGROT => class_rot%ZPLGROT,                       & !< real, dimension(:,:) : 
    ZPLSROT => class_rot%ZPLSROT,                       & !< real, dimension(:,:) : 
    ZSNLROT => class_rot%ZSNLROT,                       & !< real, dimension(:,:) : Limiting snow depth (m) 
    ZSNOROT => class_rot%ZSNOROT,                       & !< real, dimension(:,:) : 
    ALGWVROT => class_rot%ALGWVROT,                     & !< real, dimension(:,:) : 
    ALGWNROT => class_rot%ALGWNROT,                     & !< real, dimension(:,:) : 
    ALGDVROT => class_rot%ALGDVROT,                     & !< real, dimension(:,:) : 
    ALGDNROT => class_rot%ALGDNROT,                     & !< real, dimension(:,:) : 
    EMISROT => class_rot%EMISROT,                       & !< real, dimension(:,:) : 
    ALIRROT => class_rot%ALIRROT,                       & !< real, dimension(:,:) : 
    ALVSROT => class_rot%ALVSROT,                       & !< real, dimension(:,:) : 
    CDHROT => class_rot%CDHROT,                         & !< real, dimension(:,:) : 
    CDMROT => class_rot%CDMROT,                         & !< real, dimension(:,:) : 
    DRROT => class_rot%DRROT,                           & !< real, dimension(:,:) : 
    EFROT => class_rot%EFROT,                           & !< real, dimension(:,:) : 
    FLGGROT => class_rot%FLGGROT,                       & !< real, dimension(:,:) : 
    FLGSROT => class_rot%FLGSROT,                       & !< real, dimension(:,:) : 
    FLGVROT => class_rot%FLGVROT,                       & !< real, dimension(:,:) : 
    FSGGROT => class_rot%FSGGROT,                       & !< real, dimension(:,:) : 
    FSGSROT => class_rot%FSGSROT,                       & !< real, dimension(:,:) : 
    FSGVROT => class_rot%FSGVROT,                       & !< real, dimension(:,:) : 
    FSNOROT => class_rot%FSNOROT,                       & !< real, dimension(:,:) : 
    GAROT => class_rot%GAROT,                           & !< real, dimension(:,:) : 
    GTROT => class_rot%GTROT,                           & !< real, dimension(:,:) : Diagnosed effective surface black-body temperature [K] 
    HBLROT => class_rot%HBLROT,                         & !< real, dimension(:,:) : 
    HEVCROT => class_rot%HEVCROT,                       & !< real, dimension(:,:) : 
    HEVGROT => class_rot%HEVGROT,                       & !< real, dimension(:,:) : 
    HEVSROT => class_rot%HEVSROT,                       & !< real, dimension(:,:) : 
    HFSROT => class_rot%HFSROT,                         & !< real, dimension(:,:) : 
    HFSCROT => class_rot%HFSCROT,                       & !< real, dimension(:,:) : 
    HFSGROT => class_rot%HFSGROT,                       & !< real, dimension(:,:) : 
    HFSSROT => class_rot%HFSSROT,                       & !< real, dimension(:,:) : 
    HMFCROT => class_rot%HMFCROT,                       & !< real, dimension(:,:) : 
    HMFNROT => class_rot%HMFNROT,                       & !< real, dimension(:,:) : 
    HTCCROT => class_rot%HTCCROT,                       & !< real, dimension(:,:) : 
    SDEPROT => class_rot%SDEPROT,                       & !< real, dimension(:,:) : Depth to bedrock in the soil profile 
    SOCIROT => class_rot%SOCIROT,                       & !< real, dimension(:,:) : 
    HTCSROT => class_rot%HTCSROT,                       & !< real, dimension(:,:) : 
    ILMOROT => class_rot%ILMOROT,                       & !< real, dimension(:,:) : 
    PCFCROT => class_rot%PCFCROT,                       & !< real, dimension(:,:) : 
    PCLCROT => class_rot%PCLCROT,                       & !< real, dimension(:,:) : 
    PCPGROT => class_rot%PCPGROT,                       & !< real, dimension(:,:) : 
    PCPNROT => class_rot%PCPNROT,                       & !< real, dimension(:,:) : 
    PETROT => class_rot%PETROT,                         & !< real, dimension(:,:) : 
    QEVPROT => class_rot%QEVPROT,                       & !< real, dimension(:,:) : 
    QFCFROT => class_rot%QFCFROT,                       & !< real, dimension(:,:) : 
    QFCLROT => class_rot%QFCLROT,                       & !< real, dimension(:,:) : 
    QFGROT => class_rot%QFGROT,                         & !< real, dimension(:,:) : 
    QFNROT => class_rot%QFNROT,                         & !< real, dimension(:,:) : 
    QFSROT => class_rot%QFSROT,                         & !< real, dimension(:,:) : 
    QFXROT => class_rot%QFXROT,                         & !< real, dimension(:,:) : 
    QGROT => class_rot%QGROT,                           & !< real, dimension(:,:) : 
    ROFROT => class_rot%ROFROT,                         & !< real, dimension(:,:) : 
    ROFBROT => class_rot%ROFBROT,                       & !< real, dimension(:,:) : 
    ROFCROT => class_rot%ROFCROT,                       & !< real, dimension(:,:) : 
    ROFNROT => class_rot%ROFNROT,                       & !< real, dimension(:,:) : 
    ROFOROT => class_rot%ROFOROT,                       & !< real, dimension(:,:) : 
    ROFSROT => class_rot%ROFSROT,                       & !< real, dimension(:,:) : 
    ROVGROT => class_rot%ROVGROT,                       & !< real, dimension(:,:) : 
    SFCQROT => class_rot%SFCQROT,                       & !< real, dimension(:,:) : 
    SFCTROT => class_rot%SFCTROT,                       & !< real, dimension(:,:) : 
    SFCUROT => class_rot%SFCUROT,                       & !< real, dimension(:,:) : 
    SFCVROT => class_rot%SFCVROT,                       & !< real, dimension(:,:) : 
    TFXROT => class_rot%TFXROT,                         & !< real, dimension(:,:) : 
    TROBROT => class_rot%TROBROT,                       & !< real, dimension(:,:) : 
    TROFROT => class_rot%TROFROT,                       & !< real, dimension(:,:) : 
    TROOROT => class_rot%TROOROT,                       & !< real, dimension(:,:) : 
    TROSROT => class_rot%TROSROT,                       & !< real, dimension(:,:) : 
    UEROT => class_rot%UEROT,                           & !< real, dimension(:,:) : 
    WTRCROT => class_rot%WTRCROT,                       & !< real, dimension(:,:) : 
    WTRGROT => class_rot%WTRGROT,                       & !< real, dimension(:,:) : 
    WTRSROT => class_rot%WTRSROT,                       & !< real, dimension(:,:) : 
    SFRHROT => class_rot%SFRHROT,                       & !< real, dimension(:,:) : 
    wtableROT => class_rot%wtableROT,                   & !< real, dimension(:,:) : Depth of water table in soil [m] 
    maxAnnualActLyrROT => class_rot%maxAnnualActLyrROT, & !< real, dimension(:,:) : Active layer depth maximum over the e-folding period specified by parameter eftime (m). 
    groundHeatFluxROT => class_rot%groundHeatFluxROT,   & !< real, dimension(:,:) : Heat flux at soil surface \f$[W m^{-2} ]\f$ 

    ! There are allocated the dimension: 'nlat,nmos,ignd'

    ISNDROT => class_rot%ISNDROT,                       & !< integer, dimension(:,:,:) : 
    TBARROT => class_rot%TBARROT,                       & !< real, dimension(:,:,:) : 
    TCTOROT => class_rot%TCTOROT,                       & !< real, dimension(:,:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$ 
    TCBOROT => class_rot%TCBOROT,                       & !< real, dimension(:,:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$ 
    THICROT => class_rot%THICROT,                       & !< real, dimension(:,:,:) : 
    THLQROT => class_rot%THLQROT,                       & !< real, dimension(:,:,:) : 
    BIROT => class_rot%BIROT,                           & !< real, dimension(:,:,:) : 
    DLZWROT => class_rot%DLZWROT,                       & !< real, dimension(:,:,:) : Permeable thickness of soil layer [m] 
    GRKSROT => class_rot%GRKSROT,                       & !< real, dimension(:,:,:) : 
    HCPSROT => class_rot%HCPSROT,                       & !< real, dimension(:,:,:) : 
    SANDROT => class_rot%SANDROT,                       & !< real, dimension(:,:,:) : Percentage sand content of soil 
    CLAYROT => class_rot%CLAYROT,                       & !< real, dimension(:,:,:) : Percentage clay content of soil 
    ORGMROT => class_rot%ORGMROT,                       & !< real, dimension(:,:,:) : Percentage organic matter content of soil 
    PSISROT => class_rot%PSISROT,                       & !< real, dimension(:,:,:) : 
    PSIWROT => class_rot%PSIWROT,                       & !< real, dimension(:,:,:) : 
    TCSROT => class_rot%TCSROT,                         & !< real, dimension(:,:,:) : 
    THFCROT => class_rot%THFCROT,                       & !< real, dimension(:,:,:) : 
    THMROT => class_rot%THMROT,                         & !< real, dimension(:,:,:) : Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$ 
    THPROT => class_rot%THPROT,                         & !< real, dimension(:,:,:) : 
    THRROT => class_rot%THRROT,                         & !< real, dimension(:,:,:) : 
    THRAROT => class_rot%THRAROT,                       & !< real, dimension(:,:,:) : 
    ZBTWROT => class_rot%ZBTWROT,                       & !< real, dimension(:,:,:) : 
    THLWROT => class_rot%THLWROT,                       & !< real, dimension(:,:,:) : 
    GFLXROT => class_rot%GFLXROT,                       & !< real, dimension(:,:,:) : 
    HMFGROT => class_rot%HMFGROT,                       & !< real, dimension(:,:,:) : 
    HTCROT => class_rot%HTCROT,                         & !< real, dimension(:,:,:) : 
    QFCROT => class_rot%QFCROT,                         & !< real, dimension(:,:,:) : 

    ! These are allocated the dimension: 'nlat,nmos,ican'

    ACIDROT => class_rot%ACIDROT,                       & !< real, dimension(:,:,:) : 
    ACVDROT => class_rot%ACVDROT,                       & !< real, dimension(:,:,:) : 
    CMASROT => class_rot%CMASROT,                       & !< real, dimension(:,:,:) : 
    HGTDROT => class_rot%HGTDROT,                       & !< real, dimension(:,:,:) : 
    PAIDROT => class_rot%PAIDROT,                       & !< real, dimension(:,:,:) : 
    PAMNROT => class_rot%PAMNROT,                       & !< real, dimension(:,:,:) : 
    PAMXROT => class_rot%PAMXROT,                       & !< real, dimension(:,:,:) : 
    PSGAROT => class_rot%PSGAROT,                       & !< real, dimension(:,:,:) : 
    PSGBROT => class_rot%PSGBROT,                       & !< real, dimension(:,:,:) : 
    QA50ROT => class_rot%QA50ROT,                       & !< real, dimension(:,:,:) : 
    ROOTROT => class_rot%ROOTROT,                       & !< real, dimension(:,:,:) : 
    RSMNROT => class_rot%RSMNROT,                       & !< real, dimension(:,:,:) : 
    VPDAROT => class_rot%VPDAROT,                       & !< real, dimension(:,:,:) : 
    VPDBROT => class_rot%VPDBROT,                       & !< real, dimension(:,:,:) : 
    
    ! These are allocated the dimension: 'nlat,nmos,icp1'
    ALICROT => class_rot%ALICROT,                       & !< real, dimension(:,:,:) : 
    ALVCROT => class_rot%ALVCROT,                       & !< real, dimension(:,:,:) : 
    FCANROT => class_rot%FCANROT,                       & !< real, dimension(:,:,:) : 
    LNZ0ROT => class_rot%LNZ0ROT,                       & !< real, dimension(:,:,:) : 

    ! These are allocated the dimension: 'nlat,nmos,nbs'
    SALBROT => class_rot%SALBROT,                       & !< real, dimension(:,:,:)  : 
    CSALROT => class_rot%CSALROT,                       & !< real, dimension(:,:,:)  : 
    FSDBROL => class_rot%FSDBROL,                       & !< real(nlat,nmos,nbs) : Direct solar radiation in each modelled wavelength band  [W m-2]
    FSFBROL => class_rot%FSFBROL,                       & !< real(nlat,nmos,nbs) : Diffuse solar radiation in each modelled wavelength band  [W m-2]
    FSSBROL => class_rot%FSSBROL,                       & !< real(nlat,nmos,nbs) : Total solar radiation in each modelled wavelength band  [W m-2]

    ! These are allocated the dimension: 'nlat,nmos,6,50'
    ITCTROT => class_rot%ITCTROT,                       & !< integer, dimension(:,:,:,:) : 

    ! These are allocated the dimension: 'nlat,nmos,4'
    TSFSROT => class_rot%TSFSROT,                       & !< real, dimension(:,:,:)  : 

    ! Associate names with CTEM variables defined in derived types.

    spinfast          => c_switch%spinfast,             & !< integer: set this to a higher number up to 10 to spin up
    useTracer         => c_switch%useTracer,            & !< integer: Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
                                                          !!          useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
                                                          !!          the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                                          !!          useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                                          !!          useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme. 
    useStaticPeatDep  => c_switch%useStaticPeatDep,     & !< True will keep the peat depth at the sdep value. This is used for spinup to ensure
                                                          !! that the disequilibrium between the peatland and the driving climate, time period 
                                                          !! since the peatland initiation, or sub-grid heterogeneity impacts don't adversely impact the ability 
                                                          !! of the model to spinup.   
    fixedYearOBSWETF  => c_switch%fixedYearOBSWETF,     & !< integer: set the year to use for observed wetland fraction if transientOBSWETF is false. 
    IDISP             => c_switch%IDISP,                & !< integer: if idisp=0, vegetation displacement heights are ignored,
                                                          !!          because the atmospheric model considers these to be part
                                                          !!          of the "terrain".
                                                          !!          if idisp=1, vegetation displacement heights are calculated. 
    IZREF             => c_switch%IZREF,                & !< integer: if izref=1, the bottom of the atmospheric model is taken
                                                          !!          to lie at the ground surface.
                                                          !!          if izref=2, the bottom of the atmospheric model is taken
                                                          !!          to lie at the local roughness height. 
    ISLFD             => c_switch%ISLFD,                & !< integer: if islfd=0, drcoef is called for surface stability corrections
                                                          !!          and the original gcm set of screen-level diagnostic calculations
                                                          !!          is done.
                                                          !!          if islfd=1, drcoef is called for surface stability corrections
                                                          !!          and sldiag is called for screen-level diagnostic calculations. 
    IPCP              => c_switch%IPCP,                 & !< integer: if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 c.
                                                          !!          if ipcp=2, a linear partitioning of precipitation betweeen
                                                          !!          rainfall and snowfall is done between 0 c and 2 c.
                                                          !!          if ipcp=3, rainfall and snowfall are partitioned according to
                                                          !!          a polynomial curve between 0 c and 6 c. 
    IWF               => c_switch%IWF,                  & !< integer: if iwf=0, only overland flow and baseflow are modelled, and
                                                          !!          the ground surface slope is not modelled.
                                                          !!          if iwf=n (0<n<4), the watflood calculations of overland flow
                                                          !!          and interflow are performed; interflow is drawn from the top
                                                          !!          n soil layers. 
    ITC               => c_switch%ITC,                  & !< integer: itc, itcg and itg are switches to choose the iteration scheme to
                                                          !!          be used in calculating the canopy or ground surface temperature
                                                          !!          respectively.  if the switch is set to 1, a bisection method is
                                                          !!          used; if to 2, the newton-raphson method is used. 
    ITCG              => c_switch%ITCG,                 & !< integer: itc, itcg and itg are switches to choose the iteration scheme to
                                                          !!          be used in calculating the canopy or ground surface temperature
                                                          !!          respectively.  if the switch is set to 1, a bisection method is
                                                          !!          used; if to 2, the newton-raphson method is used. 
    ITG               => c_switch%ITG,                  & !< integer: itc, itcg and itg are switches to choose the iteration scheme to
                                                          !!          be used in calculating the canopy or ground surface temperature
                                                          !!          respectively.  if the switch is set to 1, a bisection method is
                                                          !!          used; if to 2, the newton-raphson method is used. 
    IPAI              => c_switch%IPAI,                 & !< integer: if ipai, ihgt, ialc, ials and ialg are zero, the values of
                                                          !!          plant area index, vegetation height, canopy albedo, snow albedo
                                                          !!          and soil albedo respectively calculated by class are used.
                                                          !!          if any of these switches is set to 1, the value of the
                                                          !!          corresponding parameter calculated by class is overridden by
                                                          !!          a user-supplied input value. 
    IHGT              => c_switch%IHGT,                 & !< integer: if ipai, ihgt, ialc, ials and ialg are zero, the values of
                                                          !!          plant area index, vegetation height, canopy albedo, snow albedo
                                                          !!          and soil albedo respectively calculated by class are used.
                                                          !!          if any of these switches is set to 1, the value of the
                                                          !!          corresponding parameter calculated by class is overridden by
                                                          !!          a user-supplied input value. 
    IALC              => c_switch%IALC,                 & !< integer: if ipai, ihgt, ialc, ials and ialg are zero, the values of
                                                          !!          plant area index, vegetation height, canopy albedo, snow albedo
                                                          !!          and soil albedo respectively calculated by class are used.
                                                          !!          if any of these switches is set to 1, the value of the
                                                          !!          corresponding parameter calculated by class is overridden by
                                                          !!          a user-supplied input value. 
    IALS              => c_switch%IALS,                 & !< integer: if ipai, ihgt, ialc, ials and ialg are zero, the values of
                                                          !!          plant area index, vegetation height, canopy albedo, snow albedo
                                                          !!          and soil albedo respectively calculated by class are used.
                                                          !!          if any of these switches is set to 1, the value of the
                                                          !!          corresponding parameter calculated by class is overridden by
                                                          !!          a user-supplied input value. 
    IALG              => c_switch%IALG,                 & !< integer: if ipai, ihgt, ialc, ials and ialg are zero, the values of
                                                          !!          plant area index, vegetation height, canopy albedo, snow albedo
                                                          !!          and soil albedo respectively calculated by class are used.
                                                          !!          if any of these switches is set to 1, the value of the
                                                          !!          corresponding parameter calculated by class is overridden by
                                                          !!          a user-supplied input value. 
    isnoalb           => c_switch%isnoalb,              & !< integer: if isnoalb is set to 0, the original two-band snow albedo algorithms are used.
                                                          !!          if it is set to 1, the new four-band routines are used.
    blackCdepon       => c_switch%blackCdepon,          & !< logical: If true, include consideration of black carbon deposition in snow processes (4-band scheme)  
    KsatScalaron      => c_switch%KsatScalaron,         & !< logical: If true, include exponential distribution scalar for saturated hydraulic conductivity (GRKSAT)
    todfrac           => vgat%todfrac,                  & !< real, dimension(:,:)  :!(ilg,icc) 
    netrad_gat        => vgat%netrad_gat,               & !< real, dimension(:)    :!(ilg) 
    preacc_gat        => vgat%preacc_gat,               & !< real, dimension(:)    :!(ilg) 
    sdepgat           => vgat%sdepgat,                  & !< real, dimension(:)    :!(ilg) 
    sandgat           => vgat%sandgat,                  & !< real, dimension(:,:)  :!(ilg,ignd) 
    claygat           => vgat%claygat,                  & !< real, dimension(:,:)  :!(ilg,ignd) 
    orgmgat           => vgat%orgmgat,                  & !< real, dimension(:,:)  :!(ilg,ignd) 
    xdiffusgat        => vgat%xdiffusgat,               & !< real, dimension(:)    :!(ilg) ! the corresponding ROW is CLASS's XDIFFUS 
    faregat           => vgat%faregat,                  & !< real, dimension(:)    :!(ilg) ! the ROT is FAREROT 

    ! Model switches:
    !agcm_classic_on   => c_switch%agcm_classic_on,      & !< logical:  Not currently used, but might be needed in the future.
    ctem_on           => c_switch%ctem_on,              & !< logical: 
    Ncycle_on         => c_switch%Ncycle_on,            & !< logical:
    dofire            => c_switch%dofire,               & !< logical: 
    PFTCompetition    => c_switch%PFTCompetition,       & !< logical: 
    lnduseon          => c_switch%lnduseon,             & !< logical: 
    doMethane         => c_switch%doMethane,            & !< logical: 
    inibioclim        => c_switch%inibioclim,           & !< logical: 
    dynamicTilingOn    => c_switch%dynamicTilingOn,     & !< logical:
    prescribedFire    => c_switch%prescribedFire,       & !< logical:
    timberHarvest     => c_switch%timberHarvest,        & !< logical
    trackTileAge      => c_switch%trackTileAge,         & !< logical

    ! ROW vars:

    pftexistrow       => vrot%pftexist,                 & !< logical, dimension(:,:,:) : 
    ccrow             => vrot%cc,                 & !< real, dimension(:,:,:) : 
    mmrow             => vrot%mm,                 & !< real, dimension(:,:,:) : 
    colddays_leaffallrow => vrot%colddays_leaffall,     & !< integer, dimension(:,:) :
    colddays_harvestrow  => vrot%colddays_harvest,      & !< integer, dimension(:,:) : 
    lfstatusrow       => vrot%lfstatus,                 & !< integer, dimension(:,:,:) : 
    pandaysrow        => vrot%pandays,                  & !< integer, dimension(:,:,:) : 

    gleafmasrow       => vrot%gleafmas,                 & !< real, dimension(:,:,:) :! 
    gleafmas_NSrow    => vrot%gleafmas_ns,              & !< real, dimension(:,:,:) :
    gleafmassrow      => vrot%gleafmas_s,               & !< real, dimension(:,:,:) : 
    stemmass_NSrow    => vrot%stemmass_ns,              & !< real, dimension(:,:,:) :
    stemmasssrow      => vrot%stemmass_s,               & !< real, dimension(:,:,:) :
    rootmass_NSrow    => vrot%rootmass_ns,              & !< real, dimension(:,:,:) :
    rootmasssrow      => vrot%rootmass_s,               & !< real, dimension(:,:,:) :
    bleafmasrow       => vrot%bleafmas,                 & !< real, dimension(:,:,:) :! 
    stemmassrow       => vrot%stemmass,                 & !< real, dimension(:,:,:) :! 
    rootmassrow       => vrot%rootmass,                 & !< real, dimension(:,:,:) :! 
    pstemmassrow      => vrot%pstemmass,                & !< real, dimension(:,:,:) :! 
    pgleafmassrow     => vrot%pgleafmass,               & !< real, dimension(:,:,:) :! 
    fcancmxrow        => vrot%fcancmx,                  & !< real, dimension(:,:,:) : 
    gavglairow        => vrot%gavglai,                  & !< real, dimension(:,:) : 
    zolncrow          => vrot%zolnc,                    & !< real, dimension(:,:,:) : 
    ailcrow           => vrot%ailc,                     & !< real, dimension(:,:,:) : 
    ailcgrow          => vrot%ailcg,                    & !< real, dimension(:,:,:) : 
    ailcgsrow         => vrot%ailcgs,                   & !< real, dimension(:,:,:) : 
    fcancsrow         => vrot%fcancs,                   & !< real, dimension(:,:,:) : 
    fcancrow          => vrot%fcanc,                    & !< real, dimension(:,:,:) : 
    co2concrow        => vrot%co2conc,                  & !< real, dimension(:,:) : 
    ch4concrow        => vrot%ch4conc,                  & !< real, dimension(:,:) : 
    co2i1cgrow        => vrot%co2i1cg,                  & !< real, dimension(:,:,:) : 
    co2i1csrow        => vrot%co2i1cs,                  & !< real, dimension(:,:,:) : 
    co2i2cgrow        => vrot%co2i2cg,                  & !< real, dimension(:,:,:) : 
    co2i2csrow        => vrot%co2i2cs,                  & !< real, dimension(:,:,:) : 
    ancsvegrow        => vrot%ancsveg,                  & !< real, dimension(:,:,:) : 
    ancgvegrow        => vrot%ancgveg,                  & !< real, dimension(:,:,:) : 
    rmlcsvegrow       => vrot%rmlcsveg,                 & !< real, dimension(:,:,:) : 
    rmlcgvegrow       => vrot%rmlcgveg,                 & !< real, dimension(:,:,:) : 
    slairow           => vrot%slai,                     & !< real, dimension(:,:,:) : 
    ailcbrow          => vrot%ailcb,                    & !< real, dimension(:,:,:) : 
    canresrow         => vrot%canres,                   & !< real, dimension(:,:) : 
    flhrlossrow       => vrot%flhrloss,                 & !< real, dimension(:,:,:) : 
    flhrloss_nsrow    => vrot%flhrloss_ns,              & !< real, dimension(:,:,:) : 
    flhrloss_srow     => vrot%flhrloss_s,               & !< real, dimension(:,:,:) : 
    grwtheffrow       => vrot%grwtheff,                 & !< real, dimension(:,:,:) : 
    lystmmasrow       => vrot%lystmmas,                 & !< real, dimension(:,:,:) : 
    lyrotmasrow       => vrot%lyrotmas,                 & !< real, dimension(:,:,:) : 
    lmaxtrow          => vrot%lmaxt,                    & !< real, dimension(:,:,:) :
    smaxtrow          => vrot%smaxt,                    & !< real, dimension(:,:,:) :
    rmaxtrow          => vrot%rmaxt,                    & !< real, dimension(:,:,:) :
    lygleafmasmaxrow  => vrot%lygleafmasmax,            & !< real, dimension(:,:,:) :
    lystemmassmaxrow  => vrot%lystemmassmax,            & !< real, dimension(:,:,:) :
    lyrootmassmaxrow  => vrot%lyrootmassmax,            & !< real, dimension(:,:,:) :
    tymaxlairow       => vrot%tymaxlai,                 & !< real, dimension(:,:,:) : 
    vgbiomasrow       => vrot%vgbiomas,                 & !< real, dimension(:,:) : 
    gavgltmsrow       => vrot%gavgltms,                 & !< real, dimension(:,:) : 
    gavgscmsrow       => vrot%gavgscms,                 & !< real, dimension(:,:) : 
    stmhrlosrow       => vrot%stmhrlos,                 & !< real, dimension(:,:,:) : 
    rmatcrow          => vrot%rmatc,                    & !< real, dimension(:,:,:,:) : 
    rmatctemrow       => vrot%rmatctem,                 & !< real, dimension(:,:,:,:) : 
    litrmassrow       => vrot%litrmass,                 & !< real, dimension(:,:,:,:) : 
    soilcmasrow       => vrot%soilcmas,                 & !< real, dimension(:,:,:,:) : 
    litresvegrow      => vrot%litresveg,                & !< real, dimension(:,:,:,:) : 
    soilcresvegrow    => vrot%soilcresveg,              & !< real, dimension(:,:,:,:) : 
    humiftrsvegrow    => vrot%humiftrsveg,              & !< real, dimension(:,:,:,:) : 
    nh4_massrow       => vrot%nh4_mass,                 & !< real, dimension(:,:,:) :
    no3_massrow       => vrot%no3_mass,                 & !< real, dimension(:,:,:) :
    ngleafmasrow      => vrot%ngleafmas,                & !< real, dimension(:,:,:) :
    ngleafmas_NSrow   => vrot%ngleafmas_ns,             & !< real, dimension(:,:,:) :
    ngleafmassrow     => vrot%ngleafmas_s,              & !< real, dimension(:,:,:) :
    nbleafmasrow      => vrot%nbleafmas,                & !< real, dimension(:,:,:) :
    nstemmassrow      => vrot%nstemmass,                & !< real, dimension(:,:,:) :
    nstemmass_NSrow   => vrot%nstemmass_ns,             & !< real, dimension(:,:,:) :
    nstemmasssrow     => vrot%nstemmass_s,              & !< real, dimension(:,:,:) :
    nrootmassrow      => vrot%nrootmass,                & !< real, dimension(:,:,:) :
    nrootmass_NSrow   => vrot%nrootmass_ns,             & !< real, dimension(:,:,:) :
    nrootmasssrow     => vrot%nrootmass_s,              & !< real, dimension(:,:,:) :
    nlitrmassrow      => vrot%nlitrmass,                & !< real, dimension(:,:,:) :
    soilnmasrow       => vrot%soilnmas,                 & !< real, dimension(:,:,:) :
    nvgbiomas_vegrow  => vrot%nvgbiomas_veg,            & !< real, dimension(:,:,:) :
    soilpHrow         => vrot%soilpH,                   & !< real, dimension(:,:) :
    nfertilrow        => vrot%nfertil,                  & !< real, dimension(:,:) :
    ndepositrow       => vrot%ndeposit,                 & !< real, dimension(:,:) :
    timharvrow        => vrot%timharvrow,                & !< real, dimension(:) :
    vgbiomas_vegrow   => vrot%vgbiomas_veg,             & !< real, dimension(:,:,:) : 

    emit_co2row       => vrot%emit_co2,                 & !< real, dimension(:,:,:) : 
    emit_corow        => vrot%emit_co,                  & !< real, dimension(:,:,:) : 
    emit_ch4row       => vrot%emit_ch4,                 & !< real, dimension(:,:,:) : 
    emit_nmhcrow      => vrot%emit_nmhc,                & !< real, dimension(:,:,:) : 
    emit_h2row        => vrot%emit_h2,                  & !< real, dimension(:,:,:) : 
    emit_noxrow       => vrot%emit_nox,                 & !< real, dimension(:,:,:) : 
    emit_n2orow       => vrot%emit_n2o,                 & !< real, dimension(:,:,:) : 
    emit_nh3row       => vrot%emit_nh3,                 & !< real, dimension(:,:,:) : 
    emit_pm25row      => vrot%emit_pm25,                & !< real, dimension(:,:,:) : 
    emit_tpmrow       => vrot%emit_tpm,                 & !< real, dimension(:,:,:) : 
    emit_tcrow        => vrot%emit_tc,                  & !< real, dimension(:,:,:) : 
    emit_ocrow        => vrot%emit_oc,                  & !< real, dimension(:,:,:) : 
    emit_bcrow        => vrot%emit_bc,                  & !< real, dimension(:,:,:) : 
    burnfracrow       => vrot%burnfrac,                 & !< real, dimension(:,:) : 
    burnvegfrow       => vrot%burnvegf,                 & !< real, dimension(:,:,:) : 
    smfuncvegrow      => vrot%smfuncveg,                & !< real, dimension(:,:,:) : 
    popdinrow         => vrot%popdin,                   & !< real, dimension(:,:) : 
    btermrow          => vrot%bterm,                    & !< real, dimension(:,:,:) : 
    ltermrow          => vrot%lterm,                    & !< real, dimension(:,:) : 
    mtermrow          => vrot%mterm,                    & !< real, dimension(:,:,:) : 

    extnprobrow       => vrot%extnprob,                 & !< real, dimension(:,:) : 
    prbfrhucrow       => vrot%prbfrhuc,                 & !< real, dimension(:,:) : 
    daylrow           => vrot%dayl,                     & !< real, dimension(:) : 
    dayl_maxrow       => vrot%dayl_max,                 & !< real, dimension(:) : 
    grclarearow       => vrot%grclarearow,              & !< real, dimension(:) : 

    bmasvegrow        => vrot%bmasveg,                  & !< real, dimension(:,:,:) : 
    cmasvegcrow       => vrot%cmasvegc,                 & !< real, dimension(:,:,:) : 
    veghghtrow        => vrot%veghght,                  & !< real, dimension(:,:,:) : 
    rootdpthrow       => vrot%rootdpth,                 & !< real, dimension(:,:,:) : 
    rmlrow            => vrot%rml,                      & !< real, dimension(:,:) : 
    rmsrow            => vrot%rms,                      & !< real, dimension(:,:) : 
    tltrleafrow       => vrot%tltrleaf,                 & !< real, dimension(:,:,:) : 
    tltrstemrow       => vrot%tltrstem,                 & !< real, dimension(:,:,:) : 
    tltrrootrow       => vrot%tltrroot,                 & !< real, dimension(:,:,:) : 
    leaflitrrow       => vrot%leaflitr,                 & !< real, dimension(:,:,:) : 
    roottemprow       => vrot%roottemp,                 & !< real, dimension(:,:,:) : 
    afrleafrow        => vrot%afrleaf,                  & !< real, dimension(:,:,:) : 
    afrstemrow        => vrot%afrstem,                  & !< real, dimension(:,:,:) : 
    afrrootrow        => vrot%afrroot,                  & !< real, dimension(:,:,:) : 
    wtstatusrow       => vrot%wtstatus,                 & !< real, dimension(:,:,:) : 
    ltstatusrow       => vrot%ltstatus,                 & !< real, dimension(:,:,:) : 
    rmrrow            => vrot%rmr,                      & !< real, dimension(:,:) : 

    slopefracrow      => vrot%slopefrac,                & !< real, dimension(:,:,:) : 
    wetfrac_presrow   => vrot%wetfrac_pres,             & !< real, dimension(:,:) : 
    ch4WetSpecrow     => vrot%ch4WetSpec,               & !< real, dimension(:,:) : 
    wetfdynrow        => vrot%wetfdyn,                  & !< real, dimension(:,:) : 
    ch4WetDynrow      => vrot%ch4WetDyn,                & !< real, dimension(:,:) : 
    ch4soillsrow      => vrot%ch4_soills,               & !< real, dimension(:,:) : 

    peatdeprow        => vrot%peatdep,                  & !< real, dimension(:,:) : 
    lucemcomrow       => vrot%lucemcom,                 & !< real, dimension(:,:) : 
    lucltrinrow       => vrot%lucltrin,                 & !< real, dimension(:,:) : 
    lucsocinrow       => vrot%lucsocin,                 & !< real, dimension(:,:) : 
    lucemcomnrow      => vrot%lucemcomn,                & !< real, dimension(:,:) : 
    lucltrinnrow      => vrot%lucltrinn,                & !< real, dimension(:,:) : 
    lucsocinnrow      => vrot%lucsocinn,                & !< real, dimension(:,:) : 

    npprow            => vrot%npp,                      & !< real, dimension(:,:) : 
    neprow            => vrot%nep,                      & !< real, dimension(:,:) : 
    nepCMIProw        => vrot%nepCMIP,                      & !< real, dimension(:,:) :
    nbprow            => vrot%nbp,                      & !< real, dimension(:,:) : 
    gpprow            => vrot%gpp,                      & !< real, dimension(:,:) : 
    hetroresrow       => vrot%hetrores,                 & !< real, dimension(:,:) : 
    autoresrow        => vrot%autores,                  & !< real, dimension(:,:) : 
    soilcresprow      => vrot%soilcresp,                & !< real, dimension(:,:) : 
    rmrow             => vrot%rm,                       & !< real, dimension(:,:) : 
    rgrow             => vrot%rg,                       & !< real, dimension(:,:) : 
    litresrow         => vrot%litres,                   & !< real, dimension(:,:) : 
    socresrow         => vrot%socres,                   & !< real, dimension(:,:) : 
    dstcemlsrow       => vrot%dstcemls,                 & !< real, dimension(:,:) : 
    litrfallrow       => vrot%litrfall,                 & !< real, dimension(:,:) : 
    humiftrsrow       => vrot%humiftrs,                 & !< real, dimension(:,:) : 

    gppvegrow         => vrot%gppveg,                   & !< real, dimension(:,:,:) : 
    leafns2srow       => vrot%leafns2s,                 & !< real, dimension(:,:,:) :
    stemns2srow       => vrot%stemns2s,                 & !< real, dimension(:,:,:) :
    rootns2srow       => vrot%rootns2s,                 & !< real, dimension(:,:,:) :
    re_alloc_s2lrow   => vrot%re_alloc_s2l,             & !< real, dimension(:,:,:) :
    re_alloc_r2lrow   => vrot%re_alloc_r2l,             & !< real, dimension(:,:,:) :
    re_alloc_sr2lrow  => vrot%re_alloc_sr2l,            & !< real, dimension(:,:,:) :
    nepvegrow         => vrot%nepveg,                   & !< real, dimension(:,:,:) : 
    nbpvegrow         => vrot%nbpveg,                   & !< real, dimension(:,:,:) : 
    nppvegrow         => vrot%nppveg,                   & !< real, dimension(:,:,:) : 
    hetroresvegrow    => vrot%hetroresveg,              & !< real, dimension(:,:,:) : 
    autoresvegrow     => vrot%autoresveg,               & !< real, dimension(:,:,:) : 
    rmlvegaccrow      => vrot%rmlvegacc,                & !< real, dimension(:,:,:) : 
    rmsvegrow         => vrot%rmsveg,                   & !< real, dimension(:,:,:) : 
    rmrvegrow         => vrot%rmrveg,                   & !< real, dimension(:,:,:) : 
    rgvegrow          => vrot%rgveg,                    & !< real, dimension(:,:,:) : 
    litrfallvegrow    => vrot%litrfallveg,              & !< real, dimension(:,:,:) : 
    rothrlosrow       => vrot%rothrlos,                 & !< real, dimension(:,:,:) : 
    pfcancmxrow       => vrot%pfcancmx,                 & !< real, dimension(:,:,:) : 
    nfcancmxrow       => vrot%nfcancmx,                 & !< real, dimension(:,:,:) : 
    alvsctmrow        => vrot%alvsctm,                  & !< real, dimension(:,:,:) : 
    paicrow           => vrot%paic,                     & !< real, dimension(:,:,:) : 
    slaicrow          => vrot%slaic,                    & !< real, dimension(:,:,:) : 
    alirctmrow        => vrot%alirctm,                  & !< real, dimension(:,:,:) : 
    cfluxcgrow        => vrot%cfluxcg,                  & !< real, dimension(:,:) : 
    cfluxcsrow        => vrot%cfluxcs,                  & !< real, dimension(:,:) : 
    vcmax0row         => vrot%vcmax0,                   & !< real, dimension(:,:,:) :
    dstcemls3row      => vrot%dstcemls3,                & !< real, dimension(:,:) : 
    anvegrow          => vrot%anveg,                    & !< real, dimension(:,:,:) : 
    rmlvegrow         => vrot%rmlveg,                   & !< real, dimension(:,:,:) : 

    twarmmrow            => vrot%twarmm,                & !< real, dimension(:,:) : 
    tcoldmrow            => vrot%tcoldm,                & !< real, dimension(:,:) : 
    gdd5row              => vrot%gdd5,                  & !< real, dimension(:,:) : 
    aridityrow           => vrot%aridity,               & !< real, dimension(:,:) : 
    srplsmonrow          => vrot%srplsmon,              & !< real, dimension(:,:) : 
    defctmonrow          => vrot%defctmon,              & !< real, dimension(:,:) : 
    anndefctrow          => vrot%anndefct,              & !< real, dimension(:,:) : 
    annsrplsrow          => vrot%annsrpls,              & !< real, dimension(:,:) : 
    annpcprow            => vrot%annpcp,                & !< real, dimension(:,:) : 
    dry_season_lengthrow => vrot%dry_season_length,     & !< real, dimension(:,:) : 

    ipeatlandrow     => vrot%ipeatland,                 & !< integer, dimension(:,:) :! This is first set in read_from_ctm. 
    anmossrow        => vrot%anmoss,                    & !< real, dimension(:,:) : 
    rmlmossrow       => vrot%rmlmoss,                   & !< real, dimension(:,:) : 
    gppmossrow       => vrot%gppmoss,                   & !< real, dimension(:,:) : 
    nppmossrow       => vrot%nppmoss,                   & !< real, dimension(:,:) : 
    armossrow        => vrot%armoss,                    & !< real, dimension(:,:) : 
    litrmsmossrow    => vrot%litrmsmoss,                & !< real, dimension(:,:) : 
    Cmossmasrow      => vrot%Cmossmas,                  & !< real, dimension(:,:) : 
    dmossrow         => vrot%dmoss,                     & !< real, dimension(:,:) : 
    peatSoilCrow     => vrot%peatSoilC,                 & !< real, dimension(:,:) : 
    pddrow           => vrot%pdd,                       & !< real, dimension(:,:) : 

    bnffreerow          => vrot%bnf_free,               & !< real, dimension(:,:,:) :
    bnfnatrow           => vrot%bnf_nat,                & !< real, dimension(:,:,:) :
    bnfantrow           => vrot%bnf_ant,                & !< real, dimension(:,:,:) :
    bnftotrow           => vrot%bnf_tot,                & !< real, dimension(:,:,:) :
    nstressrow          => vrot%nstress,                & !< real, dimension(:,:,:) :
    nitrifvegrow        => vrot%nitrifveg,              & !< real, dimension(:,:,:) :
    no_nitvegrow        => vrot%no_nitveg,              & !< real, dimension(:,:,:) :
    no_denitvegrow      => vrot%no_denitveg,            & !< real, dimension(:,:,:) :
    no_nitdenitvegrow   => vrot%no_nitdenitveg,         & !< real, dimension(:,:,:) :
    n2o_nitvegrow       => vrot%n2o_nitveg,             & !< real, dimension(:,:,:) :
    n2o_denitvegrow     => vrot%n2o_denitveg,           & !< real, dimension(:,:,:) :
    n2o_nitdenitvegrow  => vrot%n2o_nitdenitveg,        & !< real, dimension(:,:,:) :
    n2_denitvegrow      => vrot%n2_denitveg,            & !< real, dimension(:,:,:) :
    nvolvegrow          => vrot%nvolveg,                & !< real, dimension(:,:,:) :
    nleachvegrow        => vrot%nleachveg,              & !< real, dimension(:,:,:) :
    appl_fertrow        => vrot%appl_fert,              & !< real, dimension(:,:,:) :
    ndep_nh4row         => vrot%ndep_nh4,               & !< real, dimension(:,:,:) :
    ndep_no3row         => vrot%ndep_no3,               & !< real, dimension(:,:,:) :
    ndemandveg_wp_npprow=> vrot%ndemandveg_wp_npp,      & !< real, dimension(:,:,:) :
    nuptakeveg_p_nh4row => vrot%nuptakeveg_p_nh4,       & !< real, dimension(:,:,:) :
    nuptakeveg_p_no3row => vrot%nuptakeveg_p_no3,       & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_nh4row => vrot%nuptakeveg_a_actl_nh4, & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_no3row => vrot%nuptakeveg_a_actl_no3, & !< real, dimension(:,:,:) :
    nuptakevegrow      => vrot%nuptakeveg,              & !< real, dimension(:,:,:) :
    nallocveg_lrow     => vrot%nallocveg_l,             & !< real, dimension(:,:,:) :
    nallocveg_srow     => vrot%nallocveg_s,             & !< real, dimension(:,:,:) :
    nallocveg_rrow     => vrot%nallocveg_r,             & !< real, dimension(:,:,:) :
    nresorpedveg_srow  => vrot%nresorpedveg_s,          & !< real, dimension(:,:,:) :
    nresorpedveg_rrow  => vrot%nresorpedveg_r,          & !< real, dimension(:,:,:) :
    nre_allocveg_s2lrow=> vrot%nre_allocveg_s2l,        & !< real, dimension(:,:,:) :
    nre_allocveg_r2lrow=> vrot%nre_allocveg_r2l,        & !< real, dimension(:,:,:) :
    nleafns2svegrow    => vrot%nleafns2sveg,            & !< real, dimension(:,:,:) :
    nstemns2svegrow    => vrot%nstemns2sveg,            & !< real, dimension(:,:,:) :
    nrootns2svegrow    => vrot%nrootns2sveg,            & !< real, dimension(:,:,:) :
    nlitrveg_lrow      => vrot%nlitrveg_l,              & !< real, dimension(:,:,:) :
    nlitrveg_srow      => vrot%nlitrveg_s,              & !< real, dimension(:,:,:) :
    nlitrveg_rrow      => vrot%nlitrveg_r,              & !< real, dimension(:,:,:) :
    nlitrvegrow        => vrot%nlitrveg,                & !< real, dimension(:,:,:) :
    gl2bl_grass_nfluxrow => vrot%gl2bl_grass_nflux,     & !< real, dimension(:,:,:) :
    c2nveg_lrow        => vrot%c2nveg_l,                & !< real, dimension(:,:,:) :
    c2nveg_srow        => vrot%c2nveg_s,                & !< real, dimension(:,:,:) :
    c2nveg_rrow        => vrot%c2nveg_r,                & !< real, dimension(:,:,:) :
    c2nveg_wprow       => vrot%c2nveg_wp,               & !< real, dimension(:,:,:) :
    c2nveg_litrrow     => vrot%c2nveg_litr,             & !< real, dimension(:,:,:) :
    c2nveg_humusrow    => vrot%c2nveg_humus,            & !< real, dimension(:,:,:) :

    nhumtrsvegrow      => vrot%nhumtrsveg,              & !< real, dimension(:,:,:) :
    nmineralveg_litrrow=> vrot%nmineralveg_litr,        & !< real, dimension(:,:,:) :
    nmineralveg_humusrow=> vrot%nmineralveg_humus,      & !< real, dimension(:,:,:) :
    netnmineralveg_row=> vrot%netnmineralveg,           & !< real, dimension(:,:,:) :
    nimmobilveg_nh4row => vrot%nimmobilveg_nh4,         & !< real, dimension(:,:,:) :
    nimmobilveg_no3row => vrot%nimmobilveg_no3,         & !< real, dimension(:,:,:) :
    fNnetlandvegrow    => vrot%fNnetlandveg,            & !< real, dimension(:,:,:) :
    redcoeff_vcmaxrow  => vrot%redcoeff_vcmax,          & !< real, dimension(:,:,:) :


    controlVector     => vrot%controlVector,            & !< integer, dimension(nlat,nmos)
    DynTilInitializeFlag         => vrot%DynTilInitializeFlag,                & !< logical
    timharvarearow    => vrot%timharvarearow,           & !< real(nlat,nmos)
    prsfirearearow    => vrot%prsfirearearow,           & !< real(nlat,nmos)
    prsfirerow    => vrot%prsfirerow,           & !< real(nlat,nmos)
    tileAgerow        => vrot%tileAgerow,               & !< real(nlat,nmos)

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>

    ! GAT:

    pftexistgat       => vgat%pftexist,                 & !< logical, dimension(:,:) : 
    colddays_leaffallgat => vgat%colddays_leaffall,     & !< integer, dimension(:) :
    colddays_harvestgat  => vgat%colddays_harvest,      & !< integer, dimension(:) : 
    lfstatusgat       => vgat%lfstatus,                 & !< integer, dimension(:,:) : 
    pandaysgat        => vgat%pandays,                  & !< integer, dimension(:,:) : 
    lightng           => vgat%lightng,                  & !< real, dimension(:) : 

    gleafmasgat       => vgat%gleafmas,                 & !< real, dimension(:,:) :! 
    gleafmasgat_ns    => vgat%gleafmas_ns,              & !< real, dimension(:,:) :
    gleafmasgat_s     => vgat%gleafmas_s,               & !< real, dimension(:,:) :
    bleafmasgat       => vgat%bleafmas,                 & !< real, dimension(:,:) :! 
    stemmassgat       => vgat%stemmass,                 & !< real, dimension(:,:) :! 
    stemmassgat_ns    => vgat%stemmass_ns,              & !< real, dimension(:,:) :
    stemmassgat_s     => vgat%stemmass_s,               & !< real, dimension(:,:) :
    rootmassgat       => vgat%rootmass,                 & !< real, dimension(:,:) :! 
    rootmassgat_ns    => vgat%rootmass_ns,              & !< real, dimension(:,:) :! 
    rootmassgat_s     => vgat%rootmass_s,               & !< real, dimension(:,:) :! 
    pstemmassgat      => vgat%pstemmass,                & !< real, dimension(:,:) :! 
    pgleafmassgat     => vgat%pgleafmass,               & !< real, dimension(:,:) :! 
    fcancmxgat        => vgat%fcancmx,                  & !< real, dimension(:,:) : 
    gavglaigat        => vgat%gavglai,                  & !< real, dimension(:) : 
    zolncgat          => vgat%zolnc,                    & !< real, dimension(:,:) : 
    ailcgat           => vgat%ailc,                     & !< real, dimension(:,:) : 
    ailcggat          => vgat%ailcg,                    & !< real, dimension(:,:) : 
    ailcgsgat         => vgat%ailcgs,                   & !< real, dimension(:,:) : 
    fcancsgat         => vgat%fcancs,                   & !< real, dimension(:,:) : 
    fcancgat          => vgat%fcanc,                    & !< real, dimension(:,:) : 
    co2concgat        => vgat%co2conc,                  & !< real, dimension(:) : 
    ch4concgat        => vgat%ch4conc,                  & !< real, dimension(:) : 
    co2i1cggat        => vgat%co2i1cg,                  & !< real, dimension(:,:) : 
    co2i1csgat        => vgat%co2i1cs,                  & !< real, dimension(:,:) : 
    co2i2cggat        => vgat%co2i2cg,                  & !< real, dimension(:,:) : 
    co2i2csgat        => vgat%co2i2cs,                  & !< real, dimension(:,:) : 
    ancsveggat        => vgat%ancsveg,                  & !< real, dimension(:,:) : 
    ancgveggat        => vgat%ancgveg,                  & !< real, dimension(:,:) : 
    rmlcsveggat       => vgat%rmlcsveg,                 & !< real, dimension(:,:) : 
    rmlcgveggat       => vgat%rmlcgveg,                 & !< real, dimension(:,:) : 
    slaigat           => vgat%slai,                     & !< real, dimension(:,:) : 
    ailcbgat          => vgat%ailcb,                    & !< real, dimension(:,:) : 
    canresgat         => vgat%canres,                   & !< real, dimension(:) : 
    flhrlossgat       => vgat%flhrloss,                 & !< real, dimension(:,:) : 
    flhrloss_nsgat    => vgat%flhrloss_ns,              & !< real, dimension(:,:) : 
    flhrloss_sgat     => vgat%flhrloss_s,               & !< real, dimension(:,:) : 
    grwtheffgat       => vgat%grwtheff,                 & !< real, dimension(:,:) : 
    lystmmasgat       => vgat%lystmmas,                 & !< real, dimension(:,:) : 
    lyrotmasgat       => vgat%lyrotmas,                 & !< real, dimension(:,:) : 
    lmaxtgat          => vgat%lmaxt,                    & !< real, dimension(:,:) : 
    smaxtgat          => vgat%smaxt,                    & !< real, dimension(:,:) : 
    rmaxtgat          => vgat%rmaxt,                    & !< real, dimension(:,:) : 
    lygleafmasmaxgat  => vgat%lygleafmasmax,            & !< real, dimension(:,:) : 
    lystemmassmaxgat  => vgat%lystemmassmax,            & !< real, dimension(:,:) : 
    lyrootmassmaxgat  => vgat%lyrootmassmax,            & !< real, dimension(:,:) : 
    tymaxlaigat       => vgat%tymaxlai,                 & !< real, dimension(:,:) : 
    vgbiomasgat       => vgat%vgbiomas,                 & !< real, dimension(:) : 
    gavgltmsgat       => vgat%gavgltms,                 & !< real, dimension(:) : 
    gavgscmsgat       => vgat%gavgscms,                 & !< real, dimension(:) : 
    stmhrlosgat       => vgat%stmhrlos,                 & !< real, dimension(:,:) : 
    rmatcgat          => vgat%rmatc,                    & !< real, dimension(:,:,:) : 
    rmatctemgat       => vgat%rmatctem,                 & !< real, dimension(:,:,:) : 
    litrmassgat       => vgat%litrmass,                 & !< real, dimension(:,:,:) : 
    soilcmasgat       => vgat%soilcmas,                 & !< real, dimension(:,:,:) : 
    litresveggat      => vgat%litresveg,                & !< real, dimension(:,:,:) : 
    soilcresveggat    => vgat%soilcresveg,              & !< real, dimension(:,:,:) : 
    humiftrsveggat    => vgat%humiftrsveg,              & !< real, dimension(:,:,:) : 
    nh4_massgat       => vgat%nh4_mass,                 & !< real, dimension(:,:) : 
    no3_massgat       => vgat%no3_mass,                 & !< real, dimension(:,:) : 
    ngleafmasgat      => vgat%ngleafmas,                & !< real, dimension(:,:) : 
    ngleafmasgat_ns   => vgat%ngleafmas_ns,             & !< real, dimension(:,:) : 
    ngleafmasgat_s    => vgat%ngleafmas_s,              & !< real, dimension(:,:) : 
    nbleafmasgat      => vgat%nbleafmas,                & !< real, dimension(:,:) : 
    nstemmassgat      => vgat%nstemmass,                & !< real, dimension(:,:) : 
    nstemmassgat_ns   => vgat%nstemmass_ns,             & !< real, dimension(:,:) : 
    nstemmassgat_s    => vgat%nstemmass_s,              & !< real, dimension(:,:) : 
    nrootmassgat      => vgat%nrootmass,                & !< real, dimension(:,:) : 
    nrootmassgat_ns   => vgat%nrootmass_ns,             & !< real, dimension(:,:) : 
    nrootmassgat_s    => vgat%nrootmass_s,              & !< real, dimension(:,:) : 
    nlitrmassgat      => vgat%nlitrmass,                & !< real, dimension(:,:) : 
    soilnmasgat       => vgat%soilnmas,                 & !< real, dimension(:,:) : 
    nvgbiomas_veggat  => vgat%nvgbiomas_veg,            & !< real, dimension(:,:) : 
    vgbiomas_veggat   => vgat%vgbiomas_veg,             & !< real, dimension(:,:) : 

    emit_co2gat       => vgat%emit_co2,                 & !< real, dimension(:,:) : 
    emit_cogat        => vgat%emit_co,                  & !< real, dimension(:,:) : 
    emit_ch4gat       => vgat%emit_ch4,                 & !< real, dimension(:,:) : 
    emit_nmhcgat      => vgat%emit_nmhc,                & !< real, dimension(:,:) : 
    emit_h2gat        => vgat%emit_h2,                  & !< real, dimension(:,:) : 
    emit_noxgat       => vgat%emit_nox,                 & !< real, dimension(:,:) : 
    emit_n2ogat       => vgat%emit_n2o,                 & !< real, dimension(:,:) : 
    emit_nh3gat       => vgat%emit_nh3,                 & !< real, dimension(:,:) : 
    emit_pm25gat      => vgat%emit_pm25,                & !< real, dimension(:,:) : 
    emit_tpmgat       => vgat%emit_tpm,                 & !< real, dimension(:,:) : 
    emit_tcgat        => vgat%emit_tc,                  & !< real, dimension(:,:) : 
    emit_ocgat        => vgat%emit_oc,                  & !< real, dimension(:,:) : 
    emit_bcgat        => vgat%emit_bc,                  & !< real, dimension(:,:) : 
    burnfracgat       => vgat%burnfrac,                 & !< real, dimension(:) : 
    burnvegfgat       => vgat%burnvegf,                 & !< real, dimension(:,:) : 
    popdingat         => vgat%popdin,                   & !< real, dimension(:) : 
    soilPHgat         => vgat%soilPH,                   & !< real, dimension(:) : 
    nfertilgat        => vgat%nfertil,                  & !< real, dimension(:) : 
    ndepositgat       => vgat%ndeposit,                 & !< real, dimension(:) : 
    smfuncveggat      => vgat%smfuncveg,                & !< real, dimension(:,:) : 
    btermgat          => vgat%bterm,                    & !< real, dimension(:,:) : 
    ltermgat          => vgat%lterm,                    & !< real, dimension(:) : 
    mtermgat          => vgat%mterm,                    & !< real, dimension(:,:) : 
    glcaemls          => vgat%glcaemls,                 & !< real, dimension(:,:) : green leaf carbon emission disturbance losses, \f$kg c/m^2\f$ 
    blcaemls          => vgat%blcaemls,                 & !< real, dimension(:,:) : brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$ 
    rtcaemls          => vgat%rtcaemls,                 & !< real, dimension(:,:) : root carbon emission disturbance losses, \f$kg c/m^2\f$ 
    stcaemls          => vgat%stcaemls,                 & !< real, dimension(:,:) : stem carbon emission disturbance losses, \f$kg c/m^2\f$ 
    ltrcemls          => vgat%ltrcemls,                 & !< real, dimension(:,:) : litter carbon emission disturbance losses, \f$kg c/m^2\f$ 
    blfltrdt          => vgat%blfltrdt,                 & !< real, dimension(:,:) : brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$ 
    glfltrdt          => vgat%glfltrdt,                 & !< real, dimension(:,:) : green leaf litter generated due to disturbance \f$(kg c/m^2)\f$ 
    ntchlveg          => vgat%ntchlveg,                 & !< real, dimension(:,:) : fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec 
    ntchsveg          => vgat%ntchsveg,                 & !< real, dimension(:,:) : fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec 
    ntchrveg          => vgat%ntchrveg,                 & !< real, dimension(:,:) : fluxes for each pft: Net change in root biomass,
                                                          !!                        the net change is the difference between allocation and
                                                          !!                        autotrophic respiratory fluxes, u-mol CO2/m2.sec 
    extnprobgat       => vgat%extnprob,                 & !< real, dimension(:) : 
    prbfrhucgat       => vgat%prbfrhuc,                 & !< real, dimension(:) : 
    daylgat           => vgat%dayl,                     & !< real, dimension(:) : 
    dayl_maxgat       => vgat%dayl_max,                 & !< real, dimension(:) : 
    grclarea          => vgat%grclarea,                 & !< real, dimension(:) : 

    bmasveggat        => vgat%bmasveg,                  & !< real, dimension(:,:) : 
    cmasvegcgat       => vgat%cmasvegc,                 & !< real, dimension(:,:) : 
    veghghtgat        => vgat%veghght,                  & !< real, dimension(:,:) : 
    rootdpthgat       => vgat%rootdpth,                 & !< real, dimension(:,:) : 
    rmlgat            => vgat%rml,                      & !< real, dimension(:) : 
    rmsgat            => vgat%rms,                      & !< real, dimension(:) : 
    tltrleafgat       => vgat%tltrleaf,                 & !< real, dimension(:,:) : 
    tltrstemgat       => vgat%tltrstem,                 & !< real, dimension(:,:) : 
    tltrrootgat       => vgat%tltrroot,                 & !< real, dimension(:,:) : 
    leaflitrgat       => vgat%leaflitr,                 & !< real, dimension(:,:) : 
    roottempgat       => vgat%roottemp,                 & !< real, dimension(:,:) : 
    afrleafgat        => vgat%afrleaf,                  & !< real, dimension(:,:) : 
    afrstemgat        => vgat%afrstem,                  & !< real, dimension(:,:) : 
    afrrootgat        => vgat%afrroot,                  & !< real, dimension(:,:) : 
    wtstatusgat       => vgat%wtstatus,                 & !< real, dimension(:,:) : 
    ltstatusgat       => vgat%ltstatus,                 & !< real, dimension(:,:) : 
    rmrgat            => vgat%rmr,                      & !< real, dimension(:) : 

    slopefracgat      => vgat%slopefrac,                & !< real, dimension(:,:) : 
    wetfrac_presgat   => vgat%wetfrac_pres,             & !< real, dimension(:) : 
    ch4WetSpecgat     => vgat%ch4WetSpec,               & !< real, dimension(:) : 
    wetfdyngat        => vgat%wetfdyn,                  & !< real, dimension(:) : 
    ch4WetDyngat      => vgat%ch4WetDyn,                & !< real, dimension(:) : 
    ch4soillsgat      => vgat%ch4_soills,               & !< real, dimension(:) : 

    lucemcomgat       => vgat%lucemcom,                 & !< real, dimension(:) : 
    lucltringat       => vgat%lucltrin,                 & !< real, dimension(:) : 
    lucsocingat       => vgat%lucsocin,                 & !< real, dimension(:) : 
    lucemcomngat       => vgat%lucemcomn,               & !< real, dimension(:) : 
    lucltrinngat       => vgat%lucltrinn,               & !< real, dimension(:) : 
    lucsocinngat       => vgat%lucsocinn,               & !< real, dimension(:) : 

    nppgat            => vgat%npp,                      & !< real, dimension(:) : 
    nepgat            => vgat%nep,                      & !< real, dimension(:) : 
    nepCMIPgat        => vgat%nepCMIP,                  & !< real, dimension(:) :
    nbpgat            => vgat%nbp,                      & !< real, dimension(:) : 
    gppgat            => vgat%gpp,                      & !< real, dimension(:) : 
    hetroresgat       => vgat%hetrores,                 & !< real, dimension(:) : 
    autoresgat        => vgat%autores,                  & !< real, dimension(:) : 
    soilcrespgat      => vgat%soilcresp,                & !< real, dimension(:) : 
    rmgat             => vgat%rm,                       & !< real, dimension(:) : 
    rggat             => vgat%rg,                       & !< real, dimension(:) : 
    litresgat         => vgat%litres,                   & !< real, dimension(:) : 
    socresgat         => vgat%socres,                   & !< real, dimension(:) : 
    dstcemlsgat       => vgat%dstcemls,                 & !< real, dimension(:) : 
    litrfallgat       => vgat%litrfall,                 & !< real, dimension(:) : 
    humiftrsgat       => vgat%humiftrs,                 & !< real, dimension(:) : 

    gppveggat         => vgat%gppveg,                   & !< real, dimension(:,:) : 
    leafns2sgat       => vgat%leafns2s,                 & !< real, dimension(:,:) :
    stemns2sgat       => vgat%stemns2s,                 & !< real, dimension(:,:) :
    rootns2sgat       => vgat%rootns2s,                 & !< real, dimension(:,:) :
    re_alloc_s2lgat   => vgat%re_alloc_s2l,             & !< real, dimension(:,:) :
    re_alloc_r2lgat   => vgat%re_alloc_r2l,             & !< real, dimension(:,:) :
    re_alloc_sr2lgat  => vgat%re_alloc_sr2l,            & !< real, dimension(:,:) :
    nepveggat         => vgat%nepveg,                   & !< real, dimension(:,:) : 
    nbpveggat         => vgat%nbpveg,                   & !< real, dimension(:,:) : 
    nppveggat         => vgat%nppveg,                   & !< real, dimension(:,:) : 
    hetroresveggat    => vgat%hetroresveg,              & !< real, dimension(:,:) : 
    autoresveggat     => vgat%autoresveg,               & !< real, dimension(:,:) : 
    rmlvegaccgat      => vgat%rmlvegacc,                & !< real, dimension(:,:) : 
    rmsveggat         => vgat%rmsveg,                   & !< real, dimension(:,:) : 
    rmrveggat         => vgat%rmrveg,                   & !< real, dimension(:,:) : 
    rgveggat          => vgat%rgveg,                    & !< real, dimension(:,:) : 
    litrfallveggat    => vgat%litrfallveg,              & !< real, dimension(:,:) : 
    reprocost         => vgat%reprocost,                & !< real, dimension(:,:) : Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) 
    rothrlosgat       => vgat%rothrlos,                 & !< real, dimension(:,:) : 
    pfcancmxgat       => vgat%pfcancmx,                 & !< real, dimension(:,:) : 
    nfcancmxgat       => vgat%nfcancmx,                 & !< real, dimension(:,:) : 
    alvsctmgat        => vgat%alvsctm,                  & !< real, dimension(:,:) : 
    paicgat           => vgat%paic,                     & !< real, dimension(:,:) : 
    slaicgat          => vgat%slaic,                    & !< real, dimension(:,:) : 
    alirctmgat        => vgat%alirctm,                  & !< real, dimension(:,:) : 
    cfluxcggat        => vgat%cfluxcg,                  & !< real, dimension(:) : 
    cfluxcsgat        => vgat%cfluxcs,                  & !< real, dimension(:) : 
    CFLUX_GAgat       => vgat%CFLUX_GA,                 & !< real, dimension(:) : 
    vcmax0gat         => vgat%vcmax0,                   & !< real, dimension(:,:) : 
    USTARBS_GA        => vgat%USTARBS_GA,               & !< real, dimension(:) : Friction velocity to be used in nitrogen volatilization  \f$[m s^{-1} ]\f$
    ROFBGAT           => vgat%ROFB,                     & !< real, dimension(:) : Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    dstcemls3gat      => vgat%dstcemls3,                & !< real, dimension(:) : 
    anveggat          => vgat%anveg,                    & !< real, dimension(:,:) : 
    rmlveggat         => vgat%rmlveg,                   & !< real, dimension(:,:) : 

    twarmmgat            => vgat%twarmm,                & !< real, dimension(:) : 
    tcoldmgat            => vgat%tcoldm,                & !< real, dimension(:) : 
    gdd5gat              => vgat%gdd5,                  & !< real, dimension(:) : 
    ariditygat           => vgat%aridity,               & !< real, dimension(:) : 
    srplsmongat          => vgat%srplsmon,              & !< real, dimension(:) : 
    defctmongat          => vgat%defctmon,              & !< real, dimension(:) : 
    anndefctgat          => vgat%anndefct,              & !< real, dimension(:) : 
    annsrplsgat          => vgat%annsrpls,              & !< real, dimension(:) : 
    annpcpgat            => vgat%annpcp,                & !< real, dimension(:) : 
    dry_season_lengthgat => vgat%dry_season_length,     & !< real, dimension(:) : 

    bnffreegat          => vgat%bnf_free,               & !< real, dimension(:,:) :
    bnfnatgat           => vgat%bnf_nat,                & !< real, dimension(:,:) :
    bnfantgat           => vgat%bnf_ant,                & !< real, dimension(:,:) : 
    bnftotgat           => vgat%bnf_tot,                & !< real, dimension(:,:) :
    nstressgat          => vgat%nstress,                & !< real, dimension(:,:) :
    nitrifveggat        => vgat%nitrifveg,              & !< real, dimension(:,:) : 
    no_nitveggat        => vgat%no_nitveg,              & !< real, dimension(:,:) : 
    no_denitveggat      => vgat%no_denitveg,            & !< real, dimension(:,:) : 
    no_nitdenitveggat   => vgat%no_nitdenitveg,         & !< real, dimension(:,:) : 
    n2o_nitveggat       => vgat%n2o_nitveg,             & !< real, dimension(:,:) : 
    n2o_denitveggat     => vgat%n2o_denitveg,           & !< real, dimension(:,:) : 
    n2o_nitdenitveggat  => vgat%n2o_nitdenitveg,        & !< real, dimension(:,:) : 
    n2_denitveggat      => vgat%n2_denitveg,            & !< real, dimension(:,:) : 
    nvolveggat          => vgat%nvolveg,                & !< real, dimension(:,:) : 
    nleachveggat        => vgat%nleachveg,              & !< real, dimension(:,:) : 
    appl_fertgat        => vgat%appl_fert,              & !< real, dimension(:,:) : 
    ndep_nh4gat         => vgat%ndep_nh4,               & !< real, dimension(:,:) : 
    ndep_no3gat         => vgat%ndep_no3,               & !< real, dimension(:,:) : 
    ndemandveg_wp_nppgat=> vgat%ndemandveg_wp_npp,      & !< real, dimension(:,:) : 
    nuptakeveg_p_nh4gat => vgat%nuptakeveg_p_nh4,       & !< real, dimension(:,:) : 
    nuptakeveg_p_no3gat => vgat%nuptakeveg_p_no3,       & !< real, dimension(:,:) : 
    nuptakeveg_a_actl_nh4gat => vgat%nuptakeveg_a_actl_nh4, & !< real, dimension(:,:) : 
    nuptakeveg_a_actl_no3gat => vgat%nuptakeveg_a_actl_no3, & !< real, dimension(:,:) : 
    nuptakeveggat       => vgat%nuptakeveg,             & !< real, dimension(:,:) : 
    nallocveg_lgat      => vgat%nallocveg_l,            & !< real, dimension(:,:) : 
    nallocveg_sgat      => vgat%nallocveg_s,            & !< real, dimension(:,:) : 
    nallocveg_rgat      => vgat%nallocveg_r,            & !< real, dimension(:,:) : 
    nresorpedveg_sgat   => vgat%nresorpedveg_s,         & !< real, dimension(:,:) : 
    nresorpedveg_rgat   => vgat%nresorpedveg_r,         & !< real, dimension(:,:) : 
    nre_allocveg_s2lgat => vgat%nre_allocveg_s2l,       & !< real, dimension(:,:) : 
    nre_allocveg_r2lgat => vgat%nre_allocveg_r2l,       & !< real, dimension(:,:) : 
    nleafns2sveggat     => vgat%nleafns2sveg,           & !< real, dimension(:,:) : 
    nstemns2sveggat     => vgat%nstemns2sveg,           & !< real, dimension(:,:) : 
    nrootns2sveggat     => vgat%nrootns2sveg,           & !< real, dimension(:,:) : 
    nlitrveg_lgat       => vgat%nlitrveg_l,             & !< real, dimension(:,:) : 
    nlitrveg_sgat       => vgat%nlitrveg_s,             & !< real, dimension(:,:) : 
    nlitrveg_rgat       => vgat%nlitrveg_r,             & !< real, dimension(:,:) : 
    nlitrveggat         => vgat%nlitrveg,               & !< real, dimension(:,:) : 
    gl2bl_grass_nfluxgat=> vgat%gl2bl_grass_nflux,      & !< real, dimension(:,:) : 
    c2nveg_lgat         => vgat%c2nveg_l,               & !< real, dimension(:,:) : 
    c2nveg_sgat         => vgat%c2nveg_s,               & !< real, dimension(:,:) : 
    c2nveg_rgat         => vgat%c2nveg_r,               & !< real, dimension(:,:) : 
    c2nveg_wpgat        => vgat%c2nveg_wp,              & !< real, dimension(:,:) : 
    c2nveg_litrgat      => vgat%c2nveg_litr,            & !< real, dimension(:,:) : 
    c2nveg_humusgat     => vgat%c2nveg_humus,           & !< real, dimension(:,:) : 
    nhumtrsveggat       => vgat%nhumtrsveg,             & !< real, dimension(:,:) : 
    nmineralveg_litrgat => vgat%nmineralveg_litr,       & !< real, dimension(:,:) : 
    nmineralveg_humusgat=> vgat%nmineralveg_humus,      & !< real, dimension(:,:) : 
    netnmineralveg_gat  => vgat%netnmineralveg,         & !< real, dimension(:,:) : 
    nimmobilveg_nh4gat  => vgat%nimmobilveg_nh4,        & !< real, dimension(:,:) : 
    nimmobilveg_no3gat  => vgat%nimmobilveg_no3,        & !< real, dimension(:,:) : 
    fNnetlandveggat     => vgat%fNnetlandveg,           & !< real, dimension(:,:) : 
    redcoeff_vcmaxgat   => vgat%redcoeff_vcmax,         & !< real, dimension(:,:) : 

    tcurm             => vgat%tcurm,                    & !< real, dimension(:) : 
    srpcuryr          => vgat%srpcuryr,                 & !< real, dimension(:) : 
    dftcuryr          => vgat%dftcuryr,                 & !< real, dimension(:) : 
    tmonth            => vgat%tmonth,                   & !< real, dimension(:,:) : 
    anpcpcur          => vgat%anpcpcur,                 & !< real, dimension(:) : 
    anpecur           => vgat%anpecur,                  & !< real, dimension(:) : 
    gdd5cur           => vgat%gdd5cur,                  & !< real, dimension(:) : 
    surmncur          => vgat%surmncur,                 & !< real, dimension(:) : 
    defmncur          => vgat%defmncur,                 & !< real, dimension(:) : 
    srplscur          => vgat%srplscur,                 & !< real, dimension(:) : 
    defctcur          => vgat%defctcur,                 & !< real, dimension(:) : 

    geremortgat       => vgat%geremort,                 & !< real, dimension(:,:) : 
    intrmortgat       => vgat%intrmort,                 & !< real, dimension(:,:) : 
    ccgat             => vgat%cc,                       & !< real, dimension(:,:) : 
    mmgat             => vgat%mm,                       & !< real, dimension(:,:) : 

    !      Outputs
    !qevpacc_m_save    => vrot%qevpacc_m_save,          & !< real, dimension(:,:) : ! FLAG: not used? 

    !      Tile-level variables (denoted by an ending of "_t")

    fsnowacc_t        => ctem_tile%fsnowacc_t,          & !< real, dimension(:) : 
    taaccgat_t        => ctem_tile%taaccgat_t,          & !< real, dimension(:) : 
    uvaccgat_t        => ctem_tile%uvaccgat_t,          & !< real, dimension(:) : 
    vvaccgat_t        => ctem_tile%vvaccgat_t,          & !< real, dimension(:) : 
    tbaraccgat_t      => ctem_tile%tbaraccgat_t,        & !< real, dimension(:,:) : 
    thliqacc_t        => ctem_tile%thliqacc_t,          & !< real, dimension(:,:) : 
    thiceacc_t        => ctem_tile%thiceacc_t,          & !< real, dimension(:,:) : ! Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    ancgvgac_t        => ctem_tile%ancgvgac_t,          & !< real, dimension(:,:) : 
    rmlcgvga_t        => ctem_tile%rmlcgvga_t,          & !< real, dimension(:,:) : 

    CFLUX_GAacc_t     => ctem_tile%CFLUX_GAacc_t,       & !< real, dimension(:) :
    USTARBS_GAacc_t   => ctem_tile%USTARBS_GAacc_t,     & !< real, dimension(:) :
    ROFBacc_t         => ctem_tile%ROFBacc_t,           & !< real, dimension(:) :
    QFCacc_t          => ctem_tile%QFCacc_t,            & !< real, dimension(:,:) :

    anmossac_t        => ctem_tile%anmossac_t,          & !< real, dimension(:) : 
    rmlmossac_t       => ctem_tile%rmlmossac_t,         & !< real, dimension(:) : 
    gppmossac_t       => ctem_tile%gppmossac_t,         & !< real, dimension(:) : 

    ipeatlandgat     => vgat%ipeatland,                 & !< integer, dimension(:) : 
    peatdepgat       => vgat%peatdep,                   & !< real, dimension(:) : 
    anmossgat        => vgat%anmoss,                    & !< real, dimension(:) : 
    rmlmossgat       => vgat%rmlmoss,                   & !< real, dimension(:) : 
    gppmossgat       => vgat%gppmoss,                   & !< real, dimension(:) : 
    nppmossgat       => vgat%nppmoss,                   & !< real, dimension(:) : 
    armossgat        => vgat%armoss,                    & !< real, dimension(:) : 
    litrmsmossgat    => vgat%litrmsmoss,                & !< real, dimension(:) : 
    Cmossmasgat      => vgat%Cmossmas,                  & !< real, dimension(:) : 
    dmossgat         => vgat%dmoss,                     & !< real, dimension(:) : 
    peatSoilCgat     => vgat%peatSoilC,                 & !< real, dimension(:) :
    pddgat           => vgat%pdd,                       & !< real, dimension(:) : 
    ancsmoss         => vgat%ancsmoss,                  & !< real, dimension(:) : 
    angsmoss         => vgat%angsmoss,                  & !< real, dimension(:) : 
    ancmoss          => vgat%ancmoss,                   & !< real, dimension(:) : 
    angmoss          => vgat%angmoss,                   & !< real, dimension(:) : 
    rmlcsmoss        => vgat%rmlcsmoss,                 & !< real, dimension(:) : 
    rmlgsmoss        => vgat%rmlgsmoss,                 & !< real, dimension(:) : 
    rmlcmoss         => vgat%rmlcmoss,                  & !< real, dimension(:) : 
    rmlgmoss         => vgat%rmlgmoss,                  & !< real, dimension(:) : 
    
    tileAgegat       => vgat%tileAgegat,                & !< real(ilg)
    timharvareagat   => vgat%timharvareagat,            & !< real(ilg)
    prsfireareagat   => vgat%prsfireareagat,            & !< real(ilg)

    tracerMossCMassrot   => tracer%mossCMassrot,        & !< real, dimension(:,:) : Tracer mass in moss biomass, \f$kg C/m^2\f$ 
    tracerMossLitrMassrot => tracer%mossLitrMassrot,    & !< real, dimension(:,:) : Tracer mass in moss litter, \f$kg C/m^2\f$ 

    tracerGLeafMassrot   => tracer%gLeafMassrot,        & !< real, dimension(:,:,:) : Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerBLeafMassrot   => tracer%bLeafMassrot,        & !< real, dimension(:,:,:) : Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerStemMassrot    => tracer%stemMassrot,         & !< real, dimension(:,:,:) : Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerRootMassrot    => tracer%rootMassrot,         & !< real, dimension(:,:,:) : Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$ 

    ! allocated with nlat, nmos, iccp2, ignd:
    tracerLitrMassrot    => tracer%litrMassrot,         & !< real, dimension(:,:,:,:) : Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    tracerSoilCMassrot   => tracer%soilCMassrot,        & !< real, dimension(:,:,:,:) : Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    tracerCO2rot         => tracer%tracerCO2rot,        & !< real, dimension(:,:) : Atmopspheric tracer CO2 concentration (units vary) 

    ! allocated with ilg, ...:
    tracerMossCMassgat   => tracer%mossCMassgat,        & !< real, dimension(:) : Tracer mass in moss biomass, \f$kg C/m^2\f$ 
    tracerMossLitrMassgat => tracer%mossLitrMassgat,    & !< real, dimension(:) : Tracer mass in moss litter, \f$kg C/m^2\f$ 
    tracerCO2gat          => tracer%tracerCO2gat,       & !< real, dimension(:) : Atmopspheric tracer CO2 concentration (units vary) 

    ! allocated with ilg, icc:
    tracerGLeafMassgat   => tracer%gLeafMassgat,        & !< real, dimension(:,:) : Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerBLeafMassgat   => tracer%bLeafMassgat,        & !< real, dimension(:,:) : Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerStemMassgat    => tracer%stemMassgat,         & !< real, dimension(:,:) : Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$ 
    tracerRootMassgat    => tracer%rootMassgat,         & !< real, dimension(:,:) : Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$ 
    ! allocated with ilg, iccp2, ignd:
    tracerLitrMassgat    => tracer%litrMassgat,         & !< real, dimension(:,:,:) : Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    tracerSoilCMassgat   => tracer%soilCMassgat         & !< real, dimension(:,:,:) : Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$ 
    
    )

    !    =================================================================================

    ! Set the DynTilInitializeFlag to false to ensure the model isn't accidentaly re-initilized
    ! unless called for by tiledDisturbancePrep below this flag is used to re-initialize the model 
    ! after dynamic tiling runs
    DynTilInitializeFlag = .False.

    if (timberHarvest .or. dynamicTilingOn .or. prescribedFire) then
      call tiledDisturbancePrep(ctem_on, lnduseon, PFTCompetition, timberHarvest, prescribedFire, &
                            dynamicTilingOn, runyr, iday, controlVector, DynTilInitializeFlag, timharvrow, & 
                            timharvarearow, timharvareagat, prsfirerow, prsfirearearow, prsfireareagat, &
                            ALBSROT, ALICROT, ALVCROT, CLAYROT, CMASROT, Cmossmasrow, DRNROT, FAREROT, FCANROT, &
                            GROROT, LNZ0ROT, ORGMROT, PAMNROT, PAMXROT, RCANROT, RHOSROT, ROOTROT, SANDROT, SCANROT, &
                            SDEPROT, SNOROT, SOCIROT, TBARROT, TCANROT, THICROT,THLQROT, TPNDROT, TSNOROT, ZPNDROT, &
                            bleafmasrow, bnfantrow, bnfnatrow, cfluxcgrow, cfluxcsrow, co2i1cgrow, co2i1csrow, &
                            co2i2cgrow, co2i2csrow, colddays_harvestrow, colddays_leaffallrow, dmossrow, fcancmxrow, &
                            flhrloss_nsrow, flhrloss_srow, gleafmas_NSrow, gleafmassrow, grwtheffrow, ipeatlandrow, & 
                            lfstatusrow, litrmassrow, litrmsmossrow, lygleafmasmaxrow, lyrootmassmaxrow, lyrotmasrow, &
                            lystemmassmaxrow, lystmmasrow, maxAnnualActLyrROT, nbleafmasrow,ngleafmas_NSrow, ngleafmassrow, &
                            nh4_massrow, nlitrmassrow, no3_massrow, nrootmass_NSrow, nrootmasssrow, nstemmass_NSrow, &
                            nstemmasssrow, pandaysrow, peatSoilCrow, rootmass_NSrow, rootmasssrow, rothrlosrow, slopefracrow, &
                            soilcmasrow, soilnmasrow, stemmass_NSrow, stemmasssrow, stmhrlosrow, tracerBLeafMassrot, tracerGLeafMassrot, &
                            tracerLitrMassrot, tracerMossCMassrot, tracerMossLitrMassrot, tracerRootMassrot, tracerSoilCMassrot, tracerStemMassrot, &
                            tymaxlairow, RSMNROT, QA50ROT, VPDAROT, VPDBROT, PSGAROT, PSGBROT, ngleafmasrow, nstemmassrow, &
                            nrootmassrow, soilpHrow, gleafmasrow, stemmassrow, rootmassrow, flhrlossrow, TBASROT, CMAIROT, &
                            WSNOROT, ZSNLROT, TSFSROT, TACROT, QACROT, ITCTROT, DLZWROT, &
                            useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow, tileAgegat, NDAY, NCOUNT,veghghtrow)
    end if
#if defined without_agcm_
    if ((N == 0) .or. DynTilInitializeFlag) then
#endif
  
      ! This section should only execute the first time through the time loop.
      ! It also gets executed after dynamic tiling runs to re-initialized the model
      ! upon re-entry to mainCore.f90 and cut down on the number of variables that need to be transferred
  
      !> Assign soil thermal and hydraulic properties on the basis of the
      !! textural information read in for each of the soil layers.
      call soilProperties(THPROT, THRROT, THMROT, BIROT, PSISROT, GRKSROT, & ! Formerly CLASSB
                          THRAROT, HCPSROT, TCSROT, THFCROT, THLWROT, PSIWROT, &
                          DLZWROT, ZBTWROT, &
                          ALGWVROT, ALGWNROT, ALGDVROT, ALGDNROT, &
                          SANDROT, CLAYROT, ORGMROT, SOCIROT, DELZ, ZBOT, &
                          SDEPROT, ISNDROT, IGDRROT, &
                          NLAT, NMOS, 1, NLTEST, NMTEST, IGND, ipeatlandrow, &
                          KsatScalaron)

#if !defined without_agcm_
    if (N == 0) then
#endif
                          
      ! Determine the active layer depth and depth to the frozen water table.
      ! Occurs once per day, usually, but need a value to start off run with.
      call findPermafrostVars(nmtest, nltest, iday)
  
      ! ctem initializations. FLAG: EC - CLASS only runs won't work properly w/o 
      !                             at least some of what is included in the if-block.
      !                             E.g. ipeatlandgat is used in CLASS routines.
      !                             Should review what can be skipped when running CLASS only.
      !if (ctem_on) then
  
        call ctemInit(nltest, nmtest)

        if (DynTilInitializeFlag) then
          call classGatherPrep(ILMOS, JLMOS, FAREROT, & ! Formerly GATPREP
                              NML, 1, NMTEST, NMTEST, NLAT, ILG, 1, NLTEST)
        end if

        ! ctemg1 converts variables from the 'row' format (nlat, nmos, ...)
        ! to the 'gat' format (ilg, ...) which is what the model calculations
        ! are performed on. The ctemg1 subroutine is used to transform the
        ! read in state variables (which come in with the 'row' format from the
        ! various input files).
        
        call ctemg1(gleafmasgat, gleafmasgat_ns, gleafmasgat_s, & ! Out
                   stemmassgat, stemmassgat_ns, stemmassgat_s, &  ! Out
                   rootmassgat, rootmassgat_ns, rootmassgat_s, & ! Out
                   bleafmasgat, fcancmxgat, zbtwgat, & ! Out
                   dlzwgat, sdepgat, grclarea, ailcggat, & ! Out
                   ailcbgat, ailcgat, zolncgat, rmatcgat, & ! Out
                   rmatctemgat, slaigat, bmasveggat, cmasvegcgat, & ! Out
                   veghghtgat, rootdpthgat, alvsctmgat, alirctmgat, & ! Out
                   paicgat, slaicgat, faregat, & ! Out
                   ipeatlandgat, maxAnnualActLyrGAT, & ! Out
                   tracergLeafMassgat, tracerBLeafMassgat, tracerStemMassgat, & ! Out
                   tracerRootMassgat, tracerLitrMassgat, tracerSoilCMassgat, & ! Out
                   tracerMossCMassgat, tracerMossLitrMassgat, & ! Out                     
                   twarmmgat, tcoldmgat, gdd5gat, & ! Out
                   ariditygat, srplsmongat, defctmongat, anndefctgat, & ! Out
                   annsrplsgat, annpcpgat, dry_season_lengthgat, & ! Out
                   litrmsmossgat, Cmossmasgat, dmossgat, peatSoilCgat, & ! Out
                   pandaysgat, lfstatusgat, slopefracgat, pstemmassgat, & ! Out
                   pgleafmassgat,litrmassgat, soilcmasgat, grwtheffgat, & ! Out
                   lygleafmasmaxgat, lystemmassmaxgat, lyrootmassmaxgat, tileAgegat, timharvareagat, prsfireareagat, & ! Out
                   ilmos, jlmos, nml, &! In
                   gleafmasrow,  gleafmas_NSrow, gleafmassrow, & ! In
                   stemmassrow, stemmass_NSrow, stemmasssrow, & ! In
                   rootmassrow, rootmass_NSrow, rootmasssrow, &! In
                   bleafmasrow, fcancmxrow, zbtwrot, &! In 
                   dlzwrot, sdeprot, grclarearow, &! In
                   ailcgrow, ailcbrow, ailcrow, zolncrow, &! In
                   rmatcrow, rmatctemrow, slairow, bmasvegrow, &! In
                   cmasvegcrow, veghghtrow, rootdpthrow, alvsctmrow, &! In
                   alirctmrow, paicrow, slaicrow, FAREROT, &! In
                   ipeatlandrow, maxAnnualActLyrROT, &! In
                   tracergLeafMassrot, tracerBLeafMassrot, tracerStemMassrot, &! In
                   tracerRootMassrot, tracerLitrMassrot, tracerSoilCMassrot, &! In
                   tracerMossCMassrot, tracerMossLitrMassrot, & ! In
                   twarmmrow, tcoldmrow, gdd5row, & ! In
                   aridityrow, srplsmonrow, defctmonrow, anndefctrow, & ! In
                   annsrplsrow, annpcprow, dry_season_lengthrow, & ! In
                   litrmsmossrow, Cmossmasrow, dmossrow, peatSoilCrow, & ! In
                   pandaysrow, lfstatusrow, slopefracrow, pstemmassrow, & ! In
                   pgleafmassrow, litrmassrow, soilcmasrow, grwtheffrow, &
                   lygleafmasmaxrow, lystemmassmaxrow, lyrootmassmaxrow, tileAgerow, timharvarearow, prsfirearearow)! In
  
        !> Find mosaic tile (grid) average vegetation biomass, litter mass, and soil c mass.
        !! Set growth efficiency to some large number so that no growth related mortality
        !! occurs in first year. For peatlands determine the peatdepth and the peat soil 
        !! carbon amounts. Lastly find the maximum daylength for this location.
        ! NOTE: All of the non-structural/structural C pools are not 
        ! presently in use but remain passed in for future use!
        call allometry(gleafmasgat, gleafmasgat_ns, gleafmasgat_s, & ! In
                       bleafmasgat, stemmassgat, stemmassgat_ns, & ! In
                       stemmassgat_s, rootmassgat, rootmassgat_ns, & ! In
                       rootmassgat_s, 1, nml, ilg, zbtwgat, & !
                       sdepgat, fcancmxgat, & ! In
                       ipeatlandgat, maxAnnualActLyrGAT, & ! In
                       ailcggat, ailcbgat, ailcgat, zolncgat, & ! Out
                       rmatcgat, rmatctemgat, slaigat, bmasveggat, & ! Out
                       cmasvegcgat, veghghtgat, rootdpthgat, alvsctmgat, & ! Out
                       alirctmgat, paicgat, slaicgat) ! Out
  
      !end if   ! if (ctem_on)
  
      !     ctem initial preparation done

    end if ! if (N == 0)

    N = N + 1

#if defined without_agcm_
    !> atmosphericVarsCalc evaluates a series of derived atmospheric variables

    call atmosphericVarsCalc(VPDROW, TADPROW, PADRROW, RHOAROW, RHSIROW, & ! Formerly CLASSI
                             RPCPROW, TRPCROW, SPCPROW, TSPCROW, TAROW, QAROW, &
                             PREROW, RPREROW, SPREROW, PRESROW, VMODROW, &
                             IPCP, NLAT, 1, NLTEST)
#endif

    !> classGather performs the gather operation, gathering variables from their
    !> positions as mosaic tiles within the modelled areas to long vectors of mosaic tiles
    
    call classGather(TBARGAT, THLQGAT, THICGAT, TPNDGAT, ZPNDGAT, & ! Formerly CLASSG
                     TBASGAT, ALBSGAT, TSNOGAT, RHOSGAT, SNOGAT, &
                     TCANGAT, RCANGAT, SCANGAT, GROGAT, CMAIGAT, &
                     FCANGAT, LNZ0GAT, ALVCGAT, ALICGAT, PAMXGAT, &
                     PAMNGAT, CMASGAT, ROOTGAT, RSMNGAT, QA50GAT, &
                     VPDAGAT, VPDBGAT, PSGAGAT, PSGBGAT, PAIDGAT, &
                     HGTDGAT, ACVDGAT, ACIDGAT, TSFSGAT, WSNOGAT, &
                     THPGAT, THRGAT, THMGAT, BIGAT, PSISGAT, &
                     GRKSGAT, THRAGAT, HCPSGAT, TCSGAT, IGDRGAT, &
                     THFCGAT, THLWGAT, PSIWGAT, DLZWGAT, ZBTWGAT, &
                     VMODGAT, ZSNLGAT, ZPLGGAT, ZPLSGAT, TACGAT, &
                     QACGAT, DRNGAT, XSLPGAT, GRKFGAT, WFSFGAT, &
                     WFCIGAT, ALGWVGAT, ALGWNGAT, ALGDVGAT, &
                     ALGDNGAT, ASVDGAT, ASIDGAT, AGVDGAT, &
                     AGIDGAT, ISNDGAT, RADJGAT, ZBLDGAT, Z0ORGAT, &
                     ZRFMGAT, ZRFHGAT, ZDMGAT, ZDHGAT, FSVHGAT, &
                     FSIHGAT, FSDBGAT, FSFBGAT, FSSBGAT, CSZGAT, &
                     FSGGAT, FLGGAT, FDLGAT, ULGAT, VLGAT, &
                     TAGAT, QAGAT, PRESGAT, PREGAT, PADRGAT, &
                     VPDGAT, TADPGAT, RHOAGAT, RPCPGAT, TRPCGAT, &
                     SPCPGAT, TSPCGAT, RHSIGAT, FCLOGAT, DLONGAT, &
                     GGEOGAT, GUSTGAT, REFGAT, BCSNGAT, DEPBGAT, &
                     DLATGAT, maxAnnualActLyrGAT, ILMOS, JLMOS, &
                     NML, NLAT, NTLD, NMOS, ILG, IGND, ICAN, ICAN + 1, NBS, &
                     TBARROT, THLQROT, THICROT, TPNDROT, ZPNDROT, &
                     TBASROT, ALBSROT, TSNOROT, RHOSROT, SNOROT, &
                     TCANROT, RCANROT, SCANROT, GROROT, CMAIROT, &
                     FCANROT, LNZ0ROT, ALVCROT, ALICROT, PAMXROT, &
                     PAMNROT, CMASROT, ROOTROT, RSMNROT, QA50ROT, &
                     VPDAROT, VPDBROT, PSGAROT, PSGBROT, PAIDROT, &
                     HGTDROT, ACVDROT, ACIDROT, TSFSROT, WSNOROT, &
                     THPROT, THRROT, THMROT, BIROT, PSISROT, &
                     GRKSROT, THRAROT, HCPSROT, TCSROT, IGDRROT, &
                     THFCROT, THLWROT, PSIWROT, DLZWROT, ZBTWROT, &
                     VMODROW, ZSNLROT, ZPLGROT, ZPLSROT, TACROT, &
                     QACROT, DRNROT, XSLPROT, GRKFROT, WFSFROT, &
                     WFCIROT, ALGWVROT, ALGWNROT, ALGDVROT, &
                     ALGDNROT, ASVDROT, ASIDROT, AGVDROT, &
                     AGIDROT, ISNDROT, RADJROW, ZBLDROW, Z0ORROW, &
                     ZRFMROW, ZRFHROW, ZDMROW, ZDHROW, FSVHROW, &
                     FSIHROW, FSDBROL, FSFBROL, FSSBROL, CSZROW, &
                     FSGROL, FLGROL, FDLROT, ULROW, VLROW, &
                     TAROW, QAROW, PRESROW, PREROW, PADRROW, &
                     VPDROW, TADPROW, RHOAROW, RPCPROW, TRPCROW, &
                     SPCPROW, TSPCROW, RHSIROW, FCLOROW, DLONROW, &
                     GGEOROW, GUSTROL, REFROT, BCSNROW, DEPBROW, &
                     DLATROW, maxAnnualActLyrROT)

    ! Calculate snow fall flux in kg/m2/s. This is used by the CanESM snow processes subroutines.
    if (isnoalb == 1) then
      do K = 1,NML
        PCSNGAT(K) = SPCPGAT(K) * RHSIGAT(K)  ! snowfall rate * density of fresh snow
      end do
    end if

    !    * INITIALIZATION OF DIAGNOSTIC VARIABLES SPLIT OUT OF classGather
    !    * FOR CONSISTENCY WITH GCM APPLICATIONS.
    call initDiagnosticVars(nml, ilg)

    !========================================================================
    
#if defined without_agcm_
    !> energyWaterBalanceCheck does the initial calculations for the energy and water balance checks

    call energyWaterBalanceCheck(0, CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP, & ! Formerly CLASSZ
                                 WTVSTP, WTSSTP, WTGSTP, &
                                 FSGVGAT, FLGVGAT, HFSCGAT, HEVCGAT, HMFCGAT, HTCCGAT, &
                                 FSGSGAT, FLGSGAT, HFSSGAT, HEVSGAT, HMFNGAT, HTCSGAT, &
                                 FSGGGAT, FLGGGAT, HFSGGAT, HEVGGAT, HMFGGAT, HTCGAT, &
                                 PCFCGAT, PCLCGAT, QFCFGAT, QFCLGAT, ROFCGAT, WTRCGAT, &
                                 PCPNGAT, QFNGAT, ROFNGAT, WTRSGAT, PCPGGAT, QFGGAT, &
                                 QFCGAT, ROFGAT, WTRGGAT, CMAIGAT, RCANGAT, SCANGAT, &
                                 TCANGAT, SNOGAT, WSNOGAT, TSNOGAT, THLQGAT, THICGAT, &
                                 HCPSGAT, THPGAT, DLZWGAT, TBARGAT, ZPNDGAT, TPNDGAT, &
                                 DELZ, FCS, FGS, FC, FG, &
                                 1, NML, ILG, IGND, N)
#endif

    ! ctemg2 takes variables in the 'row' format (nlat, nmos, ...)
    ! and converts them to the 'gat' format (ilg, ...). At present
    ! ctemg2 is bloated with many variables that do not require
    ! gathering. This subroutine should ideally be only used for
    ! state variables that are updated from external files as
    ! the run progresses. Since the model calculations operate
    ! on the 'gat' form, any other variables need not be gathered
    ! as they will already be in the correct format from the previous
    ! model timestep.
    call ctemg2(fcancmxgat, ailcgsgat, fcancsgat, fcancgat, &
                co2concgat, co2i1cggat, co2i1csgat, co2i2cggat, &
                co2i2csgat, xdiffusgat, cfluxcggat, &
                cfluxcsgat, ancsveggat, ancgveggat, rmlcsveggat, &
                rmlcgveggat, canresgat, sdepgat, ch4concgat, &
                sandgat, claygat, orgmgat, &
                anveggat, rmlveggat, prbfrhucgat, &
                extnprobgat, pfcancmxgat, nfcancmxgat, &
                stemmassgat, stemmassgat_ns, stemmassgat_s, &
                rootmassgat, rootmassgat_ns, rootmassgat_s, &
                litrmassgat, gleafmasgat, gleafmasgat_ns, gleafmasgat_s, &
                bleafmasgat, soilcmasgat, flhrlossgat, flhrloss_nsgat, flhrloss_sgat, &
                pandaysgat, lfstatusgat, grwtheffgat, lystmmasgat, &
                lyrotmasgat, lmaxtgat, smaxtgat, rmaxtgat, &
                lygleafmasmaxgat, lystemmassmaxgat, lyrootmassmaxgat, &
                tymaxlaigat, vgbiomasgat, &
                gavgltmsgat, stmhrlosgat, colddays_leaffallgat, colddays_harvestgat, &
                rothrlosgat, gavglaigat, nppgat, &
                nepgat, nepCMIPgat, hetroresgat, autoresgat, soilcrespgat, &
                rmgat, rggat, nbpgat, litresgat, &
                socresgat, gppgat, dstcemlsgat, litrfallgat, &
                litrfallveggat, humiftrsgat, rmlgat, &
                rmsgat, rmrgat, tltrleafgat, tltrstemgat, &
                tltrrootgat, leaflitrgat, roottempgat, afrleafgat, &
                afrstemgat, afrrootgat, wtstatusgat, ltstatusgat, &
                burnfracgat, smfuncveggat, lucemcomgat, lucltringat, &
                lucsocingat, lucemcomngat, lucltrinngat, lucsocinngat, dstcemls3gat, popdingat, &
                faregat, gavgscmsgat, rmlvegaccgat, pftexistgat, ccgat, mmgat, &
                rmsveggat, rmrveggat, rgveggat, vgbiomas_veggat, &
                gppveggat, nepveggat, nppveggat, vcmax0gat, &
                emit_co2gat, emit_cogat, emit_ch4gat, emit_nmhcgat, &
                emit_h2gat, emit_noxgat, emit_n2ogat, emit_nh3gat, emit_pm25gat, &
                emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &
                btermgat, ltermgat, mtermgat, daylgat, dayl_maxgat, &
                nbpveggat, hetroresveggat, autoresveggat, litresveggat, &
                soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat, &
                ch4WetSpecgat, slopefracgat, &
                wetfdyngat, ch4WetDyngat, ch4soillsgat, &
                leafns2sgat, stemns2sgat, rootns2sgat, &
                re_alloc_s2lgat, re_alloc_r2lgat, re_alloc_sr2lgat, &
                anmossgat, rmlmossgat, gppmossgat, armossgat, nppmossgat, &
                litrmsmossgat, peatdepgat, Cmossmasgat, dmossgat, &
                peatSoilCgat, pddgat, tracerCO2gat, &
                bnffreegat, bnfnatgat, bnfantgat, bnftotgat, nstressgat, &
                soilpHgat, nfertilgat, ndepositgat, nh4_massgat, &
                no3_massgat, nitrifveggat, &
                no_nitveggat, no_denitveggat, &
                no_nitdenitveggat, &
                n2o_nitveggat, n2o_denitveggat, &
                n2o_nitdenitveggat, &
                n2_denitveggat, nvolveggat, &
                nleachveggat, appl_fertgat, ndep_nh4gat, ndep_no3gat, &
                ngleafmasgat, ngleafmasgat_ns, ngleafmasgat_s, &
                nbleafmasgat, nstemmassgat, nstemmassgat_ns, &
                nstemmassgat_s, nrootmassgat, nrootmassgat_ns, &
                nrootmassgat_s, nlitrmassgat, soilnmasgat, &
                nvgbiomas_veggat, &
                ndemandveg_wp_nppgat, &
                nuptakeveg_p_nh4gat, &
                nuptakeveg_p_no3gat, &
                nuptakeveg_a_actl_nh4gat, &
                nuptakeveg_a_actl_no3gat, &
                nuptakeveggat, &
                nallocveg_lgat, nallocveg_sgat, &
                nallocveg_rgat, &
                nresorpedveg_sgat, &
                nresorpedveg_rgat, &
                nre_allocveg_s2lgat, &
                nre_allocveg_r2lgat, nleafns2sveggat, &
                nstemns2sveggat, &
                nrootns2sveggat, nlitrveg_lgat, &
                nlitrveg_sgat, nlitrveg_rgat, &
                nlitrveggat, &
                gl2bl_grass_nfluxgat, &
                c2nveg_lgat, c2nveg_sgat, c2nveg_rgat, &
                c2nveg_wpgat, c2nveg_litrgat, c2nveg_humusgat, &
                nhumtrsveggat, &
                nmineralveg_litrgat, &
                nmineralveg_humusgat, netnmineralveg_gat, nimmobilveg_nh4gat, &
                nimmobilveg_no3gat, &
                fNnetlandveggat, redcoeff_vcmaxgat, &
                tileAgegat, timharvareagat, prsfireareagat, &
                ilmos, jlmos, &
                nml, fcancmxrow, ailcgsrow, fcancsrow, fcancrow, &
                co2concrow, co2i1cgrow, co2i1csrow, co2i2cgrow, &
                co2i2csrow, xdiffus, cfluxcgrow, &
                cfluxcsrow, ancsvegrow, ancgvegrow, rmlcsvegrow, &
                rmlcgvegrow, canresrow, SDEPROT, ch4concrow, &
                SANDROT, CLAYROT, ORGMROT, &
                anvegrow, rmlvegrow, prbfrhucrow, &
                extnprobrow, pfcancmxrow, nfcancmxrow, &
                stemmassrow, stemmass_NSrow, stemmasssrow, rootmassrow, &
                rootmass_NSrow, rootmasssrow, litrmassrow, gleafmasrow, &
                gleafmas_NSrow, gleafmassrow, bleafmasrow, soilcmasrow, flhrlossrow, flhrloss_nsrow, flhrloss_srow, &
                pandaysrow, lfstatusrow, grwtheffrow, lystmmasrow, &
                lyrotmasrow, lmaxtrow, smaxtrow, rmaxtrow, &
                lygleafmasmaxrow, lystemmassmaxrow, lyrootmassmaxrow, &
                tymaxlairow, vgbiomasrow, &
                gavgltmsrow, stmhrlosrow, colddays_leaffallrow, colddays_harvestrow, &
                rothrlosrow, gavglairow, npprow, &
                neprow, nepCMIProw, hetroresrow, autoresrow, soilcresprow, &
                rmrow, rgrow, nbprow, litresrow, &
                socresrow, gpprow, dstcemlsrow, litrfallrow, &
                litrfallvegrow, humiftrsrow, rmlrow, &
                rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                lucsocinrow, lucemcomnrow, lucltrinnrow, lucsocinnrow, dstcemls3row, popdinrow, &
                FAREROT, gavgscmsrow, rmlvegaccrow, pftexistrow, ccrow, mmrow, &
                rmsvegrow, rmrvegrow, rgvegrow, vgbiomas_vegrow, &
                gppvegrow, nepvegrow, nppvegrow, vcmax0row, &
                emit_co2row, emit_corow, emit_ch4row, emit_nmhcrow, &
                emit_h2row, emit_noxrow, emit_n2orow, emit_nh3row, emit_pm25row, &
                emit_tpmrow, emit_tcrow, emit_ocrow, emit_bcrow, &
                btermrow, ltermrow, mtermrow, daylrow, dayl_maxrow, &
                nbpvegrow, hetroresvegrow, autoresvegrow, litresvegrow, &
                soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow, &
                ch4WetSpecrow, slopefracrow, &
                wetfdynrow, ch4WetDynrow, ch4soillsrow, &
                leafns2srow, stemns2srow,  rootns2srow, &
                re_alloc_s2lrow, re_alloc_r2lrow, re_alloc_sr2lrow, &
                anmossrow, rmlmossrow, gppmossrow, armossrow, nppmossrow, &
                litrmsmossrow, peatdeprow, Cmossmasrow, dmossrow, &
                peatSoilCrow, pddrow, tracerCO2rot, &
                bnffreerow, bnfnatrow, &
                bnfantrow, bnftotrow, nstressrow, &
                soilpHrow, nfertilrow, ndepositrow, &
                nh4_massrow, no3_massrow, nitrifvegrow, &
                no_nitvegrow, no_denitvegrow, &
                no_nitdenitvegrow, &
                n2o_nitvegrow, n2o_denitvegrow, &
                n2o_nitdenitvegrow, &
                n2_denitvegrow, nvolvegrow, &
                nleachvegrow, appl_fertrow, ndep_nh4row, ndep_no3row, &
                ngleafmasrow, ngleafmas_NSrow, ngleafmassrow, &
                nbleafmasrow, nstemmassrow, nstemmass_NSrow, nstemmasssrow, &
                nrootmassrow, nrootmass_NSrow, nrootmasssrow, nlitrmassrow, &
                soilnmasrow, nvgbiomas_vegrow, &
                ndemandveg_wp_npprow, &
                nuptakeveg_p_nh4row, &
                nuptakeveg_p_no3row, &
                nuptakeveg_a_actl_nh4row, &
                nuptakeveg_a_actl_no3row, &
                nuptakevegrow, nallocveg_lrow, &
                nallocveg_srow, nallocveg_rrow, &
                nresorpedveg_srow, &
                nresorpedveg_rrow, nre_allocveg_s2lrow, &
                nre_allocveg_r2lrow, &
                nleafns2svegrow, nstemns2svegrow, &
                nrootns2svegrow, nlitrveg_lrow, &
                nlitrveg_srow, nlitrveg_rrow, &
                nlitrvegrow, &
                gl2bl_grass_nfluxrow, &
                c2nveg_lrow, c2nveg_srow, c2nveg_rrow, &
                c2nveg_wprow, c2nveg_litrrow, c2nveg_humusrow, &
                nhumtrsvegrow, &
                nmineralveg_litrrow, &
                nmineralveg_humusrow, &
                netnmineralveg_row, &
                nimmobilveg_nh4row, &
                nimmobilveg_no3row, &
                fNnetlandvegrow, redcoeff_vcmaxrow, &
                tileAgerow, timharvarearow, prsfirearearow)

    !-----------------------------------------------------------------------
    !* ALBEDO AND TRANSMISSIVITY CALCULATIONS; GENERAL VEGETATION
    !* CHARACTERISTICS.

    !     * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM by including: zolnc,
    !     * cmasvegc, alvsctm, alirctm, ipeatlandgat in the arguments.

    !> radiationDriver manages the calculation of albedos and other surface parameters

    call radiationDriver(FC, FG, FCS, FGS, ALVSCN, ALIRCN, & ! Formerly CLASSA
                         ALVSG, ALIRG, ALVSCS, ALIRCS, ALVSSN, ALIRSN, &
                         ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN, &
                         TRVSCS, TRIRCS, FSVF, FSVFS, &
                         RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC, &
                         FRAICS, FSNOCS, DISP, DISPS, ZOMLNC, ZOMLCS, &
                         ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, &
                         CHCAP, CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP, &
                         CWLCPS, CWFCPS, RC, RCS, RBCOEF, FROOT, &
                         FROOTS, ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZSNOW, &
                         WSNOGAT, ALVSGAT, ALIRGAT, HTCCGAT, HTCSGAT, HTCGAT, &
                         ALTG, ALSNO, TRSNOWC, TRSNOWG, &
                         WTRCGAT, WTRSGAT, WTRGGAT, CMAIGAT, FSNOGAT, &
                         FCANGAT, LNZ0GAT, ALVCGAT, ALICGAT, PAMXGAT, PAMNGAT, &
                         CMASGAT, ROOTGAT, RSMNGAT, QA50GAT, VPDAGAT, VPDBGAT, &
                         PSGAGAT, PSGBGAT, PAIDGAT, HGTDGAT, ACVDGAT, ACIDGAT, &
                         ASVDGAT, ASIDGAT, AGVDGAT, AGIDGAT, &
                         ALGWVGAT, ALGWNGAT, ALGDVGAT, ALGDNGAT, &
                         THLQGAT, THICGAT, TBARGAT, RCANGAT, SCANGAT, TCANGAT, &
                         GROGAT, SNOGAT, TSNOGAT, RHOSGAT, ALBSGAT, ZBLDGAT, &
                         Z0ORGAT, ZSNLGAT, ZPLGGAT, ZPLSGAT, &
                         FCLOGAT, TAGAT, VPDGAT, RHOAGAT, CSZGAT, &
                         FSDBGAT, FSFBGAT, REFGAT, BCSNGAT, &
                         FSVHGAT, RADJGAT, DLONGAT, RHSIGAT, DELZ, DLZWGAT, &
                         ZBTWGAT, THPGAT, THMGAT, PSISGAT, BIGAT, PSIWGAT, &
                         HCPSGAT, ISNDGAT, &
                         FCANCMXGAT, ICC, ctem_on, RMATCGAT, ZOLNCGAT, &
                         CMASVEGCGAT, AILCGAT, PAICGAT, NOL2PFTS, &
                         SLAICGAT, AILCGGAT, AILCGSGAT, FCANCGAT, FCANCSGAT, &
                         IDAY, ILG, 1, NML, NBS, &
                         JLAT, N, ICAN, ICAN + 1, IGND, IDISP, IZREF, &
                         IWF, IPAI, IHGT, IALC, IALS, IALG, &
                         ISNOALB, alvsctmgat, alirctmgat, ipeatlandgat, &
                         DSL, LAIPAIRatio, LAISPAISRatio)

    !-----------------------------------------------------------------------
    !          * SURFACE TEMPERATURE AND FLUX CALCULATIONS.

    !          * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM
    !          * by including in the arguments: lfstatus

    !> energyBudgetDriver calls the subroutines associated with the surface energy balance calculations


    call energyBudgetDriver(TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, & ! Formerly CLASST
                            THICEC, THICEG, HCPC, HCPG, TCTOPC, TCBOTC, TCTOPG, TCBOTG, &
                            GZEROC, GZEROG, GZROCS, GZROGS, G12C, G12G, G12CS, G12GS, &
                            G23C, G23G, G23CS, G23GS, QFREZC, QFREZG, QMELTC, QMELTG, &
                            EVAPC, EVAPCG, EVAPG, EVAPCS, EVPCSG, EVAPGS, TCANO, TCANS, &
                            RAICAN, SNOCAN, RAICNS, SNOCNS, CHCAP, CHCAPS, TPONDC, TPONDG, &
                            TPNDCS, TPNDGS, TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                            ITCTGAT, CDHGAT, CDMGAT, HFSGAT, TFXGAT, QEVPGAT, QFSGAT, &
                            PETGAT, GAGAT, EFGAT, GTGAT, QGGAT, &
                            SFCTGAT, SFCUGAT, SFCVGAT, SFCQGAT, SFRHGAT, &
                            GTBS, SFCUBS, SFCVBS, USTARBS, USTARBS_GA, &
                            FSGVGAT, FSGSGAT, FSGGGAT, FLGVGAT, FLGSGAT, FLGGGAT, &
                            HFSCGAT, HFSSGAT, HFSGGAT, HEVCGAT, HEVSGAT, HEVGGAT, HMFCGAT, HMFNGAT, &
                            HTCCGAT, HTCSGAT, HTCGAT, QFCFGAT, QFCLGAT, DRGAT, wtableGAT, ILMOGAT, &
                            UEGAT, HBLGAT, TACGAT, QACGAT, ZRFMGAT, ZRFHGAT, ZDMGAT, ZDHGAT, &
                            VPDGAT, TADPGAT, RHOAGAT, FSVHGAT, FSIHGAT, FDLGAT, ULGAT, VLGAT, &
                            TAGAT, QAGAT, PADRGAT, FC, FG, FCS, FGS, RBCOEF, &
                            FSVF, FSVFS, PRESGAT, VMODGAT, ALVSCN, ALIRCN, ALVSG, ALIRG, &
                            ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, ALIRSC, &
                            TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC, RCS, WTRGGAT, groundHeatFlux, QLWOGAT, &
                            FRAINC, FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP, DISPS, &
                            ZOMLNC, ZOELNC, ZOMLNG, ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS, &
                            TBARGAT, TCTOGAT, TCBOGAT, THLQGAT, THICGAT, TPNDGAT, ZPNDGAT, TBASGAT, TCANGAT, TSNOGAT, TSNBGAT, &
                            ZSNOW, RHOSGAT, WSNOGAT, THPGAT, THRGAT, THMGAT, THFCGAT, THLWGAT, &
                            TRSNOWC, TRSNOWG, ALSNO, FSSBGAT, FROOT, FROOTS, &
                            RADJGAT, PREGAT, HCPSGAT, TCSGAT, TSFSGAT, DELZ, DLZWGAT, ZBTWGAT, &
                            FTEMP, FVAP, RIB, ISNDGAT, &
                            AILCGGAT, AILCGSGAT, FCANCGAT, FCANCSGAT, CO2CONCGAT, CO2I1CGGAT, &
                            CO2I1CSGAT, CO2I2CGGAT, CO2I2CSGAT, CSZGAT, XDIFFUSGAT, SLAIGAT, ICC, &
                            ctem_on, RMATCTEMGAT, FCANCMXGAT, L2MAX, NOL2PFTS, CFLUXCGGAT, &
                            CFLUXCSGAT, CFLUX_GAgat, &
                            ANCSVEGGAT, ANCGVEGGAT, RMLCSVEGGAT, RMLCGVEGGAT, &
                            TCSNOW, GSNOW, ITC, ITCG, ITG, ILG, 1, NML, JLAT, N, ICAN, &
                            IGND, IZREF, ISLFD, NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI, &
                            NBS, ISNOALB, daylgat, dayl_maxgat, &
                            ipeatlandgat, ancsmoss, angsmoss, ancmoss, angmoss, &
                            rmlcsmoss, rmlgsmoss, rmlcmoss, rmlgmoss, &
                            Cmossmasgat, dmossgat, iday, pddgat, redcoeff_vcmaxgat, &
                            GLEAFMASGAT, NGLEAFMASGAT, Ncycle_on, vcmax0gat, PSISGAT, BIGAT, &
                            DSL, DSLC, RB)

    !-----------------------------------------------------------------------
    !          * WATER BUDGET CALCULATIONS.

    !> waterBudgetDriver calls the subroutines associated with the surface water balance calculations

    call waterBudgetDriver(THLQGAT, THICGAT, TBARGAT, TCANGAT, RCANGAT, SCANGAT, & ! Formerly CLASSW
                           ROFGAT, TROFGAT, SNOGAT, TSNOGAT, RHOSGAT, ALBSGAT, &
                           WSNOGAT, ZPNDGAT, TPNDGAT, GROGAT, TBASGAT, GFLXGAT, &
                           PCFCGAT, PCLCGAT, PCPNGAT, PCPGGAT, QFCFGAT, QFCLGAT, &
                           QFNGAT, QFGGAT, QFCGAT, HMFCGAT, HMFGGAT, HMFNGAT, &
                           HTCCGAT, HTCSGAT, HTCGAT, ROFCGAT, ROFNGAT, ROVGGAT, &
                           WTRSGAT, WTRGGAT, ROFOGAT, ROFSGAT, ROFBGAT, &
                           TROOGAT, TROSGAT, TROBGAT, QFSGAT, QFXGAT, RHOAGAT, &
                           TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, &
                           THICEC, THICEG, HCPC, HCPG, RPCPGAT, TRPCGAT, &
                           SPCPGAT, TSPCGAT, PREGAT, TAGAT, RHSIGAT, GGEOGAT, &
                           FC, FG, FCS, FGS, TPONDC, TPONDG, &
                           TPNDCS, TPNDGS, EVAPC, EVAPCG, EVAPG, EVAPCS, &
                           EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG, &
                           RAICAN, SNOCAN, RAICNS, SNOCNS, FSVF, FSVFS, &
                           CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO, &
                           TCANS, CHCAP, CHCAPS, CMASSC, CMASCS, ZSNOW, &
                           GZEROC, GZEROG, GZROCS, GZROGS, G12C, G12G, &
                           G12CS, G12GS, G23C, G23G, G23CS, G23GS, &
                           TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                           ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TSFSGAT, &
                           TCTOPC, TCBOTC, TCTOPG, TCBOTG, FROOT, FROOTS, &
                           THPGAT, THRGAT, THMGAT, BIGAT, PSISGAT, GRKSGAT, &
                           THRAGAT, THFCGAT, DRNGAT, HCPSGAT, DELZ, &
                           DLZWGAT, ZBTWGAT, XSLPGAT, GRKFGAT, WFSFGAT, WFCIGAT, &
                           ISNDGAT, IGDRGAT, &
                           IWF, ILG, 1, NML, N, &
                           JLAT, ICAN, IGND, IGND + 1, IGND + 2, &
                           NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI, &
                           RB, RC, RCS, FRAINC, FSNOWC, FRAICS, FSNOCS, &  
                           LAIPAIRatio, LAISPAISRatio, VMODGAT, ZOMLNS, ZRFMGAT, IZREF)  

    !========================================================================

#if defined without_agcm_
    !> energyWaterBalanceCheck completes the energy and water balance checks for the current time step

    call energyWaterBalanceCheck(1, CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP, & ! Formerly CLASSZ
                                 WTVSTP, WTSSTP, WTGSTP, &
                                 FSGVGAT, FLGVGAT, HFSCGAT, HEVCGAT, HMFCGAT, HTCCGAT, &
                                 FSGSGAT, FLGSGAT, HFSSGAT, HEVSGAT, HMFNGAT, HTCSGAT, &
                                 FSGGGAT, FLGGGAT, HFSGGAT, HEVGGAT, HMFGGAT, HTCGAT, &
                                 PCFCGAT, PCLCGAT, QFCFGAT, QFCLGAT, ROFCGAT, WTRCGAT, &
                                 PCPNGAT, QFNGAT, ROFNGAT, WTRSGAT, PCPGGAT, QFGGAT, &
                                 QFCGAT, ROFGAT, WTRGGAT, CMAIGAT, RCANGAT, SCANGAT, &
                                 TCANGAT, SNOGAT, WSNOGAT, TSNOGAT, THLQGAT, THICGAT, &
                                 HCPSGAT, THPGAT, DLZWGAT, TBARGAT, ZPNDGAT, TPNDGAT, &
                                 DELZ, FCS, FGS, FC, FG, &
                                 1, NML, ILG, IGND, N)
#endif

    if (ctem_on) then

      !> Accumulate variables not already accumulated but which are required by CTEM.
      call accumulateForCTEM(nml,ILMOS)

      if (ncount == nday) then

        ! Find daily averages of accumulated variables for CTEM
        call dayEndCTEMPreparation(nml, nday, ILMOS)

        ! Call Canadian Terrestrial Ecosystem Model which operates at a daily time step,
        ! and uses daily accumulated values of variables simulated by CLASS.
        call ctem(fsnowacc_t, sandgat, ILMOS, & ! In
                  ilg, 1, nml, iday, radjgat, &! In
                  taaccgat_t, dlzwgat, ancgvgac_t, rmlcgvga_t, & ! In
                  zbtwgat, doMethane, & ! In
                  uvaccgat_t, vvaccgat_t, lightng, tbaraccgat_t, &! In
                  sdepgat, spinfast, todfrac, & ! In
                  netrad_gat, preacc_gat, PSISGAT, &! In
                  grclarea, popdingat, isndgat, &! In
                  wetfrac_presgat, slopefracgat, BIGAT, &! In
                  THPGAT, DLATGAT, ch4concgat, &! In
                  THFCGAT, THLWGAT, thliqacc_t, thiceacc_t, &! In
                  ipeatlandgat, anmossac_t, rmlmossac_t, gppmossac_t, &! In
                  wtablegat, maxAnnualActLyrGAT, & ! In
                  PFTCompetition, dofire, lnduseon, inibioclim, & ! In
                  leapnow, useTracer, tracerCO2gat, useStaticPeatDep, &! In
                  pfcancmxgat, nfcancmxgat, & ! In
                  Ncycle_on, CFLUX_GAacc_t, soilpHgat, & ! In
                  USTARBS_GAacc_t, nfertilgat, ndepositgat, & ! In
                  ROFBacc_t, QFCacc_t, co2concgat, daylgat, dayl_maxgat, & ! In
                  prsfireareagat, prescribedFire, & ! In
                  stemmassgat, stemmassgat_ns, stemmassgat_s, & ! In/Out
                  rootmassgat, rootmassgat_ns, rootmassgat_s, litrmassgat, & ! In/Out
                  gleafmasgat, gleafmasgat_ns, gleafmasgat_s, & ! In/Out
                  bleafmasgat, soilcmasgat, ailcggat, ailcgat, & ! In/Out
                  zolncgat, rmatctemgat, rmatcgat, ailcbgat, & ! In/Out
                  flhrlossgat, flhrloss_nsgat, flhrloss_sgat, pandaysgat, lfstatusgat, grwtheffgat, & ! In/Out
                  lystmmasgat, lyrotmasgat, lmaxtgat, smaxtgat, rmaxtgat, & ! In/Out
                  lygleafmasmaxgat, lystemmassmaxgat, lyrootmassmaxgat, & ! In/Out
                  tymaxlaigat, vgbiomasgat, & ! In/Out
                  gavgltmsgat, gavgscmsgat, stmhrlosgat, slaigat, & ! In/Out
                  bmasveggat, cmasvegcgat, colddays_leaffallgat, colddays_harvestgat, rothrlosgat, & ! In/Out
                  fcangat, alvsctmgat, alirctmgat, gavglaigat, &! In/Out
                  Cmossmasgat, litrmsmossgat, peatdepgat, peatSoilCgat, fcancmxgat, &! In/Out
                  geremortgat, intrmortgat, pstemmassgat, pgleafmassgat, &! In/Out
                  tcurm, srpcuryr, dftcuryr, &! In/Out
                  tmonth, anpcpcur, anpecur, gdd5cur, &! In/Out
                  surmncur, defmncur, srplscur, defctcur, &! In/Out
                  ariditygat, srplsmongat, defctmongat, anndefctgat, &! In/Out
                  annsrplsgat, annpcpgat, dry_season_lengthgat, &! In/Out
                  pftexistgat, twarmmgat, tcoldmgat, gdd5gat, nppveggat, &! In/Out
                  tracerStemMassgat, tracerRootMassgat, tracerGLeafMassgat, tracerBLeafMassgat, & ! In/Out
                  tracerSoilCMassgat, tracerLitrMassgat, tracerMossCMassgat, tracerMossLitrMassgat, & ! In/Out
                  leafns2sgat, stemns2sgat, rootns2sgat, & ! In/Out
                  re_alloc_s2lgat, re_alloc_r2lgat, re_alloc_sr2lgat, & ! In/Out
                  bnfnatgat, bnfantgat, nstressgat, & ! In/Out
                  trackTileAge, dynamicTilingOn, tileAgegat, timberharvest, timharvareagat, timharvarearow, & ! In/Out
                  nppgat, nepgat, hetroresgat, autoresgat, &! Out (Primary)
                  soilcrespgat, rmgat, rggat, nbpgat, &! Out (Primary)
                  litresgat, socresgat, gppgat, dstcemlsgat, &! Out (Primary)
                  litrfallgat, humiftrsgat, veghghtgat, rootdpthgat, &! Out (Primary)
                  rmlgat, rmsgat, rmrgat, tltrleafgat, &! Out (Primary)
                  tltrstemgat, tltrrootgat, leaflitrgat, roottempgat, &! Out (Primary)
                  burnfracgat, lucemcomgat, lucltringat, &! Out (Primary)
                  lucsocingat, lucemcomngat, lucltrinngat, lucsocinngat, dstcemls3gat, &! Out (Primary)
                  ch4WetSpecgat, ch4WetDyngat, wetfdyngat, ch4soillsgat, &! Out (Primary)
                  paicgat, slaicgat, &! Out (Primary)
                  emit_co2gat, emit_ch4gat, reprocost, blfltrdt, glfltrdt, &! Out (Primary)
                  glcaemls, blcaemls, rtcaemls, stcaemls, ltrcemls, &  ! Out (Primary)
                  ntchlveg, ntchsveg, ntchrveg, &  ! Out (Primary)
                  emit_cogat, emit_nmhcgat, smfuncveggat, &! Out (Secondary)
                  emit_h2gat, emit_noxgat, emit_n2ogat, emit_nh3gat, emit_pm25gat, &! Out (Secondary)
                  emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &! Out (Secondary)
                  btermgat, ltermgat, mtermgat, burnvegfgat, &! Out (Secondary)
                  litrfallveggat, humiftrsveggat, ltstatusgat,  &! Out (Secondary)
                  afrleafgat, afrstemgat, afrrootgat, wtstatusgat, &! Out (Secondary)
                  rmlvegaccgat, rmsveggat, rmrveggat, rgveggat, &! Out (Secondary)
                  vgbiomas_veggat, &! Out (Secondary)
                  gppveggat, nepveggat, nbpveggat, &! Out (Secondary)
                  hetroresveggat, autoresveggat, litresveggat, soilcresveggat, &! Out (Secondary)
                  nppmossgat, armossgat, &! Out (Secondary)
                  ccgat, mmgat, &!Out
                  ngleafmasgat, ngleafmasgat_ns, ngleafmasgat_s, &! Out
                  nbleafmasgat, nstemmassgat, nstemmassgat_ns, &! Out
                  nstemmassgat_s, nrootmassgat, nrootmassgat_ns, &! Out
                  nrootmassgat_s, nlitrmassgat, soilnmasgat, &! Out
                  bnffreegat, bnftotgat, &! Out
                  nh4_massgat, no3_massgat, &! Out
                  nvolveggat, nleachveggat, &! Out
                  nitrifveggat, no_nitveggat, &! Out
                  no_denitveggat, no_nitdenitveggat, &! Out
                  n2o_nitveggat, n2o_denitveggat, &! Out
                  n2o_nitdenitveggat, &! Out
                  n2_denitveggat, &! Out
                  appl_fertgat, ndep_nh4gat, ndep_no3gat, &! Out
                  ndemandveg_wp_nppgat, &! Out
                  nuptakeveg_p_nh4gat, &! Out
                  nuptakeveg_p_no3gat, &! Out
                  nuptakeveg_a_actl_nh4gat, &! Out
                  nuptakeveg_a_actl_no3gat, &
                  nuptakeveggat, &
                  nallocveg_lgat, &! Out
                  nallocveg_sgat, nallocveg_rgat, &! Out
                  nresorpedveg_sgat, &! Out
                  nresorpedveg_rgat, nre_allocveg_s2lgat, &! Out
                  nre_allocveg_r2lgat, &! Out
                  nleafns2sveggat, nstemns2sveggat, &! Out
                  nrootns2sveggat, nlitrveg_lgat, &! Out
                  nlitrveg_sgat, nlitrveg_rgat, &! Out
                  nlitrveggat, &! Out
                  gl2bl_grass_nfluxgat, &! Out
                  c2nveg_lgat, c2nveg_sgat, &! Out
                  c2nveg_rgat, c2nveg_wpgat, &! Out
                  c2nveg_litrgat, c2nveg_humusgat, &! Out
                  nhumtrsveggat, &! Out
                  nmineralveg_litrgat, &! Out
                  nmineralveg_humusgat, netnmineralveg_gat, nimmobilveg_nh4gat, &! Out
                  nimmobilveg_no3gat, &! Out
                  nvgbiomas_veggat, &! Out
                  fNnetlandveggat, redcoeff_vcmaxgat, faregat, nepCMIPgat) ! Out

        ! Once a year,calculate the 14C lost to decay if using the 14C tracer.
        if (useTracer == 2 .and. &
            iday == lastdoy .and. ncount == nday) call decay14C(1,nml)

      end if  ! if (ncount==nday)
    end if  ! if (ctem_on)

    if (isnoalb == 1) then
      ! Calls the snow aging scheme from CanESM (from within the snowProcessesDriver)
      call snowProcessesDriver(NML, SNOGAT, TSNOGAT, DEPBGAT, ROFNGAT, & ! in
                              ZSNOW, TCSNOW, GSNOW, PCSNGAT, WSNOGAT, blackCdepon, & ! in
                              ALBSGAT, RHOSGAT, BCSNGAT, & ! in/out
                              TZSGAT, REFGAT) ! out
    end if

    !> classScatter performs the scatter operation, scattering the variables from
    !> the long vectors of mosaic tiles back onto the configuration of mosaic tiles within grid cells.

    call classScatter(TBARROT, THLQROT, THICROT, TSFSROT, TPNDROT, & ! Formerly CLASSS
                      ZPNDROT, TBASROT, ALBSROT, TSNOROT, RHOSROT, &
                      SNOROT, GTROT, TCANROT, RCANROT, SCANROT, &
                      GROROT, CMAIROT, TACROT, QACROT, WSNOROT, &
                      REFROT, BCSNROW, EMISROT, SALBROT, CSALROT, &
                      groundHeatFluxROT, &
                      ILMOS, JLMOS, NML, NLAT, NTLD, NMOS, &
                      ILG, IGND, ICAN, ICAN + 1, NBS, &
                      TBARGAT, THLQGAT, THICGAT, TSFSGAT, TPNDGAT, &
                      ZPNDGAT, TBASGAT, ALBSGAT, TSNOGAT, RHOSGAT, &
                      SNOGAT, GTGAT, TCANGAT, RCANGAT, SCANGAT, &
                      GROGAT, CMAIGAT, TACGAT, QACGAT, WSNOGAT, &
                      REFGAT, BCSNGAT, EMISGAT, SALBGAT, CSALGAT,&
                      groundHeatFlux)

    !
    !    * SCATTER OPERATION ON DIAGNOSTIC VARIABLES SPLIT OUT OF
    !    * classScatter FOR CONSISTENCY WITH GCM APPLICATIONS.
    !
    do K = 1,NML
      CDHROT (ILMOS(K),JLMOS(K)) = CDHGAT (K)
      CDMROT (ILMOS(K),JLMOS(K)) = CDMGAT (K)
      HFSROT (ILMOS(K),JLMOS(K)) = HFSGAT (K)
      TFXROT (ILMOS(K),JLMOS(K)) = TFXGAT (K)
      QEVPROT(ILMOS(K),JLMOS(K)) = QEVPGAT(K)
      QFSROT (ILMOS(K),JLMOS(K)) = QFSGAT (K)
      QFXROT (ILMOS(K),JLMOS(K)) = QFXGAT (K)
      PETROT (ILMOS(K),JLMOS(K)) = PETGAT (K)
      GAROT  (ILMOS(K),JLMOS(K)) = GAGAT  (K)
      EFROT  (ILMOS(K),JLMOS(K)) = EFGAT  (K)
      QGROT  (ILMOS(K),JLMOS(K)) = QGGAT  (K)
      ALVSROT(ILMOS(K),JLMOS(K)) = ALVSGAT(K)
      ALIRROT(ILMOS(K),JLMOS(K)) = ALIRGAT(K)
      SFCTROT(ILMOS(K),JLMOS(K)) = SFCTGAT(K)
      SFCUROT(ILMOS(K),JLMOS(K)) = SFCUGAT(K)
      SFCVROT(ILMOS(K),JLMOS(K)) = SFCVGAT(K)
      SFCQROT(ILMOS(K),JLMOS(K)) = SFCQGAT(K)
      SFRHROT(ILMOS(K),JLMOS(K)) = SFRHGAT(K)
      FSNOROT(ILMOS(K),JLMOS(K)) = FSNOGAT(K)
      FSGVROT(ILMOS(K),JLMOS(K)) = FSGVGAT(K)
      FSGSROT(ILMOS(K),JLMOS(K)) = FSGSGAT(K)
      FSGGROT(ILMOS(K),JLMOS(K)) = FSGGGAT(K)
      FLGVROT(ILMOS(K),JLMOS(K)) = FLGVGAT(K)
      FLGSROT(ILMOS(K),JLMOS(K)) = FLGSGAT(K)
      FLGGROT(ILMOS(K),JLMOS(K)) = FLGGGAT(K)
      HFSCROT(ILMOS(K),JLMOS(K)) = HFSCGAT(K)
      HFSSROT(ILMOS(K),JLMOS(K)) = HFSSGAT(K)
      HFSGROT(ILMOS(K),JLMOS(K)) = HFSGGAT(K)
      HEVCROT(ILMOS(K),JLMOS(K)) = HEVCGAT(K)
      HEVSROT(ILMOS(K),JLMOS(K)) = HEVSGAT(K)
      HEVGROT(ILMOS(K),JLMOS(K)) = HEVGGAT(K)
      HMFCROT(ILMOS(K),JLMOS(K)) = HMFCGAT(K)
      HMFNROT(ILMOS(K),JLMOS(K)) = HMFNGAT(K)
      HTCCROT(ILMOS(K),JLMOS(K)) = HTCCGAT(K)
      HTCSROT(ILMOS(K),JLMOS(K)) = HTCSGAT(K)
      PCFCROT(ILMOS(K),JLMOS(K)) = PCFCGAT(K)
      PCLCROT(ILMOS(K),JLMOS(K)) = PCLCGAT(K)
      PCPNROT(ILMOS(K),JLMOS(K)) = PCPNGAT(K)
      PCPGROT(ILMOS(K),JLMOS(K)) = PCPGGAT(K)
      QFGROT (ILMOS(K),JLMOS(K)) = QFGGAT (K)
      QFNROT (ILMOS(K),JLMOS(K)) = QFNGAT (K)
      QFCLROT(ILMOS(K),JLMOS(K)) = QFCLGAT(K)
      QFCFROT(ILMOS(K),JLMOS(K)) = QFCFGAT(K)
      ROFROT (ILMOS(K),JLMOS(K)) = ROFGAT (K)
      ROFOROT(ILMOS(K),JLMOS(K)) = ROFOGAT(K)
      ROFSROT(ILMOS(K),JLMOS(K)) = ROFSGAT(K)
      ROFBROT(ILMOS(K),JLMOS(K)) = ROFBGAT(K)
      TROFROT(ILMOS(K),JLMOS(K)) = TROFGAT(K)
      TROOROT(ILMOS(K),JLMOS(K)) = TROOGAT(K)
      TROSROT(ILMOS(K),JLMOS(K)) = TROSGAT(K)
      TROBROT(ILMOS(K),JLMOS(K)) = TROBGAT(K)
      ROFCROT(ILMOS(K),JLMOS(K)) = ROFCGAT(K)
      ROFNROT(ILMOS(K),JLMOS(K)) = ROFNGAT(K)
      ROVGROT(ILMOS(K),JLMOS(K)) = ROVGGAT(K)
      WTRCROT(ILMOS(K),JLMOS(K)) = WTRCGAT(K)
      WTRSROT(ILMOS(K),JLMOS(K)) = WTRSGAT(K)
      WTRGROT(ILMOS(K),JLMOS(K)) = WTRGGAT(K)
      DRROT  (ILMOS(K),JLMOS(K)) = DRGAT  (K)
      wtableROT(ILMOS(K),JLMOS(K)) = wtableGAT(K)
      ILMOROT(ILMOS(K),JLMOS(K)) = ILMOGAT(K)
      UEROT  (ILMOS(K),JLMOS(K)) = UEGAT(K)
      HBLROT (ILMOS(K),JLMOS(K)) = HBLGAT(K)
      TCSNROT(ILMOS(K),JLMOS(K)) = TCSNOW(K)
      TSNBROT(ILMOS(K),JLMOS(K)) = TSNBGAT(K)
    end do ! loop 380

    do L = 1,IGND
      do K = 1,NML
        HMFGROT(ILMOS(K),JLMOS(K),L) = HMFGGAT(K,L)
        HTCROT (ILMOS(K),JLMOS(K),L) = HTCGAT (K,L)
        QFCROT (ILMOS(K),JLMOS(K),L) = QFCGAT (K,L)
        GFLXROT(ILMOS(K),JLMOS(K),L) = GFLXGAT(K,L)
        TCTOROT(ILMOS(K),JLMOS(K),L) = TCTOGAT(K,L)
        TCBOROT(ILMOS(K),JLMOS(K),L) = TCBOGAT(K,L)
      end do
    end do ! loop 390

    do M = 1,50
      do L = 1,6
        do K = 1,NML
          ITCTROT(ILMOS(K),JLMOS(K),L,M) = ITCTGAT(K,L,M)
        end do ! loop 410
      end do ! loop 420
    end do ! loop 430

    ! ctems2 converts variables from the 'gat' format to the
    ! 'row' format, which is suitable for writing to output/restart
    ! files. If a variable is not written to either of those files,
    ! there is no need to scatter the variable as it will be in the
    ! correct format for model calclations ('gat').
    call ctems2(fcancmxrow, rmatcrow, zolncrow, paicrow, &
                ailcrow, ailcgrow, cmasvegcrow, slaicrow, &
                ailcgsrow, rmatctemrow, &
                co2concrow, co2i1cgrow, co2i1csrow, co2i2cgrow, &
                co2i2csrow, xdiffus, slairow, cfluxcgrow, &
                cfluxcsrow, ancsvegrow, ancgvegrow, rmlcsvegrow, &
                rmlcgvegrow, canresrow, SDEPROT, ch4concrow, &
                SANDROT, CLAYROT, ORGMROT, &
                anvegrow, rmlvegrow, prbfrhucrow, &
                extnprobrow, pfcancmxrow, nfcancmxrow, &
                stemmassrow, stemmass_NSrow, stemmasssrow, rootmassrow, &
                rootmass_NSrow, rootmasssrow, litrmassrow, gleafmasrow, &
                gleafmas_NSrow, gleafmassrow, &
                bleafmasrow, soilcmasrow, ailcbrow, flhrlossrow, flhrloss_nsrow, flhrloss_srow, &
                pandaysrow, lfstatusrow, grwtheffrow, lystmmasrow, &
                lyrotmasrow, lmaxtrow, smaxtrow, rmaxtrow, &
                lygleafmasmaxrow, lystemmassmaxrow, lyrootmassmaxrow, &
                tymaxlairow, vgbiomasrow, &
                gavgltmsrow, stmhrlosrow, bmasvegrow, &
                colddays_leaffallrow, colddays_harvestrow, rothrlosrow, &
                alvsctmrow, alirctmrow, gavglairow, npprow, &
                neprow, nepCMIProw, hetroresrow, autoresrow, soilcresprow, &
                rmrow, rgrow, nbprow, litresrow, &
                socresrow, gpprow, dstcemlsrow, litrfallrow, &
                humiftrsrow, veghghtrow, rootdpthrow, rmlrow, &
                litrfallvegrow, humiftrsvegrow, &
                rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                lucsocinrow, lucemcomnrow, lucltrinnrow, lucsocinnrow, nppvegrow, dstcemls3row, &
                FAREROT, gavgscmsrow, &
                rmlvegaccrow, rmsvegrow, rmrvegrow, rgvegrow, &
                vgbiomas_vegrow, &
                gppvegrow, vcmax0row, nepvegrow, &
                FCANROT, pftexistrow, ccrow, mmrow, &
                emit_co2row, emit_corow, emit_ch4row, emit_nmhcrow, &
                emit_h2row, emit_noxrow, emit_n2orow, emit_nh3row, emit_pm25row, &
                emit_tpmrow, emit_tcrow, emit_ocrow, emit_bcrow, &
                btermrow, ltermrow, mtermrow, &
                nbpvegrow, hetroresvegrow, autoresvegrow, litresvegrow, &
                soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow, &
                ch4WetSpecrow, wetfdynrow, ch4WetDynrow, ch4soillsrow, &
                twarmmrow, tcoldmrow, gdd5row, &
                aridityrow, srplsmonrow, defctmonrow, anndefctrow, &
                annsrplsrow, annpcprow, dry_season_lengthrow, &
                leafns2srow, stemns2srow, rootns2srow, &
                re_alloc_s2lrow, re_alloc_r2lrow, re_alloc_sr2lrow, &
                anmossrow, rmlmossrow, gppmossrow, armossrow, nppmossrow, &
                peatdeprow, litrmsmossrow, Cmossmasrow, dmossrow, &
                peatSoilCrow, pddrow, wetfrac_presrow, &
                tracergLeafMassrot, tracerBLeafMassrot, tracerStemMassrot, &
                tracerRootMassrot, tracerLitrMassrot, tracerSoilCMassrot, &
                tracerMossCMassrot, tracerMossLitrMassrot, &
                bnffreerow, bnfnatrow, bnfantrow, bnftotrow, nstressrow, &
                nh4_massrow, no3_massrow, &
                nitrifvegrow, soilpHrow, &
                no_nitvegrow, &
                no_denitvegrow, no_nitdenitvegrow, &
                n2o_nitvegrow, &
                n2o_denitvegrow, n2o_nitdenitvegrow, &
                n2_denitvegrow, &
                nvolvegrow, nleachvegrow, &
                appl_fertrow, ndep_nh4row, ndep_no3row, &
                ngleafmasrow, ngleafmas_NSrow, ngleafmassrow, &
                nbleafmasrow, nstemmassrow, nstemmass_NSrow, &
                nstemmasssrow, nrootmassrow, nrootmass_NSrow, &
                nrootmasssrow, nlitrmassrow, soilnmasrow, &
                nvgbiomas_vegrow, &
                ndemandveg_wp_npprow, &
                nuptakeveg_p_nh4row, &
                nuptakeveg_p_no3row, &
                nuptakeveg_a_actl_nh4row, &
                nuptakeveg_a_actl_no3row, &
                nuptakevegrow, &
                nallocveg_lrow, nallocveg_srow, &
                nallocveg_rrow, &
                nresorpedveg_srow, nresorpedveg_rrow, &
                nre_allocveg_s2lrow, &
                nre_allocveg_r2lrow, &
                nleafns2svegrow, nstemns2svegrow, &
                nrootns2svegrow, &
                nlitrveg_lrow, nlitrveg_srow, &
                nlitrveg_rrow, nlitrvegrow, &
                gl2bl_grass_nfluxrow, &
                c2nveg_lrow, c2nveg_srow, c2nveg_rrow, &
                c2nveg_wprow, c2nveg_litrrow, c2nveg_humusrow, &
                nhumtrsvegrow, &
                nmineralveg_litrrow, &
                nmineralveg_humusrow, &
                netnmineralveg_row, &
                nimmobilveg_nh4row, &
                nimmobilveg_no3row, &
                fNnetlandvegrow, redcoeff_vcmaxrow, &
                tileAgerow, timharvarearow, & 
                !    ----
                ilmos, jlmos, &
                nml, fcancmxgat, rmatcgat, zolncgat, paicgat, &
                ailcgat, ailcggat, cmasvegcgat, slaicgat, &
                ailcgsgat, rmatctemgat, &
                co2concgat, co2i1cggat, co2i1csgat, co2i2cggat, &
                co2i2csgat, xdiffusgat, slaigat, cfluxcggat, &
                cfluxcsgat, ancsveggat, ancgveggat, rmlcsveggat, &
                rmlcgveggat, canresgat, sdepgat, ch4concgat, &
                sandgat, claygat, orgmgat, &
                anveggat, rmlveggat, prbfrhucgat, &
                extnprobgat, pfcancmxgat, nfcancmxgat, &
                stemmassgat, stemmassgat_ns, stemmassgat_s, &
                rootmassgat, rootmassgat_ns, rootmassgat_s, &
                litrmassgat, gleafmasgat, &
                gleafmasgat_ns, gleafmasgat_s, &
                bleafmasgat, soilcmasgat, ailcbgat, flhrlossgat, flhrloss_nsgat, flhrloss_sgat, &
                pandaysgat, lfstatusgat, grwtheffgat, lystmmasgat, &
                lyrotmasgat, lmaxtgat, smaxtgat, rmaxtgat, &
                lygleafmasmaxgat, lystemmassmaxgat, lyrootmassmaxgat, &
                tymaxlaigat, vgbiomasgat, &
                gavgltmsgat, stmhrlosgat, bmasveggat, &
                colddays_leaffallgat, colddays_harvestgat, rothrlosgat, &
                alvsctmgat, alirctmgat, gavglaigat, nppgat, &
                nepgat, nepCMIPgat, hetroresgat, autoresgat, soilcrespgat, &
                rmgat, rggat, nbpgat, litresgat, &
                socresgat, gppgat, dstcemlsgat, litrfallgat, &
                humiftrsgat, veghghtgat, rootdpthgat, rmlgat, &
                litrfallveggat, humiftrsveggat, &
                rmsgat, rmrgat, tltrleafgat, tltrstemgat, &
                tltrrootgat, leaflitrgat, roottempgat, afrleafgat, &
                afrstemgat, afrrootgat, wtstatusgat, ltstatusgat, &
                burnfracgat, smfuncveggat, lucemcomgat, lucltringat, &
                lucsocingat, lucemcomngat, lucltrinngat, lucsocinngat, nppveggat, dstcemls3gat, &
                faregat, gavgscmsgat, &
                rmlvegaccgat, rmsveggat, rmrveggat, rgveggat, &
                vgbiomas_veggat, &
                gppveggat, vcmax0gat, nepveggat, &
                fcangat, pftexistgat, ccgat, mmgat, &
                emit_co2gat, emit_cogat, emit_ch4gat, emit_nmhcgat, &
                emit_h2gat, emit_noxgat, emit_n2ogat, emit_nh3gat, emit_pm25gat, &
                emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &
                btermgat, ltermgat, mtermgat, &
                nbpveggat, hetroresveggat, autoresveggat, litresveggat, &
                soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat, &
                ch4WetSpecgat, wetfdyngat, ch4WetDyngat, ch4soillsgat, &
                twarmmgat, tcoldmgat, gdd5gat, &
                ariditygat, srplsmongat, defctmongat, anndefctgat, &
                annsrplsgat, annpcpgat, dry_season_lengthgat, &
                leafns2sgat, stemns2sgat, rootns2sgat, &
                re_alloc_s2lgat, re_alloc_r2lgat, re_alloc_sr2lgat, &
                anmossgat, rmlmossgat, gppmossgat, armossgat, nppmossgat, &
                peatdepgat, litrmsmossgat, Cmossmasgat, dmossgat, &
                peatSoilCgat, pddgat, wetfrac_presgat, &
                tracergLeafMassgat, tracerBLeafMassgat, tracerStemMassgat, &
                tracerRootMassgat, tracerLitrMassgat, tracerSoilCMassgat, &
                tracerMossCMassgat, tracerMossLitrMassgat, &
                bnffreegat, bnfnatgat, bnfantgat, bnftotgat, nstressgat, &
                nh4_massgat, no3_massgat, &
                nitrifveggat, soilpHgat, &
                no_nitveggat, &
                no_denitveggat, no_nitdenitveggat, &
                n2o_nitveggat, &
                n2o_denitveggat, n2o_nitdenitveggat, &
                n2_denitveggat, &
                nvolveggat, nleachveggat, &
                appl_fertgat, ndep_nh4gat, ndep_no3gat, &
                ngleafmasgat, ngleafmasgat_ns, ngleafmasgat_s, &
                nbleafmasgat, nstemmassgat, nstemmassgat_ns, &
                nstemmassgat_s, nrootmassgat, nrootmassgat_ns, &
                nrootmassgat_s, nlitrmassgat, soilnmasgat, &
                nvgbiomas_veggat, &
                ndemandveg_wp_nppgat, &
                nuptakeveg_p_nh4gat, &
                nuptakeveg_p_no3gat, &
                nuptakeveg_a_actl_nh4gat, &
                nuptakeveg_a_actl_no3gat, &
                nuptakeveggat, &
                nallocveg_lgat, nallocveg_sgat, &
                nallocveg_rgat, &
                nresorpedveg_sgat,nresorpedveg_rgat, &
                nre_allocveg_s2lgat, &
                nre_allocveg_r2lgat, nleafns2sveggat, &
                nstemns2sveggat, &
                nrootns2sveggat, nlitrveg_lgat, &
                nlitrveg_sgat, &
                nlitrveg_rgat, nlitrveggat, &
                gl2bl_grass_nfluxgat, c2nveg_lgat, &
                c2nveg_sgat, c2nveg_rgat, c2nveg_wpgat, &
                c2nveg_litrgat, c2nveg_humusgat, &
                nhumtrsveggat, &
                nmineralveg_litrgat, &
                nmineralveg_humusgat, &
                netnmineralveg_gat, &
                nimmobilveg_nh4gat, &
                nimmobilveg_no3gat, &
                fNnetlandveggat, redcoeff_vcmaxgat, &
                tileAgegat, timharvareagat)
  
    end associate
  return 

  end subroutine main_core_driver

end module mainCore
