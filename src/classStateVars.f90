!> \file
!! Contains the physics variable type structures.
!! @author J. Melton
!!
!! 1. class_rot - CLASS's 'rot' and 'row' vars
!! 2. class_gat - CLASS's 'gat' vars
!! 3. class_out - CLASS's monthly outputs

module classStateVars

  ! G. Meyer  Aug 2022   - Include variables for the dry surface layer (DSL) parameterization
  ! S.R.C.    April 2022 - This module was modified to statically allocate these variables thus removing allocClassVars
  ! E. Chan   Nov 2020   - Declare soil temperature variables (TBARROT/TBARGAT, TBARC/G/CS/GS) with a precision 
  ! R. Harvey Nov 2018     of 64-bits to ensure proper accuracy in simulations of deep soil beyond ~10-20 m.
  !                        See Harvey & Verseghy (2016): https://doi.org/10.1007/s00382-015-2809-5
  !                      - Subroutines must also have these variables declared internally as 64-bit reals 
  !                        if they are passed in via the call.
  !                      - The 64-bit precision is specified using the named constant "real64" from the 
  !                        intrinsic module "iso_fortran_env", renamed to "r8".
  !                      - FLAG: R. Harvey also implemented a method of maintaining 64-bit versions of these
  !                              variables internal to CLASSIC, while being called from within a host model
  !                              which uses only a 32-bit version of the soil temperature TBARROT.
  !                              This has not been implemented since it is not needed when running off-line, 
  !                              but needs to be addressed when running coupled to a host model.
  !                              See Dawson et al. (2018): https://doi.org/10.1007/s00382-017-4034-x
  ! J. Melton Nov 2016

  use classicParams,       only : ican, icp1, nbs, nlat, nmos, ignd, ilg

  use, intrinsic :: iso_fortran_env, only: r8=>real64

  implicit none

  public :: resetClassMon
  public :: resetClassYr
  public :: resetAccVars
  public :: initDiagnosticVars
  public :: initRowVarsPhysics
  public :: classdump

  !=================================================================================
  !> Physics variables in the 'gather' structure
  type class_gather

    ! These will be allocated the dimension: 'ilg'

    integer, dimension(ilg) :: ILMOS     !< Index of grid cell corresponding to current element of gathered vector of land surface variables [ ]
    integer, dimension(ilg) :: JLMOS     !< Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ]
    integer, dimension(ilg) :: IGDRGAT   !< Index of soil layer in which bedrock is encountered

    real, dimension(ilg) :: TZSGAT  !< Vertical temperature gradient in a snow pack
    real, dimension(ilg) :: PCSNGAT !< Snow fall flux \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: ALBSGAT !< Snow albedo [ ]
    real, dimension(ilg) :: CMAIGAT !< Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
    real, dimension(ilg) :: GROGAT  !< Vegetation growth index [ ]
    real, dimension(ilg) :: QACGAT  !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
    real, dimension(ilg) :: RCANGAT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(ilg) :: RHOSGAT !< Density of snow \f$[kg m^{-3} ]\f$
    real, dimension(ilg) :: SCANGAT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(ilg) :: SNOGAT  !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(ilg) :: TACGAT  !< Temperature of air within vegetation canopy [K]
    real, dimension(ilg) :: TBASGAT !< Temperature of bedrock in third soil layer [K]
    real, dimension(ilg) :: TCANGAT !< Vegetation canopy temperature [K]
    real, dimension(ilg) :: TPNDGAT !< Temperature of ponded water [K]
    real, dimension(ilg) :: TSNOGAT !< Snowpack temperature [K]
    real, dimension(ilg) :: WSNOGAT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(ilg) :: maxAnnualActLyrGAT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, dimension(ilg) :: ZPNDGAT !< Depth of ponded water on surface [m]
    real, dimension(ilg) :: REFGAT  !<
    real, dimension(ilg) :: BCSNGAT !< Black carbon mixing ratio \f$[kg m^{-3}]\f$
    real, dimension(ilg) :: AGIDGAT !< Optional user-specified value of ground near-infrared albedo to override CLASS-calculated value [ ]
    real, dimension(ilg) :: AGVDGAT !< Optional user-specified value of ground visible albedo to override CLASS-calculated value [ ]
    real, dimension(ilg) :: ALGDGAT !< Reference albedo for dry soil [ ]
    real, dimension(ilg) :: ALGWGAT !< Reference albedo for saturated soil [ ]
    real, dimension(ilg) :: ASIDGAT !< Optional user-specified value of snow near-infrared albedo to override CLASS-calculated value [ ]
    real, dimension(ilg) :: ASVDGAT !< Optional user-specified value of snow visible albedo to override CLASS-calculated value [ ]
    real, dimension(ilg) :: DRNGAT  !< Drainage index at bottom of soil profile [ ]
    real, dimension(ilg) :: GRKFGAT !< WATROF parameter used when running MESH code [ ]
    real, dimension(ilg) :: WFCIGAT !< WATROF parameter used when running MESH code [ ]
    real, dimension(ilg) :: WFSFGAT !< WATROF parameter used when running MESH code [ ]
    real, dimension(ilg) :: XSLPGAT !< Surface slope (used when running MESH code) [degrees]
    real, dimension(ilg) :: ZPLGGAT !< Maximum water ponding depth for snow-free subareas (user-specified when running MESH code) [m]
    real, dimension(ilg) :: ZPLSGAT !< Maximum water ponding depth for snow-covered subareas (user-specified when running MESH code) [m]
    real, dimension(ilg) :: ZSNLGAT !< Limiting snow depth below which coverage is < 100% [m]
    real, dimension(ilg) :: ALGWVGAT !<
    real, dimension(ilg) :: ALGWNGAT !<
    real, dimension(ilg) :: ALGDVGAT !<
    real, dimension(ilg) :: ALGDNGAT !<
    real, dimension(ilg) :: EMISGAT  !<
    real, dimension(ilg) :: CSZGAT  !< Cosine of solar zenith angle [ ]
    real, dimension(ilg) :: DLONGAT !< Longitude of grid cell (east of Greenwich) [degrees]
    real, dimension(ilg) :: DLATGAT !< Latitude of grid cell [degrees]
    real, dimension(ilg) :: FCLOGAT !< Fractional cloud cover [ ]
    real, dimension(ilg) :: FDLGAT  !< Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FSIHGAT !< Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FSVHGAT !< Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: GGEOGAT !< Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: PADRGAT !< Partial pressure of dry air [Pa]
    real, dimension(ilg) :: PREGAT  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: PRESGAT !< Surface air pressure [Pa]
    real, dimension(ilg) :: QAGAT   !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, dimension(ilg) :: RADJGAT !< Latitude of grid cell (positive north of equator) [rad]
    real, dimension(ilg) :: RHOAGAT !< Density of air \f$[kg m^{-3} ]\f$
    real, dimension(ilg) :: RHSIGAT !< Density of fresh snow \f$[kg m^{-3} ]\f$
    real, dimension(ilg) :: RPCPGAT !< Rainfall rate over modelled area \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: SPCPGAT !< Snowfall rate over modelled area \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: TAGAT   !< Air temperature at reference height [K]
    real, dimension(ilg) :: TADPGAT !< Dew point temperature of air [K]
    real, dimension(ilg) :: TRPCGAT !< Rainfall temperature [K]
    real, dimension(ilg) :: TSPCGAT !< Snowfall temperature [K]
    real, dimension(ilg) :: ULGAT   !< Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: VLGAT   !< Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: VMODGAT !< Wind speed at reference height \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: VPDGAT  !< Vapour pressure deficit [mb]
    real, dimension(ilg) :: Z0ORGAT !< Orographic roughness length [m]
    real, dimension(ilg) :: ZBLDGAT !< Atmospheric blending height for surface roughness length averaging [m]
    real, dimension(ilg) :: ZDHGAT  !< User-specified height associated with diagnosed screen-level variables [m]
    real, dimension(ilg) :: ZDMGAT  !< User-specified height associated with diagnosed anemometer-level wind speed [m]
    real, dimension(ilg) :: ZRFHGAT !< Reference height associated with forcing air temperature and humidity [m]
    real, dimension(ilg) :: ZRFMGAT !< Reference height associated with forcing wind speed [m]
    real, dimension(ilg) :: FSGGAT  !<
    real, dimension(ilg) :: FLGGAT  !<
    real, dimension(ilg) :: GUSTGAT !<
    real, dimension(ilg) :: DEPBGAT !< Black carbon deposition flux \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: GTBS    !<
    real, dimension(ilg) :: SFCUBS  !<
    real, dimension(ilg) :: SFCVBS  !<
    real, dimension(ilg) :: USTARBS !<
    real, dimension(ilg) :: TCSNOW  !<
    real, dimension(ilg) :: GSNOW   !<
    real, dimension(ilg) :: ALIRGAT !< Diagnosed total near-infrared albedo of land surface [ ]
    real, dimension(ilg) :: ALVSGAT !< Diagnosed total visible albedo of land surface [ ]
    real, dimension(ilg) :: CDHGAT  !< Surface drag coefficient for heat [ ]
    real, dimension(ilg) :: CDMGAT  !< Surface drag coefficient for momentum [ ]
    real, dimension(ilg) :: DRGAT   !< Surface drag coefficient under neutral stability [ ]
    real, dimension(ilg) :: EFGAT   !< Evaporation efficiency at ground surface [ ]
    real, dimension(ilg) :: FLGGGAT !< Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FLGSGAT !< Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FLGVGAT !< Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FSGGGAT !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FSGSGAT !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FSGVGAT !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: FSNOGAT !< Diagnosed fractional snow coverage [ ]
    real, dimension(ilg) :: GAGAT   !< Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: GTGAT   !< Diagnosed effective surface black-body temperature [K]
    real, dimension(ilg) :: HBLGAT  !< Height of the atmospheric boundary layer [m]
    real, dimension(ilg) :: HEVCGAT !< Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HEVGGAT !< Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HEVSGAT !< Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HFSGAT  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HFSCGAT !< Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HFSGGAT !< Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HFSSGAT !< Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HMFCGAT !< Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HMFNGAT !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HTCCGAT !< Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: HTCSGAT !< Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: ILMOGAT !< Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    real, dimension(ilg) :: PCFCGAT !< Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: PCLCGAT !< Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: PCPGGAT !< Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: PCPNGAT !< Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: PETGAT  !< Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: QEVPGAT !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: QFCFGAT !< Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: QFCLGAT !< Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: QFGGAT  !< Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: QFNGAT  !< Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: QFSGAT  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: QFXGAT  !< Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: QGGAT   !< Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    real, dimension(ilg) :: ROFGAT  !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: ROFBGAT !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: ROFCGAT !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: ROFNGAT !< Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: ROFOGAT !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: ROFSGAT !< Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: ROVGGAT !< Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: SFCQGAT !< Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    real, dimension(ilg) :: SFCTGAT !< Diagnosed screen-level air temperature [K]
    real, dimension(ilg) :: SFCUGAT !< Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: SFCVGAT !< Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: TFXGAT  !< Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    real, dimension(ilg) :: TROBGAT !< Temperature of base flow from bottom of soil column [K]
    real, dimension(ilg) :: TROFGAT !< Temperature of total runoff [K]
    real, dimension(ilg) :: TROOGAT !< Temperature of overland flow from top of soil column [K]
    real, dimension(ilg) :: TROSGAT !< Temperature of interflow from sides of soil column [K]
    real, dimension(ilg) :: UEGAT   !< Friction velocity of air \f$[m s^{-1} ]\f$
    real, dimension(ilg) :: WTRCGAT !< Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: WTRGGAT !< Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: WTRSGAT !< Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: wtableGAT !< Depth of water table in soil [m]
    real, dimension(ilg) :: QLWOGAT !<
    real, dimension(ilg) :: SFRHGAT !<
    real, dimension(ilg) :: FTEMP   !<
    real, dimension(ilg) :: FVAP    !<
    real, dimension(ilg) :: RIB     !<
    real, dimension(ilg) :: FC      !< Subarea fractional coverage of modelled area - ground under canopy [ ]
    real, dimension(ilg) :: FG      !< Subarea fractional coverage of modelled area - bare ground [ ]
    real, dimension(ilg) :: FCS     !< Subarea fractional coverage of modelled area - snow-covered ground under canopy  [ ]
    real, dimension(ilg) :: FGS     !< Subarea fractional coverage of modelled area - snow-covered bare ground [ ]
    real, dimension(ilg) :: RBCOEF  !<
    real, dimension(ilg) :: ZSNOW   !<
    real, dimension(ilg) :: FSVF    !<
    real, dimension(ilg) :: FSVFS   !<
    real, dimension(ilg) :: ALVSCN  !<
    real, dimension(ilg) :: ALIRCN  !<
    real, dimension(ilg) :: ALVSG   !<
    real, dimension(ilg) :: ALIRG   !<
    real, dimension(ilg) :: ALVSCS  !<
    real, dimension(ilg) :: ALIRCS  !<
    real, dimension(ilg) :: ALVSSN  !<
    real, dimension(ilg) :: ALIRSN  !<
    real, dimension(ilg) :: ALVSGC  !<
    real, dimension(ilg) :: ALIRGC  !<
    real, dimension(ilg) :: ALVSSC  !<
    real, dimension(ilg) :: ALIRSC  !<
    real, dimension(ilg) :: TRVSCN  !<
    real, dimension(ilg) :: TRIRCN  !<
    real, dimension(ilg) :: TRVSCS  !<
    real, dimension(ilg) :: TRIRCS  !<
    real, dimension(ilg) :: RC      !<
    real, dimension(ilg) :: RCS     !<
    real, dimension(ilg) :: FRAINC  !<
    real, dimension(ilg) :: FSNOWC  !<
    real, dimension(ilg) :: FRAICS  !<
    real, dimension(ilg) :: FSNOCS  !<
    real, dimension(ilg) :: CMASSC  !<
    real, dimension(ilg) :: CMASCS  !<
    real, dimension(ilg) :: DISP    !<
    real, dimension(ilg) :: DISPS   !<
    real, dimension(ilg) :: ZOMLNC  !<
    real, dimension(ilg) :: ZOELNC  !<
    real, dimension(ilg) :: ZOMLNG  !<
    real, dimension(ilg) :: ZOELNG  !<
    real, dimension(ilg) :: ZOMLCS  !<
    real, dimension(ilg) :: ZOELCS  !<
    real, dimension(ilg) :: ZOMLNS  !<
    real, dimension(ilg) :: ZOELNS  !<
    real, dimension(ilg) :: TRSNOWC !<
    real, dimension(ilg) :: CHCAP   !<
    real, dimension(ilg) :: CHCAPS  !<
    real, dimension(ilg) :: GZEROC  !< Vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: GZEROG  !< Bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: GZROCS  !< Snow-covered vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: GZROGS  !<  Snow-covered bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: groundHeatFlux  !<  !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(ilg) :: G12C    !<
    real, dimension(ilg) :: G12G    !<
    real, dimension(ilg) :: G12CS   !<
    real, dimension(ilg) :: G12GS   !<
    real, dimension(ilg) :: G23C    !<
    real, dimension(ilg) :: G23G    !<
    real, dimension(ilg) :: G23CS   !<
    real, dimension(ilg) :: G23GS   !<
    real, dimension(ilg) :: QFREZC  !<
    real, dimension(ilg) :: QFREZG  !<
    real, dimension(ilg) :: QMELTC  !<
    real, dimension(ilg) :: QMELTG  !<
    real, dimension(ilg) :: EVAPC   !<
    real, dimension(ilg) :: EVAPCG  !<
    real, dimension(ilg) :: EVAPG   !<
    real, dimension(ilg) :: EVAPCS  !<
    real, dimension(ilg) :: EVPCSG  !<
    real, dimension(ilg) :: EVAPGS  !<
    real, dimension(ilg) :: TCANO   !< Temperature of canopy over ground [K]
    real, dimension(ilg) :: TCANS   !< Temperature of canopy over snow [K]
    real, dimension(ilg) :: RAICAN  !<
    real, dimension(ilg) :: SNOCAN  !<
    real, dimension(ilg) :: RAICNS  !<
    real, dimension(ilg) :: SNOCNS  !<
    real, dimension(ilg) :: CWLCAP  !<
    real, dimension(ilg) :: CWFCAP  !<
    real, dimension(ilg) :: CWLCPS  !<
    real, dimension(ilg) :: CWFCPS  !<
    real, dimension(ilg) :: TSNOCS  !<
    real, dimension(ilg) :: TSNOGS  !<
    real, dimension(ilg) :: RHOSCS  !<
    real, dimension(ilg) :: RHOSGS  !<
    real, dimension(ilg) :: WSNOCS  !<
    real, dimension(ilg) :: WSNOGS  !<
    real, dimension(ilg) :: TPONDC  !<
    real, dimension(ilg) :: TPONDG  !<
    real, dimension(ilg) :: TPNDCS  !<
    real, dimension(ilg) :: TPNDGS  !<
    real, dimension(ilg) :: ZPLMCS  !<
    real, dimension(ilg) :: ZPLMGS  !<
    real, dimension(ilg) :: ZPLIMC  !<
    real, dimension(ilg) :: ZPLIMG  !<
    real, dimension(ilg) :: DSL     !< Thickness of dry surface layer 
    real, dimension(ilg) :: DSLC    !< Thickness of dry surface layer under canopy 
    real, dimension(ilg) :: RB      !< Leaf boundary resistance of vegetation 
    real, dimension(ilg) :: LAIPAIRatio   !< LAI to PAI ratio for canopy over bare ground [] 
    real, dimension(ilg) :: LAISPAISRatio !< LAI to PAI ratio for canopy over snow [] 
    
    !
    !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
    !     * BALANCES.
    !
    real, dimension(ilg) :: CTVSTP !<
    real, dimension(ilg) :: CTSSTP !<
    real, dimension(ilg) :: CT1STP !<
    real, dimension(ilg) :: CT2STP !<
    real, dimension(ilg) :: CT3STP !<
    real, dimension(ilg) :: WTVSTP !<
    real, dimension(ilg) :: WTSSTP !<
    real, dimension(ilg) :: WTGSTP !<

    ! These will be allocated the dimension: 'ignd'
    real, dimension(ignd) :: DELZ    !< Overall thickness of soil layer [m]
    real, dimension(ignd) :: ZBOT    !< Depth of to the bottom of soil layer [m]

    ! These will be allocated the dimension: 'ilg,ignd'
    integer, dimension(ilg,ignd) :: ISNDGAT !< Integer identifier associated with sand content

    real(r8), dimension(ilg,ignd) :: TBARGAT !< Temperature of soil layers [K]

    real, dimension(ilg,ignd) :: THICGAT !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ignd) :: THLQGAT !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ignd) :: BIGAT   !< Clapp and Hornberger empirical “b” parameter [ ]
    real, dimension(ilg,ignd) :: DLZWGAT !< Permeable thickness of soil layer [m]
    real, dimension(ilg,ignd) :: GRKSGAT !< Saturated hydraulic conductivity of soil layers \f$[m s^{-1} ]\f$
    real, dimension(ilg,ignd) :: HCPSGAT !< Volumetric heat capacity of soil particles \f$[J m^{-3} ]\f$
    real, dimension(ilg,ignd) :: PSISGAT !< Soil moisture suction at saturation [m]
    real, dimension(ilg,ignd) :: PSIWGAT !< Soil moisture suction at wilting point [m]
    real, dimension(ilg,ignd) :: TCSGAT  !< Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$\
    real, dimension(ilg,ignd) :: THFCGAT !< Field capacity \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ignd) :: THMGAT  !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ignd) :: THPGAT  !< Pore volume in soil layer \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ignd) :: THRGAT  !< Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ignd) :: THRAGAT !< Fractional saturation of soil behind the wetting front [ ]
    real, dimension(ilg,ignd) :: ZBTWGAT !< Depth to permeable bottom of soil layer [m]
    real, dimension(ilg,ignd) :: THLWGAT !< Soil water content at wilting point, \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ignd) :: GFLXGAT !< Heat conduction between soil layers \f$[W m^{-2} ]\f$
    real, dimension(ilg,ignd) :: HMFGGAT !< Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    real, dimension(ilg,ignd) :: HTCGAT  !< Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, dimension(ilg,ignd) :: QFCGAT  !< Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$

    real(r8), dimension(ilg,ignd) :: TBARC  !<
    real(r8), dimension(ilg,ignd) :: TBARG  !<
    real(r8), dimension(ilg,ignd) :: TBARCS !<
    real(r8), dimension(ilg,ignd) :: TBARGS !<

    real, dimension(ilg,ignd) :: THLIQC !<
    real, dimension(ilg,ignd) :: THLIQG !<
    real, dimension(ilg,ignd) :: THICEC !<
    real, dimension(ilg,ignd) :: THICEG !<
    real, dimension(ilg,ignd) :: FROOT  !<
    real, dimension(ilg,ignd) :: HCPC   !<
    real, dimension(ilg,ignd) :: HCPG   !<
    real, dimension(ilg,ignd) :: FROOTS !<
    real, dimension(ilg,ignd) :: TCTOPC !<
    real, dimension(ilg,ignd) :: TCBOTC !<
    real, dimension(ilg,ignd) :: TCTOPG !<
    real, dimension(ilg,ignd) :: TCBOTG !<

    ! These will be allocated the dimension: 'ilg,ican'
    real, dimension(ilg,ican) :: ACIDGAT !< Optional user-specified value of canopy near-infrared albedo to override CLASS-calculated value [ ]
    real, dimension(ilg,ican) :: ACVDGAT !< Optional user-specified value of canopy visible albedo to override CLASS-calculated value [ ]
    real, dimension(ilg,ican) :: CMASGAT !< Maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$
    real, dimension(ilg,ican) :: HGTDGAT !< Optional user-specified values of height of vegetation categories to override CLASS-calculated values [m]
    real, dimension(ilg,ican) :: PAIDGAT !< Optional user-specified value of plant area indices of vegetation categories to override CLASS-calculated values [ ]
    real, dimension(ilg,ican) :: PAMNGAT !< Minimum plant area index of vegetation category [ ]
    real, dimension(ilg,ican) :: PAMXGAT !< Minimum plant area index of vegetation category [ ]
    real, dimension(ilg,ican) :: PSGAGAT !< Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, dimension(ilg,ican) :: PSGBGAT !< Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, dimension(ilg,ican) :: QA50GAT !< Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$
    real, dimension(ilg,ican) :: ROOTGAT !< Maximum rooting depth of vegetation category [m]
    real, dimension(ilg,ican) :: RSMNGAT !< Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$
    real, dimension(ilg,ican) :: VPDAGAT !< Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, dimension(ilg,ican) :: VPDBGAT !< Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]


    ! These will be allocated the dimension: 'ilg,icp1'
    real, dimension(ilg,icp1) :: ALICGAT !< Background average near-infrared albedo of vegetation category [ ]
    real, dimension(ilg,icp1) :: ALVCGAT !< Background average visible albedo of vegetation category [ ]
    real, dimension(ilg,icp1) :: FCANGAT !< Maximum fractional coverage of modelled area by vegetation category [ ]
    real, dimension(ilg,icp1) :: LNZ0GAT !< Natural logarithm of maximum roughness length of vegetation category [ ]

    ! These will be allocated the dimension: 'ilg,nbs'
    real, dimension(ilg,nbs) :: FSDBGAT !<
    real, dimension(ilg,nbs) :: FSFBGAT !<
    real, dimension(ilg,nbs) :: FSSBGAT !<
    real, dimension(ilg,nbs) :: SALBGAT !<
    real, dimension(ilg,nbs) :: CSALGAT !<
    real, dimension(ilg,nbs) :: ALTG    !<
    real, dimension(ilg,nbs) :: ALSNO   !<
    real, dimension(ilg,nbs) :: TRSNOWG !<

    ! These will be allocated the dimension: 'ilg,4'
    real, dimension(ilg,4) :: TSFSGAT !< Ground surface temperature over subarea [K]

    ! These will be allocated the dimension: 'ilg,6,50'
    integer, dimension(ilg,6,50) :: ITCTGAT !< Counter of number of iterations required to solve surface energy balance for the elements of the four subareas

  end type class_gather

  type (class_gather), save, target :: class_gat
  !$omp threadprivate(class_gat)

  ! ================================================================================
  !> Physics variables in the 'rotated' (rot) structure
  type class_rotated

    ! These will be allocated the dimension: 'nlat'

    real, dimension(nlat) :: CSZROW  !<
    real, dimension(nlat) :: DLONROW !< Longitude of grid cell (east of Greenwich) [degrees]
    real, dimension(nlat) :: DLATROW !< Latitude of grid cell [degrees]
    real, dimension(nlat) :: lonIndexROW!< Index of grid cell being run on the input files grid (longitude)
    real, dimension(nlat) :: latIndexROW!< Index of grid cell being run on the input files grid (latitude)
    real, dimension(nlat) :: FCLOROW !< Fractional cloud cover [ ]
    real, dimension(nlat) :: RHOSROW !< Density of snow \f$[kg m^{-3}]\f$
    real, dimension(nlat) :: GGEOROW !<
    real, dimension(nlat) :: PADRROW !<
    real, dimension(nlat) :: PREROW  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: PRESROW !< Surface air pressure \f$[P_a]\f$
    real, dimension(nlat) :: QAROW   !< Specific humidity at reference height \f$[kg kg^{-1}]\f$
    real, dimension(nlat) :: RADJROW !<
    real, dimension(nlat) :: RHOAROW !<
    real, dimension(nlat) :: RHSIROW !<
    real, dimension(nlat) :: RPCPROW !<
    real, dimension(nlat) :: RPREROW !< Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: SPCPROW !<
    real, dimension(nlat) :: SPREROW !< Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: TAROW   !< Air temperature at reference height [K]
    real, dimension(nlat) :: TSNOROW !< Snowpack temperature [K]
    real, dimension(nlat) :: TCANROW !< Vegetation canopy temperature [K]
    real, dimension(nlat) :: TPNDROW !< Temperature of ponded water [K]
    real, dimension(nlat) :: ZPNDROW !< Depth of ponded water [m]
    real, dimension(nlat) :: SCANROW !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat) :: RCANROW !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat) :: TADPROW !<
    real, dimension(nlat) :: TRPCROW !<
    real, dimension(nlat) :: TSPCROW !<
    real, dimension(nlat) :: ULROW   !< Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, dimension(nlat) :: VLROW   !< Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, dimension(nlat) :: VMODROW !< Wind speed at reference height \f$[m s^{-1} ]\f$
    real, dimension(nlat) :: VPDROW  !<
    real, dimension(nlat) :: ZBLDROW !< Atmospheric blending height for surface roughness length averaging [m]
    real, dimension(nlat) :: ZDHROW  !<
    real, dimension(nlat) :: ZDMROW  !<
    real, dimension(nlat) :: ZRFHROW !< Reference height associated with forcing air temperature and humidity [m]
    real, dimension(nlat) :: ZRFMROW !< Reference height associated with forcing wind speed [m]
    real, dimension(nlat) :: UVROW   !< Wind speed at reference height as read in from file \f$[m s^{-1} ]\f$
    real, dimension(nlat) :: Z0ORROW !<
    real, dimension(nlat) :: PRENROW !<
    !real, dimension(nlat) :: CLDTROW !<
    real, dimension(nlat) :: GUSTROL !<
    real, dimension(nlat) :: DEPBROW !< Black carbon deposition flux \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: BCSNROW !< Black carbon mixing ratio \f$[ kg m^{-3} ]\f$
    real, dimension(nlat) :: ALIRROW !<
    real, dimension(nlat) :: ALVSROW !<
    real, dimension(nlat) :: CDHROW  !<
    real, dimension(nlat) :: CDMROW  !<
    real, dimension(nlat) :: DRROW   !<
    real, dimension(nlat) :: EFROW   !<
    real, dimension(nlat) :: FLGGROW !<
    real, dimension(nlat) :: FLGSROW !<
    real, dimension(nlat) :: FLGVROW !<
    real, dimension(nlat) :: FSGGROW !<
    real, dimension(nlat) :: FSGSROW !<
    real, dimension(nlat) :: FSGVROW !<
    real, dimension(nlat) :: FSNOROW !<
    real, dimension(nlat) :: GAROW   !<
    real, dimension(nlat) :: GTROW   !<
    real, dimension(nlat) :: HBLROW  !<
    real, dimension(nlat) :: HEVCROW !<
    real, dimension(nlat) :: HEVGROW !<
    real, dimension(nlat) :: HEVSROW !<
    real, dimension(nlat) :: HFSROW  !<
    real, dimension(nlat) :: HFSCROW !<
    real, dimension(nlat) :: HFSGROW !<
    real, dimension(nlat) :: HFSSROW !<
    real, dimension(nlat) :: HMFCROW !<
    real, dimension(nlat) :: HMFNROW !<
    real, dimension(nlat) :: HTCCROW !<
    real, dimension(nlat) :: HTCSROW !<
    real, dimension(nlat) :: ILMOROW !<
    real, dimension(nlat) :: PCFCROW !<
    real, dimension(nlat) :: PCLCROW !<
    real, dimension(nlat) :: PCPGROW !<
    real, dimension(nlat) :: PCPNROW !<
    real, dimension(nlat) :: PETROW  !<
    real, dimension(nlat) :: QEVPROW !<
    real, dimension(nlat) :: QFCFROW !<
    real, dimension(nlat) :: QFCLROW !<
    real, dimension(nlat) :: QFGROW  !<
    real, dimension(nlat) :: QFNROW  !<
    real, dimension(nlat) :: QFSROW  !<
    real, dimension(nlat) :: QFXROW  !<
    real, dimension(nlat) :: QGROW   !<
    real, dimension(nlat) :: ROFROW  !<
    real, dimension(nlat) :: ROFBROW !<
    real, dimension(nlat) :: ROFCROW !<
    real, dimension(nlat) :: ROFNROW !<
    real, dimension(nlat) :: ROFOROW !<
    real, dimension(nlat) :: ROFSROW !<
    real, dimension(nlat) :: ROVGROW !<
    real, dimension(nlat) :: SFCQROW !<
    real, dimension(nlat) :: SFCTROW !<
    real, dimension(nlat) :: SFCUROW !<
    real, dimension(nlat) :: SFCVROW !<
    real, dimension(nlat) :: TFXROW  !<
    real, dimension(nlat) :: UEROW   !<
    real, dimension(nlat) :: WTRCROW !<
    real, dimension(nlat) :: WTRGROW !<
    real, dimension(nlat) :: WTRSROW !<
    real, dimension(nlat) :: SFRHROW !<
    real, dimension(nlat) :: SNOROW  !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, dimension(nlat) :: WSNOROW !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    integer, dimension(nlat) :: altotcntr_d     !< Used to count the number of time steps with the sun above the horizon

    real, dimension(nlat) :: FSSROW  !< Total shortwave radiation \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: FDLROW  !< Downwelling longwave sky radiation \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: fracFSFROW  !< Diffuse fraction of total shortwave radiation (if provided) \f$[ ]\f$

    ! These will be allocated the dimension: 'nlat,nmos'
    ! FLAG: These should be renamed to "rot".
    real, dimension(nlat,nmos) :: FSIHROW !< Near infrared shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: FSVHROW !< Visible shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: FDLROT  !< Downwelling longwave sky radiation (tile version of FDLROW) \f$[W m^{-2} ]\f$
    !real, dimension(nlat,nmos) :: FSFROW  !< Diffuse shortwave radiation \f$[W m^{-2} ]\f$ 
    real, dimension(nlat,nmos) :: FSGROL  !< 
    real, dimension(nlat,nmos) :: FLGROL  !< 
    real, dimension(nlat,nmos) :: XDIFFUS !<

    integer, dimension(nlat,nmos) :: IGDRROT !<
    real, dimension(nlat,nmos) :: ALBSROT !< Snow albedo [ ]
    real, dimension(nlat,nmos) :: CMAIROT !<
    real, dimension(nlat,nmos) :: GROROT  !< Vegetation growth index [ ]
    real, dimension(nlat,nmos) :: QACROT  !<
    real, dimension(nlat,nmos) :: RCANROT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat,nmos) :: RHOSROT !< Density of snow \f$[kg m^{-3}]\f$
    real, dimension(nlat,nmos) :: SCANROT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat,nmos) :: SNOROT  !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, dimension(nlat,nmos) :: TACROT  !<
    real, dimension(nlat,nmos) :: TBASROT !<
    real, dimension(nlat,nmos) :: TCANROT !< Vegetation canopy temperature [K]
    real, dimension(nlat,nmos) :: TPNDROT !< Temperature of ponded water [K]
    real, dimension(nlat,nmos) :: TSNOROT !< Snowpack temperature [K]
    real, dimension(nlat,nmos) :: WSNOROT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(nlat,nmos) :: ZPNDROT !< Depth of ponded water [m]
    real, dimension(nlat,nmos) :: REFROT  !<
    real, dimension(nlat,nmos) :: AGIDROT !<
    real, dimension(nlat,nmos) :: AGVDROT !<
    real, dimension(nlat,nmos) :: ALGDROT !<
    real, dimension(nlat,nmos) :: ALGWROT !<
    real, dimension(nlat,nmos) :: ASIDROT !<
    real, dimension(nlat,nmos) :: ASVDROT !<
    real, dimension(nlat,nmos) :: DRNROT  !<
    real, dimension(nlat,nmos) :: FAREROT !< Fractional coverage of mosaic tile on modelled area
    real, dimension(nlat,nmos) :: GRKFROT !<
    real, dimension(nlat,nmos) :: WFCIROT !<
    real, dimension(nlat,nmos) :: WFSFROT !<
    real, dimension(nlat,nmos) :: XSLPROT !<
    real, dimension(nlat,nmos) :: ZPLGROT !<
    real, dimension(nlat,nmos) :: ZPLSROT !<
    real, dimension(nlat,nmos) :: ZSNLROT !<
    real, dimension(nlat,nmos) :: ZSNOROT  !<
    real, dimension(nlat,nmos) :: ALGWVROT !<
    real, dimension(nlat,nmos) :: ALGWNROT !<
    real, dimension(nlat,nmos) :: ALGDVROT !<
    real, dimension(nlat,nmos) :: ALGDNROT !<
    real, dimension(nlat,nmos) :: EMISROT  !<
    real, dimension(nlat,nmos) :: ALIRROT !<
    real, dimension(nlat,nmos) :: ALVSROT !<
    real, dimension(nlat,nmos) :: CDHROT  !<
    real, dimension(nlat,nmos) :: CDMROT  !<
    real, dimension(nlat,nmos) :: DRROT   !<
    real, dimension(nlat,nmos) :: EFROT   !<
    real, dimension(nlat,nmos) :: FLGGROT !<
    real, dimension(nlat,nmos) :: FLGSROT !<
    real, dimension(nlat,nmos) :: FLGVROT !<
    real, dimension(nlat,nmos) :: FSGGROT !<
    real, dimension(nlat,nmos) :: FSGSROT !<
    real, dimension(nlat,nmos) :: FSGVROT !<
    real, dimension(nlat,nmos) :: FSNOROT !<
    real, dimension(nlat,nmos) :: GAROT   !<
    real, dimension(nlat,nmos) :: GTROT   !<
    real, dimension(nlat,nmos) :: HBLROT  !<
    real, dimension(nlat,nmos) :: HEVCROT !<
    real, dimension(nlat,nmos) :: HEVGROT !<
    real, dimension(nlat,nmos) :: HEVSROT !<
    real, dimension(nlat,nmos) :: HFSROT  !<
    real, dimension(nlat,nmos) :: HFSCROT !<
    real, dimension(nlat,nmos) :: HFSGROT !<
    real, dimension(nlat,nmos) :: HFSSROT !<
    real, dimension(nlat,nmos) :: HMFCROT !<
    real, dimension(nlat,nmos) :: HMFNROT !<
    real, dimension(nlat,nmos) :: HTCCROT !<
    real, dimension(nlat,nmos) :: SDEPROT !< Depth to bedrock in the soil profile
    real, dimension(nlat,nmos) :: SOCIROT !<
    real, dimension(nlat,nmos) :: HTCSROT !<
    real, dimension(nlat,nmos) :: ILMOROT !<
    real, dimension(nlat,nmos) :: PCFCROT !<
    real, dimension(nlat,nmos) :: PCLCROT !<
    real, dimension(nlat,nmos) :: PCPGROT !<
    real, dimension(nlat,nmos) :: PCPNROT !<
    real, dimension(nlat,nmos) :: PETROT  !<
    real, dimension(nlat,nmos) :: QEVPROT !<
    real, dimension(nlat,nmos) :: QFCFROT !<
    real, dimension(nlat,nmos) :: QFCLROT !<
    real, dimension(nlat,nmos) :: QFGROT  !<
    real, dimension(nlat,nmos) :: QFNROT  !<
    real, dimension(nlat,nmos) :: QFSROT  !<
    real, dimension(nlat,nmos) :: QFXROT  !<
    real, dimension(nlat,nmos) :: QGROT   !<
    real, dimension(nlat,nmos) :: ROFROT  !<
    real, dimension(nlat,nmos) :: ROFBROT !<
    real, dimension(nlat,nmos) :: ROFCROT !<
    real, dimension(nlat,nmos) :: ROFNROT !<
    real, dimension(nlat,nmos) :: ROFOROT !<
    real, dimension(nlat,nmos) :: ROFSROT !<
    real, dimension(nlat,nmos) :: ROVGROT !<
    real, dimension(nlat,nmos) :: SFCQROT !<
    real, dimension(nlat,nmos) :: SFCTROT !<
    real, dimension(nlat,nmos) :: SFCUROT !<
    real, dimension(nlat,nmos) :: SFCVROT !<
    real, dimension(nlat,nmos) :: TFXROT  !<
    real, dimension(nlat,nmos) :: TROBROT !<
    real, dimension(nlat,nmos) :: TROFROT !<
    real, dimension(nlat,nmos) :: TROOROT !<
    real, dimension(nlat,nmos) :: TROSROT !<
    real, dimension(nlat,nmos) :: UEROT   !<
    real, dimension(nlat,nmos) :: WTRCROT !<
    real, dimension(nlat,nmos) :: WTRGROT !<
    real, dimension(nlat,nmos) :: WTRSROT !<
    real, dimension(nlat,nmos) :: SFRHROT !<
    real, dimension(nlat,nmos) :: wtableROT !< Depth of water table in soil [m]
    real, dimension(nlat,nmos) :: FTABLE !< Depth to frozen water table (m)
    real, dimension(nlat,nmos) :: ACTLYR !< Active layer depth (m)
    real, dimension(nlat,nmos) :: maxAnnualActLyrROT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, dimension(nlat,nmos) :: actLyrThisYrROT !< Annual active layer depth maximum starting from summer solstice for the present year (m)
    real, dimension(nlat,nmos) :: groundHeatFluxROT !< Heat flux at soil surface \f$[W m^{-2} ]\f$

    ! There will be allocated the dimension: 'nlat,nmos,ignd'
    integer, dimension(nlat,nmos,ignd) :: ISNDROT !< Sand content flag, used to delineate non-soils.

    real(r8), dimension(nlat,nmos,ignd) :: TBARROT !< Temperature of soil layers [K]

    real, dimension(nlat,nmos,ignd) :: THICROT !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(nlat,nmos,ignd) :: THLQROT !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(nlat,nmos,ignd) :: BIROT   !<
    real, dimension(nlat,nmos,ignd) :: DLZWROT !< Permeable thickness of soil layer [m]
    real, dimension(nlat,nmos,ignd) :: GRKSROT !<
    real, dimension(nlat,nmos,ignd) :: HCPSROT !<
    real, dimension(nlat,nmos,ignd) :: SANDROT !< Percentage sand content of soil
    real, dimension(nlat,nmos,ignd) :: CLAYROT !< Percentage clay content of soil
    real, dimension(nlat,nmos,ignd) :: ORGMROT !< Percentage organic matter content of soil
    real, dimension(nlat,nmos,ignd) :: PSISROT !<
    real, dimension(nlat,nmos,ignd) :: PSIWROT !<
    real, dimension(nlat,nmos,ignd) :: TCSROT  !<
    real, dimension(nlat,nmos,ignd) :: THFCROT !<
    real, dimension(nlat,nmos,ignd) :: THMROT  !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, dimension(nlat,nmos,ignd) :: THPROT  !<
    real, dimension(nlat,nmos,ignd) :: THRROT  !<
    real, dimension(nlat,nmos,ignd) :: THRAROT !<
    real, dimension(nlat,nmos,ignd) :: ZBTWROT !<
    real, dimension(nlat,nmos,ignd) :: THLWROT !< Soil water content at wilting point, \f$[m^3 m^{-3} ]\f$
    real, dimension(nlat,nmos,ignd) :: GFLXROT !<
    real, dimension(nlat,nmos,ignd) :: HMFGROT !<
    real, dimension(nlat,nmos,ignd) :: HTCROT  !<
    real, dimension(nlat,nmos,ignd) :: QFCROT  !< Water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$

    ! allocated with nlat,nmos,ignd:
    real, dimension(nlat,nmos,ignd) :: TBARACC_M        !< Temperature of soil layers [K] (accumulated for means)
    real, dimension(nlat,nmos,ignd) :: THLQACC_M        !< Volumetric liquid water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)
    real, dimension(nlat,nmos,ignd) :: THICACC_M        !< Volumetric frozen water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)
    !real, dimension(nlat,nmos,ignd) :: tbaraccrow_m     !< Temperature of soil layers [K] (accumulated for CTEM)


    ! These will be allocated the dimension: 'nlat,nmos,ican'
    real, dimension(nlat,nmos,ican) :: ACIDROT !<
    real, dimension(nlat,nmos,ican) :: ACVDROT !<
    real, dimension(nlat,nmos,ican) :: CMASROT !<
    real, dimension(nlat,nmos,ican) :: HGTDROT !<
    real, dimension(nlat,nmos,ican) :: PAIDROT !<
    real, dimension(nlat,nmos,ican) :: PAMNROT !<
    real, dimension(nlat,nmos,ican) :: PAMXROT !<
    real, dimension(nlat,nmos,ican) :: PSGAROT !<
    real, dimension(nlat,nmos,ican) :: PSGBROT !<
    real, dimension(nlat,nmos,ican) :: QA50ROT !<
    real, dimension(nlat,nmos,ican) :: ROOTROT !<
    real, dimension(nlat,nmos,ican) :: RSMNROT !<
    real, dimension(nlat,nmos,ican) :: VPDAROT !<
    real, dimension(nlat,nmos,ican) :: VPDBROT !<

    ! These will be allocated the dimension: 'nlat,nmos,icp1'
    real, dimension(nlat,nmos,icp1) :: ALICROT !<
    real, dimension(nlat,nmos,icp1) :: ALVCROT !<
    real, dimension(nlat,nmos,icp1) :: FCANROT !<
    real, dimension(nlat,nmos,icp1) :: LNZ0ROT !<

    ! These will be allocated the dimension: 'nlat,nmos,nbs'
    real, dimension(nlat,nmos,nbs) :: SALBROT  !<
    real, dimension(nlat,nmos,nbs) :: CSALROT  !<
    real, dimension(nlat,nmos,nbs) :: FSDBROL  !< Direct solar radiation in each modelled wavelength band  [W m-2]
    real, dimension(nlat,nmos,nbs) :: FSFBROL  !< Diffuse solar radiation in each modelled wavelength band  [W m-2]
    real, dimension(nlat,nmos,nbs) :: FSSBROL  !< Total solar radiation in each modelled wavelength band  [W m-2]

    ! These will be allocated the dimension: 'nlat,ignd'

    real, dimension(nlat,ignd) :: TBARROW !< Temperature of soil layers [K]
    real, dimension(nlat,ignd) :: THALROW !< Total volumetric water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(nlat,ignd) :: THICROW !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(nlat,ignd) :: THLQROW !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(nlat,ignd) :: GFLXROW !<
    real, dimension(nlat,ignd) :: HMFGROW !<
    real, dimension(nlat,ignd) :: HTCROW  !<
    real, dimension(nlat,ignd) :: QFCROW  !< Water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$

    ! These will be allocated the dimension: 'nlat,nmos,6,50'
    integer, dimension(nlat,nmos,6,50) :: ITCTROT !<

    ! These will be allocated the dimension: 'nlat,nmos,4'
    real, dimension(nlat,nmos,4)  :: TSFSROT !<

    ! allocated with nlat,nmos:
    real, dimension(nlat,nmos) :: PREACC_M              !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: GTACC_M               !< Diagnosed effective surface black-body temperature [K]
    real, dimension(nlat,nmos) :: QEVPACC_M             !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: HFSACC_M              !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: HMFNACC_M             !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: ROFACC_M              !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: ZSNACC_M              !< Depth of snow pack \f$[ m ]\f$
    real, dimension(nlat,nmos) :: SNOACC_M              !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(nlat,nmos) :: FSNOACC_M              !< Fractional cover of snow pack \f$[fraction]\f$
    real, dimension(nlat,nmos) :: BCSNACC_M             !< Black carbon mixing ratio \f$[kg m^{-3}]\f$ 
    real, dimension(nlat,nmos) :: OVRACC_M              !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(nlat,nmos) :: ROFSACC_M              !< Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: ROFBACC_M              !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: ROFCACC_M              !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: ROFNACC_M              !< Liquid water runoff from snowpack \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: wtableACC_M             !< Depth of water table in soil [m]
    real, dimension(nlat,nmos) :: ALVSACC_M             !< Diagnosed total visible albedo of land surface [ ]
    real, dimension(nlat,nmos) :: ALIRACC_M             !< Diagnosed total near-infrared albedo of land surface [ ]
    real, dimension(nlat,nmos) :: RHOSACC_M             !< Density of snow \f$[kg m^{-3} ]\f$
    real, dimension(nlat,nmos) :: TSNOACC_M             !< Snowpack temperature [K]
    real, dimension(nlat,nmos) :: WSNOACC_M             !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(nlat,nmos) :: TCANACC_M             !< Vegetation canopy temperature [K]
    real, dimension(nlat,nmos) :: RCANACC_M             !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat,nmos) :: SCANACC_M             !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat,nmos) :: ALTOTACC_M            !< Daily broadband albedo
    real, dimension(nlat,nmos) :: ALSNOACC_M            !< Daily snow albedo [ ]
    real, dimension(nlat,nmos) :: GROACC_M              !< Vegetation growth index [ ]
    real, dimension(nlat,nmos) :: FSINACC_M             !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: FLINACC_M             !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: TAACC_M               !< Air temperature at reference height [K]
    real, dimension(nlat,nmos) :: TSURFACC_M            !< Ground surface temperature at reference height [K]
    real, dimension(nlat,nmos) :: UVACC_M               !< Wind speed \f$[m s^{-1} ]\f$
    real, dimension(nlat,nmos) :: PRESACC_M             !< Surface air pressure [Pa]
    real, dimension(nlat,nmos) :: QAACC_M               !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, dimension(nlat,nmos) :: EVAPACC_M             !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(nlat,nmos) :: GROUNDEVAP_M
    real, dimension(nlat,nmos) :: CANOPYEVAP_M
    real, dimension(nlat,nmos) :: TRANSPACC_M
    real, dimension(nlat,nmos) :: FLUTACC_M             !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, dimension(nlat,nmos) :: ZPNDACC_M             !< Depth of ponded water [m]
    real, dimension(nlat,nmos) :: ACTLYR_M              !< Active layer depth [m]

  end type class_rotated

  type (class_rotated), save, target :: class_rot
  !$omp threadprivate(class_rot)

  !=================================================================================
  !> CLASS's monthly and annual outputs
  type class_moyr_output

    !   MONTHLY OUTPUT FOR CLASS GRID-MEAN

    ! allocated with nlat:
    real, dimension(nlat) :: ALVSACC_MO   !< Diagnosed total visible albedo of land surface [ ]
    real, dimension(nlat) :: ALIRACC_MO   !< Diagnosed total near-infrared albedo of land surface [ ]
    real, dimension(nlat) :: FLUTACC_MO   !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: FSINACC_MO   !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: FLINACC_MO   !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: HFSACC_MO    !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: QEVPACC_MO   !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: groundHeatFlux_MO  !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: RCAN_MO      !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat) :: SCAN_MO      !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nlat) :: SNOACC_MO    !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(nlat) :: FSNOACC_MO    !< Fractional cover snow pack \f$[fraction]\f$
    real, dimension(nlat) :: BCSNACC_MO    !< Black carbon mixing ratio \f$[kg m^{-3}]\f$
    real, dimension(nlat) :: WSNOACC_MO   !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(nlat) :: ZSNACC_MO    !< Depth of snow pack \f$[ m ]\f$
    real, dimension(nlat) :: OVRACC_MO    !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(nlat) :: ROFBACC_MO   !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: ROFNACC_MO   !< Liquid water runoff from snowpack \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: ROFACC_MO    !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: PREACC_MO    !< Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(nlat) :: EVAPACC_MO   !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(nlat) :: TRANSPACC_MO !<
    real, dimension(nlat) :: TAACC_MO     !< Air temperature at reference height [K]
    real, dimension(nlat) :: ACTLYR_MO
    real, dimension(nlat) :: FTABLE_MO
    real, dimension(nlat) :: ACTLYR_MIN_MO
    real, dimension(nlat) :: ACTLYR_MAX_MO
    real, dimension(nlat) :: FTABLE_MIN_MO
    real, dimension(nlat) :: FTABLE_MAX_MO
    real, dimension(nlat) :: MRSO_MO      !< Total Soil Moisture Content [kg $m^{-2}$]
    real, dimension(nlat,ignd) :: MRSOL_MO     !< Total water content of soil layer [kg $m^{-2}$]
    real, dimension(nlat) :: ALTOTACC_MO  !< Broadband albedo [ ]
    real, dimension(nlat) :: ALSNOACC_MO  !< Snow albedo [ ]
    real, dimension(nlat) :: GROUNDEVAP   !< evaporation and sublimation from the ground surface (formed from QFG and QFN), kg /m/mon
    real, dimension(nlat) :: CANOPYEVAP   !< evaporation and sublimation from the canopy (formed from QFCL and QFCF), kg /m/mon
    real, dimension(nlat) :: EVSPSBL_MO   !< evaporation and sublimation from the canopy (formed from QFCL and QFCF), kg /m/mon
    real, dimension(nlat) :: TSURFACC_MO    !< Ground surface temperature [K]
    real, dimension(nlat) :: TCANACC_MO     !< Vegetation canopy temperature [K]
    integer, dimension(nlat) :: altotcntr_m !< Used to count the number of time steps with the sun above the horizon

    ! allocated with nlat,ignd:
    real, dimension(nlat,ignd) :: TBARACC_MO !< Temperature of soil layers [K] (accumulated for means)
    real, dimension(nlat,ignd) :: THLQACC_MO !< Volumetric liquid water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)
    real, dimension(nlat,ignd) :: THICACC_MO !< Volumetric frozen water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)

    !   YEARLY OUTPUT FOR CLASS GRID-MEAN

    real, dimension(nlat) :: ALVSACC_YR  !< Diagnosed total visible albedo of land surface [ ]
    real, dimension(nlat) :: ALIRACC_YR  !< Diagnosed total near-infrared albedo of land surface [ ]
    real, dimension(nlat) :: FLUTACC_YR  !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: FSINACC_YR  !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: FLINACC_YR  !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: HFSACC_YR   !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: QEVPACC_YR  !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nlat) :: ROFACC_YR   !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat) :: PREACC_YR   !< Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(nlat) :: EVAPACC_YR  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(nlat) :: TRANSPACC_YR !<
    real, dimension(nlat) :: TAACC_YR    !< Air temperature at reference height [K]
    real, dimension(nlat) :: ACTLYR_YR
    real, dimension(nlat) :: ACTLYR_MIN_YR
    real, dimension(nlat) :: ACTLYR_MAX_YR
    real, dimension(nlat) :: FTABLE_YR
    real, dimension(nlat) :: FTABLE_MIN_YR
    real, dimension(nlat) :: FTABLE_MAX_YR
    real, dimension(nlat) :: ALTOTACC_YR !< Broadband albedo
    integer, dimension(nlat) :: altotcntr_yr !< Used to count the number of time steps with the sun above the horizon

    ! allocated with nlat,ignd:
    real, dimension(nlat,ignd) :: THLQACC_YR !< Volumetric liquid water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)
    real, dimension(nlat,ignd) :: THICACC_YR !< Volumetric frozen water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)

  end type class_moyr_output

  type (class_moyr_output), save, target :: class_out

  !=================================================================================

contains

  !> \ingroup classstatevars_resetClassMon
  !! @{
  !> Resets the CLASS (physics) monthly variables in preparation for the next month
  subroutine resetClassMon (nltest)

    use classicParams,    only : ignd

    implicit none

    integer, intent(in) :: nltest

    integer :: i, j

    do i = 1,nltest
      class_out%ALVSACC_MO(I) = 0.
      class_out%ALIRACC_MO(I) = 0.
      class_out%FLUTACC_MO(I) = 0.
      class_out%FSINACC_MO(I) = 0.
      class_out%FLINACC_MO(I) = 0.
      class_out%HFSACC_MO(I) = 0.
      class_out%QEVPACC_MO(I) = 0.
      class_out%groundHeatFlux_MO(I) = 0.
      class_out%TRANSPACC_MO(I) = 0.
      class_out%SCAN_MO(I) = 0.
      class_out%RCAN_MO(I) = 0.
      class_out%SNOACC_MO(I) = 0.
      class_out%WSNOACC_MO(I) = 0.
      class_out%ZSNACC_MO(I) = 0.
      class_out%OVRACC_MO(I) = 0.
      class_out%ROFBACC_MO(I) = 0.
      class_out%ROFNACC_MO(I) = 0.
      class_out%ROFACC_MO(I) = 0.
      class_out%PREACC_MO(I) = 0.
      class_out%EVAPACC_MO(I) = 0.
      class_out%TAACC_MO(I) = 0.
      class_out%ACTLYR_MO(I) = 0.
      class_out%FTABLE_MO(I) = 0.
      class_out%ACTLYR_MIN_MO(I) = 100000.
      class_out%FTABLE_MIN_MO(I) = 100000.
      class_out%ACTLYR_MAX_MO(I) = 0.
      class_out%FTABLE_MAX_MO(I) = 0.
      class_out%CANOPYEVAP(I) = 0.
      class_out%GROUNDEVAP(I) = 0.
      class_out%EVSPSBL_MO(I) = 0.
      class_out%ALTOTACC_MO(I) = 0.
      class_out%ALSNOACC_MO(I) = 0.
      class_out%altotcntr_m(i) = 0
      class_out%MRSO_MO(i) = 0.
      class_out%TSURFACC_MO(I) = 0.
      class_out%TCANACC_MO(I) = 0.

      do J = 1,IGND
        class_out%TBARACC_MO(I,J) = 0.
        class_out%THLQACC_MO(I,J) = 0.
        class_out%THICACC_MO(I,J) = 0.
        class_out%MRSOL_MO(i,j) = 0.
      end do
    end do

  end subroutine resetClassMon
  !! @}
  !==================================================

  !> \ingroup classstatevars_resetClassYr
  !! @{
  !> Resets the CLASS (physics) annual variables in preparation for the next year
  subroutine resetClassYr (nltest)

    use classicParams,         only : ignd

    implicit none

    integer, intent(in) :: nltest

    integer :: i, j

    do i = 1,nltest
      class_out%ALVSACC_YR(I) = 0.
      class_out%ALIRACC_YR(I) = 0.
      class_out%FLUTACC_YR(I) = 0.
      class_out%FSINACC_YR(I) = 0.
      class_out%FLINACC_YR(I) = 0.
      class_out%HFSACC_YR(I) = 0.
      class_out%QEVPACC_YR(I) = 0.
      class_out%ROFACC_YR(I) = 0.
      class_out%PREACC_YR(I) = 0.
      class_out%EVAPACC_YR(I) = 0.
      class_out%TRANSPACC_YR(I) = 0.
      class_out%TAACC_YR(I) = 0.
      class_out%ALTOTACC_YR(I) = 0.
      class_out%altotcntr_yr(i) = 0
      class_out%ACTLYR_yr(I) = 0.
      class_out%FTABLE_yr(I) = 0.
      class_out%ACTLYR_MIN_yr(I) = 100000.
      class_out%FTABLE_MIN_yr(I) = 100000.
      class_out%ACTLYR_MAX_yr(I) = 0.
      class_out%FTABLE_MAX_yr(I) = 0.
      class_out%FSNOACC_MO(I) = 0.
      class_out%BCSNACC_MO(I) = 0.

      do J = 1,IGND
        class_out%THLQACC_YR(I,J) = 0.
        class_out%THICACC_YR(I,J) = 0.
      end do
    end do

  end subroutine resetClassYr
  !! @}
  !==================================================

  !> \ingroup classstatevars_resetAccVars
  !! @{
  !> Resets the CLASS (physics) aggregation variables in preparation for the next period
  subroutine resetAccVars (nltest,nmtest)

    use classicParams,         only : ignd

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i, m, j

    do I = 1,NLTEST

      class_rot%altotcntr_d(i) = 0

      do M = 1,NMTEST

        class_rot%PREACC_M(i,m) = 0.
        class_rot%GTACC_M(i,m) = 0.
        class_rot%QEVPACC_M(i,m) = 0.
        class_rot%HFSACC_M(i,m) = 0.
        class_rot%HMFNACC_M(i,m) = 0.
        class_rot%ROFACC_M(i,m) = 0.
        class_rot%ZSNACC_M(i,m) = 0.
        class_rot%SNOACC_M(i,m) = 0.
        class_rot%FSNOACC_M(i,m) = 0.
        class_rot%BCSNACC_M(i,m) = 0.
        class_rot%OVRACC_M(i,m) = 0.
        class_rot%ROFSACC_M(i,m) = 0.
        class_rot%ROFBACC_M(i,m) = 0.
        class_rot%ROFCACC_M(i,m) = 0.
        class_rot%ROFNACC_M(i,m) = 0.
        class_rot%wtableACC_M(i,m) = 0.
        class_rot%ALVSACC_M(i,m) = 0.
        class_rot%ALIRACC_M(i,m) = 0.
        class_rot%RHOSACC_M(i,m) = 0.
        class_rot%TSNOACC_M(i,m) = 0.
        class_rot%WSNOACC_M(i,m) = 0.
        class_rot%TCANACC_M(i,m) = 0.
        class_rot%RCANACC_M(i,m) = 0.
        class_rot%SCANACC_M(i,m) = 0.
        class_rot%GROACC_M(i,m) = 0.
        class_rot%FSINACC_M(i,m) = 0.
        class_rot%FLINACC_M(i,m) = 0.
        class_rot%TAACC_M(i,m) = 0.
        class_rot%TSURFACC_M(i,m) = 0.
        class_rot%UVACC_M(i,m) = 0.
        class_rot%PRESACC_M(i,m) = 0.
        class_rot%QAACC_M(i,m) = 0.
        class_rot%ALTOTACC_M(i,m) = 0.
        class_rot%ALSNOACC_M(i,m) = 0.
        class_rot%EVAPACC_M(i,m) = 0.
        class_rot%GROUNDEVAP_M(i,m) = 0.
        class_rot%CANOPYEVAP_M(i,m) = 0.
        class_rot%TRANSPACC_M(i,m) = 0.
        class_rot%ZPNDACC_M(i,m) = 0.
        class_rot%FLUTACC_M(i,m) = 0.
        class_rot%ACTLYR_M(i,m) = 0.

        do J = 1,IGND
          class_rot%TBARACC_M(I,M,J) = 0.
          class_rot%THLQACC_M(I,M,J) = 0.
          class_rot%THICACC_M(I,M,J) = 0.
        end do
      end do
    end do

  end subroutine resetAccVars
  !! @}
  !==================================================
  !> \ingroup classstatevars_initDiagnosticVars
  !! @{
  !> Initialization of diagnostic variables split out of classGather for consistency with gcm applications.
  subroutine initDiagnosticVars (nml, ilg)

    use classicParams,      only : ignd

    implicit none

    integer, intent(in) :: nml, ilg

    integer :: k, m, l

    !    * INITIALIZATION OF DIAGNOSTIC VARIABLES SPLIT OUT OF classGather
    !    * FOR CONSISTENCY WITH GCM APPLICATIONS.

    do K = 1,nml ! loop 330
      class_gat%CDHGAT (K) = 0.0
      class_gat%CDMGAT (K) = 0.0
      class_gat%HFSGAT (K) = 0.0
      class_gat%TFXGAT (K) = 0.0
      class_gat%QEVPGAT(K) = 0.0
      class_gat%QFSGAT (K) = 0.0
      class_gat%QFXGAT (K) = 0.0
      class_gat%PETGAT (K) = 0.0
      class_gat%GAGAT  (K) = 0.0
      class_gat%EFGAT  (K) = 0.0
      class_gat%GTGAT  (K) = 0.0
      class_gat%QGGAT  (K) = 0.0
      class_gat%ALVSGAT(K) = 0.0
      class_gat%ALIRGAT(K) = 0.0
      class_gat%SFCTGAT(K) = 0.0
      class_gat%SFCUGAT(K) = 0.0
      class_gat%SFCVGAT(K) = 0.0
      class_gat%SFCQGAT(K) = 0.0
      class_gat%SFRHGAT(K) = 0.0
      class_gat%FSNOGAT(K) = 0.0
      class_gat%FSGVGAT(K) = 0.0
      class_gat%FSGSGAT(K) = 0.0
      class_gat%FSGGGAT(K) = 0.0
      class_gat%FLGVGAT(K) = 0.0
      class_gat%FLGSGAT(K) = 0.0
      class_gat%FLGGGAT(K) = 0.0
      class_gat%HFSCGAT(K) = 0.0
      class_gat%HFSSGAT(K) = 0.0
      class_gat%HFSGGAT(K) = 0.0
      class_gat%HEVCGAT(K) = 0.0
      class_gat%HEVSGAT(K) = 0.0
      class_gat%HEVGGAT(K) = 0.0
      class_gat%HMFCGAT(K) = 0.0
      class_gat%HMFNGAT(K) = 0.0
      class_gat%HTCCGAT(K) = 0.0
      class_gat%HTCSGAT(K) = 0.0
      class_gat%PCFCGAT(K) = 0.0
      class_gat%PCLCGAT(K) = 0.0
      class_gat%PCPNGAT(K) = 0.0
      class_gat%PCPGGAT(K) = 0.0
      class_gat%QFGGAT (K) = 0.0
      class_gat%QFNGAT (K) = 0.0
      class_gat%QFCFGAT(K) = 0.0
      class_gat%QFCLGAT(K) = 0.0
      class_gat%ROFGAT (K) = 0.0
      class_gat%ROFOGAT(K) = 0.0
      class_gat%ROFSGAT(K) = 0.0
      class_gat%ROFBGAT(K) = 0.0
      class_gat%TROFGAT(K) = 0.0
      class_gat%TROOGAT(K) = 0.0
      class_gat%TROSGAT(K) = 0.0
      class_gat%TROBGAT(K) = 0.0
      class_gat%ROFCGAT(K) = 0.0
      class_gat%ROFNGAT(K) = 0.0
      class_gat%ROVGGAT(K) = 0.0
      class_gat%WTRCGAT(K) = 0.0
      class_gat%WTRSGAT(K) = 0.0
      class_gat%WTRGGAT(K) = 0.0
      class_gat%DRGAT  (K) = 0.0
      class_gat%DSL    (K) = 0.0
      class_gat%DSLC   (K) = 0.0
    end do ! loop 330

    do L = 1,IGND ! loop 334
      do K = 1,nml ! loop 332
        class_gat%HMFGGAT(K,L) = 0.0
        class_gat%HTCGAT (K,L) = 0.0
        class_gat%QFCGAT (K,L) = 0.0
        class_gat%GFLXGAT(K,L) = 0.0
      end do ! loop 332
    end do ! loop 334

    do M = 1,50 ! loop 340
      do L = 1,6 ! loop 338
        do K = 1,NML ! loop 336
          class_gat%ITCTGAT(K,L,M) = 0
        end do ! loop 336
      end do ! loop 338
    end do ! loop 340

  end subroutine initDiagnosticVars
  !! @}

  !> \ingroup classstatevars_initRowVarsPhysics
  !! @{
  !> Initialization of diagnostic variables split out of classGather for consistency with gcm applications.
  subroutine initRowVarsPhysics()

    implicit none

    integer :: i,m

    ! These are nml variables.
    class_rot%CDHROW = 0.
    class_rot%CDMROW = 0.
    class_rot%HFSROW = 0.
    class_rot%TFXROW = 0.
    class_rot%QEVPROW = 0.
    class_rot%QFSROW = 0.
    class_rot%QFXROW = 0.
    class_rot%PETROW = 0.
    class_rot%GAROW = 0.
    class_rot%EFROW = 0.
    class_rot%GTROW = 0.
    class_rot%QGROW = 0.
    class_rot%ALVSROW = 0.
    class_rot%ALIRROW = 0.
    class_rot%SFCTROW = 0.
    class_rot%SFCUROW = 0.
    class_rot%SFCVROW = 0.
    class_rot%SFCQROW = 0.
    class_rot%SFRHROW = 0.
    class_rot%SNOROW = 0.
    class_rot%FSNOROW = 0.
    class_rot%FSGVROW = 0.
    class_rot%FSGSROW = 0.
    class_rot%FSGGROW = 0.
    class_rot%FLGVROW = 0.
    class_rot%FLGSROW = 0.
    class_rot%FLGGROW = 0.
    class_rot%HFSCROW = 0.
    class_rot%HFSSROW = 0.
    class_rot%HFSGROW = 0.
    class_rot%HEVCROW = 0.
    class_rot%HEVSROW = 0.
    class_rot%HEVGROW = 0.
    class_rot%HMFCROW = 0.
    class_rot%HMFNROW = 0.
    class_rot%HTCCROW = 0.
    class_rot%HTCSROW = 0.
    class_rot%PCFCROW = 0.
    class_rot%PCLCROW = 0.
    class_rot%PCPNROW = 0.
    class_rot%PCPGROW = 0.
    class_rot%QFGROW = 0.
    class_rot%QFNROW = 0.
    class_rot%QFCLROW = 0.
    class_rot%QFCFROW = 0.
    class_rot%ROFROW = 0.
    class_rot%ROFOROW = 0.
    class_rot%ROFSROW = 0.
    class_rot%ROFBROW = 0.
    class_rot%ROFCROW = 0.
    class_rot%ROFNROW = 0.
    class_rot%ROVGROW = 0.
    class_rot%RHOSROW = 0.
    class_rot%WTRCROW = 0.
    class_rot%WTRSROW = 0.
    class_rot%WTRGROW = 0.
    class_rot%DRROW = 0.
    class_rot%TCANROW = 0.
    class_rot%SCANROW = 0.
    class_rot%RCANROW = 0.
    class_rot%TSNOROW = 0.
    class_rot%WSNOROW = 0.
    class_rot%TPNDROW = 0.
    class_rot%ZPNDROW = 0.

    class_rot%ILMOROW = 0.
    class_rot%UEROW = 0.
    class_rot%HBLROW = 0.
    ! class_rot%groundHeatFluxROW = 0.

    ! These are nml,ignd
    class_rot%HMFGROW = 0.
    class_rot%HTCROW = 0.
    class_rot%QFCROW = 0.
    class_rot%GFLXROW = 0.
    class_rot%TBARROW = 0.
    class_rot%THALROW = 0.
    class_rot%THICROW = 0.
    class_rot%THLQROW = 0.

    ! Initialize to 0 for the start of a run.
    do m = 1,nmos
      do i = 1,nlat
        if (class_rot%FAREROT(i,m) > 0) then
          class_rot%actLyrThisYrROT(i,m) = 0.
        end if
      end do
    end do

  end subroutine initRowVarsPhysics
  !! @}

  !> \ingroup classstatevars_classdump
  !! @{
  !> dumps all class statevars for diagnostic purposes.
  subroutine classdump()

    implicit none

    print *, 'begin class dump'
    print *, 'ILMOS ', class_gat%ILMOS
    print *, 'JLMOS ', class_gat%JLMOS
    print *, 'IGDRGAT ', class_gat%IGDRGAT
    print *, 'PCSNGAT ', class_gat%PCSNGAT
    print *, 'CMAIGAT ', class_gat%CMAIGAT
    print *, 'QACGAT ', class_gat%QACGAT
    print *, 'RHOSGAT ', class_gat%RHOSGAT
    print *, 'SNOGAT ', class_gat%SNOGAT
    print *, 'TBASGAT ', class_gat%TBASGAT
    print *, 'TPNDGAT ', class_gat%TPNDGAT
    print *, 'WSNOGAT ', class_gat%WSNOGAT
    print *, 'ZPNDGAT ', class_gat%ZPNDGAT
    print *, 'BCSNGAT ', class_gat%BCSNGAT
    print *, 'AGVDGAT ', class_gat%AGVDGAT
    print *, 'ALGWGAT ', class_gat%ALGWGAT
    print *, 'ASVDGAT ', class_gat%ASVDGAT
    print *, 'GRKFGAT ', class_gat%GRKFGAT
    print *, 'WFSFGAT ', class_gat%WFSFGAT
    print *, 'ZPLGGAT ', class_gat%ZPLGGAT
    print *, 'ZSNLGAT ', class_gat%ZSNLGAT
    print *, 'ALGWNGAT ', class_gat%ALGWNGAT
    print *, 'ALGDNGAT ', class_gat%ALGDNGAT
    print *, 'CSZGAT ', class_gat%CSZGAT
    print *, 'DLATGAT ', class_gat%DLATGAT
    print *, 'FDLGAT ', class_gat%FDLGAT
    print *, 'FSVHGAT ', class_gat%FSVHGAT
    print *, 'PADRGAT ', class_gat%PADRGAT
    print *, 'PRESGAT ', class_gat%PRESGAT
    print *, 'RADJGAT ', class_gat%RADJGAT
    print *, 'RHSIGAT ', class_gat%RHSIGAT
    print *, 'SPCPGAT ', class_gat%SPCPGAT
    print *, 'TADPGAT ', class_gat%TADPGAT
    print *, 'TSPCGAT ', class_gat%TSPCGAT
    print *, 'VLGAT ', class_gat%VLGAT
    print *, 'VPDGAT ', class_gat%VPDGAT
    print *, 'ZBLDGAT ', class_gat%ZBLDGAT
    print *, 'ZDMGAT ', class_gat%ZDMGAT
    print *, 'ZRFMGAT ', class_gat%ZRFMGAT
    print *, 'FLGGAT ', class_gat%FLGGAT
    print *, 'DEPBGAT ', class_gat%DEPBGAT
    print *, 'SFCUBS ', class_gat%SFCUBS
    print *, 'USTARBS ', class_gat%USTARBS
    print *, 'GSNOW ', class_gat%GSNOW
    print *, 'ALVSGAT ', class_gat%ALVSGAT
    print *, 'CDMGAT ', class_gat%CDMGAT
    print *, 'EFGAT ', class_gat%EFGAT
    print *, 'FLGSGAT ', class_gat%FLGSGAT
    print *, 'FSGGGAT ', class_gat%FSGGGAT
    print *, 'FSGVGAT ', class_gat%FSGVGAT
    print *, 'GAGAT ', class_gat%GAGAT
    print *, 'HBLGAT ', class_gat%HBLGAT
    print *, 'HEVGGAT ', class_gat%HEVGGAT
    print *, 'HFSGAT ', class_gat%HFSGAT
    print *, 'HFSGGAT ', class_gat%HFSGGAT
    print *, 'HMFCGAT ', class_gat%HMFCGAT
    print *, 'HTCCGAT ', class_gat%HTCCGAT
    print *, 'ILMOGAT ', class_gat%ILMOGAT
    print *, 'PCLCGAT ', class_gat%PCLCGAT
    print *, 'PCPNGAT ', class_gat%PCPNGAT
    print *, 'QEVPGAT ', class_gat%QEVPGAT
    print *, 'QFCLGAT ', class_gat%QFCLGAT
    print *, 'QFNGAT ', class_gat%QFNGAT
    print *, 'QFXGAT ', class_gat%QFXGAT
    print *, 'ROFGAT ', class_gat%ROFGAT
    print *, 'ROFCGAT ', class_gat%ROFCGAT
    print *, 'ROFOGAT ', class_gat%ROFOGAT
    print *, 'ROVGGAT ', class_gat%ROVGGAT
    print *, 'SFCTGAT ', class_gat%SFCTGAT
    print *, 'SFCVGAT ', class_gat%SFCVGAT
    print *, 'TROBGAT ', class_gat%TROBGAT
    print *, 'TROOGAT ', class_gat%TROOGAT
    print *, 'UEGAT ', class_gat%UEGAT
    print *, 'WTRGGAT ', class_gat%WTRGGAT
    print *, 'wtableGAT ', class_gat%wtableGAT
    print *, 'SFRHGAT ', class_gat%SFRHGAT
    print *, 'FVAP ', class_gat%FVAP
    print *, 'FC ', class_gat%FC
    print *, 'FCS ', class_gat%FCS
    print *, 'RBCOEF ', class_gat%RBCOEF
    print *, 'FSVF ', class_gat%FSVF
    print *, 'ALVSCN ', class_gat%ALVSCN
    print *, 'ALVSG ', class_gat%ALVSG
    print *, 'ALVSCS ', class_gat%ALVSCS
    print *, 'ALVSSN ', class_gat%ALVSSN
    print *, 'ALVSGC ', class_gat%ALVSGC
    print *, 'ALVSSC ', class_gat%ALVSSC
    print *, 'TRVSCN ', class_gat%TRVSCN
    print *, 'TRVSCS ', class_gat%TRVSCS
    print *, 'RC ', class_gat%RC
    print *, 'FRAINC ', class_gat%FRAINC
    print *, 'FRAICS ', class_gat%FRAICS
    print *, 'CMASSC ', class_gat%CMASSC
    print *, 'DISP ', class_gat%DISP
    print *, 'ZOMLNC ', class_gat%ZOMLNC
    print *, 'ZOMLNG ', class_gat%ZOMLNG
    print *, 'ZOMLCS ', class_gat%ZOMLCS
    print *, 'ZOMLNS ', class_gat%ZOMLNS
    print *, 'TRSNOWC ', class_gat%TRSNOWC
    print *, 'CHCAPS ', class_gat%CHCAPS
    print *, 'GZEROG ', class_gat%GZEROG
    print *, 'GZROGS ', class_gat%GZROGS
    print *, 'G12C ', class_gat%G12C
    print *, 'G12CS ', class_gat%G12CS
    print *, 'G23C ', class_gat%G23C
    print *, 'G23CS ', class_gat%G23CS
    print *, 'QFREZC ', class_gat%QFREZC
    print *, 'QMELTC ', class_gat%QMELTC
    print *, 'EVAPC ', class_gat%EVAPC
    print *, 'EVAPG ', class_gat%EVAPG
    print *, 'EVPCSG ', class_gat%EVPCSG
    print *, 'TCANO ', class_gat%TCANO
    print *, 'RAICAN ', class_gat%RAICAN
    print *, 'RAICNS ', class_gat%RAICNS
    print *, 'CWLCAP ', class_gat%CWLCAP
    print *, 'CWLCPS ', class_gat%CWLCPS
    print *, 'TSNOCS ', class_gat%TSNOCS
    print *, 'RHOSCS ', class_gat%RHOSCS
    print *, 'WSNOCS ', class_gat%WSNOCS
    print *, 'TPONDC ', class_gat%TPONDC
    print *, 'TPNDCS ', class_gat%TPNDCS
    print *, 'ZPLMCS ', class_gat%ZPLMCS
    print *, 'ZPLIMC ', class_gat%ZPLIMC
    print *, 'CTVSTP ', class_gat%CTVSTP
    print *, 'CT1STP ', class_gat%CT1STP
    print *, 'CT3STP ', class_gat%CT3STP
    print *, 'WTSSTP ', class_gat%WTSSTP
    print *, 'DELZ ', class_gat%DELZ
    print *, 'ISNDGAT ', class_gat%ISNDGAT
    print *, 'THICGAT ', class_gat%THICGAT
    print *, 'THLQGAT ', class_gat%THLQGAT
    print *, 'BIGAT ', class_gat%BIGAT
    print *, 'GRKSGAT ', class_gat%GRKSGAT
    print *, 'PSISGAT ', class_gat%PSISGAT
    print *, 'TCSGAT ', class_gat%TCSGAT
    print *, 'THMGAT ', class_gat%THMGAT
    print *, 'THRGAT ', class_gat%THRGAT
    print *, 'ZBTWGAT ', class_gat%ZBTWGAT
    print *, 'GFLXGAT ', class_gat%GFLXGAT
    print *, 'HTCGAT ', class_gat%HTCGAT
    print *, 'TBARC ', class_gat%TBARC
    print *, 'TBARCS ', class_gat%TBARCS
    print *, 'THLIQC ', class_gat%THLIQC
    print *, 'THICEC ', class_gat%THICEC
    print *, 'FROOT ', class_gat%FROOT
    print *, 'HCPG ', class_gat%HCPG
    print *, 'TCTOPC ', class_gat%TCTOPC
    print *, 'TCTOPG ', class_gat%TCTOPG
    print *, 'ACIDGAT ', class_gat%ACIDGAT
    print *, 'CMASGAT ', class_gat%CMASGAT
    print *, 'PAIDGAT ', class_gat%PAIDGAT
    print *, 'PAMXGAT ', class_gat%PAMXGAT
    print *, 'PSGBGAT ', class_gat%PSGBGAT
    print *, 'ROOTGAT ', class_gat%ROOTGAT
    print *, 'VPDAGAT ', class_gat%VPDAGAT
    print *, 'ALICGAT ', class_gat%ALICGAT
    print *, 'FCANGAT ', class_gat%FCANGAT
    print *, 'FSDBGAT ', class_gat%FSDBGAT
    print *, 'FSSBGAT ', class_gat%FSSBGAT
    print *, 'CSALGAT ', class_gat%CSALGAT
    print *, 'ALSNO ', class_gat%ALSNO
    print *, 'TSFSGAT ', class_gat%TSFSGAT
    print *, 'ITCTGAT ', class_gat%ITCTGAT
    print *, 'CSZROW ', class_rot%CSZROW
    print *, 'DLATROW ', class_rot%DLATROW
    print *, 'latIndexROW ', class_rot%latIndexROW
    print *, 'RHOSROW ', class_rot%RHOSROW
    print *, 'PADRROW ', class_rot%PADRROW
    print *, 'PRESROW ', class_rot%PRESROW
    print *, 'RADJROW ', class_rot%RADJROW
    print *, 'RHSIROW ', class_rot%RHSIROW
    print *, 'RPREROW ', class_rot%RPREROW
    print *, 'SPREROW ', class_rot%SPREROW
    print *, 'TSNOROW ', class_rot%TSNOROW
    print *, 'TPNDROW ', class_rot%TPNDROW
    print *, 'SCANROW ', class_rot%SCANROW
    print *, 'TADPROW ', class_rot%TADPROW
    print *, 'TSPCROW ', class_rot%TSPCROW
    print *, 'VLROW ', class_rot%VLROW
    print *, 'VPDROW ', class_rot%VPDROW
    print *, 'ZDHROW ', class_rot%ZDHROW
    print *, 'ZRFHROW ', class_rot%ZRFHROW
    print *, 'UVROW ', class_rot%UVROW
    print *, 'PRENROW ', class_rot%PRENROW
    print *, 'DEPBROW ', class_rot%DEPBROW
    print *, 'BCSNROW ', class_rot%BCSNROW
    print *, 'ALVSROW ', class_rot%ALVSROW
    print *, 'CDMROW ', class_rot%CDMROW
    print *, 'EFROW ', class_rot%EFROW
    print *, 'FLGSROW ', class_rot%FLGSROW
    print *, 'FSGGROW ', class_rot%FSGGROW
    print *, 'FSGVROW ', class_rot%FSGVROW
    print *, 'GAROW ', class_rot%GAROW
    print *, 'HBLROW ', class_rot%HBLROW
    print *, 'HEVGROW ', class_rot%HEVGROW
    print *, 'HFSROW ', class_rot%HFSROW
    print *, 'HFSGROW ', class_rot%HFSGROW
    print *, 'HMFCROW ', class_rot%HMFCROW
    print *, 'HTCCROW ', class_rot%HTCCROW
    print *, 'ILMOROW ', class_rot%ILMOROW
    print *, 'PCLCROW ', class_rot%PCLCROW
    print *, 'PCPNROW ', class_rot%PCPNROW
    print *, 'QEVPROW ', class_rot%QEVPROW
    print *, 'QFCLROW ', class_rot%QFCLROW
    print *, 'QFNROW ', class_rot%QFNROW
    print *, 'QFXROW ', class_rot%QFXROW
    print *, 'ROFROW ', class_rot%ROFROW
    print *, 'ROFCROW ', class_rot%ROFCROW
    print *, 'ROFOROW ', class_rot%ROFOROW
    print *, 'ROVGROW ', class_rot%ROVGROW
    print *, 'SFCTROW ', class_rot%SFCTROW
    print *, 'SFCVROW ', class_rot%SFCVROW
    print *, 'UEROW ', class_rot%UEROW
    print *, 'WTRGROW ', class_rot%WTRGROW
    print *, 'SFRHROW ', class_rot%SFRHROW
    print *, 'WSNOROW ', class_rot%WSNOROW
    print *, 'FSSROW ', class_rot%FSSROW
    print *, 'fracFSFROW ', class_rot%fracFSFROW
    print *, 'FSVHROW ', class_rot%FSVHROW
    print *, 'FSGROL ', class_rot%FSGROL
    print *, 'XDIFFUS ', class_rot%XDIFFUS
    print *, 'ALBSROT ', class_rot%ALBSROT
    print *, 'GROROT ', class_rot%GROROT
    print *, 'RCANROT ', class_rot%RCANROT
    print *, 'SCANROT ', class_rot%SCANROT
    print *, 'TACROT ', class_rot%TACROT
    print *, 'TCANROT ', class_rot%TCANROT
    print *, 'TSNOROT ', class_rot%TSNOROT
    print *, 'ZPNDROT ', class_rot%ZPNDROT
    print *, 'AGVDROT ', class_rot%AGVDROT
    print *, 'ALGWROT ', class_rot%ALGWROT
    print *, 'ASVDROT ', class_rot%ASVDROT
    print *, 'FAREROT ', class_rot%FAREROT
    print *, 'WFCIROT ', class_rot%WFCIROT
    print *, 'XSLPROT ', class_rot%XSLPROT
    print *, 'ZPLSROT ', class_rot%ZPLSROT
    print *, 'ZSNOROT ', class_rot%ZSNOROT
    print *, 'ALGWNROT ', class_rot%ALGWNROT
    print *, 'ALGDNROT ', class_rot%ALGDNROT
    print *, 'ALIRROT ', class_rot%ALIRROT
    print *, 'CDHROT ', class_rot%CDHROT
    print *, 'DRROT ', class_rot%DRROT
    print *, 'FLGGROT ', class_rot%FLGGROT
    print *, 'FLGVROT ', class_rot%FLGVROT
    print *, 'FSGSROT ', class_rot%FSGSROT
    print *, 'FSNOROT ', class_rot%FSNOROT
    print *, 'GTROT ', class_rot%GTROT
    print *, 'HEVCROT ', class_rot%HEVCROT
    print *, 'HEVSROT ', class_rot%HEVSROT
    print *, 'HFSCROT ', class_rot%HFSCROT
    print *, 'HFSSROT ', class_rot%HFSSROT
    print *, 'HMFNROT ', class_rot%HMFNROT
    print *, 'SDEPROT ', class_rot%SDEPROT
    print *, 'HTCSROT ', class_rot%HTCSROT
    print *, 'PCFCROT ', class_rot%PCFCROT
    print *, 'PCPGROT ', class_rot%PCPGROT
    print *, 'PETROT ', class_rot%PETROT
    print *, 'QFCFROT ', class_rot%QFCFROT
    print *, 'QFGROT ', class_rot%QFGROT
    print *, 'QFSROT ', class_rot%QFSROT
    print *, 'QGROT ', class_rot%QGROT
    print *, 'ROFBROT ', class_rot%ROFBROT
    print *, 'ROFNROT ', class_rot%ROFNROT
    print *, 'ROFSROT ', class_rot%ROFSROT
    print *, 'SFCQROT ', class_rot%SFCQROT
    print *, 'SFCUROT ', class_rot%SFCUROT
    print *, 'TFXROT ', class_rot%TFXROT
    print *, 'TROFROT ', class_rot%TROFROT
    print *, 'TROSROT ', class_rot%TROSROT
    print *, 'WTRCROT ', class_rot%WTRCROT
    print *, 'WTRSROT ', class_rot%WTRSROT
    print *, 'wtableROT ', class_rot%wtableROT
    print *, 'ACTLYR ', class_rot%ACTLYR
    print *, 'actLyrThisYrROT ', class_rot%actLyrThisYrROT
    print *, 'ISNDROT ', class_rot%ISNDROT
    print *, 'THICROT ', class_rot%THICROT
    print *, 'BIROT ', class_rot%BIROT
    print *, 'GRKSROT ', class_rot%GRKSROT
    print *, 'SANDROT ', class_rot%SANDROT
    print *, 'ORGMROT ', class_rot%ORGMROT
    print *, 'PSIWROT ', class_rot%PSIWROT
    print *, 'THFCROT ', class_rot%THFCROT
    print *, 'THPROT ', class_rot%THPROT
    print *, 'THRAROT ', class_rot%THRAROT
    print *, 'THLWROT ', class_rot%THLWROT
    print *, 'HMFGROT ', class_rot%HMFGROT
    print *, 'QFCROT ', class_rot%QFCROT
    print *, 'THLQACC_M ', class_rot%THLQACC_M
    print *, 'ACVDROT ', class_rot%ACVDROT
    print *, 'HGTDROT ', class_rot%HGTDROT
    print *, 'PAMNROT ', class_rot%PAMNROT
    print *, 'PSGAROT ', class_rot%PSGAROT
    print *, 'QA50ROT ', class_rot%QA50ROT
    print *, 'RSMNROT ', class_rot%RSMNROT
    print *, 'VPDBROT ', class_rot%VPDBROT
    print *, 'ALVCROT ', class_rot%ALVCROT
    print *, 'LNZ0ROT ', class_rot%LNZ0ROT
    print *, 'CSALROT ', class_rot%CSALROT
    print *, 'FSFBROL ', class_rot%FSFBROL
    print *, 'TBARROW ', class_rot%TBARROW
    print *, 'THICROW ', class_rot%THICROW
    print *, 'GFLXROW ', class_rot%GFLXROW
    print *, 'HTCROW ', class_rot%HTCROW
    print *, 'ITCTROT ', class_rot%ITCTROT
    print *, 'PREACC_M ', class_rot%PREACC_M
    print *, 'QEVPACC_M ', class_rot%QEVPACC_M
    print *, 'HMFNACC_M ', class_rot%HMFNACC_M
    print *, 'ZSNACC_M ', class_rot%ZSNACC_M
    print *, 'OVRACC_M ', class_rot%OVRACC_M
    print *, 'ROFBACC_M ', class_rot%ROFBACC_M
    print *, 'ROFNACC_M ', class_rot%ROFNACC_M
    print *, 'ALVSACC_M ', class_rot%ALVSACC_M
    print *, 'RHOSACC_M ', class_rot%RHOSACC_M
    print *, 'WSNOACC_M ', class_rot%WSNOACC_M
    print *, 'RCANACC_M ', class_rot%RCANACC_M
    print *, 'ALTOTACC_M ', class_rot%ALTOTACC_M
    print *, 'GROACC_M ', class_rot%GROACC_M
    print *, 'FLINACC_M ', class_rot%FLINACC_M
    print *, 'UVACC_M ', class_rot%UVACC_M
    print *, 'QAACC_M ', class_rot%QAACC_M
    print *, 'GROUNDEVAP_M ', class_rot%GROUNDEVAP_M
    print *, 'TRANSPACC_M ', class_rot%TRANSPACC_M
    print *, 'ZPNDACC_M ', class_rot%ZPNDACC_M
    print *, 'ALIRACC_MO ', class_out%ALIRACC_MO
    print *, 'FSINACC_MO ', class_out%FSINACC_MO
    print *, 'HFSACC_MO ', class_out%HFSACC_MO
    print *, 'groundHeatFlux_MO ', class_out%groundHeatFlux_MO
    print *, 'SCAN_MO ', class_out%SCAN_MO
    print *, 'WSNOACC_MO ', class_out%WSNOACC_MO
    print *, 'OVRACC_MO ', class_out%OVRACC_MO
    print *, 'ROFACC_MO ', class_out%ROFACC_MO
    print *, 'EVAPACC_MO ', class_out%EVAPACC_MO
    print *, 'TAACC_MO ', class_out%TAACC_MO
    print *, 'TSURFACC_MO ', class_out%TSURFACC_MO
    print *, 'FTABLE_MO ', class_out%FTABLE_MO
    print *, 'ACTLYR_MAX_MO ', class_out%ACTLYR_MAX_MO
    print *, 'FTABLE_MAX_MO ', class_out%FTABLE_MAX_MO
    print *, 'MRSOL_MO ', class_out%MRSOL_MO
    print *, 'ALSNOACC_MO ', class_out%ALSNOACC_MO
    print *, 'CANOPYEVAP ', class_out%CANOPYEVAP
    print *, 'EVSPSBL_MO ', class_out%EVSPSBL_MO
    print *, 'TBARACC_MO ', class_out%TBARACC_MO
    print *, 'THICACC_MO ', class_out%THICACC_MO
    print *, 'TCANACC_MO ', class_out%TCANACC_MO
    print *, 'ALIRACC_YR ', class_out%ALIRACC_YR
    print *, 'FSINACC_YR ', class_out%FSINACC_YR
    print *, 'HFSACC_YR ', class_out%HFSACC_YR
    print *, 'ROFACC_YR ', class_out%ROFACC_YR
    print *, 'EVAPACC_YR ', class_out%EVAPACC_YR
    print *, 'TAACC_YR ', class_out%TAACC_YR
    print *, 'ACTLYR_MIN_YR ', class_out%ACTLYR_MIN_YR
    print *, 'FTABLE_YR ', class_out%FTABLE_YR
    print *, 'FTABLE_MAX_YR ', class_out%FTABLE_MAX_YR
    print *, 'altotcntr_yr ', class_out%altotcntr_yr
    print *, 'THICACC_YR ', class_out%THICACC_YR
    print *, 'end class dump'

  end subroutine classdump
  !! @}

  !> \namespace classstatevars
  !!
  !! Contains the physics variable type structures.
end module classStateVars
