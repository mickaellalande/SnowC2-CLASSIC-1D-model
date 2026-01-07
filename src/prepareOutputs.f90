!> \file
!> Central module that handles all preparation and writing to output files
module prepareOutputs

  ! J. Melton Mar 30 2015

  implicit none

  ! subroutines contained in this module:

  public  :: class_hh_w           ! Prepares and writes the CLASS (physics) half hourly file
  public  :: class_daily_aw       ! Accumulates and writes the CLASS (physics) daily outputs
  public  :: class_monthly_aw     ! Accumulates and writes the CLASS (physics) monthly outputs
  public  :: class_annual_aw      ! Accumulates and writes the CLASS (physics) annual outputs
  public  :: convertUnitsCTEM     ! Converts units prior to output for CTEM (biogeochemistry) variables.
  public  :: ctem_daily_aw        ! Accumulates and writes the CTEM (biogeochemistry) daily outputs
  public  :: ctem_monthly_aw      ! Accumulates and writes the CTEM (biogeochemistry) monthly outputs
  public  :: ctem_annual_aw       ! Accumulates and writes the CTEM (biogeochemistry) annual outputs


contains

  !> \ingroup prepareoutputs_class_halfhourly_aw
  !> @{
  !> Prepares and writes the CLASS (physics) half hourly file
  !! @author J. Melton
  subroutine class_hh_w (lonLocalIndex, latLocalIndex, nltest, nmtest, ncount, nday, iday, realyr)

    use classStateVars, only : class_rot, class_gat, initRowVarsPhysics
    use ctemStateVars, only : c_switch, vrot
    use classicParams, only : ignd, icc, SBC, TFREZ, convertkgC
    use outputManager, only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: ncount
    integer, intent(in) :: nday
    integer, intent(in) :: iday
    integer, intent(in) :: realyr

    ! local variables
    real, dimension(1) :: timeStamp
    real :: ALTOT           !< Broadband albedo [-] (temporary variable)
    real :: EVAPSUM         !< Total evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$ (temporary variable)
    real :: FLSTAR
    real :: FSSTAR
    real :: TCN, TPN, TSN, TSNBOT, TSURF, ZSN
    ! real :: DSLav
    real :: an_grd, rml_grd, totvegarea
    real, dimension(nltest,nmtest,ignd) :: PSI
    real, dimension(nltest,ignd) :: PSI_grid
    real :: temp
    real, dimension(icc) :: anveggrd, rmlveggrd, fcanctot 
    integer :: i, j, m
    
    ! Associate names with variables defined in derived types.
    associate( &
    ctem_on           => c_switch%ctem_on,              & !< logical: True if this run includes the biogeochemistry parameterizations (CTEM) 
    Ncycle_on         => c_switch%Ncycle_on,            & !< logical: True if this run includes Nitrogen Cycle processes 
    dopertileoutput   => c_switch%dopertileoutput,      & !< logical: Switch for making extra output files that are at the per tile level 
    doperpftoutput    => c_switch%doperpftoutput,       & !< logical: Switch for making extra output files that are at the per pft level 
    TSFSGAT => class_gat%TSFSGAT,                       & !< real, dimension(:,:) : Ground surface temperature over subarea [K] 
    FG => class_gat%FG,                                 & !< real, dimension(:) : Subarea fractional coverage of modelled area - bare ground [ ] 
    FC => class_gat%FC,                                 & !< real, dimension(:) : Subarea fractional coverage of modelled area - ground under canopy [ ] 
    FGS => class_gat%FGS,                               & !< real, dimension(:) : Subarea fractional coverage of modelled area - snow-covered bare ground [ ] 
    FCS => class_gat%FCS,                               & !< real, dimension(:) : Subarea fractional coverage of modelled area - snow-covered ground under canopy  [ ] 

    ailcgrow => vrot%ailcg,                             & !< real, dimension(:,:,:) : Green LAI for CTEM's pfts 
    anvegrow => vrot%anveg,                             & !< real, dimension(:,:,:) : Net photosynthesis rate for each pft \f$[kg C m^{-2} s^{-1} ]\f$ 
    rmlvegrow=> vrot%rmlveg,                            & !< real, dimension(:,:,:) : Leaf maintenance respiration rate for each pft \f$[kg C m^{-2} s^{-1} ]\f$ 
    ancsvegrow  => vrot%ancsveg,                        & !< real, dimension(:,:,:) : Net photosynthetic rate for CTEM's pfts for canopy over snow subarea \f$[\mu mol CO_2 m^{-2} s^{-1} ]\f$ 
    ancgvegrow  => vrot%ancgveg,                        & !< real, dimension(:,:,:) : Net photosynthetic rate for CTEM's pfts for canopy over ground subarea \f$[\mu mol CO_2 m^{-2} s^{-1} ]\f$ 
    rmlcsvegrow => vrot%rmlcsveg,                       & !< real, dimension(:,:,:) : Leaf respiration rate for CTEM' pfts forcanopy over snow subarea \f$[\mu mol CO_2 m^{-2} s^{-1} ]\f$ 
    rmlcgvegrow => vrot%rmlcgveg,                       & !< real, dimension(:,:,:) : Leaf respiration rate for CTEM' pfts forcanopy over ground subarea \f$[\mu mol CO_2 m^{-2} s^{-1} ]\f$ 

    FAREROT=> class_rot%FAREROT,                        & !< real, dimension(:,:) : Fractional coverage of mosaic tile on modelled area 
    CDHROT => class_rot%CDHROT,                         & !< real, dimension(:,:) : Surface drag coefficient for heat [ ] 
    CDMROT => class_rot%CDMROT,                         & !< real, dimension(:,:) : Surface drag coefficient for momentum [ ] 
    HFSROT => class_rot%HFSROT,                         & !< real, dimension(:,:) : Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$ 
    TFXROT => class_rot%TFXROT,                         & !< real, dimension(:,:) : Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$ 
    QEVPROT => class_rot%QEVPROT,                       & !< real, dimension(:,:) : Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$ 
    QFSROT => class_rot%QFSROT,                         & !< real, dimension(:,:) : Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    QFXROT => class_rot%QFXROT,                         & !< real, dimension(:,:) : Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$ 
    PETROT => class_rot%PETROT,                         & !< real, dimension(:,:) : Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$ 
    GAROT => class_rot%GAROT,                           & !< real, dimension(:,:) : Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$ 
    EFROT => class_rot%EFROT,                           & !< real, dimension(:,:) : Evaporation efficiency at ground surface [ ] 
    GTROT => class_rot%GTROT,                           & !< real, dimension(:,:) : Diagnosed effective surface black-body temperature [K] 
    QGROT => class_rot%QGROT,                           & !< real, dimension(:,:) : Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$ 
    ALIRROT => class_rot%ALIRROT,                       & !< real, dimension(:,:) : Diagnosed total near-infrared albedo of land surface [ ] 
    ALVSROT => class_rot%ALVSROT,                       & !< real, dimension(:,:) : Diagnosed total visible albedo of land surface [ ] 
    SFCQROT => class_rot%SFCQROT,                       & !< real, dimension(:,:) : Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$ 
    SFCTROT => class_rot%SFCTROT,                       & !< real, dimension(:,:) : Diagnosed screen-level air temperature [K] 
    SFCUROT => class_rot%SFCUROT,                       & !< real, dimension(:,:) : Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$ 
    SFCVROT => class_rot%SFCVROT,                       & !< real, dimension(:,:) : Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$ 
    SFRHROT => class_rot%SFRHROT,                       & !< real, dimension(:,:) : 
    FSNOROT => class_rot%FSNOROT,                       & !< real, dimension(:,:) : Diagnosed fractional snow coverage [ ] 
    FLGGROT => class_rot%FLGGROT,                       & !< real, dimension(:,:) : Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$ 
    FLGSROT => class_rot%FLGSROT,                       & !< real, dimension(:,:) : Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$ 
    FLGVROT => class_rot%FLGVROT,                       & !< real, dimension(:,:) : Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$ 
    FSGGROT => class_rot%FSGGROT,                       & !< real, dimension(:,:) : Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$ 
    FSGSROT => class_rot%FSGSROT,                       & !< real, dimension(:,:) : Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$ 
    FSGVROT => class_rot%FSGVROT,                       & !< real, dimension(:,:) : Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$ 
    HEVCROT => class_rot%HEVCROT,                       & !< real, dimension(:,:) : Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$ 
    HEVGROT => class_rot%HEVGROT,                       & !< real, dimension(:,:) : Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$ 
    HEVSROT => class_rot%HEVSROT,                       & !< real, dimension(:,:) : Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$ 
    HFSCROT => class_rot%HFSCROT,                       & !< real, dimension(:,:) : Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$ 
    HFSGROT => class_rot%HFSGROT,                       & !< real, dimension(:,:) : Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$ 
    HFSSROT => class_rot%HFSSROT,                       & !< real, dimension(:,:) : Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$ 
    HMFCROT => class_rot%HMFCROT,                       & !< real, dimension(:,:) : Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$ 
    HMFNROT => class_rot%HMFNROT,                       & !< real, dimension(:,:) : Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$ 
    HTCCROT => class_rot%HTCCROT,                       & !< real, dimension(:,:) : Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$ 
    HTCSROT => class_rot%HTCSROT,                       & !< real, dimension(:,:) : Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$ 
    PCFCROT => class_rot%PCFCROT,                       & !< real, dimension(:,:) : Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    PCLCROT => class_rot%PCLCROT,                       & !< real, dimension(:,:) : Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    PCPGROT => class_rot%PCPGROT,                       & !< real, dimension(:,:) : Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$ 
    PCPNROT => class_rot%PCPNROT,                       & !< real, dimension(:,:) : Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$ 
    QFCFROT => class_rot%QFCFROT,                       & !< real, dimension(:,:) : Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    QFCLROT => class_rot%QFCLROT,                       & !< real, dimension(:,:) : Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    QFGROT => class_rot%QFGROT,                         & !< real, dimension(:,:) : Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$ 
    QFNROT => class_rot%QFNROT,                         & !< real, dimension(:,:) : Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFROT => class_rot%ROFROT,                         & !< real, dimension(:,:) : Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFBROT => class_rot%ROFBROT,                       & !< real, dimension(:,:) : Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFCROT => class_rot%ROFCROT,                       & !< real, dimension(:,:) : Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFNROT => class_rot%ROFNROT,                       & !< real, dimension(:,:) : Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFOROT => class_rot%ROFOROT,                       & !< real, dimension(:,:) : Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFSROT => class_rot%ROFSROT,                       & !< real, dimension(:,:) : Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    ROVGROT => class_rot%ROVGROT,                       & !< real, dimension(:,:) : Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$ 
    WTRCROT => class_rot%WTRCROT,                       & !< real, dimension(:,:) : Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$ 
    WTRGROT => class_rot%WTRGROT,                       & !< real, dimension(:,:) : Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$ 
    WTRSROT => class_rot%WTRSROT,                       & !< real, dimension(:,:) : Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$ 
    DRROT => class_rot%DRROT,                           & !< real, dimension(:,:) : Surface drag coefficient under neutral stability [ ] 
    ILMOROT => class_rot%ILMOROT,                       & !< real, dimension(:,:) : Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$ 
    UEROT => class_rot%UEROT,                           & !< real, dimension(:,:) : Friction velocity of air \f$[m s^{-1} ]\f$ 
    HBLROT => class_rot%HBLROT,                         & !< real, dimension(:,:) : Height of the atmospheric boundary layer [m] 
    HMFGROT => class_rot%HMFGROT,                       & !< real, dimension(:,:,:) : Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$ 
    GFLXROT => class_rot%GFLXROT,                       & !< real, dimension(:,:,:) : Heat conduction between soil layers \f$[W m^{-2} ]\f$ 
    HTCROT => class_rot%HTCROT,                         & !< real, dimension(:,:,:) : Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$ 
    QFCROT => class_rot%QFCROT,                         & !< real, dimension(:,:,:) : Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$ 
    TBARROT=> class_rot%TBARROT,                        & !< real, dimension(:,:,:) : Temperature of soil layers [K] 
    TCTOROT=> class_rot%TCTOROT,                        & !< real, dimension(:,:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$
    TCBOROT=> class_rot%TCBOROT,                        & !< real, dimension(:,:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$
    THICROT=> class_rot%THICROT,                        & !< real, dimension(:,:,:) : Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    THLQROT=> class_rot%THLQROT,                        & !< real, dimension(:,:,:) : Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    RHOSROT => class_rot%RHOSROT,                       & !< real, dimension(:,:) : Density of snow \f$[kg m^{-3}]\f$ 
    SCANROT => class_rot%SCANROT,                       & !< real, dimension(:,:) : Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$ 
    RCANROT => class_rot%RCANROT,                       & !< real, dimension(:,:) : Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$ 
    SNOROT => class_rot%SNOROT,                         & !< real, dimension(:,:) : Mass of snow pack \f$[kg m^{-2}]\f$ 
    TCSNROT => class_rot%TCSNROT,                       & !< real, dimension(:,:) : Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$ 
    WSNOROT => class_rot%WSNOROT,                       & !< real, dimension(:,:) : Liquid water content of snow pack \f$[kg m^{-2} ]\f$ 
    TCANROT=> class_rot%TCANROT,                        & !< real, dimension(:,:) : Vegetation canopy temperature [K] 
    TSNOROT=> class_rot%TSNOROT,                        & !< real, dimension(:,:) : Snowpack temperature [K] 
    TSNBROT=> class_rot%TSNBROT,                        & !< real, dimension(:,:) : Bottom snowpack temperature [K] 
    TPNDROT=> class_rot%TPNDROT,                        & !< real, dimension(:,:) : Temperature of ponded water [K] 
    ZPNDROT=> class_rot%ZPNDROT,                        & !< real, dimension(:,:) : Depth of ponded water [m] 
    dlzwrot => class_rot%dlzwrot,                       & !< real, dimension(:,:,:) : Permeable thickness of soil layer [m] 
        
    PREROW  => class_rot%PREROW,                        & !< real, dimension(:): Surface precipitation rate \f$[kg m^{-2}  s^{-1} ]\f$
    UVROW => class_rot%UVROW,                           & !< real, dimension(:): Wind speed at reference height \f$[m s^{-1} ]\f$
    TAROW  => class_rot%TAROW,                          & !< real, dimension(:): Air temperature at reference height [K]
    QAROW => class_rot%QAROW,                           & !< real, dimension(:) : Specific humidity at reference height \f$[kg kg^{-1}]\f$
    PRESROW => class_rot%PRESROW,                       & !< real, dimension(:) : Surface air pressure [Pa]

    ! FLAG some of the row variables here can probably be moved to local vars. JM Nov 2017.
    FSSROW => class_rot%FSSROW,                         & !< real, dimension(:) : Shortwave radiation \f$[W m^{-2} ]\f$
    FDLROW => class_rot%FDLROW,                         & !< real, dimension(:) : Downwelling longwave sky radiation \f$[W m^{-2} ]\f$
    RHOSROW => class_rot%RHOSROW,                       & !< real, dimension(:) : Density of snow \f$[kg m^{-3}]\f$
    SNOROW => class_rot%SNOROW,                         & !< real, dimension(:) : Mass of snow pack \f$[kg m^{-2}]\f$
    TCSNROW => class_rot%TCSNROW,                       & !< real, dimension(:) : Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
    BCSNROW => class_rot%BCSNROW,                       & !< Black carbon mixing ratio \f$[kg m^{-3}]\f$
    CDHROW => class_rot%CDHROW,                         & !< real, dimension(:) : Surface drag coefficient for heat [ ]
    CDMROW => class_rot%CDMROW,                         & !< real, dimension(:) : Surface drag coefficient for momentum [ ]
    HFSROW => class_rot%HFSROW,                         & !< real, dimension(:) : Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    TFXROW => class_rot%TFXROW,                         & !< real, dimension(:) : Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    QEVPROW => class_rot%QEVPROW,                       & !< real, dimension(:) : Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    QFSROW => class_rot%QFSROW,                         & !< real, dimension(:) : Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    QFXROW => class_rot%QFXROW,                         & !< real, dimension(:) : Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    PETROW => class_rot%PETROW,                         & !< real, dimension(:) : Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    GAROW => class_rot%GAROW,                           & !< real, dimension(:) : Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    EFROW => class_rot%EFROW,                           & !< real, dimension(:) : Evaporation efficiency at ground surface [ ]
    GTROW => class_rot%GTROW,                           & !< real, dimension(:) : Diagnosed effective surface black-body temperature [K]
    QGROW => class_rot%QGROW,                           & !< real, dimension(:) : Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    ALIRROW => class_rot%ALIRROW,                       & !< real, dimension(:) : Diagnosed total near-infrared albedo of land surface [ ]
    ALVSROW => class_rot%ALVSROW,                       & !< real, dimension(:) : Diagnosed total visible albedo of land surface [ ]
    SFCQROW => class_rot%SFCQROW,                       & !< real, dimension(:) : Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    SFCTROW => class_rot%SFCTROW,                       & !< real, dimension(:) : Diagnosed screen-level air temperature [K]
    SFCUROW => class_rot%SFCUROW,                       & !< real, dimension(:) : Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    SFCVROW => class_rot%SFCVROW,                       & !< real, dimension(:) : Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    SFRHROW => class_rot%SFRHROW,                       & !< real, dimension(:) : Diagnosed screen-level relative humidity [%]
    FSNOROW => class_rot%FSNOROW,                       & !< real, dimension(:) : Diagnosed fractional snow coverage [ ]
    FLGGROW => class_rot%FLGGROW,                       & !< real, dimension(:) : Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    FLGSROW => class_rot%FLGSROW,                       & !< real, dimension(:) : Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    FLGVROW => class_rot%FLGVROW,                       & !< real, dimension(:) : Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    FSGGROW => class_rot%FSGGROW,                       & !< real, dimension(:) : Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    FSGSROW => class_rot%FSGSROW,                       & !< real, dimension(:) : Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    FSGVROW => class_rot%FSGVROW,                       & !< real, dimension(:) : Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    HEVCROW => class_rot%HEVCROW,                       & !< real, dimension(:) : Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    HEVGROW => class_rot%HEVGROW,                       & !< real, dimension(:) : Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    HEVSROW => class_rot%HEVSROW,                       & !< real, dimension(:) : Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    HFSCROW => class_rot%HFSCROW,                       & !< real, dimension(:) : Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    HFSGROW => class_rot%HFSGROW,                       & !< real, dimension(:) : Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    HFSSROW => class_rot%HFSSROW,                       & !< real, dimension(:) : Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    HMFCROW => class_rot%HMFCROW,                       & !< real, dimension(:) : Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    HMFNROW => class_rot%HMFNROW,                       & !< real, dimension(:) : Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    HTCCROW => class_rot%HTCCROW,                       & !< real, dimension(:) : Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    HTCSROW => class_rot%HTCSROW,                       & !< real, dimension(:) : Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    PCFCROW => class_rot%PCFCROW,                       & !< real, dimension(:) : Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    PCLCROW => class_rot%PCLCROW,                       & !< real, dimension(:) : Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    PCPGROW => class_rot%PCPGROW,                       & !< real, dimension(:) : Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    PCPNROW => class_rot%PCPNROW,                       & !< real, dimension(:) : Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    QFCFROW => class_rot%QFCFROW,                       & !< real, dimension(:) : Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    QFCLROW => class_rot%QFCLROW,                       & !< real, dimension(:) : Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    QFGROW => class_rot%QFGROW,                         & !< real, dimension(:) : Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    QFNROW => class_rot%QFNROW,                         & !< real, dimension(:) : Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    ROFROW => class_rot%ROFROW,                         & !< real, dimension(:) : Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    ROFBROW => class_rot%ROFBROW,                       & !< real, dimension(:) : Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROFCROW => class_rot%ROFCROW,                       & !< real, dimension(:) : Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    ROFNROW => class_rot%ROFNROW,                       & !< real, dimension(:) : Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    ROFOROW => class_rot%ROFOROW,                       & !< real, dimension(:) : Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROFSROW => class_rot%ROFSROW,                       & !< real, dimension(:) : Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROVGROW => class_rot%ROVGROW,                       & !< real, dimension(:) : Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    WTRCROW => class_rot%WTRCROW,                       & !< real, dimension(:) : Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    WTRGROW => class_rot%WTRGROW,                       & !< real, dimension(:) : Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    WTRSROW => class_rot%WTRSROW,                       & !< real, dimension(:) : Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    DRROW => class_rot%DRROW,                           & !< real, dimension(:) : Surface drag coefficient under neutral stability [ ]
    ILMOROW => class_rot%ILMOROW,                       & !< real, dimension(:) : Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    UEROW => class_rot%UEROW,                           & !< real, dimension(:) : Friction velocity of air \f$[m s^{-1} ]\f$
    HBLROW => class_rot%HBLROW,                         & !< real, dimension(:) : Height of the atmospheric boundary layer [m]
    HMFGROW => class_rot%HMFGROW,                       & !< real, dimension(:,:) : Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    GFLXROW => class_rot%GFLXROW,                       & !< real, dimension(:,:) : Heat conduction between soil layers \f$[W m^{-2} ]\f$
    HTCROW => class_rot%HTCROW,                         & !< real, dimension(:,:) : Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    QFCROW => class_rot%QFCROW,                         & !< real, dimension(:,:) : Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    TBARROW => class_rot%TBARROW,                       & !< real, dimension(:,:) : Temperature of soil layers [K]
    TCTOROW=> class_rot%TCTOROW,                        & !< real, dimension(:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$
    TCBOROW=> class_rot%TCBOROW,                        & !< real, dimension(:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$
    THALROW => class_rot%THALROW,                       & !< real, dimension(:,:) : Total volumetric water content of soil layers \f$[m^3 m^{-3} ]\f$
    THLQROW => class_rot%THLQROW,                       & !< real, dimension(:,:) : Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    THICROW => class_rot%THICROW,                       & !< real, dimension(:,:) : Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    TCANROW => class_rot%TCANROW,                       & !< real, dimension(:) : Vegetation canopy temperature [K]
    SCANROW => class_rot%SCANROW,                       & !< real, dimension(:) : Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    RCANROW => class_rot%RCANROW,                       & !< real, dimension(:) : Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    TSNOROW => class_rot%TSNOROW,                       & !< real, dimension(:) : Snowpack temperature [K]
    TSNBROW => class_rot%TSNBROW,                       & !< real, dimension(:) : Bottom snowpack temperature [K]
    WSNOROW => class_rot%WSNOROW,                       & !< real, dimension(:) : Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    TPNDROW => class_rot%TPNDROW,                       & !< real, dimension(:) : Temperature of ponded water [K]
    ZPNDROW => class_rot%ZPNDROW,                       & !< real, dimension(:) : Depth of ponded water [m]
    fcancmxrow => vrot%fcancmx,                          & !< real, dimension(:,:,:) :
    PSISROT => class_rot%PSISROT,                       & !< Soil moisture suction at saturation [m]
    THPROT  => class_rot%THPROT,                        & !< Pore volume in soil layer \f$[m^3 m^{-3} ]\f$\
    BIROT   => class_rot%BIROT                          & !< Clapp and Hornberger empirical “b” parameter [ ]
    )

    ! Calculate grid cell average diagnostic fields.

    ! First set all to zero
    call initRowVarsPhysics

    do I = 1,NLTEST
      do M = 1,NMTEST
        CDHROW(I) = CDHROW(I) + CDHROT(I,M) * FAREROT(I,M)
        CDMROW(I) = CDMROW(I) + CDMROT(I,M) * FAREROT(I,M)
        HFSROW(I) = HFSROW(I) + HFSROT(I,M) * FAREROT(I,M)
        TFXROW(I) = TFXROW(I) + TFXROT(I,M) * FAREROT(I,M)
        QEVPROW(I) = QEVPROW(I) + QEVPROT(I,M) * FAREROT(I,M)
        QFSROW(I) = QFSROW(I) + QFSROT(I,M) * FAREROT(I,M)
        QFXROW(I) = QFXROW(I) + QFXROT(I,M) * FAREROT(I,M)
        PETROW(I) = PETROW(I) + PETROT(I,M) * FAREROT(I,M)
        RHOSROW(I) = RHOSROW(I) + RHOSROT(I,M) * FAREROT(I,M)
        GAROW(I) = GAROW(I) + GAROT(I,M) * FAREROT(I,M)
        EFROW(I) = EFROW(I) + EFROT(I,M) * FAREROT(I,M)
        GTROW(I) = GTROW(I) + GTROT(I,M) * FAREROT(I,M)
        QGROW(I) = QGROW(I) + QGROT(I,M) * FAREROT(I,M)
        ALVSROW(I) = ALVSROW(I) + ALVSROT(I,M) * FAREROT(I,M)
        ALIRROW(I) = ALIRROW(I) + ALIRROT(I,M) * FAREROT(I,M)
        SFCTROW(I) = SFCTROW(I) + SFCTROT(I,M) * FAREROT(I,M)
        SFCUROW(I) = SFCUROW(I) + SFCUROT(I,M) * FAREROT(I,M)
        SFCVROW(I) = SFCVROW(I) + SFCVROT(I,M) * FAREROT(I,M)
        SFCQROW(I) = SFCQROW(I) + SFCQROT(I,M) * FAREROT(I,M)
        SFRHROW(I) = SFRHROW(I) + SFRHROT(I,M) * FAREROT(I,M)
        SNOROW(I)  = SNOROW(I) + SNOROT(I,M) * FAREROT(I,M)
        TCSNROW(I)  = TCSNROW(I) + TCSNROT(I,M) * FAREROT(I,M)
        TSNOROW(I)  = TSNOROW(I) + TSNOROT(I,M) * FAREROT(I,M)
        TSNBROW(I)  = TSNBROW(I) + TSNBROT(I,M) * FAREROT(I,M)
        TCANROW(I)  = TCANROW(I) + TCANROT(I,M) * FAREROT(I,M)
        SCANROW(I)  = SCANROW(I) + SCANROT(I,M) * FAREROT(I,M)
        WSNOROW(I)  = WSNOROW(I) + WSNOROT(I,M) * FAREROT(I,M)
        FSNOROW(I) = FSNOROW(I) + FSNOROT(I,M) * FAREROT(I,M)
        FSGVROW(I) = FSGVROW(I) + FSGVROT(I,M) * FAREROT(I,M)
        FSGSROW(I) = FSGSROW(I) + FSGSROT(I,M) * FAREROT(I,M)
        FSGGROW(I) = FSGGROW(I) + FSGGROT(I,M) * FAREROT(I,M)
        FLGVROW(I) = FLGVROW(I) + FLGVROT(I,M) * FAREROT(I,M)
        FLGSROW(I) = FLGSROW(I) + FLGSROT(I,M) * FAREROT(I,M)
        FLGGROW(I) = FLGGROW(I) + FLGGROT(I,M) * FAREROT(I,M)
        HFSCROW(I) = HFSCROW(I) + HFSCROT(I,M) * FAREROT(I,M)
        HFSSROW(I) = HFSSROW(I) + HFSSROT(I,M) * FAREROT(I,M)
        HFSGROW(I) = HFSGROW(I) + HFSGROT(I,M) * FAREROT(I,M)
        HEVCROW(I) = HEVCROW(I) + HEVCROT(I,M) * FAREROT(I,M)
        HEVSROW(I) = HEVSROW(I) + HEVSROT(I,M) * FAREROT(I,M)
        HEVGROW(I) = HEVGROW(I) + HEVGROT(I,M) * FAREROT(I,M)
        HMFCROW(I) = HMFCROW(I) + HMFCROT(I,M) * FAREROT(I,M)
        HMFNROW(I) = HMFNROW(I) + HMFNROT(I,M) * FAREROT(I,M)
        HTCCROW(I) = HTCCROW(I) + HTCCROT(I,M) * FAREROT(I,M)
        HTCSROW(I) = HTCSROW(I) + HTCSROT(I,M) * FAREROT(I,M)
        PCFCROW(I) = PCFCROW(I) + PCFCROT(I,M) * FAREROT(I,M)
        PCLCROW(I) = PCLCROW(I) + PCLCROT(I,M) * FAREROT(I,M)
        PCPNROW(I) = PCPNROW(I) + PCPNROT(I,M) * FAREROT(I,M)
        PCPGROW(I) = PCPGROW(I) + PCPGROT(I,M) * FAREROT(I,M)
        QFGROW(I) = QFGROW(I) + QFGROT(I,M) * FAREROT(I,M)
        QFNROW(I) = QFNROW(I) + QFNROT(I,M) * FAREROT(I,M)
        QFCLROW(I) = QFCLROW(I) + QFCLROT(I,M) * FAREROT(I,M)
        QFCFROW(I) = QFCFROW(I) + QFCFROT(I,M) * FAREROT(I,M)
        ROFROW(I) = ROFROW(I) + ROFROT(I,M) * FAREROT(I,M)
        ROFOROW(I) = ROFOROW(I) + ROFOROT(I,M) * FAREROT(I,M)
        ROFSROW(I) = ROFSROW(I) + ROFSROT(I,M) * FAREROT(I,M)
        ROFBROW(I) = ROFBROW(I) + ROFBROT(I,M) * FAREROT(I,M)
        ROFCROW(I) = ROFCROW(I) + ROFCROT(I,M) * FAREROT(I,M)
        ROFNROW(I) = ROFNROW(I) + ROFNROT(I,M) * FAREROT(I,M)
        ROVGROW(I) = ROVGROW(I) + ROVGROT(I,M) * FAREROT(I,M)
        WTRCROW(I) = WTRCROW(I) + WTRCROT(I,M) * FAREROT(I,M)
        WTRSROW(I) = WTRSROW(I) + WTRSROT(I,M) * FAREROT(I,M)
        WTRGROW(I) = WTRGROW(I) + WTRGROT(I,M) * FAREROT(I,M)
        DRROW(I) = DRROW(I) + DRROT(I,M) * FAREROT(I,M)
        ILMOROW(I) = ILMOROW(I) + ILMOROT(I,M) * FAREROT(I,M)
        UEROW(I) = UEROW(I) + UEROT(I,M) * FAREROT(I,M)
        HBLROW(I) = HBLROW(I) + HBLROT(I,M) * FAREROT(I,M)
        
        do J = 1,IGND
          HMFGROW(I,J) = HMFGROW(I,J) + HMFGROT(I,M,J) * FAREROT(I,M)
          HTCROW(I,J) = HTCROW(I,J) + HTCROT(I,M,J) * FAREROT(I,M)
          QFCROW(I,J) = QFCROW(I,J) + QFCROT(I,M,J) * FAREROT(I,M)
          GFLXROW(I,J) = GFLXROW(I,J) + GFLXROT(I,M,J) * FAREROT(I,M)
          TBARROW(i,j) = TBARROW(i,j) + TBARROT(i,m,j) * FAREROT(i,m)
          TCTOROW(i,j) = TCTOROW(i,j) + TCTOROT(i,m,j) * FAREROT(i,m)
          TCBOROW(i,j) = TCBOROW(i,j) + TCBOROT(i,m,j) * FAREROT(i,m)
          THLQROW(i,j) = THLQROW(i,j) + THLQROT(i,m,j) * FAREROT(i,m)
          THICROW(i,j) = THICROW(i,j) + THICROT(i,m,j) * FAREROT(i,m)

          ! Calculate the soil matric potential for outputting
          if (THICROT(i,m,j) < THPROT(i,m,j)) then
          	PSI(I,M,J) = MAX(PSISROT(I,M,J) * (THLQROT(I,M,J) / (THPROT(I,M,J) - THICROT(i,m,j))) ** &
                		(-BIROT(I,M,J)), PSISROT(I,M,J)) 		
          else
          	PSI(I,M,J) = 10000.0 ! value used in heterotrophic respiration
          end if  
          PSI_grid(I,J) = PSI(I,M,J) * FAREROT(i,m)

        end do ! loop 550
      end do ! loop 575
    end do ! loop 600


    ! Prepare the timestamp for this timestep.
    timeStamp = consecDays + ((real(ncount) - 1.) / real(nday))

    ! Now prepare and write out the grid averaged physics variables to output files
    do I = 1,NLTEST

      ! First write out the model meteorological forcing so they can be compared to inputs
      ! or to check on those generated by the dissagregation module.
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fss_hh' ,timeStamp,'rsds', [FSSROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fdl_hh' ,timeStamp,'rlds', [FDLROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'pre_hh' ,timeStamp,'pr',   [PREROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'uv_hh' ,timeStamp,'uvas',  [UVROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ta_hh' ,timeStamp,'tas',   [TAROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qa_hh' ,timeStamp,'huss',  [QAROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'pres_hh' ,timeStamp,'ps',  [PRESROW(I)])

      ALTOT = 0.0
      if (FSSROW(I) > 0.0) ALTOT = (FSSROW(I) - (FSGVROW(I) + FSGSROW(I) + FSGGROW(I)))/FSSROW(I)
      FSSTAR = FSSROW(I) * (1.0 - ALTOT)

      call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_hh' ,timeStamp,'rss', [FSSTAR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'altot_hh' ,timeStamp,'albs', [ALTOT])

      FLSTAR = FDLROW(I) - SBC * GTROW(I) ** 4
      call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_hh' ,timeStamp,'rls', [FLSTAR])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_hh'     ,timeStamp,'hfss', [HFSROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_hh'     ,timeStamp,'hfls', [QEVPROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'snm_hh' ,timeStamp,'snm', [HMFNROW(I)])

      ZSN = 0.0
      if (RHOSROW(I) > 0.0) ZSN = SNOROW(I)/RHOSROW(I)
      call writeOutput1D(lonLocalIndex,latLocalIndex,'snd_hh'   ,timeStamp,'snd', [ZSN])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fsno_hh'   ,timeStamp,'snc', [FSNOROW(I)])

      TCN = 0.0
      if (TCANROW(I) > 0.01) TCN = TCANROW(I) - TFREZ
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tcs_hh'   ,timeStamp,'tcs', [TCN])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'scan_hh'   ,timeStamp,'scanopy', [SCANROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rcan_hh'   ,timeStamp,'rcanopy', [RCANROW(I)])

      ! Find the surface temperature for the grid average over all subareas (1-4)
      TSURF = FCS(I) * TSFSGAT(I,1) + FGS(I) * TSFSGAT(I,2) + FC(I) * TSFSGAT(I,3) + FG(I) * TSFSGAT(I,4)
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurf_hh'  ,timeStamp,'ts', [TSURF])

      ! Output the surface temperature over each subareas
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurfcs_hh'  ,timeStamp,'tscs', [TSFSGAT(I,1)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurfgs_hh'  ,timeStamp,'tsgs', [TSFSGAT(I,2)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurfc_hh'  ,timeStamp,'tsc', [TSFSGAT(I,3)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurfg_hh'  ,timeStamp,'tsg', [TSFSGAT(I,4)])

      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurf_subareas_hh'  ,timeStamp,'ts_subareas', [TSFSGAT(I,:)])

      ! Output the subarea fractional coverage of modelled area
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fcs_hh'  ,timeStamp,'fcs', [FCS(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fgs_hh'  ,timeStamp,'fgs', [FGS(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fc_hh'  ,timeStamp,'fc', [FC(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fg_hh'  ,timeStamp,'fg', [FG(I)])





      TSN = 0.0
      if (TSNOROW(I) > 0.01) TSN = TSNOROW(I) - TFREZ
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsno_hh'   ,timeStamp,'tsn', [TSN])
      TSNBOT = 0.0
      if (TSNBROW(I) > 0.01) TSNBOT = TSNBROW(I) - TFREZ
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsnb_hh'   ,timeStamp,'tsnbot', [TSNBOT])

      TPN = 0.0
      if (TPNDROW(I) > 0.01) TPN = TPNDROW(I) - TFREZ
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tpond_hh'   ,timeStamp,'tpond', [TPN])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'zpond_hh'   ,timeStamp,'zpond', [ZPNDROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gt_hh'   ,timeStamp,'tsblack', [GTROW(I) - TFREZ])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'sno_hh' ,timeStamp,'snw', [SNOROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'wsno_hh',timeStamp,'wsnw', [WSNOROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'snodens_hh',timeStamp,'snwdens', [RHOSROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tcsnow_hh',timeStamp,'tcsnow', [TCSNROW(I)])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'rof_hh',timeStamp,'mrro', [ROFROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofo_hh',timeStamp,'mrros', [ROFOROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofs_hh',timeStamp,'mrroi', [ROFSROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofb_hh',timeStamp,'mrrob', [ROFBROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofn_hh',timeStamp,'mrron', [ROFNROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap_hh',timeStamp,'evspsblveg', [QFCFROW(I)+QFCLROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap_hh',timeStamp,'evspsblsoi', [QFGROW(I)+QFNROW(I)])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'cdh_hh',timeStamp,'cdh', [CDHROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'cdm_hh',timeStamp,'cdm', [CDMROW(I)])
      !             call writeOutput1D(lonLocalIndex,latLocalIndex,'windu_hh',timeStamp,'windu', [SFCUROW(I)])  ! name
      !             call writeOutput1D(lonLocalIndex,latLocalIndex,'windv_hh',timeStamp,'windv', [SFCVROW(I)])  ! name
      !
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tbar_hh',timeStamp,'tsl', [TBARROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tctop_hh',timeStamp,'tctop', [TCTOROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tcbot_hh',timeStamp,'tcbot', [TCBOROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'thlq_hh',timeStamp,'mrsll', [THLQROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'thic_hh',timeStamp,'mrsfl', [THICROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gflx_hh',timeStamp,'gflx', [GFLXROW(I,:)])  
      call writeOutput1D(lonLocalIndex,latLocalIndex,'psi_hh',timeStamp,'psi', [PSI_grid(I,:)]) 

      ! If requested also prepare and write out the per tile physics variables
      if (dopertileoutput) then
        do m = 1,nmtest

          ALTOT = 0.0
          if (FSSROW(I) > 0.0) then ! there is no ROT for fss,as it will always be the same for all tiles.
            ALTOT = (FSSROW(I) - (FSGVROT(I,M) + FSGSROT(I,M) + FSGGROT(I,M))) / FSSROW(I)
          end if
          FSSTAR = FSSROW(I) * (1.0 - ALTOT)
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_hh_t' ,timeStamp,'rss', [FSSTAR])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'altot_hh_t' ,timeStamp,'albs', [ALTOT])

          FLSTAR = FDLROW(I) - SBC * GTROT(I,M) ** 4 ! there is no ROT for fdl,as it will always be the same for all tiles.
          call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_hh_t' ,timeStamp,'rls', [FLSTAR])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_hh_t'     ,timeStamp,'hfss', [HFSROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_hh_t'     ,timeStamp,'hfls', [QEVPROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'snm_hh_t'   ,timeStamp,'snm', [HMFNROT(I,M)])

          ZSN = 0.0
          if (RHOSROT(I,M) > 0.0) ZSN = SNOROT(I,M) / RHOSROT(I,M)
          call writeOutput1D(lonLocalIndex,latLocalIndex,'snd_hh_t'   ,timeStamp,'snd', [ZSN])

          TCN = 0.0
          if (TCANROT(I,M) > 0.01) TCN = TCANROT(I,M) - TFREZ
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tcs_hh_t'   ,timeStamp,'tcs', [TCN])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'scan_hh_t'   ,timeStamp,'scanopy', [SCANROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rcan_hh_t'   ,timeStamp,'rcanopy', [RCANROT(I,M)])

          ! TACGAT(I)-TFREZ  ! temp of air within canopy.
          TSN = 0.0
          if (TSNOROT(I,M) > 0.01) TSN = TSNOROT(I,M) - TFREZ
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tsno_hh_t'   ,timeStamp,'tsn', [TSN])
          TSNBOT = 0.0
          if (TSNBROT(I,M) > 0.01) TSNBOT = TSNBROT(I,M) - TFREZ
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tsnb_hh_t'   ,timeStamp,'tsnbot', [TSNBOT])

          TPN = 0.0
          if (TPNDROT(I,M) > 0.01) TPN = TPNDROT(I,M) - TFREZ
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tpond_hh_t'   ,timeStamp,'tpond', [TPN])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'zpond_hh_t'   ,timeStamp,'zpond', [ZPNDROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gt_hh_t'   ,timeStamp,'tsblack', [GTROT(I,M) - TFREZ])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'sno_hh_t' ,timeStamp,'snw', [SNOROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_hh_t',timeStamp,'wsnw', [WSNOROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'snodens_hh_t',timeStamp,'snwdens', [RHOSROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tcsnow_hh_t',timeStamp,'tcsnow', [TCSNROT(I,M)])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'rof_hh_t',timeStamp,'mrro', [ROFROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rofo_hh_t',timeStamp,'mrros', [ROFOROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rofs_hh_t',timeStamp,'mrroi', [ROFSROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rofb_hh_t',timeStamp,'mrrob', [ROFBROT(I,M)])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'cdh_hh_t',timeStamp,'cdh', [CDHROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'cdm_hh_t',timeStamp,'cdm', [CDMROT(I,M)])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'windu_hh_t',timeStamp,'windu', [SFCUROT(I,M)])  ! name
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'windv_hh_t',timeStamp,'windv', [SFCVROT(I,M)])  ! name

          do J = 1,IGND
            EVAPSUM = QFCFROT(I,M) + QFCLROT(I,M) + QFNROT(I,M) + QFGROT(I,M) + QFCROT(I,M,J)
          end do
          call writeOutput1D(lonLocalIndex,latLocalIndex,'evspsbl_hh_t'   ,timeStamp,'evspsbl', [EVAPSUM])

          !call writeOutput1D(lonLocalIndex,latLocalIndex,'tbar_hh_t',timeStamp,'tsl', [TBARROT(I,M,:)])
          ! Convert 64-bit precision TBARROT to default real type in writeOutput1D.
          ! FLAG, EC: This means TBARROT will be output as 32-bit reals if compiled in 32-bit mode. 
          !           Need to change if want to always output 64-bit reals instead.
          !call writeOutput1D(lonLocalIndex,latLocalIndex,'tbar_hh_t',timeStamp,'tsl', [real(TBARROT(I,M,:))])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tbar_hh_t',timeStamp,'tsl', real(TBARROT(I,M,:)))
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tctop_hh_t',timeStamp,'tctop', real(TCTOROT(I,M,:)))
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tcbot_hh_t',timeStamp,'tcbot', real(TCBOROT(I,M,:)))

          ! For mrsll and mrsfl, add in conversion from m3/m3 to kg/m2
          call writeOutput1D(lonLocalIndex,latLocalIndex,'thlq_hh_t',timeStamp,'mrsll', [THLQROT(I,M,:) * 1000. * DLZWROT(I,M,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'thic_hh_t',timeStamp,'mrsfl', [THICROT(I,M,:) * 1000. * DLZWROT(I,M,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gflx_hh_t',timeStamp,'gflx', [GFLXROT(I,M,:)])

        end do
      end if

      ! Other ones to make outputs with later.
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap',timeStamp,'evspsblsoi', [GROUNDEVAP(I)])
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap',timeStamp,'evspsblveg', [CANOPYEVAP(I)])
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_mo',timeStamp,'evspsbl', [EVAPACC_MO(I)])
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_mo',timeStamp,'tran', [TRANSPACC_MO(I)])
      ! TROFROT(I,M) !< Temperature of total runoff [K]
      ! TROOROT(I,M) !< Temperature of overland flow from top of soil column [K]
      ! TROSROT(I,M) !< Temperature of interflow from sides of soil column [K]
      ! TROBROT(I,M) !< Temperature of base flow from bottom of soil column [K]

      !                                 !
      !                                     &                   FCS(M),FGS(M),FC(M),FG(M),'
      !                                 WRITE(68,6800) IHOUR,IMIN,IDAY,realyr, &
      !                                     &                   FSGVROT(I,M),FSGSROT(I,M),FSGGROT(I,M), &
      !                                     &                   FLGVROT(I,M),FLGSROT(I,M),FLGGROT(I,M), &
      !                                     &                   HFSCROT(I,M),HFSSROT(I,M),HFSGROT(I,M), &
      !                                     &                   HEVCROT(I,M),HEVSROT(I,M),HEVGROT(I,M), &
      !                                     &                   HMFCROT(I,M),HMFNROT(I,M), &
      !                                     &                   (HMFGROT(I,M,J),J=1,3), &
      !                                     &                   HTCCROT(I,M),HTCSROT(I,M), &
      !                                     &                   (HTCROT(I,M,J),J=1,3),' TILE ',M
      !                                 WRITE(69,6900) IHOUR,IMIN,IDAY,realyr, &
      !                                     &                   PCFCROT(I,M),PCLCROT(I,M),PCPNROT(I,M), &
      !                                     &                   PCPGROT(I,M),QFCFROT(I,M),QFCLROT(I,M), &
      !                                     &                   QFNROT(I,M),QFGROT(I,M),(QFCROT(I,M,J),J=1,3), &
      !                                     &                   ROFCROT(I,M),ROFNROT(I,M),
      !                                     &                   ROFROT(I,M),WTRCROT(I,M),WTRSROT(I,M), &
      !                                     &                   WTRGROT(I,M),' TILE ',M

      ! Write half-hourly CTEM results to file
      !
      ! Net photosynthetic rates (GPP) and leaf maintenance respiration for each pft. however, if ctem_on then physyn subroutine
      ! is using storage lai while actual lai is zero. if actual lai is zero then we make anveg and rmlveg zero as well because these
      ! are imaginary just like storage lai. note that anveg and rmlveg are not passed to ctem. rather ancsveg, ancgveg, rmlcsveg, and
      ! rmlcgveg are passed.
      !
      if (ctem_on) then
        do m = 1,nmtest
          do j = 1,icc
            if (ailcgrow(i,m,j) <= 0.0) then
              anvegrow(i,m,j) = 0.0
              rmlvegrow(i,m,j) = 0.0
            else
              ! Add up the snow covered and non fluxes. Also convert from umol CO2/m2/s to kgC/m2/s.
              temp = (ancsvegrow(i,m,j) * FSNOROT(i,m) + ancgvegrow(i,m,j) * (1. - FSNOROT(i,m))) * convertkgC
              rmlvegrow(i,m,j) = (rmlcsvegrow(i,m,j) * FSNOROT(i,m) + rmlcgvegrow(i,m,j) * (1. - FSNOROT(i,m))) * convertkgC
              anvegrow(i,m,j) = temp + rmlvegrow(i,m,j) ! Add back in the rmLeaf to make it gross primary productivity.
            end if
          end do
          if (dopertileoutput) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_hh_t',timeStamp,'gpp', [anvegrow(I,M,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rml_hh_t',timeStamp,'rmLeaf', [rmlvegrow(I,M,:)])
          end if
          anveggrd = 0.0
          rmlveggrd = 0.0
          fcanctot = 0.0
          do j = 1,icc
            anveggrd(j) = anveggrd(j) + anvegrow(i,m,j) * FAREROT(i,m)
            rmlveggrd(j) = rmlveggrd(j) + rmlvegrow(i,m,j) * FAREROT(i,m)
            fcanctot(j) = fcanctot(j) + fcancmxrow(i,m,j) * FAREROT(i,m)
          end do
        end do ! m loop
        if (doperpftoutput) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_hh',timeStamp,'gpp', [anveggrd(:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rml_hh',timeStamp,'rmLeaf', [rmlveggrd(:)])
        end if
        ! Grid cell level outputs
        an_grd = 0.
        rml_grd = 0.
        totvegarea = sum(fcanctot)
        if (totvegarea > 0.) then
          do j = 1,icc
            an_grd  = an_grd + anveggrd(j)  * fcanctot(j) / totvegarea
            rml_grd = rml_grd + rmlveggrd(j) * fcanctot(j) / totvegarea          
          end do
        end if   
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_hh_g',timeStamp,'gpp', [an_grd])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rml_hh_g',timeStamp,'rmLeaf', [rml_grd])
      end if ! ctem_on
    end do ! nltest loop

    end associate
  end subroutine class_hh_w
  !< @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_class_daily_aw
  !> @{
  !> Accumlates and writes the daily physics variables. These are kept in pointer structures as
  !! this subroutine is called each physics timestep and we increment the timestep values to produce a daily value.
  !! The pointer to the daily data structures (in classStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine class_daily_aw (lonLocalIndex, latLocalIndex, iday, nltest, nmtest, ncount, nday, lastDOY, realyr)

    use classStateVars, only : class_rot, class_gat, resetAccVars
    use classicParams,  only : ignd, sbc, delt, TFREZ
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: iday
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: ncount
    integer, intent(in) :: nday
    integer, intent(in) :: lastDOY
    integer, intent(in) :: realyr

    ! local variables
    integer :: i, m, j, k
    real, dimension(1) :: timeStamp
    real :: FSSTAR, FLSTAR, TCN
    real, dimension(nltest) :: ALIRACC !< Diagnosed total near-infrared albedo of land surface [ ]
    real, dimension(nltest) :: ALVSACC !< Diagnosed total visible albedo of land surface [ ]
    real, dimension(nltest) :: EVAPACC !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: GROUNDEVAP  !< Surface water vapour flux over ground area \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: CANOPYEVAP  !< Water vapour flux from the canopy \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: TRANSPACC !<transpiration over modelled area \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: FLINACC !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nltest) :: FLUTACC !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, dimension(nltest) :: FSINACC !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, dimension(nltest) :: GROACC  !< Vegetation growth index [ ]
    real, dimension(nltest) :: GTACC   !< Diagnosed effective surface black-body temperature [K]
    real, dimension(nltest) :: HFSACC  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nltest) :: HMFNACC !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, dimension(nltest) :: OVRACC  !< Overland flow from top of soil column \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: PREACC  !< Surface precipitation rate \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: PRESACC !< Surface air pressure [Pa]
    real, dimension(nltest) :: QAACC   !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, dimension(nltest) :: QEVPACC !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, dimension(nltest) :: RCANACC !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: RHOSACC !< Density of snow \f$[kg m^{-3} ]\f$
    real, dimension(nltest) :: ROFACC  !< Total runoff from soil \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: ROFSACC !< Interflow from sides of the soil column \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: ROFBACC !< Base flow from bottom of soil column \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: ROFCACC !< liquid/frozen water runoff from vegetation \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: ROFNACC !< liquid water runoff from snowpack \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: SCANACC !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: ZSNACC  !< Depth of snow pack \f$[ m ]\f$
    real, dimension(nltest) :: SNOACC  !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: TCSNACC !< Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
    real, dimension(nltest) :: BCSNACC !< Black carbon mixing ratio \f$[kg m^{-3}]\f$
    real, dimension(nltest) :: FSNOACC !< Fractional cover of snow pack \f$[fraction]\f$
    real, dimension(nltest) :: TAACC   !< Air temperature at reference height [K]
    real, dimension(nltest) :: TCANACC !< Vegetation canopy temperature [K] (accumulated)
    real, dimension(nltest) :: TSURFACC !< Ground surface temperature [K] (accumulated)
    real, dimension(nltest) :: TSNOACC !< Snowpack temperature [K]
    real, dimension(nltest) :: TSNBACC !< Bottom snowpack temperature [K]
    real, dimension(nltest) :: UVACC   !< Wind speed \f$[m s^{-1} ]\f$
    real, dimension(nltest) :: WSNOACC !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, dimension(nltest) :: wtableACC !< Depth of water table in soil [m]
    real, dimension(nltest) :: ALTOTACC !< Broadband albedo [-]
    real, dimension(nltest) :: ALSNOACC !< Snow albedo [-]
    real, dimension(nltest,ignd) :: TBARACC  !< Temperature of soil layers [K] (accumulated)
    real, dimension(nltest,ignd) :: TCTOACC  !< Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$ (accumulated)
    real, dimension(nltest,ignd) :: TCBOACC  !< Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$ (accumulated)
    real, dimension(nltest,ignd) :: THLQACC  !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated)
    real, dimension(nltest,ignd) :: THICACC  !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated)
    real, dimension(nltest) :: ZPNDACC !<Depth of ponded water [m]
    real, dimension(nltest) :: ACTLYRACC !< Active layer depth [m]

    !!real, pointer, dimension(:,:) :: UVACC_M
    !!real, pointer, dimension(:,:) :: PRESACC_M
    !!real, pointer, dimension(:,:) :: QAACC_M

    ! Associate names with variables defined in derived types.
    associate( &
    TSFSGAT => class_gat%TSFSGAT,                       & !< real, dimension(:,:) : Ground surface temperature over subarea [K] 
    FG => class_gat%FG,                                 & !< real, dimension(:) : Subarea fractional coverage of modelled area - bare ground [ ] 
    FC => class_gat%FC,                                 & !< real, dimension(:) : Subarea fractional coverage of modelled area - ground under canopy [ ] 
    FGS => class_gat%FGS,                               & !< real, dimension(:) : Subarea fractional coverage of modelled area - snow-covered bare ground [ ] 
    FCS => class_gat%FCS,                               & !< real, dimension(:) : Subarea fractional coverage of modelled area - snow-covered ground under canopy  [ ] 
    FSSROW => class_rot%FSSROW,                         & !< real, dimension(:) : Shortwave radiation \f$[W m^{-2} ]\f$ 
    PREROW  => class_rot%PREROW,                        & !< real, dimension(:): Surface precipitation rate \f$[kg m^{-2}  s^{-1} ]\f$ 
    FSVHROW => class_rot%FSVHROW,                       & !< real, dimension(:): Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$ 
    FSIHROW => class_rot%FSIHROW,                       & !< real, dimension(:): Near infrared shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$ 
    FDLROW => class_rot%FDLROW,                         & !< real, dimension(:) : Downwelling longwave sky radiation \f$[W m^{-2} ]\f$ 
    TAROW  => class_rot%TAROW,                          & !< real, dimension(:): Air temperature at reference height [K] 
    TBARROT=> class_rot%TBARROT,                        & !< real, dimension(:,:,:) : Temperature of soil layers [K] 
    TCTOROT=> class_rot%TCTOROT,                        & !< real, dimension(:,:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$
    TCBOROT=> class_rot%TCBOROT,                        & !< real, dimension(:,:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$
    THICROT=> class_rot%THICROT,                        & !< real, dimension(:,:,:) : Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    THLQROT=> class_rot%THLQROT,                        & !< real, dimension(:,:,:) : Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    SNOROT => class_rot%SNOROT,                         & !< real, dimension(:,:) : Mass of snow pack \f$[kg m^{-2}]\f$ 
    BCSNROW => class_rot%BCSNROW,                       & !< Black carbon mixing ratio \f$[kg m^{-3}]\f$ 
    FSNOROT => class_rot%FSNOROT,                       & !< real, dimension(:,:) : Diagnosed fractional snow coverage [ ]
    HFSROT => class_rot%HFSROT,                         & !< real, dimension(:,:) : Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$ 
    QEVPROT => class_rot%QEVPROT,                       & !< real, dimension(:,:) : Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$ 
    QFSROT => class_rot%QFSROT,                         & !< real, dimension(:,:) : Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    QFCROT => class_rot%QFCROT,                         & !< real, dimension(:,:,:) : Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    QFCFROT => class_rot%QFCFROT,                       & !< real, dimension(:,:) : Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    QFNROT => class_rot%QFNROT,                         & !< real, dimension(:,:) : Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    QFCLROT => class_rot%QFCLROT,                       & !< real, dimension(:,:) : Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    QFGROT => class_rot%QFGROT,                         & !< real, dimension(:,:) : Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    GTROT => class_rot%GTROT,                           & !< real, dimension(:,:) : Diagnosed effective surface black-body temperature [K] 
    ALIRROT => class_rot%ALIRROT,                       & !< real, dimension(:,:) : Diagnosed total near-infrared albedo of land surface [ ] 
    ALVSROT => class_rot%ALVSROT,                       & !< real, dimension(:,:) : Diagnosed total visible albedo of land surface [ ] 
    ALBSROT => class_rot%ALBSROT,                       & !< real, dimension(:,:) : Snow albedo [ ] 
    FSGGROT => class_rot%FSGGROT,                       & !< real, dimension(:,:) : Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$ 
    FSGSROT => class_rot%FSGSROT,                       & !< real, dimension(:,:) : Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$ 
    FSGVROT => class_rot%FSGVROT,                       & !< real, dimension(:,:) : Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$ 
    HMFNROT => class_rot%HMFNROT,                       & !< real, dimension(:,:) : Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$ 
    ROFROT => class_rot%ROFROT,                         & !< real, dimension(:,:) : Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$ 
    GROROT => class_rot%GROROT,                         & !< real, dimension(:,:) : Vegetation growth index [ ] 
    ROFOROT => class_rot%ROFOROT,                       & !< real, dimension(:,:) : Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFBROT => class_rot%ROFBROT,                       & !< real, dimension(:,:) : Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROFCROT => class_rot%ROFCROT,                       & !< real, dimension(:,:) : Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    ROFNROT => class_rot%ROFNROT,                       & !< real, dimension(:,:) : Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    ROFSROT => class_rot%ROFSROT,                       & !< real, dimension(:,:) : Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    RHOSROT => class_rot%RHOSROT,                       & !< real, dimension(:,:) : Density of snow \f$[kg m^{-3}]\f$ 
    TCSNROT => class_rot%TCSNROT,                       & !< real, dimension(:,:) : Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$ 
    TSNOROT=> class_rot%TSNOROT,                        & !< real, dimension(:,:) : Snowpack temperature [K] 
    TSNBROT=> class_rot%TSNBROT,                        & !< real, dimension(:,:) : Bottom snowpack temperature [K] 
    WSNOROT => class_rot%WSNOROT,                       & !< real, dimension(:,:) : Liquid water content of snow pack \f$[kg m^{-2} ]\f$ 
    TCANROT=> class_rot%TCANROT,                        & !< real, dimension(:,:) : Vegetation canopy temperature [K] 
    RCANROT => class_rot%RCANROT,                       & !< real, dimension(:,:) : Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$ 
    SCANROT => class_rot%SCANROT,                       & !< real, dimension(:,:) : Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$ 
    FAREROT=> class_rot%FAREROT,                        & !< real, dimension(:,:) : Fractional coverage of mosaic tile on modelled area 
    actlyr    => class_rot%actlyr,                      & !< real, dimension(:,:)  : Active layer depth (m) 
    dlzwrot   => class_rot%dlzwrot,                     & !< real, dimension(:,:,:) : Permeable thickness of soil layer [m] 
    ZPNDROT   => class_rot%ZPNDROT,                     & !< real, dimension(:,:) : Depth of ponded water [m]
    wtableROT => class_rot%wtableROT,                   & !< real, dimension(:,:) : Depth to water table [m]
    PREACC_M  => class_rot%PREACC_M,                    & !< real, dimension(:,:) : Surface precipitation rate \f$[kg m^{-2} ]\f$ (accumulated)
    GTACC_M   => class_rot%GTACC_M,                     & !< real, dimension(:,:) : Diagnosed effective surface black-body temperature [K] (accumulated)
    QEVPACC_M => class_rot%QEVPACC_M,                   & !< real, dimension(:,:) : Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$ (accumulated)
    HFSACC_M  => class_rot%HFSACC_M,                    & !< real, dimension(:,:) : Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$ (accumulated)
    HMFNACC_M => class_rot%HMFNACC_M,                   & !< real, dimension(:,:) : Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$ (accumulated)
    ROFACC_M  => class_rot%ROFACC_M,                    & !< real, dimension(:,:) : Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ROFSACC_M => class_rot%ROFSACC_M,                   & !< real, dimension(:,:) : Interflow from sides of the soil column \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ROFBACC_M => class_rot%ROFBACC_M,                   & !< real, dimension(:,:) : Base flow from bottom of the soil column \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ROFCACC_M => class_rot%ROFCACC_M,                   & !< real, dimension(:,:) : Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ROFNACC_M => class_rot%ROFNACC_M,                   & !< real, dimension(:,:) : Liquid water runoff from snowpack \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ZSNACC_M  => class_rot%ZSNACC_M,                    & !< real, dimension(:,:) : Depth of snow pack \f$[ m ]\f$
    SNOACC_M  => class_rot%SNOACC_M,                    & !< real, dimension(:,:) : Mass of snow pack \f$[kg m^{-2}]\f$ (accumulated)
    BCSNACC_M  => class_rot%BCSNACC_M,                  & !< Black carbon mixing ratio \f$[kg m^{-3}]\f$ (accumulated)
    FSNOACC_M => class_rot%FSNOACC_M,                   & !< real, dimension(:,:) : Fractional snow coverage [ ] (accumulated)
    OVRACC_M  => class_rot%OVRACC_M,                    & !< real, dimension(:,:) : Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    wtableACC_M => class_rot%wtableACC_M,               & !< real, dimension(:,:) : Depth of water table in soil [m]
    TBARACC_M => class_rot%TBARACC_M,                   & !< real, dimension(:,:,:) : Temperature of soil layers [K] (accumulated)
    TCTOACC_M => class_rot%TCTOACC_M,                   & !< real, dimension(:,:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$ (accumulated)
    TCBOACC_M => class_rot%TCBOACC_M,                   & !< real, dimension(:,:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$ (accumulated)
    THLQACC_M => class_rot%THLQACC_M,                   & !< real, dimension(:,:,:) : Volumetric frozen water content of soil layers \f$[kg m^{-2}]\f$ (accumulated)
    THICACC_M => class_rot%THICACC_M,                   & !< real, dimension(:,:,:) : Volumetric liquid water content of soil layers \f$[kg m^{-2}]\f$ (accumulated)
    ALVSACC_M => class_rot%ALVSACC_M,                   & !< real, dimension(:,:) : Diagnosed total visible albedo of land surface [ ] (accumulated)
    ALIRACC_M => class_rot%ALIRACC_M,                   & !< real, dimension(:,:) : Diagnosed total near-infrared albedo of land surface [ ] (accumulated)
    RHOSACC_M => class_rot%RHOSACC_M,                   & !< real, dimension(:,:) : Density of snow \f$[kg m^{-3}]\f$ (accumulated)
    TCSNACC_M => class_rot%TCSNACC_M,                   & !< real, dimension(:,:) : Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
    TSNOACC_M => class_rot%TSNOACC_M,                   & !< real, dimension(:,:) : Snowpack temperature [K] (accumulated)
    TSNBACC_M => class_rot%TSNBACC_M,                   & !< real, dimension(:,:) : Bottom snowpack temperature [K] (accumulated)
    WSNOACC_M => class_rot%WSNOACC_M,                   & !< real, dimension(:,:) : Liquid water content of snow pack \f$[kg m^{-2} ]\f$ (accumulated)
    TCANACC_M => class_rot%TCANACC_M,                   & !< real, dimension(:,:) : Vegetation canopy temperature [K] (accumulated)
    RCANACC_M => class_rot%RCANACC_M,                   & !< real, dimension(:,:) : Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$ (accumulated)
    SCANACC_M => class_rot%SCANACC_M,                   & !< real, dimension(:,:) : Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$ (accumulated)
    GROACC_M  => class_rot%GROACC_M,                    & !< real, dimension(:,:) : Vegetation growth index [ ] (accumulated)
    FSINACC_M => class_rot%FSINACC_M,                   & !< real, dimension(:,:) : Shortwave radiation \f$[W m^{-2} ]\f$ (accumulated)
    FLINACC_M => class_rot%FLINACC_M,                   & !< real, dimension(:,:) : Downwelling longwave sky radiation \f$[W m^{-2} ]\f$ (accumulated)
    TAACC_M   => class_rot%TAACC_M,                     & !< real, dimension(:,:) : Air temperature at reference height [K] (accumulated)
    TSURFACC_M => class_rot%TSURFACC_M,                 & !< real, dimension(:,:) : Ground surface temperature [K] (accumulated)
    ZPNDACC_M => class_rot%ZPNDACC_M,                   & !< real, dimension(:,:) : Depth of ponded water [m] (accumulated)
    !         UVACC_M  => class_rot%UVACC_M,            & !< real, dimension(:,:) : 
    !         PRESACC_M => class_rot%PRESACC_M,         & !< real, dimension(:,:) : 
    !         QAACC_M  => class_rot%QAACC_M,            & !< real, dimension(:,:) : 
    ALTOTACC_M => class_rot%ALTOTACC_M,                 & !< real, dimension(:,:) : Broadband albedo [-] (accumulated) 
    ALSNOACC_M => class_rot%ALSNOACC_M,                 & !< real, dimension(:,:) : Snow albedo [-] (accumulated) 
    EVAPACC_M   => class_rot%EVAPACC_M,                 & !< real, dimension(:,:) : Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$ (accumulated) 
    TRANSPACC_M => class_rot%TRANSPACC_M,               & !< real, dimension(:,:) : Transpiration over modelled area \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    FLUTACC_M   => class_rot%FLUTACC_M,                 & !< real, dimension(:,:) : Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$ (accumulated) 
    GROUNDEVAP_M => class_rot%GROUNDEVAP_M,             & !< real, dimension(:,:) : Surface water vapour flux over ground area \f$[kg m^{-2} ]\f$ (accumulated)
    CANOPYEVAP_M => class_rot%CANOPYEVAP_M,             & !< real, dimension(:,:) : Water vapour flux from the canopy \f$[kg m^{-2} ]\f$ (accumulated)
    ACTLYR_M    => class_rot%ACTLYR_M,                  & !< real, dimension(:,:)  : Active layer depth (m) (accumulated)
    altotcntr_d => class_rot%altotcntr_d                & !< integer, dimension(:) : Used to count the number of time steps with the sun above the horizon 
    )

    ! Accumulate output data for diurnally averaged fields. Both grid mean and mosaic mean
    do I = 1,NLTEST
      do M = 1,NMTEST
        if (FSSROW(I) > 0.) then
          ALTOTACC_M(I,M) = ALTOTACC_M(I,M) + (FSSROW(I) - (FSGVROT(I,M) &
                            + FSGSROT(I,M) + FSGGROT(I,M)))/FSSROW(I)
          ALSNOACC_M(I,M) = ALSNOACC_M(I,M) + ALBSROT(I,M)
          ALVSACC_M(I,M) = ALVSACC_M(I,M) + ALVSROT(I,M)
          ALIRACC_M(I,M) = ALIRACC_M(I,M) + ALIRROT(I,M)
          if (m == 1) altotcntr_d(i) = altotcntr_d(i) + 1 ! s.r.c. altered if statement to prevent per-tile indexing
        end if

        PREACC_M(I,M) = PREACC_M(I,M) + PREROW(I)
        BCSNACC_M(I,M) = BCSNACC_M(I,M) + BCSNROW(I)
        GTACC_M(I,M) = GTACC_M(I,M) + GTROT(I,M)
        QEVPACC_M(I,M) = QEVPACC_M(I,M) + QEVPROT(I,M)
        EVAPACC_M(I,M) = EVAPACC_M(I,M) + QFSROT(I,M)
        HFSACC_M(I,M) = HFSACC_M(I,M) + HFSROT(I,M)
        HMFNACC_M(I,M) = HMFNACC_M(I,M) + HMFNROT(I,M)
        ROFACC_M(I,M) = ROFACC_M(I,M) + ROFROT(I,M)
        ROFSACC_M(I,M) = ROFSACC_M(I,M) + ROFSROT(I,M)
        ROFBACC_M(I,M) = ROFBACC_M(I,M) + ROFBROT(I,M)
        ROFCACC_M(I,M) = ROFCACC_M(I,M) + ROFCROT(I,M)
        ROFNACC_M(I,M) = ROFNACC_M(I,M) + ROFNROT(I,M)
        OVRACC_M(I,M) = OVRACC_M(I,M) + ROFOROT(I,M)
        FSNOACC_M(I,M) = FSNOACC_M(I,M) + FSNOROT(I,M)
        wtableACC_M(I,M)=wtableACC_M(I,M)+wtableROT(I,M)  
        do J = 1,IGND
          TRANSPACC_M(I,M) = TRANSPACC_M(I,M) + QFCROT(I,M,J)
          TBARACC_M(I,M,J) = TBARACC_M(I,M,J) + TBARROT(I,M,J)
          TCTOACC_M(I,M,J) = TCTOACC_M(I,M,J) + TCTOROT(I,M,J)
          TCBOACC_M(I,M,J) = TCBOACC_M(I,M,J) + TCBOROT(I,M,J)
          THLQACC_M(I,M,J) = THLQACC_M(I,M,J) + THLQROT(I,M,J) * 1000. * DLZWROT(I,M,J) ! converted to kg/m2
          THICACC_M(I,M,J) = THICACC_M(I,M,J) + THICROT(I,M,J) * 1000. * DLZWROT(I,M,J) ! converted to kg/m2
        end do
        !ALVSACC_M(I,M) = ALVSACC_M(I,M) + ALVSROT(I,M) * FSVHROW(I)
        !ALIRACC_M(I,M) = ALIRACC_M(I,M) + ALIRROT(I,M) * FSIHROW(I)
        if (SNOROT(I,M) > 0.0) then
          RHOSACC_M(I,M) = RHOSACC_M(I,M) + RHOSROT(I,M)
          TSNOACC_M(I,M) = TSNOACC_M(I,M) + TSNOROT(I,M) - TFREZ ! converted to °C
          TSNBACC_M(I,M) = TSNBACC_M(I,M) + TSNBROT(I,M) - TFREZ ! converted to °C
          WSNOACC_M(I,M) = WSNOACC_M(I,M) + WSNOROT(I,M)
          TCSNACC_M(I,M) = TCSNACC_M(I,M) + TCSNROT(I,M)
        end if
        if (TCANROT(I,M) > 0.5) then
          TCANACC_M(I,M) = TCANACC_M(I,M) + TCANROT(I,M)
        end if
        if (RHOSROT(I,M) > 0.0) ZSNACC_M(I,M) = ZSNACC_M(I,M) + SNOROT(I,M)/RHOSROT(I,M)
        SNOACC_M(I,M) = SNOACC_M(I,M) + SNOROT(I,M)
        GROUNDEVAP_M(I,M) = GROUNDEVAP_M(I,M) + QFGROT(I,M) + QFNROT(I,M)
        CANOPYEVAP_M(I,M) = CANOPYEVAP_M(I,M) + QFCLROT(I,M) + QFCFROT(I,M)
        RCANACC_M(I,M) = RCANACC_M(I,M) + RCANROT(I,M)
        SCANACC_M(I,M) = SCANACC_M(I,M) + SCANROT(I,M)
        GROACC_M(I,M) = GROACC_M(I,M) + GROROT(I,M)
        FSINACC_M(I,M) = FSINACC_M(I,M) + FSSROW(I)  ! not per tile
        FLINACC_M(I,M) = FLINACC_M(I,M) + FDLROW(I)  ! not per tile
        FLUTACC_M(I,M) = FLUTACC_M(I,M) + SBC * GTROT(I,M) ** 4
        TAACC_M(I,M) = TAACC_M(I,M) + TAROW(I)  ! not per tile
        ZPNDACC_M(I,M) = ZPNDACC_M(I,M) + ZPNDROT(I,M)
        ACTLYR_M(I,M) = ACTLYR_M(I,M) + ACTLYR(I,M)  
        ! Find the surface temperature for the grid average over all subareas (1-4)
        TSURFACC_M(I,M) = TSURFACC_M(I,M) + FCS(I) * TSFSGAT(I,1) + FGS(I) * TSFSGAT(I,2) + FC(I) * TSFSGAT(I,3) + FG(I) * TSFSGAT(I,4) ! not per tile
        !                UVACC_M(I,M)=UVACC_M(I,M)+UVROW(I)  ! not per tile
        !                 PRESACC_M(I,M)=PRESACC_M(I,M)+PRESROW(I)
        !                 QAACC_M(I,M)=QAACC_M(I,M)+QAROW(I)
      end do ! loop 50
    end do ! loop 75

    if (NCOUNT == NDAY) then

      ALIRACC(:) = 0.0 ; ALVSACC(:) = 0.0 ; CANOPYEVAP(:) = 0.0 ; EVAPACC(:) = 0.0
      FLINACC(:) = 0.0 ; FLUTACC(:) = 0.0 ; FSINACC(:) = 0.0 ; GROACC(:) = 0.0
      GROUNDEVAP(:) = 0.0 ; GTACC(:) = 0.0 ; HFSACC(:) = 0.0 ;  HMFNACC(:) = 0.0
      OVRACC(:) = 0.0 ;    PREACC(:) = 0.0 ; PRESACC(:) = 0.0 ;   QAACC(:) = 0.0
      QEVPACC(:) = 0.0 ;   RCANACC(:) = 0.0 ; RHOSACC(:) = 0.0 ;  ROFACC(:) = 0.0
      ROFSACC(:) = 0.0 ;  ROFBACC(:) = 0.0 ; ROFCACC(:) = 0.0 ;  ROFNACC(:) = 0.0
      SCANACC(:) = 0.0 ;  SNOACC(:) = 0.0 ; TCSNACC(:) = 0.0 ; FSNOACC(:) = 0.0 ;  TAACC(:) = 0.0   ; TCANACC(:) = 0.0 
      TRANSPACC(:) = 0.0 ; TSNOACC(:) = 0.0 ; TSNBACC(:) = 0.0 ; UVACC(:) = 0.0 ;   WSNOACC(:) = 0.0 ; BCSNACC(:) = 0.0
      wtableACC(:) = 0.0 ; ALTOTACC(:) = 0.0 ; THLQACC(:,:) = 0.0 ; THICACC(:,:) = 0.0
      TBARACC(:,:) = 0.0 ; ZPNDACC(:) = 0.0 ;   ALSNOACC(:) = 0.0 ; ZSNACC(:) = 0.0
      ACTLYRACC(:) = 0.0 ; TSURFACC(:) = 0.0 ; TCTOACC(:,:) = 0.0 ; TCBOACC(:,:) = 0.0

      do I = 1,NLTEST
        do M = 1,NMTEST
          PREACC(I) = PREACC(I) + PREACC_M(I,M) * FAREROT(I,M)
          GTACC(I) = GTACC(I) + GTACC_M(I,M) * FAREROT(I,M)
          QEVPACC(I) = QEVPACC(I) + QEVPACC_M(I,M) * FAREROT(I,M)
          EVAPACC(I) = EVAPACC(I) + EVAPACC_M(I,M) * FAREROT(I,M)
          GROUNDEVAP(I) = GROUNDEVAP(I) + GROUNDEVAP_M(I,M) * FAREROT(I,M)
          CANOPYEVAP(I) = CANOPYEVAP(I) + CANOPYEVAP_M(I,M) * FAREROT(I,M)
          TRANSPACC(I) = TRANSPACC(I) + TRANSPACC_M(I,M) * FAREROT(I,M)
          HFSACC(I) = HFSACC(I) + HFSACC_M(I,M) * FAREROT(I,M)
          HMFNACC(I) = HMFNACC(I) + HMFNACC_M(I,M) * FAREROT(I,M)
          ROFACC(I) = ROFACC(I) + ROFACC_M(I,M) * FAREROT(I,M)
          OVRACC(I) = OVRACC(I) + OVRACC_M(I,M) * FAREROT(I,M)
          ROFSACC(I) = ROFSACC(I) + ROFSACC_M(I,M) * FAREROT(I,M)
          ROFBACC(I) = ROFBACC(I) + ROFBACC_M(I,M) * FAREROT(I,M)
          ROFCACC(I) = ROFCACC(I) + ROFCACC_M(I,M) * FAREROT(I,M)
          ROFNACC(I) = ROFNACC(I) + ROFNACC_M(I,M) * FAREROT(I,M)
          wtableACC(I) = wtableACC(I) + wtableACC_M(I,M) * FAREROT(I,M)
          ALTOTACC(I) = ALTOTACC(I) + ALTOTACC_M(I,M) * FAREROT(I,M)
          ALSNOACC(I) = ALSNOACC(I) + ALSNOACC_M(I,M) * FAREROT(I,M)
          do J = 1,IGND
            TBARACC(I,J) = TBARACC(I,J) + TBARACC_M(I,M,J) * FAREROT(I,M)
            TCTOACC(I,J) = TCTOACC(I,J) + TCTOACC_M(I,M,J) * FAREROT(I,M)
            TCBOACC(I,J) = TCBOACC(I,J) + TCBOACC_M(I,M,J) * FAREROT(I,M)
            THLQACC(I,J) = THLQACC(I,J) + THLQACC_M(I,M,J) * FAREROT(I,M)
            THICACC(I,J) = THICACC(I,J) + THICACC_M(I,M,J) * FAREROT(I,M)
          end do
          ALVSACC(I) = ALVSACC(I) + ALVSACC_M(I,M) * FAREROT(I,M)
          ALIRACC(I) = ALIRACC(I) + ALIRACC_M(I,M) * FAREROT(I,M)
          RHOSACC(I) = RHOSACC(I) + RHOSACC_M(I,M) * FAREROT(I,M)
          TCSNACC(I) = TCSNACC(I) + TCSNACC_M(I,M) * FAREROT(I,M)
          TSNOACC(I) = TSNOACC(I) + TSNOACC_M(I,M) * FAREROT(I,M)
          TSNBACC(I) = TSNBACC(I) + TSNBACC_M(I,M) * FAREROT(I,M)
          WSNOACC(I) = WSNOACC(I) + WSNOACC_M(I,M) * FAREROT(I,M)
          BCSNACC(I) = BCSNACC(I) + BCSNACC_M(I,M) * FAREROT(I,M)
          TCANACC(I) = TCANACC(I) + TCANACC_M(I,M) * FAREROT(I,M)
          ZSNACC(I) = ZSNACC(I) + ZSNACC_M(I,M) * FAREROT(I,M)
          SNOACC(I) = SNOACC(I) + SNOACC_M(I,M) * FAREROT(I,M)
          FSNOACC(I) = FSNOACC(I) + FSNOACC_M(I,M) * FAREROT(I,M)
          RCANACC(I) = RCANACC(I) + RCANACC_M(I,M) * FAREROT(I,M)
          SCANACC(I) = SCANACC(I) + SCANACC_M(I,M) * FAREROT(I,M)
          GROACC(I) = GROACC(I) + GROACC_M(I,M) * FAREROT(I,M)
          FSINACC(I) = FSINACC(I) + FSINACC_M(I,M) * FAREROT(I,M)
          FLINACC(I) = FLINACC(I) + FLINACC_M(I,M) * FAREROT(I,M)
          FLUTACC(I) = FLUTACC(I) + FLUTACC_M(I,M) * FAREROT(I,M)
          TAACC(I) = TAACC(I) + TAACC_M(I,M) * FAREROT(I,M)
          ZPNDACC(I) = ZPNDACC(I) + ZPNDACC_M(I,M) * FAREROT(I,M)
          ! Find the surface temperature for the grid average over all subareas (1-4)
          ! TSURFACC(I) = TSURFACC(I) + FCS(I) * TSFSGAT(I,1) + FGS(I) * TSFSGAT(I,2) + FC(I) * TSFSGAT(I,3) + FG(I) * TSFSGAT(I,4)
          TSURFACC(I) = TSURFACC(I) + TSURFACC_M(I,M) * FAREROT(I,M)
          ACTLYRACC(I) = ACTLYRACC(I) + ACTLYR_M(I,M) * FAREROT(I,M)

          ! UVACC(I)=UVACC(I)+UVACC_M(I)*FAREROT(I,M)
          ! PRESACC(I)=PRESACC(I)+PRESROW(I)*FAREROT(I,M)
          ! QAACC(I)=QAACC(I)+QAROW(I)*FAREROT(I,M)
        end do
      end do

      ! Now write to file the grid average values

      ! Transfer the consecDays to timeStamp (since we need a size 1 array)
      timeStamp = consecDays

      do i = 1,nltest
        if (altotcntr_d(i) > 0) then
          ALTOTACC(I) = ALTOTACC(I)/real(altotcntr_d(i))
          ALSNOACC(I) = ALSNOACC(I)/real(altotcntr_d(i))
          ALVSACC(I) = ALVSACC(I)/real(altotcntr_d(i))
          ALIRACC(I) = ALIRACC(I)/real(altotcntr_d(i))
        else
          ALTOTACC(I) = 0.
          ALSNOACC(I) = 0.
          ALVSACC(I) = 0.
          ALIRACC(I) = 0.
        end if
        FSSTAR = FSINACC(I)/real(nday) * (1. - ALTOTACC(I))
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_d' ,timeStamp,'rss', [FSSTAR])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc',timeStamp,'albs',[ALTOTACC(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'alsnoacc',timeStamp,'albsn',[ALSNOACC(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'alvsacc',timeStamp,'albsvis',[ALVSACC(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'aliracc',timeStamp,'albsir',[ALIRACC(I)])
        FLSTAR = (FLINACC(I) - FLUTACC(I))/real(NDAY)
        call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_d' ,timeStamp,'rls', [FLSTAR])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_d'     ,timeStamp,'hfss', [HFSACC(I)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_d'     ,timeStamp,'hfls', [QEVPACC(I)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snm_d'    ,timeStamp,'snm', [HMFNACC(I)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tbaracc_d',timeStamp,'tsl', [(TBARACC(I,:)/real(NDAY))])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tctop_d',timeStamp,'tctop', [(TCTOACC(I,:)/real(NDAY))])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tcbot_d',timeStamp,'tcbot', [(TCBOACC(I,:)/real(NDAY))])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thlqacc_d',timeStamp,'mrsll', [THLQACC(I,:)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thicacc_d',timeStamp,'mrsfl', [THICACC(I,:)/real(NDAY)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_d',timeStamp,'tas', [(TAACC(I)/REAL(NDAY))-TFREZ])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_d',timeStamp,'tas', [(TAACC(I)/REAL(NDAY))])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snd_d',timeStamp,'snd', [ZSNACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snoacc_d',timeStamp,'snw', [SNOACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsnoacc_d',timeStamp,'snc', [FSNOACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_d',timeStamp,'wsnw', [WSNOACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tsnoacc_d',timeStamp,'tsn', [TSNOACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tsnbacc_d',timeStamp,'tsnbot', [TSNBACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rhosacc_d',timeStamp,'snwdens', [RHOSACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tcsnow_d',timeStamp,'tcsnow', [TCSNACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_d',timeStamp,'pr', [PREACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'scanacc_d',timeStamp,'scanopy', [SCANACC(I)/REAL(NDAY)])
        !call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_d',timeStamp,'evspsbl', [EVAPACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'evspsbl_d',timeStamp,'evspsbl', [(CANOPYEVAP(I) + GROUNDEVAP(I) + TRANSPACC(I))/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap_d',timeStamp,'evspsblsoi', [GROUNDEVAP(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap_d',timeStamp,'evspsblveg', [CANOPYEVAP(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_d',timeStamp,'tran', [TRANSPACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_d',timeStamp,'mrro', [ROFACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofsacc_d',timeStamp,'mrroi', [ROFSACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofbacc_d',timeStamp,'mrrob', [ROFBACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ovracc_d',timeStamp,'mrros', [OVRACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofcacc_d',timeStamp,'mrroc', [ROFCACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofnacc_d',timeStamp,'mrron', [ROFNACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'zpndacc_d',timeStamp,'zpond', [ZPNDACC(I)/REAL(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'wtdacc_d',timeStamp,'wtd', [wtableACC(I)/REAL(NDAY)])
        TCANACC(I) = TCANACC(I)/REAL(NDAY)
        if (TCANACC(I) > 0.01) then
          TCN = TCANACC(I) - TFREZ
        else
          TCN = 0.
        end if
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tcanacc_d',timeStamp,'tcs', [TCN])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'tcanacc_d',timeStamp,'tcs', [TCANACC(I)/REAL(NDAY)-TFREZ])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurfacc_d',timeStamp,'ts', [(TSURFACC(I)/REAL(NDAY))])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_d',timeStamp,'actlyr', [ACTLYRACC(I)/REAL(NDAY)])
        
        !                 !                 BEG=FSSTAR+FLSTAR-QH-QE

        ! Below are obviously not daily ones but I copy in as a reminder of what likely will be put out daily.
        ! after these writeOutput1D statements is the listed daily file variables. It would be good to do those again
        !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'gflx_hh',timeStamp,'gflx', [GFLXROW(I,:)])  ! name
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'snoacc_mo' ,timeStamp,'snw', [HMFNACC(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_mo',timeStamp,'wsnw', [WSNOACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_mo'  ,timeStamp,'tas', [TAACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap',timeStamp,'evspsblsoi', [GROUNDEVAP(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap',timeStamp,'evspsblveg', [CANOPYEVAP(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_mo' ,timeStamp,'mrro', [ROFACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_mo' ,timeStamp,'pr', [PREACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_mo',timeStamp,'evspsbl', [EVAPACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_mo',timeStamp,'tran', [TRANSPACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'alsnoacc',timeStamp,'albsn', [ALSNOACC(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_mo',timeStamp,'albs', [ALTOTACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_mo',timeStamp,'actlyr', [ACTLYR_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_max_mo',timeStamp,'actlyrmax', [ACTLYR_MAX_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_min_mo',timeStamp,'actlyrmin', [ACTLYR_MIN_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_mo',timeStamp,'ftable', [FTABLE_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_max_mo',timeStamp,'ftablemax', [FTABLE_MAX_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_min_mo',timeStamp,'ftablemin', [FTABLE_MIN_MO(I)])
        !
        ! !                             &                       BEG,GTOUT,SNOACC(I),RHOSACC(I), &
        ! !                             &                       WSNOACC(I),ALTOTACC(I),ROFACC(I),CUMSNO
        ! !                             WRITE(62,6200) IDAY,realyr,(TBARACC(I,J)-TFREZ, &
        ! !                                 &                       THLQACC(I,J),THICACC(I,J),J=1,3), &
        ! !                                 &                       TCN,RCANACC(I),SCANACC(I),TSN,ZSN, &
        ! !                                 &                       ACTLYR_G(I),FTABLE_g(I)
        ! !                         WRITE(63,6300) IDAY,realyr,FSINACC(I),FLINACC(I), &
        ! !                             &                       TAACC(I)-TFREZ,UVACC(I),PRESACC(I), &
        ! !                             &                       QAACC(I),PREACC(I),EVAPACC(I)
        ! !                                 &                    BEG,GTOUT,SNOACC_M(I,M),RHOSACC_M(I,M), &
        ! !                                 &                    WSNOACC_M(I,M),ALTOTACC_M(I,M),ROFACC_M(I,M), &
        ! !                                 &                    CUMSNO,' TILE ',M
        ! !                                     &                  TCN,RCANACC_M(I,M),SCANACC_M(I,M),TSN,ZSN, &
        ! !                                     &                  ' TILE ',M
        ! !                             END IF
        ! !                             WRITE(631,6300) IDAY,realyr,FSINACC_M(I,M),FLINACC_M(I,M), &
        ! !                                 &                  TAACC_M(I,M)-TFREZ,UVACC_M(I,M),PRESACC_M(I,M), &
        ! !                                 &                  QAACC_M(I,M),PREACC_M(I,M),EVAPACC_M(I,M), &
        ! !                                 &                  ' TILE ',M
        !
        ! !
      end do
      ! !             do I=1,NLTEST
      ! !                 PREACC(I)=PREACC(I)
      ! !                 GTACC(I)=GTACC(I)/REAL(NDAY)
      ! !                 EVAPACC(I)=EVAPACC(I)
      ! !                 HMFNACC(I)=HMFNACC(I)/REAL(NDAY)
      ! !                 ROFACC(I)=ROFACC(I)
      ! !                 OVRACC(I)=OVRACC(I)
      ! !                 do J=1,IGND
      ! !                     TBARACC(I,J)=TBARACC(I,J)/REAL(NDAY)
      ! !                     THLQACC(I,J)=THLQACC(I,J)/REAL(NDAY)
      ! !                     THICACC(I,J)=THICACC(I,J)/REAL(NDAY)
      ! !                     THALACC(I,J)=THALACC(I,J)/REAL(NDAY)
      ! ! 725                     CONTINUE
      ! !                 IF (FSINACC(I)>0.0) THEN
      ! !                     ALVSACC(I)=ALVSACC(I)/(FSINACC(I)*0.5)
      ! !                     ALIRACC(I)=ALIRACC(I)/(FSINACC(I)*0.5)
      ! !                 ELSE
      ! !                     ALVSACC(I)=0.0
      ! !                     ALIRACC(I)=0.0
      ! !                 END IF
      ! !                 IF (SNOARE(I)>0.0) THEN
      ! !                     RHOSACC(I)=RHOSACC(I)/SNOARE(I)
      ! !                     TSNOACC(I)=TSNOACC(I)/SNOARE(I)
      ! !                     WSNOACC(I)=WSNOACC(I)/SNOARE(I)
      ! !                 END IF
      ! !                 IF (CANARE(I)>0.0) THEN
      ! !                     TCANACC(I)=TCANACC(I)/CANARE(I)
      ! !                 END IF
      ! !                 SNOACC(I)=SNOACC(I)/REAL(NDAY)
      ! !                 RCANACC(I)=RCANACC(I)/REAL(NDAY)
      ! !                 SCANACC(I)=SCANACC(I)/REAL(NDAY)
      ! !                 GROACC(I)=GROACC(I)/REAL(NDAY)
      ! !                 FSINACC(I)=FSINACC(I)/REAL(NDAY)
      ! !                 FLINACC(I)=FLINACC(I)/REAL(NDAY)
      ! !                 FLUTACC(I)=FLUTACC(I)/REAL(NDAY)
      ! !                 TAACC(I)=TAACC(I)/REAL(NDAY)
      ! !                 UVACC(I)=UVACC(I)/REAL(NDAY)
      ! !                 PRESACC(I)=PRESACC(I)/REAL(NDAY)
      ! !                 QAACC(I)=QAACC(I)/REAL(NDAY)
      ! !                 IF (RHOSACC(I)>0.0) THEN
      ! !                     ZSN=SNOACC(I)/RHOSACC(I)
      ! !                 ELSE
      ! !                     ZSN=0.0
      ! !                 END IF
      ! !                 IF (TCANACC(I)>0.01) THEN
      ! !                     TCN=TCANACC(I)-TFREZ
      ! !                 ELSE
      ! !                     TCN=0.0
      ! !                 END IF
      ! !                 IF (TSNOACC(I)>0.01) THEN
      ! !                     TSN=TSNOACC(I)-TFREZ
      ! !                 ELSE
      ! !                     TSN=0.0
      ! !                 END IF
      ! !                 GTOUT=GTACC(I)-TFREZ
      ! !!
      ! !
      ! !                     end if
      ! !                 END IF
      ! !
      ! !             do I=1,NLTEST
      ! !                 do M=1,NMTEST
      ! !                     PREACC_M(I,M)=PREACC_M(I,M)     ! became [kg m-2 day-1] instead of [kg m-2 s-1]
      ! !                     GTACC_M(I,M)=GTACC_M(I,M)/REAL(NDAY)
      ! !                     QEVPACC_M(I,M)=QEVPACC_M(I,M)/REAL(NDAY)
      ! !                     EVAPACC_M(I,M)=EVAPACC_M(I,M)   ! became [kg m-2 day-1] instead of [kg m-2 s-1]
      ! !                     HFSACC_M(I,M)=HFSACC_M(I,M)/REAL(NDAY)
      ! !                     HMFNACC_M(I,M)=HMFNACC_M(I,M)/REAL(NDAY)
      ! !                     ROFACC_M(I,M)=ROFACC_M(I,M)   ! became [kg m-2 day-1] instead of [kg m-2 s-1
      ! !                     OVRACC_M(I,M)=OVRACC_M(I,M)   ! became [kg m-2 day-1] instead of [kg m-2 s-1]
      ! !
      ! !                     IF (FSINACC_M(I,M)>0.0) THEN
      ! !                         ALVSACC_M(I,M)=ALVSACC_M(I,M)/(FSINACC_M(I,M)*0.5)
      ! !                         ALIRACC_M(I,M)=ALIRACC_M(I,M)/(FSINACC_M(I,M)*0.5)
      ! !                     ELSE
      ! !                         ALVSACC_M(I,M)=0.0
      ! !                         ALIRACC_M(I,M)=0.0
      ! !                     END IF
      ! !
      ! !                     SNOACC_M(I,M)=SNOACC_M(I,M)/REAL(NDAY)
      ! !                     if (SNOARE_M(I,M) > 0.) THEN
      ! !                         RHOSACC_M(I,M)=RHOSACC_M(I,M)/SNOARE_M(I,M)
      ! !                         TSNOACC_M(I,M)=TSNOACC_M(I,M)/SNOARE_M(I,M)
      ! !                         WSNOACC_M(I,M)=WSNOACC_M(I,M)/SNOARE_M(I,M)
      ! !                     END IF
      ! !                     TCANACC_M(I,M)=TCANACC_M(I,M)/REAL(NDAY)
      ! !                     RCANACC_M(I,M)=RCANACC_M(I,M)/REAL(NDAY)
      ! !                     SCANACC_M(I,M)=SCANACC_M(I,M)/REAL(NDAY)
      ! !                     GROACC_M(I,M)=GROACC_M(I,M)/REAL(NDAY)
      ! !                     FSINACC_M(I,M)=FSINACC_M(I,M)/REAL(NDAY)
      ! !                     FLINACC_M(I,M)=FLINACC_M(I,M)/REAL(NDAY)
      ! !                     FLUTACC_M(I,M)=FLUTACC_M(I,M)/REAL(NDAY)
      ! !                     TAACC_M(I,M)=TAACC_M(I,M)/REAL(NDAY)
      ! !                     UVACC_M(I,M)=UVACC_M(I,M)/REAL(NDAY)
      ! !                     PRESACC_M(I,M)=PRESACC_M(I,M)/REAL(NDAY)
      ! !                     QAACC_M(I,M)=QAACC_M(I,M)/REAL(NDAY)
      ! !                     if (altotcntr_d(i) > 0) then ! altotcntr_d(i) could be 0
      ! !                         ALTOTACC_M(I,M)=ALTOTACC_M(I,M)/REAL(altotcntr_d(i))
      ! !                     else
      ! !                         ALTOTACC_M(I,M)=0.
      ! !                     end if
      ! !                     FSSTAR=FSINACC_M(I,M)*(1.-ALTOTACC_M(I,M))
      ! !                     FLSTAR=FLINACC_M(I,M)-FLUTACC_M(I,M)
      ! !                     QH=HFSACC_M(I,M)
      ! !                     QE=QEVPACC_M(I,M)
      ! !                     QEVPACC_M_SAVE(I,M)=QEVPACC_M(I,M)   ! FLAG ! What is the point of this? JM Apr 12015
      ! !                     BEG=FSSTAR+FLSTAR-QH-QE
      ! !                     SNOMLT=HMFNACC_M(I,M)
      ! !
      ! !                     IF (RHOSACC_M(I,M)>0.0) THEN
      ! !                         ZSN=SNOACC_M(I,M)/RHOSACC_M(I,M)
      ! !                     ELSE
      ! !                         ZSN=0.0
      ! !                     END IF
      ! !
      ! !                     IF (TCANACC_M(I,M)>0.01) THEN
      ! !                         TCN=TCANACC_M(I,M)-TFREZ
      ! !                     ELSE
      ! !                         TCN=0.0
      ! !                     END IF
      ! !
      ! !                     IF (TSNOACC_M(I,M)>0.01) THEN
      ! !                         TSN=TSNOACC_M(I,M)-TFREZ
      ! !                     ELSE
      ! !                         TSN=0.0
      ! !                     END IF
      ! !
      ! !                     GTOUT=GTACC_M(I,M)-TFREZ
      ! !
      ! !                             !         WRITE TO OUTPUT FILES
      ! !                             !
      ! !                         end if
      ! !
      ! ! 809                     CONTINUE
      ! ! 808                 CONTINUE

      ! RESET ACCUMULATOR ARRAYS (*_M).
      call resetAccVars(nltest, nmtest)

    end if ! IF (NCOUNT==NDAY)

    end associate
  end subroutine class_daily_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_class_monthly_aw
  !> @{
  !> Accumulate and write out the monthly physics outputs. These are kept in pointer structures as
  !! this subroutine is called each physics timestep and we increment the timestep values to produce a monthly value.
  !! The pointer to the monthly data structures (in classStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine class_monthly_aw (lonLocalIndex, latLocalIndex, IDAY, realyr, NCOUNT, NDAY, nltest, nmtest, lastDOY)

    use classStateVars, only : class_out, resetClassMon, class_rot, class_gat
    use classicParams,  only : nmon, monthend, nmos, ignd, SBC, DELT, TFREZ
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: IDAY
    integer, intent(in) :: realyr
    integer, intent(in) :: NCOUNT
    integer, intent(in) :: NDAY
    integer, intent(in) :: lastDOY
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    ! local
    real :: ALTOT_MO
    integer :: NT
    integer :: NDMONTH
    integer :: i, m, j
    integer :: IMONTH
    real :: tovere
    real :: ACTLYR_tmp
    real :: FTABLE_tmp
    real :: FSSTAR_MO
    real :: FLSTAR_MO
    real :: QH_MO
    real :: QE_MO
    real :: TCN
    real, dimension(1) :: timeStamp

    ! Associate names with variables defined in derived types.
    associate( &
    TSFSGAT => class_gat%TSFSGAT,                       & !< real, dimension(:,:) : Ground surface temperature over subarea [K]
    FG => class_gat%FG,                                 & !< real, dimension(:) : Subarea fractional coverage of modelled area - bare ground [ ]
    FC => class_gat%FC,                                 & !< real, dimension(:) : Subarea fractional coverage of modelled area - ground under canopy [ ]
    FGS => class_gat%FGS,                               & !< real, dimension(:) : Subarea fractional coverage of modelled area - snow-covered bare ground [ ]
    FCS => class_gat%FCS,                               & !< real, dimension(:) : Subarea fractional coverage of modelled area - snow-covered ground under canopy  [ ]
    TBARROT         => class_rot%TBARROT,               & !< real(r8), dimension(:,:,:): 
    TCTOROT         => class_rot%TCTOROT,               & !< real, dimension(:,:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$
    TCBOROT         => class_rot%TCBOROT,               & !< real, dimension(:,:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$
    THLQROT         => class_rot%THLQROT,               & !< real, dimension(:,:,:): 
    THICROT         => class_rot%THICROT,               & !< real, dimension(:,:,:): 
    QFCROT          => class_rot%QFCROT,                & !< real, dimension(:,:,:): 
    ALVSROT         => class_rot%ALVSROT,               & !< real, dimension(:,:): 
    ALBSROT         => class_rot%ALBSROT,               & !< real, dimension(:,:): 
    FAREROT         => class_rot%FAREROT,               & !< real, dimension(:,:): 
    ALIRROT         => class_rot%ALIRROT,               & !< real, dimension(:,:): 
    GTROT           => class_rot%GTROT,                 & !< real, dimension(:,:): 
    HFSROT          => class_rot%HFSROT,                & !< real, dimension(:,:): 
    QEVPROT         => class_rot%QEVPROT,               & !< real, dimension(:,:): 
    groundHeatFluxROT => class_rot%groundHeatFluxROT,   & !< real, dimension(:,:): Heat flux at soil surface \f$[W m^{-2} ]\f$ 
    SCANROT         => class_rot%SCANROT,               & !< real, dimension(:,:): Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$ 
    RCANROT         => class_rot%RCANROT,               & !< real, dimension(:,:): Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$ 
    BCSNROW         => class_rot%BCSNROW,               & !< real, dimension(:,:): 
    SNOROT          => class_rot%SNOROT,                & !< real, dimension(:,:): 
    FSNOROT          => class_rot%FSNOROT,              & !< real, dimension(:,:): 
    WSNOROT         => class_rot%WSNOROT,               & !< real, dimension(:,:): 
    RHOSROT         => class_rot%RHOSROT,               & !< real, dimension(:,:): Density of snow \f$[kg m^{-3}]\f$ 
    ROFOROT         => class_rot%ROFOROT,               & !< real, dimension(:,:): Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    ROFBROT         => class_rot%ROFBROT,               & !< real, dimension(:,:): Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    ROFROT          => class_rot%ROFROT,                & !< real, dimension(:,:): 
    ROFNROT         => class_rot%ROFNROT,               & !< real, dimension(:,:): 
    QFSROT          => class_rot%QFSROT,                & !< real, dimension(:,:): 
    QFGROT          => class_rot%QFGROT,                & !< real, dimension(:,:): 
    QFNROT          => class_rot%QFNROT,                & !< real, dimension(:,:): 
    QFCLROT         => class_rot%QFCLROT,               & !< real, dimension(:,:): 
    QFCFROT         => class_rot%QFCFROT,               & !< real, dimension(:,:): 
    FSGVROT         => class_rot%FSGVROT,               & !< real, dimension(:,:): Diagnosed net shortwave radiation on vegetation canopy 
    FSGSROT         => class_rot%FSGSROT,               & !< real, dimension(:,:): Diagnosed net shortwave radiation on ground snow surface 
    FSGGROT         => class_rot%FSGGROT,               & !< real, dimension(:,:): Diagnosed net shortwave radiation on ground surface 
    FSSROW          => class_rot%FSSROW,                & !< real, dimension(:): 
    FDLROW          => class_rot%FDLROW,                & !< real, dimension(:): 
    FSVHROW         => class_rot%FSVHROW,               & !< real, dimension(:): 
    FSIHROW         => class_rot%FSIHROW,               & !< real, dimension(:): 
    TAROW           => class_rot%TAROW,                 & !< real, dimension(:): 
    PREROW          => class_rot%PREROW,                & !< real, dimension(:): 
    ftable          => class_rot%ftable,                & !< real, dimension(:,:)  : Depth to frozen water table (m) 
    actlyr          => class_rot%actlyr,                & !< real, dimension(:,:)  : Active layer depth (m) 
    dlzwrot         => class_rot%dlzwrot,               & !< real, dimension(:,:,:) : Permeable thickness of soil layer [m] 
    TCANROT         => class_rot%TCANROT,               & !< real, dimension(:,:) : Vegetation canopy temperature [K] 
    ALVSACC_MO        => class_out%ALVSACC_MO,          & !< real, dimension(:) : 
    ALIRACC_MO        => class_out%ALIRACC_MO,          & !< real, dimension(:) : 
    FLUTACC_MO        => class_out%FLUTACC_MO,          & !< real, dimension(:) : 
    FSINACC_MO        => class_out%FSINACC_MO,          & !< real, dimension(:) : Surface Downwelling Shortwave Radiative flux in air [$W m^{-2}$] 
    FLINACC_MO        => class_out%FLINACC_MO,          & !< real, dimension(:) : 
    HFSACC_MO         => class_out%HFSACC_MO,           & !< real, dimension(:) : 
    QEVPACC_MO        => class_out%QEVPACC_MO,          & !< real, dimension(:) : 
    groundHeatFlux_MO => class_out%groundHeatFlux_MO,   & !< real, dimension(:) : Heat flux at soil surface \f$[W m^{-2} ]\f$ 
    SCAN_MO           => class_out%SCAN_MO,             & !< real, dimension(:) : Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$ 
    RCAN_MO           => class_out%RCAN_MO,             & !< real, dimension(:) : Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$ 
    BCSNACC_MO        => class_out%BCSNACC_MO,          & !< real, dimension(:) : Black carbon mixing ratio \f$[kg m^{-3} ]\f$ 
    SNOACC_MO         => class_out%SNOACC_MO,           & !< real, dimension(:) : Mass of snow pack \f$[kg m^{-2} ]\f$ 
    RHOSACC_MO         => class_out%RHOSACC_MO,           & !< real, dimension(:) : Density of snow pack \f$[kg m^{-3} ]\f$ 
    FSNOACC_MO         => class_out%FSNOACC_MO,           & !< real, dimension(:) : Fractional snow cover [ ] 
    ZSNACC_MO         => class_out%ZSNACC_MO,           & !< real, dimension(:) : Depth of snow pack \f$[ m ]\f$ 
    WSNOACC_MO        => class_out%WSNOACC_MO,          & !< real, dimension(:) : Liquid water content of snow pack \f$[kg m^{-2} ]\f$ 
    OVRACC_MO         => class_out%OVRACC_MO,           & !< real, dimension(:) : Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ROFBACC_MO        => class_out%ROFBACC_MO,          & !< real, dimension(:) : Base flow from bottom of the soil column \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ROFNACC_MO        => class_out%ROFNACC_MO,          & !< real, dimension(:) : Liquid water rubnoff from snowpack \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    ROFACC_MO         => class_out%ROFACC_MO,           & !< real, dimension(:) : Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$ 
    PREACC_MO         => class_out%PREACC_MO,           & !< real, dimension(:) : Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$ 
    EVAPACC_MO        => class_out%EVAPACC_MO,          & !< real, dimension(:) : Diagnosed total surface evaporation water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$ 
    TRANSPACC_MO      => class_out%TRANSPACC_MO,        & !< real, dimension(:) : 
    TAACC_MO          => class_out%TAACC_MO,            & !< real, dimension(:) : Air temperature at reference height [K] 
    TBARACC_MO        => class_out%TBARACC_MO,          & !< real, dimension(:,:) : 
    TCTOACC_MO        => class_out%TCTOACC_MO,          & !< real, dimension(:,:) : Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$ (accumulated)
    TCBOACC_MO        => class_out%TCBOACC_MO,          & !< real, dimension(:,:) : Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$ (accumulated)
    THLQACC_MO        => class_out%THLQACC_MO,          & !< real, dimension(:,:) : Volumetric liquid water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means) 
    THICACC_MO        => class_out%THICACC_MO,          & !< real, dimension(:,:) : Volumetric frozen water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means) 
    ACTLYR_MO         => class_out%ACTLYR_MO,           & !< real, dimension(:) : 
    FTABLE_MO         => class_out%FTABLE_MO,           & !< real, dimension(:) : 
    ACTLYR_MIN_MO     => class_out%ACTLYR_MIN_MO,       & !< real, dimension(:) : 
    FTABLE_MIN_MO     => class_out%FTABLE_MIN_MO,       & !< real, dimension(:) : 
    ACTLYR_MAX_MO     => class_out%ACTLYR_MAX_MO,       & !< real, dimension(:) : 
    FTABLE_MAX_MO     => class_out%FTABLE_MAX_MO,       & !< real, dimension(:) : 
    GROUNDEVAP        => class_out%GROUNDEVAP,          & !< real, dimension(:) : 
    CANOPYEVAP        => class_out%CANOPYEVAP,          & !< real, dimension(:) : 
    EVSPSBL_MO        => class_out%EVSPSBL_MO,          & !< real, dimension(:) : 
    ALTOTACC_MO       => class_out%ALTOTACC_MO,         & !< real, dimension(:) : 
    ALSNOACC_MO       => class_out%ALSNOACC_MO,         & !< real, dimension(:) : 
    altotcntr_m       => class_out%altotcntr_m,         & !< integer, dimension(:) : 
    MRSO_MO           => class_out%MRSO_MO,             & !< real, dimension(:) : 
    MRSOL_MO          => class_out%MRSOL_MO,            & !< real, dimension(:,:) : 
    TSURFACC_MO       => class_out%TSURFACC_MO,         & !< real, dimension(:) : Ground surface temperature[K]
    TCANACC_MO        => class_out%TCANACC_MO           & !< real, dimension(:) : Vegetation canopy temperature [K] (accumulated)
    )

    ! ------------

    !> Accumulate output data for monthly averaged fields for class grid-mean.
    !> for both parallel mode and stand alone mode

    FSSTAR_MO   = 0.0
    FLSTAR_MO   = 0.0
    QH_MO       = 0.0
    QE_MO       = 0.0
    ACTLYR_tmp  = 0.0
    FTABLE_tmp  = 0.0

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do M = 1,NMTEST

      ! These are presently not being outputted but the code is kept in place if the need arises.
      !     ALVSACC_MO(I)=ALVSACC_MO(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
      !     ALIRACC_MO(I)=ALIRACC_MO(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)
      FLUTACC_MO(I) = FLUTACC_MO(I) + SBC * GTROT(I,M) ** 4 * FAREROT(I,M)
      FSINACC_MO(I) = FSINACC_MO(I) + FSSROW(I) * FAREROT(I,M)
      FLINACC_MO(I) = FLINACC_MO(I) + FDLROW(I) * FAREROT(I,M)
      HFSACC_MO(I) = HFSACC_MO(I) + HFSROT(I,M) * FAREROT(I,M)
      QEVPACC_MO(I) = QEVPACC_MO(I) + QEVPROT(I,M) * FAREROT(I,M)
      groundHeatFlux_MO(I) = groundHeatFlux_MO(I) + groundHeatFluxROT(I,M) * FAREROT(I,M) !*()*()*()*()*()*()
      SNOACC_MO(I) = SNOACC_MO(I) + SNOROT(I,M) * FAREROT(I,M)
      RHOSACC_MO(I) = RHOSACC_MO(I) + RHOSROT(I,M) * FAREROT(I,M)
      SCAN_MO(I) = SCAN_MO(I) + SCANROT(I,M) * FAREROT(I,M)
      RCAN_MO(I) = RCAN_MO(I) + RCANROT(I,M) * FAREROT(I,M)
      FSNOACC_MO(I) = FSNOACC_MO(I) + FSNOROT(I,M) * FAREROT(I,M)
      if (RHOSROT(I,M) > 0.0) ZSNACC_MO(I) = &
         ZSNACC_MO(I) + SNOROT(I,M)/RHOSROT(I,M) * FAREROT(I,M)

      TAACC_MO(I) = TAACC_MO(I) + TAROW(I) * FAREROT(I,M)
      ACTLYR_MO(I) = ACTLYR_MO(I) + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_MO(I) = FTABLE_MO(I) + FTABLE(I,M) * FAREROT(I,M)
      ACTLYR_tmp = ACTLYR_tmp + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_tmp = FTABLE_tmp + FTABLE(I,M) * FAREROT(I,M)
      GROUNDEVAP(I) = GROUNDEVAP(I) + (QFGROT(I,M) + QFNROT(I,M)) * FAREROT(I,M) ! ground evap includes both evap and sublimation from snow
      CANOPYEVAP(I) = CANOPYEVAP(I) + (QFCLROT(I,M) + QFCFROT(I,M)) * FAREROT(I,M) ! canopy evap includes both evap and sublimation
      ! Find the surface temperature for the grid average over all subareas (1-4)
      TSURFACC_MO(I) = TSURFACC_MO(I) + FCS(I) * TSFSGAT(I,1) + FGS(I) * TSFSGAT(I,2) + FC(I) * TSFSGAT(I,3) + FG(I) * TSFSGAT(I,4)

      if (SNOROT(I,M) > 0.0) then
        WSNOACC_MO(I) = WSNOACC_MO(I) + WSNOROT(I,M) * FAREROT(I,M)
        BCSNACC_MO(I) = BCSNACC_MO(I) + BCSNROW(I) * FAREROT(I,M)
      end if

      if (TCANROT(I,M) > 0.5) then
        TCANACC_MO(I) = TCANACC_MO(I) + TCANROT(I,M) * FAREROT(I,M)
      end if

      OVRACC_MO(I) = OVRACC_MO(I) + ROFOROT(I,M) * FAREROT(I,M)
      ROFBACC_MO(I) = ROFBACC_MO(I) + ROFBROT(I,M) * FAREROT(I,M)
      ROFNACC_MO(I) = ROFNACC_MO(I) + ROFNROT(I,M) * FAREROT(I,M)
      ROFACC_MO(I) = ROFACC_MO(I) + ROFROT(I,M) * FAREROT(I,M)
      PREACC_MO(I) = PREACC_MO(I) + PREROW(I) * FAREROT(I,M)
      EVAPACC_MO(I) = EVAPACC_MO(I) + QFSROT(I,M) * FAREROT(I,M) ! Only evaporation

      if (FSSROW(I) > 0.0) then
        ALTOTACC_MO(I) = ALTOTACC_MO(I) + ( (FSSROW(I) - (FSGVROT(I,M) + FSGSROT(I,M) + FSGGROT(I,M))) &
                         /FSSROW(I) ) * FAREROT(I,M)
        ALSNOACC_MO(I) = ALSNOACC_MO(I) + ALBSROT(I,M)
        if (M == 1) altotcntr_m(i) = altotcntr_m(i) + 1 ! s.r.c. added if statement to prevent per-tile indexing
      end if

      do J = 1,IGND
        TBARACC_MO(I,J) = TBARACC_MO(I,J) + TBARROT(I,M,J) * FAREROT(I,M)
        TCTOACC_MO(I,J) = TCTOACC_MO(I,J) + TCTOROT(I,M,J) * FAREROT(I,M)
        TCBOACC_MO(I,J) = TCBOACC_MO(I,J) + TCBOROT(I,M,J) * FAREROT(I,M)
        ! Convert from m3/m3 to kg/m2
        THLQACC_MO(I,J) = THLQACC_MO(I,J) + THLQROT(I,M,J) * FAREROT(I,M) * 1000. * DLZWROT(I,M,J)
        THICACC_MO(I,J) = THICACC_MO(I,J) + THICROT(I,M,J) * FAREROT(I,M) * 1000. * DLZWROT(I,M,J)
        ! Find the total soil moisture content
        ! Add up each soil layers moisture and convert from m3/m3 to kg/m2
        MRSO_MO(I) = MRSO_MO(I) + (THLQROT(I,M,J) * FAREROT(I,M) + THICROT(I,M,J) * FAREROT(I,M)) * 1000. * DLZWROT(I,M,J)
        MRSOL_MO(I,J) = MRSOL_MO(I,J) + (THLQROT(I,M,J) * FAREROT(I,M) + THICROT(I,M,J) * FAREROT(I,M)) * 1000. * DLZWROT(I,M,J)
        TRANSPACC_MO(I) = TRANSPACC_MO(I) + QFCROT(I,M,J) * FAREROT(I,M)

      end do ! loop 823

    end do ! loop 821

    ! Check if the active layer has become more shallow or deepened.
    ACTLYR_MAX_MO(I) = max(ACTLYR_MAX_MO(I),ACTLYR_tmp)
    ACTLYR_MIN_MO(I) = min(ACTLYR_MIN_MO(I),ACTLYR_tmp)
    FTABLE_MAX_MO(I) = max(FTABLE_MAX_MO(I),FTABLE_tmp)
    FTABLE_MIN_MO(I) = min(FTABLE_MIN_MO(I),FTABLE_tmp)

    do NT = 1,NMON
      if (IDAY == monthend(NT + 1) .and. NCOUNT == NDAY) then
        IMONTH = NT
        NDMONTH = (monthend(NT + 1) - monthend(NT)) * NDAY

        ! These are presently not being outputted but the code is kept in place if the need arises.
        !             IF (FSINACC_MO(I)>0.0) THEN
        !                 ALVSACC_MO(I)=ALVSACC_MO(I)/(FSINACC_MO(I)*0.5)
        !                 ALIRACC_MO(I)=ALIRACC_MO(I)/(FSINACC_MO(I)*0.5)
        !             ELSE
        !                 ALVSACC_MO(I)=0.0
        !                 ALIRACC_MO(I)=0.0
        !             END IF

        ! Albedo is only counted when sun is above horizon so it uses its own counter.\

        if (altotcntr_m(i) > 0) then
          ALTOTACC_MO(I) = ALTOTACC_MO(I)/real(altotcntr_m(i))
          ALSNOACC_MO(I) = ALSNOACC_MO(I)/real(altotcntr_m(i))
        else
          ALTOTACC_MO(I) = 0.
          ALSNOACC_MO(I) = 0.
        end if

        FLUTACC_MO(I) = FLUTACC_MO(I)/real(NDMONTH)
        FSINACC_MO(I) = FSINACC_MO(I)/real(NDMONTH)
        FLINACC_MO(I) = FLINACC_MO(I)/real(NDMONTH)
        HFSACC_MO(I) = HFSACC_MO(I)/real(NDMONTH)
        QEVPACC_MO(I) = QEVPACC_MO(I)/real(NDMONTH)
        groundHeatFlux_MO(I) = groundHeatFlux_MO(I)/real(NDMONTH)
        SNOACC_MO(I) = SNOACC_MO(I)/real(NDMONTH)
        RHOSACC_MO(I) = RHOSACC_MO(I)/real(NDMONTH)
        SCAN_MO(I) = SCAN_MO(I)/real(NDMONTH)
        RCAN_MO(I) = RCAN_MO(I)/real(NDMONTH)
        FSNOACC_MO(I) = FSNOACC_MO(I)/real(NDMONTH)
        ZSNACC_MO(I) = ZSNACC_MO(I)/real(NDMONTH)
        WSNOACC_MO(I) = WSNOACC_MO(I)/real(NDMONTH)
        BCSNACC_MO(I) = BCSNACC_MO(I)/real(NDMONTH)
        TAACC_MO(I) = TAACC_MO(I)/real(NDMONTH)
        ACTLYR_MO(I) = ACTLYR_MO(I)/real(NDMONTH)
        FTABLE_MO(I) = FTABLE_MO(I)/real(NDMONTH)
        MRSO_MO(I) = MRSO_MO(I)/real(NDMONTH)
        OVRACC_MO(I) = OVRACC_MO(I) /real(NDMONTH)
        ROFBACC_MO(I) = ROFBACC_MO(I) /real(NDMONTH)
        ROFNACC_MO(I) = ROFNACC_MO(I) /real(NDMONTH)
        ROFACC_MO(I) = ROFACC_MO(I) /real(NDMONTH)
        PREACC_MO(I) = PREACC_MO(I) /real(NDMONTH)
        EVAPACC_MO(I) = EVAPACC_MO(I)/real(NDMONTH)
        TRANSPACC_MO(I) = TRANSPACC_MO(I)/real(NDMONTH)
        GROUNDEVAP(I) = GROUNDEVAP(I) /real(NDMONTH)
        CANOPYEVAP(I) = CANOPYEVAP(I) /real(NDMONTH)
        TSURFACC_MO(I) = TSURFACC_MO(I)/real(NDMONTH)
        TCANACC_MO(I) = TCANACC_MO(I)/real(NDMONTH)

        do J = 1,IGND
          TBARACC_MO(I,J) = TBARACC_MO(I,J)/real(NDMONTH)
          TCTOACC_MO(I,J) = TCTOACC_MO(I,J)/real(NDMONTH)
          TCBOACC_MO(I,J) = TCBOACC_MO(I,J)/real(NDMONTH)
          THLQACC_MO(I,J) = THLQACC_MO(I,J)/real(NDMONTH)
          THICACC_MO(I,J) = THICACC_MO(I,J)/real(NDMONTH)
          MRSOL_MO(I,J) = MRSOL_MO(I,J)/real(NDMONTH)
        end do

        FSSTAR_MO = FSINACC_MO(I) * (1. - ALTOTACC_MO(I))
        FLSTAR_MO = FLINACC_MO(I) - FLUTACC_MO(I)
        QH_MO = HFSACC_MO(I)
        QE_MO = QEVPACC_MO(I)

        tovere = 0.
        if (EVAPACC_MO(I) > 0.) tovere = TRANSPACC_MO(I)/EVAPACC_MO(I)

        ! Prepare the timestamp for this month (need in size 1 array)
        timeStamp = consecDays

        call writeOutput1D(lonLocalIndex,latLocalIndex,'scan_mo'   ,timeStamp,'scanopy', [SCAN_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rcan_mo'   ,timeStamp,'rcanopy', [RCAN_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snd_mo'   ,timeStamp,'snd', [ZSNACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsinacc_mo' ,timeStamp,'rsds', [FSINACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_mo' ,timeStamp,'rss', [FSSTAR_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_mo' ,timeStamp,'rls', [FLSTAR_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_mo'     ,timeStamp,'hfss', [QH_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_mo'     ,timeStamp,'hfls', [QE_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'hfg_mo'    ,timeStamp,'hfg', [groundHeatFlux_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snoacc_mo' ,timeStamp,'snw', [SNOACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rhosacc_mo' ,timeStamp,'snwdens', [RHOSACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsnoacc_mo' ,timeStamp,'snc', [FSNOACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_mo',timeStamp,'wsnw', [WSNOACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_mo'  ,timeStamp,'tas', [TAACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap',timeStamp,'evspsblsoi', [GROUNDEVAP(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap',timeStamp,'evspsblveg', [CANOPYEVAP(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'evspsbl_mo',timeStamp,'evspsbl', [CANOPYEVAP(I) + GROUNDEVAP(I) + TRANSPACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ovracc_mo' ,timeStamp,'mrros', [OVRACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofbacc_mo' ,timeStamp,'mrrob', [ROFBACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_mo' ,timeStamp,'mrro', [ROFACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_mo' ,timeStamp,'pr', [PREACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_mo',timeStamp,'evap', [EVAPACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_mo',timeStamp,'tran', [TRANSPACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'alsnoacc_mo',timeStamp,'albsn', [ALSNOACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_mo',timeStamp,'albs', [ALTOTACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tbaracc_mo',timeStamp,'tsl', [TBARACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tctopacc_mo',timeStamp,'tctop', [TCTOACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tcbotacc_mo',timeStamp,'tcbot', [TCBOACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thlqacc_mo',timeStamp,'mrsll', [THLQACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thicacc_mo',timeStamp,'mrsfl', [THICACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'mrso_mo',timeStamp,'mrso', [MRSO_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'mrsol_mo',timeStamp,'mrsol', [MRSOL_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_mo',timeStamp,'actlyr', [ACTLYR_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_max_mo',timeStamp,'actlyrmax', [ACTLYR_MAX_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_min_mo',timeStamp,'actlyrmin', [ACTLYR_MIN_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_mo',timeStamp,'ftable', [FTABLE_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_max_mo',timeStamp,'ftablemax', [FTABLE_MAX_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_min_mo',timeStamp,'ftablemin', [FTABLE_MIN_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'wue_mo' ,timeStamp,'wue', [tovere])
        ! tovere
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurf_mo'  ,timeStamp,'ts', [TSURFACC_MO(I)])
        if (TCANACC_MO(I) > 0.01) then
          TCN = TCANACC_MO(I) - TFREZ
        else
          TCN = 0.
        end if
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tcanacc_mo',timeStamp,'tcs', [TCN])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'tcanacc_mo'  ,timeStamp,'tcs', [TCANACC_MO(I)-TFREZ])

        call resetClassMon(nltest)

      end if ! if (IDAY==monthend(NT+1).AND.NCOUNT==NDAY)
    end do ! NMON

    end associate
  end subroutine class_monthly_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_class_annual_aw
  !> @{
  !> Accumulate and write out the annual physics outputs. These are kept in pointer structures as
  !! this subroutine is called each physics timestep and we increment the timestep values to produce annuals values.
  !! The pointer to the annual data structures (in classStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine class_annual_aw (lonLocalIndex, latLocalIndex, IDAY, realyr, NCOUNT, NDAY, &
                              nltest, nmtest, lastDOY)

    use classStateVars, only : class_out, resetClassYr, class_rot
    use classicParams,  only : nmon, monthend, nmos, ignd, SBC, DELT
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: IDAY
    integer, intent(in) :: realyr
    integer, intent(in) :: NCOUNT
    integer, intent(in) :: NDAY
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: lastDOY

    ! local
    integer :: i, m, j
    real :: tovere
    real :: FSSTAR_YR
    real :: FLSTAR_YR
    real :: QH_YR
    real :: QE_YR
    real :: ACTLYR_tmp
    real :: FTABLE_tmp
    real, dimension(1) :: timeStamp

    ! Associate names with variables defined in derived types.
    associate( &
    ALVSROT         => class_rot%ALVSROT,               & !< real, dimension(:,:):
    FAREROT         => class_rot%FAREROT,               & !< real, dimension(:,:):
    ALIRROT         => class_rot%ALIRROT,               & !< real, dimension(:,:):
    GTROT           => class_rot%GTROT,                 & !< real, dimension(:,:):
    HFSROT          => class_rot%HFSROT,                & !< real, dimension(:,:):
    QEVPROT         => class_rot%QEVPROT,               & !< real, dimension(:,:):
    ROFROT          => class_rot%ROFROT,                & !< real, dimension(:,:):
    QFSROT          => class_rot%QFSROT,                & !< real, dimension(:,:):
    QFCROT          => class_rot%QFCROT,                & !< real, dimension(:,:,:):
    FSGVROT         => class_rot%FSGVROT,               & !< real, dimension(:,:): Diagnosed net shortwave radiation on vegetation canopy
    FSGSROT         => class_rot%FSGSROT,               & !< real, dimension(:,:): Diagnosed net shortwave radiation on ground snow surface
    FSGGROT         => class_rot%FSGGROT,               & !< real, dimension(:,:): Diagnosed net shortwave radiation on ground surface
    FSSROW          => class_rot%FSSROW,                & !< real, dimension(:):
    FDLROW          => class_rot%FDLROW,                & !< real, dimension(:):
    FSVHROW         => class_rot%FSVHROW,               & !< real, dimension(:):
    FSIHROW         => class_rot%FSIHROW,               & !< real, dimension(:):
    TAROW           => class_rot%TAROW,                 & !< real, dimension(:):
    PREROW          => class_rot%PREROW,                & !< real, dimension(:):
    ftable          => class_rot%ftable,                & !< real, dimension(:,:)  : Depth to frozen water table (m)
    actlyr          => class_rot%actlyr,                & !< real, dimension(:,:)  : Active layer depth (m)
    dlzwrot         => class_rot%dlzwrot,               & !< real, dimension(:,:,:) : Permeable thickness of soil layer [m]
    THLQROT         => class_rot%THLQROT,               & !< real, dimension(:,:,:):
    THICROT         => class_rot%THICROT,               & !< real, dimension(:,:,:):
    ALVSACC_YR        => class_out%ALVSACC_YR,          & !< real, dimension(:) :
    ALIRACC_YR        => class_out%ALIRACC_YR,          & !< real, dimension(:) :
    FLUTACC_YR        => class_out%FLUTACC_YR,          & !< real, dimension(:) :
    FSINACC_YR        => class_out%FSINACC_YR,          & !< real, dimension(:) :
    FLINACC_YR        => class_out%FLINACC_YR,          & !< real, dimension(:) :
    HFSACC_YR         => class_out%HFSACC_YR,           & !< real, dimension(:) :
    QEVPACC_YR        => class_out%QEVPACC_YR,          & !< real, dimension(:) :
    ROFACC_YR         => class_out%ROFACC_YR,           & !< real, dimension(:) :
    PREACC_YR         => class_out%PREACC_YR,           & !< real, dimension(:) :
    EVAPACC_YR        => class_out%EVAPACC_YR,          & !< real, dimension(:) :
    TRANSPACC_YR      => class_out%TRANSPACC_YR,        & !< real, dimension(:) :
    TAACC_YR          => class_out%TAACC_YR,            & !< real, dimension(:) :
    THLQACC_YR        => class_out%THLQACC_YR,          & !< real, dimension(:,:) : Volumetric liquid water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)
    THICACC_YR        => class_out%THICACC_YR,          & !< real, dimension(:,:) : Volumetric frozen water content of soil layers \f$[kg m^{-2}]\f$ (accumulated for means)
    ACTLYR_YR         => class_out%ACTLYR_YR,           & !< real, dimension(:) :
    ACTLYR_MIN_YR     => class_out%ACTLYR_MIN_YR,       & !< real, dimension(:) :
    ACTLYR_MAX_YR     => class_out%ACTLYR_MAX_YR,       & !< real, dimension(:) :
    FTABLE_YR         => class_out%FTABLE_YR,           & !< real, dimension(:) :
    FTABLE_MIN_YR     => class_out%FTABLE_MIN_YR,       & !< real, dimension(:) :
    FTABLE_MAX_YR     => class_out%FTABLE_MAX_YR,       & !< real, dimension(:) :
    ALTOTACC_YR       => class_out%ALTOTACC_YR,         & !< real, dimension(:) :
    altotcntr_yr      => class_out%altotcntr_yr         & !< integer, dimension(:) :
    )

    !> Accumulate output data for yearly averaged fields for class grid-mean.
    !> for both parallel mode and stand alone mode
    FSSTAR_YR   = 0.0
    FLSTAR_YR   = 0.0
    QH_YR       = 0.0
    QE_YR       = 0.0
    ACTLYR_tmp  = 0.0
    FTABLE_tmp  = 0.0

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do M = 1,NMTEST

      ! These are presently not being outputted but the code is kept in place if the need arises.
      !         ALVSACC_YR(I)=ALVSACC_YR(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
      !         ALIRACC_YR(I)=ALIRACC_YR(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)

      FLUTACC_YR(I) = FLUTACC_YR(I) + SBC * GTROT(I,M) ** 4 * FAREROT(I,M)
      FSINACC_YR(I) = FSINACC_YR(I) + FSSROW(I) * FAREROT(I,M)
      FLINACC_YR(I) = FLINACC_YR(I) + FDLROW(I) * FAREROT(I,M)
      HFSACC_YR(I) = HFSACC_YR(I) + HFSROT(I,M) * FAREROT(I,M)
      QEVPACC_YR(I) = QEVPACC_YR(I) + QEVPROT(I,M) * FAREROT(I,M)
      TAACC_YR(I) = TAACC_YR(I) + TAROW(I) * FAREROT(I,M)
      ROFACC_YR(I) = ROFACC_YR(I) + ROFROT(I,M) * FAREROT(I,M)
      PREACC_YR(I) = PREACC_YR(I) + PREROW(I) * FAREROT(I,M)
      EVAPACC_YR(I) = EVAPACC_YR(I) + QFSROT(I,M) * FAREROT(I,M)
      ACTLYR_YR(I) = ACTLYR_YR(I) + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_YR(I) = FTABLE_YR(I) + FTABLE(I,M) * FAREROT(I,M)
      ACTLYR_TMP = ACTLYR_TMP + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_TMP = FTABLE_TMP + FTABLE(I,M) * FAREROT(I,M)

      do J = 1,IGND
        TRANSPACC_YR(I) = TRANSPACC_YR(I) + QFCROT(I,M,J) * FAREROT(I,M)
        ! Convert from m3/m3 to kg/m2
        THLQACC_YR(I,J) = THLQACC_YR(I,J) + THLQROT(I,M,J) * FAREROT(I,M) * 1000. * DLZWROT(I,M,J)
        THICACC_YR(I,J) = THICACC_YR(I,J) + THICROT(I,M,J) * FAREROT(I,M) * 1000. * DLZWROT(I,M,J)
      end do

      if (FSSROW(I) > 0.0) then
        ALTOTACC_YR(I) = ALTOTACC_YR(I) + ((FSSROW(I) - (FSGVROT(I,M) + FSGSROT(I,M) + FSGGROT(I,M))) &
                         /FSSROW(I) ) * FAREROT(I,M)
        if (M == 1) altotcntr_yr(i) = altotcntr_yr(i) + 1 ! s.r.c. added if statement to prevent per-tile indexing
      end if

    end do ! loop 828

    ! Check if the active layer has become more shallow or deepened.
    ACTLYR_MAX_YR(I) = max(ACTLYR_MAX_YR(I),ACTLYR_tmp)
    ACTLYR_MIN_YR(I) = min(ACTLYR_MIN_YR(I),ACTLYR_tmp)
    FTABLE_MAX_YR(I) = max(FTABLE_MAX_YR(I),FTABLE_tmp)
    FTABLE_MIN_YR(I) = min(FTABLE_MIN_YR(I),FTABLE_tmp)

    if (IDAY == lastDOY .and. NCOUNT == NDAY) then

      ! These are presently not being outputted but the code is kept in place if the need arises.
      !             IF (FSINACC_YR(I)>0.0) THEN
      !                 ALVSACC_YR(I)=ALVSACC_YR(I)/(FSINACC_YR(I)*0.5)
      !                 ALIRACC_YR(I)=ALIRACC_YR(I)/(FSINACC_YR(I)*0.5)
      !             ELSE
      !                 ALVSACC_YR(I)=0.0
      !                 ALIRACC_YR(I)=0.0
      !             END IF

      FLUTACC_YR(I) = FLUTACC_YR(I)/(real(NDAY) * real(lastDOY))
      FSINACC_YR(I) = FSINACC_YR(I)/(real(NDAY) * real(lastDOY))
      FLINACC_YR(I) = FLINACC_YR(I)/(real(NDAY) * real(lastDOY))
      HFSACC_YR(I) = HFSACC_YR(I)/(real(NDAY) * real(lastDOY))
      QEVPACC_YR(I) = QEVPACC_YR(I)/(real(NDAY) * real(lastDOY))
      ROFACC_YR(I) = ROFACC_YR(I)/(real(NDAY) * real(lastDOY))
      PREACC_YR(I) = PREACC_YR(I)/(real(NDAY) * real(lastDOY))
      EVAPACC_YR(I) = EVAPACC_YR(I)/(real(NDAY) * real(lastDOY))
      TRANSPACC_YR(I) = TRANSPACC_YR(I)/(real(NDAY) * real(lastDOY))
      TAACC_YR(I) = TAACC_YR(I)/(real(NDAY) * real(lastDOY))
      ACTLYR_YR(I) = ACTLYR_YR(I)/(real(NDAY) * real(lastDOY))
      FTABLE_YR(I) = FTABLE_YR(I)/(real(NDAY) * real(lastDOY))

      do J = 1,IGND
        THLQACC_YR(I,J) = THLQACC_YR(I,J)/(real(NDAY) * real(lastDOY))
        THICACC_YR(I,J) = THICACC_YR(I,J)/(real(NDAY) * real(lastDOY))
      end do

      ! Albedo is only counted when sun is above horizon so it uses its own counter.
      if (altotcntr_yr(i) > 0) ALTOTACC_YR(I) = ALTOTACC_YR(I)/(real(altotcntr_yr(i)))

      FSSTAR_YR = FSINACC_YR(I) * (1. - ALTOTACC_YR(I))
      FLSTAR_YR = FLINACC_YR(I) - FLUTACC_YR(I)
      QH_YR = HFSACC_YR(I)
      QE_YR = QEVPACC_YR(I)

      tovere = 0.
      if (EVAPACC_YR(I) > 0.) tovere = TRANSPACC_YR(I)/EVAPACC_YR(I)

      ! Prepare the timestamp for this year
      timeStamp = consecDays

      call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_yr' ,timeStamp,'rss', [FSSTAR_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_yr' ,timeStamp,'rls', [FLSTAR_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_yr'     ,timeStamp,'hfss', [QH_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_yr'     ,timeStamp,'hfls', [QE_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_yr' ,timeStamp,'mrro', [ROFACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_yr' ,timeStamp,'pr', [PREACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_yr' ,timeStamp,'tas', [TAACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_yr' ,timeStamp,'evspsbl', [EVAPACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_yr' ,timeStamp,'tran', [TRANSPACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_yr' ,timeStamp,'albs', [ALTOTACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_yr' ,timeStamp,'actlyr', [actlyr_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_yr' ,timeStamp,'ftable', [ftable_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_max_yr' ,timeStamp,'actlyrmax', [actlyr_max_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_max_yr' ,timeStamp,'ftablemax', [ftable_max_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_min_yr' ,timeStamp,'actlyrmin', [actlyr_min_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_min_yr' ,timeStamp,'ftablemin', [ftable_min_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'thlqacc_yr',timeStamp,'mrsll', [THLQACC_YR(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'thicacc_yr',timeStamp,'mrsfl', [THICACC_YR(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'wue_yr' ,timeStamp,'wue', [tovere])
      !> ADD INITIALIZTION FOR YEARLY ACCUMULATED ARRAYS
      call resetClassYr(nltest)

    end if !> IDAY==365/366 .AND. NDAY

    end associate
  end subroutine class_annual_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_convertUnitsCTEM
  !! @{
  !> Do some unit conversions for CTEM (biogeochemical processes) so they are ready to be written out
  subroutine convertUnitsCTEM (nltest, nmtest)

    use classicParams, only : icc, ignd, iccp1, wtCH4, convertkgC, convertkgN, iccp2
    use ctemStateVars, only : vrot, c_switch

    implicit none

    integer, intent(in) :: nltest   ! number of grid cells (offline = 1)
    integer, intent(in) :: nmtest   ! number of tiles per cell

    integer :: i, m, j, k

    associate( &
    Ncycle_on             => c_switch%Ncycle_on,              & !< logical:
    fcancmxrow        => vrot%fcancmx,                  & !< real, dimension(:,:,:) :
    gppvegrow         => vrot%gppveg,                   & !< real, dimension(:,:,:) :
    nepvegrow         => vrot%nepveg,                   & !< real, dimension(:,:,:) :
    nbpvegrow         => vrot%nbpveg,                   & !< real, dimension(:,:,:) :
    nppvegrow         => vrot%nppveg,                   & !< real, dimension(:,:,:) :
    hetroresvegrow    => vrot%hetroresveg,              & !< real, dimension(:,:,:) :
    autoresvegrow     => vrot%autoresveg,               & !< real, dimension(:,:,:) :
    rmrvegrow         => vrot%rmrveg,                   & !< real, dimension(:,:,:) :
    rgvegrow          => vrot%rgveg,                    & !< real, dimension(:,:,:) :
    litresvegrow      => vrot%litresveg,                & !< real, dimension(:,:,:,:) :
    soilcresvegrow    => vrot%soilcresveg,              & !< real, dimension(:,:,:,:) :
    litrfallvegrow    => vrot%litrfallveg,              & !< real, dimension(:,:,:,:) :litterfall for each pft
    humiftrsvegrow    => vrot%humiftrsveg,              & !< real, dimension(:,:,:,:) : Transfer of humidified litter from litter to soil C pool 
    npprow            => vrot%npp,                      & !< real, dimension(:,:) :
    neprow            => vrot%nep,                      & !< real, dimension(:,:) :
    nepCMIProw        => vrot%nepCMIP,                  & !< real, dimension(:,:) :
    nbprow            => vrot%nbp,                      & !< real, dimension(:,:) :
    gpprow            => vrot%gpp,                      & !< real, dimension(:,:) :
    lucemcomrow       => vrot%lucemcom,                 & !< real, dimension(:,:) :
    lucltrinrow       => vrot%lucltrin,                 & !< real, dimension(:,:) :
    lucsocinrow       => vrot%lucsocin,                 & !< real, dimension(:,:) :
    lucemcomnrow       => vrot%lucemcomn,                 & !< real, dimension(:,:) :
    lucltrinnrow       => vrot%lucltrinn,                 & !< real, dimension(:,:) :
    lucsocinnrow       => vrot%lucsocinn,                 & !< real, dimension(:,:) :
    hetroresrow       => vrot%hetrores,                 & !< real, dimension(:,:) :
    autoresrow        => vrot%autores,                  & !< real, dimension(:,:) :
    litresrow         => vrot%litres,                   & !< real, dimension(:,:) :
    socresrow         => vrot%socres,                   & !< real, dimension(:,:) :
    ch4WetSpecrow        => vrot%ch4WetSpec,            & !< real, dimension(:,:) :
    ch4WetDynrow        => vrot%ch4WetDyn,              & !< real, dimension(:,:) :
    ch4soillsrow      => vrot%ch4_soills,               & !< real, dimension(:,:) :
    nppmossrow         => vrot%nppmoss,                 & !< real, dimension(:,:) :
    armossrow          => vrot%armoss,                  & !< real, dimension(:,:) :
    emit_co2row       => vrot%emit_co2                  & !< real, dimension(:,:,:) :
    )

    !> Some unit conversions:
    !!
    !! We want to go from umol CO2/m2/s to kg C/m2/s so:
    !! umolCO2/m2/s * mol/10^6umol * mol C/ molCO2 * 12.01 g C / mol C * 1 kg/ 1000g = kgC/m2/s
    !! umolCO2/m2/s * 1.201E-8 = kgC/m2/s
    !! convertkgC = 1.201E-8

    do i = 1,nltest
      do m = 1,nmtest

        nppmossrow(i,m) = nppmossrow(i,m) * convertkgC
        armossrow(i,m) = armossrow(i,m) * convertkgC

        do j = 1,icc
          if (fcancmxrow(i,m,j) > 0.0) then

            gppvegrow(i,m,j) = gppvegrow(i,m,j) * convertkgC
            nppvegrow(i,m,j) = nppvegrow(i,m,j) * convertkgC
            nepvegrow(i,m,j) = nepvegrow(i,m,j) * convertkgC
            nbpvegrow(i,m,j) = nbpvegrow(i,m,j) * convertkgC
            litrfallvegrow(i,m,j) = litrfallvegrow(i,m,j) * convertkgC
            
            hetroresvegrow(i,m,j) = hetroresvegrow(i,m,j) * convertkgC
            autoresvegrow(i,m,j) = autoresvegrow(i,m,j) * convertkgC
            rmrvegrow(i,m,j) = rmrvegrow(i,m,j) * convertkgC
            rgvegrow(i,m,j) = rgvegrow(i,m,j) * convertkgC
            do k = 1,ignd
              litresvegrow(i,m,j,k)=litresvegrow(i,m,j,k)*convertkgC
              soilcresvegrow(i,m,j,k)=soilcresvegrow(i,m,j,k)*convertkgC
              humiftrsvegrow(i,m,j,k) = humiftrsvegrow(i,m,j,k) * convertkgC
            end do



            ! emit_co2,like all fire gas fluxes is in kg {species} / m2 / s. So we need
            ! to convert from kg CO2/m2/s to kg C/m2/s. Since 1 g C = 0.083 mole CO2 = 3.664 g CO2
            ! kg CO2/m2/s * 1 g C/ 3.664 g CO2 = kg C /m2/s
            emit_co2row(i,m,j) = emit_co2row(i,m,j) / 3.664

          end if
        end do ! loop 30 ! icc

        !> Now for the bare fraction of the grid cell and the LUC products pool.
        hetroresvegrow(i,m,iccp1) = hetroresvegrow(i,m,iccp1) * convertkgC
        nepvegrow(i,m,iccp1) = nepvegrow(i,m,iccp1) * convertkgC
        nbpvegrow(i,m,iccp1) = nbpvegrow(i,m,iccp1) * convertkgC
        do k = 1,ignd
          litresvegrow(i,m,iccp1:iccp2,k)=litresvegrow(i,m,iccp1:iccp2,k)*convertkgC
          soilcresvegrow(i,m,iccp1:iccp2,k)=soilcresvegrow(i,m,iccp1:iccp2,k)*convertkgC
        end do
        npprow(i,m)     = npprow(i,m) * convertkgC
        gpprow(i,m)     = gpprow(i,m) * convertkgC
        neprow(i,m)     = neprow(i,m) * convertkgC
        nepCMIProw(i,m)     = nepCMIProw(i,m) * convertkgC
        nbprow(i,m)     = nbprow(i,m) * convertkgC
        lucemcomrow(i,m) = lucemcomrow(i,m) * convertkgC
        lucltrinrow(i,m) = lucltrinrow(i,m) * convertkgC
        lucsocinrow(i,m) = lucsocinrow(i,m) * convertkgC
        hetroresrow(i,m) = hetroresrow(i,m) * convertkgC
        autoresrow(i,m) = autoresrow(i,m) * convertkgC
        litresrow(i,m)  = litresrow(i,m) * convertkgC
        socresrow(i,m)  = socresrow(i,m) * convertkgC

        if(Ncycle_on) then
          lucemcomnrow(i,m) = lucemcomnrow(i,m) * convertkgN
          lucltrinnrow(i,m) = lucltrinnrow(i,m) * convertkgN
          lucsocinnrow(i,m) = lucsocinnrow(i,m) * convertkgN
        end if

        ch4WetSpecrow(i,m) = ch4WetSpecrow(i,m) * convertkgC * wtCH4 / 12.01 ! convert from umolch4/m2/s to kg CH4/ m2 /s
        ch4WetDynrow(i,m) = ch4WetDynrow(i,m) * convertkgC * wtCH4 / 12.01 ! convert from umolch4/m2/s to kg CH4/ m2 /s
        ch4soillsrow(i,m) = ch4soillsrow(i,m) * convertkgC * wtCH4 / 12.01 ! convert from umolch4/m2/s to kg CH4/ m2 /s

      end do ! loop 20
    end do ! loop 10

    end associate
  end subroutine convertUnitsCTEM
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_ctem_daily_aw
  !> @{
  !> Accumulate and write the daily biogeochemical outputs
  !! @author J. Melton

  subroutine ctem_daily_aw (lonLocalIndex, latLocalIndex, nltest, nmtest, iday, ncount, nday, realyr, grclarea, ipeatlandrow)

    ! J. Melton Feb 2016.

    use classStateVars, only : class_rot
    use ctemStateVars,  only : vrot, c_switch !, resetdaily, ctem_grd ctem_tile,
    use classicParams,  only : icc, ignd, nmos, iccp1, wtCH4, convertkgC, iccp2, convertg2kg, convertkgN
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    integer, intent(in) :: ncount
    integer, intent(in) :: nday
    integer, intent(in) :: realyr
    real, intent(in), dimension(:) :: grclarea ! use pointer FLAG
    integer, intent(in), dimension(:,:) :: ipeatlandrow ! use pointer FLAG

    ! local vars
    real, dimension(nltest) :: gpp_g      !<
    real, dimension(nltest) :: npp_g      !<
    real, dimension(nltest) :: nbp_g      !<
    real, dimension(nltest) :: socres_g   !<
    real, dimension(nltest) :: autores_g  !<
    real, dimension(nltest) :: litres_g   !<
    real, dimension(nltest) :: dstcemls3_g !<
    real, dimension(nltest) :: litrfall_g !<
    real, dimension(nltest) :: rml_g      !<
    real, dimension(nltest) :: rms_g      !<
    real, dimension(nltest) :: rg_g       !<
    real, dimension(nltest) :: leaflitr_g !<
    real, dimension(nltest) :: tltrstem_g !<
    real, dimension(nltest) :: tltrroot_g !<
    real, dimension(nltest) :: nep_g      !<
    real, dimension(nltest) :: nepCMIP_g      !<
    real, dimension(nltest) :: hetrores_g !<
    real, dimension(nltest) :: dstcemls_g !<
    real, dimension(nltest) :: humiftrs_g !<
    real, dimension(nltest) :: rmr_g      !<
    real, dimension(nltest) :: tltrleaf_g !<
    real, dimension(nltest) :: gavgltms_g !<
    real, dimension(nltest) :: vgbiomas_g !<
    real, dimension(nltest) :: gavglai_g  !<
    real, dimension(nltest) :: gavgscms_g !<

    real, dimension(nltest) :: leafns2s_g !<
    real, dimension(nltest) :: stemns2s_g !<
    real, dimension(nltest) :: rootns2s_g !<
    real, dimension(nltest) :: re_alloc_s2l_g  !<
    real, dimension(nltest) :: re_alloc_r2l_g  !<
    real, dimension(nltest) :: re_alloc_sr2l_g !<
    real, dimension(nltest) :: gleafmas_g !<
    real, dimension(nltest) :: gleafmas_ns_g !<
    real, dimension(nltest) :: gleafmas_s_g  !<
    real, dimension(nltest) :: bleafmas_g !<
    real, dimension(nltest) :: stemmass_g !<
    real, dimension(nltest) :: stemmass_ns_g !<
    real, dimension(nltest) :: stemmass_s_g  !<
    real, dimension(nltest) :: rootmass_g !<
    real, dimension(nltest) :: rootmass_ns_g !<
    real, dimension(nltest) :: rootmass_s_g  !<
    real, dimension(nltest) :: litrmass_g !<
    real, dimension(nltest) :: soilcmas_g !<
    real, dimension(nltest) :: nh4_mass_g      !<
    real, dimension(nltest) :: no3_mass_g      !<
    real, dimension(nltest) :: ngleafmas_g     !<
    real, dimension(nltest) :: ngleafmas_ns_g  !<
    real, dimension(nltest) :: ngleafmas_s_g   !<
    real, dimension(nltest) :: nbleafmas_g     !<
    real, dimension(nltest) :: nstemmass_g     !<
    real, dimension(nltest) :: nstemmass_ns_g  !<
    real, dimension(nltest) :: nstemmass_s_g   !<
    real, dimension(nltest) :: nrootmass_g     !<
    real, dimension(nltest) :: nrootmass_ns_g  !<
    real, dimension(nltest) :: nrootmass_s_g   !<
    real, dimension(nltest) :: c2n_l_g         !<
    real, dimension(nltest) :: c2n_s_g         !<
    real, dimension(nltest) :: c2n_r_g         !<
    real, dimension(nltest) :: c2n_wp_g        !<
    real, dimension(nltest) :: c2n_litr_g      !<
    real, dimension(nltest) :: c2n_humus_g     !<
    real, dimension(nltest) :: nlitrmass_g     !<
    real, dimension(nltest) :: soilnmas_g      !<
    real, dimension(nltest) :: nvgbiomas_g     !<
    real, dimension(nltest) :: slai_g     !<
    real, dimension(nltest) :: ailcg_g    !<
    real, dimension(nltest) :: ailcb_g    !<
    real, dimension(nltest) :: veghght_g  !<
    real, dimension(nltest) :: rootdpth_g !<
    real, dimension(nltest) :: roottemp_g !<
    real, dimension(nltest) :: totcmass_g !<
    real, dimension(nltest) :: tcanoacc_out_g !<
    real, dimension(nltest) :: burnfrac_g !<
    real, dimension(nltest) :: smfuncveg_g !<
    real, dimension(nltest) :: lucemcom_g !<
    real, dimension(nltest) :: lucltrin_g !<
    real, dimension(nltest) :: lucsocin_g !<
    real, dimension(nltest) :: lucemcomn_g !<
    real, dimension(nltest) :: lucltrinn_g !<
    real, dimension(nltest) :: lucsocinn_g !<
    real, dimension(nltest) :: emit_co2_g !<
    real, dimension(nltest) :: emit_ch4_g !<
    real, dimension(nltest) :: ch4WetSpec_g  !<
    real, dimension(nltest) :: wetfdyn_g  !<
    real, dimension(nltest) :: ch4WetDyn_g  !<
    real, dimension(nltest) :: ch4soills_g   !<
    real, dimension(nltest,icc) :: afrleaf_g  !<
    real, dimension(nltest,icc) :: afrstem_g  !<
    real, dimension(nltest,icc) :: afrroot_g  !<
    real, dimension(nltest,icc) :: rmlvegrow_g !<
    real, dimension(nltest,icc) :: anvegrow_g !<
    real, dimension(nltest,ignd) :: rmatctem_g !<

    real, dimension(nltest)   :: bnf_tot_g     !<
    real, dimension(nltest)   :: bnf_free_g    !<
    real, dimension(nltest)   :: bnf_ant_g     !<
    real, dimension(nltest)   :: bnf_nat_g     !<
    real, dimension(nltest)   :: nitrif_g      !<
    real, dimension(nltest)   :: no_nit_g      !<
    real, dimension(nltest)   :: no_denit_g    !<
    real, dimension(nltest)   :: no_nitdenit_g !<
    real, dimension(nltest)   :: n2o_nit_g     !<
    real, dimension(nltest)   :: n2o_denit_g   !<
    real, dimension(nltest)   :: n2o_nitdenit_g!<
    real, dimension(nltest)   :: n2_denit_g    !<
    real, dimension(nltest)   :: nvol_g        !<
    real, dimension(nltest)   :: nleach_g      !<
    real, dimension(nltest)   :: appl_fert_g   !<
    real, dimension(nltest)   :: ndep_nh4_g    !<
    real, dimension(nltest)   :: ndep_no3_g    !<
    real, dimension(nltest)   :: ndemand_wp_npp_g     !<
    real, dimension(nltest)   :: nuptake_p_nh4_g      !<
    real, dimension(nltest)   :: nuptake_p_no3_g      !<
    real, dimension(nltest)   :: nuptake_a_actl_nh4_g !<
    real, dimension(nltest)   :: nuptake_a_actl_no3_g !<
    real, dimension(nltest)   :: nalloc_l_g      !<
    real, dimension(nltest)   :: nalloc_s_g      !<
    real, dimension(nltest)   :: nalloc_r_g      !<
    real, dimension(nltest)   :: nresorped_s_g   !<
    real, dimension(nltest)   :: nresorped_r_g   !<
    real, dimension(nltest)   :: nre_alloc_s2l_g !<
    real, dimension(nltest)   :: nre_alloc_r2l_g !<
    real, dimension(nltest)   :: nleafns2s_g     !<
    real, dimension(nltest)   :: nstemns2s_g     !<
    real, dimension(nltest)   :: nrootns2s_g     !<
    real, dimension(nltest)   :: nlitr_l_g       !<
    real, dimension(nltest)   :: nlitr_s_g       !<
    real, dimension(nltest)   :: nlitr_r_g       !<
    real, dimension(nltest)   :: gl2bl_grass_nflux_g !<
    real, dimension(nltest)   :: nhumtrs_g       !<
    real, dimension(nltest)   :: nmineral_litr_g !<
    real, dimension(nltest)   :: nmineral_humus_g!<
    real, dimension(nltest)   :: nimmobil_nh4_g  !<
    real, dimension(nltest)   :: nimmobil_no3_g  !<
    real, dimension(nltest)   :: fNnetland_g     !<

    real, dimension(nltest,nmtest) :: leaflitr_t !<
    real, dimension(nltest,nmtest) :: tltrleaf_t !<
    real, dimension(nltest,nmtest) :: tltrstem_t !<
    real, dimension(nltest,nmtest) :: tltrroot_t !<
    real, dimension(nltest,nmtest) :: ailcg_t    !<
    real, dimension(nltest,nmtest) :: ailcb_t    !<
    real, dimension(nltest,nmtest,ignd) :: rmatctem_t !< nlat,nmos,ignd
    real, dimension(nltest,nmtest) :: veghght_t  !<
    real, dimension(nltest,nmtest) :: rootdpth_t !<
    real, dimension(nltest,nmtest) :: roottemp_t !<
    real, dimension(nltest,nmtest) :: slai_t     !<
    real, dimension(nltest,nmtest) :: afrroot_t  !<
    real, dimension(nltest,nmtest) :: afrleaf_t  !<
    real, dimension(nltest,nmtest) :: afrstem_t  !<
    real, dimension(nltest,nmtest) :: laimaxg_t  !<
    real, dimension(nltest,nmtest) :: stemmass_t !<
    real, dimension(nltest,nmtest) :: stemmass_ns_t !<
    real, dimension(nltest,nmtest) :: stemmass_s_t  !<
    real, dimension(nltest,nmtest) :: rootmass_t !<
    real, dimension(nltest,nmtest) :: rootmass_ns_t !<
    real, dimension(nltest,nmtest) :: rootmass_s_t  !<
    real, dimension(nltest,nmtest,ignd) :: litrmass_t !<
    real, dimension(nltest,nmtest,ignd) :: soilcmas_t !<
    real, dimension(nltest,nmtest) :: gleafmas_t !<
    real, dimension(nltest,nmtest) :: gleafmas_ns_t !<
    real, dimension(nltest,nmtest) :: gleafmas_s_t  !<
    real, dimension(nltest,nmtest) :: bleafmas_t !<

    real, dimension(nltest,nmtest) :: nh4_mass_t       !<
    real, dimension(nltest,nmtest) :: no3_mass_t       !<


    real, dimension(nltest,nmtest) :: ngleafmas_t      !<
    real, dimension(nltest,nmtest) :: ngleafmas_ns_t   !<
    real, dimension(nltest,nmtest) :: ngleafmas_s_t    !<
    real, dimension(nltest,nmtest) :: nbleafmas_t      !<
    real, dimension(nltest,nmtest) :: nstemmass_t      !<
    real, dimension(nltest,nmtest) :: nstemmass_ns_t   !<
    real, dimension(nltest,nmtest) :: nstemmass_s_t    !<
    real, dimension(nltest,nmtest) :: nrootmass_t      !<
    real, dimension(nltest,nmtest) :: nrootmass_ns_t   !<
    real, dimension(nltest,nmtest) :: nrootmass_s_t    !<
    real, dimension(nltest,nmtest) :: nvgbiomas_t      !<
    real, dimension(nltest,nmtest) :: c2n_l_t          !<
    real, dimension(nltest,nmtest) :: c2n_s_t          !<
    real, dimension(nltest,nmtest) :: c2n_r_t          !<
    real, dimension(nltest,nmtest) :: c2n_wp_t         !<
    real, dimension(nltest,nmtest) :: c2n_litr_t       !<
    real, dimension(nltest,nmtest) :: c2n_humus_t      !<
    real, dimension(nltest,nmtest) :: nlitrmass_t      !<
    real, dimension(nltest,nmtest) :: soilnmas_t       !<
    real, dimension(nltest,nmtest) :: emit_co2_t !<
    real, dimension(nltest,nmtest) :: emit_ch4_t !<
    real, dimension(nltest,nmtest) :: smfuncveg_t !<
    real, dimension(nltest,nmtest) :: vgbiomas_t       !<
    real, dimension(nltest,nmtest) :: leafns2s_t       !<
    real, dimension(nltest,nmtest) :: stemns2s_t       !<
    real, dimension(nltest,nmtest) :: rootns2s_t       !<
    real, dimension(nltest,nmtest) :: re_alloc_s2l_t   !<
    real, dimension(nltest,nmtest) :: re_alloc_r2l_t   !<
    real, dimension(nltest,nmtest) :: re_alloc_sr2l_t  !<
    real, dimension(nltest,nmtest) :: bnf_tot_t        !<
    real, dimension(nltest,nmtest) :: bnf_free_t       !<
    real, dimension(nltest,nmtest) :: bnf_ant_t        !<
    real, dimension(nltest,nmtest) :: bnf_nat_t        !<
    real, dimension(nltest,nmtest) :: nitrif_t         !<
    real, dimension(nltest,nmtest) :: no_nit_t         !<
    real, dimension(nltest,nmtest) :: no_denit_t       !<
    real, dimension(nltest,nmtest) :: no_nitdenit_t    !<
    real, dimension(nltest,nmtest) :: n2o_nit_t        !<
    real, dimension(nltest,nmtest) :: n2o_denit_t      !<
    real, dimension(nltest,nmtest) :: n2o_nitdenit_t   !<
    real, dimension(nltest,nmtest) :: n2_denit_t       !<
    real, dimension(nltest,nmtest) :: nvol_t           !<
    real, dimension(nltest,nmtest) :: nleach_t         !<
    real, dimension(nltest,nmtest) :: appl_fert_t      !<
    real, dimension(nltest,nmtest) :: ndep_nh4_t       !<
    real, dimension(nltest,nmtest) :: ndep_no3_t       !<
    real, dimension(nltest,nmtest) :: ndemand_wp_npp_t     !<
    real, dimension(nltest,nmtest) :: nuptake_p_nh4_t      !<
    real, dimension(nltest,nmtest) :: nuptake_p_no3_t      !<
    real, dimension(nltest,nmtest) :: nuptake_a_actl_nh4_t !<
    real, dimension(nltest,nmtest) :: nuptake_a_actl_no3_t !<
    real, dimension(nltest,nmtest) :: nalloc_l_t       !<
    real, dimension(nltest,nmtest) :: nalloc_s_t       !<
    real, dimension(nltest,nmtest) :: nalloc_r_t       !<
    real, dimension(nltest,nmtest) :: nresorped_s_t    !<
    real, dimension(nltest,nmtest) :: nresorped_r_t    !<
    real, dimension(nltest,nmtest) :: nre_alloc_s2l_t  !<
    real, dimension(nltest,nmtest) :: nre_alloc_r2l_t  !<
    real, dimension(nltest,nmtest) :: nleafns2s_t      !<
    real, dimension(nltest,nmtest) :: nstemns2s_t      !<
    real, dimension(nltest,nmtest) :: nrootns2s_t      !<
    real, dimension(nltest,nmtest) :: nlitr_l_t        !<
    real, dimension(nltest,nmtest) :: nlitr_s_t        !<
    real, dimension(nltest,nmtest) :: nlitr_r_t        !<
    real, dimension(nltest,nmtest) :: gl2bl_grass_nflux_t !<
    real, dimension(nltest,nmtest) :: nhumtrs_t        !<
    real, dimension(nltest,nmtest) :: nmineral_litr_t  !<
    real, dimension(nltest,nmtest) :: nmineral_humus_t !<
    real, dimension(nltest,nmtest) :: nimmobil_nh4_t   !<
    real, dimension(nltest,nmtest) :: nimmobil_no3_t   !<
    real, dimension(nltest,nmtest) :: fNnetland_t      !<

    integer :: i, m, j, nt, k
    real :: barefrac
    real :: sumfare
    real, dimension(1) :: timeStamp

    ! Associate names with variables defined in derived types.

    associate( &
    dofire                => c_switch%dofire,                 & !< logical:
    prescribedFire        => c_switch%prescribedFire,         & !< logical:
    lnduseon              => c_switch%lnduseon,               & !< logical:
    PFTCompetition        => c_switch%PFTCompetition,         & !< logical:
    doperpftoutput        => c_switch%doperpftoutput,         & !< logical:
    dopertileoutput       => c_switch%dopertileoutput,        & !< logical:
    transientOBSWETF      => c_switch%transientOBSWETF,       & !< logical:
    fixedYearOBSWETF      => c_switch%fixedYearOBSWETF,       & !< integer:
    Ncycle_on             => c_switch%Ncycle_on,              & !< logical:

    FAREROT => class_rot%FAREROT,                             & !< real, dimension(:,:) : Fractional coverage of mosaic tile on modelled area

    fcancmxrow        => vrot%fcancmx,                        & !< real, dimension(:,:,:) :
    gppvegrow         => vrot%gppveg,                         & !< real, dimension(:,:,:) :
    leafns2srow       => vrot%leafns2s,                       & !< real, dimension(:,:,:) :
    stemns2srow       => vrot%stemns2s,                       & !< real, dimension(:,:,:) :
    rootns2srow       => vrot%rootns2s,                       & !< real, dimension(:,:,:) :
    re_alloc_s2lrow   => vrot%re_alloc_s2l,                   & !< real, dimension(:,:,:) :
    re_alloc_r2lrow   => vrot%re_alloc_r2l,                   & !< real, dimension(:,:,:) :
    re_alloc_sr2lrow  => vrot%re_alloc_sr2l,                  & !< real, dimension(:,:,:) :
    nepvegrow         => vrot%nepveg,                         & !< real, dimension(:,:,:) :
    nbpvegrow         => vrot%nbpveg,                         & !< real, dimension(:,:,:) :
    nppvegrow         => vrot%nppveg,                         & !< real, dimension(:,:,:) :
    hetroresvegrow    => vrot%hetroresveg,                    & !< real, dimension(:,:,:) :
    autoresvegrow     => vrot%autoresveg,                     & !< real, dimension(:,:,:) :
    litresvegrow      => vrot%litresveg,                      & !< real, dimension(:,:,:,:) :
    litrfallvegrow    => vrot%litrfallveg,                    & !< real, dimension(:,:,:,:) :
    soilcresvegrow    => vrot%soilcresveg,                    & !< real, dimension(:,:,:,:) :
    rmlvegaccrow      => vrot%rmlvegacc,                      & !< real, dimension(:,:,:) :
    rmsvegrow         => vrot%rmsveg,                         & !< real, dimension(:,:,:) :
    rmrvegrow         => vrot%rmrveg,                         & !< real, dimension(:,:,:) :
    rgvegrow          => vrot%rgveg,                          & !< real, dimension(:,:,:) :
    ailcgrow          => vrot%ailcg,                          & !< real, dimension(:,:,:) :
    emit_co2row       => vrot%emit_co2,                       & !< real, dimension(:,:,:) :
    emit_corow        => vrot%emit_co,                        & !< real, dimension(:,:,:) :
    emit_ch4row       => vrot%emit_ch4,                       & !< real, dimension(:,:,:) :
    emit_nmhcrow      => vrot%emit_nmhc,                      & !< real, dimension(:,:,:) :
    emit_h2row        => vrot%emit_h2,                        & !< real, dimension(:,:,:) :
    emit_pm25row      => vrot%emit_pm25,                      & !< real, dimension(:,:,:) :
    emit_tpmrow       => vrot%emit_tpm,                       & !< real, dimension(:,:,:) :
    emit_tcrow        => vrot%emit_tc,                        & !< real, dimension(:,:,:) :
    emit_ocrow        => vrot%emit_oc,                        & !< real, dimension(:,:,:) :
    emit_bcrow        => vrot%emit_bc,                        & !< real, dimension(:,:,:) :
    burnfracrow       => vrot%burnfrac,                       & !< real, dimension(:,:) :
    smfuncvegrow      => vrot%smfuncveg,                      & !< real, dimension(:,:,:) :
    btermrow          => vrot%bterm,                          & !< real, dimension(:,:,:) :
    ltermrow          => vrot%lterm,                          & !< real, dimension(:,:) :
    mtermrow          => vrot%mterm,                          & !< real, dimension(:,:,:) :
    lucemcomrow       => vrot%lucemcom,                       & !< real, dimension(:,:) :
    lucltrinrow       => vrot%lucltrin,                       & !< real, dimension(:,:) :
    lucsocinrow       => vrot%lucsocin,                       & !< real, dimension(:,:) :
    lucemcomnrow       => vrot%lucemcomn,                       & !< real, dimension(:,:) :
    lucltrinnrow       => vrot%lucltrinn,                       & !< real, dimension(:,:) :
    lucsocinnrow       => vrot%lucsocinn,                       & !< real, dimension(:,:) :
    ch4WetSpecrow        => vrot%ch4WetSpec,                  & !< real, dimension(:,:) :
    wetfdynrow        => vrot%wetfdyn,                        & !< real, dimension(:,:) :
    ch4WetDynrow        => vrot%ch4WetDyn,                    & !< real, dimension(:,:) :
    ch4soillsrow      => vrot%ch4_soills,                     & !< real, dimension(:,:) :
    litrmassrow       => vrot%litrmass,                       & !< real, dimension(:,:,:,:) :
    soilcmasrow       => vrot%soilcmas,                       & !< real, dimension(:,:,:,:) :
    nh4_massrow       => vrot%nh4_mass,                       & !< real, dimension(:,:,:) :
    no3_massrow       => vrot%no3_mass,                       & !< real, dimension(:,:,:) :
    nlitrmassrow      => vrot%nlitrmass,                      & !< real, dimension(:,:,:) :
    soilnmasrow       => vrot%soilnmas,                       & !< real, dimension(:,:,:) :
    ngleafmasrow      => vrot%ngleafmas,                      & !< real, dimension(:,:,:) :
    ngleafmas_NSrow    => vrot%ngleafmas_ns,                  & !< real, dimension(:,:,:) :
    ngleafmassrow     => vrot%ngleafmas_s,                    & !< real, dimension(:,:,:) :
    nbleafmasrow      => vrot%nbleafmas,                      & !< real, dimension(:,:,:) :
    nstemmassrow      => vrot%nstemmass,                      & !< real, dimension(:,:,:) :
    nstemmass_NSrow    => vrot%nstemmass_ns,                  & !< real, dimension(:,:,:) :
    nstemmasssrow     => vrot%nstemmass_s,                    & !< real, dimension(:,:,:) :
    nrootmassrow      => vrot%nrootmass,                      & !< real, dimension(:,:,:) :
    nrootmass_NSrow    => vrot%nrootmass_ns,                  & !< real, dimension(:,:,:) :
    nrootmasssrow     => vrot%nrootmass_s,                    & !< real, dimension(:,:,:) :
    nvgbiomasvegrow  => vrot%nvgbiomas_veg,                  & !< real, dimension(:,:,:) :
    c2nveg_lrow       => vrot%c2nveg_l,                       & !< real, dimension(:,:,:) :
    c2nveg_srow       => vrot%c2nveg_s,                       & !< real, dimension(:,:,:) :
    c2nveg_rrow       => vrot%c2nveg_r,                       & !< real, dimension(:,:,:) :
    c2nveg_wprow      => vrot%c2nveg_wp,                      & !< real, dimension(:,:,:) :
    c2nveg_litrrow    => vrot%c2nveg_litr,                    & !< real, dimension(:,:,:) :
    c2nveg_humusrow   => vrot%c2nveg_humus,                   & !< real, dimension(:,:,:) :
    vgbiomas_vegrow   => vrot%vgbiomas_veg,                   & !< real, dimension(:,:,:) :
    stemmassrow       => vrot%stemmass,                       & !< real, dimension(:,:,:) :
    stemmass_NSrow     => vrot%stemmass_ns,                   & !< real, dimension(:,:,:) :
    stemmasssrow      => vrot%stemmass_s,                     & !< real, dimension(:,:,:) :
    rootmassrow       => vrot%rootmass,                       & !< real, dimension(:,:,:) :
    rootmass_NSrow     => vrot%rootmass_ns,                   & !< real, dimension(:,:,:) :
    rootmasssrow      => vrot%rootmass_s,                     & !< real, dimension(:,:,:) :
    dstcemls3row      => vrot%dstcemls3,                      & !< real, dimension(:,:) :
    lfstatusrow       => vrot%lfstatus,                       & !< integer, dimension(:,:,:) :

    npprow            => vrot%npp,                            & !< real, dimension(:,:) :
    neprow            => vrot%nep,                            & !< real, dimension(:,:) :
    nepCMIProw        => vrot%nepCMIP,                        & !< real, dimension(:,:) :
    nbprow            => vrot%nbp,                            & !< real, dimension(:,:) :
    gpprow            => vrot%gpp,                            & !< real, dimension(:,:) :
    hetroresrow       => vrot%hetrores,                       & !< real, dimension(:,:) :
    autoresrow        => vrot%autores,                        & !< real, dimension(:,:) :
    soilcresprow      => vrot%soilcresp,                      & !< real, dimension(:,:) :
    rgrow             => vrot%rg,                             & !< real, dimension(:,:) :
    litresrow         => vrot%litres,                         & !< real, dimension(:,:) :
    socresrow         => vrot%socres,                         & !< real, dimension(:,:) :
    vgbiomasrow       => vrot%vgbiomas,                       & !< real, dimension(:,:) :
    gavgltmsrow       => vrot%gavgltms,                       & !< real, dimension(:,:) :
    gavgscmsrow       => vrot%gavgscms,                       & !< real, dimension(:,:) :

    nppmossrow         => vrot%nppmoss,                       & !< real, dimension(:,:) :
    armossrow          => vrot%armoss,                        & !< real, dimension(:,:) :

    bmasvegrow        => vrot%bmasveg,                        & !< real, dimension(:,:,:) :
    cmasvegcrow       => vrot%cmasvegc,                       & !< real, dimension(:,:,:) :
    veghghtrow        => vrot%veghght,                        & !< real, dimension(:,:,:) :
    rootdpthrow       => vrot%rootdpth,                       & !< real, dimension(:,:,:) :
    rmlrow            => vrot%rml,                            & !< real, dimension(:,:) :
    rmsrow            => vrot%rms,                            & !< real, dimension(:,:) :
    tltrleafrow       => vrot%tltrleaf,                       & !< real, dimension(:,:,:) :
    tltrstemrow       => vrot%tltrstem,                       & !< real, dimension(:,:,:) :
    tltrrootrow       => vrot%tltrroot,                       & !< real, dimension(:,:,:) :
    leaflitrrow       => vrot%leaflitr,                       & !< real, dimension(:,:,:) :
    roottemprow       => vrot%roottemp,                       & !< real, dimension(:,:,:) :
    afrleafrow        => vrot%afrleaf,                        & !< real, dimension(:,:,:) :
    afrstemrow        => vrot%afrstem,                        & !< real, dimension(:,:,:) :
    afrrootrow        => vrot%afrroot,                        & !< real, dimension(:,:,:) :
    wtstatusrow       => vrot%wtstatus,                       & !< real, dimension(:,:,:) :
    ltstatusrow       => vrot%ltstatus,                       & !< real, dimension(:,:,:) :
    rmrrow            => vrot%rmr,                            & !< real, dimension(:,:) :
    gleafmasrow       => vrot%gleafmas,                       & !< real, dimension(:,:,:) :
    gleafmas_NSrow     => vrot%gleafmas_ns,                   & !< real, dimension(:,:,:) :
    gleafmassrow      => vrot%gleafmas_s,                     & !< real, dimension(:,:,:) :
    bleafmasrow       => vrot%bleafmas,                       & !< real, dimension(:,:,:) :
    gavglairow        => vrot%gavglai,                        & !< real, dimension(:,:) :
    slairow           => vrot%slai,                           & !< real, dimension(:,:,:) :
    ailcbrow          => vrot%ailcb,                          & !< real, dimension(:,:,:) :
    flhrlossrow       => vrot%flhrloss,                       & !< real, dimension(:,:,:) :
    rmatctemrow       => vrot%rmatctem,                       & !< real, dimension(:,:,:,:) :
    dstcemlsrow       => vrot%dstcemls,                       & !< real, dimension(:,:) :
    litrfallrow       => vrot%litrfall,                       & !< real, dimension(:,:) :
    humiftrsrow       => vrot%humiftrs,                       & !< real, dimension(:,:) :

    bnftotrow        => vrot%bnf_tot,                      & !< real, dimension(:,:,:)   :
    bnffreerow       => vrot%bnf_free,                     & !< real, dimension(:,:,:)   :
    bnfantrow        => vrot%bnf_ant,                      & !< real, dimension(:,:,:)   :
    bnfnatrow        => vrot%bnf_nat,                      & !< real, dimension(:,:,:)   :
    nitrifvegrow      => vrot%nitrifveg,                      & !< real, dimension(:,:,:) :
    no_nitvegrow      => vrot%no_nitveg,                      & !< real, dimension(:,:,:) :
    no_denitvegrow    => vrot%no_denitveg,                    & !< real, dimension(:,:,:) :
    no_nitdenitvegrow => vrot%no_nitdenitveg,                 & !< real, dimension(:,:,:) :
    n2o_nitvegrow     => vrot%n2o_nitveg,                     & !< real, dimension(:,:,:) :
    n2o_denitvegrow   => vrot%n2o_denitveg,                   & !< real, dimension(:,:,:) :
    n2o_nitdenitvegrow=> vrot%n2o_nitdenitveg,                & !< real, dimension(:,:,:) :
    n2_denitvegrow    => vrot%n2_denitveg,                    & !< real, dimension(:,:,:) :
    nvolvegrow        => vrot%nvolveg,                        & !< real, dimension(:,:,:) :
    nleachvegrow      => vrot%nleachveg,                      & !< real, dimension(:,:,:) :
    appl_fertrow      => vrot%appl_fert,                      & !< real, dimension(:,:,:) :
    ndep_nh4row       => vrot%ndep_nh4,                       & !< real, dimension(:,:,:) :
    ndep_no3row       => vrot%ndep_no3,                       & !< real, dimension(:,:,:) :
    ndemandveg_wp_npprow     => vrot%ndemandveg_wp_npp,       & !< real, dimension(:,:,:) :
    nuptakeveg_p_nh4row      => vrot%nuptakeveg_p_nh4,        & !< real, dimension(:,:,:) :
    nuptakeveg_p_no3row      => vrot%nuptakeveg_p_no3,        & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_nh4row => vrot%nuptakeveg_a_actl_nh4,   & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_no3row => vrot%nuptakeveg_a_actl_no3,   & !< real, dimension(:,:,:) :
    nallocveg_lrow       => vrot%nallocveg_l,                 & !< real, dimension(:,:,:) :
    nallocveg_srow       => vrot%nallocveg_s,                 & !< real, dimension(:,:,:) :
    nallocveg_rrow       => vrot%nallocveg_r,                 & !< real, dimension(:,:,:) :
    nresorpedveg_srow    => vrot%nresorpedveg_s,              & !< real, dimension(:,:,:) :
    nresorpedveg_rrow    => vrot%nresorpedveg_r,              & !< real, dimension(:,:,:) :
    nre_allocveg_s2lrow  => vrot%nre_allocveg_s2l,            & !< real, dimension(:,:,:) :
    nre_allocveg_r2lrow  => vrot%nre_allocveg_r2l,            & !< real, dimension(:,:,:) :
    nleafns2svegrow      => vrot%nleafns2sveg,                & !< real, dimension(:,:,:) :
    nstemns2svegrow      => vrot%nstemns2sveg,                & !< real, dimension(:,:,:) :
    nrootns2svegrow      => vrot%nrootns2sveg,                & !< real, dimension(:,:,:) :
    nlitrveg_lrow        => vrot%nlitrveg_l,                  & !< real, dimension(:,:,:) :
    nlitrveg_srow        => vrot%nlitrveg_s,                  & !< real, dimension(:,:,:) :
    nlitrveg_rrow        => vrot%nlitrveg_r,                  & !< real, dimension(:,:,:) :
    gl2bl_grass_nfluxrow => vrot%gl2bl_grass_nflux,           & !< real, dimension(:,:,:) :
    nhumtrsvegrow        => vrot%nhumtrsveg,                  & !< real, dimension(:,:,:) :
    nmineralveg_litrrow  => vrot%nmineralveg_litr,            & !< real, dimension(:,:,:) :
    nmineralveg_humusrow => vrot%nmineralveg_humus,           & !< real, dimension(:,:,:) :
    nimmobilveg_nh4row   => vrot%nimmobilveg_nh4,             & !< real, dimension(:,:,:) :
    nimmobilveg_no3row   => vrot%nimmobilveg_no3,             & !< real, dimension(:,:,:) :
    fNnetlandvegrow      => vrot%fNnetlandveg                 & !< real, dimension(:,:,:) :
    )
    
    ! First set the local variables to 0.
    gpp_g(:) = 0.0 ; npp_g(:) = 0.0 ; nep_g(:) = 0.0 ; nepCMIP_g(:) = 0.0 ; nbp_g(:) = 0.0 ; autores_g(:) = 0.0
    hetrores_g(:) = 0.0 ; litres_g(:) = 0.0 ; socres_g(:) = 0.0 ; dstcemls_g(:) = 0.0
    dstcemls3_g(:) = 0.0 ; litrfall_g(:) = 0.0 ; humiftrs_g(:) = 0.0 ; rml_g(:) = 0.0
    rms_g(:) = 0.0 ; rmr_g(:) = 0.0 ; rg_g(:) = 0.0 ; vgbiomas_g(:) = 0.0 ; totcmass_g(:) = 0.0
    gavglai_g(:) = 0.0 ; gavgltms_g(:) = 0.0 ; gavgscms_g(:) = 0.0 ; ailcg_g(:) = 0.0
    ailcb_g(:) = 0.0 ; tcanoacc_out_g(:) = 0.0 ; burnfrac_g(:) = 0.0 ; smfuncveg_g(:) = 0.0
    lucemcom_g(:) = 0.0 ; lucltrin_g(:) = 0.0 ; lucsocin_g(:) = 0.0 ; 
    lucemcomn_g(:) = 0.0 ; lucltrinn_g(:) = 0.0 ; lucsocinn_g(:) = 0.0 ; emit_co2_g(:) = 0.0
    leaflitr_g(:) = 0.0 ; tltrleaf_g(:) = 0.0 ; tltrstem_g(:) = 0.0 ; tltrroot_g(:) = 0.0
    gleafmas_g(:) = 0.0 ; bleafmas_g(:) = 0.0 ; stemmass_g(:) = 0.0 ; rootmass_g(:) = 0.0
    litrmass_g(:) = 0.0 ; soilcmas_g(:) = 0.0 ; veghght_g(:) = 0.0 ; rootdpth_g(:) = 0.0
    roottemp_g(:) = 0.0 ; slai_g(:) = 0.0 ; ch4WetSpec_G(:) = 0.0
    WETFDYN_G(:) = 0.0 ; ch4WetDyn_G(:) = 0.0 ; ch4soills_g(:) = 0.0
    rmatctem_g(:,:) = 0.0
    afrleaf_g(:,:) = 0.0 ; afrstem_g(:,:) = 0.0 ; afrroot_g(:,:) = 0.0
    leaflitr_t(:,:) = 0.0 ; tltrleaf_t(:,:) = 0.0 ; tltrstem_t(:,:) = 0.0
    tltrroot_t(:,:) = 0.0 ; ailcg_t(:,:) = 0.0 ; ailcb_t(:,:) = 0.0
    afrleaf_t(:,:) = 0.0 ; afrstem_t(:,:) = 0.0 ; afrroot_t(:,:) = 0.0
    veghght_t(:,:) = 0.0 ; rootdpth_t(:,:) = 0.0 ; roottemp_t(:,:) = 0.0
    slai_t(:,:) = 0.0 ; gleafmas_t(:,:) = 0.0 ; bleafmas_t(:,:) = 0.0
    stemmass_t(:,:) = 0.0 ; rootmass_t(:,:) = 0.0
    emit_co2_t(:,:) = 0.0; rmatctem_t(:,:,:) = 0.0
    soilcmas_t(:,:,:) = 0.0 ; litrmass_t(:,:,:) = 0.0
    gleafmas_ns_g(:) = 0.0 ; gleafmas_s_g(:) = 0.0 ; stemmass_ns_g(:) = 0.0
    stemmass_s_g(:) = 0.0 ; rootmass_ns_g(:) = 0.0 ; rootmass_s_g(:) = 0.0
    leafns2s_g(:) = 0.0 ; stemns2s_g(:) = 0.0 ; rootns2s_g(:) = 0.0
    re_alloc_s2l_g(:) = 0.0 ; re_alloc_r2l_g(:) = 0.0 ; re_alloc_sr2l_g(:) = 0.0
    bnf_tot_g(:) = 0.0
    bnf_free_g(:) = 0.0 ; bnf_ant_g(:) = 0.0 ; bnf_nat_g(:) = 0.0
    nh4_mass_g(:) = 0.0 ; no3_mass_g(:) = 0.0 ; nitrif_g(:) = 0.0 ; no_nit_g(:) = 0.0 ; no_denit_g(:) = 0.0
    no_nitdenit_g(:) = 0.0 ; n2o_nit_g(:) = 0.0 ; n2o_denit_g(:) = 0.0
    n2o_nitdenit_g(:) = 0.0 ; n2_denit_g(:) = 0.0 ; nvol_g(:) = 0.0 ; nleach_g(:) = 0.0
    appl_fert_g(:) = 0.0 ; ndep_nh4_g(:) = 0.0 ; ndep_no3_g(:) = 0.0
    ngleafmas_g(:) = 0.0 ; ngleafmas_ns_g(:) = 0.0 ; ngleafmas_s_g(:)  = 0.0
    nbleafmas_g(:) = 0.0 ; nstemmass_g(:) = 0.0 ; nstemmass_ns_g(:) = 0.0
    nstemmass_s_g(:) = 0.0 ; nrootmass_g(:) = 0.0 ; nrootmass_ns_g(:) = 0.0
    nrootmass_s_g(:) = 0.0 ; nlitrmass_g(:) = 0.0 ; soilnmas_g(:) = 0.0
    ndemand_wp_npp_g(:) = 0.0 ;nuptake_p_nh4_g(:) = 0.0
    nuptake_p_no3_g(:) = 0.0 ; nuptake_a_actl_nh4_g(:) = 0.0
    nuptake_a_actl_no3_g(:) = 0.0
    nalloc_l_g(:) = 0.0 ; nalloc_s_g(:) = 0.0 ; nalloc_r_g(:) = 0.0
    nresorped_s_g(:) = 0.0 ; nresorped_r_g(:) = 0.0 ; nre_alloc_s2l_g(:) = 0.0
    nre_alloc_r2l_g(:) = 0.0; nleafns2s_g(:) = 0.0 ; nstemns2s_g(:) = 0.0
    nrootns2s_g(:) = 0.0 ; nlitr_l_g(:) = 0.0 ; nlitr_s_g(:) = 0.0
    nlitr_r_g(:) = 0.0 ; gl2bl_grass_nflux_g(:) = 0.0
    c2n_l_g(:) = 0.0 ; c2n_s_g(:) = 0.0 ; c2n_r_g(:) = 0.0
    c2n_wp_g(:) = 0.0 ; c2n_litr_g(:) = 0.0 ; c2n_humus_g(:) = 0.0
    nhumtrs_g(:) = 0.0 ; nmineral_litr_g(:) = 0.0 ; nmineral_humus_g(:) = 0.0
    nimmobil_nh4_g(:) = 0.0 ; nimmobil_no3_g(:) = 0.0 ; nvgbiomas_g(:) = 0.0
    fNnetland_g(:) = 0.0
    nh4_mass_t(:,:) = 0.0 ; no3_mass_t(:,:) = 0.0
    ngleafmas_t(:,:) = 0.0 ; ngleafmas_ns_t(:,:) = 0.0 ; ngleafmas_s_t(:,:) = 0.0
    nbleafmas_t(:,:) = 0.0 ; nstemmass_t(:,:) = 0.0 ; nstemmass_ns_t(:,:) = 0.0
    nstemmass_s_t(:,:) = 0.0 ; nrootmass_t(:,:) = 0.0 ; nrootmass_ns_t(:,:) = 0.0
    nrootmass_s_t(:,:) = 0.0 ; nlitrmass_t(:,:) = 0.0 ; soilnmas_t(:,:) = 0.0
    nvgbiomas_t(:,:) = 0.0
    c2n_l_t(:,:) = 0.0 ; c2n_s_t(:,:) = 0.0 ; c2n_r_t(:,:) = 0.0
    c2n_wp_t(:,:) = 0.0 ; c2n_litr_t(:,:) = 0.0 ; c2n_humus_t(:,:) = 0.0
    bnf_tot_t(:,:) = 0.0
    bnf_free_t(:,:) = 0.0 ; bnf_ant_t(:,:) = 0.0 ; bnf_nat_t(:,:) = 0.0
    nitrif_t(:,:) = 0.0 ; no_nit_t(:,:) = 0.0 ; no_denit_t(:,:) = 0.0
    no_nitdenit_t(:,:) = 0.0 ; n2o_nit_t(:,:) = 0.0 ; n2o_denit_t(:,:) = 0.0
    n2o_nitdenit_t(:,:) = 0.0 ; n2_denit_t(:,:) = 0.0 ; nvol_t(:,:) = 0.0 ; nleach_t(:,:) = 0.0
    appl_fert_t(:,:) = 0.0 ; ndep_nh4_t(:,:) = 0.0 ; ndep_no3_t(:,:) = 0.0
    ndemand_wp_npp_t(:,:) = 0.0 ;nuptake_p_nh4_t(:,:) = 0.0
    nuptake_p_no3_t(:,:) = 0.0 ; nuptake_a_actl_nh4_t(:,:) = 0.0
    nuptake_a_actl_no3_t(:,:) = 0.0
    nalloc_l_t(:,:) = 0.0 ; nalloc_s_t(:,:) = 0.0 ; nalloc_r_t(:,:) = 0.0
    nresorped_s_t(:,:) = 0.0 ; nresorped_r_t(:,:) = 0.0 ; nre_alloc_s2l_t(:,:) = 0.0
    nre_alloc_r2l_t(:,:) = 0.0; nleafns2s_t(:,:) = 0.0 ; nstemns2s_t(:,:) = 0.0
    nrootns2s_t(:,:) = 0.0 ; nlitr_l_t(:,:) = 0.0 ; nlitr_s_t(:,:) = 0.0
    nlitr_r_t(:,:) = 0.0 ; gl2bl_grass_nflux_t(:,:) = 0.0
    nhumtrs_t(:,:) = 0.0 ; nmineral_litr_t(:,:) = 0.0 ; nmineral_humus_t(:,:) = 0.0
    nimmobil_nh4_t(:,:) = 0.0 ; nimmobil_no3_t(:,:) = 0.0
    fNnetland_t(:,:) = 0.0
    leafns2s_t(:,:) = 0.0 ; stemns2s_t(:,:) = 0.0 
    rootns2s_t(:,:) = 0.0 ; re_alloc_s2l_t (:,:) = 0.0 ; vgbiomas_t(:,:) = 0.0
    re_alloc_r2l_t(:,:) = 0.0 ; re_alloc_sr2l_t(:,:) = 0.0; emit_ch4_t(:,:) = 0.0
    gleafmas_ns_t(:,:) = 0.0 ; gleafmas_s_t(:,:) = 0.0 ; emit_co2_t(:,:) = 0.0
    smfuncveg_t(:,:) = 0.0
    stemmass_ns_t(:,:) = 0.0 ; stemmass_s_t(:,:) = 0.0 ; emit_ch4_g(:) = 0.0
    rootmass_ns_t(:,:) = 0.0 ; rootmass_s_t(:,:) = 0.0 ; anvegrow_g = 0.0 ; laimaxg_t = 0.0 ;

    !> Aggregate to the tile avg vars:
    do i = 1,nltest
      do m = 1,nmtest
        barefrac = 1.0
        do j = 1,icc
          barefrac = barefrac - fcancmxrow(i,m,j)
          leaflitr_t(i,m) = leaflitr_t(i,m) + leaflitrrow(i,m,j) * fcancmxrow(i,m,j)
          tltrleaf_t(i,m) = tltrleaf_t(i,m) + tltrleafrow(i,m,j) * fcancmxrow(i,m,j)
          tltrstem_t(i,m) = tltrstem_t(i,m) + tltrstemrow(i,m,j) * fcancmxrow(i,m,j)
          tltrroot_t(i,m) = tltrroot_t(i,m) + tltrrootrow(i,m,j) * fcancmxrow(i,m,j)
          veghght_t(i,m) = veghght_t(i,m) + veghghtrow(i,m,j) * fcancmxrow(i,m,j)
          rootdpth_t(i,m) = rootdpth_t(i,m) + rootdpthrow(i,m,j) * fcancmxrow(i,m,j)
          roottemp_t(i,m) = roottemp_t(i,m) + roottemprow(i,m,j) * fcancmxrow(i,m,j)
          slai_t(i,m) = slai_t(i,m) + slairow(i,m,j) * fcancmxrow(i,m,j)
          afrleaf_t(i,m) = afrleaf_t(i,m) + afrleafrow(i,m,j) * fcancmxrow(i,m,j)
          afrstem_t(i,m) = afrstem_t(i,m) + afrstemrow(i,m,j) * fcancmxrow(i,m,j)
          afrroot_t(i,m) = afrroot_t(i,m) + afrrootrow(i,m,j) * fcancmxrow(i,m,j)
          ailcg_t(i,m) = ailcg_t(i,m) + ailcgrow(i,m,j) * fcancmxrow(i,m,j)
          ailcb_t(i,m) = ailcb_t(i,m) + ailcbrow(i,m,j) * fcancmxrow(i,m,j)
          leafns2s_t(i,m) = leafns2s_t(i,m) + leafns2srow(i,m,j) * fcancmxrow(i,m,j)
          stemns2s_t(i,m) = stemns2s_t(i,m) + stemns2srow(i,m,j) * fcancmxrow(i,m,j)
          rootns2s_t(i,m) = rootns2s_t(i,m) + rootns2srow(i,m,j) * fcancmxrow(i,m,j)
          re_alloc_s2l_t(i,m) = re_alloc_s2l_t(i,m) + re_alloc_s2lrow(i,m,j) * fcancmxrow(i,m,j)
          re_alloc_r2l_t(i,m) = re_alloc_r2l_t(i,m) + re_alloc_r2lrow(i,m,j) * fcancmxrow(i,m,j)
          re_alloc_sr2l_t(i,m) = re_alloc_sr2l_t(i,m) + re_alloc_sr2lrow(i,m,j) * fcancmxrow(i,m,j)
          gleafmas_t(i,m) = gleafmas_t(i,m) + gleafmasrow(i,m,j) * fcancmxrow(i,m,j)
          gleafmas_ns_t(i,m) = gleafmas_ns_t(i,m) + gleafmas_NSrow(i,m,j) * fcancmxrow(i,m,j)
          gleafmas_s_t(i,m) = gleafmas_s_t(i,m) + gleafmassrow(i,m,j) * fcancmxrow(i,m,j)
          vgbiomas_t(i,m) = vgbiomas_t(i,m) + vgbiomas_vegrow(i,m,j) * fcancmxrow(i,m,j)
          bleafmas_t(i,m) = bleafmas_t(i,m) + bleafmasrow(i,m,j) * fcancmxrow(i,m,j)
          stemmass_t(i,m) = stemmass_t(i,m) + stemmassrow(i,m,j) * fcancmxrow(i,m,j)
          stemmass_ns_t(i,m) = stemmass_ns_t(i,m) + stemmass_NSrow(i,m,j) * fcancmxrow(i,m,j)
          stemmass_s_t(i,m) = stemmass_s_t(i,m) + stemmasssrow(i,m,j) * fcancmxrow(i,m,j)
          rootmass_t(i,m) = rootmass_t(i,m) + rootmassrow(i,m,j) * fcancmxrow(i,m,j)
          rootmass_ns_t(i,m) = rootmass_ns_t(i,m) + rootmass_NSrow(i,m,j) * fcancmxrow(i,m,j)
          rootmass_s_t(i,m) = rootmass_s_t(i,m) + rootmasssrow(i,m,j) * fcancmxrow(i,m,j)
          emit_co2_t(i,m) = emit_co2_t(i,m) + emit_co2row(i,m,j) * fcancmxrow(i,m,j)
          emit_ch4_t(i,m) = emit_ch4_t(i,m) + emit_ch4row(i,m,j) * fcancmxrow(i,m,j)
          smfuncveg_t(i,m) = smfuncveg_t(i,m) + smfuncvegrow(i,m,j) * fcancmxrow(i,m,j)

          do k = 1,ignd
            rmatctem_t(i,m,k) = rmatctem_t(i,m,k) + rmatctemrow(i,m,j,k) * fcancmxrow(i,m,j)
            litrmass_t(i,m,k) = litrmass_t(i,m,k) + litrmassrow(i,m,j,k)*fcancmxrow(i,m,j)
            soilcmas_t(i,m,k) = soilcmas_t(i,m,k) + soilcmasrow(i,m,j,k)*fcancmxrow(i,m,j)
          end do

          if (Ncycle_on) then
             nh4_mass_t(i,m) = nh4_mass_t(i,m) + nh4_massrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             no3_mass_t(i,m) = no3_mass_t(i,m) + no3_massrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             appl_fert_t(i,m) = appl_fert_t(i,m) + appl_fertrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             ndep_nh4_t(i,m) = ndep_nh4_t(i,m) + ndep_nh4row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             ndep_no3_t(i,m) = ndep_no3_t(i,m) + ndep_no3row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             bnf_tot_t(i,m) =  bnf_tot_t(i,m) + bnftotrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             bnf_free_t(i,m) = bnf_free_t(i,m) + bnffreerow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             bnf_ant_t(i,m) = bnf_ant_t(i,m) + bnfantrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             bnf_nat_t(i,m) = bnf_nat_t(i,m) + bnfnatrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nitrif_t(i,m) = nitrif_t(i,m) + nitrifvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             no_nit_t(i,m) = no_nit_t(i,m) + no_nitvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             no_denit_t(i,m) = no_denit_t(i,m) + no_denitvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             no_nitdenit_t(i,m) = no_nitdenit_t(i,m) + no_nitdenitvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             n2o_nit_t(i,m) = n2o_nit_t(i,m) + n2o_nitvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             n2o_denit_t(i,m) = n2o_denit_t(i,m) + n2o_denitvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             n2o_nitdenit_t(i,m) = n2o_nitdenit_t(i,m) + n2o_nitdenitvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             n2_denit_t(i,m) = n2_denit_t(i,m) + n2_denitvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nvol_t(i,m) = nvol_t(i,m) + nvolvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nleach_t(i,m) = nleach_t(i,m) + nleachvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nuptake_p_nh4_t(i,m) = nuptake_p_nh4_t(i,m) + nuptakeveg_p_nh4row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nuptake_p_no3_t(i,m) = nuptake_p_no3_t(i,m) + nuptakeveg_p_no3row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nuptake_a_actl_nh4_t(i,m) = nuptake_a_actl_nh4_t(i,m) + nuptakeveg_a_actl_nh4row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nuptake_a_actl_no3_t(i,m) = nuptake_a_actl_no3_t(i,m) + nuptakeveg_a_actl_no3row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nalloc_l_t(i,m) = nalloc_l_t(i,m) + nallocveg_lrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nalloc_s_t(i,m) = nalloc_s_t(i,m) + nallocveg_srow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nalloc_r_t(i,m) = nalloc_r_t(i,m) + nallocveg_rrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nresorped_s_t(i,m) = nresorped_s_t(i,m) + nresorpedveg_srow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nresorped_r_t(i,m) = nresorped_r_t(i,m) + nresorpedveg_rrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nre_alloc_s2l_t(i,m) = nre_alloc_s2l_t(i,m) + nre_allocveg_s2lrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nre_alloc_r2l_t(i,m) = nre_alloc_r2l_t(i,m) + nre_allocveg_r2lrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nleafns2s_t(i,m) = nleafns2s_t(i,m) + nleafns2svegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nstemns2s_t(i,m) = nstemns2s_t(i,m) + nstemns2svegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nrootns2s_t(i,m) = nrootns2s_t(i,m) + nrootns2svegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nlitr_l_t(i,m) = nlitr_l_t(i,m) + nlitrveg_lrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nlitr_s_t(i,m) = nlitr_s_t(i,m) + nlitrveg_srow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nlitr_r_t(i,m) = nlitr_r_t(i,m) + nlitrveg_rrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             gl2bl_grass_nflux_t(i,m) = gl2bl_grass_nflux_t(i,m) + gl2bl_grass_nfluxrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nhumtrs_t(i,m) = nhumtrs_t(i,m) + nhumtrsvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nmineral_litr_t(i,m) = nmineral_litr_t(i,m) + nmineralveg_humusrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nmineral_humus_t(i,m) = nmineral_humus_t(i,m) + nmineralveg_humusrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nimmobil_nh4_t(i,m) = nimmobil_nh4_t(i,m) + nimmobilveg_nh4row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             nimmobil_no3_t(i,m) = nimmobil_no3_t(i,m) + nimmobilveg_no3row(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             fNnetland_t(i,m) = fNnetland_t(i,m) + fNnetlandvegrow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
             ngleafmas_t(i,m) = ngleafmas_t(i,m) + ngleafmasrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             ngleafmas_ns_t(i,m) = ngleafmas_ns_t(i,m) + ngleafmas_NSrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             ngleafmas_s_t(i,m) = ngleafmas_s_t(i,m) + ngleafmassrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nbleafmas_t(i,m) = nbleafmas_t(i,m) + nbleafmasrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nstemmass_t(i,m) = nstemmass_t(i,m) + nstemmassrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nstemmass_ns_t(i,m) = nstemmass_ns_t(i,m) + nstemmass_NSrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nstemmass_s_t(i,m) = nstemmass_s_t(i,m) + nstemmasssrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nrootmass_t(i,m) = nrootmass_t(i,m) + nrootmassrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nrootmass_ns_t(i,m) = nrootmass_ns_t(i,m) + nrootmass_NSrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nrootmass_s_t(i,m) = nrootmass_s_t(i,m) + nrootmasssrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nlitrmass_t(i,m) = nlitrmass_t(i,m) + nlitrmassrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             soilnmas_t(i,m) = soilnmas_t(i,m) + soilnmasrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             nvgbiomas_t(i,m) = nvgbiomas_t(i,m) + nvgbiomasvegrow(i,m,j) * fcancmxrow(i,m,j) * convertg2kg
             ndemand_wp_npp_t(i,m) = ndemand_wp_npp_t(i,m) + ndemandveg_wp_npprow(i,m,j) * fcancmxrow(i,m,j) * convertkgN
          end if !for Ncycle_on

        end do ! icc

        !> Do the bare ground also:
        do k = 1,ignd
          litrmass_t(i,m,k) = litrmass_t(i,m,k) + litrmassrow(i,m,iccp1,k)*barefrac
          soilcmas_t(i,m,k) = soilcmas_t(i,m,k) + soilcmasrow(i,m,iccp1,k)*barefrac
        end do
        if (Ncycle_on) then
          nh4_mass_t(i,m) = nh4_mass_t(i,m) + nh4_massrow(i,m,iccp1) * barefrac * convertg2kg
          no3_mass_t(i,m) = no3_mass_t(i,m) + no3_massrow(i,m,iccp1) * barefrac * convertg2kg
          nlitrmass_t(i,m) = nlitrmass_t(i,m) + nlitrmassrow(i,m,iccp1) * barefrac * convertg2kg
          soilnmas_t(i,m) = soilnmas_t(i,m) + soilnmasrow(i,m,iccp1) * barefrac * convertg2kg
          appl_fert_t(i,m) = appl_fert_t(i,m) + appl_fertrow(i,m,iccp1) * barefrac * convertkgN
          ndep_nh4_t(i,m) = ndep_nh4_t(i,m) + ndep_nh4row(i,m,iccp1) * barefrac * convertkgN
          ndep_no3_t(i,m) = ndep_no3_t(i,m) + ndep_no3row(i,m,iccp1)  * barefrac * convertkgN
          bnf_tot_t(i,m) =  bnf_tot_t(i,m) + bnftotrow(i,m,iccp1) * barefrac  * convertkgN
          bnf_free_t(i,m) = bnf_free_t(i,m) + bnffreerow(i,m,iccp1) * barefrac * convertkgN
          nitrif_t(i,m) = nitrif_t(i,m) + nitrifvegrow(i,m,iccp1) * barefrac * convertkgN
          no_nit_t(i,m) = no_nit_t(i,m) + no_nitvegrow(i,m,iccp1) * barefrac * convertkgN
          no_denit_t(i,m) = no_denit_t(i,m) + no_denitvegrow(i,m,iccp1) * barefrac * convertkgN
          no_nitdenit_t(i,m) = no_nitdenit_t(i,m) + no_nitdenitvegrow(i,m,iccp1) *  barefrac * convertkgN
          n2o_nit_t(i,m) = n2o_nit_t(i,m) + n2o_nitvegrow(i,m,iccp1) * barefrac * convertkgN
          n2o_denit_t(i,m) = n2o_denit_t(i,m) + n2o_denitvegrow(i,m,iccp1) * barefrac * convertkgN
          n2o_nitdenit_t(i,m) = n2o_nitdenit_t(i,m) + n2o_nitdenitvegrow(i,m,iccp1) * barefrac * convertkgN
          n2_denit_t(i,m) = n2_denit_t(i,m) + n2_denitvegrow(i,m,iccp1) * barefrac * convertkgN
          nvol_t(i,m) = nvol_t(i,m) + nvolvegrow(i,m,iccp1) * barefrac * convertkgN
          nleach_t(i,m) = nleach_t(i,m) + nleachvegrow(i,m,iccp1) * barefrac * convertkgN
          nhumtrs_t(i,m) = nhumtrs_t(i,m) + nhumtrsvegrow(i,m,iccp1) * barefrac * convertkgN
          nmineral_litr_t(i,m) = nmineral_litr_t(i,m) + nmineralveg_litrrow(i,m,iccp1) * barefrac * convertkgN
          nmineral_humus_t(i,m) = nmineral_humus_t(i,m) + nmineralveg_humusrow(i,m,iccp1) * barefrac * convertkgN
          nimmobil_nh4_t(i,m) = nimmobil_nh4_t(i,m) + nimmobilveg_nh4row(i,m,iccp1) * barefrac * convertkgN
          nimmobil_no3_t(i,m) = nimmobil_no3_t(i,m) + nimmobilveg_no3row(i,m,iccp1) * barefrac * convertkgN
          fNnetland_t(i,m) = fNnetland_t(i,m) + fNnetlandvegrow(i,m,iccp1) * barefrac * convertkgN
          if ((ngleafmas_t(i,m) + nbleafmas_t(i,m)) /= 0.0) c2n_l_t(i,m) = &
             (gleafmas_t(i,m) + bleafmas_t(i,m)) / (ngleafmas_t(i,m) + nbleafmas_t(i,m))
          if (nstemmass_t(i,m) /= 0.0) c2n_s_t(i,m) = stemmass_t(i,m) / nstemmass_t(i,m)
          if (nrootmass_t(i,m) /= 0.0) c2n_r_t(i,m) = rootmass_t(i,m) / nrootmass_t(i,m)
          if ((ngleafmas_t(i,m) + nbleafmas_t(i,m) + nstemmass_t(i,m) + nrootmass_t(i,m)) /= 0.0) &
             c2n_wp_t(i,m) = (gleafmas_t(i,m) + bleafmas_t(i,m) + stemmass_t(i,m) + rootmass_t(i,m)) &
                             / (ngleafmas_t(i,m) + nbleafmas_t(i,m) + nstemmass_t(i,m) + nrootmass_t(i,m))
          if (nlitrmass_t(i,m) /= 0.0) c2n_litr_t(i,m) = sum(litrmass_t(i,m,:)) / nlitrmass_t(i,m)
          if (soilnmas_t(i,m) /= 0.0) c2n_humus_t(i,m) = sum(soilcmas_t(i,m,:)) / soilnmas_t(i,m)
        end if !for Ncycle_on

        !> Calculation of grid averaged variables

        gpp_g(i) = gpp_g(i) + gpprow(i,m) * FAREROT(i,m)
        npp_g(i) = npp_g(i) + npprow(i,m) * FAREROT(i,m)
        nep_g(i) = nep_g(i) + neprow(i,m) * FAREROT(i,m)
        nepCMIP_g(i) = nepCMIP_g(i) + nepCMIProw(i,m) * FAREROT(i,m)
        ! NOTE: This NBP will include LUC product pool contributions
        ! since it is the nbprow varialbe. nbpveg variable does not
        ! include LUC contributions since they exist at the tile level
        ! not at the PFT level.
        nbp_g(i) = nbp_g(i) + nbprow(i,m) * FAREROT(i,m)
        autores_g(i) = autores_g(i) + autoresrow(i,m) * FAREROT(i,m)
        hetrores_g(i) = hetrores_g(i) + hetroresrow(i,m) * FAREROT(i,m)
        litres_g(i) = litres_g(i) + litresrow(i,m) * FAREROT(i,m)
        socres_g(i) = socres_g(i) + socresrow(i,m) * FAREROT(i,m)
        dstcemls_g(i) = dstcemls_g(i) + dstcemlsrow(i,m) * FAREROT(i,m)
        dstcemls3_g(i) = dstcemls3_g(i) + dstcemls3row(i,m) * FAREROT(i,m)
        litrfall_g(i) = litrfall_g(i) + litrfallrow(i,m) * FAREROT(i,m) !(NOTE:umol CO_2 /m2/s, litrfallveg is kgC/m2/s)
        humiftrs_g(i) = humiftrs_g(i) + humiftrsrow(i,m) * FAREROT(i,m)
        rml_g(i) = rml_g(i) + rmlrow(i,m) * FAREROT(i,m)
        rms_g(i) = rms_g(i) + rmsrow(i,m) * FAREROT(i,m)
        rmr_g(i) = rmr_g(i) + rmrrow(i,m) * FAREROT(i,m)
        rg_g(i) = rg_g(i) + rgrow(i,m) * FAREROT(i,m)
        leaflitr_g(i) = leaflitr_g(i) + leaflitr_t(i,m) * FAREROT(i,m)
        tltrleaf_g(i) = tltrleaf_g(i) + tltrleaf_t(i,m) * FAREROT(i,m)
        tltrstem_g(i) = tltrstem_g(i) + tltrstem_t(i,m) * FAREROT(i,m)
        tltrroot_g(i) = tltrroot_g(i) + tltrroot_t(i,m) * FAREROT(i,m)
        slai_g(i) = slai_g(i) + slai_t(i,m) * FAREROT(i,m)
        ailcg_g(i) = ailcg_g(i) + ailcg_t(i,m) * FAREROT(i,m)
        ailcb_g(i) = ailcb_g(i) + ailcb_t(i,m) * FAREROT(i,m)
        vgbiomas_g(i) = vgbiomas_g(i) + vgbiomasrow(i,m) * FAREROT(i,m)
        veghght_g(i) = veghght_g(i) + veghght_t(i,m) * FAREROT(i,m)
        gavglai_g(i) = gavglai_g(i) + gavglairow(i,m) * FAREROT(i,m)
        gavgltms_g(i) = gavgltms_g(i) + gavgltmsrow(i,m) * FAREROT(i,m)
        gavgscms_g(i) = gavgscms_g(i) + gavgscmsrow(i,m) * FAREROT(i,m)
        totcmass_g(i) = vgbiomas_g(i) + gavgltms_g(i) + gavgscms_g(i)
        leafns2s_g(i) = leafns2s_g(i) + leafns2s_t(i,m) * FAREROT(i,m)
        stemns2s_g(i) = stemns2s_g(i) + stemns2s_t(i,m) * FAREROT(i,m)
        rootns2s_g(i) = rootns2s_g(i) + rootns2s_t(i,m) * FAREROT(i,m)
        re_alloc_s2l_g(i) = re_alloc_s2l_g(i) + re_alloc_s2l_t(i,m) * FAREROT(i,m)
        re_alloc_r2l_g(i) = re_alloc_r2l_g(i) + re_alloc_r2l_t(i,m) * FAREROT(i,m)
        re_alloc_sr2l_g(i) = re_alloc_sr2l_g(i) + re_alloc_sr2l_t(i,m) * FAREROT(i,m)
        gleafmas_g(i) = gleafmas_g(i) + gleafmas_t(i,m) * FAREROT(i,m)
        gleafmas_ns_g(i) = gleafmas_ns_g(i) + gleafmas_ns_t(i,m) * FAREROT(i,m)
        gleafmas_s_g(i) = gleafmas_s_g(i) + gleafmas_s_t(i,m) * FAREROT(i,m)
        bleafmas_g(i) = bleafmas_g(i) + bleafmas_t(i,m) * FAREROT(i,m)
        stemmass_g(i) = stemmass_g(i) + stemmass_t(i,m) * FAREROT(i,m)
        stemmass_ns_g(i) = stemmass_ns_g(i) + stemmass_ns_t(i,m) * FAREROT(i,m)
        stemmass_s_g(i) = stemmass_s_g(i) + stemmass_s_t(i,m) * FAREROT(i,m)
        rootmass_g(i) = rootmass_g(i) + rootmass_t(i,m) * FAREROT(i,m)
        rootmass_ns_g(i) = rootmass_ns_g(i) + rootmass_ns_t(i,m) * FAREROT(i,m)
        rootmass_s_g(i) = rootmass_s_g(i) + rootmass_s_t(i,m) * FAREROT(i,m)
        rootdpth_g(i) = rootdpth_g(i) + rootdpth_t(i,m) * FAREROT(i,m)
        roottemp_g(i) = roottemp_g(i) + roottemp_t(i,m) * FAREROT(i,m)
        burnfrac_g(i) = burnfrac_g(i) + burnfracrow(i,m) * FAREROT(i,m)
        smfuncveg_g(i) = smfuncveg_g(i) + smfuncveg_t(i,m) * FAREROT(i,m)
        lucemcom_g(i) = lucemcom_g(i) + lucemcomrow(i,m) * FAREROT(i,m)
        lucltrin_g(i) = lucltrin_g(i) + lucltrinrow(i,m) * FAREROT(i,m)
        lucsocin_g(i) = lucsocin_g(i) + lucsocinrow(i,m) * FAREROT(i,m)
        if (Ncycle_on) then
          lucemcomn_g(i) = lucemcomn_g(i) + lucemcomnrow(i,m) * FAREROT(i,m)
          lucltrinn_g(i) = lucltrinn_g(i) + lucltrinnrow(i,m) * FAREROT(i,m)
          lucsocinn_g(i) = lucsocinn_g(i) + lucsocinnrow(i,m) * FAREROT(i,m)
        end if
        ch4WetSpec_g(i) = ch4WetSpec_g(i) + ch4WetSpecrow(i,m) * farerot(i,m)
        wetfdyn_g(i) = wetfdyn_g(i) + wetfdynrow(i,m) * farerot(i,m)
        ch4WetDyn_g(i) = ch4WetDyn_g(i) + ch4WetDynrow(i,m) * farerot(i,m)
        ch4soills_g(i) = ch4soills_g(i) + ch4soillsrow(i,m) * farerot(i,m)
        emit_co2_g(i) = emit_co2_g(i) + emit_co2_t(i,m) * FAREROT(i,m)
        emit_ch4_g(i) = emit_ch4_g(i) + emit_ch4_t(i,m) * FAREROT(i,m)
        ! nppmoss_g(i)  = nppmoss_g(i) +nppmossrow(i,m)*FAREROT(i,m)
        ! armoss_g(i)   = armoss_g(i) + armossrow(i,m)*FAREROT(i,m)

        do k = 1,ignd
          rmatctem_g(i,k) = rmatctem_g(i,k) + rmatctem_t(i,m,k) * FAREROT(i,m)
          ! FLAG not putting as per layer yet since it is not presently written out.
          litrmass_g(i) = litrmass_g(i) + litrmass_t(i,m,k) * FAREROT(i,m)
          soilcmas_g(i) = soilcmas_g(i) + soilcmas_t(i,m,k) * FAREROT(i,m)

        end do
        if (Ncycle_on) then
           bnf_tot_g(i) = bnf_tot_g(i) + bnf_tot_t(i,m) * FAREROT(i,m)
           bnf_free_g(i) = bnf_free_g(i) + bnf_free_t(i,m) * FAREROT(i,m)
           bnf_ant_g(i) = bnf_ant_g(i) + bnf_ant_t(i,m) * FAREROT(i,m)
           bnf_nat_g(i) = bnf_nat_g(i) + bnf_nat_t(i,m) * FAREROT(i,m)
           nh4_mass_g(i) = nh4_mass_g(i) + nh4_mass_t(i,m) * FAREROT(i,m)
           no3_mass_g(i) = no3_mass_g(i) + no3_mass_t(i,m) * FAREROT(i,m)
           nlitrmass_g(i) = nlitrmass_g(i) + nlitrmass_t(i,m) * FAREROT(i,m)
           soilnmas_g(i) = soilnmas_g(i) + soilnmas_t(i,m) * FAREROT(i,m)
           ngleafmas_g(i) = ngleafmas_g(i) + ngleafmas_t(i,m) * FAREROT(i,m)
           ngleafmas_ns_g(i) = ngleafmas_ns_g(i) + ngleafmas_ns_t(i,m) * FAREROT(i,m)
           ngleafmas_s_g(i) = ngleafmas_s_g(i) + ngleafmas_s_t(i,m) * FAREROT(i,m)
           nbleafmas_g(i) = nbleafmas_g(i) + nbleafmas_t(i,m) * FAREROT(i,m)
           nstemmass_g(i) = nstemmass_g(i) + nstemmass_t(i,m) * FAREROT(i,m)
           nstemmass_ns_g(i) = nstemmass_ns_g(i) + nstemmass_ns_t(i,m) * FAREROT(i,m)
           nstemmass_s_g(i) = nstemmass_s_g(i) + nstemmass_s_t(i,m)  * FAREROT(i,m)
           nrootmass_g(i) = nrootmass_g(i) + nrootmass_t(i,m) * FAREROT(i,m)
           nrootmass_ns_g(i) = nrootmass_ns_g(i) + nrootmass_ns_t(i,m) * FAREROT(i,m)
           nrootmass_s_g(i) = nrootmass_s_g(i) + nrootmass_s_t(i,m) * FAREROT(i,m)
           nvgbiomas_g(i) = nvgbiomas_g(i) + nvgbiomas_t(i,m) * FAREROT(i,m)
           nitrif_g(i) = nitrif_g(i) + nitrif_t(i,m) * FAREROT(i,m)
           no_nit_g(i) = no_nit_g(i) + no_nit_t(i,m) * FAREROT(i,m)
           no_denit_g(i) = no_denit_g(i) + no_denit_t(i,m) * FAREROT(i,m)
           no_nitdenit_g(i) = no_nitdenit_g(i) + no_nitdenit_t(i,m) * FAREROT(i,m)
           n2o_nit_g(i) = n2o_nit_g(i) + n2o_nit_t(i,m) * FAREROT(i,m)
           n2o_denit_g(i) = n2o_denit_g(i) + n2o_denit_t(i,m) * FAREROT(i,m)
           n2o_nitdenit_g(i) = n2o_nitdenit_g(i) + n2o_nitdenit_t(i,m) * FAREROT(i,m)
           n2_denit_g(i) = n2_denit_g(i) + n2_denit_t(i,m) * FAREROT(i,m)
           nvol_g(i) = nvol_g(i) + nvol_t(i,m) * FAREROT(i,m)
           nleach_g(i) = nleach_g(i) + nleach_t(i,m) * FAREROT(i,m)
           appl_fert_g(i) = appl_fert_g(i) + appl_fert_t(i,m) * FAREROT(i,m)
           ndep_nh4_g(i) = ndep_nh4_g(i) + ndep_nh4_t(i,m) * FAREROT(i,m)
           ndep_no3_g(i) = ndep_no3_g(i) + ndep_no3_t(i,m) * FAREROT(i,m)
           ndemand_wp_npp_g(i) = ndemand_wp_npp_g(i) + ndemand_wp_npp_t(i,m) * FAREROT(i,m)
           nuptake_p_nh4_g(i) = nuptake_p_nh4_g(i) + nuptake_p_nh4_t(i,m) * FAREROT(i,m)
           nuptake_p_no3_g(i) = nuptake_p_no3_g(i) + nuptake_p_no3_t(i,m) * FAREROT(i,m)
           nuptake_a_actl_nh4_g(i) = nuptake_a_actl_nh4_g(i) + nuptake_a_actl_nh4_t(i,m) * FAREROT(i,m)
           nuptake_a_actl_no3_g(i) = nuptake_a_actl_no3_g(i) + nuptake_a_actl_no3_t(i,m) * FAREROT(i,m)
           nalloc_l_g(i) = nalloc_l_g(i) + nalloc_l_t(i,m) * FAREROT(i,m)
           nalloc_s_g(i) = nalloc_s_g(i) + nalloc_s_t(i,m) * FAREROT(i,m)
           nalloc_r_g(i) = nalloc_r_g(i) + nalloc_r_t(i,m) * FAREROT(i,m)
           nresorped_s_g(i) = nresorped_s_g(i) + nresorped_s_t(i,m) * FAREROT(i,m)
           nresorped_r_g(i) = nresorped_r_g(i) + nresorped_r_t(i,m) * FAREROT(i,m)
           nre_alloc_s2l_g(i) = nre_alloc_s2l_g(i) + nre_alloc_s2l_t(i,m) * FAREROT(i,m)
           nre_alloc_r2l_g(i) = nre_alloc_r2l_g(i) + nre_alloc_r2l_t(i,m) * FAREROT(i,m)
           nleafns2s_g(i) = nleafns2s_g(i) + nleafns2s_t(i,m) * FAREROT(i,m)
           nstemns2s_g(i) = nstemns2s_g(i) + nstemns2s_t(i,m) * FAREROT(i,m)
           nrootns2s_g(i) = nrootns2s_g(i) + nrootns2s_t(i,m) * FAREROT(i,m)

           nlitr_l_g(i) = nlitr_l_g(i) + nlitr_l_t(i,m) * FAREROT(i,m)
           nlitr_s_g(i) = nlitr_s_g(i) + nlitr_s_t(i,m) * FAREROT(i,m)
           nlitr_r_g(i) = nlitr_r_g(i) + nlitr_r_t(i,m) * FAREROT(i,m)
           gl2bl_grass_nflux_g(i) = gl2bl_grass_nflux_g(i)+ &
                    gl2bl_grass_nflux_t(i,m) * FAREROT(i,m)
           nhumtrs_g(i) = nhumtrs_g(i) + nhumtrs_t(i,m) * FAREROT(i,m)
           nmineral_litr_g(i) = nmineral_litr_g(i) + nmineral_litr_t(i,m) * FAREROT(i,m)
           nmineral_humus_g(i) = nmineral_humus_g(i) + nmineral_humus_t(i,m) * FAREROT(i,m)
           nimmobil_nh4_g(i) = nimmobil_nh4_g(i) + nimmobil_nh4_t(i,m) * FAREROT(i,m)
           nimmobil_no3_g(i) = nimmobil_no3_g(i) + nimmobil_no3_t(i,m) * FAREROT(i,m)
           fNnetland_g(i) = fNnetland_g(i) + fNnetland_t(i,m) * FAREROT(i,m)
        end if !for Ncycle_on

      end do ! loop 70 ! nmtest

      if (Ncycle_on) then
         if ((ngleafmas_g(i) + nbleafmas_g(i)) /= 0.0) c2n_l_g(i) = &
             (gleafmas_g(i) + bleafmas_g(i)) / (ngleafmas_g(i) + nbleafmas_g(i))
         if (nstemmass_g(i) /= 0.0) c2n_s_g(i) = stemmass_g(i) / nstemmass_g(i)
         if (nrootmass_g(i) /= 0.0) c2n_r_g(i) = rootmass_g(i) / nrootmass_g(i)
         if ((ngleafmas_g(i) + nbleafmas_g(i) + nstemmass_g(i) + nrootmass_g(i)) /= 0.0) &
              c2n_wp_g(i) = (gleafmas_g(i) + bleafmas_g(i) + stemmass_g(i) + rootmass_g(i)) &
                           / (ngleafmas_g(i) + nbleafmas_g(i) + nstemmass_g(i) + nrootmass_g(i))
         if (nlitrmass_g(i) /= 0.0) c2n_litr_g(i) = litrmass_g(i) / nlitrmass_g(i)
         if (soilnmas_g(i) /= 0.0) c2n_humus_g(i) = soilcmas_g(i) / soilnmas_g(i)
      end if

    end do ! loop 60 ! nltest

    i = 1 ! offline nltest is always 1.

    ! Transfer the consecDays to timeStamp (since we need a size 1 array)
    timeStamp = consecDays

    !> Write grid average values

    call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_d_g' ,timeStamp,'gpp', [gpp_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'leafns2s_d_g' ,timeStamp,'leafns2s', [leafns2s_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'stemns2s_d_g' ,timeStamp,'stemns2s', [stemns2s_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'rootns2s_d_g' ,timeStamp,'rootns2s', [rootns2s_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_s2l_d_g' ,timeStamp,'realloc_s2l', [re_alloc_s2l_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_r2l_d_g' ,timeStamp,'realloc_r2l', [re_alloc_r2l_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_sr2l_d_g' ,timeStamp,'realloc_sr2l', [re_alloc_sr2l_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_d_g' ,timeStamp,'npp', [npp_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_d_g' ,timeStamp,'nep', [nep_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'nepCMIP_d_g' ,timeStamp,'nepCMIP', [nepCMIP_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_d_g' ,timeStamp,'nbp', [nbp_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_d_g' ,timeStamp,'ra', [autores_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_d_g' ,timeStamp,'rh', [hetrores_g(i)])
    if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_d_g' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_g(i)])
    end if
    call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_d_g' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_d_g' ,timeStamp,'soilCH4cons',[ch4soills_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_d_g' ,timeStamp,'fFireCH4',[emit_ch4_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_d_g' ,timeStamp,'wetlandFrac',[wetfdyn_g(i)])
    ! litres_g(i),
    ! socres_g(i), &
    !         (dstcemls_g(i)+dstcemls3_g(i)), &
    !         litrfall_g(i),
    ! humiftrs_g(i),' GRDAV'
    call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_d_g',timeStamp,'cLeaf',[gleafmas_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_NS_d_g',timeStamp,'cLeaf_ns',[gleafmas_ns_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmass_d_g',timeStamp,'cLeaf_s',[gleafmas_s_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'bleafmas_d_g' ,timeStamp,'bleafmas', [bleafmas_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_d_g',timeStamp,'cStem',[stemmass_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_NS_d_g',timeStamp,'cStem_ns',[stemmass_ns_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmasss_d_g',timeStamp,'cStem_s',[stemmass_s_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_d_g',timeStamp,'cRoot',[rootmass_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_NS_d_g',timeStamp,'cRoot_ns',[rootmass_ns_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmasss_d_g',timeStamp,'cRoot_s',[rootmass_s_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'ailcg_d_g',timeStamp,'lai',[ailcg_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_d_g',timeStamp,'burntFractionAll',[burnfrac_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_d_g' ,timeStamp,'vegHeight', [veghght_g(i)])
    if (Ncycle_on) then
       call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_d_g',timeStamp,'nLeaf',[ngleafmas_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_NS_d_g',timeStamp,'nLeaf_ns',[ngleafmas_ns_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmass_d_g',timeStamp,'nLeaf_s',[ngleafmas_s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_d_g',timeStamp,'nbLeaf',[nbleafmas_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_d_g',timeStamp,'nStem',[nstemmass_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_NS_d_g',timeStamp,'nStem_ns',[nstemmass_ns_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmasss_d_g',timeStamp,'nStem_s',[nstemmass_s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_d_g',timeStamp,'nRoot',[nrootmass_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_NS_d_g',timeStamp,'nRoot_ns',[nrootmass_ns_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmasss_d_g',timeStamp,'nRoot_s',[nrootmass_s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'ndemand_wp_npp_d_g' ,timeStamp,'ndemand_npp', [ndemand_wp_npp_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_d_g' ,timeStamp,'c2n_l', [c2n_l_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_s_d_g' ,timeStamp,'c2n_s', [c2n_s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_r_d_g' ,timeStamp,'c2n_r', [c2n_r_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_wp_d_g' ,timeStamp,'c2n_wp', [c2n_wp_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_litr_d_g' ,timeStamp,'c2n_litr', [c2n_litr_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_humus_d_g' ,timeStamp,'c2n_humus', [c2n_humus_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_d_g',timeStamp,'nh4',[nh4_mass_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_d_g',timeStamp,'no3',[no3_mass_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_d_g',timeStamp,'nLitter',[nlitrmass_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_d_g',timeStamp,'nSoil',[soilnmas_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_d_g',timeStamp,'nVeg',[nvgbiomas_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_d_g' ,timeStamp,'bnf_tot', [bnf_tot_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_d_g' ,timeStamp,'bnf_free', [bnf_free_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_d_g' ,timeStamp,'bnf_ant', [bnf_ant_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_d_g' ,timeStamp,'bnf_nat', [bnf_nat_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_d_g' ,timeStamp,'nitrif', [nitrif_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nit_d_g' ,timeStamp,'no_nit', [no_nit_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_d_g' ,timeStamp,'no_denit', [no_denit_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_d_g' ,timeStamp,'no', [no_nitdenit_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nit_d_g' ,timeStamp,'n2o_nit', [n2o_nit_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_d_g' ,timeStamp,'n2o_denit', [n2o_denit_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_d_g' ,timeStamp,'n2o', [n2o_nitdenit_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_d_g' ,timeStamp,'n2', [n2_denit_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_d_g' ,timeStamp,'nvol', [nvol_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_d_g' ,timeStamp,'nleach', [nleach_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'appl_fert_d_g' ,timeStamp,'nfer_nh4', [appl_fert_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_nh4_d_g' ,timeStamp,'ndep_nh4', [ndep_nh4_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_no3_d_g' ,timeStamp,'ndep_no3', [ndep_no3_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_nh4_d_g' ,timeStamp,'nuptake_p_nh4', [nuptake_p_nh4_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_no3_d_g' ,timeStamp,'nuptake_p_no3', [nuptake_p_no3_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_nh4_d_g' ,timeStamp,'nuptake_a_actl_nh4', [nuptake_a_actl_nh4_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_no3_d_g' ,timeStamp,'nuptake_a_actl_no3', [nuptake_a_actl_no3_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_l_d_g' ,timeStamp,'nalloc_l', [nalloc_l_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_s_d_g' ,timeStamp,'nalloc_s', [nalloc_s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_r_d_g' ,timeStamp,'nalloc_r', [nalloc_r_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_d_g' ,timeStamp,'nresorped_s', [nresorped_s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_d_g' ,timeStamp,'nresorped_r', [nresorped_r_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_s2l_d_g' ,timeStamp,'nre_alloc_s2l', [nre_alloc_s2l_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_r2l_d_g' ,timeStamp,'nre_alloc_r2l', [nre_alloc_r2l_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nleafns2s_d_g' ,timeStamp,'nleafns2s', [nleafns2s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemns2s_d_g' ,timeStamp,'nstemns2s', [nstemns2s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootns2s_d_g' ,timeStamp,'nrootns2s', [nrootns2s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_d_g' ,timeStamp,'nlitr_l', [nlitr_l_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_s_d_g' ,timeStamp,'nlitr_s', [nlitr_s_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_r_d_g' ,timeStamp,'nlitr_r', [nlitr_r_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'gl2bl_grass_nflux_d_g' ,timeStamp,'gl2bl_grass_nflux', [gl2bl_grass_nflux_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nhumtrs_d_g' ,timeStamp,'nhumtrs', [nhumtrs_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_litr_d_g' ,timeStamp,'nmineral_litr', [nmineral_litr_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_humus_d_g' ,timeStamp,'nmineral_humus', [nmineral_humus_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_nh4_d_g' ,timeStamp,'nimmobil_nh4', [nimmobil_nh4_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_no3_d_g' ,timeStamp,'nimmobil_no3', [nimmobil_no3_g(i)])
       call writeOutput1D(lonLocalIndex,latLocalIndex,'fNnetland_d_g' ,timeStamp,'fNnetland', [fNnetland_g(i)])
    end if

    if (doperpftoutput) then

      !                     !> File: .CT01D
      !                     write(72,8200)iday,realyr,gppvegrow(i,m,j),nppvegrow(i,m,j), &
      !                     nepvegrow(i,m,j),nbpvegrow(i,m,j),autoresvegrow(i,m,j), &
      !                     hetroresvegrow(i,m,j),litresvegrow(i,m,j),soilcresvegrow(i,m,j), &
      !                     (dstcemlsrow(i,m)+dstcemls3row(i,m)), &   ! FLAG at present dstcemls are only per tile values
      !                     litrfallrow(i,m),humiftrsrow(i,m), & ! same with litrfall and humiftrs.
      !                     ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
      !                 !> File: .CT01D
      !                 write(72,8200)iday,realyr,0.0,0.0, &
      !                 nepvegrow(i,m,iccp1),nbpvegrow(i,m,iccp1),0.0, &
      !                 hetroresvegrow(i,m,iccp1),litresvegrow(i,m,iccp1),soilcresvegrow(i,m,iccp1), &
      !                 (dstcemlsrow(i,m)+dstcemls3row(i,m)), &   ! FLAG at present dstcemls are only per tile values
      !                 litrfallrow(i,m),humiftrsrow(i,m), & ! same with litrfall and humiftrs.
      !                 ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac

            if (nmtest > 1) then
                print*,'Per PFT and per tile outputs together not implemented yet'
            else
                m = 1
                call writeOutput1D(lonLocalIndex,latLocalIndex,'lfstatus_d',timeStamp,'lfstatus',[real(lfstatusrow(i,m,:))])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_d',timeStamp,'cLeaf',[gleafmasrow(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_d',timeStamp,'cStem',[stemmassrow(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_d',timeStamp,'cRoot',[rootmassrow(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'ailcg_d',timeStamp,'lai',[ailcgrow(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_d',timeStamp,'fVegLitter',&
                  [(tltrleafrow(i,m,:) + tltrstemrow(i,m,:) + tltrrootrow(i,m,:)) * convertkgC]) ! convert from umol CO2/m2/s to kgC/m2/s
                
                call writeOutput1D(lonLocalIndex,latLocalIndex,'leaflitr_d',timeStamp,'fLeafLitter',&
                    [tltrleafrow(i,m,:) * convertkgC]) ! convert from umol CO2/m2/s to kgC/m2/s
                call writeOutput1D(lonLocalIndex,latLocalIndex,'stemlitr_d',timeStamp,'fStemLitter',&
                    [tltrstemrow(i,m,:) * convertkgC]) ! convert from umol CO2/m2/s to kgC/m2/s
                call writeOutput1D(lonLocalIndex,latLocalIndex,'rootlitr_d',timeStamp,'fRootLitter',&
                    [tltrrootrow(i,m,:) * convertkgC]) ! convert from umol CO2/m2/s to kgC/m2/s
                call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_d',timeStamp,'npp',[nppvegrow(i,m,:)])
                    
               

                !do k = 1,iccp2
                !  call writeOutput1D(lonLocalIndex,latLocalIndex,'litrMassPerLay_d',timeStamp,'cLitterperlay',[litrmassrow(i,m,k,:)])
                !  call writeOutput1D(lonLocalIndex,latLocalIndex,'soilCMasPerLay_d',timeStamp,'cSoilperlay',[soilcmasrow(i,m,k,:)])
                !end do

                if (Ncycle_on) then
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_d',timeStamp,'nLeaf',[ngleafmasrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_d',timeStamp,'nbLeaf',[nbleafmasrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_d',timeStamp,'nStem',[nstemmassrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_d',timeStamp,'nRoot',[nrootmassrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_d',timeStamp,'nh4',[nh4_massrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_d',timeStamp,'nLitter',[nlitrmassrow(i,m,1:iccp1)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_d',timeStamp,'nSoil',[soilnmasrow(i,m,1:iccp1)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_d' ,timeStamp,'nlitr_l', [nlitrveg_lrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_d' ,timeStamp,'nresorped_s', [nresorpedveg_srow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_d' ,timeStamp,'nresorped_r', [nresorpedveg_rrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_d',timeStamp,'bnf_nat',[bnfnatrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_d',timeStamp,'bnf_ant',[bnfantrow(i,m,:)])
                  call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_d',timeStamp,'bnf_free',[bnffreerow(i,m,:)])
                end if ! Ncycle_on
            end if

    end if

    if (dopertileoutput) then
      !                 !> File: .CT01D
      !                 write(72,8200)iday,realyr,gpprow(i,m),npprow(i,m), &
      !                 neprow(i,m),nbprow(i,m),autoresrow(i,m), &
      !                 hetroresrow(i,m),litresrow(i,m),socresrow(i,m), &
      !                 (dstcemlsrow(i,m)+dstcemls3row(i,m)), &
      !                 litrfallrow(i,m),humiftrsrow(i,m), &
      !                 ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)

      if (nmtest > 1) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_d_t',timeStamp,'cLeaf',[gleafmas_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_NS_d_t',timeStamp,'cLeaf_ns',[gleafmas_ns_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmass_d_t',timeStamp,'cLeaf_s',[gleafmas_s_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_d_t',timeStamp,'cStem',[stemmass_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_NS_d_t',timeStamp,'cStem_ns',[stemmass_ns_T(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmasss_d_t',timeStamp,'cStem_s',[stemmass_s_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_d_t',timeStamp,'cRoot',[rootmass_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_NS_d_t',timeStamp,'cRoot_ns',[rootmass_ns_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmasss_d_t',timeStamp,'cRoot_s',[rootmass_s_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'leafns2s_t',timeStamp,'leafns2s',[leafns2s_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemns2s_t',timeStamp,'stemns2s',[stemns2s_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootns2s_t',timeStamp,'rootns2s',[rootns2s_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_s2l_t',timeStamp,'realloc_s2l',[re_alloc_s2l_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_r2l_t',timeStamp,'realloc_r2l',[re_alloc_r2l_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_sr2l_t',timeStamp,'realloc_sr2l',[re_alloc_sr2l_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ailcg_d_t',timeStamp,'lai',[ailcg_t(i,:)])
          if (Ncycle_on) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_d_t',timeStamp,'nLeaf',[ngleafmas_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_NS_d_t',timeStamp,'nLeaf_ns',[ngleafmas_ns_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmass_d_t',timeStamp,'nLeaf_s',[ngleafmas_s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_d_t',timeStamp,'nbLeaf',[nbleafmas_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_d_t',timeStamp,'nStem',[nstemmass_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_NS_d_t',timeStamp,'nStem_ns',[nstemmass_ns_T(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmasss_d_t',timeStamp,'nStem_s',[nstemmass_s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_d_t',timeStamp,'nRoot',[nrootmass_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_NS_d_t',timeStamp,'nRoot_ns',[nrootmass_ns_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmasss_d_t',timeStamp,'nRoot_s',[nrootmass_s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_d_t',timeStamp,'nVeg',[nvgbiomas_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ndemand_wp_npp_d_t',timeStamp,'ndemand_npp',[ndemand_wp_npp_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_d_t' ,timeStamp,'c2n_l', [c2n_l_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_s_d_t' ,timeStamp,'c2n_s', [c2n_s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_r_d_t' ,timeStamp,'c2n_r', [c2n_r_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_wp_d_t' ,timeStamp,'c2n_wp', [c2n_wp_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_litr_d_t' ,timeStamp,'c2n_litr', [c2n_litr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_humus_d_t' ,timeStamp,'c2n_humus', [c2n_humus_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_d_t',timeStamp,'nLitter',[nlitrmass_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_d_t',timeStamp,'nSoil',[soilnmas_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_d_t',timeStamp,'nh4',[nh4_mass_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_d_t',timeStamp,'no3',[no3_mass_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_d_t',timeStamp,'bnf_tot',[bnf_tot_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_d_t',timeStamp,'bnf_free',[bnf_free_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_d_t',timeStamp,'bnf_ant',[bnf_ant_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_d_t',timeStamp,'bnf_nat',[bnf_nat_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_d_t',timeStamp,'nitrif',[nitrif_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nit_d_t',timeStamp,'no_nit',[no_nit_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_d_t',timeStamp,'no_denit',[no_denit_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_d_t',timeStamp,'no',[no_nitdenit_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nit_d_t',timeStamp,'n2o_nit',[n2o_nit_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_d_t',timeStamp,'n2o_denit',[n2o_denit_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_d_t',timeStamp,'n2o',[n2o_nitdenit_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_d_t',timeStamp,'n2',[n2_denit_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_d_t',timeStamp,'nvol',[nvol_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_d_t',timeStamp,'nleach',[nleach_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'appl_fert_d_t',timeStamp,'nfer_nh4',[appl_fert_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_nh4_d_t',timeStamp,'ndep_nh4',[ndep_nh4_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_no3_d_t',timeStamp,'ndep_no3',[ndep_no3_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_nh4_d_t',timeStamp,'nuptake_p_nh4',[nuptake_p_nh4_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_no3_d_t',timeStamp,'nuptake_p_no3',[nuptake_p_no3_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_nh4_d_t',timeStamp,'nuptake_a_actl_nh4',[nuptake_a_actl_nh4_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_no3_d_t',timeStamp,'nuptake_a_actl_no3',[nuptake_a_actl_no3_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_l_d_t' ,timeStamp,'nalloc_l', [nalloc_l_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_s_d_t' ,timeStamp,'nalloc_s', [nalloc_s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_r_d_t' ,timeStamp,'nalloc_r', [nalloc_r_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_d_t' ,timeStamp,'nresorped_s', [nresorped_s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_d_t' ,timeStamp,'nresorped_r', [nresorped_r_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_s2l_d_t' ,timeStamp,'nre_alloc_s2l', [nre_alloc_s2l_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_r2l_d_t' ,timeStamp,'nre_alloc_r2l', [nre_alloc_r2l_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nleafns2s_d_t' ,timeStamp,'nleafns2s', [nleafns2s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemns2s_d_t' ,timeStamp,'nstemns2s', [nstemns2s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootns2s_d_t' ,timeStamp,'nrootns2s', [nrootns2s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_d_t' ,timeStamp,'nlitr_l', [nlitr_l_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_s_d_t' ,timeStamp,'nlitr_s', [nlitr_s_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_r_d_t' ,timeStamp,'nlitr_r', [nlitr_r_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gl2bl_grass_nflux_d_t' ,timeStamp,'gl2bl_grass_nflux', [gl2bl_grass_nflux_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nhumtrs_d_t' ,timeStamp,'nhumtrs', [nhumtrs_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_litr_d_t' ,timeStamp,'nmineral_litr', [nmineral_litr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_humus_d_t' ,timeStamp,'nmineral_humus', [nmineral_humus_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_nh4_d_t' ,timeStamp,'nimmobil_nh4', [nimmobil_nh4_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_no3_d_t' ,timeStamp,'nimmobil_no3', [nimmobil_no3_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'fNnetland_d_t' ,timeStamp,'fNnetland', [fNnetland_t(i,:)])
          end if
      end if
    end if

    end associate
  end subroutine ctem_daily_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_ctem_monthly_aw
  !> @{
  !> Accumulate and write out the monthly CTEM outputs. These are kept in pointer structures as
  !! this subroutine is called daily and we increment the daily values to produce a monthly value.
  !! The pointer to the monthly data structures (in ctemStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine ctem_monthly_aw (lonLocalIndex, latLocalIndex, nltest, nmtest, iday, realyr, nday, lastDOY)

    ! J. Melton Feb 2016.

    use classStateVars, only : class_rot
    use ctemStateVars,  only : ctem_tile_mo, vrot, ctem_grd_mo, c_switch, &
                               resetMonthEnd, ctem_mo, tracer
    use classicParams,  only : icc, iccp1, nmon, mmday, monthend, monthdays, iccp2, ignd, convertkgN, convertg2kg, convertkgC
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    integer, intent(in) :: realyr
    integer, intent(in) :: nday
    integer, intent(in) :: lastDOY

    ! local
    integer :: i, m, j, nt, k
    real :: barefrac
    real :: sumfare
    integer :: NDMONTH
    integer :: imonth
    real, dimension(1) :: timeStamp
    real, dimension(icc) :: pftExist
    real :: oneOverDPM
    real, dimension(1) :: bulkLitterCarbon_mo_g !< Temporary variable used to produce the bulk soil litter carbon quantity for output at the grid level \f$[kg C m^{-2}]\f$
    real, dimension(1) :: bulkSoilCarbon_mo_g   !< Temporary variable used to produce the bulk soil carbon quantity for output at the grid level \f$[kg C m^{-2}]\f$
    real, dimension(1) :: bulkLitterResp_mo_g   !< Temporary variable used to produce the bulk litter respiration quantity for output at the grid level \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(1) :: bulkSoilResp_mo_g     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output at the grid level \f$[kg C m^{-2} s^{-1}]\f$

    real, dimension(iccp1) :: bulkLitterCarbon_mo !< Temporary variable used to produce the bulk soil litter carbon quantity for output at the monthly timestep \f$[kg C m^{-2}]\f$
    real, dimension(iccp1) :: bulkSoilCarbon_mo   !< Temporary variable used to produce the bulk soil carbon quantity for output at the monthly timestep \f$[kg C m^{-2}]\f$
    real, dimension(iccp1) :: bulkLitterResp_mo   !< Temporary variable used to produce the bulk litter respiration quantity for output at the monthly timestep \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(iccp1) :: bulkSoilResp_mo     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output at the monthly timestep \f$[kg C m^{-2} s^{-1}]\f$

    real, dimension(nmtest) :: bulkLitterCarbon_mo_t !< Temporary variable used to produce the bulk soil litter carbon quantity for output at the tile level \f$[kg C m^{-2}]\f$
    real, dimension(nmtest) :: bulkSoilCarbon_mo_t   !< Temporary variable used to produce the bulk soil carbon quantity for output at the tile level  \f$[kg C m^{-2}]\f$
    real, dimension(nmtest) :: bulkLitterResp_mo_t   !< Temporary variable used to produce the bulk litter respiration quantity for output at the tile level \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(nmtest) :: bulkSoilResp_mo_t     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output at the tile level \f$[kg C m^{-2} s^{-1}]\f$


    ! Associate names with variables defined in derived types.

    associate( &
    dofire                => c_switch%dofire,                      & !< logical:
    prescribedFire        => c_switch%prescribedFire,         & !< logical:
    lnduseon              => c_switch%lnduseon,                    & !< logical:
    Ncycle_on             => c_switch%Ncycle_on,                   & !< logical:
    PFTCompetition        => c_switch%PFTCompetition,              & !< logical:
    doperpftoutput        => c_switch%doperpftoutput,              & !< logical:
    dopertileoutput       => c_switch%dopertileoutput,             & !< logical:
    transientOBSWETF      => c_switch%transientOBSWETF,            & !< logical:
    fixedYearOBSWETF      => c_switch%fixedYearOBSWETF,            & !< integer:
    useTracer             => c_switch%useTracer,                   & !< integer: Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
                                                                     !!          useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
                                                                     !!                        the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                                                     !!          useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                                                     !!          useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    dynamicTilingOn         => c_switch%dynamicTilingOn,               & !< logical:
    trackTileAge         => c_switch%trackTileAge,               & !< logical:
    timberHarvest         => c_switch%timberHarvest,               & !< logical:

    FAREROT => class_rot%FAREROT,                                  & !< real, dimension(:,:) : Fractional coverage of mosaic tile on modelled area

    pftexistrow           => vrot%pftexist,                        & !< logical, dimension(:,:,:) :
    ccrow           => vrot%cc,                        & !< real, dimension(:,:,:) :
    mmrow           => vrot%mm,                        & !< real, dimension(:,:,:) :
    fcancmxrow            => vrot%fcancmx,                         & !< real, dimension(:,:,:) :
    laimaxg_mo            =>ctem_mo%laimaxg_mo,                    & !< real, dimension(:,:,:) :
    stemmass_mo           =>ctem_mo%stemmass_mo,                   & !< real, dimension(:,:,:) :
    stemmass_NS_mo         =>ctem_mo%stemmass_NS_mo,               & !< real, dimension(:,:,:) :
    stemmasss_mo          =>ctem_mo%stemmasss_mo,                  & !< real, dimension(:,:,:) :
    gleafmas_mo           =>ctem_mo%gleafmas_mo,                   & !< real, dimension(:,:,:) :
    gleafmas_NS_mo         =>ctem_mo%gleafmas_NS_mo,               & !< real, dimension(:,:,:) :
    gleafmass_mo          =>ctem_mo%gleafmass_mo,                  & !< real, dimension(:,:,:) :
    rootmass_mo           =>ctem_mo%rootmass_mo,                   & !< real, dimension(:,:,:) :
    rootmass_NS_mo         =>ctem_mo%rootmass_NS_mo,               & !< real, dimension(:,:,:) :
    rootmasss_mo          =>ctem_mo%rootmasss_mo,                  & !< real, dimension(:,:,:) :
    rootdpth_mo           =>ctem_mo%rootdpth_mo,                   & !< real, dimension(:,:,:) :
    litrfallveg_mo        =>ctem_mo%litrfallveg_mo,                & !< real, dimension(:,:,:) :
    humiftrsveg_mo        =>ctem_mo%humiftrsveg_mo,                & !< real, dimension(:,:,:) :
    tltrleaf_mo           =>ctem_mo%tltrleaf_mo,                   & !< real, dimension(:,:,:) :
    tltrstem_mo           =>ctem_mo%tltrstem_mo,                   & !< real, dimension(:,:,:) :
    tltrroot_mo           =>ctem_mo%tltrroot_mo,                   & !< real, dimension(:,:,:) :
    cc_mo                 =>ctem_mo%cc_mo,                        & !< real, dimension(:,:,:) :
    mm_mo                 =>ctem_mo%mm_mo,                        & !< real, dimension(:,:,:) :
    npp_mo                =>ctem_mo%npp_mo,                        & !< real, dimension(:,:,:) :
    gpp_mo                =>ctem_mo%gpp_mo,                        & !< real, dimension(:,:,:) :
    vcmax0_mo             =>ctem_mo%vcmax0_mo,                     & !< real, dimension(:,:,:) :
    leafns2s_mo           =>ctem_mo%leafns2s_mo,                   & !< real, dimension(:,:,:) :
    stemns2s_mo           =>ctem_mo%stemns2s_mo,                   & !< real, dimension(:,:,:) :
    rootns2s_mo           =>ctem_mo%rootns2s_mo,                   & !< real, dimension(:,:,:) :
    re_alloc_s2l_mo       =>ctem_mo%re_alloc_s2l_mo,               & !< real, dimension(:,:,:) :
    re_alloc_r2l_mo       =>ctem_mo%re_alloc_r2l_mo,               & !< real, dimension(:,:,:) :
    re_alloc_sr2l_mo      =>ctem_mo%re_alloc_sr2l_mo,              & !< real, dimension(:,:,:) :
    vgbiomas_mo           =>ctem_mo%vgbiomas_mo,                   & !< real, dimension(:,:,:) :
    autores_mo            =>ctem_mo%autores_mo,                    & !< real, dimension(:,:,:) :
    soilres_mo            =>ctem_mo%soilres_mo,                    & !< real, dimension(:,:,:) :
    totcmass_mo           =>ctem_mo%totcmass_mo,                   & !< real, dimension(:,:,:) :
    litrmass_mo           =>ctem_mo%litrmass_mo,                   & !< real, dimension(:,:,:,:) :
    soilcmas_mo           =>ctem_mo%soilcmas_mo,                   & !< real, dimension(:,:,:,:) :
    nh4_mass_mo           =>ctem_mo%nh4_mass_mo,                   & !< real, dimension(:,:,:) :
    no3_mass_mo           =>ctem_mo%no3_mass_mo,                   & !< real, dimension(:,:,:) :
    nlitrmass_mo          =>ctem_mo%nlitrmass_mo,                  & !< real, dimension(:,:,:) :
    soilnmas_mo           =>ctem_mo%soilnmas_mo,                   & !< real, dimension(:,:,:) :
    ngleafmas_mo          =>ctem_mo%ngleafmas_mo,                  & !< real, dimension(:,:,:) :
    ngleafmas_NS_mo        =>ctem_mo%ngleafmas_NS_mo,              & !< real, dimension(:,:,:) :
    ngleafmass_mo         =>ctem_mo%ngleafmass_mo,                 & !< real, dimension(:,:,:) :
    nbleafmas_mo          =>ctem_mo%nbleafmas_mo,                  & !< real, dimension(:,:,:) :
    nstemmass_mo          =>ctem_mo%nstemmass_mo,                  & !< real, dimension(:,:,:) :
    nstemmass_NS_mo        =>ctem_mo%nstemmass_NS_mo,              & !< real, dimension(:,:,:) :
    nstemmasss_mo         =>ctem_mo%nstemmasss_mo,                 & !< real, dimension(:,:,:) :
    nrootmass_mo          =>ctem_mo%nrootmass_mo,                  & !< real, dimension(:,:,:) :
    nrootmass_NS_mo        =>ctem_mo%nrootmass_NS_mo,              & !< real, dimension(:,:,:) :
    nrootmasss_mo         =>ctem_mo%nrootmasss_mo,                 & !< real, dimension(:,:,:) :
    c2n_l_mo              =>ctem_mo%c2n_l_mo,                      & !< real, dimension(:,:,:) :
    c2n_s_mo              =>ctem_mo%c2n_s_mo,                      & !< real, dimension(:,:,:) :
    c2n_r_mo              =>ctem_mo%c2n_r_mo,                      & !< real, dimension(:,:,:) :
    c2n_wp_mo             =>ctem_mo%c2n_wp_mo,                     & !< real, dimension(:,:,:) :
    c2n_litr_mo           =>ctem_mo%c2n_litr_mo,                   & !< real, dimension(:,:,:) :
    c2n_humus_mo          =>ctem_mo%c2n_humus_mo,                  & !< real, dimension(:,:,:) :
    nep_mo                =>ctem_mo%nep_mo,                        & !< real, dimension(:,:,:) :
    litres_mo             =>ctem_mo%litres_mo,                     & !< real, dimension(:,:,:,:) :
    soilcres_mo           =>ctem_mo%soilcres_mo,                   & !< real, dimension(:,:,:,:) :
    hetrores_mo           =>ctem_mo%hetrores_mo,                   & !< real, dimension(:,:,:) :
    nbp_mo                =>ctem_mo%nbp_mo,                        & !< real, dimension(:,:,:) :
    emit_co2_mo           =>ctem_mo%emit_co2_mo,                   & !< real, dimension(:,:,:) :
    emit_co_mo            =>ctem_mo%emit_co_mo,                    & !< real, dimension(:,:,:) :
    emit_ch4_mo           =>ctem_mo%emit_ch4_mo,                   & !< real, dimension(:,:,:) :
    emit_nmhc_mo          =>ctem_mo%emit_nmhc_mo,                  & !< real, dimension(:,:,:) :
    emit_h2_mo            =>ctem_mo%emit_h2_mo,                    & !< real, dimension(:,:,:) :
    emit_nox_mo           =>ctem_mo%emit_nox_mo,                   & !< real, dimension(:,:,:) :
    emit_n2o_mo           =>ctem_mo%emit_n2o_mo,                   & !< real, dimension(:,:,:) :
    emit_nh3_mo           =>ctem_mo%emit_nh3_mo,                   & !< real, dimension(:,:,:) :
    emit_pm25_mo          =>ctem_mo%emit_pm25_mo,                  & !< real, dimension(:,:,:) :
    emit_tpm_mo           =>ctem_mo%emit_tpm_mo,                   & !< real, dimension(:,:,:) :
    emit_tc_mo            =>ctem_mo%emit_tc_mo,                    & !< real, dimension(:,:,:) :
    emit_oc_mo            =>ctem_mo%emit_oc_mo,                    & !< real, dimension(:,:,:) :
    emit_bc_mo            =>ctem_mo%emit_bc_mo,                    & !< real, dimension(:,:,:) :
    bterm_mo              =>ctem_mo%bterm_mo,                      & !< real, dimension(:,:,:) :
    mterm_mo              =>ctem_mo%mterm_mo,                      & !< real, dimension(:,:,:) :
    burnfrac_mo           =>ctem_mo%burnfrac_mo,                   & !< real, dimension(:,:,:) :
    smfuncveg_mo          =>ctem_mo%smfuncveg_mo,                  & !< real, dimension(:,:,:) :

    bnf_tot_mo            =>ctem_mo%bnf_tot_mo,                    & !< real, dimension(:,:,:) :
    bnf_free_mo           =>ctem_mo%bnf_free_mo,                   & !< real, dimension(:,:,:) :
    bnf_ant_mo            =>ctem_mo%bnf_ant_mo,                    & !< real, dimension(:,:,:) :
    bnf_nat_mo            =>ctem_mo%bnf_nat_mo,                    & !< real, dimension(:,:,:) :
    nitrif_mo             =>ctem_mo%nitrif_mo,                     & !< real, dimension(:,:,:) :
    no_nit_mo             =>ctem_mo%no_nit_mo,                     & !< real, dimension(:,:,:) :
    no_denit_mo           =>ctem_mo%no_denit_mo,                   & !< real, dimension(:,:,:) :
    no_nitdenit_mo        =>ctem_mo%no_nitdenit_mo,                & !< real, dimension(:,:,:) :
    n2o_nit_mo            =>ctem_mo%n2o_nit_mo,                    & !< real, dimension(:,:,:) :
    n2o_denit_mo          =>ctem_mo%n2o_denit_mo,                  & !< real, dimension(:,:,:) :
    n2o_nitdenit_mo       =>ctem_mo%n2o_nitdenit_mo,               & !< real, dimension(:,:,:) :
    n2_denit_mo           =>ctem_mo%n2_denit_mo,                   & !< real, dimension(:,:,:) :
    nvol_mo               =>ctem_mo%nvol_mo,                       & !< real, dimension(:,:,:) :
    nleach_mo             =>ctem_mo%nleach_mo,                     & !< real, dimension(:,:,:) :
    appl_fert_mo          =>ctem_mo%appl_fert_mo,                  & !< real, dimension(:,:,:) :
    ndep_nh4_mo           =>ctem_mo%ndep_nh4_mo,                   & !< real, dimension(:,:,:) :
    ndep_no3_mo           =>ctem_mo%ndep_no3_mo,                   & !< real, dimension(:,:,:) :
    ndemand_wp_npp_mo     =>ctem_mo%ndemand_wp_npp_mo,             & !< real, dimension(:,:,:) :
    nuptake_p_nh4_mo      =>ctem_mo%nuptake_p_nh4_mo,              & !< real, dimension(:,:,:) :
    nuptake_p_no3_mo      =>ctem_mo%nuptake_p_no3_mo,              & !< real, dimension(:,:,:) :
    nuptake_a_actl_nh4_mo =>ctem_mo%nuptake_a_actl_nh4_mo,         & !< real, dimension(:,:,:) :
    nuptake_a_actl_no3_mo =>ctem_mo%nuptake_a_actl_no3_mo,         & !< real, dimension(:,:,:) :
    nuptake_mo            =>ctem_mo%nuptake_mo,                    & !< real, dimension(:,:,:) :
    nalloc_l_mo           =>ctem_mo%nalloc_l_mo,                   & !< real, dimension(:,:,:) :
    nalloc_s_mo           =>ctem_mo%nalloc_s_mo,                   & !< real, dimension(:,:,:) :
    nalloc_r_mo           =>ctem_mo%nalloc_r_mo,                   & !< real, dimension(:,:,:) :
    nresorped_s_mo        =>ctem_mo%nresorped_s_mo,                & !< real, dimension(:,:,:) :
    nresorped_r_mo        =>ctem_mo%nresorped_r_mo,                & !< real, dimension(:,:,:) :
    nre_alloc_s2l_mo      =>ctem_mo%nre_alloc_s2l_mo,              & !< real, dimension(:,:,:) :
    nre_alloc_r2l_mo      =>ctem_mo%nre_alloc_r2l_mo,              & !< real, dimension(:,:,:) :
    nleafns2s_mo          =>ctem_mo%nleafns2s_mo,                  & !< real, dimension(:,:,:) :
    nstemns2s_mo          =>ctem_mo%nstemns2s_mo,                  & !< real, dimension(:,:,:) :
    nrootns2s_mo          =>ctem_mo%nrootns2s_mo,                  & !< real, dimension(:,:,:) :
    nlitr_l_mo            =>ctem_mo%nlitr_l_mo,                    & !< real, dimension(:,:,:) :
    nlitr_s_mo            =>ctem_mo%nlitr_s_mo,                    & !< real, dimension(:,:,:) :
    nlitr_r_mo            =>ctem_mo%nlitr_r_mo,                    & !< real, dimension(:,:,:) :
    nlitr_mo              =>ctem_mo%nlitr_mo,                      & !< real, dimension(:,:,:) :
    gl2bl_grass_nflux_mo  =>ctem_mo%gl2bl_grass_nflux_mo,          & !< real, dimension(:,:,:) :
    nhumtrs_mo            =>ctem_mo%nhumtrs_mo,                    & !< real, dimension(:,:,:) :
    nmineral_litr_mo      =>ctem_mo%nmineral_litr_mo,              & !< real, dimension(:,:,:) :
    nmineral_humus_mo     =>ctem_mo%nmineral_humus_mo,             & !< real, dimension(:,:,:) :
    netnmineral_mo        =>ctem_mo%netnmineral_mo,                & !< real, dimension(:,:,:) :
    nimmobil_nh4_mo       =>ctem_mo%nimmobil_nh4_mo,               & !< real, dimension(:,:,:) :
    nimmobil_no3_mo       =>ctem_mo%nimmobil_no3_mo,               & !< real, dimension(:,:,:) :
    nvgbiomas_mo          =>ctem_mo%nvgbiomas_mo,                  & !< real, dimension(:,:,:) :
    nstress_mo            =>ctem_mo%nstress_mo,                    & !< real, dimension(:,:,:) :
    fNnetland_mo          =>ctem_mo%fNnetland_mo,                  & !< real, dimension(:,:,:) :

    laimaxg_mo_t          =>ctem_tile_mo%laimaxg_mo_t,             & !< real, dimension(:,:) :
    stemmass_mo_t         =>ctem_tile_mo%stemmass_mo_t,            & !< real, dimension(:,:) :
    stemmass_NS_mo_t       =>ctem_tile_mo%stemmass_NS_mo_t,        & !< real, dimension(:,:) :
    stemmasss_mo_t        =>ctem_tile_mo%stemmasss_mo_t,           & !< real, dimension(:,:) :
    gleafmas_mo_t         =>ctem_tile_mo%gleafmas_mo_t,            & !< real, dimension(:,:) :
    gleafmas_NS_mo_t       =>ctem_tile_mo%gleafmas_NS_mo_t,        & !< real, dimension(:,:) :
    gleafmass_mo_t        =>ctem_tile_mo%gleafmass_mo_t,           & !< real, dimension(:,:) :
    rootmass_mo_t         =>ctem_tile_mo%rootmass_mo_t,            & !< real, dimension(:,:) :
    rootmass_NS_mo_t       =>ctem_tile_mo%rootmass_NS_mo_t,        & !< real, dimension(:,:) :
    rootmasss_mo_t        =>ctem_tile_mo%rootmasss_mo_t,           & !< real, dimension(:,:) :
    litrfall_mo_t         =>ctem_tile_mo%litrfall_mo_t,            & !< real, dimension(:,:) :
    humiftrs_mo_t         =>ctem_tile_mo%humiftrs_mo_t,            & !< real, dimension(:,:) :
    tltrleaf_mo_t         =>ctem_tile_mo%tltrleaf_mo_t,            & !< real, dimension(:,:) :
    tltrstem_mo_t         =>ctem_tile_mo%tltrstem_mo_t,            & !< real, dimension(:,:) :
    tltrroot_mo_t         =>ctem_tile_mo%tltrroot_mo_t,            & !< real, dimension(:,:) :

    npp_mo_t              =>ctem_tile_mo%npp_mo_t,                 & !< real, dimension(:,:) :
    gpp_mo_t              =>ctem_tile_mo%gpp_mo_t,                 & !< real, dimension(:,:) :
    vcmax0_mo_t           =>ctem_tile_mo%vcmax0_mo_t,              & !< real, dimension(:,:) :
    leafns2s_mo_t         =>ctem_tile_mo%leafns2s_mo_t,            & !< real, dimension(:,:) :
    stemns2s_mo_t         =>ctem_tile_mo%stemns2s_mo_t,            & !< real, dimension(:,:) :
    rootns2s_mo_t         =>ctem_tile_mo%rootns2s_mo_t,            & !< real, dimension(:,:) :
    re_alloc_s2l_mo_t     =>ctem_tile_mo%re_alloc_s2l_mo_t,        & !< real, dimension(:,:) :
    re_alloc_r2l_mo_t     =>ctem_tile_mo%re_alloc_r2l_mo_t,        & !< real, dimension(:,:) :
    re_alloc_sr2l_mo_t    =>ctem_tile_mo%re_alloc_sr2l_mo_t,       & !< real, dimension(:,:) :
    vgbiomas_mo_t         =>ctem_tile_mo%vgbiomas_mo_t,            & !< real, dimension(:,:) :
    autores_mo_t          =>ctem_tile_mo%autores_mo_t,             & !< real, dimension(:,:) :
    soilres_mo_t          =>ctem_tile_mo%soilres_mo_t,             & !< real, dimension(:,:) :
    totcmass_mo_t         =>ctem_tile_mo%totcmass_mo_t,            & !< real, dimension(:,:) :
    litrmass_mo_t         =>ctem_tile_mo%litrmass_mo_t,            & !< real, dimension(:,:,:) :
    soilcmas_mo_t         =>ctem_tile_mo%soilcmas_mo_t,            & !< real, dimension(:,:,:) :
    nh4_mass_mo_t         =>ctem_tile_mo%nh4_mass_mo_t,            & !< real, dimension(:,:) :
    no3_mass_mo_t         =>ctem_tile_mo%no3_mass_mo_t,            & !< real, dimension(:,:) :
    ngleafmas_mo_t        =>ctem_tile_mo%ngleafmas_mo_t,           & !< real, dimension(:,:) :
    ngleafmas_NS_mo_t      =>ctem_tile_mo%ngleafmas_NS_mo_t,       & !< real, dimension(:,:) :
    ngleafmass_mo_t       =>ctem_tile_mo%ngleafmass_mo_t,          & !< real, dimension(:,:) :
    nbleafmas_mo_t        =>ctem_tile_mo%nbleafmas_mo_t,           & !< real, dimension(:,:) :
    nstemmass_mo_t        =>ctem_tile_mo%nstemmass_mo_t,           & !< real, dimension(:,:) :
    nstemmass_NS_mo_t      =>ctem_tile_mo%nstemmass_NS_mo_t,       & !< real, dimension(:,:) :
    nstemmasss_mo_t       =>ctem_tile_mo%nstemmasss_mo_t,          & !< real, dimension(:,:) :
    nrootmass_mo_t        =>ctem_tile_mo%nrootmass_mo_t,           & !< real, dimension(:,:) :
    nrootmass_NS_mo_t      =>ctem_tile_mo%nrootmass_NS_mo_t,       & !< real, dimension(:,:) :
    nrootmasss_mo_t       =>ctem_tile_mo%nrootmasss_mo_t,          & !< real, dimension(:,:) :
    c2n_l_mo_t            =>ctem_tile_mo%c2n_l_mo_t,               & !< real, dimension(:,:) :
    c2n_s_mo_t            =>ctem_tile_mo%c2n_s_mo_t,               & !< real, dimension(:,:) :
    c2n_r_mo_t            =>ctem_tile_mo%c2n_r_mo_t,               & !< real, dimension(:,:) :
    c2n_wp_mo_t           =>ctem_tile_mo%c2n_wp_mo_t,              & !< real, dimension(:,:) :
    c2n_litr_mo_t         =>ctem_tile_mo%c2n_litr_mo_t,            & !< real, dimension(:,:) :
    c2n_humus_mo_t        =>ctem_tile_mo%c2n_humus_mo_t,           & !< real, dimension(:,:) :
    nlitrmass_mo_t        =>ctem_tile_mo%nlitrmass_mo_t,           & !< real, dimension(:,:) :
    soilnmas_mo_t         =>ctem_tile_mo%soilnmas_mo_t,            & !< real, dimension(:,:) :
    nep_mo_t              =>ctem_tile_mo%nep_mo_t,                 & !< real, dimension(:,:) :
    litres_mo_t           =>ctem_tile_mo%litres_mo_t,              & !< real, dimension(:,:,:) :
    soilcres_mo_t         =>ctem_tile_mo%soilcres_mo_t,            & !< real, dimension(:,:,:) :
    hetrores_mo_t         =>ctem_tile_mo%hetrores_mo_t,            & !< real, dimension(:,:) :
    nbp_mo_t              =>ctem_tile_mo%nbp_mo_t,                 & !< real, dimension(:,:) :
    emit_co2_mo_t         =>ctem_tile_mo%emit_co2_mo_t,            & !< real, dimension(:,:) :
    emit_co_mo_t          =>ctem_tile_mo%emit_co_mo_t,             & !< real, dimension(:,:) :
    emit_ch4_mo_t         =>ctem_tile_mo%emit_ch4_mo_t,            & !< real, dimension(:,:) :
    emit_nmhc_mo_t        =>ctem_tile_mo%emit_nmhc_mo_t,           & !< real, dimension(:,:) :
    emit_h2_mo_t          =>ctem_tile_mo%emit_h2_mo_t,             & !< real, dimension(:,:) :
    emit_nox_mo_t         =>ctem_tile_mo%emit_nox_mo_t,            & !< real, dimension(:,:) :
    emit_n2o_mo_t         =>ctem_tile_mo%emit_n2o_mo_t,            & !< real, dimension(:,:) :
    emit_nh3_mo_t         =>ctem_tile_mo%emit_nh3_mo_t,            & !< real, dimension(:,:) :
    emit_pm25_mo_t        =>ctem_tile_mo%emit_pm25_mo_t,           & !< real, dimension(:,:) :
    emit_tpm_mo_t         =>ctem_tile_mo%emit_tpm_mo_t,            & !< real, dimension(:,:) :
    emit_tc_mo_t          =>ctem_tile_mo%emit_tc_mo_t,             & !< real, dimension(:,:) :
    emit_oc_mo_t          =>ctem_tile_mo%emit_oc_mo_t,             & !< real, dimension(:,:) :
    emit_bc_mo_t          =>ctem_tile_mo%emit_bc_mo_t,             & !< real, dimension(:,:) :
    burnfrac_mo_t         =>ctem_tile_mo%burnfrac_mo_t,            & !< real, dimension(:,:) :
    smfuncveg_mo_t        =>ctem_tile_mo%smfuncveg_mo_t,           & !< real, dimension(:,:) :
    bterm_mo_t            =>ctem_tile_mo%bterm_mo_t,               & !< real, dimension(:,:) :
    luc_emc_mo_t          =>ctem_tile_mo%luc_emc_mo_t,             & !< real, dimension(:,:) :
    lterm_mo_t            =>ctem_tile_mo%lterm_mo_t,               & !< real, dimension(:,:) :
    lucsocin_mo_t         =>ctem_tile_mo%lucsocin_mo_t,            & !< real, dimension(:,:) :
    mterm_mo_t            =>ctem_tile_mo%mterm_mo_t,               & !< real, dimension(:,:) :
    lucltrin_mo_t         =>ctem_tile_mo%lucltrin_mo_t,            & !< real, dimension(:,:) :
    luc_emcn_mo_t          =>ctem_tile_mo%luc_emcn_mo_t,             & !< real, dimension(:,:) :
    lucsocinn_mo_t         =>ctem_tile_mo%lucsocinn_mo_t,            & !< real, dimension(:,:) :
    lucltrinn_mo_t         =>ctem_tile_mo%lucltrinn_mo_t,            & !< real, dimension(:,:) :
    ch4WetSpec_mo_t       =>ctem_tile_mo%ch4WetSpec_mo_t,          & !< real, dimension(:,:) :
    wetfdyn_mo_t          =>ctem_tile_mo%wetfdyn_mo_t,             & !< real, dimension(:,:) :
    wetfpres_mo_t         =>ctem_tile_mo%wetfpres_mo_t,            & !< real, dimension(:,:) :
    ch4WetDyn_mo_t        =>ctem_tile_mo%ch4WetDyn_mo_t,           & !< real, dimension(:,:) :
    ch4soills_mo_t        =>ctem_tile_mo%ch4soills_mo_t,           & !< real, dimension(:,:) :
    wind_mo_t             =>ctem_tile_mo%wind_mo_t,                & !< real, dimension(:,:) :
    fProductDecomp_mo_t   =>ctem_tile_mo%fProductDecomp_mo_t,      & !< real, dimension(:,:) : Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$
    tileAge_mo_t          =>ctem_tile_mo%tileAge_mo_t,             & !< real, dimension(:,:) :

    bnf_tot_mo_t          =>ctem_tile_mo%bnf_tot_mo_t,             & !< real, dimension(:,:) :
    bnf_free_mo_t         =>ctem_tile_mo%bnf_free_mo_t,            & !< real, dimension(:,:) :
    bnf_ant_mo_t          =>ctem_tile_mo%bnf_ant_mo_t,             & !< real, dimension(:,:) :
    bnf_nat_mo_t          =>ctem_tile_mo%bnf_nat_mo_t,             & !< real, dimension(:,:) :
    nitrif_mo_t           =>ctem_tile_mo%nitrif_mo_t,              & !< real, dimension(:,:) :
    no_nit_mo_t           =>ctem_tile_mo%no_nit_mo_t,              & !< real, dimension(:,:) :
    no_denit_mo_t         =>ctem_tile_mo%no_denit_mo_t,            & !< real, dimension(:,:) :
    no_nitdenit_mo_t      =>ctem_tile_mo%no_nitdenit_mo_t,         & !< real, dimension(:,:) :
    n2o_nit_mo_t          =>ctem_tile_mo%n2o_nit_mo_t,             & !< real, dimension(:,:) :
    n2o_denit_mo_t        =>ctem_tile_mo%n2o_denit_mo_t,           & !< real, dimension(:,:) :
    n2o_nitdenit_mo_t     =>ctem_tile_mo%n2o_nitdenit_mo_t,        & !< real, dimension(:,:) :
    n2_denit_mo_t         =>ctem_tile_mo%n2_denit_mo_t,            & !< real, dimension(:,:) :
    nvol_mo_t             =>ctem_tile_mo%nvol_mo_t,                & !< real, dimension(:,:) :
    nleach_mo_t           =>ctem_tile_mo%nleach_mo_t,              & !< real, dimension(:,:) :
    appl_fert_mo_t        =>ctem_tile_mo%appl_fert_mo_t,           & !< real, dimension(:,:) :
    ndep_nh4_mo_t         =>ctem_tile_mo%ndep_nh4_mo_t,            & !< real, dimension(:,:) :
    ndep_no3_mo_t         =>ctem_tile_mo%ndep_no3_mo_t,            & !< real, dimension(:,:) :
    ndemand_wp_npp_mo_t   =>ctem_tile_mo%ndemand_wp_npp_mo_t,      & !< real, dimension(:,:) :
    nuptake_p_nh4_mo_t    =>ctem_tile_mo%nuptake_p_nh4_mo_t,       & !< real, dimension(:,:) :
    nuptake_p_no3_mo_t    =>ctem_tile_mo%nuptake_p_no3_mo_t,       & !< real, dimension(:,:) :
    nuptake_a_actl_nh4_mo_t=>ctem_tile_mo%nuptake_a_actl_nh4_mo_t, & !< real, dimension(:,:) :
    nuptake_a_actl_no3_mo_t=>ctem_tile_mo%nuptake_a_actl_no3_mo_t, & !< real, dimension(:,:) :
    nuptake_mo_t          =>ctem_tile_mo%nuptake_mo_t,             & !< real, dimension(:,:) :
    nalloc_l_mo_t         =>ctem_tile_mo%nalloc_l_mo_t,            & !< real, dimension(:,:) :
    nalloc_s_mo_t         =>ctem_tile_mo%nalloc_s_mo_t,            & !< real, dimension(:,:) :
    nalloc_r_mo_t         =>ctem_tile_mo%nalloc_r_mo_t,            & !< real, dimension(:,:) :
    nresorped_s_mo_t      =>ctem_tile_mo%nresorped_s_mo_t,         & !< real, dimension(:,:) :
    nresorped_r_mo_t      =>ctem_tile_mo%nresorped_r_mo_t,         & !< real, dimension(:,:) :
    nre_alloc_s2l_mo_t    =>ctem_tile_mo%nre_alloc_s2l_mo_t,       & !< real, dimension(:,:) :
    nre_alloc_r2l_mo_t    =>ctem_tile_mo%nre_alloc_r2l_mo_t,       & !< real, dimension(:,:) :
    nleafns2s_mo_t        =>ctem_tile_mo%nleafns2s_mo_t,           & !< real, dimension(:,:) :
    nstemns2s_mo_t        =>ctem_tile_mo%nstemns2s_mo_t,           & !< real, dimension(:,:) :
    nrootns2s_mo_t        =>ctem_tile_mo%nrootns2s_mo_t,           & !< real, dimension(:,:) :
    nlitr_l_mo_t          =>ctem_tile_mo%nlitr_l_mo_t,             & !< real, dimension(:,:) :
    nlitr_s_mo_t          =>ctem_tile_mo%nlitr_s_mo_t,             & !< real, dimension(:,:) :
    nlitr_r_mo_t          =>ctem_tile_mo%nlitr_r_mo_t,             & !< real, dimension(:,:) :
    nlitr_mo_t            =>ctem_tile_mo%nlitr_mo_t,               & !< real, dimension(:,:) :
    gl2bl_grass_nflux_mo_t=>ctem_tile_mo%gl2bl_grass_nflux_mo_t,   & !< real, dimension(:,:) :
    nhumtrs_mo_t          =>ctem_tile_mo%nhumtrs_mo_t,             & !< real, dimension(:,:) :
    nmineral_litr_mo_t    =>ctem_tile_mo%nmineral_litr_mo_t,       & !< real, dimension(:,:) :
    nmineral_humus_mo_t   =>ctem_tile_mo%nmineral_humus_mo_t,      & !< real, dimension(:,:) :
    netnmineral_mo_t      =>ctem_tile_mo%netnmineral_mo_t,         & !< real, dimension(:,:) :
    nimmobil_nh4_mo_t     =>ctem_tile_mo%nimmobil_nh4_mo_t,        & !< real, dimension(:,:) :
    nimmobil_no3_mo_t     =>ctem_tile_mo%nimmobil_no3_mo_t,        & !< real, dimension(:,:) :
    fNnetland_mo_t        =>ctem_tile_mo%fNnetland_mo_t,           & !< real, dimension(:,:) :
    nvgbiomas_mo_t        =>ctem_tile_mo%nvgbiomas_mo_t,           & !< real, dimension(:,:) :
    nstress_mo_t          =>ctem_tile_mo%nstress_mo_t,             & !< real, dimension(:,:) :

    gppvegrow         => vrot%gppveg,                              & !< real, dimension(:,:,:) :
    vcmax0row         => vrot%vcmax0,                              & !< real, dimension(:,:,:) :
    leafns2srow       => vrot%leafns2s,                            & !< real, dimension(:,:,:) :
    stemns2srow       => vrot%stemns2s,                            & !< real, dimension(:,:,:) :
    rootns2srow       => vrot%rootns2s,                            & !< real, dimension(:,:,:) :
    re_alloc_s2lrow   => vrot%re_alloc_s2l,                        & !< real, dimension(:,:,:) :
    re_alloc_r2lrow   => vrot%re_alloc_r2l,                        & !< real, dimension(:,:,:) :
    re_alloc_sr2lrow  => vrot%re_alloc_sr2l,                       & !< real, dimension(:,:,:) :
    nepvegrow         => vrot%nepveg,                              & !< real, dimension(:,:,:) :
    nbpvegrow         => vrot%nbpveg,                              & !< real, dimension(:,:,:) :
    nbprow            => vrot%nbp,                                 & !< real, dimension(:,:)   :
    nppvegrow         => vrot%nppveg,                              & !< real, dimension(:,:,:) :
    hetroresvegrow    => vrot%hetroresveg,                         & !< real, dimension(:,:,:) :
    autoresvegrow     => vrot%autoresveg,                          & !< real, dimension(:,:,:) :
    litresvegrow      => vrot%litresveg,                           & !< real, dimension(:,:,:,:) :
    soilcresvegrow    => vrot%soilcresveg,                         & !< real, dimension(:,:,:,:) :
    rmlvegaccrow      => vrot%rmlvegacc,                           & !< real, dimension(:,:,:) :
    rmsvegrow         => vrot%rmsveg,                              & !< real, dimension(:,:,:) :
    rmrvegrow         => vrot%rmrveg,                              & !< real, dimension(:,:,:) :
    rgvegrow          => vrot%rgveg,                               & !< real, dimension(:,:,:) :
    afrrootrow        => vrot%afrroot,                             & !< real, dimension(:,:,:) :
    ailcgrow          => vrot%ailcg,                               & !< real, dimension(:,:,:) :
    emit_co2row       => vrot%emit_co2,                            & !< real, dimension(:,:,:) :
    emit_corow        => vrot%emit_co,                             & !< real, dimension(:,:,:) :
    emit_ch4row       => vrot%emit_ch4,                            & !< real, dimension(:,:,:) :
    emit_nmhcrow      => vrot%emit_nmhc,                           & !< real, dimension(:,:,:) :
    emit_h2row        => vrot%emit_h2,                             & !< real, dimension(:,:,:) :
    emit_noxrow       => vrot%emit_nox,                            & !< real, dimension(:,:,:) :
    emit_n2orow       => vrot%emit_n2o,                            & !< real, dimension(:,:,:) :
    emit_nh3row       => vrot%emit_nh3,                            & !< real, dimension(:,:,:) :
    emit_pm25row      => vrot%emit_pm25,                           & !< real, dimension(:,:,:) :
    emit_tpmrow       => vrot%emit_tpm,                            & !< real, dimension(:,:,:) :
    emit_tcrow        => vrot%emit_tc,                             & !< real, dimension(:,:,:) :
    emit_ocrow        => vrot%emit_oc,                             & !< real, dimension(:,:,:) :
    emit_bcrow        => vrot%emit_bc,                             & !< real, dimension(:,:,:) :
    burnfracrow       => vrot%burnfrac,                            & !< real, dimension(:,:) :
    burnvegfrow       => vrot%burnvegf,                            & !< real, dimension(:,:,:) :
    smfuncvegrow      => vrot%smfuncveg,                           & !< real, dimension(:,:,:) :
    btermrow          => vrot%bterm,                               & !< real, dimension(:,:,:) :
    ltermrow          => vrot%lterm,                               & !< real, dimension(:,:) :
    mtermrow          => vrot%mterm,                               & !< real, dimension(:,:,:) :
    lucemcomrow       => vrot%lucemcom,                            & !< real, dimension(:,:) :
    lucltrinrow       => vrot%lucltrin,                            & !< real, dimension(:,:) :
    lucsocinrow       => vrot%lucsocin,                            & !< real, dimension(:,:) :
    lucemcomnrow       => vrot%lucemcomn,                            & !< real, dimension(:,:) :
    lucltrinnrow       => vrot%lucltrinn,                            & !< real, dimension(:,:) :
    lucsocinnrow       => vrot%lucsocinn,                            & !< real, dimension(:,:) :
    ch4WetSpecrow     => vrot%ch4WetSpec,                          & !< real, dimension(:,:) :
    wetfdynrow        => vrot%wetfdyn,                             & !< real, dimension(:,:) :
    wetfrac_presrow   => vrot%wetfrac_pres,                        & !< real, dimension(:,:) :

    ch4WetDynrow      => vrot%ch4WetDyn,                           & !< real, dimension(:,:) :
    ch4soillsrow      => vrot%ch4_soills,                          & !< real, dimension(:,:) :
    litrmassrow       => vrot%litrmass,                            & !< real, dimension(:,:,:,:) :
    soilcmasrow       => vrot%soilcmas,                            & !< real, dimension(:,:,:,:) :
    nh4_massrow       => vrot%nh4_mass,                            & !< real, dimension(:,:,:) :
    no3_massrow       => vrot%no3_mass,                            & !< real, dimension(:,:,:) :
    nlitrmassrow      => vrot%nlitrmass,                           & !< real, dimension(:,:,:) :
    soilnmasrow       => vrot%soilnmas,                            & !< real, dimension(:,:,:) :
    ngleafmasrow      => vrot%ngleafmas,                           & !< real, dimension(:,:,:) :
    ngleafmas_NSrow    => vrot%ngleafmas_ns,                       & !< real, dimension(:,:,:) :
    ngleafmassrow     => vrot%ngleafmas_s,                         & !< real, dimension(:,:,:) :
    nbleafmasrow      => vrot%nbleafmas,                           & !< real, dimension(:,:,:) :
    nstemmassrow      => vrot%nstemmass,                           & !< real, dimension(:,:,:) :
    nstemmass_NSrow    => vrot%nstemmass_ns,                       & !< real, dimension(:,:,:) :
    nstemmasssrow     => vrot%nstemmass_s,                         & !< real, dimension(:,:,:) :
    nrootmassrow      => vrot%nrootmass,                           & !< real, dimension(:,:,:) :
    nrootmass_NSrow    => vrot%nrootmass_ns,                       & !< real, dimension(:,:,:) :
    nrootmasssrow     => vrot%nrootmass_s,                         & !< real, dimension(:,:,:) :
    c2nveg_lrow       => vrot%c2nveg_l,                            & !< real, dimension(:,:,:) :
    c2nveg_srow       => vrot%c2nveg_s,                            & !< real, dimension(:,:,:) :
    c2nveg_rrow       => vrot%c2nveg_r,                            & !< real, dimension(:,:,:) :
    c2nveg_wprow      => vrot%c2nveg_wp,                           & !< real, dimension(:,:,:) :
    c2nveg_litrrow    => vrot%c2nveg_litr,                         & !< real, dimension(:,:,:) :
    c2nveg_humusrow   => vrot%c2nveg_humus,                        & !< real, dimension(:,:,:) :
    vgbiomas_vegrow   => vrot%vgbiomas_veg,                        & !< real, dimension(:,:,:) :
    gleafmasrow       => vrot%gleafmas,                            & !< real, dimension(:,:,:) :
    gleafmas_NSrow     => vrot%gleafmas_ns,                        & !< real, dimension(:,:,:) :
    gleafmassrow      => vrot%gleafmas_s,                          & !< real, dimension(:,:,:) :
    bleafmassrow      => vrot%bleafmas,                            & !< real, dimension(:,:,:) :
    stemmassrow       => vrot%stemmass,                            & !< real, dimension(:,:,:) :
    stemmass_NSrow     => vrot%stemmass_ns,                        & !< real, dimension(:,:,:) :
    stemmasssrow      => vrot%stemmass_s,                          & !< real, dimension(:,:,:) :
    rootmassrow       => vrot%rootmass,                            & !< real, dimension(:,:,:) :
    rootmass_NSrow     => vrot%rootmass_ns,                        & !< real, dimension(:,:,:) :
    rootmasssrow      => vrot%rootmass_s,                          & !< real, dimension(:,:,:) :
    rootdpthrow       => vrot%rootdpth,                            & !< real, dimension(:,:,:) :
    uvaccrow_m        => vrot%uvaccrow_m,                          & !< real, dimension(:,:) :
    vvaccrow_m        => vrot%vvaccrow_m,                          & !< real, dimension(:,:) :
    litrfallvegrow    => vrot%litrfallveg,                         & !< real, dimension(:,:,:) :
    humiftrsvegrow    => vrot%humiftrsveg,                         & !< real, dimension(:,:,:,:) :
    tltrleafrow       => vrot%tltrleaf,                            & !< real, dimension(:,:,:) :
    tltrstemrow       => vrot%tltrstem,                            & !< real, dimension(:,:,:) :
    tltrrootrow       => vrot%tltrroot,                            & !< real, dimension(:,:,:) :

    tracerGLeafMassrot   => tracer%gLeafMassrot,                   & !< real, dimension(:,:,:) : Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    tracerBLeafMassrot   => tracer%bLeafMassrot,                   & !< real, dimension(:,:,:) : Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    tracerStemMassrot    => tracer%stemMassrot,                    & !< real, dimension(:,:,:) : Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    tracerRootMassrot    => tracer%rootMassrot,                    & !< real, dimension(:,:,:) : Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    tracerLitrMassrot    => tracer%litrMassrot,                    & !< real, dimension(:,:,:,:) : Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    tracerSoilCMassrot   => tracer%soilCMassrot,                   & !< real, dimension(:,:,:,:) : Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    tracerMossCMassrot   => tracer%mossCMassrot,                   & !< real, dimension(:,:) : Tracer mass in moss biomass, \f$kg C/m^2\f$
    tracerMossLitrMassrot => tracer%mossLitrMassrot,               & !< real, dimension(:,:) : Tracer mass in moss litter, \f$kg C/m^2\f$

    bnftotrow        => vrot%bnf_tot,                           & !< real, dimension(:,:,:)   :
    bnffreerow       => vrot%bnf_free,                          & !< real, dimension(:,:,:)   :
    bnfantrow        => vrot%bnf_ant,                           & !< real, dimension(:,:,:)   :
    bnfnatrow        => vrot%bnf_nat,                           & !< real, dimension(:,:,:)   :
    nitrifvegrow      => vrot%nitrifveg,                           & !< real, dimension(:,:,:) :
    no_nitvegrow      => vrot%no_nitveg,                           & !< real, dimension(:,:,:) :
    no_denitvegrow    => vrot%no_denitveg,                         & !< real, dimension(:,:,:) :
    no_nitdenitvegrow => vrot%no_nitdenitveg,                      & !< real, dimension(:,:,:) :
    n2o_nitvegrow     => vrot%n2o_nitveg,                          & !< real, dimension(:,:,:) :
    n2o_denitvegrow   => vrot%n2o_denitveg,                        & !< real, dimension(:,:,:) :
    n2o_nitdenitvegrow=> vrot%n2o_nitdenitveg,                     & !< real, dimension(:,:,:) :
    n2_denitvegrow    => vrot%n2_denitveg,                         & !< real, dimension(:,:,:) :
    nvolvegrow        => vrot%nvolveg,                             & !< real, dimension(:,:,:) :
    nleachvegrow      => vrot%nleachveg,                           & !< real, dimension(:,:,:) :
    appl_fertrow      => vrot%appl_fert,                           & !< real, dimension(:,:,:) :
    ndep_nh4row       => vrot%ndep_nh4,                            & !< real, dimension(:,:,:) :
    ndep_no3row       => vrot%ndep_no3,                            & !< real, dimension(:,:,:) :
    ndemandveg_wp_npprow     => vrot%ndemandveg_wp_npp,            & !< real, dimension(:,:,:) :
    nuptakeveg_p_nh4row      => vrot%nuptakeveg_p_nh4,             & !< real, dimension(:,:,:) :
    nuptakeveg_p_no3row      => vrot%nuptakeveg_p_no3,             & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_nh4row => vrot%nuptakeveg_a_actl_nh4,        & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_no3row => vrot%nuptakeveg_a_actl_no3,        & !< real, dimension(:,:,:) :
    nuptakevegrow            => vrot%nuptakeveg,                   & !< real, dimension(:,:,:) :
    nallocveg_lrow           => vrot%nallocveg_l,                  & !< real, dimension(:,:,:) :
    nallocveg_srow           => vrot%nallocveg_s,                  & !< real, dimension(:,:,:) :
    nallocveg_rrow           => vrot%nallocveg_r,                  & !< real, dimension(:,:,:) :
    nresorpedveg_srow        => vrot%nresorpedveg_s,               & !< real, dimension(:,:,:) :
    nresorpedveg_rrow        => vrot%nresorpedveg_r,               & !< real, dimension(:,:,:) :
    nre_allocveg_s2lrow      => vrot%nre_allocveg_s2l,             & !< real, dimension(:,:,:) :
    nre_allocveg_r2lrow      => vrot%nre_allocveg_r2l,             & !< real, dimension(:,:,:) :
    nleafns2svegrow          => vrot%nleafns2sveg,                 & !< real, dimension(:,:,:) :
    nstemns2svegrow          => vrot%nstemns2sveg,                 & !< real, dimension(:,:,:) :
    nrootns2svegrow          => vrot%nrootns2sveg,                 & !< real, dimension(:,:,:) :
    nlitrveg_lrow            => vrot%nlitrveg_l,                   & !< real, dimension(:,:,:) :
    nlitrveg_srow            => vrot%nlitrveg_s,                   & !< real, dimension(:,:,:) :
    nlitrveg_rrow            => vrot%nlitrveg_r,                   & !< real, dimension(:,:,:) :
    nlitrvegrow              => vrot%nlitrveg,                     & !< real, dimension(:,:,:) :
    gl2bl_grass_nfluxrow     => vrot%gl2bl_grass_nflux,            & !< real, dimension(:,:,:) :
    nhumtrsvegrow            => vrot%nhumtrsveg,                   & !< real, dimension(:,:,:) :
    nmineralveg_litrrow      => vrot%nmineralveg_litr,             & !< real, dimension(:,:,:) :
    nmineralveg_humusrow     => vrot%nmineralveg_humus,            & !< real, dimension(:,:,:) :
    netnmineralveg_row       => vrot%netnmineralveg,               & !< real, dimension(:,:,:) :
    nimmobilveg_nh4row       => vrot%nimmobilveg_nh4,              & !< real, dimension(:,:,:) :
    nimmobilveg_no3row       => vrot%nimmobilveg_no3,              & !< real, dimension(:,:,:) :
    fNnetlandvegrow          => vrot%fNnetlandveg,                 & !< real, dimension(:,:,:) :
    nvgbiomas_vegrow         => vrot%nvgbiomas_veg,                & !< real, dimension(:,:,:) :
    nstressrow               => vrot%nstress,                      & !< real, dimension(:,:,:)   :
    tileAgerot               => vrot%tileAgerow,                   & !< real, dimension(:,:)   :

    laimaxg_mo_g        =>ctem_grd_mo%laimaxg_mo_g,                & !< real, dimension(:) :
    stemmass_mo_g       =>ctem_grd_mo%stemmass_mo_g,               & !< real, dimension(:) :
    stemmass_NS_mo_g     =>ctem_grd_mo%stemmass_NS_mo_g,           & !< real, dimension(:) :
    stemmasss_mo_g      =>ctem_grd_mo%stemmasss_mo_g,              & !< real, dimension(:) :
    gleafmas_mo_g       =>ctem_grd_mo%gleafmas_mo_g,               & !< real, dimension(:) :
    gleafmas_NS_mo_g     =>ctem_grd_mo%gleafmas_NS_mo_g,           & !< real, dimension(:) :
    gleafmass_mo_g      =>ctem_grd_mo%gleafmass_mo_g,              & !< real, dimension(:) :
    rootmass_mo_g       =>ctem_grd_mo%rootmass_mo_g,               & !< real, dimension(:) :
    rootmass_NS_mo_g     =>ctem_grd_mo%rootmass_NS_mo_g,           & !< real, dimension(:) :
    rootmasss_mo_g      =>ctem_grd_mo%rootmasss_mo_g,              & !< real, dimension(:) :
    litrmass_mo_g       =>ctem_grd_mo%litrmass_mo_g,               & !< real, dimension(:,:) :
    soilcmas_mo_g       =>ctem_grd_mo%soilcmas_mo_g,               & !< real, dimension(:,:) :
    nh4_mass_mo_g       =>ctem_grd_mo%nh4_mass_mo_g,               & !< real, dimension(:) :
    no3_mass_mo_g       =>ctem_grd_mo%no3_mass_mo_g,               & !< real, dimension(:) :
    ngleafmas_mo_g      =>ctem_grd_mo%ngleafmas_mo_g,              & !< real, dimension(:) :
    ngleafmas_NS_mo_g    =>ctem_grd_mo%ngleafmas_NS_mo_g,          & !< real, dimension(:) :
    ngleafmass_mo_g     =>ctem_grd_mo%ngleafmass_mo_g,             & !< real, dimension(:) :
    nbleafmas_mo_g      =>ctem_grd_mo%nbleafmas_mo_g,              & !< real, dimension(:) :
    nstemmass_mo_g      =>ctem_grd_mo%nstemmass_mo_g,              & !< real, dimension(:) :
    nstemmass_NS_mo_g    =>ctem_grd_mo%nstemmass_NS_mo_g,          & !< real, dimension(:) :
    nstemmasss_mo_g     =>ctem_grd_mo%nstemmasss_mo_g,             & !< real, dimension(:) :
    nrootmass_mo_g      =>ctem_grd_mo%nrootmass_mo_g,              & !< real, dimension(:) :
    nrootmass_NS_mo_g    =>ctem_grd_mo%nrootmass_NS_mo_g,          & !< real, dimension(:) :
    nrootmasss_mo_g     =>ctem_grd_mo%nrootmasss_mo_g,             & !< real, dimension(:) :
    c2n_l_mo_g          =>ctem_grd_mo%c2n_l_mo_g,                  & !< real, dimension(:) :
    c2n_s_mo_g          =>ctem_grd_mo%c2n_s_mo_g,                  & !< real, dimension(:) :
    c2n_r_mo_g          =>ctem_grd_mo%c2n_r_mo_g,                  & !< real, dimension(:) :
    c2n_wp_mo_g         =>ctem_grd_mo%c2n_wp_mo_g,                 & !< real, dimension(:) :
    c2n_litr_mo_g       =>ctem_grd_mo%c2n_litr_mo_g,               & !< real, dimension(:) :
    c2n_humus_mo_g      =>ctem_grd_mo%c2n_humus_mo_g,              & !< real, dimension(:) :
    nlitrmass_mo_g      =>ctem_grd_mo%nlitrmass_mo_g,              & !< real, dimension(:) :
    soilnmas_mo_g       =>ctem_grd_mo%soilnmas_mo_g,               & !< real, dimension(:) :
    litrfall_mo_g       =>ctem_grd_mo%litrfall_mo_g,               & !< real, dimension(:) :
    humiftrs_mo_g       =>ctem_grd_mo%humiftrs_mo_g,               & !< real, dimension(:) :
    tltrleaf_mo_g       =>ctem_grd_mo%tltrleaf_mo_g,               & !< real, dimension(:) :
    tltrstem_mo_g       =>ctem_grd_mo%tltrstem_mo_g,               & !< real, dimension(:) :
    tltrroot_mo_g       =>ctem_grd_mo%tltrroot_mo_g,               & !< real, dimension(:) :

    npp_mo_g            =>ctem_grd_mo%npp_mo_g,                    & !< real, dimension(:) :
    gpp_mo_g            =>ctem_grd_mo%gpp_mo_g,                    & !< real, dimension(:) :
    vcmax0_mo_g         =>ctem_grd_mo%vcmax0_mo_g,                 & !< real, dimension(:) :
    leafns2s_mo_g       =>ctem_grd_mo%leafns2s_mo_g,               & !< real, dimension(:) :
    stemns2s_mo_g       =>ctem_grd_mo%stemns2s_mo_g,               & !< real, dimension(:) :
    rootns2s_mo_g       =>ctem_grd_mo%rootns2s_mo_g,               & !< real, dimension(:) :
    re_alloc_s2l_mo_g   =>ctem_grd_mo%re_alloc_s2l_mo_g,           & !< real, dimension(:) :
    re_alloc_r2l_mo_g   =>ctem_grd_mo%re_alloc_r2l_mo_g,           & !< real, dimension(:) :
    re_alloc_sr2l_mo_g  =>ctem_grd_mo%re_alloc_sr2l_mo_g,          & !< real, dimension(:) :
    nep_mo_g            =>ctem_grd_mo%nep_mo_g,                    & !< real, dimension(:) :
    nbp_mo_g            =>ctem_grd_mo%nbp_mo_g,                    & !< real, dimension(:) :
    hetrores_mo_g       =>ctem_grd_mo%hetrores_mo_g,               & !< real, dimension(:) :
    autores_mo_g        =>ctem_grd_mo%autores_mo_g,                & !< real, dimension(:) :
    soilres_mo_g        =>ctem_grd_mo%soilres_mo_g,                & !< real, dimension(:) :
    litres_mo_g         =>ctem_grd_mo%litres_mo_g,                 & !< real, dimension(:,:) :
    soilcres_mo_g       =>ctem_grd_mo%soilcres_mo_g,               & !< real, dimension(:,:) :
    vgbiomas_mo_g       =>ctem_grd_mo%vgbiomas_mo_g,               & !< real, dimension(:) :
    totcmass_mo_g       =>ctem_grd_mo%totcmass_mo_g,               & !< real, dimension(:) :
    emit_co2_mo_g       =>ctem_grd_mo%emit_co2_mo_g,               & !< real, dimension(:) :
    emit_co_mo_g        =>ctem_grd_mo%emit_co_mo_g,                & !< real, dimension(:) :
    emit_ch4_mo_g       =>ctem_grd_mo%emit_ch4_mo_g,               & !< real, dimension(:) :
    emit_nmhc_mo_g      =>ctem_grd_mo%emit_nmhc_mo_g,              & !< real, dimension(:) :
    emit_h2_mo_g        =>ctem_grd_mo%emit_h2_mo_g,                & !< real, dimension(:) :
    emit_nox_mo_g       =>ctem_grd_mo%emit_nox_mo_g,               & !< real, dimension(:) :
    emit_n2o_mo_g       =>ctem_grd_mo%emit_n2o_mo_g,               & !< real, dimension(:) :
    emit_nh3_mo_g       =>ctem_grd_mo%emit_nh3_mo_g,               & !< real, dimension(:) :
    emit_pm25_mo_g      =>ctem_grd_mo%emit_pm25_mo_g,              & !< real, dimension(:) :
    emit_tpm_mo_g       =>ctem_grd_mo%emit_tpm_mo_g,               & !< real, dimension(:) :
    emit_tc_mo_g        =>ctem_grd_mo%emit_tc_mo_g,                & !< real, dimension(:) :
    emit_oc_mo_g        =>ctem_grd_mo%emit_oc_mo_g,                & !< real, dimension(:) :
    emit_bc_mo_g        =>ctem_grd_mo%emit_bc_mo_g,                & !< real, dimension(:) :
    smfuncveg_mo_g      =>ctem_grd_mo%smfuncveg_mo_g,              & !< real, dimension(:) :
    luc_emc_mo_g        =>ctem_grd_mo%luc_emc_mo_g,                & !< real, dimension(:) :
    lucltrin_mo_g       =>ctem_grd_mo%lucltrin_mo_g,               & !< real, dimension(:) :
    lucsocin_mo_g       =>ctem_grd_mo%lucsocin_mo_g,               & !< real, dimension(:) :
    luc_emcn_mo_g        =>ctem_grd_mo%luc_emcn_mo_g,                & !< real, dimension(:) :
    lucltrinn_mo_g       =>ctem_grd_mo%lucltrinn_mo_g,               & !< real, dimension(:) :
    lucsocinn_mo_g       =>ctem_grd_mo%lucsocinn_mo_g,               & !< real, dimension(:) :
    burnfrac_mo_g       =>ctem_grd_mo%burnfrac_mo_g,               & !< real, dimension(:) :
    bterm_mo_g          =>ctem_grd_mo%bterm_mo_g,                  & !< real, dimension(:) :
    lterm_mo_g          =>ctem_grd_mo%lterm_mo_g,                  & !< real, dimension(:) :
    mterm_mo_g          =>ctem_grd_mo%mterm_mo_g,                  & !< real, dimension(:) :
    ch4WetSpec_mo_g     =>ctem_grd_mo%ch4WetSpec_mo_g,             & !< real, dimension(:) :
    wetfdyn_mo_g        =>ctem_grd_mo%wetfdyn_mo_g,                & !< real, dimension(:) :
    wetfpres_mo_g       =>ctem_grd_mo%wetfpres_mo_g,               & !< real, dimension(:) :
    ch4WetDyn_mo_g      =>ctem_grd_mo%ch4WetDyn_mo_g,              & !< real, dimension(:) :
    ch4soills_mo_g      =>ctem_grd_mo%ch4soills_mo_g,              & !< real, dimension(:) :
    cProduct_mo_g       =>ctem_grd_mo%cProduct_mo_g,               & !< real, dimension(:) : Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    nProduct_mo_g       =>ctem_grd_mo%nProduct_mo_g,               & !< real, dimension(:) : Nitrogen in the LUC product pools (litter and soil N iccp2 position) \f$[g N m^{-2}]\f$
    fProductDecomp_mo_g =>ctem_grd_mo%fProductDecomp_mo_g,         & !< real, dimension(:) : Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-
    tileAge_mo_g        =>ctem_grd_mo%tileAge_mo_g,                & !< real, dimension(:)

    bnf_tot_mo_g        =>ctem_grd_mo%bnf_tot_mo_g,                & !< real, dimension(:) :
    bnf_free_mo_g       =>ctem_grd_mo%bnf_free_mo_g,               & !< real, dimension(:) :
    bnf_ant_mo_g        =>ctem_grd_mo%bnf_ant_mo_g,                & !< real, dimension(:) :
    bnf_nat_mo_g        =>ctem_grd_mo%bnf_nat_mo_g,                & !< real, dimension(:) :
    nitrif_mo_g         =>ctem_grd_mo%nitrif_mo_g,                 & !< real, dimension(:) :
    no_nit_mo_g         =>ctem_grd_mo%no_nit_mo_g,                 & !< real, dimension(:) :
    no_denit_mo_g       =>ctem_grd_mo%no_denit_mo_g,               & !< real, dimension(:) :
    no_nitdenit_mo_g    =>ctem_grd_mo%no_nitdenit_mo_g,            & !< real, dimension(:) :
    n2o_nit_mo_g        =>ctem_grd_mo%n2o_nit_mo_g,                & !< real, dimension(:) :
    n2o_denit_mo_g      =>ctem_grd_mo%n2o_denit_mo_g,              & !< real, dimension(:) :
    n2o_nitdenit_mo_g   =>ctem_grd_mo%n2o_nitdenit_mo_g,           & !< real, dimension(:) :
    n2_denit_mo_g       =>ctem_grd_mo%n2_denit_mo_g,               & !< real, dimension(:) :
    nvol_mo_g           =>ctem_grd_mo%nvol_mo_g,                   & !< real, dimension(:) :
    nleach_mo_g         =>ctem_grd_mo%nleach_mo_g,                 & !< real, dimension(:) :
    appl_fert_mo_g      =>ctem_grd_mo%appl_fert_mo_g,              & !< real, dimension(:) :
    ndep_nh4_mo_g       =>ctem_grd_mo%ndep_nh4_mo_g,               & !< real, dimension(:) :
    ndep_no3_mo_g       =>ctem_grd_mo%ndep_no3_mo_g,               & !< real, dimension(:) :
    ndemand_wp_npp_mo_g =>ctem_grd_mo%ndemand_wp_npp_mo_g,         & !< real, dimension(:) :
    nuptake_p_nh4_mo_g  =>ctem_grd_mo%nuptake_p_nh4_mo_g,          & !< real, dimension(:) :
    nuptake_p_no3_mo_g  =>ctem_grd_mo%nuptake_p_no3_mo_g,          & !< real, dimension(:) :
    nuptake_a_actl_nh4_mo_g=>ctem_grd_mo%nuptake_a_actl_nh4_mo_g,  & !< real, dimension(:) :
    nuptake_a_actl_no3_mo_g=>ctem_grd_mo%nuptake_a_actl_no3_mo_g,  & !< real, dimension(:) :
    nuptake_mo_g        =>ctem_grd_mo%nuptake_mo_g,                & !< real, dimension(:) :
    nalloc_l_mo_g       =>ctem_grd_mo%nalloc_l_mo_g,               & !< real, dimension(:) :
    nalloc_s_mo_g       =>ctem_grd_mo%nalloc_s_mo_g,               & !< real, dimension(:) :
    nalloc_r_mo_g       =>ctem_grd_mo%nalloc_r_mo_g,               & !< real, dimension(:) :
    nresorped_s_mo_g    =>ctem_grd_mo%nresorped_s_mo_g,            & !< real, dimension(:) :
    nresorped_r_mo_g    =>ctem_grd_mo%nresorped_r_mo_g,            & !< real, dimension(:) :
    nre_alloc_s2l_mo_g  =>ctem_grd_mo%nre_alloc_s2l_mo_g,          & !< real, dimension(:) :
    nre_alloc_r2l_mo_g  =>ctem_grd_mo%nre_alloc_r2l_mo_g,          & !< real, dimension(:) :
    nleafns2s_mo_g      =>ctem_grd_mo%nleafns2s_mo_g,              & !< real, dimension(:) :
    nstemns2s_mo_g      =>ctem_grd_mo%nstemns2s_mo_g,              & !< real, dimension(:) :
    nrootns2s_mo_g      =>ctem_grd_mo%nrootns2s_mo_g,              & !< real, dimension(:) :
    nlitr_l_mo_g        =>ctem_grd_mo%nlitr_l_mo_g,                & !< real, dimension(:) :
    nlitr_s_mo_g        =>ctem_grd_mo%nlitr_s_mo_g,                & !< real, dimension(:) :
    nlitr_r_mo_g        =>ctem_grd_mo%nlitr_r_mo_g,                & !< real, dimension(:) :
    nlitr_mo_g          =>ctem_grd_mo%nlitr_mo_g,                  & !< real, dimension(:) :
    gl2bl_grass_nflux_mo_g=>ctem_grd_mo%gl2bl_grass_nflux_mo_g,    & !< real, dimension(:) :
    nhumtrs_mo_g        =>ctem_grd_mo%nhumtrs_mo_g,                & !< real, dimension(:) :
    nmineral_litr_mo_g  =>ctem_grd_mo%nmineral_litr_mo_g,          & !< real, dimension(:) :
    nmineral_humus_mo_g =>ctem_grd_mo%nmineral_humus_mo_g,         & !< real, dimension(:) :
    netnmineral_mo_g    =>ctem_grd_mo%netnmineral_mo_g,            & !< real, dimension(:) :
    nimmobil_nh4_mo_g   =>ctem_grd_mo%nimmobil_nh4_mo_g,           & !< real, dimension(:) :
    nimmobil_no3_mo_g   =>ctem_grd_mo%nimmobil_no3_mo_g,           & !< real, dimension(:) :
    fNnetland_mo_g      =>ctem_grd_mo%fNnetland_mo_g,              & !< real, dimension(:) :
    nvgbiomas_mo_g      =>ctem_grd_mo%nvgbiomas_mo_g,               & !< real, dimension(:) :
    nstress_mo_g        =>ctem_grd_mo%nstress_mo_g                & !< real, dimension(:) :
    )


    


    !> ------------

    !> Accumulate monthly outputs

    ! Find which month you are in so we know the number of days per month.
    do nt = 1,nmon
      if (iday <= monthend(nt + 1)) then
        oneOverDPM = 1./real(monthdays(nt))
        exit
      else
        cycle
      end if
    end do

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do m = 1,nmtest
      do j = 1,icc

        !> Accumulate monthly outputs at the per PFT level.
        if (ailcgrow(i,m,j) > laimaxg_mo(i,m,j)) then
          laimaxg_mo(i,m,j) = ailcgrow(i,m,j)
        end if

        cc_mo(i,m,j) = cc_mo(i,m,j) + ccrow(i,m,j) * oneOverDPM
        mm_mo(i,m,j) = mm_mo(i,m,j) + mmrow(i,m,j) * oneOverDPM
        npp_mo(i,m,j) = npp_mo(i,m,j) + nppvegrow(i,m,j) * oneOverDPM
        gpp_mo(i,m,j) = gpp_mo(i,m,j) + gppvegrow(i,m,j) * oneOverDPM
        leafns2s_mo(i,m,j) = leafns2s_mo(i,m,j) + leafns2srow(i,m,j) * oneOverDPM
        stemns2s_mo(i,m,j) = stemns2s_mo(i,m,j) + stemns2srow(i,m,j) * oneOverDPM
        rootns2s_mo(i,m,j) = rootns2s_mo(i,m,j) + rootns2srow(i,m,j) * oneOverDPM
        re_alloc_s2l_mo(i,m,j) = re_alloc_s2l_mo(i,m,j) + re_alloc_s2lrow(i,m,j) * oneOverDPM
        re_alloc_r2l_mo(i,m,j) = re_alloc_r2l_mo(i,m,j) + re_alloc_r2lrow(i,m,j) * oneOverDPM
        re_alloc_sr2l_mo(i,m,j) = re_alloc_sr2l_mo(i,m,j) + re_alloc_sr2lrow(i,m,j) * oneOverDPM
        nep_mo(i,m,j) = nep_mo(i,m,j) + nepvegrow(i,m,j) * oneOverDPM
        ! NOTE: This NBP does not include LUC product pool contributions since they are
        ! not per PFT but rather per tile
        litrfallveg_mo(i,m,j) = litrfallveg_mo(i,m,j) + litrfallvegrow(i,m,j) * oneOverDPM
        tltrleaf_mo(i,m,j) = tltrleaf_mo(i,m,j) + tltrleafrow(i,m,j) * oneOverDPM * convertkgC ! convert from umol CO2/m2/s to kgC/m2/s
        tltrstem_mo(i,m,j) = tltrstem_mo(i,m,j) + tltrstemrow(i,m,j) * oneOverDPM * convertkgC ! convert from umol CO2/m2/s to kgC/m2/s
        tltrroot_mo(i,m,j) = tltrroot_mo(i,m,j) + tltrrootrow(i,m,j) * oneOverDPM * convertkgC ! convert from umol CO2/m2/s to kgC/m2/s

        nbp_mo(i,m,j) = nbp_mo(i,m,j) + nbpvegrow(i,m,j) * oneOverDPM
        hetrores_mo(i,m,j) = hetrores_mo(i,m,j) + hetroresvegrow(i,m,j) * oneOverDPM
        autores_mo(i,m,j) = autores_mo(i,m,j) + autoresvegrow(i,m,j) * oneOverDPM
        ! Calculate the total soil respiration (Rh + root contributions of Rm + Rg)
        soilres_mo(i,m,j) = soilres_mo(i,m,j) + (hetroresvegrow(i,m,j) + rmrvegrow(i,m,j) &
                                              + rgvegrow(i,m,j) * afrrootrow(i,m,j)) * oneOverDPM
        do k = 1,ignd
          litres_mo(i,m,j,k) = litres_mo(i,m,j,k) + litresvegrow(i,m,j,k)*oneOverDPM
          soilcres_mo(i,m,j,k) = soilcres_mo(i,m,j,k) + soilcresvegrow(i,m,j,k)*oneOverDPM
          humiftrsveg_mo(i,m,j) = humiftrsveg_mo(i,m,j) + humiftrsvegrow(i,m,j,k)*oneOverDPM
        end do
        emit_co2_mo(i,m,j) = emit_co2_mo(i,m,j) + emit_co2row(i,m,j) * oneOverDPM
        emit_co_mo(i,m,j) = emit_co_mo(i,m,j) + emit_corow(i,m,j) * oneOverDPM
        emit_ch4_mo(i,m,j) = emit_ch4_mo(i,m,j) + emit_ch4row(i,m,j) * oneOverDPM
        emit_nmhc_mo(i,m,j) = emit_nmhc_mo(i,m,j) + emit_nmhcrow(i,m,j) * oneOverDPM
        emit_h2_mo(i,m,j) = emit_h2_mo(i,m,j) + emit_h2row(i,m,j) * oneOverDPM
        emit_nox_mo(i,m,j) = emit_nox_mo(i,m,j) + emit_noxrow(i,m,j) * oneOverDPM
        emit_n2o_mo(i,m,j) = emit_n2o_mo(i,m,j) + emit_n2orow(i,m,j) * oneOverDPM
        emit_nh3_mo(i,m,j) = emit_nh3_mo(i,m,j) + emit_nh3row(i,m,j) * oneOverDPM
        emit_pm25_mo(i,m,j) = emit_pm25_mo(i,m,j) + emit_pm25row(i,m,j) * oneOverDPM
        emit_tpm_mo(i,m,j) = emit_tpm_mo(i,m,j) + emit_tpmrow(i,m,j) * oneOverDPM
        emit_tc_mo(i,m,j) = emit_tc_mo(i,m,j) + emit_tcrow(i,m,j) * oneOverDPM
        emit_oc_mo(i,m,j) = emit_oc_mo(i,m,j) + emit_ocrow(i,m,j) * oneOverDPM
        emit_bc_mo(i,m,j) = emit_bc_mo(i,m,j) + emit_bcrow(i,m,j) * oneOverDPM
        bterm_mo(i,m,j) = bterm_mo(i,m,j) + btermrow(i,m,j) * oneOverDPM
        mterm_mo(i,m,j) = mterm_mo(i,m,j) + mtermrow(i,m,j) * oneOverDPM
        smfuncveg_mo(i,m,j) = smfuncveg_mo(i,m,j) + smfuncvegrow(i,m,j) * oneOverDPM

        ! Let accumulate,not fluxes nor meant to be mean values.
        burnfrac_mo(i,m,j) = burnfrac_mo(i,m,j) + burnvegfrow(i,m,j)

        vcmax0_mo(i,m,j) = vcmax0_mo(i,m,j) + vcmax0row(i,m,j)*oneOverDPM

        ! Note that the following are recorded in units of g N m^-2 mo^-1:
        if (Ncycle_on) then
           bnf_free_mo(i,m,j) = bnf_free_mo(i,m,j) + bnffreerow(i,m,j) * oneOverDPM * convertkgN
           bnf_ant_mo(i,m,j) = bnf_ant_mo(i,m,j) + bnfantrow(i,m,j) * oneOverDPM * convertkgN
           bnf_nat_mo(i,m,j) = bnf_nat_mo(i,m,j) + bnfnatrow(i,m,j) * oneOverDPM * convertkgN
           bnf_tot_mo(i,m,j) = bnf_tot_mo(i,m,j) + bnftotrow(i,m,j) * oneOverDPM * convertkgN
           nstress_mo(i,m,j) = nstress_mo(i,m,j) + nstressrow(i,m,j) * oneOverDPM * convertkgN
           nitrif_mo(i,m,j) = nitrif_mo(i,m,j) + nitrifvegrow(i,m,j) * oneOverDPM * convertkgN
           no_nit_mo(i,m,j) = no_nit_mo(i,m,j) + no_nitvegrow(i,m,j) * oneOverDPM * convertkgN
           no_denit_mo(i,m,j) = no_denit_mo(i,m,j) + no_denitvegrow(i,m,j) * oneOverDPM * convertkgN
           no_nitdenit_mo(i,m,j) = no_nitdenit_mo(i,m,j) + no_nitdenitvegrow(i,m,j) * oneOverDPM * convertkgN
           n2o_nit_mo(i,m,j) = n2o_nit_mo(i,m,j) + n2o_nitvegrow(i,m,j) * oneOverDPM * convertkgN
           n2o_denit_mo(i,m,j) = n2o_denit_mo(i,m,j) + n2o_denitvegrow(i,m,j) * oneOverDPM * convertkgN
           n2o_nitdenit_mo(i,m,j) = n2o_nitdenit_mo(i,m,j) + n2o_nitdenitvegrow(i,m,j) * oneOverDPM * convertkgN
           n2_denit_mo(i,m,j) = n2_denit_mo(i,m,j) + n2_denitvegrow(i,m,j) * oneOverDPM * convertkgN
           nvol_mo(i,m,j) = nvol_mo(i,m,j) + nvolvegrow(i,m,j) * oneOverDPM * convertkgN
           nleach_mo(i,m,j) = nleach_mo(i,m,j) + nleachvegrow(i,m,j) * oneOverDPM * convertkgN
           appl_fert_mo(i,m,j) = appl_fert_mo(i,m,j) + appl_fertrow(i,m,j) * oneOverDPM * convertkgN
           ndep_nh4_mo(i,m,j) = ndep_nh4_mo(i,m,j) + ndep_nh4row(i,m,j) * oneOverDPM * convertkgN
           ndep_no3_mo(i,m,j) = ndep_no3_mo(i,m,j) + ndep_no3row(i,m,j) * oneOverDPM * convertkgN
           ndemand_wp_npp_mo(i,m,j) = ndemand_wp_npp_mo(i,m,j) + ndemandveg_wp_npprow(i,m,j) * oneOverDPM * convertkgN
           nuptake_p_nh4_mo(i,m,j) = nuptake_p_nh4_mo(i,m,j) + nuptakeveg_p_nh4row(i,m,j) * oneOverDPM * convertkgN
           nuptake_p_no3_mo(i,m,j) = nuptake_p_no3_mo(i,m,j) + nuptakeveg_p_no3row(i,m,j) * oneOverDPM * convertkgN
           nuptake_a_actl_nh4_mo(i,m,j) = nuptake_a_actl_nh4_mo(i,m,j) + nuptakeveg_a_actl_nh4row(i,m,j) * oneOverDPM * convertkgN
           nuptake_a_actl_no3_mo(i,m,j) = nuptake_a_actl_no3_mo(i,m,j) + nuptakeveg_a_actl_no3row(i,m,j) * oneOverDPM * convertkgN
           nuptake_mo(i,m,j) = nuptake_mo(i,m,j) + nuptakevegrow(i,m,j) * oneOverDPM * convertkgN
           nalloc_l_mo(i,m,j) = nalloc_l_mo(i,m,j) + nallocveg_lrow(i,m,j) * oneOverDPM * convertkgN
           nalloc_s_mo(i,m,j) = nalloc_s_mo(i,m,j) + nallocveg_srow(i,m,j) * oneOverDPM * convertkgN
           nalloc_r_mo(i,m,j) = nalloc_r_mo(i,m,j) + nallocveg_rrow(i,m,j) * oneOverDPM * convertkgN
           nresorped_s_mo(i,m,j) = nresorped_s_mo(i,m,j) + nresorpedveg_srow(i,m,j) * oneOverDPM * convertkgN
           nresorped_r_mo(i,m,j) = nresorped_r_mo(i,m,j) + nresorpedveg_rrow(i,m,j) * oneOverDPM * convertkgN
           nre_alloc_s2l_mo(i,m,j) = nre_alloc_s2l_mo(i,m,j) + nre_allocveg_s2lrow(i,m,j) * oneOverDPM * convertkgN
           nre_alloc_r2l_mo(i,m,j) = nre_alloc_r2l_mo(i,m,j) + nre_allocveg_r2lrow(i,m,j) * oneOverDPM * convertkgN
           nleafns2s_mo(i,m,j) = nleafns2s_mo(i,m,j) + nleafns2svegrow(i,m,j) * oneOverDPM * convertkgN
           nstemns2s_mo(i,m,j) = nstemns2s_mo(i,m,j) + nstemns2svegrow(i,m,j) * oneOverDPM * convertkgN
           nrootns2s_mo(i,m,j) = nrootns2s_mo(i,m,j) + nrootns2svegrow(i,m,j) * oneOverDPM * convertkgN
           nlitr_l_mo(i,m,j) = nlitr_l_mo(i,m,j) + nlitrveg_lrow(i,m,j) * oneOverDPM * convertkgN
           nlitr_s_mo(i,m,j) = nlitr_s_mo(i,m,j) + nlitrveg_srow(i,m,j) * oneOverDPM * convertkgN
           nlitr_r_mo(i,m,j) = nlitr_r_mo(i,m,j) + nlitrveg_rrow(i,m,j) * oneOverDPM * convertkgN
           nlitr_mo(i,m,j) = nlitr_mo(i,m,j) + nlitrvegrow(i,m,j) * oneOverDPM * convertkgN
           gl2bl_grass_nflux_mo(i,m,j) = gl2bl_grass_nflux_mo(i,m,j) + gl2bl_grass_nfluxrow(i,m,j) * oneOverDPM * convertkgN
           nhumtrs_mo(i,m,j) = nhumtrs_mo(i,m,j) + nhumtrsvegrow(i,m,j) * oneOverDPM * convertkgN
           nmineral_litr_mo(i,m,j) = nmineral_litr_mo(i,m,j) + nmineralveg_litrrow(i,m,j) * oneOverDPM * convertkgN
           nmineral_humus_mo(i,m,j) = nmineral_humus_mo(i,m,j) + nmineralveg_humusrow(i,m,j) * oneOverDPM * convertkgN
           netnmineral_mo(i,m,j) = netnmineral_mo(i,m,j) + netnmineralveg_row(i,m,j) * oneOverDPM * convertkgN
           nimmobil_nh4_mo(i,m,j) = nimmobil_nh4_mo(i,m,j) + nimmobilveg_nh4row(i,m,j) * oneOverDPM * convertkgN
           nimmobil_no3_mo(i,m,j) = nimmobil_no3_mo(i,m,j) + nimmobilveg_no3row(i,m,j) * oneOverDPM * convertkgN
           fNnetland_mo(i,m,j) = fNnetland_mo(i,m,j) + fNnetlandvegrow(i,m,j) * oneOverDPM * convertkgN
        end if

      end do ! j

      !> Also do the bare ground
      nep_mo(i,m,iccp1) = nep_mo(i,m,iccp1) + nepvegrow(i,m,iccp1) * oneOverDPM
      nbp_mo(i,m,iccp1) = nbp_mo(i,m,iccp1) + nbpvegrow(i,m,iccp1) * oneOverDPM
      hetrores_mo(i,m,iccp1) = hetrores_mo(i,m,iccp1) + hetroresvegrow(i,m,iccp1) * oneOverDPM
      do k = 1,ignd
        litres_mo(i,m,iccp1,k) = litres_mo(i,m,iccp1,k) + litresvegrow(i,m,iccp1,k) * oneOverDPM
        soilcres_mo(i,m,iccp1,k) = soilcres_mo(i,m,iccp1,k) &
                                 + soilcresvegrow(i,m,iccp1,k) * oneOverDPM
      end do
      if (Ncycle_on) then
        bnf_free_mo(i,m,iccp1) = bnf_free_mo(i,m,iccp1) + bnffreerow(i,m,iccp1) * oneOverDPM * convertkgN
        bnf_tot_mo(i,m,iccp1) = bnf_tot_mo(i,m,iccp1) + bnftotrow(i,m,iccp1) * oneOverDPM * convertkgN
        nitrif_mo(i,m,iccp1) = nitrif_mo(i,m,iccp1) + nitrifvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        no_nit_mo(i,m,iccp1) = no_nit_mo(i,m,iccp1) + no_nitvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        no_denit_mo(i,m,iccp1) = no_denit_mo(i,m,iccp1) + no_denitvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        no_nitdenit_mo(i,m,iccp1) = no_nitdenit_mo(i,m,iccp1) + no_nitdenitvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        n2o_nit_mo(i,m,iccp1) = n2o_nit_mo(i,m,iccp1) + n2o_nitvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        n2o_denit_mo(i,m,iccp1) = n2o_denit_mo(i,m,iccp1) + n2o_denitvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        n2o_nitdenit_mo(i,m,iccp1) = n2o_nitdenit_mo(i,m,iccp1) + n2o_nitdenitvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        n2_denit_mo(i,m,iccp1) = n2_denit_mo(i,m,iccp1) + n2_denitvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        nvol_mo(i,m,iccp1) = nvol_mo(i,m,iccp1) + nvolvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        nleach_mo(i,m,iccp1) = nleach_mo(i,m,iccp1) + nleachvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        appl_fert_mo(i,m,iccp1) = appl_fert_mo(i,m,iccp1) + appl_fertrow(i,m,iccp1) * oneOverDPM * convertkgN
        ndep_nh4_mo(i,m,iccp1) = ndep_nh4_mo(i,m,iccp1) + ndep_nh4row(i,m,iccp1) * oneOverDPM * convertkgN
        ndep_no3_mo(i,m,iccp1) = ndep_no3_mo(i,m,iccp1) + ndep_no3row(i,m,iccp1) * oneOverDPM * convertkgN
        nhumtrs_mo(i,m,iccp1) = nhumtrs_mo(i,m,iccp1) + nhumtrsvegrow(i,m,iccp1) * oneOverDPM * convertkgN
        nmineral_litr_mo(i,m,iccp1) = nmineral_litr_mo(i,m,iccp1) + nmineralveg_litrrow(i,m,iccp1) * oneOverDPM * convertkgN
        nmineral_humus_mo(i,m,iccp1) = nmineral_humus_mo(i,m,iccp1) + nmineralveg_humusrow(i,m,iccp1) * oneOverDPM * convertkgN
        netnmineral_mo(i,m,iccp1) = netnmineral_mo(i,m,iccp1) + netnmineralveg_row(i,m,iccp1) * oneOverDPM * convertkgN
        nimmobil_nh4_mo(i,m,iccp1) = nimmobil_nh4_mo(i,m,iccp1) + nimmobilveg_nh4row(i,m,iccp1) * oneOverDPM * convertkgN
        nimmobil_no3_mo(i,m,iccp1) = nimmobil_no3_mo(i,m,iccp1) + nimmobilveg_no3row(i,m,iccp1) * oneOverDPM * convertkgN
        fNnetland_mo(i,m,iccp1) = fNnetland_mo(i,m,iccp1) + fNnetlandvegrow(i,m,iccp1) * oneOverDPM * convertkgN
      end if ! Ncycle_on

      !> Accumulate monthly outputs at the per tile level.
      luc_emc_mo_t(i,m) = luc_emc_mo_t(i,m) + lucemcomrow(i,m) * oneOverDPM
      lucsocin_mo_t(i,m) = lucsocin_mo_t(i,m) + lucsocinrow(i,m) * oneOverDPM
      lucltrin_mo_t(i,m) = lucltrin_mo_t(i,m) + lucltrinrow(i,m) * oneOverDPM
      tileAge_mo_t(i,m) = tileAge_mo_t(i,m) + tileAgerot(i,m) * oneOverDPM
      if (Ncycle_on) then
        luc_emcn_mo_t(i,m) = luc_emcn_mo_t(i,m) + lucemcomnrow(i,m) * oneOverDPM
        lucsocinn_mo_t(i,m) = lucsocinn_mo_t(i,m) + lucsocinnrow(i,m) * oneOverDPM
        lucltrinn_mo_t(i,m) = lucltrinn_mo_t(i,m) + lucltrinnrow(i,m) * oneOverDPM
      end if
      ch4WetSpec_mo_t(i,m) = ch4WetSpec_mo_t(i,m) + ch4WetSpecrow(i,m) * oneOverDPM
      wetfdyn_mo_t(i,m) = wetfdyn_mo_t(i,m) + wetfdynrow(i,m) * oneOverDPM
      wetfpres_mo_t(i,m) = wetfpres_mo_t(i,m) + wetfrac_presrow(i,m) * oneOverDPM
      ch4WetDyn_mo_t(i,m) = ch4WetDyn_mo_t(i,m) + ch4WetDynrow(i,m) * oneOverDPM
      ch4soills_mo_t(i,m) = ch4soills_mo_t(i,m) + ch4soillsrow(i,m) * oneOverDPM
      lterm_mo_t(i,m) = lterm_mo_t(i,m) + ltermrow(i,m) * oneOverDPM
      ! wind_mo_t(i,m) = wind_mo_t(i,m) + (sqrt(uvaccrow_m(i,m)**2.0 + vvaccrow_m(i,m)**2.0))*3.6 !> take mean wind speed and convert to km/h
      ! NOTE: LUC product pools are kept in layer 1.
      fProductDecomp_mo_t(i,m) = fProductDecomp_mo_t(i,m) &
                                + (soilcresvegrow(i,m,iccp2,1) + litresvegrow(i,m,iccp2,1)) * oneOverDPM

      ! NOTE: NBP is a special case here. The LUC product pool contributions are not
      ! per PFT, they exist uniformly across a tile, so they are not inclued in the
      ! nbp_mo calculation. Instead we need to use the nbp, not nbpveg variable
      ! for per tile and per gridcell outputting.
      nbp_mo_t(i,m) = nbp_mo_t(i,m) + nbprow(i,m) * oneOverDPM

    end do ! loop 863 ! m

    do nt = 1,nmon

      if (iday == mmday(nt)) then

        !> Do the mid-month variables (these are not accumulated,we just keep the mid month value for printing in the monthly file)

        do m = 1,nmtest
          do j = 1,icc

            vgbiomas_mo(i,m,j) = vgbiomas_vegrow(i,m,j)
            stemmass_mo(i,m,j) = stemmassrow(i,m,j)
            stemmass_NS_mo(i,m,j) = stemmass_NSrow(i,m,j)
            stemmasss_mo(i,m,j) = stemmasssrow(i,m,j)
            gleafmas_mo(i,m,j) = gleafmasrow(i,m,j)
            gleafmas_NS_mo(i,m,j) = gleafmas_NSrow(i,m,j)
            gleafmass_mo(i,m,j) = gleafmassrow(i,m,j)
            rootmass_mo(i,m,j) = rootmassrow(i,m,j)
            rootmass_NS_mo(i,m,j) = rootmass_NSrow(i,m,j)
            rootmasss_mo(i,m,j) = rootmasssrow(i,m,j)
            rootdpth_mo(i,m,j) = rootdpthrow(i,m,j)
            totcmass_mo(i,m,j) = vgbiomas_vegrow(i,m,j)
            do k = 1,ignd
              litrmass_mo(i,m,j,k) = litrmassrow(i,m,j,k)
              soilcmas_mo(i,m,j,k) = soilcmasrow(i,m,j,k)
              totcmass_mo(i,m,j) = totcmass_mo(i,m,j) + litrmassrow(i,m,j,k) + soilcmasrow(i,m,j,k)
            end do
            if (Ncycle_on) then
               nh4_mass_mo(i,m,j) = nh4_massrow(i,m,j) * convertg2kg
               no3_mass_mo(i,m,j) = no3_massrow(i,m,j) * convertg2kg
               ngleafmas_mo(i,m,j) = ngleafmasrow(i,m,j) * convertg2kg
               ngleafmas_NS_mo(i,m,j) = ngleafmas_NSrow(i,m,j) * convertg2kg
               ngleafmass_mo(i,m,j) = ngleafmassrow(i,m,j) * convertg2kg
               nbleafmas_mo(i,m,j) = nbleafmasrow(i,m,j) * convertg2kg
               nstemmass_mo(i,m,j) = nstemmassrow(i,m,j) * convertg2kg
               nstemmass_NS_mo(i,m,j) = nstemmass_NSrow(i,m,j) * convertg2kg
               nstemmasss_mo(i,m,j) = nstemmasssrow(i,m,j) * convertg2kg
               nrootmass_mo(i,m,j) = nrootmassrow(i,m,j) * convertg2kg
               nrootmass_NS_mo(i,m,j) = nrootmass_NSrow(i,m,j) * convertg2kg
               nrootmasss_mo(i,m,j) = nrootmasssrow(i,m,j) * convertg2kg
               nlitrmass_mo(i,m,j) = nlitrmassrow(i,m,j) * convertg2kg
               soilnmas_mo(i,m,j) = soilnmasrow(i,m,j) * convertg2kg
               nvgbiomas_mo(i,m,j) = nvgbiomas_vegrow(i,m,j) * convertg2kg
               if (ngleafmas_mo(i,m,j) /= 0.0) c2n_l_mo(i,m,j) = gleafmas_mo(i,m,j) / ngleafmas_mo(i,m,j)
               if (nstemmass_mo(i,m,j) /= 0.0) c2n_s_mo(i,m,j) = stemmass_mo(i,m,j) / nstemmass_mo(i,m,j)
               if (nrootmass_mo(i,m,j) /= 0.0) c2n_r_mo(i,m,j) = rootmass_mo(i,m,j) / nrootmass_mo(i,m,j)
               if ((ngleafmas_mo(i,m,j) + nstemmass_mo(i,m,j) + nrootmass_mo(i,m,j)) /= 0.0) &
                    c2n_wp_mo(i,m,j) = (gleafmas_mo(i,m,j) + stemmass_mo(i,m,j) + rootmass_mo(i,m,j)) &
                                 / (ngleafmas_mo(i,m,j) + nstemmass_mo(i,m,j) + nrootmass_mo(i,m,j))
               if (nlitrmass_mo(i,m,j) /= 0.0) c2n_litr_mo(i,m,j) = sum(litrmass_mo(i,m,j,:)) / nlitrmass_mo(i,m,j)
               if (soilnmas_mo(i,m,j) /= 0.0) c2n_humus_mo(i,m,j) = sum(soilcmas_mo(i,m,j,:)) / soilnmas_mo(i,m,j)
            end if

          end do ! j

          !> Do the bare fraction too
          do k = 1,ignd
            litrmass_mo(i,m,iccp1,k) = litrmassrow(i,m,iccp1,k)
            soilcmas_mo(i,m,iccp1,k) = soilcmasrow(i,m,iccp1,k)
            totcmass_mo(i,m,iccp1) = totcmass_mo(i,m,iccp1) &
                                    + soilcmasrow(i,m,iccp1,k) + litrmassrow(i,m,iccp1,k)
          end do
          if (Ncycle_on) then
            nh4_mass_mo(i,m,iccp1) = nh4_massrow(i,m,iccp1) * convertg2kg
            no3_mass_mo(i,m,iccp1) = no3_massrow(i,m,iccp1) * convertg2kg
            nlitrmass_mo(i,m,iccp1) = nlitrmassrow(i,m,iccp1) * convertg2kg
            soilnmas_mo(i,m,iccp1) = soilnmasrow(i,m,iccp1) * convertg2kg
            if (nlitrmass_mo(i,m,iccp1) /= 0.0) c2n_litr_mo(i,m,iccp1) = sum(litrmass_mo(i,m,iccp1,:)) / nlitrmass_mo(i,m,iccp1)
            if (soilnmas_mo(i,m,iccp1) /= 0.0) c2n_humus_mo(i,m,iccp1) = sum(soilcmas_mo(i,m,iccp1,:)) / soilnmas_mo(i,m,iccp1)
          end if

          barefrac = 1.0

          !> Now find the per tile values:
          do j = 1,icc
            vgbiomas_mo_t(i,m) = vgbiomas_mo_t(i,m) + vgbiomas_mo(i,m,j) * fcancmxrow(i,m,j)
            do k = 1,ignd
              litrmass_mo_t(i,m,k) = litrmass_mo_t(i,m,k) + litrmass_mo(i,m,j,k) * fcancmxrow(i,m,j)
              soilcmas_mo_t(i,m,k) = soilcmas_mo_t(i,m,k) + soilcmas_mo(i,m,j,k) * fcancmxrow(i,m,j)
            end do
            gleafmas_mo_t(i,m) = gleafmas_mo_t(i,m) + gleafmas_mo(i,m,j) * fcancmxrow(i,m,j)
            gleafmas_NS_mo_t(i,m) = gleafmas_NS_mo_t(i,m) + gleafmas_NS_mo(i,m,j) * fcancmxrow(i,m,j)
            gleafmass_mo_t(i,m) = gleafmass_mo_t(i,m) + gleafmass_mo(i,m,j) * fcancmxrow(i,m,j)
            stemmass_mo_t(i,m) = stemmass_mo_t(i,m) + stemmass_mo(i,m,j) * fcancmxrow(i,m,j)
            stemmass_NS_mo_t(i,m) = stemmass_NS_mo_t(i,m) + stemmass_NS_mo(i,m,j) * fcancmxrow(i,m,j)
            stemmasss_mo_t(i,m) = stemmasss_mo_t(i,m) + stemmasss_mo(i,m,j) * fcancmxrow(i,m,j)
            rootmass_mo_t(i,m) = rootmass_mo_t(i,m) + rootmass_mo(i,m,j) * fcancmxrow(i,m,j)
            rootmass_NS_mo_t(i,m) = rootmass_NS_mo_t(i,m) + rootmass_NS_mo(i,m,j) * fcancmxrow(i,m,j)
            rootmasss_mo_t(i,m) = rootmasss_mo_t(i,m) + rootmasss_mo(i,m,j) * fcancmxrow(i,m,j)
            totcmass_mo_t(i,m) = totcmass_mo_t(i,m) + totcmass_mo(i,m,j) * fcancmxrow(i,m,j)
            if (Ncycle_on) then
              nh4_mass_mo_t(i,m) = nh4_mass_mo_t(i,m) + nh4_mass_mo(i,m,j) * fcancmxrow(i,m,j)
              no3_mass_mo_t(i,m) = no3_mass_mo_t(i,m) + no3_mass_mo(i,m,j) * fcancmxrow(i,m,j)
              nlitrmass_mo_t(i,m) = nlitrmass_mo_t(i,m) + nlitrmass_mo(i,m,j) * fcancmxrow(i,m,j)
              soilnmas_mo_t(i,m) = soilnmas_mo_t(i,m) + soilnmas_mo(i,m,j) * fcancmxrow(i,m,j)
              ngleafmas_mo_t(i,m) = ngleafmas_mo_t(i,m) + ngleafmas_mo(i,m,j) * fcancmxrow(i,m,j)
              ngleafmas_NS_mo_t(i,m) = ngleafmas_NS_mo_t(i,m) + ngleafmas_NS_mo(i,m,j) * fcancmxrow(i,m,j)
              ngleafmass_mo_t(i,m) = ngleafmass_mo_t(i,m) + ngleafmass_mo(i,m,j) * fcancmxrow(i,m,j)
              nbleafmas_mo_t(i,m) = nbleafmas_mo_t(i,m) + nbleafmas_mo(i,m,j) * fcancmxrow(i,m,j)
              nstemmass_mo_t(i,m) = nstemmass_mo_t(i,m) + nstemmass_mo(i,m,j) * fcancmxrow(i,m,j)
              nstemmass_NS_mo_t(i,m) = nstemmass_NS_mo_t(i,m) + nstemmass_NS_mo(i,m,j) * fcancmxrow(i,m,j)
              nstemmasss_mo_t(i,m) = nstemmasss_mo_t(i,m) + nstemmasss_mo(i,m,j) * fcancmxrow(i,m,j)
              nrootmass_mo_t(i,m) = nrootmass_mo_t(i,m) + nrootmass_mo(i,m,j) * fcancmxrow(i,m,j)
              nrootmass_NS_mo_t(i,m) = nrootmass_NS_mo_t(i,m) + nrootmass_NS_mo(i,m,j) * fcancmxrow(i,m,j)
              nrootmasss_mo_t(i,m) = nrootmasss_mo_t(i,m) + nrootmasss_mo(i,m,j) * fcancmxrow(i,m,j)
              nvgbiomas_mo_t(i,m) = nvgbiomas_mo_t(i,m) + nvgbiomas_mo(i,m,j) * fcancmxrow(i,m,j)
            end if !Ncycle_on
            barefrac = barefrac - fcancmxrow(i,m,j)
          end do !j
          if(Ncycle_on) then
            if (ngleafmas_mo_t(i,m) /= 0.0) c2n_l_mo_t(i,m) = gleafmas_mo_t(i,m) / ngleafmas_mo_t(i,m)
            if (nstemmass_mo_t(i,m) /= 0.0) c2n_s_mo_t(i,m) = stemmass_mo_t(i,m) / nstemmass_mo_t(i,m)
            if (nrootmass_mo_t(i,m) /= 0.0) c2n_r_mo_t(i,m) = rootmass_mo_t(i,m) / nrootmass_mo_t(i,m)
            if ((ngleafmas_mo_t(i,m) + nstemmass_mo_t(i,m) + nrootmass_mo_t(i,m)) /= 0.0) &
                 c2n_wp_mo_t(i,m) = (gleafmas_mo_t(i,m) + stemmass_mo_t(i,m) + rootmass_mo_t(i,m)) &
                            / (ngleafmas_mo_t(i,m) + nstemmass_mo_t(i,m) + nrootmass_mo_t(i,m))
          end if

          !> Also add in the bare fraction contributions.
          do k = 1,ignd
            litrmass_mo_t(i,m,k) = litrmass_mo_t(i,m,k) + litrmass_mo(i,m,iccp1,k) * barefrac
            soilcmas_mo_t(i,m,k) = soilcmas_mo_t(i,m,k) + soilcmas_mo(i,m,iccp1,k) * barefrac
            totcmass_mo_t(i,m) = totcmass_mo_t(i,m) &
                                + (litrmass_mo(i,m,iccp1,k) + soilcmas_mo(i,m,iccp1,k)) * barefrac
          end do
          if (Ncycle_on) then
            nh4_mass_mo_t(i,m) = nh4_mass_mo_t(i,m) + nh4_mass_mo(i,m,iccp1) * barefrac
            no3_mass_mo_t(i,m) = no3_mass_mo_t(i,m) + no3_mass_mo(i,m,iccp1) * barefrac
            nlitrmass_mo_t(i,m)= nlitrmass_mo_t(i,m)+ nlitrmass_mo(i,m,iccp1) * barefrac
            soilnmas_mo_t(i,m) = soilnmas_mo_t(i,m) + soilnmas_mo(i,m,iccp1) * barefrac
            if (nlitrmass_mo_t(i,m) /= 0.0) c2n_litr_mo_t(i,m) = sum(litrmass_mo_t(i,m,:)) / nlitrmass_mo_t(i,m)
            if (soilnmas_mo_t(i,m) /= 0.0) c2n_humus_mo_t(i,m) = sum(soilcmas_mo_t(i,m,:)) / soilnmas_mo_t(i,m)
          end if

          !> Now find the gridcell level values:
          vgbiomas_mo_g(i) = vgbiomas_mo_g(i) + vgbiomas_mo_t(i,m) * FAREROT(i,m)
          do k = 1,ignd
            litrmass_mo_g(i,k) = litrmass_mo_g(i,k) + litrmass_mo_t(i,m,k) * FAREROT(i,m)
            soilcmas_mo_g(i,k) = soilcmas_mo_g(i,k) + soilcmas_mo_t(i,m,k) * FAREROT(i,m)
          end do
          gleafmas_mo_g(i) = gleafmas_mo_g(i) + gleafmas_mo_t(i,m) * FAREROT(i,m)
          gleafmas_NS_mo_g(i) = gleafmas_NS_mo_g(i) + gleafmas_NS_mo_t(i,m) * FAREROT(i,m)
          gleafmass_mo_g(i) = gleafmass_mo_g(i) + gleafmass_mo_t(i,m) * FAREROT(i,m)
          stemmass_mo_g(i) = stemmass_mo_g(i) + stemmass_mo_t(i,m) * FAREROT(i,m)
          stemmass_NS_mo_g(i) = stemmass_NS_mo_g(i) + stemmass_NS_mo_t(i,m) * FAREROT(i,m)
          stemmasss_mo_g(i) = stemmasss_mo_g(i) + stemmasss_mo_t(i,m) * FAREROT(i,m)
          rootmass_mo_g(i) = rootmass_mo_g(i) + rootmass_mo_t(i,m) * FAREROT(i,m)
          rootmass_NS_mo_g(i) = rootmass_NS_mo_g(i) + rootmass_NS_mo_t(i,m) * FAREROT(i,m)
          rootmasss_mo_g(i) = rootmasss_mo_g(i) + rootmasss_mo_t(i,m) * FAREROT(i,m)
          totcmass_mo_g(i) = totcmass_mo_g(i) + totcmass_mo_t(i,m) * FAREROT(i,m)

          ! Including the LUC product pools. They are per tile values and
          ! are assumed to occupy the whole tile. Only kept in layer 1.
          cProduct_mo_g(i) = cProduct_mo_g(i) &
                             + (litrmassrow(i,m,iccp2,1) &
                             + soilcmasrow(i,m,iccp2,1)) * FAREROT(i,m)
          if (Ncycle_on) then
            nh4_mass_mo_g(i) = nh4_mass_mo_g(i) + nh4_mass_mo_t(i,m) * FAREROT(i,m)
            no3_mass_mo_g(i) = no3_mass_mo_g(i) + no3_mass_mo_t(i,m) * FAREROT(i,m)
            nlitrmass_mo_g(i) = nlitrmass_mo_g(i) + nlitrmass_mo_t(i,m) * FAREROT(i,m)
            soilnmas_mo_g(i) = soilnmas_mo_g(i) + soilnmas_mo_t(i,m) * FAREROT(i,m)
            ngleafmas_mo_g(i) = ngleafmas_mo_g(i) + ngleafmas_mo_t(i,m) * FAREROT(i,m)
            ngleafmas_NS_mo_g(i) = ngleafmas_NS_mo_g(i) + ngleafmas_NS_mo_t(i,m) * FAREROT(i,m)
            ngleafmass_mo_g(i) = ngleafmass_mo_g(i) + ngleafmass_mo_t(i,m) * FAREROT(i,m)
            nbleafmas_mo_g(i) = nbleafmas_mo_g(i) + nbleafmas_mo_t(i,m) * FAREROT(i,m)
            nstemmass_mo_g(i) = nstemmass_mo_g(i) + nstemmass_mo_t(i,m) * FAREROT(i,m)
            nstemmass_NS_mo_g(i) = nstemmass_NS_mo_g(i) + nstemmass_NS_mo_t(i,m) * FAREROT(i,m)
            nstemmasss_mo_g(i) = nstemmasss_mo_g(i) + nstemmasss_mo_t(i,m) * FAREROT(i,m)
            nrootmass_mo_g(i) = nrootmass_mo_g(i) + nrootmass_mo_t(i,m) * FAREROT(i,m)
            nrootmass_NS_mo_g(i) = nrootmass_NS_mo_g(i) + nrootmass_NS_mo_t(i,m) * FAREROT(i,m)
            nrootmasss_mo_g(i) = nrootmasss_mo_g(i) + nrootmasss_mo_t(i,m) * FAREROT(i,m)
            nvgbiomas_mo_g(i) = nvgbiomas_mo_g(i) + nvgbiomas_mo_t(i,m) * FAREROT(i,m)
            nProduct_mo_g(i) = nProduct_mo_g(i) &
                               + (nlitrmassrow(i,m,iccp2) &
                               + soilnmasrow(i,m,iccp2)) * convertg2kg * FAREROT(i,m)
          end if !for Ncycle_on

        end do ! m 
        if(Ncycle_on) then
          if (ngleafmas_mo_g(i) /= 0.0) c2n_l_mo_g(i) = gleafmas_mo_g(i) / ngleafmas_mo_g(i)
          if (nstemmass_mo_g(i) /= 0.0) c2n_s_mo_g(i) = stemmass_mo_g(i) / nstemmass_mo_g(i)
          if (nrootmass_mo_g(i) /= 0.0) c2n_r_mo_g(i) = rootmass_mo_g(i) / nrootmass_mo_g(i)
          if ((ngleafmas_mo_g(i) + nstemmass_mo_g(i) + nrootmass_mo_g(i)) /= 0.0) &
              c2n_wp_mo_g(i) = (gleafmas_mo_g(i) + stemmass_mo_g(i) + rootmass_mo_g(i)) &
                           / (ngleafmas_mo_g(i) + nstemmass_mo_g(i) + nrootmass_mo_g(i))
          if (nlitrmass_mo_g(i) /= 0.0) c2n_litr_mo_g(i) = sum(litrmass_mo_g(i,:)) / nlitrmass_mo_g(i)
          if (soilnmas_mo_g(i) /= 0.0) c2n_humus_mo_g(i) = sum(soilcmas_mo_g(i,:)) / soilnmas_mo_g(i)
        end if

      end if ! mmday (mid-month instantaneous value)

      if (iday == monthend(nt + 1)) then

        !> Do the end of month variables
        ndmonth = (monthend(nt + 1) - monthend(nt)) * nday

        do m = 1,nmtest


            !> Convert some quantities into per day values
            !wetfdyn_mo_t(i,m)=wetfdyn_mo_t(i,m)*(1./real(monthdays(nt)))
            !wetfpres_mo_t(i,m)=wetfpres_mo_t(i,m)*(1./real(monthdays(nt)))
            lterm_mo_t(i,m)=lterm_mo_t(i,m)*(1./real(monthdays(nt)))
            !wind_mo_t(i,m) = wind_mo_t(i,m)*(1./real(monthdays(nt)))

          barefrac = 1.0

          do j = 1,icc

            !> Find the monthly outputs at the per tile level from the outputs at the per PFT level
            npp_mo_t(i,m) = npp_mo_t(i,m) + npp_mo(i,m,j) * fcancmxrow(i,m,j)
            gpp_mo_t(i,m) = gpp_mo_t(i,m) + gpp_mo(i,m,j) * fcancmxrow(i,m,j)
            leafns2s_mo_t(i,m) = leafns2s_mo_t(i,m) + leafns2s_mo(i,m,j) * fcancmxrow(i,m,j)
            stemns2s_mo_t(i,m) = stemns2s_mo_t(i,m) + stemns2s_mo(i,m,j) * fcancmxrow(i,m,j)
            rootns2s_mo_t(i,m) = rootns2s_mo_t(i,m) + rootns2s_mo(i,m,j) * fcancmxrow(i,m,j)
            re_alloc_s2l_mo_t(i,m) = re_alloc_s2l_mo_t(i,m) + re_alloc_s2l_mo(i,m,j) * fcancmxrow(i,m,j)
            re_alloc_r2l_mo_t(i,m) = re_alloc_r2l_mo_t(i,m) + re_alloc_r2l_mo(i,m,j) * fcancmxrow(i,m,j)
            re_alloc_sr2l_mo_t(i,m) = re_alloc_sr2l_mo_t(i,m)+ re_alloc_sr2l_mo(i,m,j)* fcancmxrow(i,m,j)
            nep_mo_t(i,m) = nep_mo_t(i,m) + nep_mo(i,m,j) * fcancmxrow(i,m,j)
            ! nbp_mo_t(i,m)=nbp_mo_t(i,m)+nbp_mo(i,m,j)*fcancmxrow(i,m,j)
            hetrores_mo_t(i,m) = hetrores_mo_t(i,m) + hetrores_mo(i,m,j) * fcancmxrow(i,m,j)
            autores_mo_t(i,m) = autores_mo_t(i,m) + autores_mo(i,m,j) * fcancmxrow(i,m,j)
            soilres_mo_t(i,m) = soilres_mo_t(i,m) + soilres_mo(i,m,j) * fcancmxrow(i,m,j)

            do k = 1,ignd
              litres_mo_t(i,m,k) = litres_mo_t(i,m,k) + litres_mo(i,m,j,k) * fcancmxrow(i,m,j)
              soilcres_mo_t(i,m,k) = soilcres_mo_t(i,m,k) + soilcres_mo(i,m,j,k) * fcancmxrow(i,m,j)
            end do
            emit_co2_mo_t(i,m) = emit_co2_mo_t(i,m) + emit_co2_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_co_mo_t(i,m) = emit_co_mo_t(i,m) + emit_co_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_ch4_mo_t(i,m) = emit_ch4_mo_t(i,m) + emit_ch4_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_nmhc_mo_t(i,m) = emit_nmhc_mo_t(i,m) + emit_nmhc_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_h2_mo_t(i,m) = emit_h2_mo_t(i,m) + emit_h2_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_nox_mo_t(i,m) = emit_nox_mo_t(i,m) + emit_nox_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_n2o_mo_t(i,m) = emit_n2o_mo_t(i,m) + emit_n2o_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_nh3_mo_t(i,m) = emit_nh3_mo_t(i,m) + emit_nh3_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_pm25_mo_t(i,m) = emit_pm25_mo_t(i,m) + emit_pm25_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_tpm_mo_t(i,m) = emit_tpm_mo_t(i,m) + emit_tpm_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_tc_mo_t(i,m) = emit_tc_mo_t(i,m) + emit_tc_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_oc_mo_t(i,m) = emit_oc_mo_t(i,m) + emit_oc_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_bc_mo_t(i,m) = emit_bc_mo_t(i,m) + emit_bc_mo(i,m,j) * fcancmxrow(i,m,j)
            bterm_mo_t(i,m) = bterm_mo_t(i,m) + bterm_mo(i,m,j) * fcancmxrow(i,m,j)
            mterm_mo_t(i,m) = mterm_mo_t(i,m) + mterm_mo(i,m,j) * fcancmxrow(i,m,j)
            smfuncveg_mo_t(i,m) = smfuncveg_mo_t(i,m) + smfuncveg_mo(i,m,j) * fcancmxrow(i,m,j)
            burnfrac_mo_t(i,m) = burnfrac_mo_t(i,m) + burnfrac_mo(i,m,j) * fcancmxrow(i,m,j)
            laimaxg_mo_t(i,m) = laimaxg_mo_t(i,m) + laimaxg_mo(i,m,j) * fcancmxrow(i,m,j)
            litrfall_mo_t(i,m) = litrfall_mo_t(i,m) + litrfallveg_mo(i,m,j) * fcancmxrow(i,m,j)
            humiftrs_mo_t(i,m) = humiftrs_mo_t(i,m) + humiftrsveg_mo(i,m,j) * fcancmxrow(i,m,j)
            vcmax0_mo_t(i,m) = vcmax0_mo_t(i,m) + vcmax0_mo(i,m,j) * fcancmxrow(i,m,j)
            tltrleaf_mo_t(i,m) = tltrleaf_mo_t(i,m) + tltrleaf_mo(i,m,j) * fcancmxrow(i,m,j)
            tltrstem_mo_t(i,m) = tltrstem_mo_t(i,m) + tltrstem_mo(i,m,j) * fcancmxrow(i,m,j)
            tltrroot_mo_t(i,m) = tltrroot_mo_t(i,m) + tltrroot_mo(i,m,j) * fcancmxrow(i,m,j)

            if (Ncycle_on) then
                bnf_free_mo_t(i,m) = bnf_free_mo_t(i,m) + bnf_free_mo(i,m,j) * fcancmxrow(i,m,j)
                bnf_ant_mo_t(i,m) = bnf_ant_mo_t(i,m) + bnf_ant_mo(i,m,j) * fcancmxrow(i,m,j)
                bnf_nat_mo_t(i,m) = bnf_nat_mo_t(i,m) + bnf_nat_mo(i,m,j) * fcancmxrow(i,m,j)
                bnf_tot_mo_t(i,m) = bnf_tot_mo_t(i,m) + bnf_tot_mo(i,m,j) * fcancmxrow(i,m,j)
                nstress_mo_t(i,m) = nstress_mo_t(i,m) + nstress_mo(i,m,j) * fcancmxrow(i,m,j)
                nitrif_mo_t(i,m) = nitrif_mo_t(i,m) + nitrif_mo(i,m,j) * fcancmxrow(i,m,j)
                no_nit_mo_t(i,m) = no_nit_mo_t(i,m) + no_nit_mo(i,m,j) * fcancmxrow(i,m,j)
                no_denit_mo_t(i,m) = no_denit_mo_t(i,m) + no_denit_mo(i,m,j) * fcancmxrow(i,m,j)
                no_nitdenit_mo_t(i,m) = no_nitdenit_mo_t(i,m) &
                                        + no_nitdenit_mo(i,m,j) * fcancmxrow(i,m,j)
                n2o_nit_mo_t(i,m) = n2o_nit_mo_t(i,m) + n2o_nit_mo(i,m,j) * fcancmxrow(i,m,j)
                n2o_denit_mo_t(i,m) = n2o_denit_mo_t(i,m) &
                                      + n2o_denit_mo(i,m,j) * fcancmxrow(i,m,j)
                n2o_nitdenit_mo_t(i,m) = n2o_nitdenit_mo_t(i,m) &
                                         + n2o_nitdenit_mo(i,m,j) * fcancmxrow(i,m,j)
                n2_denit_mo_t(i,m) = n2_denit_mo_t(i,m) + n2_denit_mo(i,m,j) * fcancmxrow(i,m,j)
                nvol_mo_t(i,m) = nvol_mo_t(i,m) + nvol_mo(i,m,j) * fcancmxrow(i,m,j)
                nleach_mo_t(i,m) = nleach_mo_t(i,m) + nleach_mo(i,m,j) * fcancmxrow(i,m,j)
                appl_fert_mo_t(i,m) = appl_fert_mo_t(i,m) &
                                      + appl_fert_mo(i,m,j) * fcancmxrow(i,m,j)
                ndep_nh4_mo_t(i,m) = ndep_nh4_mo_t(i,m) + ndep_nh4_mo(i,m,j) * fcancmxrow(i,m,j)
                ndep_no3_mo_t(i,m) = ndep_no3_mo_t(i,m) + ndep_no3_mo(i,m,j) * fcancmxrow(i,m,j)
                ndemand_wp_npp_mo_t(i,m) = ndemand_wp_npp_mo_t(i,m) + ndemand_wp_npp_mo(i,m,j) * fcancmxrow(i,m,j)
                nuptake_p_nh4_mo_t(i,m) = nuptake_p_nh4_mo_t(i,m) &
                                          + nuptake_p_nh4_mo(i,m,j) * fcancmxrow(i,m,j)
                nuptake_p_no3_mo_t(i,m) = nuptake_p_no3_mo_t(i,m) &
                                          + nuptake_p_no3_mo(i,m,j) * fcancmxrow(i,m,j)
                nuptake_a_actl_nh4_mo_t(i,m) = nuptake_a_actl_nh4_mo_t(i,m) &
                                               + nuptake_a_actl_nh4_mo(i,m,j) * fcancmxrow(i,m,j)
                nuptake_a_actl_no3_mo_t(i,m) = nuptake_a_actl_no3_mo_t(i,m) &
                                               + nuptake_a_actl_no3_mo(i,m,j) * fcancmxrow(i,m,j)
                nuptake_mo_t(i,m) = nuptake_mo_t(i,m) + nuptake_mo(i,m,j) * fcancmxrow(i,m,j)
                nalloc_l_mo_t(i,m) = nalloc_l_mo_t(i,m) + nalloc_l_mo(i,m,j) * fcancmxrow(i,m,j)
                nalloc_s_mo_t(i,m) = nalloc_s_mo_t(i,m) + nalloc_s_mo(i,m,j) * fcancmxrow(i,m,j)
                nalloc_r_mo_t(i,m) = nalloc_r_mo_t(i,m) + nalloc_r_mo(i,m,j) * fcancmxrow(i,m,j)
                nresorped_s_mo_t(i,m) = nresorped_s_mo_t(i,m) &
                                        + nresorped_s_mo(i,m,j) * fcancmxrow(i,m,j)
                nresorped_r_mo_t(i,m) = nresorped_r_mo_t(i,m) &
                                        + nresorped_r_mo(i,m,j) * fcancmxrow(i,m,j)
                nre_alloc_s2l_mo_t(i,m) = nre_alloc_s2l_mo_t(i,m) &
                                          + nre_alloc_s2l_mo(i,m,j) * fcancmxrow(i,m,j)
                nre_alloc_r2l_mo_t(i,m) = nre_alloc_r2l_mo_t(i,m) &
                                          + nre_alloc_r2l_mo(i,m,j) * fcancmxrow(i,m,j)
                nleafns2s_mo_t(i,m) = nleafns2s_mo_t(i,m) &
                                      + nleafns2s_mo(i,m,j) * fcancmxrow(i,m,j)
                nstemns2s_mo_t(i,m) = nstemns2s_mo_t(i,m) &
                                      + nstemns2s_mo(i,m,j) * fcancmxrow(i,m,j)
                nrootns2s_mo_t(i,m) = nrootns2s_mo_t(i,m) &
                                      + nrootns2s_mo(i,m,j) * fcancmxrow(i,m,j)
                nlitr_l_mo_t(i,m) = nlitr_l_mo_t(i,m) + nlitr_l_mo(i,m,j) * fcancmxrow(i,m,j)
                nlitr_s_mo_t(i,m) = nlitr_s_mo_t(i,m) + nlitr_s_mo(i,m,j) * fcancmxrow(i,m,j)
                nlitr_r_mo_t(i,m) = nlitr_r_mo_t(i,m) + nlitr_r_mo(i,m,j) * fcancmxrow(i,m,j)
                nlitr_mo_t(i,m) = nlitr_mo_t(i,m) + nlitr_mo(i,m,j) * fcancmxrow(i,m,j)
                gl2bl_grass_nflux_mo_t(i,m) = gl2bl_grass_nflux_mo_t(i,m) &
                                              + gl2bl_grass_nflux_mo(i,m,j) * fcancmxrow(i,m,j)
                nhumtrs_mo_t(i,m) = nhumtrs_mo_t(i,m) + nhumtrs_mo(i,m,j) * fcancmxrow(i,m,j)
                nmineral_litr_mo_t(i,m) = nmineral_litr_mo_t(i,m) &
                                          + nmineral_litr_mo(i,m,j) * fcancmxrow(i,m,j)
                nmineral_humus_mo_t(i,m) = nmineral_humus_mo_t(i,m) &
                                           + nmineral_humus_mo(i,m,j) * fcancmxrow(i,m,j)
                netnmineral_mo_t(i,m) = netnmineral_mo_t(i,m) &
                                           + netnmineral_mo(i,m,j) * fcancmxrow(i,m,j)
                nimmobil_nh4_mo_t(i,m) = nimmobil_nh4_mo_t(i,m) &
                                         + nimmobil_nh4_mo(i,m,j) * fcancmxrow(i,m,j)
                nimmobil_no3_mo_t(i,m) = nimmobil_no3_mo_t(i,m) &
                                         + nimmobil_no3_mo(i,m,j) * fcancmxrow(i,m,j)
                fNnetland_mo_t(i,m) = fNnetland_mo_t(i,m) &
                                      + fNnetland_mo(i,m,j) * fcancmxrow(i,m,j)
            end if !for Ncycle_on
            barefrac = barefrac - fcancmxrow(i,m,j)

          end do ! j
          if (Ncycle_on .and. sum(fcancmxrow(i,m,:)) > 0.0) nstress_mo_t(i,m) = nstress_mo_t(i,m)/sum(fcancmxrow(i,m,:))

          nep_mo_t(i,m) = nep_mo_t(i,m) + nep_mo(i,m,iccp1) * barefrac
          ! nbp_mo_t(i,m)=nbp_mo_t(i,m)+nbp_mo(i,m,iccp1)*barefrac
          hetrores_mo_t(i,m) = hetrores_mo_t(i,m) + hetrores_mo(i,m,iccp1) * barefrac
          humiftrs_mo_t(i,m) = humiftrs_mo_t(i,m) + humiftrsveg_mo(i,m,iccp1) * barefrac
          do k = 1,ignd
            litres_mo_t(i,m,k) = litres_mo_t(i,m,k) + litres_mo(i,m,iccp1,k) * barefrac
            soilcres_mo_t(i,m,k) = soilcres_mo_t(i,m,k) + soilcres_mo(i,m,iccp1,k) * barefrac
          end do

          if (Ncycle_on) then
            bnf_free_mo_t(i,m) = bnf_free_mo_t(i,m) + bnf_free_mo(i,m,iccp1) * barefrac
            bnf_tot_mo_t(i,m) = bnf_tot_mo_t(i,m) + bnf_tot_mo(i,m,iccp1) * barefrac
            nitrif_mo_t(i,m) = nitrif_mo_t(i,m) + nitrif_mo(i,m,iccp1) * barefrac
            no_nit_mo_t(i,m) = no_nit_mo_t(i,m) + no_nit_mo(i,m,iccp1) * barefrac
            no_denit_mo_t(i,m) = no_denit_mo_t(i,m) + no_denit_mo(i,m,iccp1) * barefrac
            no_nitdenit_mo_t(i,m) = no_nitdenit_mo_t(i,m) + no_nitdenit_mo(i,m,iccp1) * barefrac
            n2o_nit_mo_t(i,m) = n2o_nit_mo_t(i,m) + n2o_nit_mo(i,m,iccp1) * barefrac
            n2o_denit_mo_t(i,m) = n2o_denit_mo_t(i,m) + n2o_denit_mo(i,m,iccp1) * barefrac
            n2o_nitdenit_mo_t(i,m) = n2o_nitdenit_mo_t(i,m) &
                                     + n2o_nitdenit_mo(i,m,iccp1) * barefrac
            n2_denit_mo_t(i,m) = n2_denit_mo_t(i,m) + n2_denit_mo(i,m,iccp1) * barefrac
            nvol_mo_t(i,m) = nvol_mo_t(i,m) + nvol_mo(i,m,iccp1) * barefrac
            nleach_mo_t(i,m) = nleach_mo_t(i,m) + nleach_mo(i,m,iccp1) * barefrac
            appl_fert_mo_t(i,m) = appl_fert_mo_t(i,m) + appl_fert_mo(i,m,iccp1) * barefrac
            ndep_nh4_mo_t(i,m) = ndep_nh4_mo_t(i,m) + ndep_nh4_mo(i,m,iccp1) * barefrac
            ndep_no3_mo_t(i,m) = ndep_no3_mo_t(i,m) + ndep_no3_mo(i,m,iccp1) * barefrac
            nhumtrs_mo_t(i,m) = nhumtrs_mo_t(i,m) + nhumtrs_mo(i,m,iccp1) * barefrac
            nmineral_litr_mo_t(i,m) = nmineral_litr_mo_t(i,m) &
                                      + nmineral_litr_mo(i,m,iccp1) * barefrac
            nmineral_humus_mo_t(i,m) = nmineral_humus_mo_t(i,m) &
                                       + nmineral_humus_mo(i,m,iccp1)* barefrac
            netnmineral_mo_t(i,m) = netnmineral_mo_t(i,m) &
                                       + netnmineral_mo(i,m,iccp1)* barefrac
            nimmobil_nh4_mo_t(i,m) = nimmobil_nh4_mo_t(i,m) &
                                     + nimmobil_nh4_mo(i,m,iccp1) * barefrac
            nimmobil_no3_mo_t(i,m) = nimmobil_no3_mo_t(i,m) &
                                     + nimmobil_no3_mo(i,m,iccp1) * barefrac
            fNnetland_mo_t(i,m) = fNnetland_mo_t(i,m) + fNnetland_mo(i,m,iccp1) * barefrac
          end if ! Ncycle_on

          !> Find the monthly outputs at the per grid cell level from the outputs at the per tile level
          npp_mo_g(i) = npp_mo_g(i) + npp_mo_t(i,m) * FAREROT(i,m)
          gpp_mo_g(i) = gpp_mo_g(i) + gpp_mo_t(i,m) * FAREROT(i,m)
          leafns2s_mo_g(i) = leafns2s_mo_g(i) + leafns2s_mo_t(i,m) * FAREROT(i,m)
          stemns2s_mo_g(i) = stemns2s_mo_g(i) + stemns2s_mo_t(i,m) * FAREROT(i,m)
          rootns2s_mo_g(i) = rootns2s_mo_g(i) + rootns2s_mo_t(i,m) * FAREROT(i,m)
          re_alloc_s2l_mo_g(i) = re_alloc_s2l_mo_g(i) + re_alloc_s2l_mo_t(i,m) * FAREROT(i,m)
          re_alloc_r2l_mo_g(i) = re_alloc_r2l_mo_g(i) + re_alloc_r2l_mo_t(i,m) * FAREROT(i,m)
          re_alloc_sr2l_mo_g(i) = re_alloc_sr2l_mo_g(i) + re_alloc_sr2l_mo_t(i,m) * FAREROT(i,m)
          nep_mo_g(i) = nep_mo_g(i) + nep_mo_t(i,m) * FAREROT(i,m)
          nbp_mo_g(i) = nbp_mo_g(i) + nbp_mo_t(i,m) * FAREROT(i,m)
          hetrores_mo_g(i) = hetrores_mo_g(i) + hetrores_mo_t(i,m) * FAREROT(i,m)
          autores_mo_g(i) = autores_mo_g(i) + autores_mo_t(i,m) * FAREROT(i,m)
          soilres_mo_g(i) = soilres_mo_g(i) + soilres_mo_t(i,m) * FAREROT(i,m)
          do k = 1,ignd
            litres_mo_g(i,k) = litres_mo_g(i,k) + litres_mo_t(i,m,k) * FAREROT(i,m)
            soilcres_mo_g(i,k) = soilcres_mo_g(i,k) + soilcres_mo_t(i,m,k) * FAREROT(i,m)
          end do
          laimaxg_mo_g(i) = laimaxg_mo_g(i) + laimaxg_mo_t(i,m) * FAREROT(i,m)
          emit_co2_mo_g(i) = emit_co2_mo_g(i) + emit_co2_mo_t(i,m) * FAREROT(i,m)
          emit_co_mo_g(i) = emit_co_mo_g(i) + emit_co_mo_t(i,m) * FAREROT(i,m)
          emit_ch4_mo_g(i) = emit_ch4_mo_g(i) + emit_ch4_mo_t(i,m) * FAREROT(i,m)
          emit_nmhc_mo_g(i) = emit_nmhc_mo_g(i) + emit_nmhc_mo_t(i,m) * FAREROT(i,m)
          emit_h2_mo_g(i) = emit_h2_mo_g(i) + emit_h2_mo_t(i,m) * FAREROT(i,m)
          emit_nox_mo_g(i) = emit_nox_mo_g(i) + emit_nox_mo_t(i,m) * FAREROT(i,m)
          emit_n2o_mo_g(i) = emit_n2o_mo_g(i) + emit_n2o_mo_t(i,m) * FAREROT(i,m)
          emit_nh3_mo_g(i) = emit_nh3_mo_g(i) + emit_nh3_mo_t(i,m) * FAREROT(i,m)
          emit_pm25_mo_g(i) = emit_pm25_mo_g(i) + emit_pm25_mo_t(i,m) * FAREROT(i,m)
          emit_tpm_mo_g(i) = emit_tpm_mo_g(i) + emit_tpm_mo_t(i,m) * FAREROT(i,m)
          emit_tc_mo_g(i) = emit_tc_mo_g(i) + emit_tc_mo_t(i,m) * FAREROT(i,m)
          emit_oc_mo_g(i) = emit_oc_mo_g(i) + emit_oc_mo_t(i,m) * FAREROT(i,m)
          emit_bc_mo_g(i) = emit_bc_mo_g(i) + emit_bc_mo_t(i,m) * FAREROT(i,m)
          burnfrac_mo_g(i) = burnfrac_mo_g(i) + burnfrac_mo_t(i,m) * FAREROT(i,m)
          luc_emc_mo_g(i) = luc_emc_mo_g(i) + luc_emc_mo_t(i,m) * FAREROT(i,m)
          lucsocin_mo_g(i) = lucsocin_mo_g(i) + lucsocin_mo_t(i,m) * FAREROT(i,m)
          lucltrin_mo_g(i) = lucltrin_mo_g(i) + lucltrin_mo_t(i,m) * FAREROT(i,m)
          tileAge_mo_g(i) = tileAge_mo_g(i) + tileAge_mo_t(i,m) * FAREROT(i,m)
          luc_emcn_mo_g(i) = luc_emcn_mo_g(i) + luc_emcn_mo_t(i,m) * FAREROT(i,m)
          lucsocinn_mo_g(i) = lucsocinn_mo_g(i) + lucsocinn_mo_t(i,m) * FAREROT(i,m)
          lucltrinn_mo_g(i) = lucltrinn_mo_g(i) + lucltrinn_mo_t(i,m) * FAREROT(i,m)
          ch4WetSpec_mo_g(i) = ch4WetSpec_mo_g(i) + ch4WetSpec_mo_t(i,m) * FAREROT(i,m)
          wetfdyn_mo_g(i) = wetfdyn_mo_g(i) + wetfdyn_mo_t(i,m) * FAREROT(i,m)
          wetfpres_mo_g(i) = wetfpres_mo_g(i) + wetfpres_mo_t(i,m) * FAREROT(i,m)
          ch4WetDyn_mo_g(i) = ch4WetDyn_mo_g(i) + ch4WetDyn_mo_t(i,m) * FAREROT(i,m)
          ch4soills_mo_g(i) = ch4soills_mo_g(i) + ch4soills_mo_t(i,m) * FAREROT(i,m)
          smfuncveg_mo_g(i) = smfuncveg_mo_g(i) + smfuncveg_mo_t(i,m) * FAREROT(i,m)
          bterm_mo_g(i) = bterm_mo_g(i) + bterm_mo_t(i,m) * FAREROT(i,m)
          lterm_mo_g(i) = lterm_mo_g(i) + lterm_mo_t(i,m) * FAREROT(i,m)
          mterm_mo_g(i) = mterm_mo_g(i) + mterm_mo_t(i,m) * FAREROT(i,m)
          litrfall_mo_g(i) = litrfall_mo_g(i) + litrfall_mo_t(i,m) * FAREROT(i,m)
          humiftrs_mo_g(i) = humiftrs_mo_g(i) + humiftrs_mo_t(i,m) * FAREROT(i,m)
          vcmax0_mo_g(i) = vcmax0_mo_g(i) + vcmax0_mo_t(i,m) * FAREROT(i,m)
          tltrleaf_mo_g(i) = tltrleaf_mo_g(i) + tltrleaf_mo_t(i,m) * FAREROT(i,m)
          tltrstem_mo_g(i) = tltrstem_mo_g(i) + tltrstem_mo_t(i,m) * FAREROT(i,m)
          tltrroot_mo_g(i) = tltrroot_mo_g(i) + tltrroot_mo_t(i,m) * FAREROT(i,m)
          if (Ncycle_on) then
            bnf_tot_mo_g(i) = bnf_tot_mo_g(i) + bnf_tot_mo_t(i,m) * FAREROT(i,m)
            bnf_free_mo_g(i) = bnf_free_mo_g(i) + bnf_free_mo_t(i,m) * FAREROT(i,m)
            bnf_ant_mo_g(i) = bnf_ant_mo_g(i) + bnf_ant_mo_t(i,m) * FAREROT(i,m)
            bnf_nat_mo_g(i) = bnf_nat_mo_g(i) + bnf_nat_mo_t(i,m) * FAREROT(i,m)
            nstress_mo_g(i) = nstress_mo_g(i) + nstress_mo_t(i,m) * FAREROT(i,m)
            nitrif_mo_g(i) = nitrif_mo_g(i) + nitrif_mo_t(i,m) * FAREROT(i,m)
            no_nit_mo_g(i) = no_nit_mo_g(i) + no_nit_mo_t(i,m) * FAREROT(i,m)
            no_denit_mo_g(i) = no_denit_mo_g(i) + no_denit_mo_t(i,m) * FAREROT(i,m)
            no_nitdenit_mo_g(i) = no_nitdenit_mo_g(i) + no_nitdenit_mo_t(i,m) * FAREROT(i,m)
            n2o_nit_mo_g(i) = n2o_nit_mo_g(i) + n2o_nit_mo_t(i,m) * FAREROT(i,m)
            n2o_denit_mo_g(i) = n2o_denit_mo_g(i) + n2o_denit_mo_t(i,m) * FAREROT(i,m)
            n2o_nitdenit_mo_g(i) = n2o_nitdenit_mo_g(i) + n2o_nitdenit_mo_t(i,m)* FAREROT(i,m)
            n2_denit_mo_g(i) = n2_denit_mo_g(i) + n2_denit_mo_t(i,m) * FAREROT(i,m)
            nvol_mo_g(i) = nvol_mo_g(i) + nvol_mo_t(i,m) * FAREROT(i,m)
            nleach_mo_g(i) = nleach_mo_g(i) + nleach_mo_t(i,m) * FAREROT(i,m)
            appl_fert_mo_g(i) = appl_fert_mo_g(i) + appl_fert_mo_t(i,m) * FAREROT(i,m)
            ndep_nh4_mo_g(i) = ndep_nh4_mo_g(i) + ndep_nh4_mo_t(i,m) * FAREROT(i,m)
            ndep_no3_mo_g(i) = ndep_no3_mo_g(i) + ndep_no3_mo_t(i,m) * FAREROT(i,m)
            ndemand_wp_npp_mo_g(i) = ndemand_wp_npp_mo_g(i) &
                                     + ndemand_wp_npp_mo_t(i,m) * FAREROT(i,m)
            nuptake_p_nh4_mo_g(i) = nuptake_p_nh4_mo_g(i) &
                                    + nuptake_p_nh4_mo_t(i,m) * FAREROT(i,m)
            nuptake_p_no3_mo_g(i) = nuptake_p_no3_mo_g(i) &
                                    + nuptake_p_no3_mo_t(i,m) * FAREROT(i,m)
            nuptake_a_actl_nh4_mo_g(i) = nuptake_a_actl_nh4_mo_g(i) &
                                         + nuptake_a_actl_nh4_mo_t(i,m) * FAREROT(i,m)
            nuptake_a_actl_no3_mo_g(i) = nuptake_a_actl_no3_mo_g(i) &
                                         + nuptake_a_actl_no3_mo_t(i,m) * FAREROT(i,m)
            nuptake_mo_g(i) = nuptake_mo_g(i) + nuptake_mo_t(i,m) * FAREROT(i,m)
            nalloc_l_mo_g(i) = nalloc_l_mo_g(i) + nalloc_l_mo_t(i,m) * FAREROT(i,m)
            nalloc_s_mo_g(i) = nalloc_s_mo_g(i) + nalloc_s_mo_t(i,m) * FAREROT(i,m)
            nalloc_r_mo_g(i) = nalloc_r_mo_g(i) + nalloc_r_mo_t(i,m) * FAREROT(i,m)
            nresorped_s_mo_g(i) = nresorped_s_mo_g(i) + nresorped_s_mo_t(i,m) * FAREROT(i,m)
            nresorped_r_mo_g(i) = nresorped_r_mo_g(i) + nresorped_r_mo_t(i,m) * FAREROT(i,m)
            nre_alloc_s2l_mo_g(i) = nre_alloc_s2l_mo_g(i) + nre_alloc_s2l_mo_t(i,m) * FAREROT(i,m)
            nre_alloc_r2l_mo_g(i) = nre_alloc_r2l_mo_g(i) + nre_alloc_r2l_mo_t(i,m) * FAREROT(i,m)
            nleafns2s_mo_g(i) = nleafns2s_mo_g(i) + nleafns2s_mo_t(i,m) * FAREROT(i,m)
            nstemns2s_mo_g(i) = nstemns2s_mo_g(i) + nstemns2s_mo_t(i,m) * FAREROT(i,m)
            nrootns2s_mo_g(i) = nrootns2s_mo_g(i) + nrootns2s_mo_t(i,m) * FAREROT(i,m)
            nlitr_l_mo_g(i) = nlitr_l_mo_g(i) + nlitr_l_mo_t(i,m) * FAREROT(i,m)
            nlitr_s_mo_g(i) = nlitr_s_mo_g(i) + nlitr_s_mo_t(i,m) * FAREROT(i,m)
            nlitr_r_mo_g(i) = nlitr_r_mo_g(i) + nlitr_r_mo_t(i,m) * FAREROT(i,m)
            nlitr_mo_g(i) = nlitr_mo_g(i) + nlitr_mo_t(i,m) * FAREROT(i,m)
            gl2bl_grass_nflux_mo_g(i) = gl2bl_grass_nflux_mo_g(i) &
                                        + gl2bl_grass_nflux_mo_t(i,m)* FAREROT(i,m)
            nhumtrs_mo_g(i) = nhumtrs_mo_g(i) + nhumtrs_mo_t(i,m) * FAREROT(i,m)
            nmineral_litr_mo_g(i) = nmineral_litr_mo_g(i)+ nmineral_litr_mo_t(i,m)* FAREROT(i,m)
            nmineral_humus_mo_g(i) = nmineral_humus_mo_g(i)+ nmineral_humus_mo_t(i,m)* FAREROT(i,m)
            netnmineral_mo_g(i) = netnmineral_mo_g(i)+ netnmineral_mo_t(i,m)* FAREROT(i,m)
            nimmobil_nh4_mo_g(i) = nimmobil_nh4_mo_g(i) + nimmobil_nh4_mo_t(i,m) * FAREROT(i,m)
            nimmobil_no3_mo_g(i) = nimmobil_no3_mo_g(i) + nimmobil_no3_mo_t(i,m) * FAREROT(i,m)
            fNnetland_mo_g(i) = fNnetland_mo_g(i) + fNnetland_mo_t(i,m) * FAREROT(i,m)
          end if !for Ncycle_on

          fProductDecomp_mo_g(i) = fProductDecomp_mo_g(i) + fProductDecomp_mo_t(i,m) * FAREROT(i,m)

        end do !m
        if (Ncycle_on .and. sum(FAREROT(i,:)) > 0.0) nstress_mo_g(i) = nstress_mo_g(i)/sum(FAREROT(i,:))

        imonth = nt

        ! Prepare the timestamp for this month

        ! Transfer the timestamp (need in size 1 array)
        timeStamp(1) = consecDays

        call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo_g' ,timeStamp,'lai', [laimaxg_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo_g',timeStamp,'cVeg',[vgbiomas_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo_g'     ,timeStamp,'npp',[npp_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo_g'     ,timeStamp,'gpp',[gpp_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo_g'     ,timeStamp,'nep',[nep_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo_g'     ,timeStamp,'nbp',[nbp_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo_g',timeStamp,'rh',[hetrores_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo_g' ,timeStamp,'ra',[autores_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilres_mo_g' ,timeStamp,'rSoil',[soilres_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_mo_g' ,timeStamp,'fVegLitter',[litrfall_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrs_mo_g' ,timeStamp,'fLitterSoil',[humiftrs_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'leaflitr_mo_g' ,timeStamp,'fLeafLitter',[tltrleaf_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'stemlitr_mo_g' ,timeStamp,'fStemLitter',[tltrstem_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rootlitr_mo_g' ,timeStamp,'fRootLitter',[tltrroot_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_mo_g' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_mo_g' ,timeStamp,'wetlandFrac',[wetfdyn_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_mo_g' ,timeStamp,'soilCH4cons',[ch4soills_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_mo_g',timeStamp,'cLeaf',[gleafmas_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_NS_mo_g',timeStamp,'cLeaf_ns',[gleafmas_NS_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmass_mo_g',timeStamp,'cLeaf_s',[gleafmass_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_mo_g',timeStamp,'cStem',[stemmass_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_NS_mo_g',timeStamp,'cStem_ns',[stemmass_NS_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmasss_mo_g',timeStamp,'cStem_s',[stemmasss_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_mo_g',timeStamp,'cRoot',[rootmass_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_NS_mo_g',timeStamp,'cRoot_ns',[rootmass_NS_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmasss_mo_g',timeStamp,'cRoot_s',[rootmasss_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'leafns2s_mo_g',timeStamp,'leafns2s',[leafns2s_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'stemns2s_mo_g',timeStamp,'stemns2s',[stemns2s_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rootns2s_mo_g',timeStamp,'rootns2s',[rootns2s_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_s2l_mo_g',timeStamp,'realloc_s2l',[re_alloc_s2l_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_r2l_mo_g',timeStamp,'realloc_r2l',[re_alloc_r2l_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_sr2l_mo_g',timeStamp,'realloc_sr2l',[re_alloc_sr2l_mo_g(i)])

        call writeOutput1D(lonLocalIndex,latLocalIndex,'vcmax0_mo_g',timeStamp,'vcmax0',[vcmax0_mo_g(i)])

        ! Make the bulk litter and soil C pool and respiration temporary variables:
        bulkLitterCarbon_mo_g(1) = sum(litrmass_mo_g(i,:))
        bulkSoilCarbon_mo_g(1) = sum(soilcmas_mo_g(i,:))
        bulkLitterResp_mo_g(1) = sum(litres_mo_g(i,:))
        bulkSoilResp_mo_g(1) = sum(soilcres_mo_g(i,:))
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_g',timeStamp,'cLitter',[bulkLitterCarbon_mo_g])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_g',timeStamp,'cSoil',[bulkSoilCarbon_mo_g])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_g'  ,timeStamp,'rhLitter',[bulkLitterResp_mo_g])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_g',timeStamp,'rhSoil',[bulkSoilResp_mo_g])

        if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_mo_g' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfpres_mo_g' ,timeStamp,'wetlandFracPresc',[wetfpres_mo_g(i)])
        end if
        if (Ncycle_on) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_mo_g',timeStamp,'bnf_tot',[bnf_tot_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_mo_g',timeStamp,'bnf_free',[bnf_free_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_mo_g',timeStamp,'bnf_ant',[bnf_ant_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_mo_g',timeStamp,'bnf_nat',[bnf_nat_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_mo_g',timeStamp,'nh4',[nh4_mass_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_mo_g',timeStamp,'no3',[no3_mass_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_mo_g',timeStamp,'nLeaf',[ngleafmas_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_NS_mo_g',timeStamp,'nLeaf_ns',[ngleafmas_NS_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmass_mo_g',timeStamp,'nLeaf_s',[ngleafmass_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_mo_g',timeStamp,'nbLeaf',[nbleafmas_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_mo_g',timeStamp,'nStem',[nstemmass_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_NS_mo_g',timeStamp,'nStem_ns',[nstemmass_NS_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmasss_mo_g',timeStamp,'nStem_s',[nstemmasss_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_mo_g',timeStamp,'nRoot',[nrootmass_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_NS_mo_g',timeStamp,'nRoot_ns',[nrootmass_NS_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmasss_mo_g',timeStamp,'nRoot_s',[nrootmasss_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_mo_g',timeStamp,'nLitter',[nlitrmass_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_mo_g',timeStamp,'nSoil',[soilnmas_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_mo_g',timeStamp,'nitrif',[nitrif_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nit_mo_g',timeStamp,'no_nit',[no_nit_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_mo_g',timeStamp,'no_denit',[no_denit_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_mo_g',timeStamp,'no',[no_nitdenit_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nit_mo_g',timeStamp,'n2o_nit',[n2o_nit_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_mo_g',timeStamp,'n2o_denit',[n2o_denit_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_mo_g',timeStamp,'n2o',[n2o_nitdenit_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_mo_g',timeStamp,'n2',[n2_denit_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_mo_g',timeStamp,'nvol',[nvol_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_mo_g',timeStamp,'nleach',[nleach_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'appl_fert_mo_g',timeStamp,'nfer_nh4',[appl_fert_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_nh4_mo_g',timeStamp,'ndep_nh4',[ndep_nh4_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_no3_mo_g',timeStamp,'ndep_no3',[ndep_no3_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ndemand_wp_npp_mo_g',timeStamp,'ndemand_npp',[ndemand_wp_npp_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_nh4_mo_g',timeStamp,'nuptake_p_nh4',[nuptake_p_nh4_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_no3_mo_g',timeStamp,'nuptake_p_no3',[nuptake_p_no3_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_nh4_mo_g',timeStamp,'nuptake_a_actl_nh4',[nuptake_a_actl_nh4_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_no3_mo_g',timeStamp,'nuptake_a_actl_no3',[nuptake_a_actl_no3_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_mo_g',timeStamp,'nuptake',[nuptake_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_l_mo_g' ,timeStamp,'nalloc_l', [nalloc_l_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_s_mo_g' ,timeStamp,'nalloc_s', [nalloc_s_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_r_mo_g' ,timeStamp,'nalloc_r', [nalloc_r_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_mo_g' ,timeStamp,'nresorped_s', [nresorped_s_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_mo_g' ,timeStamp,'nresorped_r', [nresorped_r_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_s2l_mo_g' ,timeStamp,'nre_alloc_s2l', [nre_alloc_s2l_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_r2l_mo_g' ,timeStamp,'nre_alloc_r2l', [nre_alloc_r2l_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nleafns2s_mo_g' ,timeStamp,'nleafns2s', [nleafns2s_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemns2s_mo_g' ,timeStamp,'nstemns2s', [nstemns2s_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootns2s_mo_g' ,timeStamp,'nrootns2s', [nrootns2s_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_mo_g' ,timeStamp,'nlitr_l', [nlitr_l_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_s_mo_g' ,timeStamp,'nlitr_s', [nlitr_s_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_r_mo_g' ,timeStamp,'nlitr_r', [nlitr_r_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_mo_g' ,timeStamp,'nlitr', [nlitr_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gl2bl_grass_nflux_mo_g' ,timeStamp,'gl2bl_grass_nflux', [gl2bl_grass_nflux_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_mo_g' ,timeStamp,'c2n_l', [c2n_l_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_s_mo_g' ,timeStamp,'c2n_s', [c2n_s_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_r_mo_g' ,timeStamp,'c2n_r', [c2n_r_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_wp_mo_g' ,timeStamp,'c2n_wp', [c2n_wp_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_litr_mo_g' ,timeStamp,'c2n_litr', [c2n_litr_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_humus_mo_g' ,timeStamp,'c2n_humus', [c2n_humus_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nhumtrs_mo_g',timeStamp,'nhumtrs',[nhumtrs_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_litr_mo_g',timeStamp,'nmineral_litr',[nmineral_litr_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_humus_mo_g',timeStamp,'nmineral_humus',[nmineral_humus_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'netnmineral_mo_g',timeStamp,'netnmineral',[netnmineral_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_nh4_mo_g',timeStamp,'nimmobil_nh4',[nimmobil_nh4_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_no3_mo_g',timeStamp,'nimmobil_no3',[nimmobil_no3_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_mo_g',timeStamp,'nVeg',[nvgbiomas_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nstress_mo_g',timeStamp,'nstress',[nstress_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fNnetland_mo_g',timeStamp,'fNnetland',[fNnetland_mo_g(i)])
        end if ! Ncycle_on

        do m = 1,nmtest
          sumfare = 0.0
          do j = 1,icc
            sumfare = sumfare + fcancmxrow(i,m,j)
          end do ! j
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_mo_g' ,timeStamp,'landCoverFrac',[fcancmxrow(i,m,1:icc),1. - sumfare])
        end do ! m

        if (dofire .or. prescribedFire .or. lnduseon) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_mo_g' ,timeStamp,'fFireCH4',[emit_ch4_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co_mo_g' ,timeStamp,'fFireCO',[emit_co_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo_g' ,timeStamp,'fFire',[emit_co2_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_bc_mo_g' ,timeStamp,'fFireBC',[emit_bc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_oc_mo_g' ,timeStamp,'fFireOC',[emit_oc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_nox_mo_g' ,timeStamp,'fFireNOX',[emit_nox_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_n2o_mo_g' ,timeStamp,'fFireN2O',[emit_n2o_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_nh3_mo_g' ,timeStamp,'fFireNH3',[emit_nh3_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_nmhc_mo_g' ,timeStamp,'fFireNMHC',[emit_nmhc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo_g' ,timeStamp,'burntFractionAll',[burnfrac_mo_g(i)])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'bterm_mo_g' ,timeStamp,'bterm',[bterm_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'mterm_mo_g' ,timeStamp,'mterm',[mterm_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'lterm_mo_g' ,timeStamp,'lterm',[lterm_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'smfuncveg_mo_g' ,timeStamp,'smfuncveg',[smfuncveg_mo_g(i)])
        end if

        if(dynamicTilingOn .or. trackTileAge) call writeOutput1D(lonLocalIndex,latLocalIndex,'tileAge_mo_g' ,timeStamp,'tileAge',[tileAge_mo_g(i)])

        if (lnduseon .or. timberHarvest) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_mo_g' ,timeStamp,'fDeforestToAtmos',[luc_emc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_mo_g' ,timeStamp,'fDeforestToLitter',[lucltrin_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_mo_g' ,timeStamp,'fDeforestToSoil',[lucsocin_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_mo_g' ,timeStamp,'fDeforestTotal', &
                             [lucsocin_mo_g(i) + lucltrin_mo_g(i) + luc_emc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fProductDecomp_mo_g' ,timeStamp,'fProductDecomp',[fProductDecomp_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'cProduct_mo_g' ,timeStamp,'cProduct',[cProduct_mo_g(i)])
          if (Ncycle_on) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nProduct_mo_g' ,timeStamp,'nProduct',[nProduct_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emcn_mo_g' ,timeStamp,'fDeforestToAtmosN',[luc_emcn_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrinn_mo_g' ,timeStamp,'fDeforestToLitterN',[lucltrinn_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocinn_mo_g' ,timeStamp,'fDeforestToSoilN',[lucsocinn_mo_g(i)])
          end if
        end if

        if (PFTCompetition) then
          m = 1
          pftExist = 0.0
          do j = 1,icc
            if (pftexistrow(i,1,j)) pftExist(j) = 1.0
          end do
          call writeOutput1D(lonLocalIndex,latLocalIndex,'pftexistrow_mo_g' ,timeStamp,'landCoverExist',[pftExist]) ! flag only set up for one tile !
          if (doperpftoutput) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'cc_mo' ,timeStamp,'cc',[cc_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'mm_mo' ,timeStamp,'mm',[mm_mo(i,m,:)])
          end if
        end if

        if (doperpftoutput) then
          if (nmtest > 1) then
            print * ,'Per PFT and per tile outputs together not implemented yet'
          else
            m = 1
            call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo' ,timeStamp,'lai', [laimaxg_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_mo',timeStamp,'cLeaf',[gleafmas_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_mo',timeStamp,'cStem',[stemmass_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_mo',timeStamp,'cRoot',[rootmass_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootdpth_mo',timeStamp,'rootdpth',[rootdpth_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo',timeStamp,'cVeg',[vgbiomas_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo'     ,timeStamp,'npp',[npp_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo'     ,timeStamp,'gpp',[gpp_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo'     ,timeStamp,'nep',[nep_mo(i,m,:)])
            ! NOTE: This NBP does not include LUC product pool contributions since they are
            ! not per PFT but rather per tile
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo'     ,timeStamp,'nbp',[nbp_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo',timeStamp,'rh',[hetrores_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo' ,timeStamp,'ra',[autores_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilres_mo' ,timeStamp,'rSoil',[soilres_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfallveg_mo' ,timeStamp,'fVegLitter',[litrfallveg_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrsveg_mo' ,timeStamp,'fLitterSoil',[humiftrsveg_mo(i,m,1:iccp1)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'leaflitr_mo' ,timeStamp,'fLeafLitter',[tltrleaf_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'stemlitr_mo' ,timeStamp,'fStemLitter',[tltrstem_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootlitr_mo' ,timeStamp,'fRootLitter',[tltrroot_mo(i,m,:)])

            ! Make the bulk litter and soil C pool and respiration temporary variables:
            do k = 1,iccp1
              bulkLitterCarbon_mo(k) = sum(litrmass_mo(i,m,k,:))
              bulkSoilCarbon_mo(k) = sum(soilcmas_mo(i,m,k,:))
              bulkLitterResp_mo(k) = sum(litres_mo(i,m,k,:))
              bulkSoilResp_mo(k) = sum(soilcres_mo(i,m,k,:))
            end do
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo',timeStamp,'cLitter',[bulkLitterCarbon_mo])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo',timeStamp,'cSoil',[bulkSoilCarbon_mo])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo'  ,timeStamp,'rhLitter',[bulkLitterResp_mo])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo',timeStamp,'rhSoil',[bulkSoilResp_mo])

            call writeOutput1D(lonLocalIndex,latLocalIndex,'vcmax0_mo',timeStamp,'vcmax0',[vcmax0_mo(i,m,:)])

            if (Ncycle_on) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_mo',timeStamp,'nLeaf',[ngleafmas_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_mo',timeStamp,'nbLeaf',[nbleafmas_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_mo',timeStamp,'nStem',[nstemmass_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_mo',timeStamp,'nRoot',[nrootmass_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_mo',timeStamp,'nLitter',[nlitrmass_mo(i,m,1:iccp1)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_mo',timeStamp,'nSoil',[soilnmas_mo(i,m,1:iccp1)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_mo',timeStamp,'nh4',[nh4_mass_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_mo',timeStamp,'no3',[no3_mass_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_mo',timeStamp,'bnf_tot',[bnf_tot_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_mo',timeStamp,'bnf_nat',[bnf_nat_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_mo',timeStamp,'bnf_ant',[bnf_ant_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_mo',timeStamp,'bnf_free',[bnf_free_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_mo',timeStamp,'nVeg',[nvgbiomas_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nstress_mo',timeStamp,'nstress',[nstress_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_mo',timeStamp,'n2o',[n2o_nitdenit_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_mo',timeStamp,'no',[no_nitdenit_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_mo',timeStamp,'nitrif',[nitrif_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'netnmineral_mo',timeStamp,'netnmineral',[netnmineral_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_mo',timeStamp,'nuptake',[nuptake_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_mo',timeStamp,'c2n_l',[c2n_l_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_mo',timeStamp,'nlitr_l',[nlitr_l_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_mo',timeStamp,'no_denit',[no_denit_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_mo',timeStamp,'n2o_denit',[n2o_denit_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_mo',timeStamp,'n2',[n2_denit_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_mo',timeStamp,'nvol',[nvol_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_mo',timeStamp,'nleach',[nleach_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'appl_fert_mo',timeStamp,'nfer_nh4',[appl_fert_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_mo',timeStamp,'nresorped_s',[nresorped_s_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_mo',timeStamp,'nresorped_r',[nresorped_r_mo(i,m,:)])
            end if ! Ncycle_on

            if (dofire .or. lnduseon .or. timberHarvest .or. prescribedFire) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_mo' ,timeStamp,'fFireCH4',[emit_ch4_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo' ,timeStamp,'fFire',[emit_co2_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo' ,timeStamp,'burntFractionAll',[burnfrac_mo(i,m,:)])
            end if

            if (useTracer > 0) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerGLeafMass' ,timeStamp,'cLeafTracer',[tracerGLeafMassrot(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerStemMass' ,timeStamp,'cStemTracer',[tracerStemMassrot(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerRootMass' ,timeStamp,'cRootTracer',[tracerRootMassrot(i,m,:)])
              ! call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerLitrMass' ,timeStamp,'cLitterTracer',[tracerLitrMassrot(i,m,:)])
              ! call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerSoilCMass' ,timeStamp,'cSoilTracer',[tracerSoilCMassrot(i,m,:)])
            end if

          end if
        end if

        if (dopertileoutput) then
          if (nmtest > 1) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo_t' ,timeStamp,'lai', [laimaxg_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo_t',timeStamp,'cVeg',[vgbiomas_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_mo_t',timeStamp,'cLeaf',[gleafmas_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_NS_mo_t',timeStamp,'cLeaf_ns',[gleafmas_NS_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmass_mo_t',timeStamp,'cLeaf_s',[gleafmass_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_mo_t',timeStamp,'cStem',[stemmass_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_NS_mo_t',timeStamp,'cStem_ns',[stemmass_NS_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmasss_mo_t',timeStamp,'cStem_s',[stemmasss_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_mo_t',timeStamp,'cRoot',[rootmass_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_NS_mo_t',timeStamp,'cRoot_ns',[rootmass_NS_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmasss_mo_t',timeStamp,'cRoot_s',[rootmasss_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'leafns2s_mo_t',timeStamp,'leafns2s',[leafns2s_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'stemns2s_mo_t',timeStamp,'stemns2s',[stemns2s_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootns2s_mo_t',timeStamp,'rootns2s',[rootns2s_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_s2l_mo_t' ,timeStamp,'realloc_s2l',[re_alloc_s2l_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_r2l_mo_t' ,timeStamp,'realloc_r2l',[re_alloc_r2l_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_sr2l_mo_t',timeStamp,'realloc_sr2l',[re_alloc_sr2l_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo_t'     ,timeStamp,'npp',[npp_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo_t'     ,timeStamp,'gpp',[gpp_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo_t'     ,timeStamp,'nep',[nep_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo_t'     ,timeStamp,'nbp',[nbp_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo_t',timeStamp,'rh',[hetrores_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo_t' ,timeStamp,'ra',[autores_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilres_mo_t' ,timeStamp,'rSoil',[soilres_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_mo_t' ,timeStamp,'fVegLitter',[litrfall_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrs_mo_t' ,timeStamp,'fLitterSoil',[humiftrs_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'leaflitr_mo_t' ,timeStamp,'fLeafLitter',[tltrleaf_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'stemlitr_mo_t' ,timeStamp,'fStemLitter',[tltrstem_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rootlitr_mo_t' ,timeStamp,'fRootLitter',[tltrroot_mo_t(i,:)])
            !call writeOutput1D(lonLocalIndex,latLocalIndex,'leafmass_mo_t',timeStamp,'cLeaf',[leafmass_mo_t(i,:)])

            ! Make the bulk litter and soil C pool and respiration temporary variables:

            do m = 1,nmtest
              bulkLitterCarbon_mo_t(m) = sum(litrmass_mo_t(i,m,:))
              bulkSoilCarbon_mo_t(m) = sum(soilcmas_mo_t(i,m,:))
              bulkLitterResp_mo_t(m) = sum(litres_mo_t(i,m,:))
              bulkSoilResp_mo_t(m) = sum(soilcres_mo_t(i,m,:))
            end do

            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_t',timeStamp,'cLitter',[bulkLitterCarbon_mo_t])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_t',timeStamp,'cSoil',[bulkSoilCarbon_mo_t])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_t'  ,timeStamp,'rhLitter',[bulkLitterResp_mo_t])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_t',timeStamp,'rhSoil',[bulkSoilResp_mo_t])

            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_mo_t' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_mo_t' ,timeStamp,'wetlandFrac',[wetfdyn_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_mo_t' ,timeStamp,'soilCH4cons',[ch4soills_mo_t(i,:)])

            if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfpres_mo_t' ,timeStamp,'wetlandFracPresc',[wetfpres_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_mo_t' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_mo_t(i,:)])
            end if

            call writeOutput1D(lonLocalIndex,latLocalIndex,'vcmax0_mo_t',timeStamp,'vcmax0',[vcmax0_mo_t(i,:)])
            if (Ncycle_on) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_mo_t',timeStamp,'nLeaf',[ngleafmas_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_NS_mo_t',timeStamp,'nLeaf_ns',[ngleafmas_NS_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmass_mo_t',timeStamp,'nLeaf_s',[ngleafmass_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_mo_t',timeStamp,'nbLeaf',[nbleafmas_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_mo_t',timeStamp,'nStem',[nstemmass_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_NS_mo_t',timeStamp,'nStem_ns',[nstemmass_NS_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmasss_mo_t',timeStamp,'nStem_s',[nstemmasss_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_mo_t',timeStamp,'nRoot',[nrootmass_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_NS_mo_t',timeStamp,'nRoot_ns',[nrootmass_NS_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmasss_mo_t',timeStamp,'nRoot_s',[nrootmasss_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_mo_t' ,timeStamp,'c2n_l', [c2n_l_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_s_mo_t' ,timeStamp,'c2n_s', [c2n_s_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_r_mo_t' ,timeStamp,'c2n_r', [c2n_r_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_wp_mo_t' ,timeStamp,'c2n_wp', [c2n_wp_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_litr_mo_t' ,timeStamp,'c2n_litr', [c2n_litr_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_humus_mo_t' ,timeStamp,'c2n_humus', [c2n_humus_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_mo_t',timeStamp,'nLitter',[nlitrmass_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_mo_t',timeStamp,'nSoil',[soilnmas_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_mo_t',timeStamp,'bnf_tot',[bnf_tot_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_mo_t',timeStamp,'bnf_free',[bnf_free_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_mo_t',timeStamp,'bnf_ant',[bnf_ant_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_mo_t',timeStamp,'bnf_nat',[bnf_nat_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_mo_t',timeStamp,'nh4',[nh4_mass_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_mo_t',timeStamp,'no3',[no3_mass_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_mo_t',timeStamp,'nitrif',[nitrif_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nit_mo_t',timeStamp,'no_nit',[no_nit_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_mo_t',timeStamp,'no_denit',[no_denit_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_mo_t',timeStamp,'no',[no_nitdenit_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nit_mo_t',timeStamp,'n2o_nit',[n2o_nit_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_mo_t',timeStamp,'n2o_denit',[n2o_denit_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_mo_t',timeStamp,'n2o',[n2o_nitdenit_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_mo_t',timeStamp,'n2',[n2_denit_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_mo_t',timeStamp,'nvol',[nvol_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_mo_t',timeStamp,'nleach',[nleach_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'appl_fert_mo_t',timeStamp,'nfer_nh4',[appl_fert_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_nh4_mo_t',timeStamp,'ndep_nh4',[ndep_nh4_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_no3_mo_t',timeStamp,'ndep_no3',[ndep_no3_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ndemand_wp_npp_mo_t',timeStamp,'ndemand_npp',[ndemand_wp_npp_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_nh4_mo_t',timeStamp,'nuptake_p_nh4',[nuptake_p_nh4_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_no3_mo_t',timeStamp,'nuptake_p_no3',[nuptake_p_no3_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_nh4_mo_t',timeStamp,'nuptake_a_actl_nh4',[nuptake_a_actl_nh4_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_no3_mo_t',timeStamp,'nuptake_a_actl_no3',[nuptake_a_actl_no3_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_mo_t',timeStamp,'nuptake',[nuptake_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_l_mo_t' ,timeStamp,'nalloc_l', [nalloc_l_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_s_mo_t' ,timeStamp,'nalloc_s', [nalloc_s_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_r_mo_t' ,timeStamp,'nalloc_r', [nalloc_r_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_mo_t' ,timeStamp,'nresorped_s', [nresorped_s_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_mo_t' ,timeStamp,'nresorped_r', [nresorped_r_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_s2l_mo_t' ,timeStamp,'nre_alloc_s2l', [nre_alloc_s2l_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_r2l_mo_t' ,timeStamp,'nre_alloc_r2l', [nre_alloc_r2l_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nleafns2s_mo_t' ,timeStamp,'nleafns2s', [nleafns2s_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemns2s_mo_t' ,timeStamp,'nstemns2s', [nstemns2s_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootns2s_mo_t' ,timeStamp,'nrootns2s', [nrootns2s_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_mo_t' ,timeStamp,'nlitr_l', [nlitr_l_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_s_mo_t' ,timeStamp,'nlitr_s', [nlitr_s_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_r_mo_t' ,timeStamp,'nlitr_r', [nlitr_r_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_mo_t' ,timeStamp,'nlitr', [nlitr_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'gl2bl_grass_nflux_mo_t' ,timeStamp,'gl2bl_grass_nflux', [gl2bl_grass_nflux_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nhumtrs_mo_t' ,timeStamp,'nhumtrs', [nhumtrs_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_litr_mo_t',timeStamp,'nmineral_litr',[nmineral_litr_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_humus_mo_t',timeStamp,'nmineral_humus',[nmineral_humus_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'netnmineral_mo_t',timeStamp,'netnmineral',[netnmineral_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_nh4_mo_t',timeStamp,'nimmobil_nh4',[nimmobil_nh4_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_no3_mo_t',timeStamp,'nimmobil_no3',[nimmobil_no3_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'fNnetland_mo_t',timeStamp,'fNnetland',[fNnetland_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_mo_t',timeStamp,'nVeg',[nvgbiomas_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'nstress_mo_t',timeStamp,'nstress',[nstress_mo_t(i,:)])
             end if ! Ncycle
             
            if (dofire .or. prescribedFire .or. timberHarvest .or. lnduseon) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_mo_t' ,timeStamp,'fFireCH4',[emit_ch4_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo_t' ,timeStamp,'fFire',[emit_co2_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo_t' ,timeStamp,'burntFractionAll',[burnfrac_mo_t(i,:)])
            end if

            if(dynamicTilingOn) call writeOutput1D(lonLocalIndex,latLocalIndex,'FAREROT_mo_t' ,timeStamp,'FARE',[FAREROT(i,:)]) !the final value of FARE for each month can be written directly (FARE shouldnt change within months and agreegating it across tiles isn't informative)
            if(dynamicTilingOn .or. trackTileAge) call writeOutput1D(lonLocalIndex,latLocalIndex,'tileAge_mo_t' ,timeStamp,'tileAge',[tileAge_mo_t(i,:)])

            if (lnduseon .or. timberHarvest) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_mo_t' ,timeStamp,'fDeforestToAtmos',[luc_emc_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_mo_t' ,timeStamp,'fDeforestToLitter',[lucltrin_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_mo_t' ,timeStamp,'fDeforestToSoil',[lucsocin_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_mo_t' ,timeStamp,'fDeforestTotal', &
                                 [lucsocin_mo_t(i,:) + lucltrin_mo_t(i,:) + luc_emc_mo_t(i,:)])
              if (Ncycle_on) then
                call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emcn_mo_t' ,timeStamp,'fDeforestToAtmosN',[luc_emcn_mo_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrinn_mo_t' ,timeStamp,'fDeforestToLitterN',[lucltrinn_mo_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocinn_mo_t' ,timeStamp,'fDeforestToSoilN',[lucsocinn_mo_t(i,:)])  
              end if
            end if
          end if
        end if

        !> Reset all end of month accumulated arrays
        call resetMonthEnd(nltest ,nmtest)

      end if ! end of month
    end do ! loop 865 ! nmon

    end associate
  end subroutine ctem_monthly_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_ctem_annual_aw
  !> @{
  !> Accumulate and write out the annual biogeochemical (CTEM) outputs. These are kept in pointer structures as
  !! this subroutine is called daily and we increment the daily values to produce annual values.
  !! The pointer to the annual data structures (in ctemStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine ctem_annual_aw (lonLocalIndex, latLocalIndex, iday, realyr, nltest, nmtest, lastDOY)

    use classStateVars, only : class_rot
    use ctemStateVars,  only : ctem_tile_yr, vrot, ctem_grd_yr, c_switch, ctem_yr, &
                               resetYearEnd, vgat
    use classicParams,  only : icc, iccp1, iccp2, ignd, nlat, nmos, convertkgN, convertg2kg
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    integer, intent(in) :: realyr
    integer, intent(in) :: lastDOY

    ! local
    integer :: i, m, j, nt, k
    real :: barefrac
    real :: sumfare
    real, dimension(1) :: timeStamp
    real, dimension(icc) :: pftExist
    real, dimension(icc) :: fcancmxNoSeed

    real, dimension(1) :: bulkLitterCarbon_yr_g !< Temporary variable used to produce the bulk soil litter carbon quantity for output at the grid annual scale \f$[kg C m^{-2}]\f$
    real, dimension(1) :: bulkSoilCarbon_yr_g   !< Temporary variable used to produce the bulk soil carbon quantity for output at the grid annual scale \f$[kg C m^{-2}]\f$
    real, dimension(1) :: bulkLitterResp_yr_g   !< Temporary variable used to produce the bulk litter respiration quantity for output at the grid annual scale \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(1) :: bulkSoilResp_yr_g     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output at the grid annual scale \f$[kg C m^{-2} s^{-1}]\f$

    real, dimension(iccp1) :: bulkLitterCarbon_yr !< Temporary variable used to produce the bulk soil litter carbon quantity for output at the annual scale \f$[kg C m^{-2}]\f$
    real, dimension(iccp1) :: bulkSoilCarbon_yr   !< Temporary variable used to produce the bulk soil carbon quantity for output at the annual scale  \f$[kg C m^{-2}]\f$
    real, dimension(iccp1) :: bulkLitterResp_yr   !< Temporary variable used to produce the bulk litter respiration quantity for output at the annual scale  \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(iccp1) :: bulkSoilResp_yr     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output at the annual scale  \f$[kg C m^{-2} s^{-1}]\f$

    real, dimension(nmtest) :: bulkLitterCarbon_yr_t !< Temporary variable used to produce the bulk soil litter carbon quantity for output at the annual tile scale \f$[kg C m^{-2}]\f$
    real, dimension(nmtest) :: bulkSoilCarbon_yr_t   !< Temporary variable used to produce the bulk soil carbon quantity for output at the annual tile scale  \f$[kg C m^{-2}]\f$
    real, dimension(nmtest) :: bulkLitterResp_yr_t   !< Temporary variable used to produce the bulk litter respiration quantity for output at the annual tile scale  \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(nmtest) :: bulkSoilResp_yr_t     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output at the annual tile scale  \f$[kg C m^{-2} s^{-1}]\f$

    real :: oneOverDPY

    real :: coverLow, coverHigh

    ! Associate names with variables defined in derived types.

    associate( &
    dofire                => c_switch%dofire,                      & !< logical:
    prescribedFire        => c_switch%prescribedFire,         & !< logical:
    lnduseon              => c_switch%lnduseon,                    & !< logical:
    PFTCompetition        => c_switch%PFTCompetition,              & !< logical:
    doperpftoutput        => c_switch%doperpftoutput,              & !< logical:
    dopertileoutput       => c_switch%dopertileoutput,             & !< logical:
    transientOBSWETF      => c_switch%transientOBSWETF,            & !< logical:
    fixedYearOBSWETF      => c_switch%fixedYearOBSWETF,            & !< integer:
    Ncycle_on             => c_switch%Ncycle_on,                   & !< logical:
    dynamicTilingOn         => c_switch%dynamicTilingOn,               & !< logical:
    trackTileAge         => c_switch%trackTileAge,               & !< logical:
    timberHarvest         => c_switch%timberHarvest,               & !< logical:

    FAREROT => class_rot%FAREROT,                                  & !< real, dimension(:,:) : Fractional coverage of mosaic tile on modelled area

    laimaxg_yr          =>ctem_yr%laimaxg_yr,                      & !< real, dimension(:,:,:) :
    c2n_l_yr            =>ctem_yr%c2n_l_yr,                        & !< real, dimension(:,:,:) :
    c2n_s_yr            =>ctem_yr%c2n_s_yr,                        & !< real, dimension(:,:,:) :
    c2n_r_yr            =>ctem_yr%c2n_r_yr,                        & !< real, dimension(:,:,:) :
    c2n_wp_yr           =>ctem_yr%c2n_wp_yr,                       & !< real, dimension(:,:,:) :
    c2n_litr_yr         =>ctem_yr%c2n_litr_yr,                     & !< real, dimension(:,:,:) :
    c2n_humus_yr        =>ctem_yr%c2n_humus_yr,                    & !< real, dimension(:,:,:) :
    gleafmas_yr         =>ctem_yr%gleafmas_yr,                     & !< real, dimension(:,:,:) :
    gleafmas_NS_yr       =>ctem_yr%gleafmas_NS_yr,                 & !< real, dimension(:,:,:) :
    gleafmass_yr        =>ctem_yr%gleafmass_yr,                    & !< real, dimension(:,:,:) :
    bleafmas_yr         =>ctem_yr%bleafmas_yr,                     & !< real, dimension(:,:,:) :
    stemmass_yr         =>ctem_yr%stemmass_yr,                     & !< real, dimension(:,:,:) :
    stemmass_NS_yr       =>ctem_yr%stemmass_NS_yr,                 & !< real, dimension(:,:,:) :
    stemmasss_yr        =>ctem_yr%stemmasss_yr,                    & !< real, dimension(:,:,:) :
    rootmass_yr         =>ctem_yr%rootmass_yr,                     & !< real, dimension(:,:,:) :
    rootmass_NS_yr       =>ctem_yr%rootmass_NS_yr,                 & !< real, dimension(:,:,:) :
    rootmasss_yr        =>ctem_yr%rootmasss_yr,                    & !< real, dimension(:,:,:) :
    cc_yr               =>ctem_yr%cc_yr,                           & !< real, dimension(:,:) :
    mm_yr               =>ctem_yr%mm_yr,                           & !< real, dimension(:,:) :
    npp_yr              =>ctem_yr%npp_yr,                          & !< real, dimension(:,:,:) :
    gpp_yr              =>ctem_yr%gpp_yr,                          & !< real, dimension(:,:,:) :
    vcmax0_yr           =>ctem_yr%vcmax0_yr,                       & !< real, dimension(:,:,:) :
    leafns2s_yr         =>ctem_yr%leafns2s_yr,                     & !< real, dimension(:,:,:) :
    stemns2s_yr         =>ctem_yr%stemns2s_yr,                     & !< real, dimension(:,:,:) :
    rootns2s_yr         =>ctem_yr%rootns2s_yr,                     & !< real, dimension(:,:,:) :
    re_alloc_s2l_yr     =>ctem_yr%re_alloc_s2l_yr,                 & !< real, dimension(:,:,:) :
    re_alloc_r2l_yr     =>ctem_yr%re_alloc_r2l_yr,                 & !< real, dimension(:,:,:) :
    re_alloc_sr2l_yr    =>ctem_yr%re_alloc_sr2l_yr,                & !< real, dimension(:,:,:) :
    vgbiomas_yr         =>ctem_yr%vgbiomas_yr,                     & !< real, dimension(:,:,:) :
    autores_yr          =>ctem_yr%autores_yr,                      & !< real, dimension(:,:,:) :
    rmrveg_yr           =>ctem_yr%rmrveg_yr,                       & !< real, dimension(:,:,:) :
    totcmass_yr         =>ctem_yr%totcmass_yr,                     & !< real, dimension(:,:,:) :
    litrmass_yr         =>ctem_yr%litrmass_yr,                     & !< real, dimension(:,:,:,:) :
    soilcmas_yr         =>ctem_yr%soilcmas_yr,                     & !< real, dimension(:,:,:,:) :
    nh4_mass_yr         =>ctem_yr%nh4_mass_yr,                     & !< real, dimension(:,:,:) :
    no3_mass_yr         =>ctem_yr%no3_mass_yr,                     & !< real, dimension(:,:,:) :
    ngleafmas_yr        =>ctem_yr%ngleafmas_yr,                    & !< real, dimension(:,:,:) :
    ngleafmas_NS_yr      =>ctem_yr%ngleafmas_NS_yr,                & !< real, dimension(:,:,:) :
    ngleafmass_yr       =>ctem_yr%ngleafmass_yr,                   & !< real, dimension(:,:,:) :
    nbleafmas_yr        =>ctem_yr%nbleafmas_yr,                    & !< real, dimension(:,:,:) :
    nstemmass_yr        =>ctem_yr%nstemmass_yr,                    & !< real, dimension(:,:,:) :
    nstemmass_NS_yr      =>ctem_yr%nstemmass_NS_yr,                & !< real, dimension(:,:,:) :
    nstemmasss_yr       =>ctem_yr%nstemmasss_yr,                   & !< real, dimension(:,:,:) :
    nrootmass_yr        =>ctem_yr%nrootmass_yr,                    & !< real, dimension(:,:,:) :
    nrootmass_NS_yr      =>ctem_yr%nrootmass_NS_yr,                & !< real, dimension(:,:,:) :
    nrootmasss_yr       =>ctem_yr%nrootmasss_yr,                   & !< real, dimension(:,:,:) :
    nlitrmass_yr        =>ctem_yr%nlitrmass_yr,                    & !< real, dimension(:,:,:) :
    soilnmas_yr         =>ctem_yr%soilnmas_yr,                     & !< real, dimension(:,:,:) :
    nvgbiomas_yr        =>ctem_yr%nvgbiomas_yr,                    & !< real, dimension(:,:,:) :
    nep_yr              =>ctem_yr%nep_yr,                          & !< real, dimension(:,:,:) :
    litres_yr           =>ctem_yr%litres_yr,                       & !< real, dimension(:,:,:,:) :
    soilcres_yr         =>ctem_yr%soilcres_yr,                     & !< real, dimension(:,:,:,:) :
    litrfall_yr         =>ctem_yr%litrfall_yr,                     & !< real, dimension(:,:,:) :
    hetrores_yr         =>ctem_yr%hetrores_yr,                     & !< real, dimension(:,:,:) :
    nbp_yr              =>ctem_yr%nbp_yr,                          & !< real, dimension(:,:,:) :
    emit_co2_yr         =>ctem_yr%emit_co2_yr,                     & !< real, dimension(:,:,:) :
    emit_co_yr          =>ctem_yr%emit_co_yr,                      & !< real, dimension(:,:,:) :
    emit_ch4_yr         =>ctem_yr%emit_ch4_yr,                     & !< real, dimension(:,:,:) :
    emit_nmhc_yr        =>ctem_yr%emit_nmhc_yr,                    & !< real, dimension(:,:,:) :
    emit_h2_yr          =>ctem_yr%emit_h2_yr,                      & !< real, dimension(:,:,:) :
    emit_nox_yr         =>ctem_yr%emit_nox_yr,                     & !< real, dimension(:,:,:) :
    emit_n2o_yr         =>ctem_yr%emit_n2o_yr,                     & !< real, dimension(:,:,:) :
    emit_nh3_yr         =>ctem_yr%emit_nh3_yr,                     & !< real, dimension(:,:,:) :
    emit_pm25_yr        =>ctem_yr%emit_pm25_yr,                    & !< real, dimension(:,:,:) :
    emit_tpm_yr         =>ctem_yr%emit_tpm_yr,                     & !< real, dimension(:,:,:) :
    emit_tc_yr          =>ctem_yr%emit_tc_yr,                      & !< real, dimension(:,:,:) :
    emit_oc_yr          =>ctem_yr%emit_oc_yr,                      & !< real, dimension(:,:,:) :
    emit_bc_yr          =>ctem_yr%emit_bc_yr,                      & !< real, dimension(:,:,:) :
    bterm_yr            =>ctem_yr%bterm_yr,                        & !< real, dimension(:,:,:) :
    mterm_yr            =>ctem_yr%mterm_yr,                        & !< real, dimension(:,:,:) :
    burnfrac_yr         =>ctem_yr%burnfrac_yr,                     & !< real, dimension(:,:,:) :
    smfuncveg_yr        =>ctem_yr%smfuncveg_yr,                    & !< real, dimension(:,:,:) :
    veghght_yr          =>ctem_yr%veghght_yr,                      & !< real, dimension(:,:,:) :
    bnf_tot_yr          =>ctem_yr%bnf_tot_yr,                      & !< real, dimension(:,:,:) :
    bnf_free_yr         =>ctem_yr%bnf_free_yr,                     & !< real, dimension(:,:,:) :
    bnf_ant_yr          =>ctem_yr%bnf_ant_yr,                      & !< real, dimension(:,:,:) :
    bnf_nat_yr          =>ctem_yr%bnf_nat_yr,                      & !< real, dimension(:,:,:) :
    nitrif_yr           =>ctem_yr%nitrif_yr,                       & !< real, dimension(:,:,:) :
    no_nit_yr           =>ctem_yr%no_nit_yr,                       & !< real, dimension(:,:,:) :
    no_denit_yr         =>ctem_yr%no_denit_yr,                     & !< real, dimension(:,:,:) :
    no_nitdenit_yr      =>ctem_yr%no_nitdenit_yr,                  & !< real, dimension(:,:,:) :
    n2o_nit_yr          =>ctem_yr%n2o_nit_yr,                      & !< real, dimension(:,:,:) :
    n2o_denit_yr        =>ctem_yr%n2o_denit_yr,                    & !< real, dimension(:,:,:) :
    n2o_nitdenit_yr     =>ctem_yr%n2o_nitdenit_yr,                 & !< real, dimension(:,:,:) :
    n2_denit_yr         =>ctem_yr%n2_denit_yr,                     & !< real, dimension(:,:,:) :
    nvol_yr             =>ctem_yr%nvol_yr,                         & !< real, dimension(:,:,:) :
    nleach_yr           =>ctem_yr%nleach_yr,                       & !< real, dimension(:,:,:) :
    appl_fert_yr        =>ctem_yr%appl_fert_yr,                    & !< real, dimension(:,:,:) :
    ndep_nh4_yr         =>ctem_yr%ndep_nh4_yr,                     & !< real, dimension(:,:,:) :
    ndep_no3_yr         =>ctem_yr%ndep_no3_yr,                     & !< real, dimension(:,:,:) :
    ndemand_wp_npp_yr   =>ctem_yr%ndemand_wp_npp_yr,               & !< real, dimension(:,:,:) :
    nuptake_p_nh4_yr    =>ctem_yr%nuptake_p_nh4_yr,                & !< real, dimension(:,:,:) :
    nuptake_p_no3_yr    =>ctem_yr%nuptake_p_no3_yr,                & !< real, dimension(:,:,:) :
    nuptake_a_actl_nh4_yr=>ctem_yr%nuptake_a_actl_nh4_yr,          & !< real, dimension(:,:,:) :
    nuptake_a_actl_no3_yr=>ctem_yr%nuptake_a_actl_no3_yr,          & !< real, dimension(:,:,:) :
    nuptake_yr          =>ctem_yr%nuptake_yr,                      & !< real, dimension(:,:,:) :
    nalloc_l_yr         =>ctem_yr%nalloc_l_yr,                     & !< real, dimension(:,:,:) :
    nalloc_s_yr         =>ctem_yr%nalloc_s_yr,                     & !< real, dimension(:,:,:) :
    nalloc_r_yr         =>ctem_yr%nalloc_r_yr,                     & !< real, dimension(:,:,:) :
    nresorped_s_yr      =>ctem_yr%nresorped_s_yr,                  & !< real, dimension(:,:,:) :
    nresorped_r_yr      =>ctem_yr%nresorped_r_yr,                  & !< real, dimension(:,:,:) :
    nre_alloc_s2l_yr    =>ctem_yr%nre_alloc_s2l_yr,                & !< real, dimension(:,:,:) :
    nre_alloc_r2l_yr    =>ctem_yr%nre_alloc_r2l_yr,                & !< real, dimension(:,:,:) :
    nleafns2s_yr        =>ctem_yr%nleafns2s_yr,                    & !< real, dimension(:,:,:) :
    nstemns2s_yr        =>ctem_yr%nstemns2s_yr,                    & !< real, dimension(:,:,:) :
    nrootns2s_yr        =>ctem_yr%nrootns2s_yr,                    & !< real, dimension(:,:,:) :
    nlitr_l_yr          =>ctem_yr%nlitr_l_yr,                      & !< real, dimension(:,:,:) :
    nlitr_s_yr          =>ctem_yr%nlitr_s_yr,                      & !< real, dimension(:,:,:) :
    nlitr_r_yr          =>ctem_yr%nlitr_r_yr,                      & !< real, dimension(:,:,:) :
    nlitr_yr          =>ctem_yr%nlitr_yr,                      & !< real, dimension(:,:,:) :
    gl2bl_grass_nflux_yr=>ctem_yr%gl2bl_grass_nflux_yr,            & !< real, dimension(:,:,:) :
    nhumtrs_yr          =>ctem_yr%nhumtrs_yr,                      & !< real, dimension(:,:,:) :
    nmineral_litr_yr    =>ctem_yr%nmineral_litr_yr,                & !< real, dimension(:,:,:) :
    nmineral_humus_yr   =>ctem_yr%nmineral_humus_yr,               & !< real, dimension(:,:,:) :
    netnmineral_yr      =>ctem_yr%netnmineral_yr,                  & !< real, dimension(:,:,:) :
    nimmobil_nh4_yr     =>ctem_yr%nimmobil_nh4_yr,                 & !< real, dimension(:,:,:) :
    nimmobil_no3_yr     =>ctem_yr%nimmobil_no3_yr,                 & !< real, dimension(:,:,:) :
    fNnetland_yr        =>ctem_yr%fNnetland_yr,                    & !< real, dimension(:,:,:) :

    laimaxg_yr_t          =>ctem_tile_yr%laimaxg_yr_t,             & !< real, dimension(:,:) :
    c2n_l_yr_t            =>ctem_tile_yr%c2n_l_yr_t,               & !< real, dimension(:,:) :
    c2n_s_yr_t            =>ctem_tile_yr%c2n_s_yr_t,               & !< real, dimension(:,:) :
    c2n_r_yr_t            =>ctem_tile_yr%c2n_r_yr_t,               & !< real, dimension(:,:) :
    c2n_wp_yr_t           =>ctem_tile_yr%c2n_wp_yr_t,              & !< real, dimension(:,:) :
    c2n_litr_yr_t         =>ctem_tile_yr%c2n_litr_yr_t,            & !< real, dimension(:,:) :
    c2n_humus_yr_t        =>ctem_tile_yr%c2n_humus_yr_t,           & !< real, dimension(:,:) :
    gleafmas_yr_t         =>ctem_tile_yr%gleafmas_yr_t,            & !< real, dimension(:,:) :
    gleafmas_NS_yr_t       =>ctem_tile_yr%gleafmas_NS_yr_t,        & !< real, dimension(:,:) :
    gleafmass_yr_t        =>ctem_tile_yr%gleafmass_yr_t,           & !< real, dimension(:,:) :
    bleafmas_yr_t         =>ctem_tile_yr%bleafmas_yr_t,            & !< real, dimension(:,:) :
    stemmass_yr_t         =>ctem_tile_yr%stemmass_yr_t,            & !< real, dimension(:,:) :
    stemmass_NS_yr_t       =>ctem_tile_yr%stemmass_NS_yr_t,        & !< real, dimension(:,:) :
    stemmasss_yr_t        =>ctem_tile_yr%stemmasss_yr_t,           & !< real, dimension(:,:) :
    rootmass_yr_t         =>ctem_tile_yr%rootmass_yr_t,            & !< real, dimension(:,:) :
    rootmass_NS_yr_t       =>ctem_tile_yr%rootmass_NS_yr_t,        & !< real, dimension(:,:) :
    rootmasss_yr_t        =>ctem_tile_yr%rootmasss_yr_t,           & !< real, dimension(:,:) :
    npp_yr_t              =>ctem_tile_yr%npp_yr_t,                 & !< real, dimension(:,:) :
    gpp_yr_t              =>ctem_tile_yr%gpp_yr_t,                 & !< real, dimension(:,:) :
    vcmax0_yr_t           =>ctem_tile_yr%vcmax0_yr_t,              & !< real, dimension(:,:) :
    leafns2s_yr_t         =>ctem_tile_yr%leafns2s_yr_t,            & !< real, dimension(:,:) :
    stemns2s_yr_t         =>ctem_tile_yr%stemns2s_yr_t,            & !< real, dimension(:,:) :
    rootns2s_yr_t         =>ctem_tile_yr%rootns2s_yr_t,            & !< real, dimension(:,:) :
    re_alloc_s2l_yr_t     =>ctem_tile_yr%re_alloc_s2l_yr_t,        & !< real, dimension(:,:) :
    re_alloc_r2l_yr_t     =>ctem_tile_yr%re_alloc_r2l_yr_t,        & !< real, dimension(:,:) :
    re_alloc_sr2l_yr_t    =>ctem_tile_yr%re_alloc_sr2l_yr_t,       & !< real, dimension(:,:) :
    vgbiomas_yr_t         =>ctem_tile_yr%vgbiomas_yr_t,            & !< real, dimension(:,:) :
    autores_yr_t          =>ctem_tile_yr%autores_yr_t,             & !< real, dimension(:,:) :
    rmrveg_yr_t           =>ctem_tile_yr%rmrveg_yr_t,              & !< real, dimension(:,:) :
    totcmass_yr_t         =>ctem_tile_yr%totcmass_yr_t,            & !< real, dimension(:,:) :
    litrmass_yr_t         =>ctem_tile_yr%litrmass_yr_t,            & !< real, dimension(:,:,:) :
    soilcmas_yr_t         =>ctem_tile_yr%soilcmas_yr_t,            & !< real, dimension(:,:,:) :
    nh4_mass_yr_t         =>ctem_tile_yr%nh4_mass_yr_t,            & !< real, dimension(:,:) :
    no3_mass_yr_t         =>ctem_tile_yr%no3_mass_yr_t,            & !< real, dimension(:,:) :
    ngleafmas_yr_t        =>ctem_tile_yr%ngleafmas_yr_t,           & !< real, dimension(:,:) :
    ngleafmas_NS_yr_t      =>ctem_tile_yr%ngleafmas_NS_yr_t,       & !< real, dimension(:,:) :
    ngleafmass_yr_t       =>ctem_tile_yr%ngleafmass_yr_t,          & !< real, dimension(:,:) :
    nbleafmas_yr_t        =>ctem_tile_yr%nbleafmas_yr_t,           & !< real, dimension(:,:) :
    nstemmass_yr_t        =>ctem_tile_yr%nstemmass_yr_t,           & !< real, dimension(:,:) :
    nstemmass_NS_yr_t      =>ctem_tile_yr%nstemmass_NS_yr_t,       & !< real, dimension(:,:) :
    nstemmasss_yr_t       =>ctem_tile_yr%nstemmasss_yr_t,          & !< real, dimension(:,:) :
    nrootmass_yr_t        =>ctem_tile_yr%nrootmass_yr_t,           & !< real, dimension(:,:) :
    nrootmass_NS_yr_t      =>ctem_tile_yr%nrootmass_NS_yr_t,       & !< real, dimension(:,:) :
    nrootmasss_yr_t       =>ctem_tile_yr%nrootmasss_yr_t,          & !< real, dimension(:,:) :
    nlitrmass_yr_t        =>ctem_tile_yr%nlitrmass_yr_t,           & !< real, dimension(:,:) :
    soilnmas_yr_t         =>ctem_tile_yr%soilnmas_yr_t,            & !< real, dimension(:,:) :
    nvgbiomas_yr_t        =>ctem_tile_yr%nvgbiomas_yr_t,           & !< real, dimension(:,:) :
    nep_yr_t              =>ctem_tile_yr%nep_yr_t,                 & !< real, dimension(:,:) :
    litres_yr_t           =>ctem_tile_yr%litres_yr_t,              & !< real, dimension(:,:,:) :
    soilcres_yr_t         =>ctem_tile_yr%soilcres_yr_t,            & !< real, dimension(:,:,:) :
    litrfall_yr_t         =>ctem_tile_yr%litrfall_yr_t,            & !< real, dimension(:,:) :
    hetrores_yr_t         =>ctem_tile_yr%hetrores_yr_t,            & !< real, dimension(:,:) :
    nbp_yr_t              =>ctem_tile_yr%nbp_yr_t,                 & !< real, dimension(:,:) :
    emit_co2_yr_t         =>ctem_tile_yr%emit_co2_yr_t,            & !< real, dimension(:,:) :
    emit_co_yr_t          =>ctem_tile_yr%emit_co_yr_t,             & !< real, dimension(:,:) :
    emit_ch4_yr_t         =>ctem_tile_yr%emit_ch4_yr_t,            & !< real, dimension(:,:) :
    emit_nmhc_yr_t        =>ctem_tile_yr%emit_nmhc_yr_t,           & !< real, dimension(:,:) :
    emit_h2_yr_t          =>ctem_tile_yr%emit_h2_yr_t,             & !< real, dimension(:,:) :
    emit_nox_yr_t         =>ctem_tile_yr%emit_nox_yr_t,            & !< real, dimension(:,:) :
    emit_n2o_yr_t         =>ctem_tile_yr%emit_n2o_yr_t,            & !< real, dimension(:,:) :
    emit_nh3_yr_t         =>ctem_tile_yr%emit_nh3_yr_t,            & !< real, dimension(:,:) :
    emit_pm25_yr_t        =>ctem_tile_yr%emit_pm25_yr_t,           & !< real, dimension(:,:) :
    emit_tpm_yr_t         =>ctem_tile_yr%emit_tpm_yr_t,            & !< real, dimension(:,:) :
    emit_tc_yr_t          =>ctem_tile_yr%emit_tc_yr_t,             & !< real, dimension(:,:) :
    emit_oc_yr_t          =>ctem_tile_yr%emit_oc_yr_t,             & !< real, dimension(:,:) :
    emit_bc_yr_t          =>ctem_tile_yr%emit_bc_yr_t,             & !< real, dimension(:,:) :
    burnfrac_yr_t         =>ctem_tile_yr%burnfrac_yr_t,            & !< real, dimension(:,:) :
    smfuncveg_yr_t        =>ctem_tile_yr%smfuncveg_yr_t,           & !< real, dimension(:,:) :
    bterm_yr_t            =>ctem_tile_yr%bterm_yr_t,               & !< real, dimension(:,:) :
    luc_emc_yr_t          =>ctem_tile_yr%luc_emc_yr_t,             & !< real, dimension(:,:) :
    lterm_yr_t            =>ctem_tile_yr%lterm_yr_t,               & !< real, dimension(:,:) :
    lucsocin_yr_t         =>ctem_tile_yr%lucsocin_yr_t,            & !< real, dimension(:,:) :
    mterm_yr_t            =>ctem_tile_yr%mterm_yr_t,               & !< real, dimension(:,:) :
    lucltrin_yr_t         =>ctem_tile_yr%lucltrin_yr_t,            & !< real, dimension(:,:) :
    tileAge_yr_t          =>ctem_tile_yr%tileAge_yr_t,            & !< real, dimension(:,:) :
    luc_emcn_yr_t          =>ctem_tile_yr%luc_emcn_yr_t,             & !< real, dimension(:,:) :
    lucsocinn_yr_t         =>ctem_tile_yr%lucsocinn_yr_t,            & !< real, dimension(:,:) :
    lucltrinn_yr_t         =>ctem_tile_yr%lucltrinn_yr_t,            & !< real, dimension(:,:) :
    ch4WetSpec_yr_t       =>ctem_tile_yr%ch4WetSpec_yr_t,          & !< real, dimension(:,:) :
    wetfdyn_yr_t          =>ctem_tile_yr%wetfdyn_yr_t,             & !< real, dimension(:,:) :
    ch4WetDyn_yr_t        =>ctem_tile_yr%ch4WetDyn_yr_t,           & !< real, dimension(:,:) :
    ch4soills_yr_t        =>ctem_tile_yr%ch4soills_yr_t,           & !< real, dimension(:,:) :
    veghght_yr_t          =>ctem_tile_yr%veghght_yr_t,             & !< real, dimension(:,:) :
    peatdep_yr_t          =>ctem_tile_yr%peatdep_yr_t,             & !< real, dimension(:,:) :
    peatSoilC_yr_t        =>ctem_tile_yr%peatSoilC_yr_t,           & !< real, dimension(:,:) :
    fProductDecomp_yr_t   =>ctem_tile_yr%fProductDecomp_yr_t,      & !< real, dimension(:,:) : Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    bnf_tot_yr_t        =>ctem_tile_yr%bnf_tot_yr_t,               & !< real, dimension(:,:) :
    bnf_free_yr_t       =>ctem_tile_yr%bnf_free_yr_t,              & !< real, dimension(:,:) :
    bnf_ant_yr_t        =>ctem_tile_yr%bnf_ant_yr_t,               & !< real, dimension(:,:) :
    bnf_nat_yr_t        =>ctem_tile_yr%bnf_nat_yr_t,               & !< real, dimension(:,:) :
    nitrif_yr_t         =>ctem_tile_yr%nitrif_yr_t,                & !< real, dimension(:,:) :
    no_nit_yr_t         =>ctem_tile_yr%no_nit_yr_t,                & !< real, dimension(:,:) :
    no_denit_yr_t       =>ctem_tile_yr%no_denit_yr_t,              & !< real, dimension(:,:) :
    no_nitdenit_yr_t    =>ctem_tile_yr%no_nitdenit_yr_t,           & !< real, dimension(:,:) :
    n2o_nit_yr_t        =>ctem_tile_yr%n2o_nit_yr_t,               & !< real, dimension(:,:) :
    n2o_denit_yr_t      =>ctem_tile_yr%n2o_denit_yr_t,             & !< real, dimension(:,:) :
    n2o_nitdenit_yr_t   =>ctem_tile_yr%n2o_nitdenit_yr_t,          & !< real, dimension(:,:) :
    n2_denit_yr_t       =>ctem_tile_yr%n2_denit_yr_t,              & !< real, dimension(:,:) :
    nvol_yr_t           =>ctem_tile_yr%nvol_yr_t,                  & !< real, dimension(:,:) :
    nleach_yr_t         =>ctem_tile_yr%nleach_yr_t,                & !< real, dimension(:,:) :
    appl_fert_yr_t      =>ctem_tile_yr%appl_fert_yr_t,             & !< real, dimension(:,:) :
    ndep_nh4_yr_t       =>ctem_tile_yr%ndep_nh4_yr_t,              & !< real, dimension(:,:) :
    ndep_no3_yr_t       =>ctem_tile_yr%ndep_no3_yr_t,              & !< real, dimension(:,:) :
    ndemand_wp_npp_yr_t =>ctem_tile_yr%ndemand_wp_npp_yr_t,        & !< real, dimension(:,:) :
    nuptake_p_nh4_yr_t  =>ctem_tile_yr%nuptake_p_nh4_yr_t,         & !< real, dimension(:,:) :
    nuptake_p_no3_yr_t  =>ctem_tile_yr%nuptake_p_no3_yr_t,         & !< real, dimension(:,:) :
    nuptake_a_actl_nh4_yr_t=>ctem_tile_yr%nuptake_a_actl_nh4_yr_t, & !< real, dimension(:,:) :
    nuptake_a_actl_no3_yr_t=>ctem_tile_yr%nuptake_a_actl_no3_yr_t, & !< real, dimension(:,:) :
    nuptake_yr_t        =>ctem_tile_yr%nuptake_yr_t,               & !< real, dimension(:,:) :
    nalloc_l_yr_t       =>ctem_tile_yr%nalloc_l_yr_t,              & !< real, dimension(:,:) :
    nalloc_s_yr_t       =>ctem_tile_yr%nalloc_s_yr_t,              & !< real, dimension(:,:) :
    nalloc_r_yr_t       =>ctem_tile_yr%nalloc_r_yr_t,              & !< real, dimension(:,:) :
    nresorped_s_yr_t    =>ctem_tile_yr%nresorped_s_yr_t,           & !< real, dimension(:,:) :
    nresorped_r_yr_t    =>ctem_tile_yr%nresorped_r_yr_t,           & !< real, dimension(:,:) :
    nre_alloc_s2l_yr_t  =>ctem_tile_yr%nre_alloc_s2l_yr_t,         & !< real, dimension(:,:) :
    nre_alloc_r2l_yr_t  =>ctem_tile_yr%nre_alloc_r2l_yr_t,         & !< real, dimension(:,:) :
    nleafns2s_yr_t      =>ctem_tile_yr%nleafns2s_yr_t,             & !< real, dimension(:,:) :
    nstemns2s_yr_t      =>ctem_tile_yr%nstemns2s_yr_t,             & !< real, dimension(:,:) :
    nrootns2s_yr_t      =>ctem_tile_yr%nrootns2s_yr_t,             & !< real, dimension(:,:) :
    nlitr_l_yr_t        =>ctem_tile_yr%nlitr_l_yr_t,               & !< real, dimension(:,:) :
    nlitr_s_yr_t        =>ctem_tile_yr%nlitr_s_yr_t,               & !< real, dimension(:,:) :
    nlitr_r_yr_t        =>ctem_tile_yr%nlitr_r_yr_t,               & !< real, dimension(:,:) :
    nlitr_yr_t          =>ctem_tile_yr%nlitr_yr_t,                 & !< real, dimension(:,:) :
    gl2bl_grass_nflux_yr_t=>ctem_tile_yr%gl2bl_grass_nflux_yr_t,   & !< real, dimension(:,:) :
    nhumtrs_yr_t        =>ctem_tile_yr%nhumtrs_yr_t,               & !< real, dimension(:,:) :
    nmineral_litr_yr_t  =>ctem_tile_yr%nmineral_litr_yr_t,         & !< real, dimension(:,:) :
    nmineral_humus_yr_t =>ctem_tile_yr%nmineral_humus_yr_t,        & !< real, dimension(:,:) :
    netnmineral_yr_t    =>ctem_tile_yr%netnmineral_yr_t,           & !< real, dimension(:,:) :
    nimmobil_nh4_yr_t   =>ctem_tile_yr%nimmobil_nh4_yr_t,          & !< real, dimension(:,:) :
    nimmobil_no3_yr_t   =>ctem_tile_yr%nimmobil_no3_yr_t,          & !< real, dimension(:,:) :
    fNnetland_yr_t      =>ctem_tile_yr%fNnetland_yr_t,             & !< real, dimension(:,:) :
    timharvarea_yr_t    =>ctem_tile_yr%timharvarea_yr_t,           & !< real, dimension(nlat,nmos)

    pftexistrow       => vrot%pftexist,                            & !< logical, dimension(:,:,:) :
    ccrow       => vrot%cc,                            & !< real, dimension(:,:,:) :
    mmrow       => vrot%mm,                            & !< real, dimension(:,:,:) :
    twarmmrow         => vrot%twarmm,                              & !< real, dimension(:,:) :
    tcoldmrow         => vrot%tcoldm,                              & !< real, dimension(:,:) :
    gdd5row         => vrot%gdd5,                              & !< real, dimension(:,:) :
    aridityrow         => vrot%aridity,                              & !< real, dimension(:,:) :
    srplsmonrow         => vrot%srplsmon,                              & !< real, dimension(:,:) :
    defctmonrow         => vrot%defctmon,                              & !< real, dimension(:,:) :
    anndefctrow         => vrot%anndefct,                              & !< real, dimension(:,:) :
    annsrplsrow         => vrot%annsrpls,                              & !< real, dimension(:,:) :
    annpcprow         => vrot%annpcp,                              & !< real, dimension(:,:) :
    dry_season_lengthrow         => vrot%dry_season_length,                              & !< real, dimension(:,:) :
    gppvegrow         => vrot%gppveg,                              & !< real, dimension(:,:,:) :
    vcmax0row         => vrot%vcmax0,                              & !< real, dimension(:,:,:) :
    leafns2srow       => vrot%leafns2s,                            & !< real, dimension(:,:,:) :
    stemns2srow       => vrot%stemns2s,                            & !< real, dimension(:,:,:) :
    rootns2srow       => vrot%rootns2s,                            & !< real, dimension(:,:,:) :
    re_alloc_s2lrow   => vrot%re_alloc_s2l,                        & !< real, dimension(:,:,:) :
    re_alloc_r2lrow   => vrot%re_alloc_r2l,                        & !< real, dimension(:,:,:) :
    re_alloc_sr2lrow  => vrot%re_alloc_sr2l,                       & !< real, dimension(:,:,:) :
    nepvegrow         => vrot%nepveg,                              & !< real, dimension(:,:,:) :
    nbpvegrow         => vrot%nbpveg,                              & !< real, dimension(:,:,:) :
    nbprow            => vrot%nbp,                                 & !< real, dimension(:,:)   :
    nppvegrow         => vrot%nppveg,                              & !< real, dimension(:,:,:) :
    hetroresvegrow    => vrot%hetroresveg,                         & !< real, dimension(:,:,:) :
    autoresvegrow     => vrot%autoresveg,                          & !< real, dimension(:,:,:) :
    litresvegrow      => vrot%litresveg,                           & !< real, dimension(:,:,:,:) :
    soilcresvegrow    => vrot%soilcresveg,                         & !< real, dimension(:,:,:,:) :
    litrfallvegrow    => vrot%litrfallveg,                         & !< real, dimension(:,:,:) :
    rmlvegaccrow      => vrot%rmlvegacc,                           & !< real, dimension(:,:,:) :
    rmsvegrow         => vrot%rmsveg,                              & !< real, dimension(:,:,:) :
    rmrvegrow         => vrot%rmrveg,                              & !< real, dimension(:,:,:) :
    rgvegrow          => vrot%rgveg,                               & !< real, dimension(:,:,:) :
    ailcgrow          => vrot%ailcg,                               & !< real, dimension(:,:,:) :
    emit_co2row       => vrot%emit_co2,                            & !< real, dimension(:,:,:) :
    emit_corow        => vrot%emit_co,                             & !< real, dimension(:,:,:) :
    emit_ch4row       => vrot%emit_ch4,                            & !< real, dimension(:,:,:) :
    emit_nmhcrow      => vrot%emit_nmhc,                           & !< real, dimension(:,:,:) :
    emit_h2row        => vrot%emit_h2,                             & !< real, dimension(:,:,:) :
    emit_noxrow       => vrot%emit_nox,                            & !< real, dimension(:,:,:) :
    emit_n2orow       => vrot%emit_n2o,                            & !< real, dimension(:,:,:) :
    emit_nh3row       => vrot%emit_nh3,                            & !< real, dimension(:,:,:) :
    emit_pm25row      => vrot%emit_pm25,                           & !< real, dimension(:,:,:) :
    emit_tpmrow       => vrot%emit_tpm,                            & !< real, dimension(:,:,:) :
    emit_tcrow        => vrot%emit_tc,                             & !< real, dimension(:,:,:) :
    emit_ocrow        => vrot%emit_oc,                             & !< real, dimension(:,:,:) :
    emit_bcrow        => vrot%emit_bc,                             & !< real, dimension(:,:,:) :
    burnfracrow       => vrot%burnfrac,                            & !< real, dimension(:,:) :
    burnvegfrow       => vrot%burnvegf,                            & !< real, dimension(:,:,:) :
    smfuncvegrow      => vrot%smfuncveg,                           & !< real, dimension(:,:,:) :
    btermrow          => vrot%bterm,                               & !< real, dimension(:,:,:) :
    ltermrow          => vrot%lterm,                               & !< real, dimension(:,:) :
    mtermrow          => vrot%mterm,                               & !< real, dimension(:,:,:) :
    lucemcomrow       => vrot%lucemcom,                            & !< real, dimension(:,:) :
    lucltrinrow       => vrot%lucltrin,                            & !< real, dimension(:,:) :
    tileAgerow       => vrot%tileAgerow,                            & !< real, dimension(:,:) :
    lucsocinrow       => vrot%lucsocin,                            & !< real, dimension(:,:) :
    lucemcomnrow       => vrot%lucemcomn,                            & !< real, dimension(:,:) :
    lucltrinnrow       => vrot%lucltrinn,                            & !< real, dimension(:,:) :
    lucsocinnrow       => vrot%lucsocinn,                            & !< real, dimension(:,:) :
    timharvarearot    => vrot%timharvarearow,                       & !< real, dimension(:,:) :
    ch4WetSpecrow     => vrot%ch4WetSpec,                          & !< real, dimension(:,:) :
    wetfdynrow        => vrot%wetfdyn,                             & !< real, dimension(:,:) :
    ch4WetDynrow      => vrot%ch4WetDyn,                           & !< real, dimension(:,:) :
    ch4soillsrow      => vrot%ch4_soills,                          & !< real, dimension(:,:) :
    litrmassrow       => vrot%litrmass,                            & !< real, dimension(:,:,:,:) :
    soilcmasrow       => vrot%soilcmas,                            & !< real, dimension(:,:,:,:) :
    nh4_massrow       => vrot%nh4_mass,                            & !< real, dimension(:,:,:) :
    no3_massrow       => vrot%no3_mass,                            & !< real, dimension(:,:,:) :
    nlitrmassrow      => vrot%nlitrmass,                           & !< real, dimension(:,:,:) :
    soilnmasrow       => vrot%soilnmas,                            & !< real, dimension(:,:,:) :
    ngleafmasrow      => vrot%ngleafmas,                           & !< real, dimension(:,:,:) :
    ngleafmas_NSrow    => vrot%ngleafmas_ns,                       & !< real, dimension(:,:,:) :
    ngleafmassrow     => vrot%ngleafmas_s,                         & !< real, dimension(:,:,:) :
    nbleafmasrow      => vrot%nbleafmas,                           & !< real, dimension(:,:,:) :
    nstemmassrow      => vrot%nstemmass,                           & !< real, dimension(:,:,:) :
    nstemmass_NSrow    => vrot%nstemmass_ns,                       & !< real, dimension(:,:,:) :
    nstemmasssrow     => vrot%nstemmass_s,                         & !< real, dimension(:,:,:) :
    nrootmassrow      => vrot%nrootmass,                           & !< real, dimension(:,:,:) :
    nrootmass_NSrow    => vrot%nrootmass_ns,                       & !< real, dimension(:,:,:) :
    nrootmasssrow     => vrot%nrootmass_s,                         & !< real, dimension(:,:,:) :
    vgbiomas_vegrow   => vrot%vgbiomas_veg,                        & !< real, dimension(:,:,:) :
    gleafmasrow       => vrot%gleafmas,                            & !< real, dimension(:,:,:) :
    gleafmas_NSrow     => vrot%gleafmas_ns,                        & !< real, dimension(:,:,:) :
    gleafmassrow      => vrot%gleafmas_s,                          & !< real, dimension(:,:,:) :
    stemmassrow       => vrot%stemmass,                            & !< real, dimension(:,:,:) :
    stemmass_NSrow     => vrot%stemmass_ns,                        & !< real, dimension(:,:,:) :
    stemmasssrow      => vrot%stemmass_s,                          & !< real, dimension(:,:,:) :
    bleafmasrow       => vrot%bleafmas,                            & !< real, dimension(:,:,:) :
    rootmassrow       => vrot%rootmass,                            & !< real, dimension(:,:,:) :
    rootmass_NSrow     => vrot%rootmass_ns,                        & !< real, dimension(:,:,:) :
    rootmasssrow      => vrot%rootmass_s,                          & !< real, dimension(:,:,:) :
    fcancmxrow        => vrot%fcancmx,                             & !< real, dimension(:,:,:) :
    veghghtrow        => vrot%veghght,                             & !< real, dimension(:,:,:) :
    peatdeprow        => vrot%peatdep,                             & !< real, dimension(:,:) :
    peatSoilC        => vrot%peatSoilC,                            & !< real, dimension(:,:) :
    tileAgerot               => vrot%tileAgerow,                   & !< real, dimension(:,:)   :

    bnftotrow        => vrot%bnf_tot,                           & !< real, dimension(:,:,:)   :
    bnffreerow       => vrot%bnf_free,                          & !< real, dimension(:,:,:)   :
    bnfantrow        => vrot%bnf_ant,                           & !< real, dimension(:,:,:)   :
    bnfnatrow        => vrot%bnf_nat,                           & !< real, dimension(:,:,:)   :
    nitrifvegrow      => vrot%nitrifveg,                           & !< real, dimension(:,:,:) :
    no_nitvegrow      => vrot%no_nitveg,                           & !< real, dimension(:,:,:) :
    no_denitvegrow    => vrot%no_denitveg,                         & !< real, dimension(:,:,:) :
    no_nitdenitvegrow => vrot%no_nitdenitveg,                      & !< real, dimension(:,:,:) :
    n2o_nitvegrow     => vrot%n2o_nitveg,                          & !< real, dimension(:,:,:) :
    n2o_denitvegrow   => vrot%n2o_denitveg,                        & !< real, dimension(:,:,:) :
    n2o_nitdenitvegrow=> vrot%n2o_nitdenitveg,                     & !< real, dimension(:,:,:) :
    n2_denitvegrow    => vrot%n2_denitveg,                         & !< real, dimension(:,:,:) :
    nvolvegrow        => vrot%nvolveg,                             & !< real, dimension(:,:,:) :
    nleachvegrow      => vrot%nleachveg,                           & !< real, dimension(:,:,:) :
    appl_fertrow      => vrot%appl_fert,                           & !< real, dimension(:,:,:) :
    ndep_nh4row       => vrot%ndep_nh4,                            & !< real, dimension(:,:,:) :
    ndep_no3row       => vrot%ndep_no3,                            & !< real, dimension(:,:,:) :
    ndemandveg_wp_npprow    => vrot%ndemandveg_wp_npp,             & !< real, dimension(:,:,:) :
    nuptakeveg_p_nh4row     => vrot%nuptakeveg_p_nh4,              & !< real, dimension(:,:,:) :
    nuptakeveg_p_no3row     => vrot%nuptakeveg_p_no3,              & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_nh4row=> vrot%nuptakeveg_a_actl_nh4,         & !< real, dimension(:,:,:) :
    nuptakeveg_a_actl_no3row=> vrot%nuptakeveg_a_actl_no3,         & !< real, dimension(:,:,:) :
    nuptakevegrow     => vrot%nuptakeveg,                          & !< real, dimension(:,:,:) :
    nallocveg_lrow    => vrot%nallocveg_l,                         & !< real, dimension(:,:,:) :
    nallocveg_srow    => vrot%nallocveg_s,                         & !< real, dimension(:,:,:) :
    nallocveg_rrow    => vrot%nallocveg_r,                         & !< real, dimension(:,:,:) :
    nresorpedveg_srow => vrot%nresorpedveg_s,                      & !< real, dimension(:,:,:) :
    nresorpedveg_rrow => vrot%nresorpedveg_r,                      & !< real, dimension(:,:,:) :
    nre_allocveg_s2lrow=> vrot%nre_allocveg_s2l,                   & !< real, dimension(:,:,:) :
    nre_allocveg_r2lrow=> vrot%nre_allocveg_r2l,                   & !< real, dimension(:,:,:) :
    nleafns2svegrow   => vrot%nleafns2sveg,                        & !< real, dimension(:,:,:) :
    nstemns2svegrow   => vrot%nstemns2sveg,                        & !< real, dimension(:,:,:) :
    nrootns2svegrow   => vrot%nrootns2sveg,                        & !< real, dimension(:,:,:) :
    nlitrveg_lrow     => vrot%nlitrveg_l,                          & !< real, dimension(:,:,:) :
    nlitrveg_srow     => vrot%nlitrveg_s,                          & !< real, dimension(:,:,:) :
    nlitrveg_rrow     => vrot%nlitrveg_r,                          & !< real, dimension(:,:,:) :
    nlitrvegrow       => vrot%nlitrveg,                            & !< real, dimension(:,:,:) :
    gl2bl_grass_nfluxrow=> vrot%gl2bl_grass_nflux,                 & !< real, dimension(:,:,:) :
    nhumtrsvegrow     => vrot%nhumtrsveg,                          & !< real, dimension(:,:,:) :
    nmineralveg_litrrow=> vrot%nmineralveg_litr,                   & !< real, dimension(:,:,:) :
    nmineralveg_humusrow=> vrot%nmineralveg_humus,                 & !< real, dimension(:,:,:) :
    netnmineralveg_row=> vrot%netnmineralveg,                      & !< real, dimension(:,:,:) :
    nimmobilveg_nh4row=> vrot%nimmobilveg_nh4,                     & !< real, dimension(:,:,:) :
    nimmobilveg_no3row=> vrot%nimmobilveg_no3,                     & !< real, dimension(:,:,:) :
    fNnetlandvegrow   => vrot%fNnetlandveg,                        & !< real, dimension(:,:,:) :
    c2nveg_lrow       => vrot%c2nveg_l,                            & !< real, dimension(:,:,:) :
    c2nveg_srow       => vrot%c2nveg_s,                            & !< real, dimension(:,:,:) :
    c2nveg_rrow       => vrot%c2nveg_r,                            & !< real, dimension(:,:,:) :
    c2nveg_wprow      => vrot%c2nveg_wp,                           & !< real, dimension(:,:,:) :
    c2nveg_litrrow    => vrot%c2nveg_litr,                         & !< real, dimension(:,:,:) :
    c2nveg_humusrow   => vrot%c2nveg_humus,                        & !< real, dimension(:,:,:) :
    nvgbiomas_vegrow  => vrot%nvgbiomas_veg,                       & !< real, dimension(:,:,:) :

    laimaxg_yr_g          =>ctem_grd_yr%laimaxg_yr_g,              & !< real, dimension(:) :
    c2n_l_yr_g            =>ctem_grd_yr%c2n_l_yr_g,                & !< real, dimension(:) :
    c2n_s_yr_g            =>ctem_grd_yr%c2n_s_yr_g,                & !< real, dimension(:) :
    c2n_r_yr_g            =>ctem_grd_yr%c2n_r_yr_g,                & !< real, dimension(:) :
    c2n_wp_yr_g           =>ctem_grd_yr%c2n_wp_yr_g,               & !< real, dimension(:) :
    c2n_litr_yr_g         =>ctem_grd_yr%c2n_litr_yr_g,             & !< real, dimension(:) :
    c2n_humus_yr_g        =>ctem_grd_yr%c2n_humus_yr_g,            & !< real, dimension(:) :
    gleafmas_yr_g         =>ctem_grd_yr%gleafmas_yr_g,             & !< real, dimension(:) :
    gleafmas_NS_yr_g       =>ctem_grd_yr%gleafmas_NS_yr_g,         & !< real, dimension(:) :
    gleafmass_yr_g        =>ctem_grd_yr%gleafmass_yr_g,            & !< real, dimension(:) :
    stemmass_yr_g         =>ctem_grd_yr%stemmass_yr_g,             & !< real, dimension(:) :
    stemmass_NS_yr_g       =>ctem_grd_yr%stemmass_NS_yr_g,         & !< real, dimension(:) :
    stemmasss_yr_g        =>ctem_grd_yr%stemmasss_yr_g,            & !< real, dimension(:) :
    bleafmas_yr_g         =>ctem_grd_yr%bleafmas_yr_g,             & !< real, dimension(:) :
    rootmass_yr_g         =>ctem_grd_yr%rootmass_yr_g,             & !< real, dimension(:) :
    rootmass_NS_yr_g       =>ctem_grd_yr%rootmass_NS_yr_g,         & !< real, dimension(:) :
    rootmasss_yr_g        =>ctem_grd_yr%rootmasss_yr_g,            & !< real, dimension(:) :
    litrmass_yr_g         =>ctem_grd_yr%litrmass_yr_g,             & !< real, dimension(:,:) :
    soilcmas_yr_g         =>ctem_grd_yr%soilcmas_yr_g,             & !< real, dimension(:,:) :
    nh4_mass_yr_g         =>ctem_grd_yr%nh4_mass_yr_g,             & !< real, dimension(:) :
    no3_mass_yr_g         =>ctem_grd_yr%no3_mass_yr_g,             & !< real, dimension(:) :
    ngleafmas_yr_g        =>ctem_grd_yr%ngleafmas_yr_g,            & !< real, dimension(:) :
    ngleafmas_NS_yr_g      =>ctem_grd_yr%ngleafmas_NS_yr_g,        & !< real, dimension(:) :
    ngleafmass_yr_g       =>ctem_grd_yr%ngleafmass_yr_g,           & !< real, dimension(:) :
    nbleafmas_yr_g        =>ctem_grd_yr%nbleafmas_yr_g,            & !< real, dimension(:) :
    nstemmass_yr_g        =>ctem_grd_yr%nstemmass_yr_g,            & !< real, dimension(:) :
    nstemmass_NS_yr_g      =>ctem_grd_yr%nstemmass_NS_yr_g,        & !< real, dimension(:) :
    nstemmasss_yr_g       =>ctem_grd_yr%nstemmasss_yr_g,           & !< real, dimension(:) :
    nrootmass_yr_g        =>ctem_grd_yr%nrootmass_yr_g,            & !< real, dimension(:) :
    nrootmass_NS_yr_g      =>ctem_grd_yr%nrootmass_NS_yr_g,        & !< real, dimension(:) :
    nrootmasss_yr_g       =>ctem_grd_yr%nrootmasss_yr_g,           & !< real, dimension(:) :
    nlitrmass_yr_g        =>ctem_grd_yr%nlitrmass_yr_g,            & !< real, dimension(:) :
    soilnmas_yr_g         =>ctem_grd_yr%soilnmas_yr_g,             & !< real, dimension(:) :
    npp_yr_g              =>ctem_grd_yr%npp_yr_g,                  & !< real, dimension(:) :
    gpp_yr_g              =>ctem_grd_yr%gpp_yr_g,                  & !< real, dimension(:) :
    vcmax0_yr_g           =>ctem_grd_yr%vcmax0_yr_g,               & !< real, dimension(:) :
    leafns2s_yr_g         =>ctem_grd_yr%leafns2s_yr_g,             & !< real, dimension(:) :
    stemns2s_yr_g         =>ctem_grd_yr%stemns2s_yr_g,             & !< real, dimension(:) :
    rootns2s_yr_g         =>ctem_grd_yr%rootns2s_yr_g,             & !< real, dimension(:) :
    re_alloc_s2l_yr_g     =>ctem_grd_yr%re_alloc_s2l_yr_g,         & !< real, dimension(:) :
    re_alloc_r2l_yr_g     =>ctem_grd_yr%re_alloc_r2l_yr_g,         & !< real, dimension(:) :
    re_alloc_sr2l_yr_g    =>ctem_grd_yr%re_alloc_sr2l_yr_g,        & !< real, dimension(:) :
    nep_yr_g              =>ctem_grd_yr%nep_yr_g,                  & !< real, dimension(:) :
    nbp_yr_g              =>ctem_grd_yr%nbp_yr_g,                  & !< real, dimension(:) :
    hetrores_yr_g         =>ctem_grd_yr%hetrores_yr_g,             & !< real, dimension(:) :
    autores_yr_g          =>ctem_grd_yr%autores_yr_g,              & !< real, dimension(:) :
    rmrveg_yr_g           =>ctem_grd_yr%rmrveg_yr_g,               & !< real, dimension(:) :
    litres_yr_g           =>ctem_grd_yr%litres_yr_g,               & !< real, dimension(:,:) :
    soilcres_yr_g         =>ctem_grd_yr%soilcres_yr_g,             & !< real, dimension(:,:) :
    litrfall_yr_g         =>ctem_grd_yr%litrfall_yr_g,             & !< real, dimension(:) :
    vgbiomas_yr_g         =>ctem_grd_yr%vgbiomas_yr_g,             & !< real, dimension(:) :
    totcmass_yr_g         =>ctem_grd_yr%totcmass_yr_g,             & !< real, dimension(:) :
    emit_co2_yr_g         =>ctem_grd_yr%emit_co2_yr_g,             & !< real, dimension(:) :
    emit_co_yr_g          =>ctem_grd_yr%emit_co_yr_g,              & !< real, dimension(:) :
    emit_ch4_yr_g         =>ctem_grd_yr%emit_ch4_yr_g,             & !< real, dimension(:) :
    emit_nmhc_yr_g        =>ctem_grd_yr%emit_nmhc_yr_g,            & !< real, dimension(:) :
    emit_h2_yr_g          =>ctem_grd_yr%emit_h2_yr_g,              & !< real, dimension(:) :
    emit_nox_yr_g         =>ctem_grd_yr%emit_nox_yr_g,             & !< real, dimension(:) :
    emit_n2o_yr_g         =>ctem_grd_yr%emit_n2o_yr_g,             & !< real, dimension(:) :
    emit_nh3_yr_g         =>ctem_grd_yr%emit_nh3_yr_g,             & !< real, dimension(:) :
    emit_pm25_yr_g        =>ctem_grd_yr%emit_pm25_yr_g,            & !< real, dimension(:) :
    emit_tpm_yr_g         =>ctem_grd_yr%emit_tpm_yr_g,             & !< real, dimension(:) :
    emit_tc_yr_g          =>ctem_grd_yr%emit_tc_yr_g,              & !< real, dimension(:) :
    emit_oc_yr_g          =>ctem_grd_yr%emit_oc_yr_g,              & !< real, dimension(:) :
    emit_bc_yr_g          =>ctem_grd_yr%emit_bc_yr_g,              & !< real, dimension(:) :
    smfuncveg_yr_g        =>ctem_grd_yr%smfuncveg_yr_g,            & !< real, dimension(:) :
    luc_emc_yr_g          =>ctem_grd_yr%luc_emc_yr_g,              & !< real, dimension(:) :
    lucltrin_yr_g         =>ctem_grd_yr%lucltrin_yr_g,             & !< real, dimension(:) :
    tileAge_yr_g         =>ctem_grd_yr%tileAge_yr_g,             & !< real, dimension(:) :
    lucsocin_yr_g         =>ctem_grd_yr%lucsocin_yr_g,             & !< real, dimension(:) :
    luc_emcn_yr_g          =>ctem_grd_yr%luc_emcn_yr_g,              & !< real, dimension(:) :
    lucltrinn_yr_g         =>ctem_grd_yr%lucltrinn_yr_g,             & !< real, dimension(:) :
    lucsocinn_yr_g         =>ctem_grd_yr%lucsocinn_yr_g,             & !< real, dimension(:) :
    burnfrac_yr_g         =>ctem_grd_yr%burnfrac_yr_g,             & !< real, dimension(:) :
    lterm_yr_g            =>ctem_grd_yr%lterm_yr_g,                & !< real, dimension(:) :
    bterm_yr_g            =>ctem_grd_yr%bterm_yr_g,                & !< real, dimension(:) :
    mterm_yr_g            =>ctem_grd_yr%mterm_yr_g,                & !< real, dimension(:) :
    ch4WetSpec_yr_g       =>ctem_grd_yr%ch4WetSpec_yr_g,           & !< real, dimension(:) :
    wetfdyn_yr_g          =>ctem_grd_yr%wetfdyn_yr_g,              & !< real, dimension(:) :
    ch4WetDyn_yr_g        =>ctem_grd_yr%ch4WetDyn_yr_g,            & !< real, dimension(:) :
    ch4soills_yr_g        =>ctem_grd_yr%ch4soills_yr_g,            & !< real, dimension(:) :
    veghght_yr_g          =>ctem_grd_yr%veghght_yr_g,              & !< real, dimension(:) :
    peatdep_yr_g          =>ctem_grd_yr%peatdep_yr_g,              & !< real, dimension(:) :
    peatSoilC_yr_g        =>ctem_grd_yr%peatSoilC_yr_g,            & !< real, dimension(:) :
    cProduct_yr_g         =>ctem_grd_yr%cProduct_yr_g,             & !< real, dimension(:) : Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    nProduct_yr_g         =>ctem_grd_yr%nProduct_yr_g,             & !< real, dimension(:) : Nitrogen in the LUC product pools (litter and soil N iccp2 position) \f$[g N m^{-2}]\f$
    fProductDecomp_yr_g   =>ctem_grd_yr%fProductDecomp_yr_g,       & !< real, dimension(:) : Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$
    bnf_tot_yr_g          =>ctem_grd_yr%bnf_tot_yr_g,              & !< real, dimension(:) :
    bnf_free_yr_g         =>ctem_grd_yr%bnf_free_yr_g,             & !< real, dimension(:) :
    bnf_ant_yr_g          =>ctem_grd_yr%bnf_ant_yr_g,              & !< real, dimension(:) :
    bnf_nat_yr_g          =>ctem_grd_yr%bnf_nat_yr_g,              & !< real, dimension(:) :
    nitrif_yr_g           =>ctem_grd_yr%nitrif_yr_g,               & !< real, dimension(:) :
    no_nit_yr_g           =>ctem_grd_yr%no_nit_yr_g,               & !< real, dimension(:) :
    no_denit_yr_g         =>ctem_grd_yr%no_denit_yr_g,             & !< real, dimension(:) :
    no_nitdenit_yr_g      =>ctem_grd_yr%no_nitdenit_yr_g,          & !< real, dimension(:) :
    n2o_nit_yr_g          =>ctem_grd_yr%n2o_nit_yr_g,              & !< real, dimension(:) :
    n2o_denit_yr_g        =>ctem_grd_yr%n2o_denit_yr_g,            & !< real, dimension(:) :
    n2o_nitdenit_yr_g     =>ctem_grd_yr%n2o_nitdenit_yr_g,         & !< real, dimension(:) :
    n2_denit_yr_g         =>ctem_grd_yr%n2_denit_yr_g,             & !< real, dimension(:) :
    nvol_yr_g             =>ctem_grd_yr%nvol_yr_g,                 & !< real, dimension(:) :
    nleach_yr_g           =>ctem_grd_yr%nleach_yr_g,               & !< real, dimension(:) :
    appl_fert_yr_g        =>ctem_grd_yr%appl_fert_yr_g,            & !< real, dimension(:) :
    ndep_nh4_yr_g         =>ctem_grd_yr%ndep_nh4_yr_g,             & !< real, dimension(:) :
    ndep_no3_yr_g         =>ctem_grd_yr%ndep_no3_yr_g,             & !< real, dimension(:) :
    ndemand_wp_npp_yr_g   =>ctem_grd_yr%ndemand_wp_npp_yr_g,       & !< real, dimension(:) :
    nuptake_p_nh4_yr_g    =>ctem_grd_yr%nuptake_p_nh4_yr_g,        & !< real, dimension(:) :
    nuptake_p_no3_yr_g    =>ctem_grd_yr%nuptake_p_no3_yr_g,        & !< real, dimension(:) :
    nuptake_a_actl_nh4_yr_g=>ctem_grd_yr%nuptake_a_actl_nh4_yr_g,  & !< real, dimension(:) :
    nuptake_a_actl_no3_yr_g=>ctem_grd_yr%nuptake_a_actl_no3_yr_g,  & !< real, dimension(:) :
    nuptake_yr_g          =>ctem_grd_yr%nuptake_yr_g,              & !< real, dimension(:) :
    nalloc_l_yr_g         =>ctem_grd_yr%nalloc_l_yr_g,             & !< real, dimension(:) :
    nalloc_s_yr_g         =>ctem_grd_yr%nalloc_s_yr_g,             & !< real, dimension(:) :
    nalloc_r_yr_g         =>ctem_grd_yr%nalloc_r_yr_g,             & !< real, dimension(:) :
    nresorped_s_yr_g      =>ctem_grd_yr%nresorped_s_yr_g,          & !< real, dimension(:) :
    nresorped_r_yr_g      =>ctem_grd_yr%nresorped_r_yr_g,          & !< real, dimension(:) :
    nre_alloc_s2l_yr_g    =>ctem_grd_yr%nre_alloc_s2l_yr_g,        & !< real, dimension(:) :
    nre_alloc_r2l_yr_g    =>ctem_grd_yr%nre_alloc_r2l_yr_g,        & !< real, dimension(:) :
    nleafns2s_yr_g        =>ctem_grd_yr%nleafns2s_yr_g,            & !< real, dimension(:) :
    nstemns2s_yr_g        =>ctem_grd_yr%nstemns2s_yr_g,            & !< real, dimension(:) :
    nrootns2s_yr_g        =>ctem_grd_yr%nrootns2s_yr_g,            & !< real, dimension(:) :
    nlitr_l_yr_g          =>ctem_grd_yr%nlitr_l_yr_g,              & !< real, dimension(:) :
    nlitr_s_yr_g          =>ctem_grd_yr%nlitr_s_yr_g,              & !< real, dimension(:) :
    nlitr_r_yr_g          =>ctem_grd_yr%nlitr_r_yr_g,              & !< real, dimension(:) :
    nlitr_yr_g            =>ctem_grd_yr%nlitr_yr_g,                & !< real, dimension(:) :
    gl2bl_grass_nflux_yr_g=>ctem_grd_yr%gl2bl_grass_nflux_yr_g,    & !< real, dimension(:) :
    nhumtrs_yr_g          =>ctem_grd_yr%nhumtrs_yr_g,              & !< real, dimension(:) :
    nmineral_litr_yr_g    =>ctem_grd_yr%nmineral_litr_yr_g,        & !< real, dimension(:) :
    nmineral_humus_yr_g   =>ctem_grd_yr%nmineral_humus_yr_g,       & !< real, dimension(:) :
    netnmineral_yr_g      =>ctem_grd_yr%netnmineral_yr_g,          & !< real, dimension(:) :
    nimmobil_nh4_yr_g     =>ctem_grd_yr%nimmobil_nh4_yr_g,         & !< real, dimension(:) :
    nimmobil_no3_yr_g     =>ctem_grd_yr%nimmobil_no3_yr_g,         & !< real, dimension(:) :
    nvgbiomas_yr_g        =>ctem_grd_yr%nvgbiomas_yr_g,            & !< real, dimension(:) :
    fNnetland_yr_g        =>ctem_grd_yr%fNnetland_yr_g             & !< real, dimension(:) :
    )

    !------------
    !> Accumulate yearly outputs
    oneOverDPY = 1./real(lastDOY)
    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do m = 1,nmtest
      do j = 1,icc

        !> Accumulate the variables at the per PFT level

        if (ailcgrow(i,m,j) > laimaxg_yr(i,m,j)) then
          laimaxg_yr(i,m,j) = ailcgrow(i,m,j)
        end if
        cc_yr(i,m,j) = cc_yr(i,m,j) + ccrow(i,m,j) * oneOverDPY
        mm_yr(i,m,j) = mm_yr(i,m,j) + mmrow(i,m,j) * oneOverDPY
        npp_yr(i,m,j) = npp_yr(i,m,j) + nppvegrow(i,m,j) * oneOverDPY
        gpp_yr(i,m,j) = gpp_yr(i,m,j) + gppvegrow(i,m,j) * oneOverDPY
        vcmax0_yr(i,m,j) = vcmax0_yr(i,m,j) + vcmax0row(i,m,j) * oneOverDPY
        leafns2s_yr(i,m,j) = leafns2s_yr(i,m,j) + leafns2srow(i,m,j) * oneOverDPY
        stemns2s_yr(i,m,j) = stemns2s_yr(i,m,j) + stemns2srow(i,m,j) * oneOverDPY
        rootns2s_yr(i,m,j) = rootns2s_yr(i,m,j) + rootns2srow(i,m,j) * oneOverDPY
        re_alloc_s2l_yr(i,m,j) = re_alloc_s2l_yr(i,m,j) + re_alloc_s2lrow(i,m,j) * oneOverDPY
        re_alloc_r2l_yr(i,m,j) = re_alloc_r2l_yr(i,m,j) + re_alloc_r2lrow(i,m,j) * oneOverDPY
        re_alloc_sr2l_yr(i,m,j) = re_alloc_sr2l_yr(i,m,j) + re_alloc_sr2lrow(i,m,j) * oneOverDPY
        nep_yr(i,m,j) = nep_yr(i,m,j) + nepvegrow(i,m,j) * oneOverDPY
        ! NOTE: This NBP does not include LUC product pool contributions since they are
        ! not per PFT but rather per tile
        nbp_yr(i,m,j) = nbp_yr(i,m,j) + nbpvegrow(i,m,j) * oneOverDPY
        litrfall_yr(i,m,j) = litrfall_yr(i,m,j) + litrfallvegrow(i,m,j) * oneOverDPY
        emit_co2_yr(i,m,j) = emit_co2_yr(i,m,j) + emit_co2row(i,m,j) * oneOverDPY
        emit_co_yr(i,m,j) = emit_co_yr(i,m,j) + emit_corow(i,m,j) * oneOverDPY
        emit_ch4_yr(i,m,j) = emit_ch4_yr(i,m,j) + emit_ch4row(i,m,j) * oneOverDPY
        emit_nmhc_yr(i,m,j) = emit_nmhc_yr(i,m,j) + emit_nmhcrow(i,m,j) * oneOverDPY
        emit_h2_yr(i,m,j) = emit_h2_yr(i,m,j) + emit_h2row(i,m,j) * oneOverDPY
        emit_nox_yr(i,m,j) = emit_nox_yr(i,m,j) + emit_noxrow(i,m,j) * oneOverDPY
        emit_n2o_yr(i,m,j) = emit_n2o_yr(i,m,j) + emit_n2orow(i,m,j) * oneOverDPY
        emit_nh3_yr(i,m,j) = emit_nh3_yr(i,m,j) + emit_nh3row(i,m,j) * oneOverDPY
        emit_pm25_yr(i,m,j) = emit_pm25_yr(i,m,j) + emit_pm25row(i,m,j) * oneOverDPY
        emit_tpm_yr(i,m,j) = emit_tpm_yr(i,m,j) + emit_tpmrow(i,m,j) * oneOverDPY
        emit_tc_yr(i,m,j) = emit_tc_yr(i,m,j) + emit_tcrow(i,m,j) * oneOverDPY
        emit_oc_yr(i,m,j) = emit_oc_yr(i,m,j) + emit_ocrow(i,m,j) * oneOverDPY
        emit_bc_yr(i,m,j) = emit_bc_yr(i,m,j) + emit_bcrow(i,m,j) * oneOverDPY

        bterm_yr(i,m,j) = bterm_yr(i,m,j) + btermrow(i,m,j) * oneOverDPY
        mterm_yr(i,m,j) = mterm_yr(i,m,j) + mtermrow(i,m,j) * oneOverDPY
        smfuncveg_yr(i,m,j) = smfuncveg_yr(i,m,j) + smfuncvegrow(i,m,j) * oneOverDPY
        hetrores_yr(i,m,j) = hetrores_yr(i,m,j) + hetroresvegrow(i,m,j) * oneOverDPY
        autores_yr(i,m,j) = autores_yr(i,m,j) + autoresvegrow(i,m,j) * oneOverDPY
        rmrveg_yr(i,m,j) = rmrveg_yr(i,m,j) + rmrvegrow(i,m,j) * oneOverDPY
        do k = 1,ignd
          litres_yr(i,m,j,k) = litres_yr(i,m,j,k) + litresvegrow(i,m,j,k) * oneOverDPY
          soilcres_yr(i,m,j,k) = soilcres_yr(i,m,j,k) + soilcresvegrow(i,m,j,k) * oneOverDPY
        end do
        ! Let accumulate, not a flux or a mean value.
        burnfrac_yr(i,m,j) = burnfrac_yr(i,m,j) + burnvegfrow(i,m,j)

        ! Note that the following are recorded in units of g N m^-2 yr^-1:
        if (Ncycle_on) then
          appl_fert_yr(i,m,j) = appl_fert_yr(i,m,j) + appl_fertrow(i,m,j) * oneOverDPY * convertkgN
          ndep_nh4_yr(i,m,j) = ndep_nh4_yr(i,m,j) + ndep_nh4row(i,m,j) * oneOverDPY * convertkgN
          ndep_no3_yr(i,m,j) = ndep_no3_yr(i,m,j) + ndep_no3row(i,m,j) * oneOverDPY * convertkgN
          bnf_free_yr(i,m,j) = bnf_free_yr(i,m,j) + bnffreerow(i,m,j) * oneOverDPY * convertkgN
          bnf_ant_yr(i,m,j) = bnf_ant_yr(i,m,j) + bnfantrow(i,m,j) * oneOverDPY * convertkgN
          bnf_nat_yr(i,m,j) = bnf_nat_yr(i,m,j) + bnfnatrow(i,m,j) * oneOverDPY * convertkgN
          bnf_tot_yr(i,m,j) = bnf_tot_yr(i,m,j) + bnftotrow(i,m,j) * oneOverDPY * convertkgN
          nitrif_yr(i,m,j) = nitrif_yr(i,m,j) + nitrifvegrow(i,m,j) * oneOverDPY * convertkgN
          no_nit_yr(i,m,j) = no_nit_yr(i,m,j) + no_nitvegrow(i,m,j) * oneOverDPY * convertkgN
          no_denit_yr(i,m,j) = no_denit_yr(i,m,j) + no_denitvegrow(i,m,j) * oneOverDPY * convertkgN
          no_nitdenit_yr(i,m,j) = no_nitdenit_yr(i,m,j) + no_nitdenitvegrow(i,m,j) * oneOverDPY * convertkgN
          n2o_nit_yr(i,m,j) = n2o_nit_yr(i,m,j) + n2o_nitvegrow(i,m,j) * oneOverDPY * convertkgN
          n2o_denit_yr(i,m,j) = n2o_denit_yr(i,m,j) + n2o_denitvegrow(i,m,j) * oneOverDPY * convertkgN
          n2o_nitdenit_yr(i,m,j) = n2o_nitdenit_yr(i,m,j) + n2o_nitdenitvegrow(i,m,j) * oneOverDPY * convertkgN
          n2_denit_yr(i,m,j) = n2_denit_yr(i,m,j) + n2_denitvegrow(i,m,j) * oneOverDPY * convertkgN
          nvol_yr(i,m,j) = nvol_yr(i,m,j) + nvolvegrow(i,m,j) * oneOverDPY * convertkgN
          nleach_yr(i,m,j) = nleach_yr(i,m,j) + nleachvegrow(i,m,j) * oneOverDPY * convertkgN
          ndemand_wp_npp_yr(i,m,j) = ndemand_wp_npp_yr(i,m,j) + ndemandveg_wp_npprow(i,m,j) * oneOverDPY * convertkgN
          nuptake_p_nh4_yr(i,m,j) = nuptake_p_nh4_yr(i,m,j) + nuptakeveg_p_nh4row(i,m,j) * oneOverDPY * convertkgN
          nuptake_p_no3_yr(i,m,j) = nuptake_p_no3_yr(i,m,j) + nuptakeveg_p_no3row(i,m,j) * oneOverDPY * convertkgN
          nuptake_a_actl_nh4_yr(i,m,j) = nuptake_a_actl_nh4_yr(i,m,j) &
                                         + nuptakeveg_a_actl_nh4row(i,m,j) * oneOverDPY * convertkgN
          nuptake_a_actl_no3_yr(i,m,j) = nuptake_a_actl_no3_yr(i,m,j) &
                                         + nuptakeveg_a_actl_no3row(i,m,j) * oneOverDPY * convertkgN
          nuptake_yr(i,m,j) = nuptake_yr(i,m,j) + nuptakevegrow(i,m,j) * oneOverDPY * convertkgN
          nalloc_l_yr(i,m,j) = nalloc_l_yr(i,m,j) + nallocveg_lrow(i,m,j) * oneOverDPY * convertkgN
          nalloc_s_yr(i,m,j) = nalloc_s_yr(i,m,j) + nallocveg_srow(i,m,j) * oneOverDPY * convertkgN
          nalloc_r_yr(i,m,j) = nalloc_r_yr(i,m,j) + nallocveg_rrow(i,m,j) * oneOverDPY * convertkgN
          nresorped_s_yr(i,m,j) = nresorped_s_yr(i,m,j) + nresorpedveg_srow(i,m,j) * oneOverDPY * convertkgN
          nresorped_r_yr(i,m,j) = nresorped_r_yr(i,m,j) + nresorpedveg_rrow(i,m,j) * oneOverDPY * convertkgN
          nre_alloc_s2l_yr(i,m,j) = nre_alloc_s2l_yr(i,m,j) + nre_allocveg_s2lrow(i,m,j) * oneOverDPY * convertkgN
          nre_alloc_r2l_yr(i,m,j) = nre_alloc_r2l_yr(i,m,j) + nre_allocveg_r2lrow(i,m,j) * oneOverDPY * convertkgN
          nleafns2s_yr(i,m,j) = nleafns2s_yr(i,m,j) + nleafns2svegrow(i,m,j) * oneOverDPY * convertkgN
          nstemns2s_yr(i,m,j) = nstemns2s_yr(i,m,j) + nstemns2svegrow(i,m,j) * oneOverDPY * convertkgN
          nrootns2s_yr(i,m,j) = nrootns2s_yr(i,m,j) + nrootns2svegrow(i,m,j) * oneOverDPY * convertkgN
          nlitr_l_yr(i,m,j) = nlitr_l_yr(i,m,j) + nlitrveg_lrow(i,m,j) * oneOverDPY * convertkgN
          nlitr_s_yr(i,m,j) = nlitr_s_yr(i,m,j) + nlitrveg_srow(i,m,j) * oneOverDPY * convertkgN
          nlitr_r_yr(i,m,j) = nlitr_r_yr(i,m,j) + nlitrveg_rrow(i,m,j) * oneOverDPY * convertkgN
          nlitr_yr(i,m,j) = nlitr_yr(i,m,j) + nlitrvegrow(i,m,j) * oneOverDPY * convertkgN
          gl2bl_grass_nflux_yr(i,m,j) = gl2bl_grass_nflux_yr(i,m,j) + gl2bl_grass_nfluxrow(i,m,j) * oneOverDPY * convertkgN
          nhumtrs_yr(i,m,j) = nhumtrs_yr(i,m,j) + nhumtrsvegrow(i,m,j) * oneOverDPY * convertkgN
          nmineral_litr_yr(i,m,j) = nmineral_litr_yr(i,m,j) + nmineralveg_litrrow(i,m,j) * oneOverDPY * convertkgN
          nmineral_humus_yr(i,m,j) = nmineral_humus_yr(i,m,j) + nmineralveg_humusrow(i,m,j) * oneOverDPY * convertkgN
          netnmineral_yr(i,m,j) = netnmineral_yr(i,m,j) + netnmineralveg_row(i,m,j) * oneOverDPY * convertkgN
          nimmobil_nh4_yr(i,m,j) = nimmobil_nh4_yr(i,m,j) + nimmobilveg_nh4row(i,m,j) * oneOverDPY * convertkgN
          nimmobil_no3_yr(i,m,j) = nimmobil_no3_yr(i,m,j) + nimmobilveg_no3row(i,m,j) * oneOverDPY * convertkgN
          fNnetland_yr(i,m,j) = fNnetland_yr(i,m,j) + fNnetlandvegrow(i,m,j) * oneOverDPY * convertkgN
        end if !Ncycle_on

      end do ! loop 884

      !>   Also do the bare fraction amounts
      hetrores_yr(i,m,iccp1) = hetrores_yr(i,m,iccp1) + hetroresvegrow(i,m,iccp1) * oneOverDPY
      do k = 1,ignd
        litres_yr(i,m,iccp1,k) = litres_yr(i,m,iccp1,k) + litresvegrow(i,m,iccp1,k) * oneOverDPY
        soilcres_yr(i,m,iccp1,k) = soilcres_yr(i,m,iccp1,k) &
                                   + soilcresvegrow(i,m,iccp1,k) * oneOverDPY
      end do
      nep_yr(i,m,iccp1) = nep_yr(i,m,iccp1) + nepvegrow(i,m,iccp1) * oneOverDPY
      nbp_yr(i,m,iccp1) = nbp_yr(i,m,iccp1) + nbpvegrow(i,m,iccp1) * oneOverDPY

      if (Ncycle_on) then
          appl_fert_yr(i,m,iccp1) = appl_fert_yr(i,m,iccp1) + appl_fertrow(i,m,iccp1) * oneOverDPY * convertkgN
          ndep_nh4_yr(i,m,iccp1) = ndep_nh4_yr(i,m,iccp1) + ndep_nh4row(i,m,iccp1) * oneOverDPY * convertkgN
          ndep_no3_yr(i,m,iccp1) = ndep_no3_yr(i,m,iccp1) + ndep_no3row(i,m,iccp1) * oneOverDPY * convertkgN
          bnf_free_yr(i,m,iccp1) = bnf_free_yr(i,m,iccp1) + bnffreerow(i,m,iccp1) * oneOverDPY * convertkgN
          bnf_tot_yr(i,m,iccp1) = bnf_tot_yr(i,m,iccp1) + bnftotrow(i,m,iccp1) * oneOverDPY * convertkgN
          nitrif_yr(i,m,iccp1) = nitrif_yr(i,m,iccp1) + nitrifvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          no_nit_yr(i,m,iccp1) = no_nit_yr(i,m,iccp1) + no_nitvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          no_denit_yr(i,m,iccp1) = no_denit_yr(i,m,iccp1) + no_denitvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          no_nitdenit_yr(i,m,iccp1) = no_nitdenit_yr(i,m,iccp1) + no_nitdenitvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          n2o_nit_yr(i,m,iccp1) = n2o_nit_yr(i,m,iccp1) + n2o_nitvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          n2o_denit_yr(i,m,iccp1) = n2o_denit_yr(i,m,iccp1) + n2o_denitvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          n2o_nitdenit_yr(i,m,iccp1) = n2o_nitdenit_yr(i,m,iccp1) + n2o_nitdenitvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          n2_denit_yr(i,m,iccp1) = n2_denit_yr(i,m,iccp1) + n2_denitvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          nvol_yr(i,m,iccp1) = nvol_yr(i,m,iccp1) + nvolvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          nleach_yr(i,m,iccp1) = nleach_yr(i,m,iccp1) + nleachvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          nhumtrs_yr(i,m,iccp1) = nhumtrs_yr(i,m,iccp1) + nhumtrsvegrow(i,m,iccp1) * oneOverDPY * convertkgN
          nmineral_litr_yr(i,m,iccp1) = nmineral_litr_yr(i,m,iccp1) + nmineralveg_litrrow(i,m,iccp1) * oneOverDPY * convertkgN
          nmineral_humus_yr(i,m,iccp1) = nmineral_humus_yr(i,m,iccp1) + nmineralveg_humusrow(i,m,iccp1) * oneOverDPY * convertkgN
          netnmineral_yr(i,m,iccp1) = netnmineral_yr(i,m,iccp1) + netnmineralveg_row(i,m,iccp1) * oneOverDPY * convertkgN
          nimmobil_nh4_yr(i,m,iccp1) = nimmobil_nh4_yr(i,m,iccp1) + nimmobilveg_nh4row(i,m,iccp1) * oneOverDPY * convertkgN
          nimmobil_no3_yr(i,m,iccp1) = nimmobil_no3_yr(i,m,iccp1) + nimmobilveg_no3row(i,m,iccp1) * oneOverDPY * convertkgN
          fNnetland_yr(i,m,iccp1) = fNnetland_yr(i,m,iccp1) + fNnetlandvegrow(i,m,iccp1) * oneOverDPY * convertkgN
        end if

      peatdep_yr_t(i,m) = peatdeprow(i,m)      ! YW September 04,2015
      peatSoilC_yr_t(i,m) = peatSoilC(i,m)     !s.r.c added peat carbon output
      timharvarea_yr_t(i,m) = timharvarea_yr_t(i,m) + timharvarearot(i,m)

      !> Accumulate the variables at the per tile level
      lterm_yr_t(i,m) = lterm_yr_t(i,m) + ltermrow(i,m) * oneOverDPY
      wetfdyn_yr_t(i,m) = wetfdyn_yr_t(i,m) + wetfdynrow(i,m) * oneOverDPY
      luc_emc_yr_t(i,m) = luc_emc_yr_t(i,m) + lucemcomrow(i,m) * oneOverDPY
      lucsocin_yr_t(i,m) = lucsocin_yr_t(i,m) + lucsocinrow(i,m) * oneOverDPY
      lucltrin_yr_t(i,m) = lucltrin_yr_t(i,m) + lucltrinrow(i,m) * oneOverDPY
      tileAge_yr_t(i,m) = tileAge_yr_t(i,m) + tileAgerot(i,m) * oneOverDPY
      if (Ncycle_on) then
        luc_emcn_yr_t(i,m) = luc_emcn_yr_t(i,m) + lucemcomnrow(i,m) * oneOverDPY
        lucsocinn_yr_t(i,m) = lucsocinn_yr_t(i,m) + lucsocinnrow(i,m) * oneOverDPY
        lucltrinn_yr_t(i,m) = lucltrinn_yr_t(i,m) + lucltrinnrow(i,m) * oneOverDPY
      end if
      ch4WetSpec_yr_t(i,m) = ch4WetSpec_yr_t(i,m) + ch4WetSpecrow(i,m) * oneOverDPY
      ch4WetDyn_yr_t(i,m) = ch4WetDyn_yr_t(i,m) + ch4WetDynrow(i,m) * oneOverDPY
      ch4soills_yr_t(i,m) = ch4soills_yr_t(i,m) + ch4soillsrow(i,m) * oneOverDPY

      ! NOTE: LUC product pools are only in layer 1.
      fProductDecomp_yr_t(i,m) = fProductDecomp_yr_t(i,m) &
                                + (soilcresvegrow(i,m,iccp2,1) &
                                + litresvegrow(i,m,iccp2,1)) * oneOverDPY

      ! NOTE: NBP is a special case here. The LUC product pool contributions are not
      ! per PFT, they exist uniformly across a tile, so they are not included in the
      ! nbp_yr calculation. Instead we need to use the nbp, not nbpveg variable
      ! for per tile and per gridcell outputting.
      nbp_yr_t(i,m) = nbp_yr_t(i,m) + nbprow(i,m) * oneOverDPY

    end do ! loop 883 ! m

    if (iday == lastDOY) then
      do m = 1,nmtest
        do j = 1,icc

          !> The pools are looked at just at the end of the year.
          gleafmas_yr(i,m,j) = gleafmasrow(i,m,j)
          gleafmas_NS_yr(i,m,j) = gleafmas_NSrow(i,m,j)
          gleafmass_yr(i,m,j) = gleafmassrow(i,m,j)
          bleafmas_yr(i,m,j) = bleafmasrow(i,m,j)
          stemmass_yr(i,m,j) = stemmassrow(i,m,j)
          stemmass_NS_yr(i,m,j) = stemmass_NSrow(i,m,j)
          stemmasss_yr(i,m,j) = stemmasssrow(i,m,j)
          rootmass_yr(i,m,j) = rootmassrow(i,m,j)
          rootmass_NS_yr(i,m,j) = rootmass_NSrow(i,m,j)
          rootmasss_yr(i,m,j) = rootmasssrow(i,m,j)
          veghght_yr(i,m,j) = veghghtrow(i,m,j)
          vgbiomas_yr(i,m,j) = vgbiomas_vegrow(i,m,j)
          totcmass_yr(i,m,j) = vgbiomas_yr(i,m,j)
          do k = 1,ignd
            litrmass_yr(i,m,j,k) = litrmassrow(i,m,j,k)
            soilcmas_yr(i,m,j,k) = soilcmasrow(i,m,j,k)
            totcmass_yr(i,m,j) = totcmass_yr(i,m,j) + litrmass_yr(i,m,j,k) &
                                                    + soilcmas_yr(i,m,j,k)
          end do

          if (Ncycle_on) then
              nh4_mass_yr(i,m,j) = nh4_massrow(i,m,j) * convertg2kg
              no3_mass_yr(i,m,j) = no3_massrow(i,m,j) * convertg2kg
              ngleafmas_yr(i,m,j) = ngleafmasrow(i,m,j) * convertg2kg
              ngleafmas_NS_yr(i,m,j) = ngleafmas_NSrow(i,m,j) * convertg2kg
              ngleafmass_yr(i,m,j) = ngleafmassrow(i,m,j) * convertg2kg
              nbleafmas_yr(i,m,j) = nbleafmasrow(i,m,j) * convertg2kg
              nstemmass_yr(i,m,j) = nstemmassrow(i,m,j) * convertg2kg
              nstemmass_NS_yr(i,m,j) = nstemmass_NSrow(i,m,j) * convertg2kg
              nstemmasss_yr(i,m,j) = nstemmasssrow(i,m,j) * convertg2kg
              nrootmass_yr(i,m,j) = nrootmassrow(i,m,j) * convertg2kg
              nrootmass_NS_yr(i,m,j) = nrootmass_NSrow(i,m,j) * convertg2kg
              nrootmasss_yr(i,m,j) = nrootmasssrow(i,m,j) * convertg2kg
              nlitrmass_yr(i,m,j) = nlitrmassrow(i,m,j) * convertg2kg
              soilnmas_yr(i,m,j) = soilnmasrow(i,m,j) * convertg2kg
              nvgbiomas_yr(i,m,j) = nvgbiomas_vegrow(i,m,j) * convertg2kg
              if ((ngleafmas_yr(i,m,j) + nbleafmas_yr(i,m,j)) /= 0.0) c2n_l_yr(i,m,j) = &
                  (gleafmas_yr(i,m,j) + bleafmas_yr(i,m,j)) / (ngleafmas_yr(i,m,j) + nbleafmas_yr(i,m,j))
              if (nstemmass_yr(i,m,j) /= 0.0) c2n_s_yr(i,m,j) = stemmass_yr(i,m,j) / nstemmass_yr(i,m,j)
              if (nrootmass_yr(i,m,j) /= 0.0) c2n_r_yr(i,m,j) = rootmass_yr(i,m,j) / nrootmass_yr(i,m,j)
              if ((ngleafmas_yr(i,m,j) + nbleafmas_yr(i,m,j) + nstemmass_yr(i,m,j) + nrootmass_yr(i,m,j)) /= 0.0) &
                c2n_wp_yr(i,m,j) = (gleafmas_yr(i,m,j) + bleafmas_yr(i,m,j) + stemmass_yr(i,m,j) + rootmass_yr(i,m,j)) &
                                 / (ngleafmas_yr(i,m,j) + nbleafmas_yr(i,m,j) + nstemmass_yr(i,m,j) + nrootmass_yr(i,m,j))
              if (nlitrmass_yr(i,m,j) /= 0.0) c2n_litr_yr(i,m,j) = sum(litrmass_yr(i,m,j,:)) / nlitrmass_yr(i,m,j)
              if (soilnmas_yr(i,m,j) /= 0.0) c2n_humus_yr(i,m,j) = sum(soilcmas_yr(i,m,j,:)) / soilnmas_yr(i,m,j)
            end if ! Ncycle_on
        end do ! j 

        peatdep_yr_g(i) = peatdep_yr_g(i) + peatdep_yr_t(i,m) * farerot(i,m)    ! YW September 04,2015
        peatSoilC_yr_g(i) = peatSoilC_yr_g(i) + peatSoilC_yr_t(i,m) * farerot(i,m)    !s.r.c added peat carbon output

        do k = 1,ignd
          litrmass_yr(i,m,iccp1,k) = litrmassrow(i,m,iccp1,k)
          soilcmas_yr(i,m,iccp1,k) = soilcmasrow(i,m,iccp1,k)
          totcmass_yr(i,m,iccp1) = totcmass_yr(i,m,iccp1) + litrmassrow(i,m,iccp1,k) + soilcmasrow(i,m,iccp1,k)
        end do
        barefrac = 1.0
        if (Ncycle_on) then
          nh4_mass_yr(i,m,iccp1) = nh4_massrow(i,m,iccp1) * convertg2kg
          no3_mass_yr(i,m,iccp1) = no3_massrow(i,m,iccp1) * convertg2kg
          nlitrmass_yr(i,m,iccp1)= nlitrmassrow(i,m,iccp1) * convertg2kg
          soilnmas_yr(i,m,iccp1) = soilnmasrow(i,m,iccp1) * convertg2kg
          if (nlitrmass_yr(i,m,iccp1) /= 0.0) c2n_litr_yr(i,m,iccp1) = &
              sum(litrmass_yr(i,m,iccp1,:)) / nlitrmass_yr(i,m,iccp1)
          if (soilnmas_yr(i,m,iccp1) /= 0.0) c2n_humus_yr(i,m,iccp1) = &
              sum(soilcmas_yr(i,m,iccp1,:)) / soilnmas_yr(i,m,iccp1)
        end if ! Ncycle_on

        !> Add values to the per tile vars
        ! NOTE: This implictly assumes that the fcancmx is only changing annually !
        do j = 1,icc

          laimaxg_yr_t(i,m) = laimaxg_yr_t(i,m) + laimaxg_yr(i,m,j) * fcancmxrow(i,m,j)
          gleafmas_yr_t(i,m) = gleafmas_yr_t(i,m) + gleafmas_yr(i,m,j) * fcancmxrow(i,m,j)
          gleafmas_NS_yr_t(i,m) = gleafmas_NS_yr_t(i,m) + gleafmas_NS_yr(i,m,j) * fcancmxrow(i,m,j)
          gleafmass_yr_t(i,m) = gleafmass_yr_t(i,m) + gleafmass_yr(i,m,j) * fcancmxrow(i,m,j)
          bleafmas_yr_t(i,m) = bleafmas_yr_t(i,m) + bleafmas_yr(i,m,j) * fcancmxrow(i,m,j)
          stemmass_yr_t(i,m) = stemmass_yr_t(i,m) + stemmass_yr(i,m,j) * fcancmxrow(i,m,j)
          stemmass_NS_yr_t(i,m) = stemmass_NS_yr_t(i,m) + stemmass_NS_yr(i,m,j) * fcancmxrow(i,m,j)
          stemmasss_yr_t(i,m) = stemmasss_yr_t(i,m) + stemmasss_yr(i,m,j) * fcancmxrow(i,m,j)
          rootmass_yr_t(i,m) = rootmass_yr_t(i,m) + rootmass_yr(i,m,j) * fcancmxrow(i,m,j)
          rootmass_NS_yr_t(i,m) = rootmass_NS_yr_t(i,m) + rootmass_NS_yr(i,m,j) * fcancmxrow(i,m,j)
          rootmasss_yr_t(i,m) = rootmasss_yr_t(i,m) + rootmasss_yr(i,m,j) * fcancmxrow(i,m,j)

          if (Ncycle_on) then
            nh4_mass_yr_t(i,m) = nh4_mass_yr_t(i,m) + nh4_mass_yr(i,m,j) * fcancmxrow(i,m,j)
            no3_mass_yr_t(i,m) = no3_mass_yr_t(i,m) + no3_mass_yr(i,m,j) * fcancmxrow(i,m,j)
            ngleafmas_yr_t(i,m) = ngleafmas_yr_t(i,m) + ngleafmas_yr(i,m,j) * fcancmxrow(i,m,j)
            ngleafmas_NS_yr_t(i,m) = ngleafmas_NS_yr_t(i,m) + ngleafmas_NS_yr(i,m,j) * fcancmxrow(i,m,j)
            ngleafmass_yr_t(i,m) = ngleafmass_yr_t(i,m) + ngleafmass_yr(i,m,j) * fcancmxrow(i,m,j)
            nbleafmas_yr_t(i,m) = nbleafmas_yr_t(i,m) + nbleafmas_yr(i,m,j) * fcancmxrow(i,m,j)
            nstemmass_yr_t(i,m) = nstemmass_yr_t(i,m) + nstemmass_yr(i,m,j) * fcancmxrow(i,m,j)
            nstemmass_NS_yr_t(i,m) = nstemmass_NS_yr_t(i,m) + nstemmass_NS_yr(i,m,j) * fcancmxrow(i,m,j)
            nstemmasss_yr_t(i,m) = nstemmasss_yr_t(i,m) + nstemmasss_yr(i,m,j) * fcancmxrow(i,m,j)
            nrootmass_yr_t(i,m) = nrootmass_yr_t(i,m) + nrootmass_yr(i,m,j) * fcancmxrow(i,m,j)
            nrootmass_NS_yr_t(i,m) = nrootmass_NS_yr_t(i,m) + nrootmass_NS_yr(i,m,j) * fcancmxrow(i,m,j)
            nrootmasss_yr_t(i,m) = nrootmasss_yr_t(i,m) + nrootmasss_yr(i,m,j) * fcancmxrow(i,m,j)
            nlitrmass_yr_t(i,m) = nlitrmass_yr_t(i,m) + nlitrmass_yr(i,m,j) * fcancmxrow(i,m,j)
            soilnmas_yr_t(i,m) = soilnmas_yr_t(i,m) + soilnmas_yr(i,m,j) * fcancmxrow(i,m,j)
            appl_fert_yr_t(i,m) = appl_fert_yr_t(i,m) + appl_fert_yr(i,m,j) * fcancmxrow(i,m,j)
            ndep_nh4_yr_t(i,m) = ndep_nh4_yr_t(i,m) + ndep_nh4_yr(i,m,j) * fcancmxrow(i,m,j)
            ndep_no3_yr_t(i,m) = ndep_no3_yr_t(i,m) + ndep_no3_yr(i,m,j) * fcancmxrow(i,m,j)
            bnf_free_yr_t(i,m) = bnf_free_yr_t(i,m) + bnf_free_yr(i,m,j) * fcancmxrow(i,m,j)
            bnf_ant_yr_t(i,m) = bnf_ant_yr_t(i,m) + bnf_ant_yr(i,m,j) * fcancmxrow(i,m,j)
            bnf_nat_yr_t(i,m) = bnf_nat_yr_t(i,m) + bnf_nat_yr(i,m,j) * fcancmxrow(i,m,j)
            bnf_tot_yr_t(i,m) = bnf_tot_yr_t(i,m) + bnf_tot_yr(i,m,j) * fcancmxrow(i,m,j)
            nitrif_yr_t(i,m) = nitrif_yr_t(i,m) + nitrif_yr(i,m,j) * fcancmxrow(i,m,j)
            no_nit_yr_t(i,m) = no_nit_yr_t(i,m) + no_nit_yr(i,m,j) * fcancmxrow(i,m,j)
            no_denit_yr_t(i,m) = no_denit_yr_t(i,m) + no_denit_yr(i,m,j) * fcancmxrow(i,m,j)
            no_nitdenit_yr_t(i,m) = no_nitdenit_yr_t(i,m) &
                                    + no_nitdenit_yr(i,m,j) * fcancmxrow(i,m,j)
            n2o_nit_yr_t(i,m) = n2o_nit_yr_t(i,m) + n2o_nit_yr(i,m,j) * fcancmxrow(i,m,j)
            n2o_denit_yr_t(i,m) = n2o_denit_yr_t(i,m) + n2o_denit_yr(i,m,j) * fcancmxrow(i,m,j)
            n2o_nitdenit_yr_t(i,m) = n2o_nitdenit_yr_t(i,m) &
                                     + n2o_nitdenit_yr(i,m,j) * fcancmxrow(i,m,j)
            n2_denit_yr_t(i,m) = n2_denit_yr_t(i,m) + n2_denit_yr(i,m,j) * fcancmxrow(i,m,j)
            nvol_yr_t(i,m) = nvol_yr_t(i,m) + nvol_yr(i,m,j) * fcancmxrow(i,m,j)
            nleach_yr_t(i,m) = nleach_yr_t(i,m) + nleach_yr(i,m,j) * fcancmxrow(i,m,j)
            ndemand_wp_npp_yr_t(i,m) = ndemand_wp_npp_yr_t(i,m) &
                                       + ndemand_wp_npp_yr(i,m,j) * fcancmxrow(i,m,j)
            nuptake_p_nh4_yr_t(i,m) = nuptake_p_nh4_yr_t(i,m) &
                                      + nuptake_p_nh4_yr(i,m,j) * fcancmxrow(i,m,j)
            nuptake_p_no3_yr_t(i,m) = nuptake_p_no3_yr_t(i,m) &
                                      + nuptake_p_no3_yr(i,m,j) * fcancmxrow(i,m,j)
            nuptake_a_actl_nh4_yr_t(i,m) = nuptake_a_actl_nh4_yr_t(i,m) &
                                           + nuptake_a_actl_nh4_yr(i,m,j) * fcancmxrow(i,m,j)
            nuptake_a_actl_no3_yr_t(i,m) = nuptake_a_actl_no3_yr_t(i,m) &
                                           + nuptake_a_actl_no3_yr(i,m,j) * fcancmxrow(i,m,j)
            nuptake_yr_t(i,m) = nuptake_yr_t(i,m) + nuptake_yr(i,m,j) * fcancmxrow(i,m,j)
            nalloc_l_yr_t(i,m) = nalloc_l_yr_t(i,m) + nalloc_l_yr(i,m,j) * fcancmxrow(i,m,j)
            nalloc_s_yr_t(i,m) = nalloc_s_yr_t(i,m) + nalloc_s_yr(i,m,j) * fcancmxrow(i,m,j)
            nalloc_r_yr_t(i,m) = nalloc_r_yr_t(i,m) + nalloc_r_yr(i,m,j) * fcancmxrow(i,m,j)
            nresorped_s_yr_t(i,m) = nresorped_s_yr_t(i,m) &
                                    + nresorped_s_yr(i,m,j) * fcancmxrow(i,m,j)
            nresorped_r_yr_t(i,m) = nresorped_r_yr_t(i,m) &
                                    + nresorped_r_yr(i,m,j) * fcancmxrow(i,m,j)
            nre_alloc_s2l_yr_t(i,m) = nre_alloc_s2l_yr_t(i,m) &
                                      + nre_alloc_s2l_yr(i,m,j) * fcancmxrow(i,m,j)
            nre_alloc_r2l_yr_t(i,m) = nre_alloc_r2l_yr_t(i,m) &
                                      + nre_alloc_r2l_yr(i,m,j) * fcancmxrow(i,m,j)
            nleafns2s_yr_t(i,m) = nleafns2s_yr_t(i,m) + nleafns2s_yr(i,m,j) * fcancmxrow(i,m,j)
            nstemns2s_yr_t(i,m) = nstemns2s_yr_t(i,m) + nstemns2s_yr(i,m,j) * fcancmxrow(i,m,j)
            nrootns2s_yr_t(i,m) = nrootns2s_yr_t(i,m) + nrootns2s_yr(i,m,j) * fcancmxrow(i,m,j)
            nlitr_l_yr_t(i,m) = nlitr_l_yr_t(i,m) + nlitr_l_yr(i,m,j) * fcancmxrow(i,m,j)
            nlitr_s_yr_t(i,m) = nlitr_s_yr_t(i,m) + nlitr_s_yr(i,m,j) * fcancmxrow(i,m,j)
            nlitr_r_yr_t(i,m) = nlitr_r_yr_t(i,m) + nlitr_r_yr(i,m,j) * fcancmxrow(i,m,j)
            nlitr_yr_t(i,m) = nlitr_yr_t(i,m) + nlitr_yr(i,m,j) * fcancmxrow(i,m,j)
            gl2bl_grass_nflux_yr_t(i,m) = gl2bl_grass_nflux_yr_t(i,m) &
                                          + gl2bl_grass_nflux_yr(i,m,j) * fcancmxrow(i,m,j)
            nhumtrs_yr_t(i,m) = nhumtrs_yr_t(i,m) + nhumtrs_yr(i,m,j) * fcancmxrow(i,m,j)
            nmineral_litr_yr_t(i,m) = nmineral_litr_yr_t(i,m) &
                                      + nmineral_litr_yr(i,m,j) * fcancmxrow(i,m,j)
            nmineral_humus_yr_t(i,m) = nmineral_humus_yr_t(i,m) &
                                       + nmineral_humus_yr(i,m,j) * fcancmxrow(i,m,j)
            netnmineral_yr_t(i,m) = netnmineral_yr_t(i,m) &
                                       + netnmineral_yr(i,m,j) * fcancmxrow(i,m,j)
            nimmobil_nh4_yr_t(i,m) = nimmobil_nh4_yr_t(i,m) &
                                     + nimmobil_nh4_yr(i,m,j) * fcancmxrow(i,m,j)
            nimmobil_no3_yr_t(i,m) = nimmobil_no3_yr_t(i,m) &
                                     + nimmobil_no3_yr(i,m,j) * fcancmxrow(i,m,j)
            nvgbiomas_yr_t(i,m) = nvgbiomas_yr_t(i,m) + nvgbiomas_yr(i,m,j) * fcancmxrow(i,m,j)
            fNnetland_yr_t(i,m) = fNnetland_yr_t(i,m) + fNnetland_yr(i,m,j) * fcancmxrow(i,m,j)
          end if !for Ncycle_on
          vgbiomas_yr_t(i,m) = vgbiomas_yr_t(i,m) + vgbiomas_yr(i,m,j) * fcancmxrow(i,m,j)
          totcmass_yr_t(i,m) = totcmass_yr_t(i,m) + totcmass_yr(i,m,j) * fcancmxrow(i,m,j)
          npp_yr_t(i,m) = npp_yr_t(i,m) + npp_yr(i,m,j) * fcancmxrow(i,m,j)
          gpp_yr_t(i,m) = gpp_yr_t(i,m) + gpp_yr(i,m,j) * fcancmxrow(i,m,j)
          vcmax0_yr_t(i,m) = vcmax0_yr_t(i,m) + vcmax0_yr(i,m,j) * fcancmxrow(i,m,j)
          leafns2s_yr_t(i,m) = leafns2s_yr_t(i,m) + leafns2s_yr(i,m,j) * fcancmxrow(i,m,j)
          stemns2s_yr_t(i,m) = stemns2s_yr_t(i,m) + stemns2s_yr(i,m,j) * fcancmxrow(i,m,j)
          rootns2s_yr_t(i,m) = rootns2s_yr_t(i,m) + rootns2s_yr(i,m,j) * fcancmxrow(i,m,j)
          re_alloc_s2l_yr_t(i,m) = re_alloc_s2l_yr_t(i,m) + re_alloc_s2l_yr(i,m,j) * fcancmxrow(i,m,j)
          re_alloc_r2l_yr_t(i,m) = re_alloc_r2l_yr_t(i,m) + re_alloc_r2l_yr(i,m,j) * fcancmxrow(i,m,j)
          re_alloc_sr2l_yr_t(i,m) = re_alloc_sr2l_yr_t(i,m) + re_alloc_sr2l_yr(i,m,j) * fcancmxrow(i,m,j)
          nep_yr_t(i,m) = nep_yr_t(i,m) + nep_yr(i,m,j) * fcancmxrow(i,m,j)
          ! nbp_yr_t(i,m)=nbp_yr_t(i,m)+nbp_yr(i,m,j)*fcancmxrow(i,m,j)
          emit_co2_yr_t(i,m) = emit_co2_yr_t(i,m) + emit_co2_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_co_yr_t(i,m) = emit_co_yr_t(i,m) + emit_co_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_ch4_yr_t(i,m) = emit_ch4_yr_t(i,m) + emit_ch4_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_nmhc_yr_t(i,m) = emit_nmhc_yr_t(i,m) + emit_nmhc_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_h2_yr_t(i,m) = emit_h2_yr_t(i,m) + emit_h2_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_nox_yr_t(i,m) = emit_nox_yr_t(i,m) + emit_nox_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_n2o_yr_t(i,m) = emit_n2o_yr_t(i,m) + emit_n2o_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_nh3_yr_t(i,m) = emit_nh3_yr_t(i,m) + emit_nh3_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_pm25_yr_t(i,m) = emit_pm25_yr_t(i,m) + emit_pm25_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_tpm_yr_t(i,m) = emit_tpm_yr_t(i,m) + emit_tpm_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_tc_yr_t(i,m) = emit_tc_yr_t(i,m) + emit_tc_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_oc_yr_t(i,m) = emit_oc_yr_t(i,m) + emit_oc_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_bc_yr_t(i,m) = emit_bc_yr_t(i,m) + emit_bc_yr(i,m,j) * fcancmxrow(i,m,j)
          bterm_yr_t(i,m) = bterm_yr_t(i,m) + bterm_yr(i,m,j) * fcancmxrow(i,m,j)
          mterm_yr_t(i,m) = mterm_yr_t(i,m) + mterm_yr(i,m,j) * fcancmxrow(i,m,j)
          smfuncveg_yr_t(i,m) = smfuncveg_yr_t(i,m) + smfuncveg_yr(i,m,j) * fcancmxrow(i,m,j)
          hetrores_yr_t(i,m) = hetrores_yr_t(i,m) + hetrores_yr(i,m,j) * fcancmxrow(i,m,j)
          autores_yr_t(i,m) = autores_yr_t(i,m) + autores_yr(i,m,j) * fcancmxrow(i,m,j)
          rmrveg_yr_t(i,m) = rmrveg_yr_t(i,m) + rmrveg_yr(i,m,j) * fcancmxrow(i,m,j)
          litrfall_yr_t(i,m) = litrfall_yr_t(i,m) + litrfall_yr(i,m,j) * fcancmxrow(i,m,j)
          burnfrac_yr_t(i,m) = burnfrac_yr_t(i,m) + burnfrac_yr(i,m,j) * fcancmxrow(i,m,j)
          veghght_yr_t(i,m) = veghght_yr_t(i,m) + veghght_yr(i,m,j) * fcancmxrow(i,m,j)

          barefrac = barefrac - fcancmxrow(i,m,j)
          do k = 1,ignd
            litrmass_yr_t(i,m,k) = litrmass_yr_t(i,m,k) + litrmass_yr(i,m,j,k) * fcancmxrow(i,m,j)
            soilcmas_yr_t(i,m,k) = soilcmas_yr_t(i,m,k) + soilcmas_yr(i,m,j,k) * fcancmxrow(i,m,j)
            litres_yr_t(i,m,k) = litres_yr_t(i,m,k) + litres_yr(i,m,j,k) * fcancmxrow(i,m,j)
            soilcres_yr_t(i,m,k) = soilcres_yr_t(i,m,k) + soilcres_yr(i,m,j,k) * fcancmxrow(i,m,j)
          end do

        end do ! j
        if(Ncycle_on) then
          if ((ngleafmas_yr_t(i,m) + nbleafmas_yr_t(i,m)) /= 0.0) c2n_l_yr_t(i,m) = &
              (gleafmas_yr_t(i,m) + bleafmas_yr_t(i,m)) / (ngleafmas_yr_t(i,m) + nbleafmas_yr_t(i,m))
          if (nstemmass_yr_t(i,m) /= 0.0) c2n_s_yr_t(i,m) = stemmass_yr_t(i,m) / nstemmass_yr_t(i,m)
          if (nrootmass_yr_t(i,m) /= 0.0) c2n_r_yr_t(i,m) = rootmass_yr_t(i,m) / nrootmass_yr_t(i,m)
          if ((ngleafmas_yr_t(i,m) + nbleafmas_yr_t(i,m) + nstemmass_yr_t(i,m) + nrootmass_yr_t(i,m)) /= 0.0) &
              c2n_wp_yr_t(i,m) = (gleafmas_yr_t(i,m) + bleafmas_yr_t(i,m) + stemmass_yr_t(i,m) + rootmass_yr_t(i,m)) &
                                / (ngleafmas_yr_t(i,m) + nbleafmas_yr_t(i,m) + nstemmass_yr_t(i,m) + nrootmass_yr_t(i,m))
        end if

        do k = 1,ignd
          litrmass_yr_t(i,m,k) = litrmass_yr_t(i,m,k) + litrmass_yr(i,m,iccp1,k) * barefrac
          soilcmas_yr_t(i,m,k) = soilcmas_yr_t(i,m,k) + soilcmas_yr(i,m,iccp1,k) * barefrac
          litres_yr_t(i,m,k) = litres_yr_t(i,m,k) + litres_yr(i,m,iccp1,k) * barefrac
          soilcres_yr_t(i,m,k) = soilcres_yr_t(i,m,k) + soilcres_yr(i,m,iccp1,k) * barefrac
        totcmass_yr_t(i,m) = totcmass_yr_t(i,m) &
                            + (litrmass_yr(i,m,iccp1,k) + soilcmas_yr(i,m,iccp1,k)) * barefrac
        end do

        hetrores_yr_t(i,m) = hetrores_yr_t(i,m) + hetrores_yr(i,m,iccp1) * barefrac
        nep_yr_t(i,m) = nep_yr_t(i,m) + nep_yr(i,m,iccp1) * barefrac
        ! nbp_yr_t(i,m)=nbp_yr_t(i,m)+nbp_yr(i,m,iccp1)*barefrac

        if (Ncycle_on) then
          nh4_mass_yr_t(i,m) = nh4_mass_yr_t(i,m) + nh4_mass_yr(i,m,iccp1) * barefrac
          no3_mass_yr_t(i,m) = no3_mass_yr_t(i,m) + no3_mass_yr(i,m,iccp1) * barefrac
          nlitrmass_yr_t(i,m) = nlitrmass_yr_t(i,m) + nlitrmass_yr(i,m,iccp1) * barefrac
          soilnmas_yr_t(i,m) = soilnmas_yr_t(i,m) + soilnmas_yr(i,m,iccp1) * barefrac
          appl_fert_yr_t(i,m) = appl_fert_yr_t(i,m) + appl_fert_yr(i,m,iccp1) * barefrac
          ndep_nh4_yr_t(i,m) = ndep_nh4_yr_t(i,m) + ndep_nh4_yr(i,m,iccp1) * barefrac
          ndep_no3_yr_t(i,m) = ndep_no3_yr_t(i,m) + ndep_no3_yr(i,m,iccp1) * barefrac
          bnf_free_yr_t(i,m) = bnf_free_yr_t(i,m) + bnf_free_yr(i,m,iccp1) * barefrac
          bnf_tot_yr_t(i,m) = bnf_tot_yr_t(i,m) + bnf_tot_yr(i,m,iccp1) * barefrac
          nitrif_yr_t(i,m) = nitrif_yr_t(i,m) + nitrif_yr(i,m,iccp1) * barefrac
          no_nit_yr_t(i,m) = no_nit_yr_t(i,m) + no_nit_yr(i,m,iccp1) * barefrac
          no_denit_yr_t(i,m) = no_denit_yr_t(i,m) + no_denit_yr(i,m,iccp1) * barefrac
          no_nitdenit_yr_t(i,m) = no_nitdenit_yr_t(i,m) + no_nitdenit_yr(i,m,iccp1) * barefrac
          n2o_nit_yr_t(i,m) = n2o_nit_yr_t(i,m) + n2o_nit_yr(i,m,iccp1) * barefrac
          n2o_denit_yr_t(i,m) = n2o_denit_yr_t(i,m) + n2o_denit_yr(i,m,iccp1) * barefrac
          n2o_nitdenit_yr_t(i,m) = n2o_nitdenit_yr_t(i,m) + n2o_nitdenit_yr(i,m,iccp1) * barefrac
          n2_denit_yr_t(i,m) = n2_denit_yr_t(i,m) + n2_denit_yr(i,m,iccp1) * barefrac
          nvol_yr_t(i,m) = nvol_yr_t(i,m) + nvol_yr(i,m,iccp1) * barefrac
          nleach_yr_t(i,m) = nleach_yr_t(i,m) + nleach_yr(i,m,iccp1) * barefrac
          nhumtrs_yr_t(i,m) = nhumtrs_yr_t(i,m) + nhumtrs_yr(i,m,iccp1) * barefrac
          nmineral_litr_yr_t(i,m) = nmineral_litr_yr_t(i,m) &
                                    + nmineral_litr_yr(i,m,iccp1) * barefrac
          nmineral_humus_yr_t(i,m) = nmineral_humus_yr_t(i,m) &
                                     + nmineral_humus_yr(i,m,iccp1) * barefrac
          netnmineral_yr_t(i,m) = netnmineral_yr_t(i,m) &
                                     + netnmineral_yr(i,m,iccp1) * barefrac
          nimmobil_nh4_yr_t(i,m) = nimmobil_nh4_yr_t(i,m) + nimmobil_nh4_yr(i,m,iccp1) * barefrac
          nimmobil_no3_yr_t(i,m) = nimmobil_no3_yr_t(i,m) + nimmobil_no3_yr(i,m,iccp1) * barefrac
          fNnetland_yr_t(i,m) = fNnetland_yr_t(i,m) + fNnetland_yr(i,m,iccp1) * barefrac
          if (nlitrmass_yr_t(i,m) /= 0.0) c2n_litr_yr_t(i,m) = sum(litrmass_yr_t(i,m,:)) / nlitrmass_yr_t(i,m)
          if (soilnmas_yr_t(i,m) /= 0.0) c2n_humus_yr_t(i,m) = sum(soilcmas_yr_t(i,m,:)) / soilnmas_yr_t(i,m)
        end if !for Ncycle_on

        !> Add values to the per gridcell vars
        laimaxg_yr_g(i) = laimaxg_yr_g(i) + laimaxg_yr_t(i,m) * FAREROT(i,m)
        gleafmas_yr_g(i) = gleafmas_yr_g(i) + gleafmas_yr_t(i,m) * FAREROT(i,m)
        gleafmas_NS_yr_g(i) = gleafmas_NS_yr_g(i) + gleafmas_NS_yr_t(i,m) * FAREROT(i,m)
        gleafmass_yr_g(i) = gleafmass_yr_g(i) + gleafmass_yr_t(i,m) * FAREROT(i,m)
        bleafmas_yr_g(i) = bleafmas_yr_g(i) + bleafmas_yr_t(i,m) * FAREROT(i,m)
        stemmass_yr_g(i) = stemmass_yr_g(i) + stemmass_yr_t(i,m) * FAREROT(i,m)
        stemmass_NS_yr_g(i) = stemmass_NS_yr_g(i) + stemmass_NS_yr_t(i,m) * FAREROT(i,m)
        stemmasss_yr_g(i) = stemmasss_yr_g(i) + stemmasss_yr_t(i,m) * FAREROT(i,m)
        rootmass_yr_g(i) = rootmass_yr_g(i) + rootmass_yr_t(i,m) * FAREROT(i,m)
        rootmass_NS_yr_g(i) = rootmass_NS_yr_g(i) + rootmass_NS_yr_t(i,m) * FAREROT(i,m)
        rootmasss_yr_g(i) = rootmasss_yr_g(i) + rootmasss_yr_t(i,m) * FAREROT(i,m)
        do k = 1,ignd
          litrmass_yr_g(i,k) = litrmass_yr_g(i,k) + litrmass_yr_t(i,m,k) * FAREROT(i,m)
          soilcmas_yr_g(i,k) = soilcmas_yr_g(i,k) + soilcmas_yr_t(i,m,k) * FAREROT(i,m)
          litres_yr_g(i,k) = litres_yr_g(i,k) + litres_yr_t(i,m,k) * FAREROT(i,m)
          soilcres_yr_g(i,k) = soilcres_yr_g(i,k) + soilcres_yr_t(i,m,k) * FAREROT(i,m)
        end do
        vgbiomas_yr_g(i) = vgbiomas_yr_g(i) + vgbiomas_yr_t(i,m) * FAREROT(i,m)
        totcmass_yr_g(i) = totcmass_yr_g(i) + totcmass_yr_t(i,m) * FAREROT(i,m)
        npp_yr_g(i) = npp_yr_g(i) + npp_yr_t(i,m) * FAREROT(i,m)
        gpp_yr_g(i) = gpp_yr_g(i) + gpp_yr_t(i,m) * FAREROT(i,m)
        vcmax0_yr_g(i) = vcmax0_yr_g(i) + vcmax0_yr_t(i,m) * FAREROT(i,m)
        leafns2s_yr_g(i) = leafns2s_yr_g(i) + leafns2s_yr_t(i,m) * FAREROT(i,m)
        stemns2s_yr_g(i) = stemns2s_yr_g(i) + stemns2s_yr_t(i,m) * FAREROT(i,m)
        rootns2s_yr_g(i) = rootns2s_yr_g(i) + rootns2s_yr_t(i,m) * FAREROT(i,m)
        re_alloc_s2l_yr_g(i) = re_alloc_s2l_yr_g(i) + re_alloc_s2l_yr_t(i,m) * FAREROT(i,m)
        re_alloc_r2l_yr_g(i) = re_alloc_r2l_yr_g(i) + re_alloc_r2l_yr_t(i,m) * FAREROT(i,m)
        re_alloc_sr2l_yr_g(i) = re_alloc_sr2l_yr_g(i) + re_alloc_sr2l_yr_t(i,m) * FAREROT(i,m)
        nep_yr_g(i) = nep_yr_g(i) + nep_yr_t(i,m) * FAREROT(i,m)
        nbp_yr_g(i) = nbp_yr_g(i) + nbp_yr_t(i,m) * FAREROT(i,m)
        emit_co2_yr_g(i) = emit_co2_yr_g(i) + emit_co2_yr_t(i,m) * FAREROT(i,m)
        emit_co_yr_g(i) = emit_co_yr_g(i) + emit_co_yr_t(i,m) * FAREROT(i,m)
        emit_ch4_yr_g(i) = emit_ch4_yr_g(i) + emit_ch4_yr_t(i,m) * FAREROT(i,m)
        emit_nmhc_yr_g(i) = emit_nmhc_yr_g(i) + emit_nmhc_yr_t(i,m) * FAREROT(i,m)
        emit_h2_yr_g(i) = emit_h2_yr_g(i) + emit_h2_yr_t(i,m) * FAREROT(i,m)
        emit_nox_yr_g(i) = emit_nox_yr_g(i) + emit_nox_yr_t(i,m) * FAREROT(i,m)
        emit_n2o_yr_g(i) = emit_n2o_yr_g(i) + emit_n2o_yr_t(i,m) * FAREROT(i,m)
        emit_nh3_yr_g(i) = emit_nh3_yr_g(i) + emit_nh3_yr_t(i,m) * FAREROT(i,m)
        emit_pm25_yr_g(i) = emit_pm25_yr_g(i) + emit_pm25_yr_t(i,m) * FAREROT(i,m)
        emit_tpm_yr_g(i) = emit_tpm_yr_g(i) + emit_tpm_yr_t(i,m) * FAREROT(i,m)
        emit_tc_yr_g(i) = emit_tc_yr_g(i) + emit_tc_yr_t(i,m) * FAREROT(i,m)
        emit_oc_yr_g(i) = emit_oc_yr_g(i) + emit_oc_yr_t(i,m) * FAREROT(i,m)
        emit_bc_yr_g(i) = emit_bc_yr_g(i) + emit_bc_yr_t(i,m) * FAREROT(i,m)
        hetrores_yr_g(i) = hetrores_yr_g(i) + hetrores_yr_t(i,m) * FAREROT(i,m)
        autores_yr_g(i) = autores_yr_g(i) + autores_yr_t(i,m) * FAREROT(i,m)
        rmrveg_yr_g(i) = rmrveg_yr_g(i) + rmrveg_yr_t(i,m) * FAREROT(i,m)
        litrfall_yr_g(i) = litrfall_yr_g(i) + litrfall_yr_t(i,m) * FAREROT(i,m)
        burnfrac_yr_g(i) = burnfrac_yr_g(i) + burnfrac_yr_t(i,m) * FAREROT(i,m)
        smfuncveg_yr_g(i) = smfuncveg_yr_g(i) + smfuncveg_yr_t(i,m) * FAREROT(i,m)
        bterm_yr_g(i) = bterm_yr_g(i) + bterm_yr_t(i,m) * FAREROT(i,m)
        lterm_yr_g(i) = lterm_yr_g(i) + lterm_yr_t(i,m) * FAREROT(i,m)
        mterm_yr_g(i) = mterm_yr_g(i) + mterm_yr_t(i,m) * FAREROT(i,m)
        luc_emc_yr_g(i) = luc_emc_yr_g(i) + luc_emc_yr_t(i,m) * FAREROT(i,m)
        lucsocin_yr_g(i) = lucsocin_yr_g(i) + lucsocin_yr_t(i,m) * FAREROT(i,m)
        lucltrin_yr_g(i) = lucltrin_yr_g(i) + lucltrin_yr_t(i,m) * FAREROT(i,m)
        tileAge_yr_g(i) = tileAge_yr_g(i) + tileAge_yr_t(i,m) * FAREROT(i,m)
        luc_emcn_yr_g(i) = luc_emcn_yr_g(i) + luc_emcn_yr_t(i,m) * FAREROT(i,m)
        lucsocinn_yr_g(i) = lucsocinn_yr_g(i) + lucsocinn_yr_t(i,m) * FAREROT(i,m)
        lucltrinn_yr_g(i) = lucltrinn_yr_g(i) + lucltrinn_yr_t(i,m) * FAREROT(i,m)
        ch4WetSpec_yr_g(i) = ch4WetSpec_yr_g(i) + ch4WetSpec_yr_t(i,m) * FAREROT(i,m)
        wetfdyn_yr_g(i) = wetfdyn_yr_g(i) + wetfdyn_yr_t(i,m) * FAREROT(i,m)
        ch4WetDyn_yr_g(i) = ch4WetDyn_yr_g(i) + ch4WetDyn_yr_t(i,m) * FAREROT(i,m)
        ch4soills_yr_g(i) = ch4soills_yr_g(i) + ch4soills_yr_t(i,m) * FAREROT(i,m)
        veghght_yr_g(i) = veghght_yr_g(i) + veghght_yr_t(i,m) * FAREROT(i,m)

        ! Including the LUC product pools. They are per tile values and
        ! are assumed to occupy the whole tile.
        fProductDecomp_yr_g(i) = fProductDecomp_yr_g(i) &
                                 + fProductDecomp_yr_t(i,m) * FAREROT(i,m)
        cProduct_yr_g(i) = cProduct_yr_g(i) &
                           + (litrmassrow(i,m,iccp2,1) &
                           + soilcmasrow(i,m,iccp2,1)) * FAREROT(i,m)
        if (Ncycle_on) then
          nh4_mass_yr_g(i) = nh4_mass_yr_g(i) + nh4_mass_yr_t(i,m) * FAREROT(i,m)
          no3_mass_yr_g(i) = no3_mass_yr_g(i) + no3_mass_yr_t(i,m) * FAREROT(i,m)
          nlitrmass_yr_g(i) = nlitrmass_yr_g(i) + nlitrmass_yr_t(i,m) * FAREROT(i,m)
          soilnmas_yr_g(i) = soilnmas_yr_g(i) + soilnmas_yr_t(i,m) * FAREROT(i,m)
          ngleafmas_yr_g(i) = ngleafmas_yr_g(i) + ngleafmas_yr_t(i,m) * FAREROT(i,m)
          ngleafmas_NS_yr_g(i) = ngleafmas_NS_yr_g(i) + ngleafmas_NS_yr_t(i,m) * FAREROT(i,m)
          ngleafmass_yr_g(i) = ngleafmass_yr_g(i) + ngleafmass_yr_t(i,m) * FAREROT(i,m)
          nbleafmas_yr_g(i) = nbleafmas_yr_g(i) + nbleafmas_yr_t(i,m) * FAREROT(i,m)
          nstemmass_yr_g(i) = nstemmass_yr_g(i) + nstemmass_yr_t(i,m) * FAREROT(i,m)
          nstemmass_NS_yr_g(i) = nstemmass_NS_yr_g(i) + nstemmass_NS_yr_t(i,m) * FAREROT(i,m)
          nstemmasss_yr_g(i) = nstemmasss_yr_g(i) + nstemmasss_yr_t(i,m) * FAREROT(i,m)
          nrootmass_yr_g(i) = nrootmass_yr_g(i) + nrootmass_yr_t(i,m) * FAREROT(i,m)
          nrootmass_NS_yr_g(i) = nrootmass_NS_yr_g(i) + nrootmass_NS_yr_t(i,m) * FAREROT(i,m)
          nrootmasss_yr_g(i) = nrootmasss_yr_g(i) + nrootmasss_yr_t(i,m) * FAREROT(i,m)
          nvgbiomas_yr_g(i) = nvgbiomas_yr_g(i) + nvgbiomas_yr_t(i,m) * FAREROT(i,m)
          bnf_tot_yr_g(i) = bnf_tot_yr_g(i) + bnf_tot_yr_t(i,m) * FAREROT(i,m)
          bnf_free_yr_g(i) = bnf_free_yr_g(i) + bnf_free_yr_t(i,m) * FAREROT(i,m)
          bnf_ant_yr_g(i) = bnf_ant_yr_g(i) + bnf_ant_yr_t(i,m) * FAREROT(i,m)
          bnf_nat_yr_g(i) = bnf_nat_yr_g(i) + bnf_nat_yr_t(i,m) * FAREROT(i,m)
          nitrif_yr_g(i) = nitrif_yr_g(i) + nitrif_yr_t(i,m) * FAREROT(i,m)
          no_nit_yr_g(i) = no_nit_yr_g(i) + no_nit_yr_t(i,m) * FAREROT(i,m)
          no_denit_yr_g(i) = no_denit_yr_g(i) + no_denit_yr_t(i,m) * FAREROT(i,m)
          no_nitdenit_yr_g(i) = no_nitdenit_yr_g(i) + no_nitdenit_yr_t(i,m) * FAREROT(i,m)
          n2o_nit_yr_g(i) = n2o_nit_yr_g(i) + n2o_nit_yr_t(i,m) * FAREROT(i,m)
          n2o_denit_yr_g(i) = n2o_denit_yr_g(i) + n2o_denit_yr_t(i,m) * FAREROT(i,m)
          n2o_nitdenit_yr_g(i) = n2o_nitdenit_yr_g(i) + n2o_nitdenit_yr_t(i,m) * FAREROT(i,m)
          n2_denit_yr_g(i) = n2_denit_yr_g(i) + n2_denit_yr_t(i,m) * FAREROT(i,m)
          nvol_yr_g(i) = nvol_yr_g(i) + nvol_yr_t(i,m) * FAREROT(i,m)
          nleach_yr_g(i) = nleach_yr_g(i) + nleach_yr_t(i,m) * FAREROT(i,m)
          appl_fert_yr_g(i) = appl_fert_yr_g(i) + appl_fert_yr_t(i,m) * FAREROT(i,m)
          ndep_nh4_yr_g(i) = ndep_nh4_yr_g(i) + ndep_nh4_yr_t(i,m) * FAREROT(i,m)
          ndep_no3_yr_g(i) = ndep_no3_yr_g(i) + ndep_no3_yr_t(i,m) * FAREROT(i,m)
          ndemand_wp_npp_yr_g(i) = ndemand_wp_npp_yr_g(i) + ndemand_wp_npp_yr_t(i,m) * FAREROT(i,m)
          nuptake_p_nh4_yr_g(i) = nuptake_p_nh4_yr_g(i) + nuptake_p_nh4_yr_t(i,m) * FAREROT(i,m)
          nuptake_p_no3_yr_g(i) = nuptake_p_no3_yr_g(i) + nuptake_p_no3_yr_t(i,m) * FAREROT(i,m)
          nuptake_a_actl_nh4_yr_g(i) = nuptake_a_actl_nh4_yr_g(i)+ nuptake_a_actl_nh4_yr_t(i,m) * FAREROT(i,m)
          nuptake_a_actl_no3_yr_g(i) = nuptake_a_actl_no3_yr_g(i)+ nuptake_a_actl_no3_yr_t(i,m) * FAREROT(i,m)
          nuptake_yr_g(i) = nuptake_yr_g(i)+ nuptake_yr_t(i,m) * FAREROT(i,m)
          nalloc_l_yr_g(i) = nalloc_l_yr_g(i) + nalloc_l_yr_t(i,m) * FAREROT(i,m)
          nalloc_s_yr_g(i) = nalloc_s_yr_g(i) + nalloc_s_yr_t(i,m) * FAREROT(i,m)
          nalloc_r_yr_g(i) = nalloc_r_yr_g(i) + nalloc_r_yr_t(i,m) * FAREROT(i,m)
          nresorped_s_yr_g(i) = nresorped_s_yr_g(i) + nresorped_s_yr_t(i,m) * FAREROT(i,m)
          nresorped_r_yr_g(i) = nresorped_r_yr_g(i) + nresorped_r_yr_t(i,m) * FAREROT(i,m)
          nre_alloc_s2l_yr_g(i) = nre_alloc_s2l_yr_g(i)+ nre_alloc_s2l_yr_t(i,m) * FAREROT(i,m)
          nre_alloc_r2l_yr_g(i) = nre_alloc_r2l_yr_g(i)+ nre_alloc_r2l_yr_t(i,m) * FAREROT(i,m)
          nleafns2s_yr_g(i) = nleafns2s_yr_g(i) + nleafns2s_yr_t(i,m) * FAREROT(i,m)
          nstemns2s_yr_g(i) = nstemns2s_yr_g(i) + nstemns2s_yr_t(i,m) * FAREROT(i,m)
          nrootns2s_yr_g(i) = nrootns2s_yr_g(i) + nrootns2s_yr_t(i,m) * FAREROT(i,m)
          nlitr_l_yr_g(i) = nlitr_l_yr_g(i) + nlitr_l_yr_t(i,m) * FAREROT(i,m)
          nlitr_s_yr_g(i) = nlitr_s_yr_g(i) + nlitr_s_yr_t(i,m) * FAREROT(i,m)
          nlitr_r_yr_g(i) = nlitr_r_yr_g(i) + nlitr_r_yr_t(i,m) * FAREROT(i,m)
          nlitr_yr_g(i) = nlitr_yr_g(i) + nlitr_yr_t(i,m) * FAREROT(i,m)
          gl2bl_grass_nflux_yr_g(i) = gl2bl_grass_nflux_yr_g(i) &
                                      + gl2bl_grass_nflux_yr_t(i,m) * FAREROT(i,m)
          nhumtrs_yr_g(i) = nhumtrs_yr_g(i) + nhumtrs_yr_t(i,m) * FAREROT(i,m)
          nmineral_litr_yr_g(i) = nmineral_litr_yr_g(i) + nmineral_litr_yr_t(i,m) * FAREROT(i,m)
          nmineral_humus_yr_g(i) = nmineral_humus_yr_g(i) + nmineral_humus_yr_t(i,m) * FAREROT(i,m)
          netnmineral_yr_g(i) = netnmineral_yr_g(i) + netnmineral_yr_t(i,m) * FAREROT(i,m)
          nimmobil_nh4_yr_g(i) = nimmobil_nh4_yr_g(i) + nimmobil_nh4_yr_t(i,m) * FAREROT(i,m)
          nimmobil_no3_yr_g(i) = nimmobil_no3_yr_g(i) + nimmobil_no3_yr_t(i,m) * FAREROT(i,m)
          fNnetland_yr_g(i) = fNnetland_yr_g(i) + fNnetland_yr_t(i,m) * FAREROT(i,m)
          nProduct_yr_g(i) = nProduct_yr_g(i) &
                             + (nlitrmassrow(i,m,iccp2) &
                             + soilnmasrow(i,m,iccp2)) * convertg2kg * FAREROT(i,m)
        end if !for Ncycle_on

      end do ! loop 900 ! m

      if (Ncycle_on) then
        if ((ngleafmas_yr_g(i) + nbleafmas_yr_g(i)) /= 0.0) c2n_l_yr_g(i) = &
            (gleafmas_yr_g(i) + bleafmas_yr_g(i)) / (ngleafmas_yr_g(i) + nbleafmas_yr_g(i))
        if (nstemmass_yr_g(i) /= 0.0) c2n_s_yr_g(i) = stemmass_yr_g(i) / nstemmass_yr_g(i)
        if (nrootmass_yr_g(i) /= 0.0) c2n_r_yr_g(i) = rootmass_yr_g(i) / nrootmass_yr_g(i)
        if ((ngleafmas_yr_g(i) + nbleafmas_yr_g(i) + nstemmass_yr_g(i) + nrootmass_yr_g(i)) /= 0.0) &
            c2n_wp_yr_g(i) = (gleafmas_yr_g(i) + bleafmas_yr_g(i) + stemmass_yr_g(i) + rootmass_yr_g(i)) &
                              / (ngleafmas_yr_g(i) + nbleafmas_yr_g(i) + nstemmass_yr_g(i) + nrootmass_yr_g(i))
        if (nlitrmass_yr_g(i) /= 0.0) c2n_litr_yr_g(i) = sum(litrmass_yr_g(i,:)) / nlitrmass_yr_g(i)
        if (soilnmas_yr_g(i) /= 0.0) c2n_humus_yr_g(i) = sum(soilcmas_yr_g(i,:)) / soilnmas_yr_g(i)
      end if

      !> Write to annual output files:

      ! Prepare the timestamp for this year
      timeStamp = consecDays

      !> First write out the per gridcell values
      call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr_g' ,timeStamp,'lai', [laimaxg_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr_g',timeStamp,'cVeg',[vgbiomas_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_yr_g',timeStamp,'cLeaf',[gleafmas_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_NS_yr_g',timeStamp,'cLeaf_ns',[gleafmas_NS_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmass_yr_g',timeStamp,'cLeaf_s',[gleafmass_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr_g',timeStamp,'cStem',[stemmass_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_NS_yr_g',timeStamp,'cStem_ns',[stemmass_NS_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmasss_yr_g',timeStamp,'cStem_s',[stemmasss_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr_g',timeStamp,'cRoot',[rootmass_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_NS_yr_g',timeStamp,'cRoot_ns',[rootmass_NS_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmasss_yr_g',timeStamp,'cRoot_s',[rootmasss_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr_g',timeStamp,'cLand',[totcmass_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr_g'     ,timeStamp,'npp',[npp_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr_g'     ,timeStamp,'gpp',[gpp_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'leafns2s_yr_g',timeStamp,'leafns2s',[leafns2s_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'stemns2s_yr_g',timeStamp,'stemns2s',[stemns2s_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rootns2s_yr_g',timeStamp,'rootns2s',[rootns2s_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_s2l_yr_g',timeStamp,'realloc_s2l',[re_alloc_s2l_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_r2l_yr_g',timeStamp,'realloc_r2l',[re_alloc_r2l_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_sr2l_yr_g',timeStamp,'realloc_sr2l',[re_alloc_sr2l_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr_g'     ,timeStamp,'nep',[nep_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr_g'     ,timeStamp,'nbp',[nbp_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr_g',timeStamp,'rh',[hetrores_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr_g' ,timeStamp,'ra',[autores_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rmrveg_yr_g' ,timeStamp,'rmrveg',[rmrveg_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_yr_g',timeStamp,'fVegLitter',[litrfall_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr_g' ,timeStamp,'vegHeight',[veghght_yr_g(i)])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_yr_g' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_yr_g' ,timeStamp,'wetlandFrac',[wetfdyn_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_yr_g' ,timeStamp,'soilCH4cons',[ch4soills_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'vcmax0_yr_g',timeStamp,'vcmax0',[vcmax0_yr_g(i)])
      if (Ncycle_on) then
        call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_yr_g',timeStamp,'bnf_tot',[bnf_tot_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_yr_g',timeStamp,'bnf_free',[bnf_free_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_yr_g',timeStamp,'bnf_ant',[bnf_ant_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_yr_g',timeStamp,'bnf_nat',[bnf_nat_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_yr_g',timeStamp,'nLeaf',[ngleafmas_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_NS_yr_g',timeStamp,'nLeaf_ns',[ngleafmas_NS_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmass_yr_g',timeStamp,'nLeaf_s',[ngleafmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_yr_g',timeStamp,'nbLeaf',[nbleafmas_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_yr_g',timeStamp,'nStem',[nstemmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_NS_yr_g',timeStamp,'nStem_ns',[nstemmass_NS_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmasss_yr_g',timeStamp,'nStem_s',[nstemmasss_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_yr_g',timeStamp,'nRoot',[nrootmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_NS_yr_g',timeStamp,'nRoot_ns',[nrootmass_NS_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmasss_yr_g',timeStamp,'nRoot_s',[nrootmasss_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ndemand_wp_npp_yr_g',timeStamp,'ndemand_npp',[ndemand_wp_npp_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_yr_g',timeStamp,'c2n_l', [c2n_l_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_s_yr_g',timeStamp,'c2n_s', [c2n_s_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_r_yr_g',timeStamp,'c2n_r', [c2n_r_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_wp_yr_g',timeStamp,'c2n_wp', [c2n_wp_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_litr_yr_g',timeStamp,'c2n_litr', [c2n_litr_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_humus_yr_g',timeStamp,'c2n_humus', [c2n_humus_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_yr_g',timeStamp,'nLitter',[nlitrmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_yr_g',timeStamp,'nSoil',[soilnmas_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_yr_g',timeStamp,'nVeg',[nvgbiomas_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_yr_g',timeStamp,'nh4',[nh4_mass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_yr_g',timeStamp,'no3',[no3_mass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_yr_g',timeStamp,'nitrif',[nitrif_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nit_yr_g',timeStamp,'no_nit',[no_nit_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_yr_g',timeStamp,'no_denit',[no_denit_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_yr_g',timeStamp,'no',[no_nitdenit_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nit_yr_g',timeStamp,'n2o_nit',[n2o_nit_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_yr_g',timeStamp,'n2o_denit',[n2o_denit_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_yr_g',timeStamp,'n2o',[n2o_nitdenit_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_yr_g',timeStamp,'n2',[n2_denit_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_yr_g',timeStamp,'nvol',[nvol_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_yr_g',timeStamp,'nleach',[nleach_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'appl_fert_yr_g',timeStamp,'nfer_nh4',[appl_fert_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_nh4_yr_g',timeStamp,'ndep_nh4',[ndep_nh4_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_no3_yr_g',timeStamp,'ndep_no3',[ndep_no3_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_nh4_yr_g',timeStamp,'nuptake_p_nh4',[nuptake_p_nh4_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_no3_yr_g',timeStamp,'nuptake_p_no3',[nuptake_p_no3_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_nh4_yr_g',timeStamp,'nuptake_a_actl_nh4',[nuptake_a_actl_nh4_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_no3_yr_g',timeStamp,'nuptake_a_actl_no3',[nuptake_a_actl_no3_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_yr_g',timeStamp,'nuptake',[nuptake_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_l_yr_g',timeStamp,'nalloc_l', [nalloc_l_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_s_yr_g',timeStamp,'nalloc_s', [nalloc_s_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_r_yr_g',timeStamp,'nalloc_r', [nalloc_r_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_yr_g',timeStamp,'nresorped_s', [nresorped_s_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_yr_g',timeStamp,'nresorped_r', [nresorped_r_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_s2l_yr_g',timeStamp,'nre_alloc_s2l', [nre_alloc_s2l_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_r2l_yr_g',timeStamp,'nre_alloc_r2l', [nre_alloc_r2l_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nleafns2s_yr_g',timeStamp,'nleafns2s', [nleafns2s_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemns2s_yr_g',timeStamp,'nstemns2s', [nstemns2s_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootns2s_yr_g',timeStamp,'nrootns2s', [nrootns2s_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_yr_g',timeStamp,'nlitr_l', [nlitr_l_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_s_yr_g',timeStamp,'nlitr_s', [nlitr_s_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_r_yr_g',timeStamp,'nlitr_r', [nlitr_r_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_yr_g',timeStamp,'nlitr', [nlitr_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gl2bl_grass_nflux_yr_g' ,timeStamp,'gl2bl_grass_nflux', [gl2bl_grass_nflux_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nhumtrs_yr_g',timeStamp,'nhumtrs', [nhumtrs_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_litr_yr_g',timeStamp,'nmineral_litr',[nmineral_litr_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_humus_yr_g',timeStamp,'nmineral_humus',[nmineral_humus_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'netnmineral_yr_g',timeStamp,'netnmineral',[netnmineral_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_nh4_yr_g',timeStamp,'nimmobil_nh4',[nimmobil_nh4_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_no3_yr_g',timeStamp,'nimmobil_no3',[nimmobil_no3_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fNnetland_yr_g',timeStamp,'fNnetland',[fNnetland_yr_g(i)])
      end if ! Ncycle

      call writeOutput1D(lonLocalIndex,latLocalIndex,'peatdep_yr_g' ,timeStamp,'peatdep', [peatdep_yr_g(i)]) !s.r.c added peat depth output
      call writeOutput1D(lonLocalIndex,latLocalIndex,'peatSoilC_yr_g' ,timeStamp,'peatSoilC', [peatSoilC_yr_g(i)]) !s.r.c added peat carbon output
      
      call writeOutput1D(lonLocalIndex,latLocalIndex,'litrMassPerLay_yr_g',timeStamp,'cLitterperlay',[litrmass_yr_g(i,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'soilCMasPerLay_yr_g',timeStamp,'cSoilperlay',[soilcmas_yr_g(i,:)])

      ! Make the bulk litter and soil C pool and respiration temporary variables:
      bulkLitterCarbon_yr_g(1) = sum(litrmass_yr_g(i,:))
      bulkSoilCarbon_yr_g(1) = sum(soilcmas_yr_g(i,:))
      bulkLitterResp_yr_g(1) = sum(litres_yr_g(i,:))
      bulkSoilResp_yr_g(1) = sum(soilcres_yr_g(i,:))
      call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_g',timeStamp,'cLitter',[bulkLitterCarbon_yr_g])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_g',timeStamp,'cSoil',[bulkSoilCarbon_yr_g])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_g'  ,timeStamp,'rhLitter',[bulkLitterResp_yr_g])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_g',timeStamp,'rhSoil',[bulkSoilResp_yr_g])

      if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_yr_g' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_yr_g(i)])
      end if

      if (PFTCompetition) then
        do m = 1,nmtest
          sumfare = 0.0
          fcancmxNoSeed = 0.0
          pftExist = 0.0
          do j = 1,icc
            if (pftexistrow(i,1,j)) then
              pftExist(j) = 1.0
            end if
              sumfare = sumfare + fcancmxrow(i,m,j)
              fcancmxNoSeed(j) = fcancmxrow(i,m,j)
          end do
          call writeOutput1D(lonLocalIndex,latLocalIndex,'pftexistrow_yr_g' ,timeStamp,'landCoverExist',[pftExist])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_yr_g' ,timeStamp,'landCoverFrac',[fcancmxNoSeed(1:icc),1 - sumfare])
          if (doperpftoutput) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'cc_yr' ,timeStamp,'cc',[cc_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'mm_yr' ,timeStamp,'mm',[mm_yr(i,m,:)])
          end if  
        end do
        call writeOutput1D(lonLocalIndex,latLocalIndex,'twarmm_yr_g' ,timeStamp,'twarmm', [twarmmrow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tcoldm_yr_g' ,timeStamp,'tcoldm', [tcoldmrow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gdd5_yr_g' ,timeStamp,'gdd5', [gdd5row(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'aridity_yr_g' ,timeStamp,'aridity', [aridityrow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'srplsmon_yr_g' ,timeStamp,'srplsmon', [srplsmonrow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'defctmon_yr_g' ,timeStamp,'defctmon', [defctmonrow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'anndefct_yr_g' ,timeStamp,'anndefct', [anndefctrow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'annsrpls_yr_g' ,timeStamp,'annsrpls', [annsrplsrow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'annpcp_yr_g' ,timeStamp,'annpcp', [annpcprow(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'dry_season_length_yr_g' ,timeStamp,'dry_season_length', [dry_season_lengthrow(i)])
      else
        do m = 1,nmtest
          sumfare = 0.0
          do j = 1,icc
            sumfare = sumfare + fcancmxrow(i,m,j)
          end do ! j
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_yr_g' ,timeStamp,'landCoverFrac',[fcancmxrow(i,m,1:icc),1 - sumfare])
        end do ! m
      end if

      if (dofire .or. lnduseon .or. timberHarvest .or. prescribedFire) then
        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_yr_g' ,timeStamp,'fFireCH4',[emit_ch4_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr_g' ,timeStamp,'fFire',[emit_co2_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_nox_yr_g' ,timeStamp,'fFireNOX',[emit_nox_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_nh3_yr_g' ,timeStamp,'fFireNH3',[emit_nh3_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_n2o_yr_g' ,timeStamp,'fFireN2O',[emit_n2o_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr_g' ,timeStamp,'burntFractionAll',[burnfrac_yr_g(i)])

        call writeOutput1D(lonLocalIndex,latLocalIndex,'bterm_yr_g' ,timeStamp,'bterm',[bterm_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'mterm_yr_g' ,timeStamp,'mterm',[mterm_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'lterm_yr_g' ,timeStamp,'lterm',[lterm_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'smfuncveg_yr_g' ,timeStamp,'smfuncveg',[smfuncveg_yr_g(i)])
      end if

      if (dynamicTilingOn .or. trackTileAge) call writeOutput1D(lonLocalIndex,latLocalIndex,'tileAge_yr_g' ,timeStamp,'tileAge',[tileAge_yr_g(i)])
      if (dynamicTilingOn) call writeOutput1D(lonLocalIndex,latLocalIndex,'FAREROT_yr' ,timeStamp,'FARE',[sum(FAREROT(i,:))]) !useful diagnostic
      if (timberHarvest) call writeOutput1D(lonLocalIndex,latLocalIndex,'timharvarea_yr_g' ,timeStamp,'timharvarea',[sum(timharvarea_yr_t(i,:)*FAREROT(i,:))]) !timber harvest area can be written directly

      if (lnduseon .or. timberHarvest) then
        call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_yr_g' ,timeStamp,'fDeforestToAtmos',[luc_emc_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_yr_g' ,timeStamp,'fDeforestToLitter',[lucltrin_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_yr_g' ,timeStamp,'fDeforestToSoil',[lucsocin_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_yr_g' ,timeStamp,'fDeforestTotal', &
                           [lucsocin_yr_g(i) + lucltrin_yr_g(i) + luc_emc_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fProductDecomp_yr_g' ,timeStamp,'fProductDecomp',[fProductDecomp_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'cProduct_yr_g' ,timeStamp,'cProduct',[cProduct_yr_g(i)])
        if (Ncycle_on) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nProduct_yr_g' ,timeStamp,'nProduct',[nProduct_yr_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emcn_yr_g' ,timeStamp,'fDeforestToAtmosN',[luc_emcn_yr_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrinn_yr_g' ,timeStamp,'fDeforestToLitterN',[lucltrinn_yr_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocinn_yr_g' ,timeStamp,'fDeforestToSoilN',[lucsocinn_yr_g(i)])
        end if
      end if

      if (doperpftoutput) then
        if (nmtest > 1) then
          print * ,'Per PFT and per tile outputs not implemented yet'
        else

          m = 1 ! FLAG only implemented for composite mode,tiles == 1 !! !

          call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr' ,timeStamp,'lai', [laimaxg_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr',timeStamp,'cVeg',[vgbiomas_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_yr',timeStamp,'cLeaf',[gleafmas_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_NS_yr',timeStamp,'cLeaf_ns',[gleafmas_NS_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmass_yr',timeStamp,'cLeaf_s',[gleafmass_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr',timeStamp,'cStem',[stemmass_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_NS_yr',timeStamp,'cStem_ns',[stemmass_NS_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmasss_yr',timeStamp,'cStem_s',[stemmasss_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr',timeStamp,'cRoot',[rootmass_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_NS_yr',timeStamp,'cRoot_ns',[rootmass_NS_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmasss_yr',timeStamp,'cRoot_s',[rootmasss_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr',timeStamp,'cLand',[totcmass_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr'     ,timeStamp,'npp',[npp_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr'     ,timeStamp,'gpp',[gpp_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr'     ,timeStamp,'nep',[nep_yr(i,m,:)])
          ! NOTE: This NBP does not include LUC product pool contributions since they are
          ! not per PFT but rather per tile
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr'     ,timeStamp,'nbp',[nbp_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr',timeStamp,'rh',[hetrores_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr' ,timeStamp,'ra',[autores_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_yr',timeStamp,'fVegLitter',[litrfall_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr' ,timeStamp,'vegHeight',[veghght_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_s2l_yr',timeStamp,'realloc_s2l',[re_alloc_s2l_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_r2l_yr',timeStamp,'realloc_r2l',[re_alloc_r2l_yr(i,m,:)])

          ! For FireMIP calc the low vs high cover.
          coverLow = 0.
          coverHigh = 0.
          do j = 1, icc
            if (fcancmxrow(i,m,j) > 0.) then
              if (veghght_yr(i,m,j) < 5.) then
                coverLow = coverLow + fcancmxrow(i,m,j)
              else
                coverHigh = coverHigh + fcancmxrow(i,m,j)
              end if
            end if
          end do
          call writeOutput1D(lonLocalIndex,latLocalIndex,'coverLow' ,timeStamp,'coverLow',[coverLow])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'coverHigh' ,timeStamp,'coverHigh',[coverHigh])


          ! Make the bulk litter and soil C pool and respiration temporary variables:
          do k = 1,iccp1
            bulkLitterCarbon_yr(k) = sum(litrmass_yr(i,m,k,:))
            bulkSoilCarbon_yr(k) = sum(soilcmas_yr(i,m,k,:))
            bulkLitterResp_yr(k) = sum(litres_yr(i,m,k,:))
            bulkSoilResp_yr(k) = sum(soilcres_yr(i,m,k,:))
          end do

          call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr',timeStamp,'cLitter',[bulkLitterCarbon_yr])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr',timeStamp,'cSoil',[bulkSoilCarbon_yr])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr'  ,timeStamp,'rhLitter',[bulkLitterResp_yr])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr',timeStamp,'rhSoil',[bulkSoilResp_yr])
        
          call writeOutput1D(lonLocalIndex,latLocalIndex,'vcmax0_yr'  ,timeStamp,'vcmax0',[vcmax0_yr(i,m,:)])

          if (Ncycle_on) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_yr',timeStamp,'nVeg',[nvgbiomas_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_yr',timeStamp,'nh4',[nh4_mass_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_yr',timeStamp,'no3',[no3_mass_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_yr',timeStamp,'nLeaf',[ngleafmas_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_NS_yr',timeStamp,'nLeaf_ns',[ngleafmas_NS_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmass_yr',timeStamp,'nLeaf_s',[ngleafmass_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_yr',timeStamp,'nbLeaf',[nbleafmas_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_yr',timeStamp,'nStem',[nstemmass_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_NS_yr',timeStamp,'nStem_ns',[nstemmass_NS_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmasss_yr',timeStamp,'nStem_s',[nstemmasss_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_yr',timeStamp,'nRoot',[nrootmass_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_NS_yr',timeStamp,'nRoot_ns',[nrootmass_NS_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmasss_yr',timeStamp,'nRoot_s',[nrootmasss_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_yr',timeStamp,'nLitter',[nlitrmass_yr(i,m,1:iccp1)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_yr',timeStamp,'nSoil',[soilnmas_yr(i,m,1:iccp1)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_yr',timeStamp,'bnf_tot',[bnf_tot_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_yr',timeStamp,'bnf_nat',[bnf_nat_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_yr',timeStamp,'bnf_ant',[bnf_ant_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_yr',timeStamp,'bnf_free',[bnf_free_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_nh4_yr',timeStamp,'nuptake_p_nh4',[nuptake_p_nh4_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_no3_yr',timeStamp,'nuptake_p_no3',[nuptake_p_no3_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_nh4_yr',timeStamp,'nuptake_a_actl_nh4',[nuptake_a_actl_nh4_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_no3_yr',timeStamp,'nuptake_a_actl_no3',[nuptake_a_actl_no3_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_yr',timeStamp,'n2o',[n2o_nitdenit_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_yr',timeStamp,'no',[no_nitdenit_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_yr',timeStamp,'nitrif',[nitrif_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'netnmineral_yr',timeStamp,'netnmineral',[netnmineral_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_nh4_yr',timeStamp,'nimmobil_nh4',[nimmobil_nh4_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_no3_yr',timeStamp,'nimmobil_no3',[nimmobil_no3_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_litr_yr',timeStamp,'nmineral_litr',[nmineral_litr_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_humus_yr',timeStamp,'nmineral_humus',[nmineral_humus_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_yr',timeStamp,'nuptake',[nuptake_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_yr',timeStamp,'c2n_l',[c2n_l_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ndemand_wp_npp_yr',timeStamp,'ndemand_npp',[ndemand_wp_npp_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_yr',timeStamp,'nlitr_l',[nlitr_l_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_yr',timeStamp,'no_denit',[no_denit_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_yr',timeStamp,'n2o_denit',[n2o_denit_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_yr',timeStamp,'n2',[n2_denit_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_yr',timeStamp,'nvol',[nvol_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_yr',timeStamp,'nleach',[nleach_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_yr',timeStamp,'nlitr', [nlitr_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nhumtrs_yr',timeStamp,'nhumtrs', [nhumtrs_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_yr',timeStamp,'nresorped_s',[nresorped_s_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_yr',timeStamp,'nresorped_r',[nresorped_r_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_s2l_yr',timeStamp,'nre_alloc_s2l',[nre_alloc_s2l_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_r2l_yr',timeStamp,'nre_alloc_r2l',[nre_alloc_r2l_yr(i,m,:)])
          end if ! Ncycle

          if (dofire .or. lnduseon .or. timberHarvest .or. prescribedFire) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_yr' ,timeStamp,'fFireCH4',[emit_ch4_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr' ,timeStamp,'fFire',[emit_co2_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr' ,timeStamp,'burntFractionAll',[burnfrac_yr(i,m,:)])
          end if

        end if
      end if

      if (dopertileoutput) then

        if (nmtest == 1) then
          print * ,'Switch selected for per tile output but number of tiles is only one.'
        else

          !> Write out the per tile values
          call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr_t' ,timeStamp,'lai', [laimaxg_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr_t',timeStamp,'cVeg',[vgbiomas_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_yr_t',timeStamp,'cLeaf',[gleafmas_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmas_NS_yr_t',timeStamp,'cLeaf_ns',[gleafmas_NS_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gleafmass_yr_t',timeStamp,'cLeaf_s',[gleafmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr_t',timeStamp,'cStem',[stemmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_NS_yr_t',timeStamp,'cStem_ns',[stemmass_NS_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmasss_yr_t',timeStamp,'cStem_s',[stemmasss_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr_t',timeStamp,'cRoot',[rootmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_NS_yr_t',timeStamp,'cRoot_ns',[rootmass_NS_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmasss_yr_t',timeStamp,'cRoot_s',[rootmasss_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr_t',timeStamp,'cLand',[totcmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr_t'     ,timeStamp,'npp',[npp_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr_t'     ,timeStamp,'gpp',[gpp_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'vcmax0_yr_t',timeStamp,'vcmax0',[vcmax0_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'leafns2s_yr_t',timeStamp,'leafns2s',[leafns2s_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemns2s_yr_t',timeStamp,'stemns2s',[stemns2s_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootns2s_yr_t',timeStamp,'rootns2s',[rootns2s_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_s2l_yr_t',timeStamp,'realloc_s2l',[re_alloc_s2l_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_r2l_yr_t',timeStamp,'realloc_r2l',[re_alloc_r2l_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'re_alloc_sr2l_yr_t',timeStamp,'realloc_sr2l',[re_alloc_sr2l_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr_t'     ,timeStamp,'nep',[nep_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr_t'     ,timeStamp,'nbp',[nbp_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr_t',timeStamp,'rh',[hetrores_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr_t' ,timeStamp,'ra',[autores_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rmrveg_yr_t' ,timeStamp,'rmrveg',[rmrveg_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_yr_t',timeStamp,'fVegLitter',[litrfall_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr_t' ,timeStamp,'vegHeight',[veghght_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_yr_t' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_yr_t' ,timeStamp,'wetlandFrac',[wetfdyn_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_yr_t' ,timeStamp,'soilCH4cons',[ch4soills_yr_t(i,:)])

          ! Make the bulk litter and soil C pool and respiration temporary variables:
          do m = 1,nmtest
            bulkLitterCarbon_yr_t(m) = sum(litrmass_yr_t(i,m,:))
            bulkSoilCarbon_yr_t(m) = sum(soilcmas_yr_t(i,m,:))
            bulkLitterResp_yr_t(m) = sum(litres_yr_t(i,m,:))
            bulkSoilResp_yr_t(m) = sum(soilcres_yr_t(i,m,:))
          end do

          call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_t',timeStamp,'cLitter',[bulkLitterCarbon_yr_t])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_t',timeStamp,'cSoil',[bulkSoilCarbon_yr_t])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_t'  ,timeStamp,'rhLitter',[bulkLitterResp_yr_t])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_t',timeStamp,'rhSoil',[bulkSoilResp_yr_t])

          if (Ncycle_on) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_tot_yr_t',timeStamp,'bnf_tot',[bnf_tot_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_free_yr_t',timeStamp,'bnf_free',[bnf_free_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_ant_yr_t',timeStamp,'bnf_ant',[bnf_ant_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'bnf_nat_yr_t',timeStamp,'bnf_nat',[bnf_nat_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nh4_mass_yr_t',timeStamp,'nh4',[nh4_mass_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no3_mass_yr_t',timeStamp,'no3',[no3_mass_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_yr_t',timeStamp,'nLeaf',[ngleafmas_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmas_NS_yr_t',timeStamp,'nLeaf_ns',[ngleafmas_NS_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ngleafmass_yr_t',timeStamp,'nLeaf_s',[ngleafmass_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbleafmas_yr_t',timeStamp,'nbLeaf',[nbleafmas_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_yr_t',timeStamp,'nStem',[nstemmass_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmass_NS_yr_t',timeStamp,'nStem_ns',[nstemmass_NS_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemmasss_yr_t',timeStamp,'nStem_s',[nstemmasss_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_yr_t',timeStamp,'nRoot',[nrootmass_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmass_NS_yr_t',timeStamp,'nRoot_ns',[nrootmass_NS_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootmasss_yr_t',timeStamp,'nRoot_s',[nrootmasss_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ndemand_wp_npp_yr_t',timeStamp,'ndemand_npp',[ndemand_wp_npp_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_l_yr_t',timeStamp,'c2n_l', [c2n_l_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_s_yr_t',timeStamp,'c2n_s', [c2n_s_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_r_yr_t',timeStamp,'c2n_r', [c2n_r_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_wp_yr_t',timeStamp,'c2n_wp', [c2n_wp_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_litr_yr_t',timeStamp,'c2n_litr', [c2n_litr_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'c2n_humus_yr_t',timeStamp,'c2n_humus', [c2n_humus_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitrmass_yr_t',timeStamp,'nLitter',[nlitrmass_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilnmas_yr_t',timeStamp,'nSoil',[soilnmas_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nvgbiomas_yr_t',timeStamp,'nVeg',[nvgbiomas_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nitrif_yr_t',timeStamp,'nitrif',[nitrif_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nit_yr_t',timeStamp,'no_nit',[no_nit_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_denit_yr_t',timeStamp,'no_denit',[no_denit_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'no_nitdenit_yr_t',timeStamp,'no',[no_nitdenit_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nit_yr_t',timeStamp,'n2o_nit',[n2o_nit_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_denit_yr_t',timeStamp,'n2o_denit',[n2o_denit_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2o_nitdenit_yr_t',timeStamp,'n2o',[n2o_nitdenit_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'n2_denit_yr_t',timeStamp,'n2',[n2_denit_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nvol_yr_t',timeStamp,'nvol',[nvol_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nleach_yr_t',timeStamp,'nleach',[nleach_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'appl_fert_yr_t',timeStamp,'nfer_nh4',[appl_fert_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_nh4_yr_t',timeStamp,'ndep_nh4',[ndep_nh4_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ndep_no3_yr_t',timeStamp,'ndep_no3',[ndep_no3_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_nh4_yr_t',timeStamp,'nuptake_p_nh4',[nuptake_p_nh4_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_p_no3_yr_t',timeStamp,'nuptake_p_no3',[nuptake_p_no3_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_nh4_yr_t',timeStamp,'nuptake_a_actl_nh4',[nuptake_a_actl_nh4_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_a_actl_no3_yr_t',timeStamp,'nuptake_a_actl_no3',[nuptake_a_actl_no3_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nuptake_yr_t',timeStamp,'nuptake',[nuptake_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_l_yr_t',timeStamp,'nalloc_l', [nalloc_l_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_s_yr_t',timeStamp,'nalloc_s', [nalloc_s_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nalloc_r_yr_t',timeStamp,'nalloc_r', [nalloc_r_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_s_yr_t',timeStamp,'nresorped_s', [nresorped_s_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nresorped_r_yr_t',timeStamp,'nresorped_r', [nresorped_r_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_s2l_yr_t',timeStamp,'nre_alloc_s2l', [nre_alloc_s2l_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nre_alloc_r2l_yr_t',timeStamp,'nre_alloc_r2l', [nre_alloc_r2l_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nleafns2s_yr_t',timeStamp,'nleafns2s', [nleafns2s_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nstemns2s_yr_t',timeStamp,'nstemns2s', [nstemns2s_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nrootns2s_yr_t',timeStamp,'nrootns2s', [nrootns2s_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_l_yr_t',timeStamp,'nlitr_l', [nlitr_l_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_s_yr_t',timeStamp,'nlitr_s', [nlitr_s_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_r_yr_t',timeStamp,'nlitr_r', [nlitr_r_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nlitr_yr_t',timeStamp,'nlitr', [nlitr_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gl2bl_grass_nflux_yr_t',timeStamp,'gl2bl_grass_nflux', [gl2bl_grass_nflux_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nhumtrs_yr_t',timeStamp,'nhumtrs', [nhumtrs_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_litr_yr_t',timeStamp,'nmineral_litr',[nmineral_litr_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nmineral_humus_yr_t',timeStamp,'nmineral_humus',[nmineral_humus_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'netnmineral_yr_t',timeStamp,'netnmineral',[netnmineral_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_nh4_yr_t',timeStamp,'nimmobil_nh4',[nimmobil_nh4_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nimmobil_no3_yr_t',timeStamp,'nimmobil_no3',[nimmobil_no3_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'fNnetland_yr_t',timeStamp,'fNnetland',[fNnetland_yr_t(i,:)])
          end if ! Ncycle

          call writeOutput1D(lonLocalIndex,latLocalIndex,'peatdep_yr_t' ,timeStamp,'peatdep', [peatdep_yr_t(i,:)]) !s.r.c added peat depth output
          call writeOutput1D(lonLocalIndex,latLocalIndex,'peatSoilC_yr_t' ,timeStamp,'peatSoilC', [peatSoilC_yr_t(i,:)]) !s.r.c added peat carbon output

          if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_yr_t' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_yr_t(i,:)])
          end if
          if (dofire .or. lnduseon .or. timberHarvest .or. prescribedFire) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_yr_t' ,timeStamp,'fFireCH4',[emit_ch4_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr_t' ,timeStamp,'fFire',[emit_co2_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr_t' ,timeStamp,'burntFractionAll',[burnfrac_yr_t(i,:)])
          end if

          if (dynamicTilingOn .or. trackTileAge) call writeOutput1D(lonLocalIndex,latLocalIndex,'tileAge_yr_t' ,timeStamp,'tileAge',[tileAge_yr_t(i,:)])
          if (dynamicTilingOn) call writeOutput1D(lonLocalIndex,latLocalIndex,'FAREROT_yr_t' ,timeStamp,'FARE',[FAREROT(i,:)]) !the final value of FARE can be written directly (averaging FARE throughout the year or across tiles isn't informative)
          if (timberHarvest) call writeOutput1D(lonLocalIndex,latLocalIndex,'timharvarea_yr_t' ,timeStamp,'timharvarea',[timharvarea_yr_t(i,:)]) !timber harvest area can be written directly

          if (lnduseon .or. timberHarvest) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_yr_t' ,timeStamp,'fDeforestToAtmos',[luc_emc_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_yr_t' ,timeStamp,'fDeforestToLitter',[lucltrin_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_yr_t' ,timeStamp,'fDeforestToSoil',[lucsocin_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_yr_t' ,timeStamp,'fDeforestTotal', &
                               [lucsocin_yr_t(i,:) + lucltrin_yr_t(i,:) + luc_emc_yr_t(i,:)])
            if (Ncycle_on) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emcn_yr_t' ,timeStamp,'fDeforestToAtmosN',[luc_emcn_yr_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrinn_yr_t' ,timeStamp,'fDeforestToLitterN',[lucltrinn_yr_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocinn_yr_t' ,timeStamp,'fDeforestToSoilN',[lucsocinn_yr_t(i,:)])
            end if
          end if
        end if
      end if

      !> Reset all annual vars in preparation for next year
      call resetYearEnd(nltest, nmtest)

    end if ! if iday=365/366
  
    end associate
  end subroutine ctem_annual_aw
  !! @}
  !> \namespace prepareoutputs
  !> Central module that handles all CTEM preparation and writing of output files

end module prepareOutputs
