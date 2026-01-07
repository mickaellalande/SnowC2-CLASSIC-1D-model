!> \file
!> Wind sublimation losses calculation for the snowpack on the ground.
!
module snowWindSublimation
  !
  !     * JUL 04/24 - M.LALANDE. IMPLEMENT WIND SUBLIMATION LOSSES FROM GORDON ET AL. (2006).
  !
  use generalutils,  only : calcEsat

  implicit none

  public :: GD06

contains

  !---------------------------------------------------------------------------------------
  !> \ingroup snowWindSublimation_GD06
  !! @{
  !> Compute the blowing snow sublimation losses rate, \f$Q_s\f$ [kg m\f$^{-2}\f$ s\f$^{-1}\f$], following Gordon et al. (2006)'s parameterization:
  !!
  !! $$Q_s=A\left(\frac{T_0}{T_a}\right)^\gamma v_t \rho_a q_{s i}\left(1-R h_i\right)\left(v_a / v_t\right)^B \text {, for } v_a>v_t{,}$$
  !!
  !! where \f$T_0\f$ [K] is the surface temperature, \f$T_a\f$ [K] is the air temperature at the reference height,
  !! \f$v_t\f$ [m s\f$^{-1}\f$] is the threshold wind speed for the initiation of blowing snow,
  !! \f$v_a\f$ [m s\f$^{-1}\f$] is the wind speed at the reference height, \f$\rho_a\f$ [kg m\f$^{-3}\f$] is the density of air, 
  !! \f$q_{s i}\f$ [kg kg\f$^{-1}\f$] is the saturation specific humidity, \f$R h_i\f$ [-] is the relative humidity with respect to ice,
  !! and the dimensionless constants \f$\gamma = 4\f$, \f$A = 0.0018\f$, and \f$B = 3.6\f$.
  !!
  !! The threshold wind speed at which blowing snow occurs over flat homogeneous terrain without protruding vegetation is calculated as 
  !! Li and Pomeroy (1997):
  !!
  !! $$v_t=v_{t^*}+0.0033\left(T_a-245.88\right)^2,$$
  !!
  !! where where \f$v_{t^*} = 6.98 \f$ m s\f$^{–1}\f$ is the minimum threshold velocity. Using the criteria of Déry and Yau (1999), 
  !! if the wind speed exceeds the threshold velocity, and the temperature is below 0°C, a blowing snow event is presumed to have occurred. 
  !! These criteria ignore physical properties of the snowpack, such as snow aging, wind hardening, or hardening due to melting and refreezing, 
  !! which would decrease the amount of blowing snow, and blowing snow sublimation. The effects of redistribution due to snow transport 
  !! are also ignored. It is simply assumed that the net flux of blowing snow through the grid square boundary is zero.
  !!
  !! References:
  !! - Déry, S. J., & Yau, M. K. (1999). A climatology of adverse winter-type weather events. Journal of 
  !! Geophysical Research: Atmospheres, 104(D14), 16657‑16672. https://doi.org/10.1029/1999JD900158
  !! - Gordon, M., Simon, K., & Taylor, P. A. (2006). On snow depth predictions with the Canadian land 
  !! surface scheme including a parametrization of blowing snow sublimation. Atmosphere-Ocean, 44(3), 
  !! 239‑255. https://doi.org/10.3137/ao.440303
  !! - Li, L., & Pomeroy, J. W. (1997). Probability of occurrence of blowing snow. Journal of Geophysical 
  !! Research: Atmospheres, 102(D18), 21955‑21964. https://doi.org/10.1029/97JD01522
  !!
  !! @author M. LALANDE
  !!
  real function GD06(TA, VA10, RHOAIR, QA, PADRY)

    use classicParams, only : TFREZ

    implicit none
    !
    !     * INPUT ARRAYS.
    !
    real, intent(in) :: TA      !< Air temperature at reference height [K] (\f$T_a\f$)
    real, intent(in) :: VA10    !< Wind speed at 10 m height [m s\f$^{-1}\f$]
    real, intent(in) :: RHOAIR  !< Density of air [kg m\f$^{-3}\f$] (\f$\rho_a\f$)
    real, intent(in) :: QA      !< Specific humidity at reference height [kg kg\f$^{-1}\f$] (\f$q_a\f$)
    real, intent(in) :: PADRY   !< Partial pressure of dry air [Pa] (\f$p_{dry}\f$)
    !
    !     * TEMPORARY VARIABLES.
    !
    real :: VT      !< Threshold wind speed at which blowing snow occurs over flat homogeneous terrain 
                    !< without protruding vegetation is calculated as (Li and Pomeroy, 1997) [m s\f$^{-1}\f$]
    real :: ESAT    !< Saturated vapor pressure [Pa]
    real :: QASAT   !< Saturation specific humidity at reference height [kg kg\f$^{-1}\f$]
    !
    !-----------------------------------------------------------------------
    ! Compute the threshold wind speed at which blowing snow occurs
    VT = 6.98 + 0.0033*(TA - 245.88)**2.0

    ! If the wind speed is higher than the threshold and the air temperature is below 0°C, compute the blowing snow sublimation loss rate
    if (VA10 > VT .and. TA < TFREZ) then
      ESAT = calcEsat(TA)
      QASAT = 0.622 * ESAT / (PADRY + 0.622 * ESAT)
      ! PRINT '(A9 F12.5)', 'TZERO = ', TZERO
      ! PRINT '(A9 F12.5)', 'TA = ', TA
      ! PRINT '(A9 F12.5)', 'QA = ', QA
      ! PRINT '(A9 F12.5)', 'QASAT = ', QASAT
      ! PRINT '(A9 F12.5)', 'VA = ', VA
      ! PRINT '(A9 F12.5)', 'VT = ', VT
      ! GD06 = MAX(0.0, 0.0018 * (TZERO / TA)**4 * VT * RHOAIR * QASAT * (1 - QA / QASAT) * (VA / VT)**3.6)
      GD06 = MAX(0.0, 0.0018 * (TFREZ / TA)**4 * VT * RHOAIR * QASAT * (1 - QA / QASAT) * (VA10 / VT)**3.6)
      ! PRINT '(A9 F14.7)', 'GD06 = ', GD06
    else
      GD06 = 0.0
    end if
  end function GD06
end module snowWindSublimation
!> \file
!!
!! @author M. LALANDE
!!
!! Compute the wind sublimation losses for the snowpack on the ground according to the following parameterizations:
!! - Gordon et al. (2006): GD06 (https://doi.org/10.3137/ao.440303)
!! 
!! See the corresponding functions for details of each parameterization.
