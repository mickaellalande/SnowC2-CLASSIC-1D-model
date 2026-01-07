!> \file
!! Updates liquid and frozen water stores on canopy and in soil in
!! response to calculated sublimation, evaporation and transpiration
!! rates.
!! @author D. Verseghy, M. Lazare, G. Meyer
!
subroutine canopyWaterUpdate (EVAP, SUBL, RAICAN, SNOCAN, TCAN, THLIQ, TBAR, ZSNOW, & ! Formerly CANVAP
                              WLOST, CHCAP, QFCF, QFCL, QFN, QFC, HTCC, HTCS, HTC, &
                              FI, CMASS, TSNOW, HCPSNO, RHOSNO, FROOT, THPOR, &
                              THLMIN, DELZW, EVLOST, RLOST, IROOT, &
                              IG, ILG, IL1, IL2, JL, N, &
                              RB, RC, FRAINC, FSNOWC, LAIPAIRatio) 

  !     * APR 08/22 - G.MEYER.    Make canopy evaporation - transpiration partitioning
  !     *                         dependent on wet and dry canopy fractions and the 
  !     *                         ratio of rb to rc (EcTResistancePartSwitch; loop 300). 
  !     * SEP 15/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3.
  !     * SEP 13/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUN 20/02 - D.VERSEGHY. TIDY UP SUBROUTINE CALL; SHORTENED
  !     *                         CLASS4 COMMON BLOCK.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE SOIL
  !     *                         PERMEABLE DEPTH.
  !     * DEC 30/96 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         BUGFIXES IN CALCULATION OF QFN AND
  !     *                         QFC.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         RATIONALIZE CALCULATION OF WLOST
  !     *                         REFINE CALCULATION OF QFCL.
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         ADDITIONAL DIAGNOSTIC CALCULATIONS -
  !     *                         HTCC AND HTC.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !                                        NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CALCULATE ACTUAL EVAPORATION,
  !     *                         SUBLIMATION AND TRANSPIRATION FROM
  !     *                         VEGETATION CANOPY.
  !
  use classicParams,       only : DELT, TFREZ, HCPW, SPHW, SPHICE, SPHVEG, &
                                  RHOW, CLHMLT, CLHVAP

  use, intrinsic :: iso_fortran_env, only: r8=>real64

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IG, ILG, IL1, IL2, JL, N
  integer             :: I, J
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real(r8), intent(inout) :: TBAR(ILG,IG) !< Temperature of soil layer \f$[K] (T_g)\f$

  real, intent(inout) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
                                          !! layer \f$[m^3 m^{-3}] (\theta_l)\f$
  real, intent(inout) :: QFC   (ILG,IG)   !< Transpired water removed from soil layer \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: HTC   (ILG,IG)   !< Internal energy change of soil layer due to
                                          !! conduction and/or change in mass \f$[W m^{-2}] (I_g)\f$
  !
  real, intent(inout) :: EVAP  (ILG)  !< Evapotranspiration rate from vegetation canopy \f$[m s^{-1}]\f$
  real, intent(inout) :: SUBL  (ILG)  !< Calculated sublimation rate from vegetation canopy \f$[m s^{-1}]\f$
  real, intent(inout) :: RAICAN(ILG)  !< Intercepted liquid water stored on the canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: SNOCAN(ILG)  !< Intercepted frozen water stored on the canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: TCAN  (ILG)  !< Temperature of vegetation canopy \f$[K] (T_c)\f$
  real, intent(inout) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_g)\f$
  real, intent(inout) :: WLOST (ILG)  !< Residual amount of water that cannot be
  !! supplied by surface stores \f$[kg m^{-2}]\f$
  real, intent(inout) :: CHCAP (ILG)  !< Heat capacity of vegetation canopy \f$[J m^{-2} K^{-1}] (C_c)\f$
  real, intent(inout) :: QFCF  (ILG)  !< Sublimation from frozen water in canopy
  !! interception store \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFCL  (ILG)  !< Evaporation from liquid water in canopy
  !! interception store \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(inout) :: QFN   (ILG)  !< Sublimation from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: HTCC  (ILG)  !< Internal energy change of canopy due to changes
  !! in temperature and/or mass \f$[W m^{-2}] (I_c)\f$
  real, intent(inout) :: HTCS  (ILG)  !< Internal energy change of snow pack due to
  !! conduction and/or change in mass \f$[W m^{-2}] (I_s)\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FROOT (ILG,IG)   !< Fractional contribution of soil layer to
  !! transpiration [ ]
  real, intent(in) :: THPOR(ILG,IG)    !< Pore volume in soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !! remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable depth of soil layer \f$[m] (\Delta z_{g,w})\f$
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on
  !! modelled area \f$[ ] (X_i)\f$
  real, intent(in) :: CMASS (ILG)  !< Mass of vegetation canopy \f$[kg m^{-2}]\f$
  real, intent(in) :: TSNOW (ILG)  !< Temperature of the snow pack [C]
  real, intent(in) :: HCPSNO(ILG)  !< Heat capacity of snow pack \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(in) :: RHOSNO(ILG)  !< Density of snow pack \f$[kg m^{-3}]\f$
  real, intent(in) :: RC    (ILG) !< Stomatal resistance of vegetation over bare ground \f$[s m^{-1} ]\f$ 
  real, intent(in) :: RB    (ILG) !< Leaf boundary resistance of vegetation \f$[s m^{-1}] (r_b)\f$ 
  real, intent(in) :: FRAINC(ILG) !< Fractional coverage of canopy by liquid water over snow-free subarea [ ]
  real, intent(in) :: FSNOWC(ILG) !< Fractional coverage of canopy by frozen water over snow-free subarea [ ]
  real, intent(in) :: LAIPAIRatio(ILG)  !< LAI to PAI ratio for canopy over bare ground or snow [] 
  !
  !     * WORK ARRAYS.
  !
  real, intent(inout)    :: EVLOST(ILG), RLOST(ILG)
  !
  integer, intent(inout) :: IROOT(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: SLOST, THTRAN, THLLIM, fWet(ILG), fDry(ILG), RBRCperc(ILG) !added fWet, fDry and RBRCperc 
  ! fWet determines the contribution of the fraction of the canopy that is wet (covered by liquid or frozen water); used to
  ! determine the contributions of canopy evaporation and transpiration; RBRCperc determines the contribution of the leaf 
  ! boundary layer resistance to the leaf boundary and stomatal resistance
  logical :: EcTResistancePartSwitch
  
  EcTResistancePartSwitch = .True.  ! Switch determining whether to use the leaf boundary layer and stomatal resistances to determine 
  ! the partitioning into canopy evaporation and transpiration or the previous standard approach, where transpiration only occurs, when
  ! there is no intercepted water on the canopy; set to false for original calculation

  ! C-----------------------------------------------------------------------
  !>
  !! The calculated fluxes of liquid and frozen water from the canopy
  !! to the overlying air, obtained as outputs of waterCalcPrep.f90, are applied
  !! to the liquid and frozen intercepted water stores on the
  !! vegetation canopy, and to the liquid and frozen moisture stores
  !! at the surface and in the soil. Since there may not be sufficient
  !! water in one or more of these stores to sustain the calculated
  !! rate over the whole time step, a hierarchy of operations is
  !! followed as described below. The change of internal energy HTC in
  !! the canopy, snow and soil layers as a result of these processes
  !! is calculated as the difference in HTC between the beginning and
  !! end of the subroutine:
  !! \f[
  !! \Delta I_c = X_i \Delta (C_c T_c)/ \Delta t
  !! \f]\f[
  !! \Delta I_s = X_i \Delta (C_s T_s z_s)/ \Delta t
  !! \f]
  !! where the C terms represent volumetric heat capacities and the T
  !! terms temperatures of the canopy and snow pack, \f$\Delta t\f$ is the
  !! length of the time step, \f$z_s\f$ the snow depth, and \f$X_i\f$ the fractional
  !! coverage of the subarea under consideration relative to the
  !! modelled area. For the soil layers, since only the liquid water
  !! content is affected by these calculations, the change in internal
  !! energy of each layer is calculated from the change in liquid
  !! water content \f$\theta_l\f$ as:
  !! \f[ \Delta I_g = X_i C_w \Delta z_{g, w} \Delta (T_g \theta_l)/ \Delta t \f]
  !!


  !     * INITIALIZE ARRAYS.
  !     * (THE WORK ARRAY "IROOT" INDICATES POINTS WHERE TRANSPIRATION
  !     * CAN OCCUR.)
  !
  do I = IL1,IL2 ! loop 50
    if (FI(I) > 0.) then
      RLOST (I) = 0.0
      EVLOST(I) = 0.0
      IROOT (I) = 0
      HTCC  (I) = HTCC(I) - FI(I) * TCAN(I) * CHCAP(I) / DELT
      HTCS(I) = HTCS(I) - FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
    end if
  end do ! loop 50
  !
  do J = 1,IG ! loop 100
    do I = IL1,IL2
      if (FI(I) > 0.) then
        HTC (I,J) = HTC(I,J) - FI(I) * (TBAR(I,J) + TFREZ) * THLIQ(I,J) * &
                    HCPW * DELZW(I,J) / DELT
        if (FROOT(I,J) > 1.0E-5) IROOT(I) = 1
      end if
    end do
  end do ! loop 100
  !>
  !! Sublimation is addressed first. The predicted mass of sublimated
  !! water SLOST is calculated and compared to the frozen water in the
  !! canopy interception store, SNOCAN. If SLOST \f$\leq\f$ SNOCAN, all of the
  !! sublimated water is subtracted from SNOCAN. If not, the excess
  !! sublimation is calculated as SLOST-SNOCAN, QFCF is corrected for
  !! the canopy sublimation difference, and SNOCAN is set to zero.
  !! Next, the new value of SLOST is compared to the snowpack mass,
  !! calculated as ZSNOW*RHOSNO. If SLOST \f$\leq\f$ ZSNOW*RHOSNO, all of the
  !! remaining sublimated water is taken from the snow pack, and QFN
  !! is modified to reflect this loss. Otherwise, the excess
  !! sublimation is calculated as SLOST - ZSNOW*RHOSNO, QFN is
  !! adjusted accordingly, and ZSNOW is set to zero. There now remain
  !! no further frozen moisture stores from which sublimated water can
  !! be taken (frozen water in the soil is assumed to be immobile), so
  !! the remaining energy that had been assigned to sublimation is
  !! assigned to canopy evaporation instead, and QFCL is duly
  !! recalculated. This means, however, that a small imbalance will
  !! arise in the water budget owing to the difference between the
  !! latent heats of sublimation and evaporation. This imbalance is
  !! assigned to the housekeeping variable WLOST.
  !
  !     * SUBLIMATION CASE.  IF SNOW ON CANOPY IS INSUFFICIENT TO SUPPLY
  !     * DEMAND, RESIDUAL IS TAKEN FIRST FROM SNOW UNDERLYING CANOPY AND
  !     * THEN FROM LIQUID WATER ON CANOPY.
  !
  do I = IL1,IL2 ! loop 200
    if (FI(I) > 0. .and. SUBL(I) > 0.) then
      SLOST = SUBL(I) * DELT * RHOW
      if (SLOST <= SNOCAN(I)) then
        SNOCAN(I) = SNOCAN(I) - SLOST
        SUBL(I) = 0.0
      else
        SLOST = SLOST - SNOCAN(I)
        QFCF(I) = QFCF(I) - FI(I) * SLOST / DELT
        SNOCAN(I) = 0.0
        if (SLOST <= ZSNOW(I) * RHOSNO(I)) then
          ZSNOW(I) = ZSNOW(I) - SLOST / RHOSNO(I)
          SUBL(I) = 0.0
          QFN(I) = QFN(I) + FI(I) * SLOST / DELT
        else
          SLOST = SLOST - ZSNOW(I) * RHOSNO(I)
          QFN(I) = QFN(I) + FI(I) * ZSNOW(I) * RHOSNO(I) / DELT
          ZSNOW(I) = 0.0
          WLOST(I) = WLOST(I) - SLOST * CLHMLT / CLHVAP
          EVAP(I) = EVAP(I) + SLOST * (CLHMLT + CLHVAP) / &
                    (CLHVAP * DELT * RHOW)
          QFCL(I) = QFCL(I) + FI(I) * SLOST * (CLHMLT + CLHVAP) / &
                    (CLHVAP * DELT)
        end if
      end if
    end if
  end do ! loop 200
  !>
  !! Now canopy evaporation is addressed. It is assumed that all
  !! intercepted liquid water evaporates before transpiration begins,
  !! since there is a canopy stomatal resistance associated with
  !! transpiration and there is none associated with evaporation. The
  !! predicted mass of evaporated water RLOST is calculated and
  !! compared to the liquid water in the canopy interception store,
  !! RAICAN. If RLOST \f$\leq\f$ RAICAN, all of the evaporated water is
  !! subtracted from RAICAN. If not, the excess evaporation is
  !! calculated as RLOST - RAICAN, and QFCL is corrected for the
  !! canopy evaporation difference. This excess evaporation is now
  !! treated as transpiration. An initial check is done by referring
  !! to the diagnostic flag IROOT, which was set to 1 at the beginning
  !! of the subroutine if there was water available for transpiration
  !! in any of the soil layers. If IROOT is zero, no transpiration can
  !! occur and the excess evaporation is stored in the temporary
  !! variable EVLOST.
  !!
  !! Modifications to canopy evaporation and transpiration (in loop 300 when EcTResistancePartSwitch is set to True): \n
  !! In CLASSIC, for snow-covered canopies the evaporative flux from the canopy (ET\f$_{c}\f$) is first assigned to sublimation, as liquid water, if present, 
  !! is assumed to be within or underneath the snow and no $T$ is expected to occur. When there is only liquid water present on the canopy, we modified the original 
  !! formulation so that \f$T\f$ is allowed to occur from a partially-wet canopy instead of from a completely dry canopy only. To separate the calculated amount of 
  !! evapotranspired water (ET\f$_{c}\f$) into \f$E_c\f$ and \f$T\f$, we determine the wet (\f$f_{wet}\f$) and dry (\f$f_{dry}\f$) fractions of the canopy 
  !! similar to Fan et al. (2019) \cite Fan2019-ap as 
  !! \f[f_{wet} = \left\{\begin{array}{l l} F_l \qquad \text {for $F_l$ $\geq$ 0.01 and $F_l$ $\leq$ 0.99} \\
  !! 0 \qquad \text {for $F_l$ $<$ 0.01} \\
  !! 1 \qquad \text {for $F_l$ $>$ 0.99} \end{array} \right.
  !! \f] 
  !! \f[f_{dry} = \left\{\begin{array}{l l} (1 - f_{wet}) \frac{\text{LAI}}{\text{PAI}}  \qquad \text {for $f_{wet}$ $\geq$ 0.01 and $f_{wet}$ $\leq$ 0.99} \\
  !! 1 \qquad \text {for $f_{wet} <$ 0.01} \\
  !! 0 \qquad \text {for $f_{wet} >$ 0.99} \end{array} \right.
  !! \f]
  !! In general, only the leaves and not stems of a canopy can transpire, \f$f_{dry}\f$ is adjusted by the LAI to PAI ratio. \f$F_l\f$ is the fractional coverage of 
  !! the canopy covered by liquid water (unitless) determined as
  !! \f[F_l = \left\{\begin{array}{l l} \text{min($W_l$ / $W_{l,max}$, 1)} \qquad \text {for $W_{l,max} >$ 0} \\
  !! 0 \qquad \text {for $W_{l,max}$ = 0} \end{array} \right. 
  !! \f]
  !! where \f$W_l\f$ (\f$kg ~m^{-2}\f$) is the amount of liquid water stored on the canopy and \f$W_{l,max}\f$ (\f$kg ~m^{-2}\f$) is the storage capacity of the canopy 
  !! for liquid water, which is calculated as 
  !! \f$W_{l,max} = p_l \times \text{PAI}\f$
  !! with the maximum storage of liquid water \f$p_l\f$ set as 0.20 \f$kg ~m^{-2}\f$ Bartlett et al. (2006) \cite Bartlett2006-xp .
  !! \f$W_l\f$ is calculated as the sum of \f$W_l\f$ of the previous time step and the rainfall intercepted by the canopy during the current time step
  !! \f$W_{l,t} = \min(W_{l,t-1} + \Delta t \rho_w (P - \chi P), W_{l,max}),\f$
  !! where \f$P\f$ is the rainfall rate (\f$m s^{-1}\f$), \f$\chi\f$ is the canopy gap fraction (unitless), \f$\Delta t\f$ is the model physics timestep (s) and 
  !! \f$\rho_w\f$ (\f$kg ~m^{-3}\f$) the density of liquid water. To determine the canopy fractional coverage of liquid water exposed to the air, \f$F_l\f$ is decreased 
  !! by the fractional snow coverage (\f$F_s\f$).
  !! \f$F_l = \text{max}(0,\text{min}(F_l - F_s, 1))\f$
  !! and, similar to \f$F_l\f$, \f$F_s\f$ is found by
  !! \f[F_s = \left\{\begin{array}{l l} \text{min($W_f$ / $W_{f,max}$,1)} \qquad \text {for $W_{f,max} >$ 0}) \\
  !! 0 \qquad \text {for $W_{f,max}$ = 0}), \end{array} \right.
  !! \f]
  !! where \f$W_f\f$ (\f$kg ~m^{-2}\f$) is the amount of frozen water stored on the canopy and \f$W_{f,max}\f$ (\f$kg ~m^{-2}\f$) is the storage capacity of the canopy 
  !! for frozen water. If there is no plant  available water in the root zone, the wet canopy fraction is set to 1, as \f$T\f$ is not allowed to occur. 
  !!
  !! The predicted mass of water evapotranspired from the canopy (\f$W_E; ~kg ~m^{-2}\f$; RLOST in the code) is calculated as 
  !! \f$W_E = \text{ET}_c \times \rho_w \Delta t,\f$
  !! where ET\f$_{c}\f$ is the evapotranspiration rate from the canopy (\f$m ~s^{-1}\f$). The wet and dry canopy fractions as well as \f$F_{RbRc}\f$ determine 
  !! the fractions of \f$W_E\f$ coming from \f$E_c\f$ and \f$T\f$, respectively. The amount of \f$W_l\f$ is adjusted by the amount of water evaporated from the wet 
  !! canopy fraction as
  !! \f$W_l = \text{$W_l$ - (1 - $f_{dry}$) (1 - $F_{RbRc}$) $W_E$}~\text{for $W_E$ (1 - $f_{dry}$)(1 - $F_{RbRc}$) $\le$ $W_l$}\f$
  !! and \f$W_E\f$ is reduced by the amount being evaporated
  !! \f[W_E = \left\{\begin{array}{l l} \text{$W_E$ ($F_{RbRc}$ + $f_{dry}$ - $f_{dry}$ $F_{RbRc}$)} \qquad \text{for $W_E$ (1 - $f_{dry}$) (1 - $F_{RbRc}$) $\le$ $W_l$} \\
  !! \text{$W_E$ - $W_l$} \qquad \text{for $W_E$ (1 - $f_{dry}$) (1 - $F_{RbRc}$) $> W_l$} \end{array} \right.
  !! \f]
  !! The contribution of the leaf boundary layer resistance (\f$r_b; s ~m^{-1}\f$) to the total resistance, the sum of \f$r_b\f$ and the stomatal resistance (\f$r_c\f$ 
  !! or \f$1/g_c; s ~m^{-1}\f$), is calculated as a proportion of the total resistance from the leaf boundary layer and stomata to determine when canopy evaporation is 
  !! dominant and when \f$T\f$ can occur, as
  !! \f$F_{RbRc} = r_b / (r_b + r_c).\f$
  !! In the second case, where \f$W_l\f$ could not meet the calculated amount of water to be evaporated, \f$W_l\f$ is then set to zero. If the 
  !! predicted mass of water evapotranspired from the vegetation (\f$W_E\f$) after considering evaporation from wet leaves is greater than zero, it is treated as 
  !! \f$T\f$. If there is enough water available in the root zone and \f$T\f$ can occur, the soil water content removed by \f$T\f$ and the \f$T\f$ flux are calculated 
  !! for each soil layer and the liquid water content of each soil layer containing roots is updated. 
  !
  !     * EVAPORATION.  IF WATER ON CANOPY IS INSUFFICIENT TO SUPPLY
  !     * DEMAND, ASSIGN RESIDUAL TO TRANSPIRATION.
  !
  if (EcTResistancePartSwitch) then
    do I = IL1,IL2 ! loop 300 (version using Ec and T partitioning dependent on resistances; G. Meyer, July2021)  
      RBRCperc(I) = RB(I) / (RB(I) + RC(I)) 
      if (RBRCperc(I) < 0.01 .or. IROOT(I) == 0) RBRCperc(I) = 0. 
      fWet(I) = FRAINC(I)
      fDry(I) = (1 - fWet(I)) * LAIPAIRatio(I)
      if (fWet(I) < 0.01) then
      fWet(I) = 0.  
      fDry(I) = 1.
      end if  
      if (fWet(I) > 0.99 .or. IROOT(I) == 0) then
        fWet(I) = 1.   
        fDry(I) = 0.
      end if
      if (FI(I) > 0. .and. EVAP(I) > 0.) then
        RLOST(I) = EVAP(I) * RHOW * DELT
        if (RLOST(I) * (1 - fDry(I)) * (1.0 - RBRCperc(I)) <= RAICAN(I)) then
          RAICAN(I) = RAICAN(I) - RLOST(I) * (1 - fDry(I)) * (1.0 - RBRCperc(I))
          EVAP  (I) = 0.
          RLOST (I) = RLOST(I) * (RBRCperc(I) + fDry(I) - fDry(I) * RBRCperc(I))
          if (RLOST(I) > 0.) QFCL(I) = QFCL(I) - FI(I) * RLOST(I) / DELT
        else
          RLOST(I) = RLOST(I) - RAICAN(I)
          QFCL(I) = QFCL(I) - FI(I) * RLOST(I) / DELT
          if (IROOT(I) == 0) EVLOST(I) = RLOST(I)
          EVAP  (I) = 0.
          RAICAN(I) = 0.
          end if
        end if
    end do ! loop 300
  else
    do I = IL1,IL2 ! loop 300 (original way)
      if (FI(I) > 0. .and. EVAP(I) > 0.) then
        RLOST(I) = EVAP(I) * RHOW * DELT
        if (RLOST(I) <= RAICAN(I)) then
          RAICAN(I) = RAICAN(I) - RLOST(I)
          EVAP  (I) = 0.
          RLOST (I) = 0.
        else
          RLOST(I) = RLOST(I) - RAICAN(I)
          QFCL(I) = QFCL(I) - FI(I) * RLOST(I) / DELT
          if (IROOT(I) == 0) EVLOST(I) = RLOST(I)
          EVAP  (I) = 0.
          RAICAN(I) = 0.
        end if
      end if
    end do ! loop 300
  end if
  !>
  !! The next loop (loop 400) is performed if IROOT=1, i.e. if transpiration is
  !! possible. For each soil layer, the volumetric water content that
  !! is removed by transpiration, THTRAN, is calculated from RLOST
  !! (converted to a volumetric content by dividing by the density of
  !! water and the permeable thickness of the soil layer), and the
  !! fractional contribution of the soil layer FROOT. If there is
  !! enough liquid water in the soil layer to supply THTRAN, the
  !! diagnostic transpiration flux QFC for the layer is updated using
  !! THTRAN, and the liquid water content of the layer is updated as
  !! THLIQ-THTRAN. If not, QFC is updated using the available water in
  !! the soil layer, THLIQ is set to THLMIN, and the residual
  !! untranspired water is added to EVLOST.
  !
  !     * TRANSPIRATION.
  !
  do J = 1,IG ! loop 400
    do I = IL1,IL2
      if (FI(I) > 0. .and. IROOT(I) > 0) then
        if (DELZW(I,J) > 0.0) then
          THTRAN = RLOST(I) * FROOT(I,J) / (RHOW * DELZW(I,J))
        else
          THTRAN = 0.0
        end if
        if (THPOR(I,J) < THLMIN(I,J)) then
          THLLIM = THPOR(I,J)
        else
          THLLIM = THLMIN(I,J)
        end if
        if (THTRAN <= (THLIQ(I,J) - THLLIM)) then
          QFC  (I,J) = QFC(I,J) + FI(I) * RLOST(I) * FROOT(I,J) / DELT
          THLIQ(I,J) = THLIQ(I,J) - THTRAN
        else
          QFC  (I,J) = QFC(I,J) + FI(I) * (THLIQ(I,J) - THLLIM) * RHOW * &
                       DELZW(I,J) / DELT
          EVLOST (I) = EVLOST(I) + (THTRAN + THLLIM - THLIQ(I,J)) * RHOW * &
                       DELZW(I,J)
          THLIQ(I,J) = THLLIM
        end if
      end if
    end do
  end do ! loop 400
  !>
  !! In the final cleanup (loop 500), the canopy heat capacity is recalculated,
  !! the contents of EVLOST are added to WLOST, and the remaining
  !! internal energy calculations are completed.
  !
  !     * CLEANUP.
  !
  do I = IL1,IL2 ! loop 500
    if (FI(I) > 0.) then
      CHCAP(I) = RAICAN(I) * SPHW + SNOCAN(I) * SPHICE + CMASS(I) * SPHVEG
      WLOST(I) = WLOST(I) + EVLOST(I)
      HTCC  (I) = HTCC(I) + FI(I) * TCAN(I) * CHCAP(I) / DELT
      HTCS(I) = HTCS(I) + FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
    end if
  end do ! loop 500
  !
  do J = 1,IG ! loop 550
    do I = IL1,IL2
      if (FI(I) > 0.) then
        HTC (I,J) = HTC(I,J) + FI(I) * (TBAR(I,J) + TFREZ) * THLIQ(I,J) * &
                    HCPW * DELZW(I,J) / DELT
      end if
    end do
  end do ! loop 550
  !
  return
end subroutine canopyWaterUpdate
