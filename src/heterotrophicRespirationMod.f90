!> \file
!> Central module for all heterotrophic respiration-related operations
module heterotrophicRespirationMod

  implicit none

  ! Subroutines contained in this module:
  public  :: heterotrophicRespiration
  public  :: soilCResPeat
  public  :: updatePoolsHetResp

contains

  !> \ingroup heterotrophicrespirationmod_heterotrophicRespiration
  !! @{
  !> Heterotrophic respiration subroutine
  !> Calculates litter and soil carbon respiration for mineral and 
  !! peatland soils. Vegetated areas are treated separately than 
  !! bare ground. As well the LUC product pools are also considered.
  !> @author Vivek Arora, Joe Melton, Yuanqiao Wu
  subroutine heterotrophicRespiration (il1, il2, ilg, ipeatland, fcancmx, & ! In
                                       litrmass, soilcmas, delzw, thpor, tbar, & ! In
                                       psisat, thliq, sort, bi, isand, thice, & ! In
                                       fg, litrmsmoss, peatdep, wtable, zbotw, & ! In
                                       useTracer, tracerLitrMass, tracerSoilCMass, tracerMossLitrMass, & ! In
                                       ltresveg, scresveg, litresmoss, socres_peat, & ! Out
                                       resoxic, resanoxic, & ! Out
                                       ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer) ! Out
    ! ------
    use classicParams,         only : icc, ignd, zero, tanhq10, &
                                      bsratelt_g, bsratesc_g, r_depthredu, &
                                      tcrit, frozered, kk, bsratelt, bsratesc, &
                                      tfrez, bsrateltms, iccp2, iccp1, bsratesc_peat

    implicit none

    integer, intent(in) :: ilg   !< number of gridcells/tiles in a latitude circle
    integer, intent(in) :: il1   !< il1=1
    integer, intent(in) :: il2   !< il2=ilg
    integer, intent(in) :: isand(:,:) !< flag for soil/bedrock/ice/glacier
    integer, intent(in) ::  ipeatland(:)  !< Peatland flag. 0 non-peatlands, 1 = bog, 2 = fen.
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
                                  !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the 
                                  !! tracer values in the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                  !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                  !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    real, intent(in) :: fcancmx(:,:)      !< max. fractional coverage of CTEM's pfts, but this can be
                                          !< modified by land-use change, and competition between pfts
    real, intent(in) :: litrmass(:,:,:)   !< litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(in) :: soilcmas(:,:,:)   !< soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(in) :: delzw(:,:)         !< thicknesses of the permeable soil layers (m)
    real, intent(in) :: thpor(:,:)         !< Soil total porosity per layer  \f$(cm^3 cm^{-3})\f$ -
    real, intent(in) :: tbar(:,:)          !< Soil temperature, K
    real, intent(in) :: psisat(:,:)        !< Saturated soil matric potential (m)
    real, intent(in) :: thliq(:,:)         !< liquid moisture content of soil layers
    
    integer, intent(in) :: sort(:)         !< index for correspondence between biogeochemical (CTEM) pfts 
                                           !! and size of parameters vectors
    real, intent(in) :: bi(:,:)            !< Brooks and Corey/Clapp and Hornberger b term
    real, intent(in) :: thice(:,:)         !< frozen moisture content of soil layers
    real, intent(in) :: fg(:)              !< fraction of bare ground []  
    real, intent(in) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$    
    real, intent(in) :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(in) :: wtable(:)          !< water table (m)
    real, intent(in) :: zbotw(:,:)         !< bottom of soil layers
    real, intent(in) :: peatdep(:)         !< peat depth (m)
    real, intent(in) :: litrmsmoss(:)   !< moss litter C pool (kgC/m2)
    
    real, intent(out) :: ltresveg(ilg,iccp2,ignd)     !< litter respiration for the given vegetated sub-area [ \f$u-mol co2/m2.sec\f$ ]
    real, intent(out) :: scresveg(ilg,iccp2,ignd)     !< soil carbon respiration for the given vegetated sub-area [ \f$u-mol co2/m2.sec\f$ ]
    real, intent(out) :: litresmoss(ilg)    !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: socres_peat(ilg)   !< heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: resoxic(ilg)       !< oxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: resanoxic(ilg)     !< anoxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: ltResTracer(ilg,iccp2,ignd)  !< Tracer litter respiration over the given unvegetated sub-area in umol co2/m2.s TODO: Units 
    real, intent(out) :: sCResTracer(ilg,iccp2,ignd)  !< Tracer soil C respiration over the given unvegetated sub-area in umol co2/m2.s TODO: Units
    real, intent(out) :: litResMossTracer(ilg)    !< Tracer moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) TODO: Units
    real, intent(out) :: soCResPeatTracer(ilg)    !< Tracer heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) TODO: Units

    ! Local
    real :: litrq10               !<
    real :: soilcq10              !<
    real :: q10funcLitr(ilg,ignd) !<
    real :: q10funcSoilC(ilg,ignd)!<
    real :: socmoscl(ilg,ignd)    !< soil moisture scalar for soil carbon decomposition
    real :: ltrmoscl(ilg,ignd)    !< soil moisture scalar for litter decomposition
    real :: psi(ilg,ignd)         !< soil matric potential (m)
    real :: reduceatdepth(ilg,ignd) !<
    integer :: i, j, k

    ! -------------------------------------------------------------------------

    ! initialize required arrays to zero

    socmoscl = 0.0   ! soil moisture scalar for soil carbon decomposition
    ltrmoscl = 0.0   ! soil moisture scalar for litter decomposition
    ltresveg = 0.0   ! litter resp. rate for each pft
    scresveg = 0.0   ! soil c resp. rate for each pft
    ltResTracer = 0.
    sCResTracer = 0.
    litresmoss = 0.0
    socres_peat = 0.0
    litResMossTracer = 0.
    soCResPeatTracer = 0.

    ! initialization ends
    do i = il1,il2
      ! Quick warning check to make sure only fens and bogs are considered as presently 
      ! formulated. If other systems are considered (such as moss on uplands), this scheme 
      ! will need to be adapted so stop the user now.
      if (ipeatland(i) > 2) then 
        print*,'Error: heterotrophicRespiration - not yet setup for ipeatland > 2'
        call errorHandler('heterotrophicRespiration', - 1)
      end if 
    end do

    !> Find moisture scalar for litter and soil C decomposition
    !! this is modelled as function of logarithm of matric potential.

    ! First find the matric potential:
    do j = 1,ignd
      do i = il1,il2
        if (isand(i,j) /= -3 .and. isand(i,j) /= -4) then ! for mineral or peat soils
          if (ipeatland(i) == 0) then ! mineral soils

            ! We don't place a lower limit on psi as it is only used here and the
            ! value of psi <= psisat is just used to select a scomotrm value.
            
            if (thice(i,j) <= thpor(i,j)) then ! there is some remaining pore space 
              psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-bi(i,j))
            else ! if the whole pore space is ice then suction is assumed to be very high.                
              psi(i,j) = 10000.0
            end if
            
          else ! peatlands
            
            if ((thliq(i,j) + thice(i,j) + 0.01 < thpor(i,j) .and. &
                tbar(i,j) < 273.16) & !freezing/frozen soils with some empty pore space 
                .or. thice(i,j) > thpor(i,j)) then  ! or  saturated frozen soils 
              psi(i,j) = 0.001 ! set to saturation  !FLAG: This is opposite mineral soils, intended? JM
            else
              psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-bi(i,j))
            end if
          end if 
          
          ! Now find the moisture scalars.       
          !> The difference between moisture scalar for litter and soil C
          !! is that the litter decomposition is not constrained by high
          !! soil moisture (assuming that litter is always exposed to air).
          !! In addition, we use moisture content of the top soil layer
          !! as a surrogate for litter moisture content. 

          if (psi(i,j) >= 10000.0) then ! very dry 
            ltrmoscl(i,j) = 0.2            
            socmoscl(i,j) = 0.2
          else if (psi(i,j) < 10000.0 .and.  psi(i,j) > 6.0) then ! dry to moist 
            ltrmoscl(i,j) = 1.0 - 0.8 * ( (log10(psi(i,j)) - log10(6.0)) &
                                        / (log10(10000.0) - log10(6.0)) )
            socmoscl(i,j) = ltrmoscl(i,j)
          else if (psi(i,j) <= 6.0 .and. psi(i,j) >= 4.0) then ! optimal 
            ltrmoscl(i,j) = 1.0
            socmoscl(i,j) = 1.0
          else if (psi(i,j) < 4.0 .and. psi(i,j) > psisat(i,j)) then ! too wet (for soil C)
            socmoscl(i,j) = 1.0 - 0.5 * ((log10(4.0) - log10(psi(i,j))) &
                                       / (log10(4.0) - log10(psisat(i,j))) )
            ! Now we distinguish the response of mineral soils from peat soils at high 
            ! soil moistures:
            if (ipeatland(i) == 0) then ! mineral soils 
              if (j ==1) then  ! only the litter at the surface (layer 1) are not impeded by too moist conditions.
                ltrmoscl(i,j) = 1.0
              else ! otherwise assume same environment as soil C (hence impeded)
                ltrmoscl(i,j) = socmoscl(i,j)
              end if 
            else ! peatlands ! FLAG Keeping depth same as surface for now. JM
              ltrmoscl(i,j) = 1.0 - 0.99*((log10(4.0) - log10(psi(i,j))) &
                                       / (log10(4.0) - log10(psisat(i,j))))
            end if 
          else if (psi(i,j) <= psisat(i,j)) then ! saturated conditions (impeded respiration)
            socmoscl(i,j) = 0.5
            ! Again we have a different response for litter in peatlands vs. mineral soils:
            if (ipeatland(i) == 0) then !mineral soils 
              if (j ==1) then  ! only the litter at the surface is not impeded at all. 
                ltrmoscl(i,j) = 1.0
              else ! otherwise assume same environment as soil C
                ltrmoscl(i,j) = socmoscl(i,j)
              end if
            else ! peatlands 
              ! Assume the litter respiration is almost shutdown.
              ltrmoscl(i,j) = 0.01  !FLAG Keeping depth same as surface for now.  JM
            end if  
          end if
        else  ! bedrock or ice 
          socmoscl(i,j) = 0.2
          ltrmoscl(i,j) = 0.2
        end if

        socmoscl(i,j) = max(0.2,min(socmoscl(i,j),1.0))
        ltrmoscl(i,j) = max(0.2,min(ltrmoscl(i,j),1.0))

        !! Next find the q10 response function to scale base respiration
        !! rate from 15 C to current temperature, we do litter first
        litrq10 = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - (tbar(i,j) - tfrez))))
        
        ! Soil has the same values so copy the litter Q10 values.
        soilcq10 = litrq10
        
        !! Apply a step reduction in q10func when the soil layer freezes. This is to reflect the
        !! lost mobility of the microbial population due to less liquid water
        if (tbar(i,j) - tfrez > tcrit) then ! unfrozen soils 
          q10funcLitr(i,j)  =  litrq10**(0.1 * (tbar(i,j) - tfrez - 15.0))
          q10funcSoilC(i,j) = soilcq10**(0.1 * (tbar(i,j) - tfrez - 15.0))
        else !frozen/freezing soils 
          q10funcLitr(i,j)  =  litrq10**(0.1 * (tbar(i,j) - tfrez - 15.0)) * frozered
          q10funcSoilC(i,j) = soilcq10**(0.1 * (tbar(i,j) - tfrez - 15.0)) * frozered
        end if

        !> Reduce the respiration at depth due to unresolved depth dependent processes including
        !! soil microbial population dynamics, pore-scale oxygen availability, mineral sorption,
        !! priming effects, and other unrepresented processes. This is following Lawrence et al.
        !! Enviro Res Lett 2015 \cite Lawrence2015-tj. We apply this for all soils.
        reduceatdepth(i,j) = 1.0 ! exp(-zbotw(i,j) / r_depthredu) FLAG turned off for now.
        
      end do !il1, il2
    end do !ignd

  !> Use temperature of the litter and soil C pools,and their soil
  !! moisture scalars to find respiration rates from these pools
  !! Here we consider the mineral soil, peatland (peat), and LUC products pool 
  !! soil carbon and litter respiration fluxes.
  !!
  !! Vegetated areas are calculated the same for peatland and mineral soil tiles.
  !! Bareground respiration only applies to mineral soils. Peatlands have peat and thus for them 
  !! we calculate peat and moss respiration separately. LUC product pools are calculated the same as 
  !! mineral soil bareground. LUC product pools cannot exist on peatland tiles. This would imply that 
  !! a peatland has been drained and converted to croplands. Presently we can't account for the impacts 
  !! of peat drainage so all LUC in a grid cell must occur on the upland tiles.

    do i = il1,il2
      ! Bare ground (non-peatland) and the LUC products pool 
      do j = iccp1, iccp2 
        do k = 1,ignd
          ! Bare ground calculated only if there is some and it is not a 
          ! peatland tile. LUC product pools are only considered looking at mineral tiles
          ! and don't need a fractional area.
          if (ipeatland(i) == 0 .and. (fg(i) > zero .or. j == iccp2)) then 
            ! 2.64 converts bsratelt_g from kg c/kg c.year to u-mol co2/kg c.s
            ltresveg(i,j,k) = ltrmoscl(i,k) * litrmass(i,j,k) * bsratelt_g  &
                              * 2.64 * q10funcLitr(i,k) * reduceatdepth(i,k)
          
            if (useTracer > 0) ltResTracer(i,j,k) = ltrmoscl(i,k) * tracerLitrMass(i,j,k) &
                                                  * bsratelt_g * 2.64 * q10funcLitr(i,k) * reduceatdepth(i,k)
            
            scresveg(i,j,k) = socmoscl(i,k) * soilcmas(i,j,k) * bsratesc_g &
                              * 2.64 * q10funcSoilC(i,k) * reduceatdepth(i,k)

            if (useTracer > 0) sCResTracer(i,j,k) = socmoscl(i,k) * tracerSoilCMass(i,j,k) &
                                                  * bsratesc_g * 2.64 * q10funcSoilC(i,k) * reduceatdepth(i,k)
          end if
        end do       
      end do 
      
      !> Peatlands do not have depth dependent soil and litter carbon pools. Instead the 
      !! peatland litter from mosses is assumed to exist only in the first layer (mosses 
      !! are 3 - 4 cm thick) and their humification leads to peat soil carbon. The peat 
      !! soil carbon pool is of known depth but is not tracked per layer like mineral 
      !! soil carbon. If the peatland tiles is vegetated, 
      if (ipeatland(i) > 0) then 
        
       !> Moss litter respiration for tiles with moss cover. Only the first layer is 
       !! used in this calculation. 
       !! Convert from kg c/kg c.year to u-mol co2/kg c.s using 2.64
        litresmoss(i) =  ltrmoscl(i,1) * litrmsmoss(i) * bsrateltms * 2.64 * q10funcLitr(i,1)
        
        ! FLAG litResMossTracer not yet dealt with!
      
        !! When peatlands are simulated find the peat heterotrophic respiration.
        ! Called regardless to set outputs to 0.
        call soilCResPeat(il1, il2, ilg, ipeatland, & ! In
                        isand, peatdep, wtable, tbar,  & ! In
                        zbotw, delzw, useTracer, & ! In
                        socres_peat, resoxic, resanoxic, & ! Out
                        soCResPeatTracer) ! Out
      
        ! FLAG soCResPeatTracer not dealt with!

      end if 
      
      ! Vegetated areas on both mineral and peatland tiles. 
      do j = 1,icc
        do k= 1,ignd
          if (fcancmx(i,j) > 0.) then
            ! 2.64 converts bsratelt from kg c/kg c.year to u-mol co2/kg c.s
            ltresveg(i,j,k) = ltrmoscl(i,k) * litrmass(i,j,k) * bsratelt(sort(j)) * 2.64 &
                                            * q10funcLitr(i,k) * reduceatdepth(i,k)
            if (useTracer > 0) ltResTracer(i,j,k) = ltrmoscl(i,k) * tracerLitrMass(i,j,k)&
                                                    * bsratelt(sort(j)) * 2.64 * q10funcLitr(i,k) &
                                                    * reduceatdepth(i,k)
            if (ipeatland(i) == 0) then !uplands
              scresveg(i,j,k) = socmoscl(i,k) * soilcmas(i,j,k) * bsratesc(sort(j)) * 2.64 &
                              * q10funcSoilC(i,k) * reduceatdepth(i,k)
              if (useTracer > 0) sCResTracer(i,j,k) = socmoscl(i,k) * tracerSoilCMass(i,j,k) &
                                                  * bsratesc(sort(j)) * 2.64 * q10funcSoilC(i,k) * reduceatdepth(i,k)
            else ! peatlands
              scresveg(i,j,k) = socmoscl(i,k) * soilcmas(i,j,k) * bsratesc_peat(sort(j)) * 2.64 &
                              * q10funcSoilC(i,k) * reduceatdepth(i,k)
              if (useTracer > 0) sCResTracer(i,j,k) = socmoscl(i,k) * tracerSoilCMass(i,j,k) &
                                                  * bsratesc_peat(sort(j)) * 2.64 * q10funcSoilC(i,k) * reduceatdepth(i,k)
            end if
          end if
        end do
      end do 
    end do 

  return

  end subroutine heterotrophicRespiration
  !! @}
  ! ! ---------------------------------------------------------------
  !> \ingroup heterotrophicrespirationmod_soilCResPeat
  !! @{
  !> Grid average peat soil heterotrophic respiration subroutine (equations are in module-level description)
  !> @author Yuanqiao Wu, J. Melton
  subroutine  soilCResPeat (il1, il2, ilg, ipeatland, & ! In
                          isand, peatdep, wtable, tbar, & ! In
                          zbotw, delzw, useTracer, & ! In
                          socres_peat, resoxic, resanoxic, & ! Out
                          soCResPeatTracer) ! Out

    use classicParams, only :icc, ignd, tanhq10, dctmin, dcbaset, TFREZ
    use peatlandsMod,  only : peatStorage
    
    
    implicit none

    integer, intent(in) :: il1, il2, ilg
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
                                    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
                                    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    integer, dimension(ilg,ignd), intent(in) :: isand   !<
    integer, dimension(ilg), intent(in) :: ipeatland    !< peatland flag, 0 = not peatland, 1 = bog, 2 = fen
    real, dimension(ilg), intent(in) :: peatdep         !< peat depth (m)
    real, dimension(ilg,ignd), intent(in) :: zbotw      !< bottom of soil layers (m)
    real, dimension(ilg,ignd), intent(in) :: delzw      !< thicknesses of the soil layers (m)
    real, dimension(ilg), intent(in) :: wtable          !< water table (m)
    real, dimension(ilg,ignd), intent(in) :: tbar       !< soil temperature (K)
    real, dimension(ilg), intent(inout) :: socres_peat         !< soil C respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: resoxic         !< respiration rate of the oxic compartment (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: resanoxic       !< respiration rate of the oxic compartment (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: soCResPeatTracer(ilg)           !< Tracer heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    !     internal variables
    integer :: i,j
    real :: Cso (ilg)            !< carbon mass in the oxic compartment (\f$kg C m^{-2}\f$)
    real :: Csa (ilg)            !< carbon mass in the anxic compartment (\f$kg C m^{-2}\f$)
    ! real :: tracerCso (ilg)      !< Tracer carbon mass in the oxic compartment (\f$kg C m^{-2}\f$)
    ! real :: tracerCsa (ilg)      !< Tracer carbon mass in the anxic compartment (\f$kg C m^{-2}\f$)
    real :: fto(ilg)             !< temperature factor of the oxic soil respiration
    real :: fta(ilg)             !< temperature factor of the anoxic soil respiration
    real :: tsoilo(ilg)          !< average temperature in the oxic compartment(C)
    real :: tsoila(ilg)          !< average temperature of the anoxic compartment (C)
    real :: ewtable(ilg)         !< effective water table depth (m)
    real :: ratescpo(ilg)        !< oxic respiration rate constant(\f$\mu mol CO_2 [kg C]^{-2} s^{-1}\f$)
    real :: ratescpa(ilg)        !< anoxic respiration rate constant(\f$\mu mol CO_2 [kg C]^{-2} s^{-1}\f$)
    real :: soilq10o(ilg)        !< q10 coefficient of oxic soil
    real :: soilq10a(ilg)        !< q10 coefficient of anoxic soil
    integer :: lewtable(ilg)     !< layer index of the water table layer
    integer :: indexLastPermLayer !index of the last permable layer in the ground column 
    !    ------------------------------------------------------------------
    
    ! Initialization
    resoxic(:) = 0.0
    resanoxic(:) = 0.0
    tsoilo(:) = 0.0
    tsoila(:) = 0.0
    soCResPeatTracer = 0.

    do i = il1,il2
      if (ipeatland(i) == 1 .or. ipeatland(i) == 2) then 

        !> Calculate soil respiration from peat

        !> Find the effective water table depth 
        ewtable(i) = wtable(i)
        
        ! And the layer index of the water table to divide the peat soil into two compartments
        if (ewtable(i) <= 0.0) then ! flooded
          lewtable(i) = 0
        else if (ewtable(i) > 0.0 .and. ewtable(i) <= 0.1) then ! below ground surface but in the moss layer
          lewtable(i) = 1
        else if (ewtable(i) > 999.) then 
          write(6,*),'Uninitialized water table detected'
          call errorHandler('soilCResPeat', - 8)
        else ! deeper in the soil column.
          do j = 1,ignd
            if (ewtable(i) < zbotw(i,j)) exit
            lewtable(i) = j + 1
          end do
        end if

        ! Find in the index of the last permeable layer in the ground column
        do j = 2,ignd 
          indexLastPermLayer = j - 1
          if (delzw(i,j) < 0.0001) exit           
        end do

        !> Find the temperature in litter, oxic soil and anoxic soil in kelvin
        !! lewtable is the layer index of the water table layer, lewtable = 0
        !! indicates WTD is above the ground surface.
        !! Set the oxic layer temperature to dctmin (minimum soil respiration
        !! temperature) when the entire soil is in the anoxic zone.

        if (lewtable(i) == 0) then ! WT is at or above the surface

          do j = 1,indexLastPermLayer ! so find the temp of the total soil column.
            ! FLAG JM - this might not be appropriate for runs with the newer deeper soils and 20 layers. There 
            ! likely should be weighting by where in the column the temperatures are coming from, rather than thickness. So deeper have 
            ! less impact since respiration is reduced at depth anyway. Could use the reduction with depth from the 
            ! multilayer soil C scheme. March 24 2020.
            tsoila(i) = tsoila(i) + tbar(i,j) * delzw(i,j) 
          end do
          tsoila(i) = tsoila(i) / zbotw(i,indexLastPermLayer)
          tsoilo(i) = dctmin 

        else if (lewtable(i) == 1) then ! WT is in the moss layer

          tsoilo(i) = tbar(i,1) ! oxic temperature is just the temperature of the moss layer.
          do j = 2, indexLastPermLayer
            tsoila(i) = tsoila(i) + tbar(i,j) * delzw(i,j)
          end do
          ! Assume that, while the true water table depth is somewhere within the moss 
          ! layer, we can ignore that small contribution to the tsoila calculation.
          tsoila(i) = tsoila(i) / (zbotw(i,indexLastPermLayer) - delzw(i,1))
        
        else ! WT is below layer 1

          do j = 1, lewtable(i) - 1
            tsoilo(i) = tsoilo(i) + tbar(i,j) * delzw(i,j)
          end do
          ! add in the partial layer contribution.
          tsoilo(i) = tsoilo(i) + tbar(i,lewtable(i)) * (ewtable(i) - zbotw(i,lewtable(i) - 1))
          tsoilo(i) = tsoilo(i) / ewtable(i)
          
          if (zbotw(i,indexLastPermLayer) - ewtable(i) > 0.) then ! water table is within soil permeable depth 
            do j = indexLastPermLayer,lewtable(i) + 1, - 1
              tsoila(i) = tsoila(i) + tbar(i,j) * delzw(i,j)
            end do
            tsoila(i) = tsoila(i) + tbar(i,lewtable(i)) * (zbotw(i,lewtable(i)) -         (i)) 
            tsoila(i) = tsoila(i) / (zbotw(i,indexLastPermLayer) - ewtable(i))
          end if 

        end if

        !> Calculate the temperature multiplier (ftsocres) for oxic and anoxic
        !! soil compartments

        soilq10o(i) = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - (tsoilo(i) - tfrez))))
        soilq10a(i) = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - (tsoila(i) - tfrez))))

        fto(i) = soilq10o(i) ** (0.1 * (tsoilo(i) - tfrez - 15.0))
        fta(i) = soilq10a(i) ** (0.1 * (tsoila(i) - tfrez - 15.0))

        !> Find the heterotrophic respiration rate constant in the oxic and
        !! anoxic (unit in yr-1), based on Fig.2b in Frolking 2001

        if (ipeatland(i) == 1) then ! bogs
          if (ewtable(i) < 0.0) then ! flooded

            ratescpo(i) = 0.0
            ratescpa(i) = -0.183 * exp(-18.0 * peatdep(i)) + 0.03 * peatdep(i) + 0.0134

          else if (ewtable(i) >= 0.0 .and. ewtable(i) < 0.30) then ! within the first 30 cm of surface

            ratescpo(i) = 0.009 * (1. - exp(-20. * ewtable(i))) + 0.015 * ewtable(i)
            ratescpa(i) = abs(0.009 * exp(-20. * ewtable(i)) - 0.183 * exp(-18. * peatdep(i)) - 0.015 * ewtable(i) + 0.0044) ! FLAG was missing the abs. JM Oct 2020.

          else if (ewtable(i) >= 0.30) then ! deeper in the soil column

            ratescpo(i) = 0.0134 - 0.183 * exp(-18. * ewtable(i)) + 0.003 * ewtable(i)
            if (zbotw(i,indexLastPermLayer) - ewtable(i) > 0.) then ! water table is within soil permeable depth 
              ratescpa(i) = -0.183 * exp(-18. * peatdep(i)) + 0.003 * (peatdep(i) - wtable(i)) &
                          + 0.183 * exp(-18. * ewtable(i))
            else 
              ratescpa(i) = 0.
            end if 
            ! ratescpa(i)=-0.183*exp(-18*peatdep(i))+0.003*(peatdep(i)-wtable(i))+0.183*exp(-18*ewtable(i))-0.004504   ! for continuity
          end if
        else if (ipeatland(i) == 2) then ! fens
          if (ewtable(i) < 0.0) then ! flooded

            ratescpo(i) = 0.0
            ratescpa(i) = 0.01512 - 1.12 * exp(-25. * peatdep(i))

          else if (ewtable(i) >= 0.0 .and. ewtable(i) < 0.30) then ! within the first 30 cm of surface

            ratescpo(i) = -0.01 * exp(-40. * ewtable(i)) + 0.015 * ewtable(i) + 0.01
            ratescpa(i) = abs(-0.01 * exp(-40. * ewtable(i)) - 1.12 * exp(-25. * &
                          peatdep(i)) + 0.015 * ewtable(i) + 0.005119)

          else if (ewtable(i) >= 0.30) then ! deeper in the soil column

            ratescpo(i) = 0.01512 - 1.12 * exp(-25 * ewtable(i))
            if (zbotw(i,indexLastPermLayer) - ewtable(i) > 0.) then ! water table is within soil permeable depth 
              ratescpa(i) = -1.12 * (exp(-25. * peatdep(i)) - exp(-25. * ewtable(i)))
            else 
              ratescpa(i) = 0.
            end if 
          end if
        end if

        !>  Convert respiration rates from kg c/kg c.year to u-mol co2/kgC/s
        ratescpo(i) = 2.64 * ratescpo(i)
        ratescpa(i) = 2.64 * ratescpa(i)

        !> Find the carbon storage in oxic and anoxic compartments (Cso. Csa)
        !! The water table depth delineates the oxic and anoxic compartments.
        !! functions (R**2 = 0.9999) determines the carbon content of each
        !! compartment from a peat bulk density profile based on unpulished
        !! data from P.J.H. Richard (described in fig. 1, Frokling et al.(2001)
        !! conversion of peat into carbon with 48.7% (Mer Bleue unpublished data,
        !! Moore)
        !! The oxic portion is limited by the total peat depth in case the water 
        !! table is below the peat and in the permeable soil depth.

        Cso(i) = peatStorage(min(ewtable(i),peatdep(i)))
        Csa(i) = peatStorage(peatdep(i)) - Cso(i)
        
        !> Find the soil respiration rate in Cso and Csa umol/m2/s.
        !! Moisture multiplier (0.025) indicates rate reduction in decomposition due
        !! to anoxia (Frolking et al. 2001), only applied to anoxic layer

        resoxic(i)   = ratescpo(i) * Cso(i) * fto(i)
        resanoxic(i) = ratescpa(i) * Csa(i) * fta(i) * 0.025
        socres_peat(i)   = resoxic(i) + resanoxic(i)

      end if 
    end do ! i loop
    return

  end subroutine soilCResPeat
  !! @}
  ! ---------------------------------------------------------------

  !> \ingroup heterotrophicrespirationmod_updatePoolsHetResp
  !! @{ Find vegetation and tile averaged litter and soil C respiration rates
  !! using values from canopy over ground and canopy over snow subareas.
  !! Also adds the moss and peat soil respiration to the tile level quantities.
  !! Next the litter and soil C pools are updated based on litter and soil C respiration rates.
  !! The humidified litter is then transferred to the soil C pool.
  !! Soil respiration is estimated as the sum of heterotrophic respiration
  !! and root maintenance respiration. For peatlands, we additionally add 
  !! moss values to the grid (litter respiration
  !! and moss root respiration).
  !> @author Vivek Arora, Joe Melton
  subroutine updatePoolsHetResp (il1, il2, ilg, fcancmx, ltresveg, scresveg, & ! In
                                 ipeatland, fg, litresmoss, socres_peat, & ! In
                                 sort, spinfast, rmrveg, rmr, leapnow, & ! In
                                 useTracer, ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer, & ! In
                                 litrmass, soilcmas, Cmossmas, litrmsmoss, peatSoilC, & ! In / Out
                                 tracerLitrMass, tracerSoilCMass, tracerMossCMass, & ! In / Out
                                 hetrsveg, litres, socres, hetrores, humtrsvg, soilresp, & ! Out
                                 humiftrs, litrfallmoss, ltrestepmoss, & ! Out
                                 humstepmoss, socrestep) ! Out

    use classicParams, only : icc, iccp1, iccp2, ignd, zero, humicfac, humicfac_bg, &
                              deltat, humicfacmoss, rmortmoss

    implicit none

    ! Inputs
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    integer, intent(in) ::  spinfast        !< spinup factor for soil carbon whose default value is 1. as this factor increases the
                                            !< soil c pool will come into equilibrium faster. reasonable value for spinfast is
                                            !< between 5 and 10. when spinfast/=1 then the balcar subroutine is not run.
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
                                    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
                                    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
                                    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    integer, intent(in) :: ipeatland(ilg)       !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    logical, intent(in) :: leapnow            !< true if this year is a leap year. Only used if the switch 'leap' is true.
    real, intent(in) :: fcancmx(ilg,icc)      !< max. fractional coverage of CTEM's pfts, but this can be
                                              !< modified by land-use change, and competition between pfts
    real, intent(in) :: ltresveg(ilg,iccp2,ignd)     !< fluxes for each pft: litter respiration for each pft + bare fraction
    real, intent(in) :: scresveg(ilg,iccp2,ignd)     !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    real, intent(in) :: fg(ilg)             !< Fraction of grid cell that is bare ground.
    integer, intent(in) :: sort(icc)
    real, intent(in) :: litresmoss(ilg)    !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: socres_peat(ilg)   !< heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: rmrveg(ilg,icc)    !< Maintenance respiration for root for the CTEM pfts in u mol co2/m2. sec
    real, intent(in) :: rmr(ilg)           !< Root maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: ltResTracer(ilg,iccp2,ignd)  !< Tracer fluxes for each pft: litter respiration for each pft + bare fraction
    real, intent(in) :: sCResTracer(ilg,iccp2,ignd)  !< Tracer soil carbon respiration for the given sub-area in umol co2/m2.s,for ctem's pfts
    real, intent(in) :: litResMossTracer(ilg)    !< Tracer moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) TODO: Units
    real, intent(in) :: soCResPeatTracer(ilg)    !< Tracer heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    ! Updates
    real, intent(inout) :: litrmass (ilg,iccp2,ignd)  !< Litter mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(inout) :: soilcmas(ilg,iccp2,ignd)   !< Soil carbon mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(inout) :: Cmossmas(ilg)              !< Moss biomass C pool (kgC/m2)
    real, intent(inout) :: litrmsmoss(ilg)           !< Moss litter C (kgC/m2)
    real, intent(inout) :: peatSoilC(ilg)              !< Peat soil C pool (kgC/m2)
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout) :: tracerMossCMass(:)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$

    ! Outputs
    real, intent(out) :: litrfallmoss(ilg)  !< moss litter fall (kgC/m2/timestep)
    real, intent(out) :: ltrestepmoss(ilg)  !< litter respiration from moss (kgC/m2/timestep)
    real, intent(out) :: humstepmoss(ilg)   !< moss humification (kgC/m2/timestep)
    real, intent(out) :: socrestep(ilg)     !< heterotrophic respiration from soil (kgC/m2/timestep)
    real, intent(out) :: hetrsveg(ilg,iccp1) !< Vegetation averaged litter and soil C respiration rates (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: litres(ilg)      !< Litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: socres(ilg)      !< Soil carbon respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: hetrores(ilg)    !< Heterotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: humtrsvg(ilg,iccp2,ignd)     !< transfer of humidified litter from litter to soil C pool per PFT. (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: soilresp(ilg)    !< Soil respiration. This includes root respiration and respiration from
                                          !! litter and soil carbon pools. Note that soilresp is different from
                                          !! socres,which is respiration from the soil C pool.(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: humiftrs(ilg)    !< Transfer of humidified litter from litter to soil C pool (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    ! Local
    integer :: k, j, i
    real :: ltrestep !<
    real :: screstep !<
    real :: hutrstep(ilg,iccp2,ignd) !<
    real :: tracerhutrstep !<
    real :: soilrsvg(ilg,iccp2) !<

    ! -----------
    ! initialize all outputs to zero:
    hetrsveg = 0.0 
    soilresp = 0.0
    humiftrs = 0.0
    litrfallmoss = 0.0
    ltrestepmoss = 0.0
    humstepmoss = 0.0
    socrestep = 0.0
    litres = 0.0
    socres = 0.0
    hetrores = 0.0
    hutrstep = 0.0
    humtrsvg = 0.0
    soilrsvg = 0.0

    do i = il1,il2
      ! Quick warning check to make sure only fens and bogs are considered as presently 
      ! formulated. If other systems are considered (such as moss on uplands), this scheme 
      ! will need to be adapted so stop the user now.
      if (ipeatland(i) > 2) then 
        print*,'Error: updatePoolsHetResp - not yet setup for ipeatland > 2'
        call errorHandler('updatePoolsHetResp', - 1)
      end if 
    end do

    !> We have litter and soil C respiration from heterotrophicRespiration, if we
    !! sum them we get hetrsveg (PFT-level + bareground), also adding root respiration 
    !! we can get soil respiration.
    do j = 1,iccp1! loop 340
      do i = il1,il2
        if (j < iccp1) then
          if (fcancmx(i,j) > zero) then ! vegetated 
            do k = 1,ignd
              ! hetrsveg is kept per PFT and tile (not per layer) at the moment.
              ! This calculation will include both mineral and peatland tiles.
              ! For peatland tiles, iccp1 (bareground) will have a value of zero. 
              hetrsveg(i,j) =  hetrsveg(i,j) + ltresveg(i,j,k) + scresveg(i,j,k)            
            end do            
            soilrsvg(i,j) = hetrsveg(i,j) + rmrveg(i,j)
          end if 
        else 
          if (fg(i) > zero) then ! bareground 
            do k = 1,ignd
              ! hetrsveg is kept per PFT and tile (not per layer) at the moment.
              ! This calculation will include both mineral and peatland tiles.
              ! For peatland tiles, iccp1 (bareground) will have a value of zero. 
              hetrsveg(i,j) =  hetrsveg(i,j) + ltresveg(i,j,k) + scresveg(i,j,k)            
            end do
            !bare ground has no roots respiration.        
            soilrsvg(i,j) = hetrsveg(i,j)
          end if           
        end if
      end do ! loop 350
    end do ! loop 340


    !> Find tile averaged litter and soil C respiration rates. 
    !> We add over the vegetated areas and the bareground (iccp1) values
    !! to the grid sum if it's not peatland.
    !! If it is a peatland, we also add litresmoss and socres_peat to the grid sum but 
    !! no bareground values as we assume peatlands have no bareground (and thus no 
    !! flux was calculated). That happens below this do loop.
    !!  As well,we don't add the LUC product pools as they 
    !! are considered a LUC flux, and thus not a heterotrophic respiration flux.
    do j = 1,iccp1 ! loop 360
      do i = il1,il2
        if (j < iccp1) then 
          soilresp(i) = soilresp(i) + fcancmx(i,j) * soilrsvg(i,j)
          do k = 1,ignd
            litres(i) = litres(i) + fcancmx(i,j) * ltresveg(i,j,k)
            socres(i) = socres(i) + fcancmx(i,j) * scresveg(i,j,k)            
          end do
        else ! bareground areas                          
          soilresp(i) = soilresp(i) + fg(i) * soilrsvg(i,j)     
          do k = 1,ignd
            litres(i) = litres(i) + fg(i) * ltresveg(i,j,k)
            socres(i) = socres(i) + fg(i) * scresveg(i,j,k)
            end do
        end if 
      end do ! loop 370
    end do ! loop 360

    ! Adjust for peatland contributions and also calculate moss litterfall.    
    do i = il1,il2
      if (ipeatland(i) > 0) then
    
        ! Both peatland fluxes (moss litter respiration and peat soil carbon respiration)
        ! are not tracked per layer so we just can add the tile values.
        litres(i) = litres(i) + litresmoss(i) ! add the moss litter, which is assumed to cover whole tile.
        
        ! Calculate moss timestep C fluxes, '/365*deltat' converts per year
        ! to per timestep, 'deltat/963.62' converts umol CO2/m2/s to kgC/m2/deltat.
        ltrestepmoss(i) = litresmoss(i) * (1.0 / 963.62) * deltat   ! kgC/m2/dt

        ! Calculate moss timestep C fluxes, '/365*deltat' converts per year
        ! to per timestep, 'deltat/963.62' converts umol CO2/m2/s to kgC/m2/deltat.
        if (leapnow) then
          litrfallmoss(i) = Cmossmas(i) * rmortmoss / 366. * deltat ! kgC/m2/day(dt)
          ! if (useTracer > 0) FLAG,not connected up !
        else
          litrfallmoss(i) = Cmossmas(i) * rmortmoss / 365. * deltat ! kgC/m2/day(dt)
          ! if (useTracer > 0) FLAG,not connected up !
        end if

        ! we now add the peat value to any from the vegetation.
        socres(i) = socres(i) + socres_peat(i) 
                    
        ! Soil respiration at the tile level is similarly updated with the contribution 
        ! from moss litter and peat respiration.
        soilresp(i) = soilresp(i) + litresmoss(i) + socres_peat(i)

      end if 
    end do 

    ! Find the tile sum hetrores based on the litres and socres values
    ! and adjust the units for some output variables.
    do i = il1,il2 ! loop 380
      hetrores(i) = litres(i) + socres(i)
      socrestep(i) = socres(i) * (1.0 / 963.62) * deltat     ! kgC/m2/dt
      soilresp(i) = soilresp(i) * (1.0 / 963.62) * deltat   ! kgC/m2/dt
    end do ! loop 380

    !> Update the litter and soil C pools based on litter and soil C respiration rates
    !! found above. Also transfer humidified litter to the soil C pool.
    do j = 1,iccp2 ! loop 420
      do i = il1,il2 ! loop 435
        do k = 1,ignd
  
          !> Convert u mol co2/m2.sec -> \f$(kg C/m^2)\f$ respired over the model time step
          ltrestep = ltresveg(i,j,k) * deltat / 963.62
          screstep = scresveg(i,j,k) * deltat / 963.62
  
          !> Update litter and soil C pools
          if (j < iccp1) then ! vegetated areas 
            litrmass(i,j,k) = litrmass(i,j,k) - (ltrestep * (1.0 + humicfac(sort(j))))
            hutrstep(i,j,k) = humicfac(sort(j)) * ltrestep
          else ! bare ground pools 
            !> Adjust bareground and LUC pool litter mass and humification (non-peatlands).
            litrmass(i,j,k) = litrmass(i,j,k) - (ltrestep * (1.0 + humicfac_bg))
            hutrstep(i,j,k) = humicfac_bg * ltrestep              
          end if
          
          ! Find the litter transfered (humified) to soil C pool in u-mol co2/m2.sec
          humtrsvg(i,j,k) = hutrstep(i,j,k) * (963.62 / deltat)
  
          ! Update the soil C pool for the inputs from humification and losses to respiration 
          soilcmas(i,j,k) = soilcmas(i,j,k) + real(spinfast) * (hutrstep(i,j,k) -  screstep)
  
          ! Don't allow negative values 
          if (litrmass(i,j,k) < zero) litrmass(i,j,k) = 0.0
          if (soilcmas(i,j,k) < zero) soilcmas(i,j,k) = 0.0
  
          ! Now perform the same calculations to the tracer pools
          if (useTracer > 0) then
            if (j < iccp1) then ! vegetated
              tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) - (ltResTracer(i,j,k) * deltat / 963.62 &
                                                              * (1.0 + humicfac(sort(j))))
              tracerhutrstep = humicfac(sort(j)) * ltResTracer(i,j,k) * deltat / 963.62
            else ! bare ground
              if (ipeatland(i) == 0) then ! non-peatlands
                tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) - (ltResTracer(i,j,k) * deltat / 963.62 &
                                                               * (1.0 + humicfac_bg))
                tracerhutrstep = humicfac_bg * ltResTracer(i,j,k) * deltat / 963.62
              end if
            end if
  
            tracerSoilCMass(i,j,k) = tracerSoilCMass(i,j,k) + real(spinfast) &
                                            * (tracerhutrstep -  sCResTracer(i,j,k) * deltat / 963.62)
  
            if (litrmass(i,j,k) < zero) tracerLitrMass(i,j,k) = 0.0
            if (soilcmas(i,j,k) < zero) tracerSoilCMass(i,j,k) = 0.0
          end if
        end do
      end do ! loop 430
    end do ! loop 420

    !> Find tile averaged humification and soil respiration rates. Also update
    !! the moss and peat soil C poools based on their humification and respiration 
    !! rates.    
    do i = il1,il2 ! loop 470
      do j = 1,iccp1  
        if (j < iccp1) then ! vegetated areas
          do k = 1,ignd
            humiftrs(i) = humiftrs(i) + fcancmx(i,j) * humtrsvg(i,j,k)
          end do
        else !> Now the bare ground areas for mineral soils. 
          do k = 1,ignd
            humiftrs(i) = humiftrs(i) + fg(i) * humtrsvg(i,j,k)
          end do    
        end if 
      end do 
      ! For peatlands, we additionally add moss values to the grid (litter respiration and moss root respiration). 
      if (ipeatland(i) > 0) then 
        humstepmoss(i) = humicfacmoss * ltrestepmoss(i)    ! kgC/m2/dt
        humiftrs(i)  = humiftrs(i) + humstepmoss(i) * (963.62 / deltat)! umol/m2/s    
        
        !> Update the moss litter and peat soil C pools based on their respiration rates
        litrmsmoss(i) = litrmsmoss(i) + litrfallmoss(i) - ltrestepmoss(i) - humstepmoss(i)     ! kg/m2
        
        ! change the units of socres_peat so it is in kgC/m2/timestep.
        peatSoilC(i)  = peatSoilC(i)  + real(spinfast) * (humstepmoss(i) &
                                                       - socres_peat(i) * (1.0 / 963.62) * deltat)
      end if 
      
    end do 
    
  end subroutine updatePoolsHetResp
  !! @}
  ! ------------------------------------------------------------------------------
  !> \namespace heterotrophicrespirationmod
  !! Central module for all heterotrophic respiration-related operations
  !!
  !! # Heterotrophic Respiration (Vegetated Ground)
  !!
  !! Heterotrophic respiration, \f$R_\mathrm{h}\f$ (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$), in CLASSIC is
  !! based on respiration from the litter (which includes contributions from the stem, leaf
  !! and root components), \f$R_{h, D}\f$, and soil carbon, \f$R_{h, H}\f$, pools,
  !!
  !! \f[ \label{hetres_all} R_\mathrm{h}=R_{h, D}+R_{h, H}. \hspace{10pt}[Eqn 1] \f]
  !!
  !! Heterotrophic respiration is regulated by soil temperature and moisture and is
  !! calculated on a daily time step. The original heterotrophic respiration scheme is
  !! described in Arora (2003) \cite Arora2003-3b7 while the modified parametrization used in CLASSIC
  !! is detailed in Melton and Arora (2014) \cite Melton2014-xy.
  !! Respiration from the litter and soil carbon pools takes the following basic form
  !!
  !! \f[ R_{\mathrm{h}, i} = 2.64 \times 10^{-6}\, \varsigma_i C_i f_{15}(Q_{10}) f(\Psi)_i f(z),
  !! \nonumber \\ i = \mathrm{D}, \mathrm{H}. \hspace{10pt}[Eqn 2] \f]
  !!
  !! The soil carbon and litter respiration depends on the amount of carbon in these components
  !! (\f$C_\mathrm{H}\f$ and \f$C_\mathrm{D}\f$; \f$kg\, C\, m^{-2}\f$) and on a PFT-dependent
  !! respiration rate specified at \f$15\, ^{\circ}{C}\f$ (\f$\varsigma_\mathrm{H}\f$ and
  !! \f$\varsigma_\mathrm{D}\f$; \f$kg\, C\, (kg\, C)^{-1}\, yr^{-1}\f$; see also
  !! classicParams.f90). The constant \f$2.64 \times 10^{-6}\f$ converts units from
  !! \f$kg\, C\, m^{-2}\, yr^{-1}\f$ to \f$mol\, CO_2\, m^{-2}\, s^{-1}\f$.
  !!
  !! The effect of soil moisture is accounted for via dependence on soil matric
  !! potential (\f$f(\Psi)\f$), described later. The temperature dependency of
  !! microbial soil respiration rates has been estimated by several different
  !! formulations, ranging from simple \f$Q_{10}\f$ (exponential) to Arrhenius-type
  !! formulations (see review by Lloyd and Taylor (1994) \cite Lloyd1994-ct). In CLASSIC, soil temperature
  !! influences heterotrophic respiration through a temperature-dependent
  !! \f$Q_{10}\f$ function (\f$f_{15}(Q_{10})\f$). The value of \f$Q_{10}\f$
  !! itself is assumed to be a function of temperature following a hyperbolic
  !! tan function:
  !!
  !! \f[ Q_{10} = 1.44 + 0.56\, \tanh[0.075 (46.0 - T_i)], \nonumber\\ i
  !! = \mathrm{D}, \mathrm{H}, \hspace{10pt}[Eqn 3]\f]
  !!
  !! where \f$T_{\{D, H\}}\f$ is the temperature of either the litter or soil
  !! carbon pool (\f$C\f$), respectively. The parametrization is a compromise
  !! between the temperature-independent \f$Q_{10}\f$ commonly found in many
  !! terrestrial ecosystem models (Cox, 2011) \cite Cox2001-am and the temperature-dependent
  !! \f$Q_{10}\f$ of Kirschbaum (1995) \cite Kirschbaum1995-db. While a constant \f$Q_{10}\f$ yields
  !! an indefinitely increasing respiration rate with increasing temperature, the
  !! formulation of Kirschbaum (1995) \cite Kirschbaum1995-db gives a continuously increasing
  !! \f$Q_{10}\f$ under decreasing temperature, which leads to unreasonably high
  !! soil and litter carbon pools at high latitudes. The CLASSIC
  !! parametrization avoids these issues with a \f$Q_{10}\f$ value of about 2.0
  !! for temperatures less than \f$20\, ^{\circ}C\f$, while a decreasing value of
  !! \f$Q_{10}\f$ at temperatures above \f$20\, ^{\circ}C\f$ ensures that the
  !! respiration rate does not increase indefinitely. As the soil temperature decreases below
  !! \f$ T_{crit} \f$ (typically \f$1\, ^{\circ}C\f$, see classicParams.f90) a step function is applied to the \f$f_{15}(Q_{10})\f$
  !! function to reflect lost mobility of the microbial populations due to less liquid water as:
  !! \f[ f_{15}(Q_{10, (T_i < T_{crit})}) = f_{15}(Q_{10}) * 0.1 \f]
  !!
  !!   \image html "Q10_response_sm.png" "Q10 response"
  !!
  !! The soil detrital pools are explictly tracked per soil layer.
  !!
  !! The response of heterotrophic respiration to soil moisture is formulated through
  !! soil matric potential (\f$\Psi\f$; \f$MPa\f$). While soil matric potential values
  !! are usually negative, the formulation uses absolute values to allow its logarithm
  !! to be taken. Absolute values of soil matric potential are high when soil is dry
  !! and low when it is wet. The primary premise of soil moisture control on heterotrophic
  !! respiration is that heterotrophic respiration is constrained both when the soils
  !! are dry (due to reduced microbial activity) and when they are wet (due to impeded
  !! oxygen supply to microbes) with optimum conditions in-between. The exception is the
  !! respiration from the litter component of the first soil layer, which is assumed to be continually exposed
  !! to air, and thus never oxygen deprived, even when soil moisture content is high
  !! (\f$0.04 > \vert \Psi \vert \geq \vert \Psi_{sat} \vert\f$, where \f$\Psi_{sat}\f$
  !! is the soil matric potential at saturation). The soil moisture dependence for each
  !! soil layer thus varies between 0 and 1 with matric potential as follows:
  !!
  !! for \f$0.04 > \vert\Psi_i\vert \geq \vert\Psi_{sat, i}\vert\f$
  !!
  !! \f[ f(\Psi_i)_\mathrm{H, (D, i>1)} = 1 - 0.5  \frac{\log(0.04)-\log\vert\Psi_i\vert}
  !! {\log(0.04)-\log\vert\Psi_{sat, i}\vert} \hspace{10pt}[Eqn 4]\f]
  !!
  !! \f[f(\Psi_i)_{D, i=1} = 1\nonumber; \hspace{10pt}[Eqn 5]\f]
  !!
  !! for \f$0.06 \geq \vert\Psi_i\vert \geq 0.04\f$
  !! \f[ f(\Psi_i)_{D, H} = 1; \hspace{10pt}[Eqn 6]\f]
  !!
  !! for \f$100.0 \gt \vert\Psi_i\vert > 0.06\f$
  !! \f[ f(\Psi_i)_{D, H} = 1 - 0.8\frac{\log\vert\Psi_i\vert-\log(0.06)}{\log(100)-\log(0.06)}; \hspace{10pt}[Eqn 7]\f]
  !!
  !! for \f$\vert\Psi_i\vert >= 100.0\f$
  !! \f[ f(\Psi_i)_{D, H}=0.2. \hspace{10pt}[Eqn 8]\f]
  !!
  !! Respiration also is reduced at depth in soil (\f$f(z)\f$) following Lawrence et al. (2015)
  !! \cite Lawrence2015-tj. This term is meant to represent unresolved depth dependent soil
  !! processes (such as oxygen availability, microbial community changes, etc.). The reduction
  !! in respiration with depth per soil layer is dependent upon the layer depth and
  !! a term, \f$z_t\f$, which is given a value of 10.0 (see classicParams.f90) as,
  !!
  !! \f[ f(z_i) =\exp (-z_i / z_t) \hspace{10pt}[Eqn 9]\f]
  !!
  !!   \image html "decr_resp_wit_depth_sm.png" "Decrease in respiration with depth"
  !!
  !! # Heterotrophic respiration for bare ground  
  !!
  !! The carbon contributions to the bare ground litter and soil carbon pools come via processes
  !! such as creation of bare ground due to fire, competition between PFTs and land use
  !! change. The heterotrophic respiration is sensitive to temperature and moisture in
  !! the same manner as vegetated areas using Eqs. (2)--(8). The
  !! base respiration rates of \f$\varsigma_{D, bare}\f$ and \f$\varsigma_{H, bare}\f$ are
  !! set to 0.5605 and 0.02258 kg C (kg C)\f$^{-1}\f$ yr\f$^{-1}\f$, respectively.
  !!
  !! The amount of humidified litter, which is transferred from the litter to the soil
  !! carbon pool (\f$C_{\mathrm{D} \rightarrow \mathrm{H}}\f$) is modelled as a fraction
  !! of litter respiration (\f$R_{h, D}\f$) as
  !!
  !! \f[  C_{\mathrm{D} \rightarrow \mathrm{H}} = \chi\, R_{h, D} \hspace{10pt}[Eqn 10] \f]
  !!
  !! where \f$\chi\f$ (see also classicParams.f90) is the PFT-dependent humification factor
  !! and varies between 0.4 and 0.5. For crops, \f$\chi\f$ is set to 0.1 to account for
  !! reduced transfer of humidified litter to the soil carbon pool which leads to loss in
  !! soil carbon when natural vegetation is converted to croplands. Over the bare ground
  !! fraction \f$\chi\f$ is set to 0.45.
  !!
  !! With heterotrophic respiration known, net ecosystem productivity (\f$NEP\f$) is
  !! calculated as
  !! \f[ NEP = G_{canopy} - R_\mathrm{m} - R_\mathrm{g} - R_\mathrm{h}. \hspace{10pt}[Eqn 11] \f]
  !!
  !! # Peatland soil heterotrophic respiration
  !!
  !!   \image html "soilCDiagram.png" "Relation between mineral soil and peatland soil detrital pools"
  !!
  !! Over the non-peatland fraction, HR is calculated
  !! as the sum of the respiration from litter and soil carbon pools as described above. 
  !! In peatlands a large amount of humic soil is generally
  !! located in the permanently saturated zone and the bulk density increases
  !! with soil depth (Loisel and Garneau, 2010) \cite Loisel2010-dj. Thus, the assumption of
  !! exponentially decreasing distribution of C content with increasing soil
  !! depth is not valid in peatlands (as was done in prior versions of CTEM). We used a quadratic equation to calculate
  !! the distribution of soil C content over depth based on an empirically
  !! determined bulk density profile (Frolking et al., 2001) \cite Frolking2001d42.
  !!
  !! HR over the peatland fraction of a grid cell is modelled using a two-pool
  !! approach with a flexible boundary between the pools that depends on the
  !! depth of the water table:
  !!
  !! \f$ R_{o}=C_{SOM, o}k_{o}f_{T, {o}} \hspace{10pt}[Eqn 12] \f$ 
  !!
  !! \f$ R_{a}=C_{SOM, a}k_{a}f_{T, {a}}f_{anoxic}  \hspace{10pt}[Eqn 13] \f$  
  !!
  !! where \f$o\f$ and \f$a\f$ denote the oxic and anoxic portions of the soil C pool respectively. The respiration rate \f$R\f$ (unit:
  !! \f$\mu\f$ mol C m\f$^{-2}\f$/s) is obtained from the respiration rate
  !! coefficient \f$k\f$ (\f$\mu\f$ mol C / kg C / s), 
  !!
  !! \f$ k_{o}= 0 \f$ for \f$ z_{wt}<0  \hspace{10pt}[Eqn 14]\f$
  !!
  !!  \f$ k_{o}= k_{1}\left( 1-e^{k_{2}z_{wt}} \right)+k_{3}z_{wt} \f$ for \f$0.3 > z_t \ge 0  \hspace{10pt}[Eqn 15]\f$
  !!
  !!  \f$ k_{o}= k_{4}e^{k_{5}z_{wt}}+k_{6}z_{wt}+k_{7} \f$ for \f$ z_{wt}\ge 0.3 \hspace{10pt}[Eqn 16]\f$
  !!
  !! \f$ k_{a}= k_{4}e^{k_{5}z_{p}}+{10k}_{6}z_{p}+k_{7} \f$ for \f$ z_{wt}<0 \hspace{10pt}[Eqn 17]\f$
  !!
  !! \f$ k_{a} = \f$ | \f$ k_1 e^{ k_2 z_{wt} } - k_4 e^{k_5 z_p} -k_3 z_{wt}+k_8 | \f$ for \f$0.3>z_{wt} \ge 0 \hspace{10pt}[Eqn 18]\f$
  !!
  !! \f$ k_{a}=  k_{4}{(e}^{k_{5}z_{P-}}e^{k_{5}z_{wt}})+k_{6}\left( z_p-z_{wt} \right)\f$ for \f$ z_{wt}\ge 0.3 \hspace{10pt}[Eqn 19]\f$
  !!
  !! where the values of \f$k_{1}\f$ to \f$k_{8}\f$ differ between fens and bogs. The variation of \f$k_{o}\f$ and \f$k_{a}\f$ with water table depth for
  !! bogs and fens is shown below (Fig. 2 of (Wu et al. 2016) \cite Wu2016-zt).
  !!
  !!   \image html "peatlandDecompDepth_Wu2016Fig2.png" "Variation of peatland respiration rate coefficients with water table depth."
  !!
  !! It will be noted that there is a sharp
  !! transition in decomposition rate at a depth of 0.3 m, reflecting the work of
  !! Frolking et al. (2001) \cite Frolking2001d42. As noted above, this value is widely accepted
  !! as a representative estimate of the depth dividing the acrotelm and catotelm.
  !! In reality, of course, this depth will vary among peatlands.
  !!
  !! The water table depth \f$z_{wt}\f$ is deduced by searching for a soil layer
  !! below, which the soil is saturated and above which the soil moisture is at or
  !! below the retention capacity with respect to gravitational drainage. Within
  !! this soil layer \f$j\f$, \f$z_{wt}\f$ is calculated in the generalutils.f90 function findWaterTable as
  !!
  !! \f$ z_{wt}=z_{{b}, j}-\Delta z\left[
  !! \frac{\theta_{{l}, j}+\theta_{i, j}-\theta
  !! _{{ret}, j}}{\theta_{{p}, j}-\theta_{{ret}, j}} \right] \hspace{10pt}[Eqn 20]\f$
  !!
  !! where \f$\Delta z \f$ is the thickness of soil layer (unit: m),
  !! \f$\theta_l \f$ and \f$\theta_i\f$ are the liquid and frozen water contents
  !! (unit, m\f$^3\f$ / m\f$^3\f$), \f$\theta_{ret}\f$ and \f$\theta_p\f$ are the
  !! water retention capacity and the porosity, and \f$z_b\f$ (unit: m) is the
  !! bottom depth of the soil layer.
  !!
  !! Also affecting the respiration rate (Eqns 12 and 13) is the temperature functions \f$f_T\f$, 
  !!
  !!  \f$ f_{T, o} =Q_{10, o}^{\left(\int\limits_{0}^{z_{wt}} T_{j} -15\right)/10}  \hspace{10pt}[Eqn 21]\f$
  !!
  !!  \f$ f_{T, {a}} =Q_{10, {a}}^{\left(\int\limits_{z_{wt}}^{z_\mathrm{p}} T_{j} -15\right)/10}  \hspace{10pt}[Eqn 22]\f$
  !!
  !! \f$Q_{10}\f$ is calculated using
  !! a hyperbolic tan function of the soil temperatures (\f$T_s)\f$ of the oxic
  !! and anoxic zones similar to the mineral soil calculation (Melton and Arora, 2016) \cite Melton2016-zx, which are in turn functions of
  !! water table depth. The \f$Q_{10}\f$ values of the anoxic and the oxic
  !! zones of the soil are indicated as \f$Q_{10, a}\f$ and \f$Q_{10, o}\f$. The anoxic and oxic soil temperatures are found via
  !!
  !! \f$ T_{s, o}=\int\limits_{0}^{z_{wt}} T_{j} \, /(z_{wt})  \hspace{10pt}[Eqn 23]\f$
  !!
  !! \f$ T_{s, a}=\int\limits_{z_{wt}}^{z_{p}} T_{j} /(z_{p}-z_{wt})  \hspace{10pt}[Eqn 24]\f$
  !!
  !! The final terms of the respiration rate (Eqns 12 and 13) are the soil C mass \f$C_{SOM}\f$ (kg),
  !! which is calculated in the peatStorage function of peatlandsMod.f90 as,
  !!
  !! \f$ C_{SOM, o}= 0.487\ast (k_{9}z_{wt}^{2}+k_{10}z_{wt})  \hspace{10pt}[Eqn 25]\f$
  !!
  !! where 0.487 is a parameter that converts from soil mass to soil C content.
  !!
  !!  \f$ C_{SOM, a}= C_{SOM}-C_{SOM, o} \hspace{10pt}[Eqn 26]\f$
  !!
  !! and a scaling factor
  !! \f$f_{anoxic}\f$ after Frolking et al.~(2010, 2001) \cite Frolking2001d42 \cite Frolking2010-jq, which represents the
  !! inhibition of microbial respiration under anoxic conditions. The value of
  !! this parameter is uncertain, varying in those two papers between 0.001, 0.025
  !! and 0.1. Based on calibration runs using two of the data sets described in Wu et al. (2016) \cite Wu2016-zt
  !! (MB-Bog and AB-Fen), we adopted a value of 0.025. 
  !! The values of \f$k\f$, \f$f_T\f$, and \f$C_{SOM}\f$ are updated along with the
  !! water table depth (\f$z_{wt}\f$, unit: m, positive downward) and the peat
  !! depth (\f$z_p\f$, unit: m) at each biogeochemistry time step. The equations for \f$k\f$
  !! and \f$C_{SOM}\f$ are derived from Fig. 2 in Frolking et al. (2001) \cite Frolking2001d42, and
  !! parameterized differently for fens and bogs.
  !!
  !! As only organic soil is considered in peatlands, the peat soil C is updated
  !! from the humification (C\f$_{hum}\f$, kg C/ m\f$^2\f$ / day) and soil
  !! respiration from the oxic (\f$R_o\f$ in kg C/ m\f$^2\f$ / day) and
  !! anoxic (\f$R_a\f$ in kg C/ m\f$^2\f$ / day) components during the
  !! time step:
  !!
  !! \f$ \frac{{dC}_{SOM}}{dt} = C_{hum}-R_{o}-R_{a}\hspace{10pt}[Eqn 27] \f$
  !!
  !! \f$C_{hum}\f$ is calculated as a PFT-dependent fraction of the decomposition
  !! rate. At the end of each time step, the peat depth (i.e. the depth
  !! of the organic soil) \f$z_p\f$ is updated from the updated peat C mass
  !! (C\f$_{SOM}\f$ in kg) by solving the quadratic equation
  !!
  !! \f$ z_{p}=
  !! \frac{{-k}_{10}+\sqrt{k_{10}+\frac{4k_{9}{C}_{SOM}}{0.487}}.
  !! }{2k_{9}} \hspace{10pt}[Eqn 28]\f$
  !!
  !! in the peatDepth function of peatlandsMod.f90
end module heterotrophicRespirationMod
