!> \file
!> Peatland specific parameterizations (photosynthesis and heterotrophic respiration)
module peatlandsMod

  ! J. Melton. Sep 26, 2016

  implicit none

  ! Subroutines contained in this module:
  public  :: mossPht
  public  :: updateMossC
  public  :: peatDayEnd
  public  :: peatDepth

contains

  ! ------------------------------------------------------------------

  !> \ingroup peatlandsmod_mossPht
  !! @{
  !> Moss photosynthesis subroutine (equations are in module-level description)
  !> @author Yuanqiao Wu, J. Melton
  subroutine mossPht (il1, il2, iday, qswnv, thliq, co2conc, tsurfk, zsnow, &
                      pres, Cmossmas, dmoss, anmoss, rmlmoss, cevapmoss, &
                      ievapmoss, ipeatland, daylength, pdd)

    use classicParams,  only : rmlmoss25, tau25m, ektau, gasc, kc25, ko25, ec, ej, eo, evc, sj, &
                               hj, alpha_moss, thpmoss, thmmoss, ilg, ignd, TFREZ, RHOW, mossMaxSnow, &
                               mossMinTemp, pi

    implicit none

    ! arguments:

    integer, intent(in) ::  il1                     !<
    integer, intent(in) ::  il2                     !<
    integer, intent(in) ::  iday                    !< model day of year
    integer, intent(in) ::  ipeatland(ilg)          !< Flag for peatland tiles
    real, dimension(ilg), intent(in) :: qswnv       !< visible short wave radiation = qswnv in energBalNoVegSolve
                                                    !! and qswnvg in energBalVegSolve (W/m2)
    real, dimension(ilg,ignd), intent(in) :: thliq  !<
    real, dimension(ilg), intent(in) :: zsnow       !< Snow depth (m)
    real, dimension(ilg), intent(in) :: co2conc     !<
    real, dimension(ilg), intent(in) :: daylength   !<
    real, dimension(ilg), intent(in) :: tsurfk      !< grid average ground surface temprature in K
    real, dimension(ilg), intent(in) :: pres        !<
    real, dimension(ilg), intent(in) :: Cmossmas    !< unit kg moss C updated in ctem
    real, dimension(ilg), intent(in) :: dmoss       !< unit m, depth of living moss. assume = 2 cm
                                                    !! can be related to mmoss as a variable
    real, dimension(ilg), intent(inout) :: pdd(ilg) !< peatland degree days above 0 deg C.
    integer, dimension(ilg), intent(out) ::  ievapmoss   !< Value is 0 is no evaporation from moss, 1 otherwise
    real, dimension(ilg), intent(out) :: anmoss     !< net photosynthesis (umol CO2/m2/s)
    real, dimension(ilg), intent(out) :: rmlmoss    !< moss autotrophic respiration (umol CO2/m2/s)
    real, dimension(ilg), intent(out) :: cevapmoss  !< evaporation coefficent for moss surface

    ! Local variables

    integer :: i               !>
    integer :: j               !>
    real :: mmoss(ilg)       !<  dry moss biomasss (kg dry moss biomass)
    real :: parm(ilg)        !< par at the ground (moss layer) surface umol/m2/s
    real :: tsurf(ilg)       !< grid average ground surface temperature in C
    real :: wmoss(ilg)       !< water content extraporated from the surface
                             !! humidity qg and thliq of the first soil layer
                             !! unit kg water/ kg dw
    real :: wmosmin(ilg)     !< residual water content kg water /kg moss
    real :: wmosmax(ilg)     !< maximum water content kg water /kg moss
    real :: fwmoss(ilg)      !< relative water content of mosses in g fw /g dw
    real :: dsmoss(ilg)      !< degree of moss saturation = relative water
                             !! content/maximum relative water content
    real :: g_moss(ilg)      !< moss conductance umol CO2/m2/s (based on
                             !! Williams and Flanagan, 1998 for Sphagnum)
    real :: mwce (ilg)       !< moisture function of dark respiration of moss
    real :: tmoss(ilg)       !< moss temperature extraporated from the tbar 1
                             !! and grid averaged ground surface temperature, tsurf
    real :: tmossk(ilg)      !< moss temperature in K
    real :: q10rmlmos(ilg)   !< temperature function of the moss dark respiration
    real :: gamma(ilg)       !< compensation point for gross photosynthesis (Pa)
    real :: o2(ilg)          !< partial presure of oxygen (Pa)
    real :: co2a(ilg)        !< partical pressure of co2 (pa) same as in PHTSYN
    real :: tau(ilg)         !< arrhenius funciton of temperature
    real :: kc(ilg)          !< kinetic coeffficient of CO2 for photosynthesis
    real :: ko(ilg)          !< kinetic coeffficient of O2 for photosynthesis
    real :: bc(ilg)          !< coefficient used for wc
    real :: vcmax25(ilg)     !< seasonal varied maximum carboylation at 25
                             !! sphagnum (fig. 6, Williams and Flanagan, 1998)
    real :: vcmax(ilg)       !< max carboxylation rate (umol/m2/s)
    real :: jmax25(ilg)      !< maximum electorn transport rate at 25 degrees (umol/m2/s)
    real :: jmax(ilg)        !< maximum electorn transport rate (umol/m2/s)
    real :: wj(ilg)          !< net co2 assimilation rate limited by
                             !! electron transport (umol/m2/s)=jE in PHTSYN
    real :: wc(ilg)          !< net co2 assimilation rate limited by
    real :: ws(ilg)          !< net co2 assimilation rate limited by
                             !! sucrose availability (umol/m2/s)=JE IN PHTSYN
    real :: photon(ilg)      !< electron transport rate (umol/m2/s)
    real :: term1(ilg)       !< temporary terms for photosynthesis calculations
    real :: term2(ilg)       !< temporary terms for photosynthesis calculations
    real :: term3(ilg)       !< temporary terms for photosynthesis calculations
    real :: psna(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: psnb(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: psne(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: mI(ilg)          !< coefficients of the solutions for net psn
    real :: mII(ilg)         !< coefficients of the solutions for net psn

    ! Local parameters:
    real, parameter :: tref = 298.16    !< unit K

    !     PHOTOSYNTHESIS COUPLING OR CURVATURE COEFFICIENTS
    ! real, parameter :: BETA1 = 0.950
    ! real, parameter :: BETA2 = 0.990

    ! ...........................
    ! Begin calculations:

    do i = il1,il2 ! loop 100
      wj(i)     = 0.0
      ws(i)     = 0.0
      wc(i)     = 0.0
      anmoss(i) = 0.0
      rmlmoss(i) = 0.0
    end do ! loop 100

    do i =  il1,il2 ! loop 200

      if (ipeatland(i) > 0) then ! only do for regions with moss.
        
        tsurf(i) = tsurfk(i) - tfrez
        tmossk(i) = tsurfk(i)
        tmoss(i) = tsurf(i)

        !>    phenology water factor on mosses, grow when temperature > mossMinTemp
        !!    and snowpack < mossMaxSnow
        
        if (zsnow(i) <= mossMaxSnow .or. tsurf (i) >= mossMinTemp) then

          !> Find the light level (parm) at the ground surface for moss photosynthesis,
          !! and a scaling factor degree of saturation of the moss layer phenology
          !! parm is in umol/m2/s and converted from qswnv in W/m2
          parm(i) = qswnv(i) * 4.6
          o2(i)  = 20.9/100.0 * pres(i)
          co2a(i) = co2conc(i)/1000000.0 * pres(i)

          !> Water content used for the living moss (depth dmoss)
          !! dmoss is an input and site specific. Preferably make dmoss a function
          !! of Cmoss and ipeatland (different species in fens and bogs)
          !! observed range of wmoss: 5 to 40 in Robrek (2007, 2009), 5 to 25
          !! (Flanagen and Williams 1998)
          !! dmoss is between 2.5 to 5cm based on the species (Lamberty et al. 2006)

          ! Convert from moss C to moss dry biomass using a standard conversion of 0.46 kg C 
          ! per kg dry biomass
          mmoss(i) = Cmossmas(i)/0.46
          
          ! Find the water content of the moss layer based on the total layer content and 
          ! the depth of the living moss. It is constrained to be within the minumum and 
          ! maximum water holding capacity of the moss.
          wmoss(i) = thliq(i,1) * rhow/(mmoss(i)/dmoss(i))
          wmosmax(i) = min(45.0,thpmoss * dmoss(i) * rhow/mmoss(i)) 
          wmosmin(i) = max(5.0,thmmoss * dmoss(i) * rhow/mmoss(i))
          wmoss(i) = min(wmosmax(i),max(wmosmin(i),wmoss(i)))
          
          ! Lastly convert to the relative water content.
          fwmoss(i) = wmoss(i) + 1.     ! g fresh weight /g dry weight

          !> Find moss conductance (g_moss) in umol CO2/m2/s
          !! (Williams and Flanagan, 1998 for Sphagnum). follow MWM, fwmoss is
          !! the mosswat_fd in MWM. Empirical equation is only valid up to
          !! fwmoss=13, above 13 apply a linear extension to the equation.

          if (fwmoss(i) <= 13.0) then
            g_moss(i) = - 0.195 + 0.134 * fwmoss(i) - 0.0256 * (fwmoss(i)) &
                        ** 2 + 0.00228 * (fwmoss(i)) ** 3 - 0.0000984 * &
                        (fwmoss(i)) ** 4 + 0.00000168 * (fwmoss(i)) ** 5
          else
            g_moss (i) = - 0.000447 * fwmoss(i) + 0.0489
          end if

          g_moss (i) = g_moss(i) * 1000000.0
          g_moss(i) = max(0.0,g_moss(i))

          !> Find moss surface evaporation coefficient
          !! controled by the degree of saturation in moss, pass to energBalVegSolve and energBalNoVegSolve
          !! apply a similar equation of soil surface cevap in energyBudgetPrep
          !! CEVAP = 0.25*[1 – cos(THLIQ*pi/THFC)]^2

          ! Fist get degree of moss saturation
          dsmoss(i) = (wmoss(i) - wmosmin(i))/(wmosmax(i) - wmosmin(i))
          
          ! then dsmoss is used to set cevap to 0, 1, or somewhere in between.
          if (dsmoss(i) <   0.001) then
            ievapmoss(i) = 0
            cevapmoss(i) = 0.
          else if (dsmoss(i) >= 1.0) then
            ievapmoss(i) = 1
            cevapmoss(i) = 1.0
          else
            ievapmoss(i) = 1
            cevapmoss(i) = 0.25 * (1.0 - cos(pi * dsmoss(i))) ** 2
          end if


          !> Calculate the moss water content effect on dark respiration
          !! in MWM and PDM an optimal wmoss is at 5.8 gw/gdw(fig. 2e, Frolking et al., 1996)
          !! Recent studies show weak but significant increases of sphagnum dark respiration
          !! with moss water content above 5.8 gw/gdw (Adkinson and Humphreys, 2011 and ref.)
          !! this change has improved the ER simulation greatly
          if (wmoss(i) < 0.4) then
            mwce (i) = 0.0
          else if (wmoss(i) < 5.8 .and. wmoss(i) > 0.4) then
            mwce(i) = 0.35 * wmoss(i) ** (2.0/3.0) - 0.14
          else
            mwce(i) = 0.01 * wmoss(i) + 0.942
          end if

          !> Moss dark respiration
          !! observed range of rmlmoss 0.60 to 1.60 umol/m2/s (e.g. Adkinson 2006)

          q10rmlmos(i) = (3.22 - (0.046 * tmoss(i))) ** ((tmoss(i) - 25.0)/10.0)
          rmlmoss(i) = rmlmoss25 * mwce(i) * q10rmlmos(i)

          !> Moss photosynthesis
          !! calculate bc (coefficient used for Wc, limited by Rubisco)

          tau(i) = tau25m * exp((tmossk(i) - tref) * ektau/(tref * gasc * tmossk(i)))
          gamma(i) = 0.5 * o2(i)/ tau(i)
          kc(i) = kc25 * exp((tmossk(i) - tref) * ec/(tref * gasc * tmossk(i)))
          ko(i) = ko25 * exp((tmossk(i) - tref) * eo/(tref * gasc * tmossk(i)))
          bc(i)  = kc(i) * (1.0 + (o2(i)/ ko(i)))

          !    seasonal change of Vcmax, sphagnum (fig. 6, Williams and Flanagan, 1998)
          !    from May 1st to september 1st Vmax is maximum
          !         if (iday < 121)        then
          !              vcmax25(i) = 6.0
          !         else if (iday > 245)   then
          !              vcmax25(i) = 7.0
          !         else
          !              vcmax25(i) = 13.5
          !         end if

          !!     use a function that is also valid for the southern hemisphere
          if (daylength(i) > 14.0 .and. pdd(i) > 200. .and. pdd(i) < 2000. ) then
            vcmax25(i) = 14.0
          else
            vcmax25(i) = 6.5
          end if

          vcmax(i) = vcmax25(i) * exp((tmossk(i) - tref) * evc/(tref * gasc * tmossk(i)))

          !> Calculate ws (phototysnthesis rate limited by transport capacity)
          !! = js in photosynCanopyConduct
          
          if (qswnv(i) > 0.) then
            ws(i) = 0.5 * vcmax(i)
          end if


          !> Calculate the maximum electron transport rate Jmax (umol/m2/s)
          !! 1.67 = vcmax25m/jmax25m ratio

          jmax25(i) = 1.67 * vcmax25(i)
          term1(i) = exp(((tmossk(i)/tref) - 1.) * ej/(gasc * tmossk(i)))
          term2(i) = 1. + exp(((tref * sj) - hj)/(tref * gasc))
          term3(i) = 1. + exp(((sj * tmossk(i)) - hj)/(gasc * tmossk(i)))
          jmax(i) = jmax25(i) * term1(i) * term2(i) * term3(i)
           
          ! find the  electron trasport rate in mosses
          if (jmax(i) > 0.0) then
            photon(i) = alpha_moss * parm(i)/sqrt(1.0 + (alpha_moss ** 2 * parm(i) ** 2/(jmax(i) ** 2)))
          else
            photon(i) = 0.0    
          end if

          !> Calculate Wj, Wc (Farquhar and Caemmerer 1982)
          !! wj = light limited, = je in photosynCanopyConduct

          wj(i) = photon(i) * (co2a(i) - gamma(i))/(4. * co2a(i) + (8. * gamma(i)))

          !> Carboxylase(rubisco) limitation = jc in photosynCanopyConduct
          wc(i) = vcmax(i) * (co2a(i) - gamma(i))/(co2a(i) + bc(i))

          !> Choose the minimum of Wj and Wj both having the form:
          !! W = (a Ci - ad) / (e Ci + b)
          !! Then set psna, psnb, psnd and psne for the quadratic solution for net photosynthesis.

          if (wj(i) < wc(i)) then
            psnb(i) = 8. * gamma(i)
            psna(i) = photon(i)
            psne(i) = 4.0
          else if (wc(i) < wj(i)) then
            psnb(i) = bc(i)
            psna(i) = vcmax(i)
            psne(i) = 1.0
          end if

          !> Calculate net and gross photosynthesis by solve the quadratic equation
          !! first root of solution is net photosynthesis An= min(Wj,Wc) - Rd
          !! gross photosynthesis GPP = min(Wc,Wj) = An + Rd

          if (psna(i) > 0.0) then
            mI(i) = rmlmoss(i) - (psnb(i) * g_moss(i)/pres(i)/psne(i)) &
                    - (co2a(i) * g_moss(i)/pres(i)) - (psna(i)/psne(i))

            mII(i) = (psna(i) * co2a(i) * g_moss(i)/pres(i)/psne(i)) - &
                     (rmlmoss(i) * co2a(i) * g_moss(i)/pres(i)) - (rmlmoss(i) * &
                     psnb(i) * g_moss(i)/pres(i)/psne(i)) - (psna(i) * gamma(i) * &
                     g_moss(i)/pres(i)/psne(i))
          else
            mI(i) = 0.0
            mII(i) = 0.0
          end if
          anmoss(i) = ( - mI(i) - (mI(i) * mI(i) - 4 * mII(i)) ** 0.5)/2
          anmoss(i) = min(anmoss(i),ws(i))
        end if !check if the phenology permits photosynthesis
      end if ! ipeatland.
    end do ! loop 200
    return
  end subroutine mossPht
  !! @}
  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatDayEnd
  !! @{
  !> At the end of the day update the moss carbon pool
  !> @author Joe Melton 
  subroutine updateMossC(Cmossmas,nppmosstep,litrfallmoss)
  
  real, intent(in) :: nppmosstep(:) !< moss npp \f$(kg C/m^2/timestep)\f$
  real, intent(in) :: litrfallmoss(:) !< moss litter fall \f$(kg C/m^2/timestep)\f$
  real, intent(inout) :: Cmossmas(:) !< C in moss biomass, \f$kg C/m^2\f$
  
  Cmossmas = Cmossmas + nppmosstep - litrfallmoss
  
end subroutine updateMossC
!! @}
! ---------------------------------------------------------------------------------------------------
  
  !> \ingroup peatlandsmod_peatDayEnd
  !! @{
  !> At the end of the day update the degree days for moss photosynthesis and the peat bottom layer depth
  !> @author Yuanqiao Wu
  subroutine peatDayEnd (nml)

    use ctemStateVars, only : ctem_tile, vgat
    use classStateVars, only : class_gat
    use classicParams, only : ignd, TFREZ

    implicit none

    integer, intent(in) :: nml

    integer :: i, k
    integer :: botlyr

    associate( &
    taaccgat_t       => ctem_tile%taaccgat_t,           & !< real, dimension(:)    : Daily mean air temperature [K] 
    pddgat           => vgat%pdd,                       & !< real, dimension(:)    : peatland degree days above 0 deg C. 
    ipeatlandgat     => vgat%ipeatland,                 & !< integer, dimension(:) : peatland flag, 0 = not peatland, 1 = bog, 2 = fen 
    dlzwgat          => class_gat%dlzwgat,              & !< real, dimension(:,:)  : Permeable thickness of soil layer [m] 
    peatdepgat       => vgat%peatdep,                   & !< real, dimension(:)    : Depth of peat column [m] 
    sdepgat          => vgat%sdepgat,                   & !< real, dimension(:)    : Depth to bedrock in the soil profile [m] 
    zbot             => class_gat%zbot                  & !< real, dimension(:)    : Bottom of soil layers (m) 
    )

    !> Calculate degree days for mossPht Vmax seasonality (only once per day)
    do i = 1,nml
      if (taaccgat_t(i) > tfrez) then
        pddgat(i) = pddgat(i) + taaccgat_t(i) - tfrez
      end if

      !> Update peatland bottom layer depth for fens and bogs 
      if (ipeatlandgat(i) == 1 .or. ipeatlandgat(i) == 2) then
        botlyr = 1
        do k = 1,ignd
          if (peatdepgat(i) < zbot(k)) exit
          botlyr = k
        end do
        
        ! Small error checking, don't let it set dlzwgat to less than the second layer (generally 20 cm).
        ! if 1 is kept it will cause cause an indexing issue below.
        if (botlyr == 1) botlyr = 2
        
        ! The permeable thickness of the soil layer just above the peat is set to the 
        ! difference between the peat depth and the bottom of the layer above.
        dlzwgat(i,botlyr) = peatdepgat(i) - zbot(botlyr-1)  
        
        ! NOTE: The soil permeable depth is just the peat depth as we assume that the 
        ! total soil column in peatlands is peat. If useStaticPeatDep is true then this
        ! will just set sdepgat to itself, since peatdepgat was set to sdepgat up in ctemDriver.
        sdepgat(i) = peatdepgat(i)  
      end if
    end do

    end associate
  end subroutine peatDayEnd

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatDepth
  !! @{
  !> Calculate the peat depth based on equation 18 in Wu, Verseghy, Melton 2016 GMD.
  !! @author Y. Wu, J. Melton
  !!
  real function peatDepth (peatSoilC)

    implicit none

    real, intent(in) :: peatSoilC  !< Peat Soil C mass, \f$kg C/m^2\f$

    ! Calculate the peat depth based on equation 18 in Wu, Verseghy, Melton 2016 GMD.
    ! Since peatlands are run over tiles this uses the tile average soil mass (so only over
    ! the peat tile). 
    !NOTE peatSoilC is, itself, based on peat depth (see peatStorage below), so this is just 
    ! used to update the peat depth as the peat soil C pool changes

    peatDepth = ( -72067.0 + sqrt((72067.0 **2) - (4.0 * 4056.6 &
                * ( -peatSoilC * 1000 / 0.487)))) / (2. * 4056.6)
  end function peatDepth
  !! @}

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatStorage
  !! @{ Finds the carbon storage in the peat based on depth (or oxic and anoxic compartments)
  !! The water table depth delineates the oxic and anoxic compartments.
  !! functions (R**2 = 0.9999) determines the carbon content of each
  !! compartment from a peat bulk density profile based on unpulished
  !! data from P.J.H. Richard (described in fig. 1, Frokling et al.(2001)
  !! conversion of peat into carbon with 48.7% (Mer Bleue unpublished data,
  !! Moore)
  !>
  !! @author Y. Wu, J. Melton
  !!
  real function peatStorage (depth)

    implicit none

    real, intent(in) :: depth  !< Peat compartment depth (m). Either total column or oxic/anoxic.

    peatStorage = (4056.6 * depth**2 + 72067.0 * depth) * 0.487 / 1000.0

  end function peatStorage
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \namespace peatlandsmod
  !! Peatland specific processes 
  !! @author Y. Wu, J. Melton
  !!
  !! The peatland module is published in Geoscientific Model Development (Wu et al. 2016) \cite Wu2016-zt.
  !!
  !! To account for the eco-hydrological and biogeochemical interactions among
  !! vegetation, atmosphere and soil in peatlands, the following modifications
  !! were made to the coupled CLASS3.6--CTEM2.0 modelling framework:
  !!
  !! 1. The top soil layer was characterized as a moss layer with a higher heat
  !! and hydraulic capacity than a mineral soil layer. The moss layer buffers the
  !! exchange of energy and water at the soil surface and regulates the soil
  !! temperature and moisture (Turetsky et al. 2012) \cite Turetsky2012-qh.
  !!
  !! 2. Three peatland vascular PFTs (evergreen shrubs, deciduous shrubs and
  !! sedges) as well as mosses were added to the existing nine CTEM PFTs. These
  !! peatland-specific PFTs are adapted to cold climate and inundated soil with
  !! optimized plant structure (shoot/root ratio, rooting depth), growth strategy
  !! and metabolic acclimations to light, water and temperature.
  !!
  !! 3. We considered the soil inundation stress on microbial respiration in the
  !! litter C pool. The original CTEM assumed that litter respiration was not
  !! affected by oxygen deficit as a result of flooding, since litter was always
  !! assumed to have access to air. This assumption does not hold for peatlands
  !! where high water table positions occur routinely.
  !!
  !! 4. To provide the framework for future runs coupled to the global earth system model, we separated the soil C balance and heterotrophic respiration
  !! (HR) calculations for peatland and non-peatland fractions for each grid cell
  !! in the global model. Over the non-peatland fraction, we use the original CTEM
  !! approach that aggregates the HR from each PFT weighted by the fractional
  !! cover. Over the peatland fraction the soil C pool and decomposition are
  !! controlled by the water table position, following the two-compartment
  !! approach used in the MWM (St. Hilaire et al. 2010) \cite St-Hilaire2010-5e9.
  !!
  !! The standard configuration of soil layers in CLASS consists of three layers with
  !! thickness of 0.10, 0.25, and 3.75, m. Organic soil in CLASS was parameterized
  !! by Letts et al. (2000) \cite Letts2000-pg as fibric, hemic and sapric peat in the three soil
  !! layers respectively, representing fresh, moderately decomposed and highly
  !! decomposed organic matter. Tests of CLASS on peatlands revealed improved
  !! performance in the energy simulations for fens and bogs with this organic
  !! soil parameterization. However, the model overestimated energy and water
  !! fluxes at bog surfaces during dry periods due to the neglect of the moss
  !! cover (Comer et al., 2000) \cite Comer2000-mz.
  !!
  !! To take into account the interaction amongst the moss and the soil layers and
  !! the overlying atmosphere for energy and water transfer, we added a new soil
  !! layer 0.10 m thick above the fibric organic soil to represent living and dead
  !! peatland bryophytes, such as Sphagnum mosses and true mosses
  !! (Bryopsida). The physical characteristics of mosses differ from those of
  !! either the shoots or the roots of vascular plants (Rice et al., 2008). In
  !! particular, mosses can hold more than 30 g of water per gram of biomass
  !! (Robroek et al., 2009). More than 90 \% of the moss leaf volume is
  !! occupied by the water-holding hyaline cells (Rice et al., 2008), which retain
  !! water even when the water table depth declines to 1--10 m below the surface
  !! (Hayward and Clymo, 1982).
  !!
  !! The parameter values of the moss layer for water and energy properties were
  !! derived from a number of recent experiments measuring the hydraulic
  !! properties of mosses (Price et al., 2008 \cite Price2008-fr; Price and Whittington, 2010
  !! McCarter and Price, 2012). Living mosses range from 2--3 to over
  !! 5 cm in height (Rice et al., 2008) and have lower values of dry bulk density
  !! and field capacity than fibric peat (Price et al., 2008) \cite Price2008-fr. Compared to fibric
  !! peat, the saturated hydraulic conductivity of living moss is higher by orders
  !! of magnitude (Price et al., 2008) \cite Price2008-fr and the thermal conductivity is more
  !! affected by the water content (O'Donnell et al., 2009). To fully account for
  !! the effect of mosses, we set the depth of the living moss (\f$z_{m}\f$)
  !! within the top soil (i.e. moss) layer to 3 cm for fens and 4 cm for bogs,
  !! and interpolated its water content \f$w_m\f$ (kg water /(kg dry mass)) from the water content of the overall layer
  !! \f$\theta_{l, 1}\f$ (m\f$^3\f$ water /(m soil)\f$^3\f$) and the depth of the
  !! living moss:
  !!
  !! \f$w_m =\frac{z_m\theta_{l, 1}\rho_w}{B_m} \f$
  !!
  !! where the dry moss biomass (\f$B_m) \f$ is converted from moss C
  !! (C\f$_m \f$) using the standard conversion factor of 0.46 kg C per kg dry
  !! biomass, \f$\theta _{l, 1} \f$ (m\f$^3\f$\, m\f$^{-3}\f$) is the liquid water
  !! content of the top soil layer, and \f$\rho_w \f$ is the density of water
  !! (1000 kg m\f$^{-3} \f$). The maximum and minimum moss water contents were
  !! estimated from a number of observed moss water contents (e.g. Williams and Flanagan, 1998; Robroek et al., 2009). In CLASS, evaporation at the soil
  !! surface is controlled by a soil evaporation efficiency coefficient \f$\beta \f$
  !! (Verseghy, 2012). This parameter is calculated from the liquid water content
  !! and the field capacity of the first soil layer following Lee and
  !! Pielke (1992). For peatlands, \f$\beta \f$ was assumed to be regulated by the
  !! relative moisture of the living moss rather than the ratio of relative liquid
  !! water content of the first soil layer:
  !!
  !! \f$ \beta = 0.25 [ 1- \cos \left( \frac{w_m
  !! -w_{m, min}}{w_m-w_{m, max}} \right)]^{2} \f$
  !!
  !! where \f$w_m \f$, \f$w_{m, max} \f$, and \f$w_{m, min} \f$ are the water content
  !! and the maximum and minimum water contents of the living moss in kg water / (kg dry mass).
  !!
  !!
  !! Moss photosynthesis subroutine
  !!
  !! Mosses are an important contributor to the primary production and the C
  !! sequestration in peatlands, owing to the low decomposability of the moss
  !! tissue. Sphagnum in peatlands grows at
  !! 20--1600 g biomass m\f$^{-2} \f$/yr and accounts for about 50 \% of
  !! the total peat volume (Turetsky, 2003). We have modified CTEM to include a
  !! moss C pool and moss litter pool along with the related C fluxes, i.e.
  !! photosynthesis, autotrophic respiration, heterotrophic respiration, and
  !! humification. The net photosynthesis of moss (\f$G_m)\f$ is calculated
  !! from the gross photosynthesis (\f$G_{0, m})\f$ and dark respiration
  !! (\f$R_{d, m}\f$):
  !!
  !! \f$ G_m=G_{0, m}-R_{d, m} \f$.
  !!
  !! The moss photosynthesis and dark respiration are calculated using the
  !! Farquhar~(1989) biochemical approach following the MWM (St-Hilaire et al.,
  !! 2010) \cite St-Hilaire2010-5e9 and CTEM (Melton and Arora, 2016) \cite Melton2016-zx, with modifications for integration
  !! with CLASS--CTEM and moss phenology. The leaf-level gross photosynthesis rate
  !! \f$G_{0, m} \f$ (\f$\mu \f$ mol CO\f$_2 \f$ /m \f$^2 \f$/s) is obtained as
  !! the minimum of the transportation limited photosynthesis rates (\f$J_s)\f$ and
  !! the first root of the quadratic solution of the light-limited rate (\f$J_e)\f$
  !! and the Rubisco limited rate (\f$J_c)\f$. A logistic factor (\f$\varsigma\f$) is
  !! added with values 0 or 1 to introduce a seasonal control of moss
  !! photosynthesis. In the MWM, spring photosynthesis starts when the snow depth
  !! is below 0.05 m and the soil temperature at 5 cm depth goes above
  !! 0.5 C (Moore et al., 2006). Since in our case CLASS sets the
  !! minimum depth for melting, discontinuous snow to 0.10 m, this limits the
  !! spring photosynthesis to starting only once the snow is completely melted.
  !!
  !! \f$ G_{0, m} = \quad \varsigma \min \left(J_{s}, \frac{(J_{c}+J_{e})\pm
  !! \sqrt{{(J_c +J_ e)}^2 - 4(J_c +J_e)} }{2}\right) \f$
  !!
  !! The dark respiration in mosses (\f$R_{d, m})\f$ is calculated as a function
  !! of the base dark respiration rate (\f$R_{d, m, 0})\f$, which has a value of
  !! 1.1 (\f$\mu \f$ mol CO\f$_2 \f$ /m \f$^2 \f$/s) (Adkinson and Humphreys, 2011) scaled
  !! by the moss moisture (\f$f_{m, rd}\f$) and soil temperature functions
  !! (\f$f_{T, rd})\f$. The moss moisture function is based on the volumetric
  !! water content of the moss, \f$\theta_m\f$ (kg water /(kg dry
  !! mass)). The MWM models the relation between water content in mosses and dark
  !! respiration with optimal water content at 5.8 g water per g dry weight,
  !! following the approach in Frolking et al. (1996). We modified the relation
  !! for water content above the optimal water content, based on a recent
  !! discovery of a weak linear positive relation between the dark respiration
  !! rate and the water content above the optimal water content during the late
  !! summer and fall (Adkinson and Humphreys, 2011):
  !!
  !! \f$R_{d, m} = R_{d, m, 0} f_{m, rd} f_{T, {rd}}\f$
  !!
  !! \f$ f_{T, {rd}} = (3.22-(0.046 \cdot
  !! T_{moss})^{(T_{moss}-25/10)} \f$
  !!
  !! \f$ f_{m, rd} = 0\f$ for\f$ \theta_{m} <0.4 \f$
  !!
  !! \f$ f_{m, rd} = 0.35 \theta_{m}^{2/3}-0.14\f$ for\f$ 0.4\le \theta_m<5.8 \f$
  !!
  !! \f$ f_{m, rd} = 0.01 \theta_{m} +0.942\f$ for \f$5.8<\theta _{m} \f$
  !!
  !! Photosynthetic photon flux density (PPFD) is measured by the
  !! photosynthetically active radiation (PAR), which is defined as the solar
  !! radiation between 0.4 to 0.7 \f$\mu \f$ mol that can be used by plants via
  !! photosynthesis. In the coupled CLASS--CTEM system, the PAR received by the
  !! moss (PAR\f$_m\f$, unit \f$\mu \f$mol photons m\f$^-2\f$/s is
  !! converted from the visible short-wave radiation reaching the ground (\f$K_{\ast g}\f$,
  !! unit W/m\f$^2\f$) in CLASS by a factor of
  !! 4.6 \f$\mu\f$ mol /m\f$^2\f$/s per W/m\f$^2\f$ (McCree, 1972).
  !! \f$K_{\ast g}\f$ is a function of the incoming short-wave radiation
  !! (\f$K\downarrow \f$, unit: W/m\f$^2\f$), the surface albedo (\f$\alpha_g\f$),
  !! and the canopy transmissivity (\f$\tau _c\f$):
  !!
  !! \f$K_{\ast g} = K \downarrow \tau_{c}\left( 1-\alpha_{g}\right)\f$
  !!
  !! The energy uptake by the moss layer is thus a function of the total incoming
  !! short-wave radiation, the aggregated LAI of the PFTs present, the snow depth,
  !! the fractional vegetation cover, and the soil water content (Verseghy, 2012).
  !! In peatland C models that do not consider vegetation dynamics, the
  !! transmissivity of the vegetation canopy is usually assumed to be constant
  !! (e.g. St-Hilaire et al., 2010) \cite St-Hilaire2010-5e9. Compared with such models, CLASS enables a
  !! more detailed representation of light incident on the moss surface since it
  !! includes partitioning of direct/diffuse and visible/near-IR radiation,
  !! PFT-specific transmissivities, and time-varying LAI and fractional PFT
  !! coverages (Verseghy, 2012) \cite Verseghy2012-c0e.
  !!
  !> \file
end module peatlandsMod
