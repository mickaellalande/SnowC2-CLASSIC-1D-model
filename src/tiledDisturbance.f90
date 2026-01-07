!> \file
!> Module paired with dynamic tiling to simulate disturbance
module tiledDisturbance

  ! S.R. Curasi May, 2022

  use classicParams,  only : nlat, nmos, ilg, ican, ignd, icp1, icc, iccp2, iccp1, convertg2kg, dynTilingTolrance, paper, furniture, combust

  implicit none

  ! Subroutines contained in this module:
    public  :: harvestTile
    public  :: tiledDisturbancePrep
    public  :: initialMassBalance
    public  :: finalMassBalance

contains

  ! ------------------------------------------------------------------

  !> \ingroup tiledDisturbance_harvestTile
  !! @{
  !> Tile harvesting subroutine which removes material from the tile and transfers it to the LUC pools
  !> @author S.R. Curasi
  subroutine harvestTile(Ncycle_on, il1, il2, timharvareagat, timharvarearow, tileAgegat, &
                  gleafmas, gleafmas_ns, gleafmas_s, &
                  bleafmas, stemmass, stemmass_ns, stemmass_s, &
                  rootmass, rootmass_ns, rootmass_s, &
                  litrmass, soilcmas, vgbiomas, gavgltms, &
                  gavgscms, fcancmx, ngleafmas, ngleafmas_ns, &
                  ngleafmas_s, nbleafmas, nstemmass, nstemmass_ns, &
                  nstemmass_s, nrootmass, nrootmass_ns, nrootmass_s, &
                  nlitrmass, soilnmas, &
                  lucemcom, lucltrin, lucsocin, lucemcomn, lucltrinn, lucsocinn, faregat, rmatctem, &
                  trackTileAge, dynamicTilingOn, lfstatus, pandays, flhrloss, flhrloss_ns, flhrloss_s) 

    use classicParams,    only : icc, ican, iccp1, iccp2, ignd, icp1, ilg, deltat, reindexPFTs, classpfts, bmasthrs

    implicit none

    integer, intent(in) :: il1   !< il1=1
    integer, intent(in) :: il2   !< il2=ilg

    logical, intent(in) :: Ncycle_on
    logical, intent(in) :: trackTileAge
    logical, intent(in) :: dynamicTilingOn

    real, intent(in) :: faregat(ilg)
    real, intent(in), dimension(ilg,icc) :: fcancmx          !< max. fractional coverages of ctem's 9 pfts.

    real, intent(in), dimension(ilg) :: timharvareagat       !< the same as timharvarea, but in CTEM's 'gat' format 
    real, intent(in), dimension(nlat,nmos) :: timharvarearow !< the fractional area of the cell where timber will be harvested
                                                             !! (0-1,this is the actual area harvested per tile by the model, which is determined 
                                                             !! from timharvrow by the harvest subroutines)
    real, intent(inout), dimension(ilg) :: tileAgegat        !< the age of the tile since the start of the run in months
                                                             !! this is reset by harvest and fire and incremented at the CTEM timestep by indexTileAge
    real, intent(in), dimension(ilg,icc,ignd) :: rmatctem    !< fraction of live roots in each soil layer for each of ctem's 9 pfts

    real, intent(inout), dimension(ilg,icc) :: gleafmas      !< green or live leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: gleafmas_ns   !< non-structural green or live leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: gleafmas_s    !< structural green or live leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: bleafmas      !< brown or dead leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: stemmass      !< stem biomass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: stemmass_ns   !< non-structural stem biomass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: stemmass_s    !< structurla stem biomass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: rootmass      !< root biomass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: rootmass_ns   !< non-structural root biomass in kg c/m2, for the 9 pfts
    real, intent(inout), dimension(ilg,icc) :: rootmass_s    !< structural root biomass in kg c/m2, for the 9 pfts 
    real, intent(inout), dimension(ilg,icc) :: ngleafmas     !< green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: ngleafmas_ns  !< non-structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: ngleafmas_s   !< structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: nbleafmas     !< brown leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: nstemmass     !< stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: nstemmass_ns  !< non-structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: nstemmass_s   !< structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: nrootmass     !< root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: nrootmass_ns  !< non-structural root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,icc) :: nrootmass_s   !< structural root nitrogen mass for individual PFTs (\f$g N/m^2\f)
    real, intent(inout), dimension(ilg) :: vgbiomas          !< grid averaged vegetation biomass, kg c/m2
    real, intent(inout), dimension(ilg,iccp2) :: soilnmas    !< soil organic nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,iccp2) :: nlitrmass   !< litter nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout), dimension(ilg,iccp2,ignd) :: soilcmas !< soil c mass in kg c/m2, for the 9 pfts + bare
    real, intent(inout), dimension(ilg,iccp2,ignd) :: litrmass !< litter mass in kg c/m2, for the 9 pfts + bare
    integer, intent(inout), dimension(ilg,icc) :: lfstatus 
    integer, intent(inout), dimension(ilg,icc) :: pandays    !< counter for positive net photosynthesis (an) days for initiating leaf onset
    real, intent(inout), dimension(ilg,icc) :: flhrloss      !< fall & harvest loss for bdl dcd plants and crops, respectively, \f$kg c/m^2\f$.
    real, intent(inout), dimension(ilg,icc) :: flhrloss_ns   !< non-structural fall & harvest loss for bdl dcd plants and crops, respectively, \f$kg c/m^2\f$.
    real, intent(inout), dimension(ilg,icc) :: flhrloss_s    !< structural fall & harvest loss for bdl dcd plants and crops, respectively, \f$kg c/m^2\f$.

    real, intent(inout), dimension(ilg) :: gavgltms       !< grid averaged litter mass including the LUC product pool, kg c/m2
    real, intent(inout), dimension(ilg) :: gavgscms       !< grid averaged soil c mass including the LUC product pool, kg c/m2

    real, intent(out), dimension(ilg) :: lucemcom !< luc related carbon emission losses from combustion u-mol co2/m2.sec
    real, intent(out), dimension(ilg) :: lucltrin !< luc related input to litter pool, u-mol co2/m2.sec
    real, intent(out), dimension(ilg) :: lucsocin !< luc related input to soil carbon pool, u-mol co2/m2.sec
    real, intent(out), dimension(ilg) :: lucemcomn !< luc related carbon emission losses from combustion of N (g N m^-2 day^-1)
    real, intent(out), dimension(ilg) :: lucltrinn !< luc related input to litter pool of N (g N m^-2 day^-1)
    real, intent(out), dimension(ilg) :: lucsocinn !< luc related input to soil nitrogen pool of N (g N m^-2 day^-1)

    !mass balance check variables
    real, dimension(ilg) :: vegCMass_initial !< temporary storage variable for the total mass of carbon in vegetation
    real, dimension(ilg) :: soilCMass_initial !< temporary storage variable for the total mass of carbon in soil and litter
    real, dimension(ilg) :: peatCMass_initial !< temporary storage variable for the total mass of carbon in the peatland pools
    real, dimension(ilg) :: lucMass_initial !< temporary storage variable for the total mass of carbon in the LUC pools
    real, dimension(ilg) :: totalCMass_initial !< temporary storage variable for the total mass of carbon in all pools

    real, dimension(ilg) :: vegNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, dimension(ilg) :: soilNMass_initial !< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, dimension(ilg) :: lucNmass_initial !< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, dimension(ilg) :: totalNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    integer, dimension(ilg,icc) :: treatind !< treatment index for combust, paper, & furniture
    real, dimension(ilg,icc) :: abvgmass   !< above-ground biomass (kg c/m^2) used to detemine the treatment thats applied for timber harvest

    integer :: i, j, k, m

    real, dimension(ilg) :: barefrac = 1. !< initialized to 1.0
    real :: abvpooltots        !< temp pool used for aboveground biomass, \f$kg c/m^2\f$.
    real :: Nabvpooltots       !< temp pool used for aboveground N mass, (\f$g N/m^2\f$).

    if (any(timharvareagat > 0.0)) then
        call initialMassBalance(fcancmx, gleafmas, bleafmas, stemmass, rootmass, &
                                soilcmas, litrmass, nlitrmass, soilnmas, &
                                ngleafmas, nbleafmas ,nstemmass, nrootmass, & 
                                vegCMass_initial, soilCMass_initial, lucMass_initial, totalCMass_initial, &
                                vegNMass_initial, soilNMass_initial, lucNMass_initial, totalNMass_initial, &
                                Ncycle_on, faregat)

        abvgmass(:,:) = 0.0                        

        ! Find the barefrac
        do i = il1,il2
          do j = 1,icc
            barefrac(i) = barefrac(i) - fcancmx(i,j)
          end do

          !< get the treatment indexes, which are the same values used in land use change (i.e. Arora & Boer 2010 table 1)
          !< they're used to fraction biomass into the combust, paper and wood products pools based upon PFT and biomass thresholds
          do j = 1,ican
            do m = reindexPFTs(j,1),reindexPFTs(j,2)
              abvgmass(i,m) = gleafmas(i,m) + bleafmas(i,m) + stemmass(i,m)
              select case (classpfts(j))
              case ('Crops', 'Grass') ! Crops and Grass
                treatind(i,m) = 3
              case ('NdlTr','BdlTr')
                if (abvgmass(i,m) >= bmasthrs(1)) then ! forest
                  treatind(i,m) = 1
                else if (abvgmass(i,m) <= bmasthrs(2)) then ! bush
                  treatind(i,m) = 3
                else  ! shrubland
                  treatind(i,m) = 2
                end if
              case ('BdlSh')
                if (abvgmass(i,m) >= bmasthrs(1)) then ! forest
                  treatind(i,m) = 1
                else if (abvgmass(i,m) <= bmasthrs(2)) then ! bush
                  treatind(i,m) = 3
                else  ! shrubland
                  treatind(i,m) = 2
                end if
              case default
                 print * ,'Unknown CLASS PFT in LUC ',classpfts(j)
                call errorHandler('harvestTile', - 1)
              end select
            end do
          end do
        end do

        lucemcom(:) = 0.0
        lucltrin(:) = 0.0
        lucsocin(:) = 0.0
        lucemcomn(:) = 0.0
        lucltrinn(:) = 0.0
        lucsocinn(:) = 0.0

        do i = il1,il2
          if (timharvareagat(i) > 0.0) then !if there is going to be harvest right now
            do j = 1,icc
              
                !< change the leaf status and reset phenological variables 
                !< in line with the vegetation having no leaves
                if(dynamicTilingOn) then
                  lfstatus(i,j) = 0
                  pandays(i,j) = 0
                  flhrloss(i,j) = 0.0
                  flhrloss_ns(i,j) = 0.0
                  flhrloss_s(i,j)  = 0.0
                end if

                !< distribute above ground biomass to the appropriate luc pools
                abvpooltots = gleafmas_ns(i,j) + gleafmas_s(i,j) + bleafmas(i,j) + stemmass_ns(i,j) + stemmass_s(i,j)
                litrmass(i,iccp2,1) = litrmass(i,iccp2,1) + (abvpooltots * timharvareagat(i) * fcancmx(i,j) * paper(treatind(i,j)))
                soilcmas(i,iccp2,1) = soilcmas(i,iccp2,1) + (abvpooltots * timharvareagat(i) * fcancmx(i,j) * furniture(treatind(i,j)))
                
                !< and fluxes
                lucemcom(i) = lucemcom(i) + (abvpooltots * timharvareagat(i) * fcancmx(i,j) * combust(treatind(i,j))* 963.62)
                lucltrin(i) = lucltrin(i) + (abvpooltots * timharvareagat(i) * fcancmx(i,j) * paper(treatind(i,j))* 963.62)
                lucsocin(i) = lucsocin(i) + (abvpooltots * timharvareagat(i) * fcancmx(i,j) * furniture(treatind(i,j))* 963.62)

                !<distrbute roots as litter into the soil profile
                !! Both roots and litter a tracked on a per-pft basis unlike the aggregated land use change pools or fluxes 
                !! so we transfer the mass directly without multiplying by fcancmx
                do k = 1,ignd
                  litrmass(i,j,k) = litrmass(i,j,k) + ((rootmass_ns(i,j) + rootmass_s(i,j)) * timharvareagat(i) * rmatctem(i,j,k))
                end do 

                !<finalize the re-distribution of the material
                gleafmas_ns(i,j) = gleafmas_ns(i,j) * (1.0 - timharvareagat(i))
                gleafmas_s(i,j) = gleafmas_s(i,j) * (1.0 - timharvareagat(i))
                bleafmas(i,j) = bleafmas(i,j) * (1.0 - timharvareagat(i))
                stemmass_ns(i,j) = stemmass_ns(i,j) * (1.0 - timharvareagat(i))
                stemmass_s(i,j) = stemmass_s(i,j) * (1.0 - timharvareagat(i))
                rootmass_ns(i,j) = rootmass_ns(i,j) * (1.0 - timharvareagat(i))
                rootmass_s(i,j) = rootmass_s(i,j) * (1.0 - timharvareagat(i))

                if (Ncycle_on) then

                  !< distribute above ground nitrogen to the appropriate luc pools
                  Nabvpooltots = ngleafmas_ns(i,j) + ngleafmas_s(i,j) + nbleafmas(i,j) + nstemmass_ns(i,j) + nstemmass_s(i,j)
                  nlitrmass(i,iccp2) = nlitrmass(i,iccp2) + (Nabvpooltots * timharvareagat(i) * fcancmx(i,j) * paper(treatind(i,j)))
                  soilnmas(i,iccp2) = soilnmas(i,iccp2) + (Nabvpooltots * timharvareagat(i) * fcancmx(i,j) * furniture(treatind(i,j)))
                  
                  !< and fluxes
                  lucemcomn(i) = lucemcomn(i) + (Nabvpooltots * timharvareagat(i) * fcancmx(i,j) * combust(treatind(i,j)))
                  lucltrinn(i) = lucltrinn(i) + (Nabvpooltots * timharvareagat(i) * fcancmx(i,j) * paper(treatind(i,j)))
                  lucsocinn(i) = lucsocinn(i) + (Nabvpooltots * timharvareagat(i) * fcancmx(i,j) * furniture(treatind(i,j)))

                  !<distrbute root N into the correct N litter pools
                  !! Both roots and litter a tracked on a per-pft basis unlike the aggregated land use change pools or fluxes 
                  !! so we transfer the mass directly without multiplying by fcancmx
                  nlitrmass(i,j) = nlitrmass(i,j) + ((nrootmass_ns(i,j) + nrootmass_s(i,j)) * timharvareagat(i))
   
                  !<finalize the re-distribution of the material
                  ngleafmas_ns(i,j) = ngleafmas_ns(i,j) * (1.0 - timharvareagat(i))
                  ngleafmas_s(i,j) = ngleafmas_s(i,j) * (1.0 - timharvareagat(i))
                  nbleafmas(i,j) = nbleafmas(i,j) * (1.0 - timharvareagat(i))
                  nrootmass_ns(i,j) = nrootmass_ns(i,j) * (1.0 - timharvareagat(i))
                  nrootmass_s(i,j) = nrootmass_s(i,j) * (1.0 - timharvareagat(i))
                  nstemmass_ns(i,j) = nstemmass_ns(i,j) * (1.0 - timharvareagat(i))
                  nstemmass_s(i,j) = nstemmass_s(i,j) * (1.0 - timharvareagat(i))
                end if
            end do

            !> update the cumulative vegetation mass pools
            gleafmas(i,:) = gleafmas_ns(i,:) + gleafmas_s(i,:) 
            stemmass(i,:) = stemmass_ns(i,:) + stemmass_s(i,:) 
            rootmass(i,:) = rootmass_ns(i,:) + rootmass_s(i,:) 

            !> and n pools
            if (Ncycle_on) then
              ngleafmas(i,:) = ngleafmas_ns(i,:) + ngleafmas_s(i,:)
              nstemmass(i,:) = nstemmass_ns(i,:) + nstemmass_s(i,:)
              nrootmass(i,:) = nrootmass_ns(i,:) + nrootmass_s(i,:)
            end if

            !> update grid averaged vegetation biomass,and litter and soil c densities
            vgbiomas(i) = 0.0
            gavgltms(i) = 0.0
            gavgscms(i) = 0.0

            do j = 1,icc 
              vgbiomas(i) = vgbiomas(i) + fcancmx(i,j) * (gleafmas(i,j) + &
                            bleafmas(i,j) + stemmass(i,j) + rootmass(i,j))
              do k = 1,ignd
                gavgltms(i)=gavgltms(i) + fcancmx(i,j) * litrmass(i,j,k)
                gavgscms(i)=gavgscms(i) + fcancmx(i,j) * soilcmas(i,j,k)
              end do
            end do
            do k = 1,ignd
                gavgltms(i)=gavgltms(i) + ( barefrac(i) * litrmass(i,iccp1,k) )
                gavgscms(i)=gavgscms(i) + ( barefrac(i) * soilcmas(i,iccp1,k) )
            end do
          end if
        end do

        call finalMassBalance(fcancmx, gleafmas, bleafmas, stemmass, rootmass, &
                                soilcmas, litrmass, lucemcom, nlitrmass, soilnmas, &
                                ngleafmas, nbleafmas ,nstemmass, nrootmass, lucemcomn, & 
                                vegCMass_initial, soilCMass_initial, lucMass_initial, totalCMass_initial, &
                                vegNMass_initial, soilNMass_initial, lucNMass_initial, totalNMass_initial, &
                                Ncycle_on, faregat)
                                    
        !< alter the tile age based upon the harvest area
        ! timharvareagat ranges from zero to one, with one representing an entire tile being harvested (i.e. when dynamic tiling on and there is a dedicate tile
        ! to represent the harvested fraction of the gridcell) and values less than one representing some fraction of a tile being harvested (i.e. when dynamic 
        ! tiling is off and harvest imapcts the average indvidual). This equation calculates the new age of tile based upon a weighted average of the ages of the 
        ! harvested and non-harvested areas. In it's non-simplified form new_age=((age*(1-timharvareagat))+(0*timharvareagat))/((1-timharvareagat)+timharvareagat)
        ! by the above definition ((1-timharvareagat)+timharvareagat)=1
        if(trackTileAge .or. dynamicTilingOn) tileAgegat(:) = tileAgegat(:)*(1-timharvareagat(:))
    end if

 

  end subroutine harvestTile
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup tiledDisturbance_tiledDisturbancePrep
  !! @{
  !> Subroutine that prepares information for the tiled disturbance subroutines 
  subroutine tiledDisturbancePrep(ctem_on, lnduseon, PFTCompetition, timberHarvest, prescribedFire, &
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
    
    use, intrinsic :: iso_fortran_env, only: r8=>real64
    use modelStateDrivers, only : updateInput
    use dynamicTiling, only : moveSplitCopyTiles, cleanTiles, indexSort
    use classicParams, only : DELT, crop, dynTilingMaxDisturb

    implicit none

    !flags and some control values from main
    logical, intent(in) :: ctem_on
    logical, intent(in) :: lnduseon
    logical, intent(in) :: PFTCompetition
    logical, intent(in) :: timberHarvest
    logical, intent(in) :: dynamicTilingOn
    logical, intent(in) :: prescribedFire
    integer, intent(in) :: runyr         !< Year of the model run (counts up starting with readMetStartYear continously, even if metLoop > 1)
    integer, intent(in) :: IDAY       !< Julian day of the year (computed in updateMet, which is not done for CanESM, so must be passed in).
    integer, intent(in) :: NDAY          !< Number of short (physics) timesteps in one day. e.g., if physics timestep is 15 min this is 48.
    integer, intent(in) :: NCOUNT        !< Counter for daily averaging

    real, intent(inout), dimension(nlat) :: prsfirerow                !< the annual fractional area where fire burned (0-1, read in from the external forcing file)
    real, intent(inout), dimension(nlat,nmos) :: prsfirearearow       !< the fractional area of the cell where fire burned 
                                                                      !!(0-1,this is the actual area burned per tile by the model, which is determined from prsfirerow by the
                                                                      !! prescribed fire subroutines)
    real, intent(inout), dimension(ilg) :: prsfireareagat             !< the same as prsfirearea, but in CTEM's 'gat' format 
    real, intent(inout), dimension(nlat) :: timharvrow                !< the annual fractional area where timber was harvested (0-1,read in from the external forcing file)
    real, intent(inout), dimension(nlat,nmos) :: timharvarearow       !< the fractional area of the cell where timber will be harvested
                                                                      !! (0-1,this is the actual area harvested per tile by the model, which is determined 
                                                                      !! from timharvrow by the harvest subroutines)
    real, intent(inout), dimension(ilg) :: timharvareagat             !< the same as timharvarea, but in CTEM's 'gat' format 
    real, intent(inout), dimension(ilg) :: tileAgegat                 !< the age of the tile since the start of the run in months 
                                                                      !! this is reset by harvest and fire and incremented at the CTEM timestep by indexTileAge

    !setup the intputs for class state vars
    real, intent(inout), dimension(nlat,nmos,ignd) :: DLZWROT !< Permeable thickness of soil layer [m]
    real, intent(in), dimension(ignd) :: DELZ    !< Overall thickness of soil layer [m]
    real, intent(inout), dimension(nlat,nmos) :: ALBSROT !< Snow albedo [ ]
    real, intent(inout), dimension(nlat,nmos) :: CMAIROT !<
    real, intent(inout), dimension(nlat,nmos) :: GROROT  !< Vegetation growth index [ ]
    real, intent(inout), dimension(nlat,nmos) :: QACROT  !<
    real, intent(inout), dimension(nlat,nmos) :: RCANROT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(inout), dimension(nlat,nmos) :: RHOSROT !< Density of snow \f$[kg m^{-3}]\f$
    real, intent(inout), dimension(nlat,nmos) :: SCANROT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(inout), dimension(nlat,nmos) :: SNOROT  !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, intent(inout), dimension(nlat,nmos) :: TACROT  !<
    real, intent(inout), dimension(nlat,nmos) :: TBASROT !<
    real, intent(inout), dimension(nlat,nmos) :: TCANROT !< Vegetation canopy temperature [K]
    real, intent(inout), dimension(nlat,nmos) :: TPNDROT !< Temperature of ponded water [K]
    real, intent(inout), dimension(nlat,nmos) :: TSNOROT !< Snowpack temperature [K]
    real, intent(inout), dimension(nlat,nmos) :: WSNOROT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, intent(inout), dimension(nlat,nmos) :: ZPNDROT !< Depth of ponded water [m]
    real, intent(inout), dimension(nlat,nmos) :: DRNROT  !<
    real, intent(inout), dimension(nlat,nmos) :: FAREROT !< Fractional coverage of mosaic tile on modelled area
    real, intent(inout), dimension(nlat,nmos) :: ZSNLROT !<
    real, intent(inout), dimension(nlat,nmos) :: SDEPROT !< Depth to bedrock in the soil profile
    real, intent(inout), dimension(nlat,nmos) :: SOCIROT !<
    real(r8), intent(inout), dimension(nlat,nmos,ignd) :: TBARROT !< Temperature of soil layers [K]
    real, intent(inout), dimension(nlat,nmos,ignd) :: THICROT !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, intent(inout), dimension(nlat,nmos,ignd) :: THLQROT !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, intent(inout), dimension(nlat,nmos,ignd) :: SANDROT !< Percentage sand content of soil
    real, intent(inout), dimension(nlat,nmos,ignd) :: CLAYROT !< Percentage clay content of soil
    real, intent(inout), dimension(nlat,nmos,ignd) :: ORGMROT !< Percentage organic matter content of soil
    real, intent(inout), dimension(nlat,nmos,ican) :: CMASROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: PAMNROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: PAMXROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: PSGAROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: PSGBROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: QA50ROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: ROOTROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: RSMNROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: VPDAROT !<
    real, intent(inout), dimension(nlat,nmos,ican) :: VPDBROT !<
    real, intent(inout), dimension(nlat,nmos,icp1) :: ALICROT !<
    real, intent(inout), dimension(nlat,nmos,icp1) :: ALVCROT !<
    real, intent(inout), dimension(nlat,nmos,icp1) :: FCANROT !<
    real, intent(inout), dimension(nlat,nmos,icp1) :: LNZ0ROT !<
    integer, intent(inout), dimension(nlat,nmos,6,50) :: ITCTROT !<
    real, intent(inout), dimension(nlat,nmos,4)  :: TSFSROT !<
    real, intent(inout), dimension(nlat,nmos) :: maxAnnualActLyrROT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, intent(inout), dimension(nlat,nmos,ignd) :: THPROT
    real, intent(inout), dimension(nlat,nmos,ignd) :: HCPSROT

    !setup the intputs for ctem state vars
    real, intent(inout), dimension(nlat,nmos) :: litrmsmossrow
    real, intent(inout), dimension(nlat,nmos) :: Cmossmasrow
    real, intent(inout), dimension(nlat,nmos) :: dmossrow
    real, intent(inout), dimension(nlat,nmos) :: peatSoilCrow         !< peat soil C mass, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,8)  :: slopefracrow          !< prescribed fraction of wetlands based on slope
    integer, intent(inout), dimension(nlat,nmos,icc) :: lfstatusrow
    integer, intent(inout), dimension(nlat,nmos,icc) :: pandaysrow
    real, intent(inout), dimension(nlat,nmos,icc) :: gleafmassrow
    real, intent(inout), dimension(nlat,nmos,icc) :: gleafmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: gleafmas_NSrow
    real, intent(inout), dimension(nlat,nmos,icc) :: bleafmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: stemmassrow
    real, intent(inout), dimension(nlat,nmos,icc) :: stemmasssrow
    real, intent(inout), dimension(nlat,nmos,icc) :: stemmass_NSrow
    real, intent(inout), dimension(nlat,nmos,icc) :: rootmassrow
    real, intent(inout), dimension(nlat,nmos,icc) :: rootmasssrow
    real, intent(inout), dimension(nlat,nmos,icc) :: rootmass_NSrow
    real, intent(inout), dimension(nlat,nmos,icc) :: fcancmxrow
    real, intent(inout), dimension(nlat,nmos,icc) :: ngleafmassrow
    real, intent(inout), dimension(nlat,nmos,icc) :: ngleafmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: ngleafmas_NSrow
    real, intent(inout), dimension(nlat,nmos,icc) :: nbleafmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: nstemmassrow
    real, intent(inout), dimension(nlat,nmos,icc) :: nstemmasssrow
    real, intent(inout), dimension(nlat,nmos,icc) :: nstemmass_NSrow
    real, intent(inout), dimension(nlat,nmos,icc) :: nrootmassrow
    real, intent(inout), dimension(nlat,nmos,icc) :: nrootmasssrow
    real, intent(inout), dimension(nlat,nmos,icc) :: nrootmass_NSrow
    real, intent(inout), dimension(nlat,nmos,icc) :: co2i1cgrow
    real, intent(inout), dimension(nlat,nmos,icc) :: co2i1csrow
    real, intent(inout), dimension(nlat,nmos,icc) :: co2i2cgrow
    real, intent(inout), dimension(nlat,nmos,icc) :: co2i2csrow
    real, intent(inout), dimension(nlat,nmos,icc) :: flhrlossrow
    real, intent(inout), dimension(nlat,nmos,icc) :: flhrloss_nsrow
    real, intent(inout), dimension(nlat,nmos,icc) :: flhrloss_srow
    real, intent(inout), dimension(nlat,nmos,icc) :: grwtheffrow
    real, intent(inout), dimension(nlat,nmos,icc) :: lystmmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: lyrotmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: tymaxlairow
    real, intent(inout), dimension(nlat,nmos,icc) :: stmhrlosrow
    real, intent(inout), dimension(nlat,nmos,icc) :: rothrlosrow
    real, intent(inout), dimension(nlat,nmos,icc) :: bnfantrow
    real, intent(inout), dimension(nlat,nmos,icc) :: bnfnatrow
    real, intent(inout), dimension(nlat,nmos) :: soilpHrow
    real, intent(inout), dimension(nlat,nmos) :: cfluxcgrow
    real, intent(inout), dimension(nlat,nmos) :: cfluxcsrow
    integer, intent(inout), dimension(nlat,nmos) :: ipeatlandrow
    integer, intent(inout), dimension(nlat,nmos) :: colddays_leaffallrow
    integer, intent(inout), dimension(nlat,nmos) :: colddays_harvestrow
    real, intent(inout), dimension(nlat,nmos,iccp2,ignd) :: litrmassrow
    real, intent(inout), dimension(nlat,nmos,iccp2,ignd) :: soilcmasrow
    real, intent(inout), dimension(nlat,nmos,iccp1) :: nh4_massrow
    real, intent(inout), dimension(nlat,nmos,iccp1) :: no3_massrow
    real, intent(inout), dimension(nlat,nmos,iccp2) :: nlitrmassrow
    real, intent(inout), dimension(nlat,nmos,iccp2) :: soilnmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: lygleafmasmaxrow
    real, intent(inout), dimension(nlat,nmos,icc) :: lystemmassmaxrow
    real, intent(inout), dimension(nlat,nmos,icc) :: lyrootmassmaxrow
    real, intent(inout), dimension(nlat,nmos) :: tileAgerow            !< the age of the tile since the start of the run in months
    real, dimension(nlat,nmos,icc) :: veghghtrow
    integer, dimension(nmos) :: ageSort            !< the tile indicies sorted by age
    integer, dimension(nmos) :: btermSort  !< the tile indicies sorted by the bterm value
    real, dimension(nlat,nmos) :: bterm_veg !< the vegetation biomass fire term used to rank tiles to be burned

    real :: areaTemp            !< the age of the tile since the start of the run in months

    !declare the tracer pools variables
    real, intent(inout), dimension(nlat,nmos) :: tracerMossCMassrot     !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos) :: tracerMossLitrMassrot   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerGLeafMassrot      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerBLeafMassrot      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerStemMassrot       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerRootMassrot       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(inout), dimension(nlat,nmos,iccp2,ignd) :: tracerLitrMassrot       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout), dimension(nlat,nmos,iccp2,ignd) :: tracerSoilCMassrot      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    !these are the switched
    logical, intent(in) :: Ncycle_on
    integer, intent(in) :: useTracer

    !declare input variables that control the tile copying subroutine
    integer, intent(in), dimension(nlat,nmos)     :: controlVector !< the control vector specified in the initializtion file
    logical, intent(inout)  :: DynTilInitializeFlag !< a flag used to re-initialize the model after dynamic tiling runs
    integer, dimension(nlat,nmos)                :: inputVector !< the input vector (1 = operate on tile, 0 = do not operate on tile) 
    character(len=5), dimension(nlat)        :: mode !< a character string that specifies the operation (i.e. move, split, copy, or skip) 
    real, dimension(nlat)                    :: FAREAROUT !< a real value greater than zero and less than one which specifies the desired size of the output tile
    CHARACTER(LEN=24) :: title
    integer, dimension(nlat) :: outputIndex !< the index of the tile where the output of movecopy tiles will be written
    real, dimension(nlat) :: tneeded !< the number of tiles needed to continue running the model
    real, dimension(nlat) :: thave !< the number of tiles we have to continue running the model
    real, dimension(nlat) :: fcancmxcrop !< the fractional area of crops for use in calculating area burned

    integer :: i, j, m
  
    !perform a check and stop the model if CTEM is off, competition is on or agricultural LUC is on
    if ((.not.(ctem_on)) .or. lnduseon .or. PFTCompetition) then
        write(6, * ) 'tiledDisturbancePrep: tiled Disturbance not avaliable with competition, agricultural LUC or when CTEM is off'
        call errorHandler('tiledDisturbancePrep', -1)
    end if

    !set FAREAROUT to 0.0. moveSplitCopyTiles will error out with a value of 0.0
    !so the value must be re-set based on the file below
    FAREAROUT = 0.0

    !if timber harvest is on and this is the first day of the year execute the timber harvest case
    if ((iday == 1) .and. (ncount == nday) .and. (timberHarvest .or. prescribedFire)) then

      !print *,"timharvrow = ", timharvrow
      !print *,"prsfirerow = ", prsfirerow

      !decide if harvest need to occur and select the appropriate case, probably the simplified dyn cases are not needed
      if (any(timharvrow(:) > 0.0) .and. dynamicTilingOn .and. timberHarvest .and. .not.(prescribedFire)) then
        title = 'timhar'  ! Timber harvest
      else if (any(prsfirerow(:) > 0.0) .and. dynamicTilingOn .and. prescribedFire .and. .not.(timberHarvest) ) then
        title = 'prefire'  ! Prescribed fire
      else if ((any(timharvrow(:) > 0.0) .or. any(prsfirerow(:) > 0.0)) .and. dynamicTilingOn .and. prescribedFire .and. timberHarvest) then
        title = 'firehar'  ! Both timber harvest and prescribed fire
      else if ((any(prsfirerow(:) > 0.0) .or. any(timharvrow(:) > 0.0)) .and. .not.(dynamicTilingOn)) then
        title = 'notdyn'  ! Prescribed fire or timber harvest but not with dynamic tiles
      else if (any(timharvrow(:) < 0.0) .or. any(prsfirerow(:) < 0.0)) then
        write(6, * ) 'tiledDisturbancePrep: timber harvest or burned area is less than zero'
        call errorHandler('tiledDisturbancePrep', -2)
      else if (dynamicTilingOn .and. (timberHarvest .or. prescribedFire)) then
       title = 'cleandyn' ! run cleanTiles to pre-emptively join tiles in years with no harvest or fire prescribed
      else 
        title = 'noaction' ! The defeault return, nothing needs to happen
      end if 
    else if (ncount == nday .and. (timberHarvest .or. prescribedFire)) then
       title = 'reset' ! On the second day of the year leave the tiles as is but make sure that timharvarearow is set to zero so extra biomass isn't harvested
    else
      title = 'noaction' ! The defeault return, nothing needs to happen
    end if

    select case (title)

    !The case structure here is intended to make things easier later on
    !and make easier to create special cases

    case ('timhar') ! Timber harvest -----------------------

        ! Find the number of tiles needed
        tneeded(:) = 0.0
        do i=1,nlat 
          if((timharvrow(i) > 0.0) .and. (.not. any((controlVector(i,:) /= 0 .and. FAREROT(i,:) == 0.0)))) then
            tneeded(i) = 1
          end if
        end do

        ! call routine which simplifies tiles based on a vegetation height threshold and/or when space is needed. 
        call cleanTiles(controlVector,tneeded,DynTilInitializeFlag,outputIndex, & !these are the main intputs that control the routine
                  ALBSROT, ALICROT, ALVCROT, CLAYROT, CMASROT, Cmossmasrow, DRNROT, FAREROT, FCANROT, & !these are just other major variables
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
                  useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow, veghghtrow)
        
        tneeded(:) = 0  !set to 0 now that we have made room.

        do i=1,nlat 
          if(timharvrow(i) > 0.0) then
            !take the area needed from the oldest active tile)
            inputVector(i,:) = 0
            ageSort = indexSort(tileAgerow(i,:),nmos)

            areaTemp = 0.0
            m = 1

            !keep setting valid tiles as input tiles until there's enough input tile area to cary out a split
            wloop4: do while(areaTemp < timharvrow(i))
              if(m > nmos) then
                write(6, * ) 'tiledDisturbancePrep: not enough tile area'
                call errorHandler('tiledDisturbancePrep', -3)
                !this exit statement is important or else the code can get stuck here
                exit wloop4
              end if
              if(controlVector(i,ageSort(m)) /= 0 .and. FAREROT(i,ageSort(m)) > 0) then 
                inputVector(i,ageSort(m)) = 1
                areaTemp = areaTemp + FAREROT(i,ageSort(m))
              end if 
              m = m + 1
            end do wloop4

            FAREAROUT(i) = min(dynTilingMaxDisturb, timharvrow(i))
            mode(i) = "split"
          else
            mode(i) = "skip" !this skip option is in here for nlat > 1, but it should never be used offline
          end if
        end do

        !now we do the actual splits
        call moveSplitCopyTiles(controlVector,inputVector,FAREAROUT,mode,DynTilInitializeFlag,outputIndex, &
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
                            useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow)

        do i=1,nlat 
          if(timharvrow(i) > 0.0) then
            !now we pass forward information about what harvesting will actually occur later on
            timharvarearow(i,:) = 0.0
            timharvarearow(i,outputIndex(i)) = 1.0
            prsfirearearow(i,:) = 0.0 ! because we defined timhar as having no prescribed fire
          else
            !if these no disturbance we make certain timharvarearow is set to zero 
            timharvarearow(i,:) = 0.0
            prsfirearearow(i,:) = 0.0 ! because we defined timhar as having no prescribed fire
          end if 
        end do

    case ('prefire') ! Prescribed fire --------------------------------

        tneeded(:) = 0
        do i=1,nlat 
          if((prsfirerow(i) > 0.0) .and. (.not. any((controlVector(i,:) /= 0 .and. FAREROT(i,:) == 0.0)))) then
            tneeded(i) = 1
          end if
        end do

        call cleanTiles(controlVector,tneeded,DynTilInitializeFlag,outputIndex, & !these are the main intputs that control the routine
                  ALBSROT, ALICROT, ALVCROT, CLAYROT, CMASROT, Cmossmasrow, DRNROT, FAREROT, FCANROT, & !these are just other major variables
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
                  useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow, veghghtrow)
        
        tneeded(:) = 0  !set to 0 now that we have made room.

        do i=1,nlat 
          if(prsfirerow(i) > 0.0) then
            !this just adds a small amount of area to account for the models representation of area burned and the fact that crops don't burn
            !this calculation appears a few times but with subtle differences depending on the job options selected
            !when dynamic tiling is on this has to be done using characteristics from the head tile (i.e. FINDLOC(controlVector(i,:),VALUE = -1, DIM = 1))
            !because harvest alone is being applied here there's not need to further checks after FAREAROUT is determined
            fcancmxcrop(i) = 0.0
            do j = 1,icc
              if (crop(j)) then
                fcancmxcrop(i) = fcancmxcrop(i) + fcancmxrow(i,FINDLOC(controlVector(i,:),VALUE = -1, DIM = 1),j)
              end if
            end do
            FAREAROUT(i) = min(dynTilingMaxDisturb, (prsfirerow(i) / (-fcancmxcrop(i) + sum(fcancmxrow(i,FINDLOC(controlVector(i,:),VALUE = -1, DIM = 1),:))))) 
            mode(i) = "split" 
            
            !<select and rank tiles to be burned using the bterm from CLASSIC's prognostic fire scheme
            call fireRank(gleafmasrow,bleafmasrow,stemmassrow,litrmassrow,fcancmxrow,bterm_veg) 
            btermSort = indexSort(bterm_veg(i,:),nmos)

            !< first mark all the tiles with a bterm of one
            inputVector(i,:) = 0
            areaTemp = 0.0
            do m=1,nmos
              if(bterm_veg(i,m) >= 1.0 .and. controlVector(i,m) /= 0 .and. FAREROT(i,m) > 0.0) then
                inputVector(i,m) = 1
                areaTemp = areaTemp + FAREROT(i,m)
              end if
            end do

            !< second add more tiles ranked using bterm if there's not enough area to burn
            !! with the already marked tiles.
            m = 1

            wloop6: do while(areaTemp < FAREAROUT(i))
              if(m > nmos) then
                write(6, * ) 'tiledDisturbancePrep: not enough tile area'
                write(6, * ) 'FAREAROUT = ', FAREAROUT(i)
                write(6, * ) 'areaTemp = ', areaTemp
                call errorHandler('tiledDisturbancePrep', -6)
                !this exit statement is important or else the code can get stuck here
                exit wloop6
              end if
              if(controlVector(i,btermSort(m)) /= 0 .and. FAREROT(i,btermSort(m)) > 0.0 .and. inputVector(i,btermSort(m)) /= 1) then 
                inputVector(i,btermSort(m)) = 1
                areaTemp = areaTemp + FAREROT(i,btermSort(m))
              end if 
              m = m + 1
            end do wloop6

          else
            mode(i) = "skip" !this skip option is in here for nlat > 1, but it should never be used offline
            FAREAROUT(i) = 0.0
          end if
        end do

        !now we do the actual splits
        call moveSplitCopyTiles(controlVector,inputVector,FAREAROUT,mode,DynTilInitializeFlag,outputIndex, &
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
                            useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow)

        do i=1,nlat 
          if(prsfirerow(i) > 0.0) then
            !now we pass forward information about what burning will actually occur later on
            prsfirearearow(i,:) = 0.0
            prsfirearearow(i,outputIndex(i)) = 1.0
            timharvarearow(i,:) = 0.0 ! because we defined this as only a prescribed fire operation, not timber harvest allowed.
          else
            !if these no disturbance we make certain prsfirearearow is set to zero 
            prsfirearearow(i,:) = 0.0
            timharvarearow(i,:) = 0.0 ! because we defined this as only a prescribed fire operation, not timber harvest allowed.
          end if 
        end do

    case ('firehar')  ! Prescribed fire and timber harvest ------------------------------

        do i=1,nlat 
          if ((prsfirerow(i) + timharvrow(i)) > 1.0) then
              write(6, * ) 'tiledDisturbancePrep: distrubances exceed avaliable tile area'
              call errorHandler('tiledDisturbancePrep', -4)
          end if
        end do

        ! Find the number of tiles needed for the upcoming timber harvest and prescribed fire.
        tneeded(:) = 0
        thave(:) = 0
        do i=1,nlat 

          if((prsfirerow(i) > 0.0) .and. (count((controlVector(i,:) /= 0 .and. FAREROT(i,:) == 0.0)) < 1)) then
            tneeded(i) = tneeded(i) + 1
          else if ((prsfirerow(i) > 0.) .and. (count((controlVector(i,:) /= 0 .and. FAREROT(i,:) == 0.0)) >= 1)) then
            thave(i) = thave(i) + 1
          end if

          if ((timharvrow(i) > 0.0) .and. ((count((controlVector(i,:) /= 0 .and. FAREROT(i,:) == 0.0)) - thave(i)) < 1)) then
            tneeded(i) = tneeded(i) + 1
          end if 

        end do

        call cleanTiles(controlVector,tneeded,DynTilInitializeFlag,outputIndex, & !these are the main intputs that control the routine
                  ALBSROT, ALICROT, ALVCROT, CLAYROT, CMASROT, Cmossmasrow, DRNROT, FAREROT, FCANROT, & !these are just other major variables
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
                  useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow, veghghtrow)
        
        tneeded(:) = 0   !set to 0 now that we have made room.

        do i=1,nlat 
          if (prsfirerow(i) > 0.0) then
            !this just adds a small amount of area to account for the models representation of area burned and the fact that crops don't burn
            !this calculation appears a few times but with subtle differences depending on the job options selected
            !when dynamic tiling is on this has to be done using characteristics from the head tile (i.e. FINDLOC(controlVector(i,:),VALUE = -1, DIM = 1))
            !because harvest and fire are being applied here a bunch of checks are run on FAREAROUT
            fcancmxcrop(i) = 0.0
            do j = 1,icc
              if (crop(j)) then
                fcancmxcrop(i) = fcancmxcrop(i) + fcancmxrow(i,FINDLOC(controlVector(i,:),VALUE = -1, DIM = 1),j)
              end if
            end do
            FAREAROUT(i) = min(dynTilingMaxDisturb, (prsfirerow(i) / (-fcancmxcrop(i) + sum(fcancmxrow(i,FINDLOC(controlVector(i,:),VALUE = -1, DIM = 1),:)))))
            mode(i) = "split" 

            if(FAREAROUT(i) + timharvrow(i) > dynTilingMaxDisturb) then

              !If this trips it is likely because the burned fraction requested in the input file exceeds 
              !the about to be harvested/non-vegetated/non-crop area avaliable in the grid cell
              !as a result fire is downscaled and priority is given to harvest because it is stand replacing

              write(6, * ) '########################################'
              write(6, * ) 'warning in tiledDisturbancePrep: disturbance was corrected as to not exceed avaliable land area'
              write(6, * )  "prior requested total area = ", timharvrow(i) + FAREAROUT(i) 
              write(6, * )  "harvested area = ", timharvrow(i)
              write(6, * )  "burned area = ", FAREAROUT(i) 

              if (timharvrow(i) < 0.0) then
                timharvrow(i) = 0.0
              else if (timharvrow(i) > dynTilingMaxDisturb) then
                timharvrow(i) = dynTilingMaxDisturb
              end if

              FAREAROUT(i) = (dynTilingMaxDisturb - timharvrow(i))

              if (FAREAROUT(i) < 0.0) then
                FAREAROUT(i) = 0.0
              else if (FAREAROUT(i) > dynTilingMaxDisturb) then
                FAREAROUT(i) = dynTilingMaxDisturb
              end if

              write(6, * )  "burned area used = ", FAREAROUT(i)
              write(6, * )  "harvested area used = ", timharvrow(i)
              write(6, * )  "new requested total area = ", FAREAROUT(i) + timharvrow(i)
              write(6, * ) '########################################'

              if ((FAREAROUT(i) + timharvrow(i)) > dynTilingMaxDisturb) then
                  write(6, * ) 'tiledDisturbancePrep: distrubance rescaling failed'
                  call errorHandler('tiledDisturbancePrep', -5)
              end if

            end if 

            !<select and rank tiles to be burned using the bterm
            call fireRank(gleafmasrow, bleafmasrow,stemmassrow,litrmassrow,fcancmxrow,bterm_veg)
            btermSort = indexSort(bterm_veg(i,:),nmos)

            !< first mark all the tiles with a bterm above the threshold of one
            inputVector(i,:) = 0
            areaTemp = 0.0
            do m=1,nmos
              if(bterm_veg(i,m) >= 1.0 .and. controlVector(i,m) /= 0 .and. FAREROT(i,m) > 0.0) then
                inputVector(i,m) = 1
                areaTemp = areaTemp + FAREROT(i,m)
              end if
            end do

            !< second add more tiles ranked using bterm if there's not enough area to burn, this must skip tiles already marked above
            m = 1
            wloop7: do while(areaTemp < FAREAROUT(i))
              if(m > nmos) then
                write(6, * ) 'test result = ', areaTemp < FAREAROUT(i)
                write(6, * ) 'remainder = ', areaTemp-FAREAROUT(i)
                write(6, * ) 'tiledDisturbancePrep: not enough tile area'
                write(6, * ) 'FAREAROUT = ', FAREAROUT(i)
                write(6, * ) 'areaTemp = ', areaTemp
                write(6, * ) 'tile itteration = ', m
                write(6, * ) 'nmos = ', nmos
                !write(6, * ) 'FAREROT = ', FAREROT(i,:)
                !write(6, * ) 'controlVector = ', controlVector(i,:)
                !write(6, * ) 'inputVector = ', inputVector(i,:)
                call errorHandler('tiledDisturbancePrep', -8)
                !this exit statement is important or else the code can get stuck here
                exit wloop7
              end if
              if(controlVector(i,btermSort(m)) /= 0 .and. FAREROT(i,btermSort(m)) > 0.0 .and. inputVector(i,btermSort(m)) /= 1) then 
                inputVector(i,btermSort(m)) = 1
                areaTemp = areaTemp + FAREROT(i,btermSort(m))
              end if 
              m = m + 1
            end do wloop7

          else
            mode(i) = "skip" !this skip option is in here for nlat > 1, but it should never be used offline
            FAREAROUT(i) = 0.0
          end if
        end do

        !now we do the actual splits
        call moveSplitCopyTiles(controlVector,inputVector,FAREAROUT,mode,DynTilInitializeFlag,outputIndex, &
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
                            useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow)

        do i=1,nlat 
          if(prsfirerow(i) > 0.0) then
            !now we pass forward information about what burning will actually occur later on
            prsfirearearow(i,:) = 0.0
            prsfirearearow(i,outputIndex(i)) = 1.0
          else
            !if these no disturbance we make certain prsfirearearow is set to zero 
            prsfirearearow(i,:) = 0.0
          end if 
        end do

        do i=1,nlat 
          if (timharvrow(i) > 0.0) then
            !now we take the area for harvest from the oldest active tile excluding the last output tile
            inputVector(i,:) = 0
            ageSort = indexSort(tileAgerow(i,:),nmos)

            areaTemp = 0.0
            m = 1

            !keep setting valid tiles as input tiles until there's enough input tile area to cary out a split
            wloop5: do while(areaTemp < timharvrow(i))
              if(m > nmos) then
                write(6, * ) 'tiledDisturbancePrep: not enough tile area'
                call errorHandler('tiledDisturbancePrep', -5)
                !this exit statement is important or else the code can get stuck here
                exit wloop5
              end if
              if(controlVector(i,ageSort(m)) /= 0 .and. FAREROT(i,ageSort(m)) > 0.0 .and. m /= outputIndex(i)) then 
                inputVector(i,ageSort(m)) = 1
                areaTemp = areaTemp + FAREROT(i,ageSort(m))
              end if 
              m = m + 1
            end do wloop5

            FAREAROUT(i) = min(dynTilingMaxDisturb, timharvrow(i))
            mode(i) = "split"

          else
            mode(i) = "skip" !this skip option is in here for nlat > 1, but it should never be used offline
            FAREAROUT(i) = 0
          end if
        end do

        !now we do the actual splits
        call moveSplitCopyTiles(controlVector,inputVector,FAREAROUT,mode,DynTilInitializeFlag,outputIndex, &
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
                            useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow)

        do i=1,nlat 
          if(timharvrow(i) > 0.0) then
            !now we pass forward information about what harvesting will actually occur later on
            timharvarearow(i,:) = 0.0
            timharvarearow(i,outputIndex(i)) = 1.0
          else
            !if these no disturbance we make certain timharvarearow is set to zero 
            timharvarearow(i,:) = 0.0
          end if 
        end do

    case('notdyn')  ! Prescribed fire or timber harvest but without dyanmic tiling -----------------

        do i=1,nlat 
          if ((prsfirerow(i) + timharvrow(i)) > 1.0) then
            write(6, * ) 'tiledDisturbancePrep: distrubances exceed avaliable tile area'
            call errorHandler('tiledDisturbancePrep', -4)
          end if
        end do

        do i=1,nlat 
            timharvarearow(i,:) = 0.0
            prsfirearearow(i,:) = 0.0
            timharvarearow(i,1) = min(dynTilingMaxDisturb, timharvrow(i))
            !this just adds a small amount of area to account for the models representation of area burned and the fact that crops don't burn
            !this calculation appears a few times but with subtle differences depending on the job options selected
            !when dynamic tiling is off there is only a single tile involved and the value is applied directly to prsfirearearow so the
            !"average individual" can be harvested because harvest and fire are being applied here a bunch of checks are then run on prsfirearearow
            fcancmxcrop(i) = 0.0
            do j = 1,icc
              if (crop(j)) then
                fcancmxcrop(i) = fcancmxcrop(i) + fcancmxrow(i,1,j)
              end if
            end do
            prsfirearearow(i,1) = min(dynTilingMaxDisturb, (prsfirerow(i) / (-fcancmxcrop(i) + sum(fcancmxrow(i,1,:)))))

            if((prsfirearearow(i,1) + timharvarearow(i,1)) > dynTilingMaxDisturb) then

              !If this trips it is likely becuase the burned fraction requested in the input file exceeds 
              !the not about to be harvested/non-vegetated/non-crop area avaliable in the grid cell
              !as a result fire is downscaled and priority is given to harvest because it is stand replacing

              write(6, * ) '########################################'
              write(6, * ) 'warning in tiledDisturbancePrep: disturbance was corrected as to not exceed avaliable land area'
              write(6, * ) "prior requested total area = ", timharvarearow(i,1) + prsfirearearow(i,1)
              write(6, * ) "harvested area = ", timharvarearow(i,1)
              write(6, * )  "burned area = ", prsfirearearow(i,1)
              
              if (timharvarearow(i,1) < 0.0) then
                timharvarearow(i,1) = 0.0
              else if (timharvarearow(i,1) > dynTilingMaxDisturb) then
                timharvarearow(i,1) = dynTilingMaxDisturb
              end if

              prsfirearearow(i,1) = (dynTilingMaxDisturb - timharvarearow(i,1))

              if (prsfirearearow(i,1) < 0.0) then
                prsfirearearow(i,1) = 0.0
              else if (prsfirearearow(i,1) > dynTilingMaxDisturb) then
                prsfirearearow(i,1) = dynTilingMaxDisturb
              end if

              write(6, * )  "burned area used = ", prsfirearearow(i,1)
              write(6, * )  "harvested area used = ", timharvarearow(i,1)
              write(6, * )  "new requested total area = ", timharvarearow(i,1) + prsfirearearow(i,1)
              write(6, * ) '########################################'

              if ((prsfirearearow(i,1) + timharvarearow(i,1)) > dynTilingMaxDisturb) then
                write(6, * ) 'tiledDisturbancePrep: distrubance rescaling failed'
                call errorHandler('tiledDisturbancePrep', -5)
              end if
            end if
        end do 

    case ('cleandyn')  ! run cleanTiles to pre-emptively join tiles in years with no harvest or fire prescribed

        tneeded(:) = 0

        call cleanTiles(controlVector,tneeded,DynTilInitializeFlag,outputIndex, & !these are the main intputs that control the routine
                  ALBSROT, ALICROT, ALVCROT, CLAYROT, CMASROT, Cmossmasrow, DRNROT, FAREROT, FCANROT, & !these are just other major variables
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
                  useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow, veghghtrow)

    case ('noaction') ! The defeault return, nothing needs to happen

    case ('reset') ! On the second day of the year leave the tiles as is but make sure that timharvarearow is set to zero so extra biomass isn't harvested
        timharvarearow(:,:) = 0.0
        timharvareagat(:) = 0.0
        prsfirearearow(:,:) = 0.0
        prsfireareagat(:) = 0.0
    case default
        write(6, * ) 'tiledDisturbancePrep: something strange happened while selecting the case'
        call errorHandler('tiledDisturbancePrep', -7)
    end select

  end subroutine tiledDisturbancePrep
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup tiledDisturbance_initialMassBalance
  !! @{
  !> This function calulates the intial total mass of carbon, water, etc in a grid cell under consideration
  !> @author S.R. Curasi
  subroutine initialMassBalance(fcancmx, gleafmas, bleafmas, stemmass, rootmass, &
                                soilcmas, litrmass, nlitrmass, soilnmas, &
                                ngleafmas, nbleafmas ,nstemmass, nrootmass, & 
                                vegCMass_initial, soilCMass_initial, lucMass_initial, totalCMass_initial, &
                                vegNMass_initial, soilNMass_initial, lucNMass_initial, totalNMass_initial, &
                                Ncycle_on, faregat)

    use, intrinsic :: iso_fortran_env, only: r8=>real64

    implicit none 

    !these are used or calcualted within the function
    integer :: i, k, j, l
    real :: barefrac(ilg)

    !these are important descriptors of the land surface taken by the function
    real, intent(in), dimension(ilg,icc) :: fcancmx
    logical, intent(in) :: Ncycle_on
    real, intent(in), dimension(ilg) :: faregat

    !vegetation biomass related variables taken by mass balance checks
    real, intent(in), dimension(ilg,icc) :: gleafmas
    real, intent(in), dimension(ilg,icc) :: bleafmas
    real, intent(in), dimension(ilg,icc) :: stemmass
    real, intent(in), dimension(ilg,icc) :: rootmass

    !soil and litter variables taken by mass balance checks
    real, intent(in), dimension(ilg,iccp2,ignd) :: soilcmas
    real, intent(in), dimension(ilg,iccp2,ignd) :: litrmass

    !ncycling related variables taken by mass balance checks
    real, intent(in), dimension(ilg,iccp2) :: nlitrmass
    real, intent(in), dimension(ilg,iccp2) :: soilnmas
    real, intent(in), dimension(ilg,icc) :: ngleafmas
    real, intent(in), dimension(ilg,icc) :: nbleafmas
    real, intent(in), dimension(ilg,icc) :: nstemmass
    real, intent(in), dimension(ilg,icc) :: nrootmass

    !aggregated outputs from mass balance checks
    real, intent(inout), dimension(ilg) :: vegCMass_initial !< temporary storage variable for the total mass of carbon in vegetation
    real, intent(inout), dimension(ilg) :: soilCMass_initial !< temporary storage variable for the total mass of carbon in soil and litter
    real, intent(inout), dimension(ilg) :: lucMass_initial !< temporary storage variable for the total mass of carbon in the LUC pools
    real, intent(inout), dimension(ilg) :: totalCMass_initial !< temporary storage variable for the total mass of carbon in all pools

    real, intent(inout), dimension(ilg) :: vegNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, intent(inout), dimension(ilg) :: soilNMass_initial !< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, intent(inout), dimension(ilg) :: lucNMass_initial !< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, intent(inout), dimension(ilg) :: totalNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation

    vegCMass_initial = 0.0
    soilCMass_initial = 0.0
    lucMass_initial = 0.0
    totalCMass_initial = 0.0

    vegNMass_initial = 0.0
    soilNMass_initial = 0.0
    lucNMass_initial = 0.0
    totalNMass_initial = 0.0

    do i=1,ilg
      if (faregat(i) > 0.0) then
        barefrac(i) = 1.0 - sum(fcancmx(i,:))
        
        do k = 1,icc
          vegCMass_initial(i) =  vegCMass_initial(i) + &
                              (gleafmas(i,k) * fcancmx(i,k)) + &
                              (bleafmas(i,k) * fcancmx(i,k)) + &
                              (stemmass(i,k) * fcancmx(i,k)) + &
                              (rootmass(i,k) * fcancmx(i,k))

          if (Ncycle_on) then
            vegNMass_initial(i) =  vegNMass_initial(i) + &
                                (ngleafmas(i,k) * fcancmx(i,k) * convertg2kg) + &
                                (nbleafmas(i,k) * fcancmx(i,k) * convertg2kg) + &
                                (nstemmass(i,k) * fcancmx(i,k) * convertg2kg) + &
                                (nrootmass(i,k) * fcancmx(i,k) * convertg2kg)
          end if
        end do
        
        do k = 1,ignd
          do l=1,icc 
            soilCMass_initial(i) = soilCMass_initial(i) + &
                                (soilcmas(i,l,k) * fcancmx(i,l)) + &
                                (litrmass(i,l,k) * fcancmx(i,l))                   
          end do
          
          soilCMass_initial(i) = soilCMass_initial(i) + &
                                (soilcmas(i,iccp1,k) * barefrac(i)) + &
                                (litrmass(i,iccp1,k) * barefrac(i))  
        end do

        lucMass_initial(i) = lucMass_initial(i) + &
                            (litrmass(i,iccp2,1)) + &
                            (soilcmas(i,iccp2,1))  
        
        if (Ncycle_on) then 
          do l=1,icc
            soilNMass_initial(i) = soilNMass_initial(i) + &
                                (soilnmas(i,l) * fcancmx(i,l) *  convertg2kg) + &
                                (nlitrmass(i,l) * fcancmx(i,l) *  convertg2kg)
          end do
          
          soilNMass_initial(i) = soilNMass_initial(i) + &
                                (soilnmas(i,iccp1) * barefrac(i)  * convertg2kg) + &
                                (nlitrmass(i,iccp1) * barefrac(i)  * convertg2kg)
          lucNMass_initial(i) = lucMass_initial(i) + &
                                (soilnmas(i,iccp2)  * convertg2kg) + &
                                (nlitrmass(i,iccp2)  * convertg2kg) 
        end if
      end if
    end do

    totalCMass_initial = vegCMass_initial + soilCMass_initial + lucMass_initial
    totalNMass_initial = soilNMass_initial + vegNMass_initial + lucNmass_initial

  end subroutine initialMassBalance
  !! @}

  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup tiledDisturbance_finalMassBalance
  !! @{
  !> This function calulates the final total mass of carbon, water, etc in a grid cell under consideration
  !> @author S.R. Curasi
  subroutine finalMassBalance(fcancmx, gleafmas, bleafmas, stemmass, rootmass, &
                                soilcmas, litrmass, lucemcom, nlitrmass, soilnmas, &
                                ngleafmas, nbleafmas ,nstemmass, nrootmass, lucemcomn, & 
                                vegCMass_initial, soilCMass_initial, lucMass_initial, totalCMass_initial, &
                                vegNMass_initial, soilNMass_initial, lucNMass_initial, totalNMass_initial, &
                                Ncycle_on, faregat)

    use, intrinsic :: iso_fortran_env, only: r8=>real64

    implicit none 

    !these are used or calcualted within the function
    integer :: i, k, j, l
    real :: barefrac(ilg)

    !these are important descriptors of the land surface taken by the function
    real, intent(in), dimension(ilg,icc) :: fcancmx
    logical, intent(in) :: Ncycle_on
    real, intent(in), dimension(ilg) :: faregat

    !vegetation biomass related variables taken by mass balance checks
    real, intent(in), dimension(ilg,icc) :: gleafmas
    real, intent(in), dimension(ilg,icc) :: bleafmas
    real, intent(in), dimension(ilg,icc) :: stemmass
    real, intent(in), dimension(ilg,icc) :: rootmass

    !soil and litter variables taken by mass balance checks
    real, intent(in), dimension(ilg,iccp2,ignd) :: soilcmas
    real, intent(in), dimension(ilg,iccp2,ignd) :: litrmass

    !fire flux taken bas mass balance checks
    real, intent(in), dimension(ilg) :: lucemcom !< luc related carbon emission losses from combustion (u-mol CO2-C/m2.sec)

    !ncycling related variables taken by mass balance checks
    real, intent(in), dimension(ilg,iccp2) :: nlitrmass
    real, intent(in), dimension(ilg,iccp2) :: soilnmas
    real, intent(in), dimension(ilg,icc) :: ngleafmas
    real, intent(in), dimension(ilg,icc) :: nbleafmas
    real, intent(in), dimension(ilg,icc) :: nstemmass
    real, intent(in), dimension(ilg,icc) :: nrootmass

    !n cycling realted fire flux taken in by mass balance checks
    real, intent(in), dimension(ilg) :: lucemcomn !< luc related carbon emission losses from combustion of N (g N m^-2 day^-1)

    !aggregated outputs from mass balance checks
    real, intent(in), dimension(ilg) :: vegCMass_initial !< temporary storage variable for the total mass of carbon in vegetation
    real, intent(in), dimension(ilg) :: soilCMass_initial !< temporary storage variable for the total mass of carbon in soil and litter
    real, intent(in), dimension(ilg) :: lucMass_initial !< temporary storage variable for the total mass of carbon in the LUC pools
    real, intent(in), dimension(ilg) :: totalCMass_initial !< temporary storage variable for the total mass of carbon in all pools

    real, intent(in), dimension(ilg) :: vegNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, intent(in), dimension(ilg) :: soilNMass_initial !< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, intent(in), dimension(ilg) :: lucNMass_initial !< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, intent(in), dimension(ilg) :: totalNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation

    !aggregated outputs from mass balance checks
    real, dimension(ilg) :: vegCMass_final !< temporary storage variable for the total mass of carbon in vegetation
    real, dimension(ilg) :: soilCMass_final !< temporary storage variable for the total mass of carbon in soil and litter
    real, dimension(ilg) :: lucMass_final !< temporary storage variable for the total mass of carbon in the LUC pools
    real, dimension(ilg) :: totalCMass_final !< temporary storage variable for the total mass of carbon in all pools

    real, dimension(ilg) :: vegNMass_final !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, dimension(ilg) :: soilNMass_final !< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, dimension(ilg) :: lucNMass_final !< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, dimension(ilg) :: totalNMass_final !< temporary storage variable for the total mass of nitrogen in the vegetation

    real, dimension(ilg) :: lucFireFlux
    real, dimension(ilg) :: lucNFireFlux

    vegCMass_final = 0.0
    soilCMass_final = 0.0
    lucMass_final = 0.0
    totalCMass_final = 0.0

    vegNMass_final = 0.0
    soilNMass_final = 0.0
    lucNMass_final = 0.0
    totalNMass_final = 0.0

    lucFireFlux = 0.0
    lucNFireFlux = 0.0

    do i=1,ilg
      if(faregat(i) > 0.0) then
        barefrac(i) = 1.0 - sum(fcancmx(i,:))
        
        do k = 1,icc
          vegCMass_final(i) =  vegCMass_final(i) + &
                              (gleafmas(i,k) * fcancmx(i,k)) + &
                              (bleafmas(i,k) * fcancmx(i,k)) + &
                              (stemmass(i,k) * fcancmx(i,k)) + &
                              (rootmass(i,k) * fcancmx(i,k))

          if (Ncycle_on) then
            vegNMass_final(i) =  vegNMass_final(i) + &
                                (ngleafmas(i,k) * fcancmx(i,k) * convertg2kg) + &
                                (nbleafmas(i,k) * fcancmx(i,k) * convertg2kg) + &
                                (nstemmass(i,k) * fcancmx(i,k) * convertg2kg) + &
                                (nrootmass(i,k) * fcancmx(i,k) * convertg2kg)
          end if
        end do
        
        do k = 1,ignd
          do l=1,icc
            soilCMass_final(i) = soilCMass_final(i) + &
                                (soilcmas(i,l,k) * fcancmx(i,l)) + &
                                (litrmass(i,l,k) * fcancmx(i,l))                   
          end do
          
          soilCMass_final(i) = soilCMass_final(i) + &
                                (soilcmas(i,iccp1,k) * barefrac(i)) + &
                                (litrmass(i,iccp1,k) * barefrac(i))  

        end do

        lucMass_final(i) = lucMass_final(i) + &
                          (litrmass(i,iccp2,1)) + &
                          (soilcmas(i,iccp2,1))  
        
        lucFireFlux(i) = lucemcom(i)/963.62
        
        if (Ncycle_on) then 
          do l=1,icc
            soilNMass_final(i) = soilNMass_final(i) + &
                                (soilnmas(i,l) * fcancmx(i,l) *  convertg2kg) + &
                                (nlitrmass(i,l) * fcancmx(i,l) *  convertg2kg)
          end do
          
          soilNMass_final(i) = soilNMass_final(i) + &
                                (soilnmas(i,iccp1) * barefrac(i)  * convertg2kg) + &
                                (nlitrmass(i,iccp1) * barefrac(i)  * convertg2kg)
          lucNMass_final(i) = lucMass_final(i) + &
                                (soilnmas(i,iccp2)  * convertg2kg) + &
                                (nlitrmass(i,iccp2)  * convertg2kg) 
          
          lucNFireFlux(i) = lucemcomn(i)

        end if
      end if
    end do

    totalCMass_final = vegCMass_final + soilCMass_final + lucMass_final
    totalNMass_final = soilNMass_final + vegNMass_final + lucNmass_final

    do i=1,ilg
      if(faregat(i) > 0.0) then
        if ((abs((totalCMass_final(i) + lucFireFlux(i)) - totalCMass_initial(i)) >= dynTilingTolrance) .or. &
            (abs((totalNMass_final(i) + lucNFireFlux(i)) - totalNMass_initial(i)) >= dynTilingTolrance)) then
            ! Note: currently C and N use the same dynTilingTolrance without issue. Both pools are tightly conserved in dynamic tiling
            ! operations with any differences stemming from numerical precision. As dynamic tiling with N is tested more extensively it may
            ! make sense to have two separate thresholds  
          print *, "error in finalMassBalance checks"
          print *, "---mass balence check results---"
          print *, "ilg = ",i
          print *, "---carbon"
          print *, 'vegCMass_initial =', vegCMass_initial(i), 'vegCMass_final =', vegCMass_final(i)
          print *, 'soilCMass_initial =', soilCMass_initial(i), 'soilCMass_final =', soilCMass_final(i)
          print *, 'lucMass_initial =', lucMass_initial(i), 'lucMass_final =', lucMass_final(i)
          print *, 'lucFireFlux = ', lucFireFlux(i)
          print *, 'totalCMass_initial = ', totalCMass_initial(i), ' totalCMass_final = ', totalCMass_final(i), 'totalCMassFire_final = ', (totalCMass_final(i) + lucFireFlux(i)),'diff = ', (totalCMass_final(i) + lucFireFlux(i)) - totalCMass_initial(i)
          if (Ncycle_on) then 
            print *, "---nitrogen"
            print *, 'vegNMass_initial =', vegNMass_initial(i), 'vegNMass_final =', vegNMass_final(i)
            print *, 'soilNMass_initial =', soilNMass_initial(i), 'soilNMass_final =', soilNMass_final(i)
            print *, 'lucNMass_initial =', lucNMass_initial(i), 'lucNMass_final =', lucNMass_final(i)
            print *, 'lucNFireFlux = ', lucNFireFlux(i)
            print *, 'totalNMass_initial = ', totalNMass_initial(i), ' totalNMass_final = ', totalNMass_final(i), 'totalNMassFire_final = ', (totalNMass_final(i) + lucNFireFlux(i)),'diff = ', (totalNMass_final(i) + lucNFireFlux(i)) - totalNMass_initial(i)
            print *, "--------------------------------"
          end if
          call errorHandler('finalMassBalance', - 1)
        end if
      end if
    end do

  end subroutine finalMassBalance
  !! @}

  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup tiledDisturbance_fireRank
  !! @{
  !> This subroutine uses the total biomass term [bterm_veg] from CLASSIC's prognostic fire model averaged 
  !! for all the PFTs in the gridcell to rank tiles to be burned
  !!
  !> @author S.R. Curasi
  subroutine fireRank(gleafmasrow, bleafmasrow,stemmassrow,litrmassrow,fcancmxrow,bterm_veg)

    use classicParams, only : bmasthrs_fire, crop

    implicit none 

    real, intent(in), dimension(nlat,nmos,icc) :: gleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: bleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: stemmassrow

    real, intent(in), dimension(nlat,nmos,icc,ignd) :: litrmassrow

    real, intent(in), dimension(nlat,nmos,icc) :: fcancmxrow

    real, dimension(nlat,nmos,icc) :: biomass = 0.0
    real, dimension(nlat,nmos,icc) :: bterm_veg_pp = 0.0
    
    real, intent(out), dimension(nlat,nmos) :: bterm_veg

    integer :: i, j, k

    !> 1. Total biomass term

    bterm_veg = 0.0

    do i=1,nlat 
      do j=1,nmos
        do k=1,icc
          if (.not. crop(k)) then
            biomass(i,j,k) = gleafmasrow(i,j,k) + bleafmasrow(i,j,k) + stemmassrow(i,j,k) + sum(litrmassrow(i,j,k,:)) !biomass
          end if
          bterm_veg_pp(i,j,k) = min(1.0,max(0.0,(biomass(i,j,k) - bmasthrs_fire(1))/(bmasthrs_fire(2) - bmasthrs_fire(1)))) !actual bterm calculation
          bterm_veg(i,j) = bterm_veg(i,j) + bterm_veg_pp(i,j,k) * fcancmxrow(i,j,k)
        end do
        bterm_veg(i,j) = bterm_veg(i,j)/sum(fcancmxrow(i,j,:))
      end do
    end do
    
  end subroutine fireRank
  !! @}

  ! ---------------------------------------------------------------------------------------------------

  !> \namespace tiledDisturbance
  !! Module paired with dynamic tiling to simulate disturbance
  !! @author S.R. Curasi
  !!
  !! Overview
  !!
  !! These subroutines are designed to allow for prescribed fire and timber harvest with or without dynamic tiling to represent the resulting subgrid-scale heterogeneity. 
  !! These operations allow biomass to be harvested or burned basses upon inputs read in from a file. They call the subroutines within dynmaicTiling to construct
  !! the tiled version of the simulation. The burning of grid cells for prescribed fire is carried out via some minor modifications to the subroutines in disturbance.f90
  !!
  !! Constituent functions and subroutine
  !!
  !! 1. harvestTile: This subroutine harvest biomass from a tile-based upon inputs read in from a file. It is called in ctemDriver.f90. It takes as its input timharvareagat
  !! which is the fraction area (i.e. 1 - 0) of a given tile that should be subject to harvest. If dynamic tiling is on harvest usually impacts the entire area of the desired tile. However if dynamic
  !! tiling is off just a portion of the single grid cell is harvested. After the tile is harvested it also alters its age [tileAgegat] as needed. It used the same parameters
  !! and storage pools as the default land use change module (i.e. Arora and Boer 2010).  
  !!
  !! 2. tiledDisturbancePrep: The subroutine preps the simulation of disturbance. If dynamic tiling is off if just read in information from the disturbance files and 
  !! prepares that information for later use by subroutines. If dynamic tiling is on it manages the tiles used by the simulation. It reads in the information from the disturbance
  !! files. Determines how many tiles the simulation needs and then runs cleanTiles and moveSplitCopyTiles as needed. It then sets up the information needed by the subroutines
  !! that carry out the disturbance simulation and calls for the model to be re-initialized.
  !!
  !! 3. initialMassBalance: The subroutine calculated the initial mass balance of the harvested tile for use later on.
  !!
  !! 4. finalMassBalance: The subroutine calculates the final mass balance of the harvested tile to ensure that no material has been lost when simulating the harvest.
  !!
  !! 5. fireRank: This subroutine uses the total biomass term [bterm_veg] for the standard fire model to rank tiles to be burned
  !!
  !!
  !> \file
end module tiledDisturbance