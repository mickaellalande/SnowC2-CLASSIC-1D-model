!> \file
!> Central module for all land use change operations
module landuseChange

  use classicParams

  implicit none

  ! subroutines contained in this module:
  public  :: initializeLandCover
  public  :: luc
  private :: adjustLucFracs
  private :: adjustFracsComp

contains

  !-------------------------------------------------------------------------------------------------------------
  !> \ingroup landusechange_initializeLandCover
  !> @{
  !> Initializes and checks over the PFT fractional cover.
  !> @author Joe Melton

  subroutine initializeLandCover
    !
    !    6   Oct  2017  - Rewrote now to only initialize and not read in the LUC since that
    !     J. Melton       is now done via netcdfs files. Also removed onetile_perPFT since we
    !                     no longer support that code option.
    !
    !    11  Jul  2016  - Further bug fixes for competing for space within a tile
    !     J. Melton
    !
    !     9  Mar  2016  - Adapt for tiling where we compete for space within a tile
    !     J. Melton
    !
    !     3  Feb  2016  - Remove mosaic flag and replace it with the onetile_perPFT one.
    !     J. Melton

    !     9  Jan. 2013  - this subroutine takes in luc information from
    !     J. Melton       a luc file and adapts them for runclassctem
    !                     it is run once to set up the luc info before the
    !                     model timestepping begins.
    !
    !     7  Feb. 2014  - Adapt it to work with competition and start_bare
    !     J. Melton
    !
    !     28 Mar. 2014  - Add in check to ensure we don't make a negative bare
    !     J. Melton       fraction when we add in seed fractions.
    !

    use classicParams,   only : nmos, nlat, icc, ican, icp1, seed, crop, numcrops, minbare, &
                                modelpft, l2max, nol2pfts, reindexPFTs
    use ctemStateVars,   only : c_switch, vrot
    use classStateVars,  only : class_rot

    implicit none

    ! local variables:
    real, dimension(nlat,nmos) :: barfm  !< bare fraction, used in competition runs to ensure that the bare ground fraction is at least seed.
    integer, dimension(1) :: bigpft
    real, dimension(nlat,nmos,icc - numcrops)    :: pftarrays ! temp variable
    integer, dimension(nlat,nmos,icc - numcrops) :: indexposj ! temp var
    integer, dimension(nlat,nmos,icc - numcrops) :: indexposm ! temp var
    integer :: j, m, i, n, k
    ! integer :: k2,k1

    associate( &
    start_bare        => c_switch%start_bare,           & !< logical: 
    PFTCompetition    => c_switch%PFTCompetition,       & !< logical: 
    FCANROT           => class_rot%FCANROT,             & !< real, dimension(:,:,:): 
    pfcancmxrow       => vrot%pfcancmx,                 & !< real, dimension(:,:,:): 
    nfcancmxrow       => vrot%nfcancmx,                 & !< real, dimension(:,:,:): 
    fcancmxrow        => vrot%fcancmx                   & !< real, dimension(:,:,:): 
    )

    !> -------------------------
    ! Some initilizations
    FCANROT = 0.0
    barfm = 1.0
    pftarrays = 0.0

    ! Set nfcancmxrow to fcancmxrow at start. fcancmxrow is what is initialized from the initial conditions file or
    ! from the LUC file.
    nfcancmxrow(:,:,:) = fcancmxrow(:,:,:)

    if (PFTCompetition) then
      do i = 1,nlat
        do m = 1,nmos
          k = 1  ! reset k for next tile
          do j = 1,icc
            if (.not. crop(j)) then
              if (start_bare) then
                nfcancmxrow(i,m,j) = seed
              else !> not starting bare,but still make sure you have at least seed
                nfcancmxrow(i,m,j) = max(seed,fcancmxrow(i,m,j))
              end if

              !> Keep track of the non-crop nfcancmx for use in loop below.
              !> pftarrays keeps track of the nfcancmxrow for all non-crops
              !> indexposj and indexposm store the index values of the non-crops
              !> in a continuous array for use later. k is then the indexes used by
              !> these arrays.
              pftarrays(i,m,k) = nfcancmxrow(i,m,j)
              indexposj(i,m,k) = j
              k = k + 1

            end if ! crops
            ! Reduce the bare fraction by the new PFT fractions
            barfm(i,m) = barfm(i,m) - nfcancmxrow(i,m,j)
          end do ! icc
        end do ! nmos
      end do ! nlat

      !> check that in making these seed fraction we haven't made our total fraction
      !> more than 1.0.
      do i = 1,nlat
        do m = 1,nmos
          if (barfm(i,m) < 0.) then

            !> Find out which of the non-crop PFTs covers the largest area.
            bigpft = maxloc(pftarrays(i,m,:))
            
            !> j is then the nmos index and m is the icc index of the PFT with the largest area
            j = indexposj(i,m,bigpft(1))
      
            !>
            !> Reduce the most dominant PFT by barf and minbare. The extra
            !! amount is to ensure we don't have trouble later with an extremely
            !! small bare fraction. barf is a negative value.
            !!
            nfcancmxrow(i,m,j) = nfcancmxrow(i,m,j) + barfm(i,m) - minbare

            if (nfcancmxrow(i,m,j) < seed) then 
              ! If now it is too small, add seed back and see if that was sufficient
              ! to make it none negative.
              nfcancmxrow(i,m,j) = nfcancmxrow(i,m,j) + seed
            
              ! now recheck to make sure it isn't negative. If so kill the run 
              ! and reassess the inputs fractions.
              if (nfcancmxrow(i,m,j) < 0.) then
                write(6, * )' lat/pft/tile',i,m,j,' is negative after adding seed to all. nfcancmx:',nfcancmxrow(i,m,j),' Check if initial values given are reasonable.'
                call errorHandler('initializeLandCover', - 1)
          end if
            end if 
          end if
        end do ! nmos
      end do ! nlat
    end if  ! PFTCompetition

    !> get fcans for use by class using the nfcancmxs just read in
    do j = 1,ican
      do n = reindexPFTs(j,1),reindexPFTs(j,2)
        do i = 1,nlat
          do m = 1,nmos
            FCANROT(i,m,j) = FCANROT(i,m,j) + nfcancmxrow(i,m,n)
          end do
        end do
      end do ! loop 998
    end do ! loop 997

    !> assign the present pft fractions from those just read in
    fcancmxrow(:,:,:) = nfcancmxrow(:,:,:)
    pfcancmxrow(:,:,:) = nfcancmxrow(:,:,:)

    end associate
    return

  end subroutine initializeLandCover
  !! @}
  !=======================================================================
  !> \ingroup landusechange_luc
  !> @{
  !> Deals with the changes in the land cover and estimates land use change (LUC)
  !! related carbon emissions. based on whether there is conversion of forests/grassland to crop area,
  !! or croplands abandonment, this subroutine reallocates carbon from live vegetation to litter
  !! and soil c components. the decomposition from the litter and soil c pools thus implicitly models luc
  !! related carbon emissions. set of rules are followed to determine the fate of carbon that
  !! results from deforestation or replacement of grasslands by crops.
  !> @author Vivek Arora, A. Asaadi
  subroutine luc (il1, il2, nilg, PFTCompetition, leapnow, useTracer, & ! In
                  grclarea, iday, todfrac, yesfrac, interpol, & ! In
                  pfcancmx, nfcancmx, gleafmas, gleafmas_ns, gleafmas_s, & ! In / Out
                  bleafmas, stemmass, stemmass_ns, stemmass_s, & ! In / Out
                  rootmass, rootmass_ns, rootmass_s, & ! In / Out
                  litrmass, soilcmas, vgbiomas, nvgbiomas, gavgltms, nlitrmass_g, & ! In / Out
                  gavgscms, soilnmas_g, nh4_mass_g, no3_mass_g, fcancmx, fcanmx, ngleafmas, ngleafmas_ns, & ! In / Out
                  ngleafmas_s, nbleafmas, nstemmass, nstemmass_ns, & ! In / Out
                  nstemmass_s, nrootmass, nrootmass_ns, nrootmass_s, & ! In / Out
                  nlitrmass, soilnmas, nh4_mass, no3_mass, & ! In / Out
                  tracerLitrMass, tracerSoilCMass, & ! In/Out
                  tracerGLeafMass, tracerBLeafMass, tracerStemMass, tracerRootMass, & ! In / Out
                  lucemcom, lucltrin, lucsocin, lucemcomn, lucltrinn, lucsocinn)                ! Out
    !
    !     ----------------------------------------------------------------
    !
    !     24  July 2019  - Add in N cycle related parts
    !     A. Asaadi
    !
    !      8  Feb 2016  - Adapted subroutine for multilayer soilc and litter (fast decaying)
    !     J. Melton       carbon pools
    !
    !     19  Jan 2016  - Implemented new LUC litter and soil C pools
    !     J. Melton

    !     31  Jan 2014  - Moved parameters to global file (classicParams.f90)
    !     J. Melton
    !
    !     18  Apr. 2013 - made it so that you will exit luc if the grid cell has
    !     J. Melton       no actual luc in this timestep. removed some extraneous checks.
    !
    !     02  Jan. 2004 - this subroutine deals with the changes in the land
    !     V. Arora        cover and estimates land use change (luc)
    !                     related carbon emissions. based on whether there
    !                     is conversion of forests/grassland to crop area,
    !                     or croplands abandonment, this subroutine
    !                     reallocates carbon from live vegetation to litter
    !                     and soil c components. the decomposition from the
    !                     litter and soil c pools thus implicitly models luc
    !                     related carbon emissions. set of rules are
    !                     followed to determine the fate of carbon that
    !                     results from deforestation or replacement of
    !                     grasslands by crops.
    !
    use classicParams,    only : icc, ican, zero, km2tom2, iccp1, &
                                 combust, paper, furniture, bmasthrs, &
                                 tolrnce1, tolrance, crop, numcrops, &
                                 minbare, iccp2, ignd, classpfts, icp1, &
                                 nol2pfts, reindexPFTs, &
                                 rel_tolerance_amount, rel_tolerance_density
    use ctemStatevars,     only : c_switch

    implicit none

    integer, intent(in) :: il1   !< il1=1
    integer, intent(in) :: il2   !< il2=nilg
    integer, intent(in) :: nilg  !< no. of grid cells in latitude circle(this is passed in as
    !! either ilg or nlat depending on comp/mos)
    integer, intent(in) :: iday       !< day of year
    logical, intent(in) :: leapnow    !< true if this year is a leap year. Only used if the switch 'leap' is true.
    logical, intent(in) ::  interpol  !< if todfrac & yesfrac are provided then interpol must be set to false so that
    !< this subroutine doesn't do its own interpolation using pfcancmx and nfcancmx
    !< which are year end values
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    real, intent(in) :: grclarea(nilg)       !< gcm grid cell area, km2
    logical, intent(in) :: PFTCompetition   !< true if the competition subroutine is on.
    real, intent(in) :: todfrac(nilg,icc)    !< today's fractional coverage of all pfts
    real, intent(in) :: yesfrac(nilg,icc)    !< yesterday's fractional coverage of all pfts

    real, intent(inout) :: gleafmas(nilg,icc)   !< green or live leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout) :: gleafmas_ns(nilg,icc) !< non-structural green or live leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout) :: gleafmas_s(nilg,icc)  !< structural green or live leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout) :: bleafmas(nilg,icc)   !< brown or dead leaf mass in kg c/m2, for the 9 pfts
    real, intent(inout) :: stemmass(nilg,icc)   !< stem biomass in kg c/m2, for the 9 pfts
    real, intent(inout) :: stemmass_ns(nilg,icc) !< non-structural stem biomass in kg c/m2, for the 9 pfts
    real, intent(inout) :: stemmass_s(nilg,icc)  !< structurla stem biomass in kg c/m2, for the 9 pfts
    real, intent(inout) :: rootmass(nilg,icc)   !< root biomass in kg c/m2, for the 9 pfts
    real, intent(inout) :: rootmass_ns(nilg,icc)   !< non-structural root biomass in kg c/m2, for the 9 pfts
    real, intent(inout) :: rootmass_s(nilg,icc)    !< structural root biomass in kg c/m2, for the 9 pfts 
    real, intent(inout) :: ngleafmas(nilg,icc)     !< green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: ngleafmas_ns(nilg,icc)  !< non-structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: ngleafmas_s(nilg,icc)   !< structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nbleafmas(nilg,icc)     !< brown leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nstemmass(nilg,icc)     !< stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nstemmass_ns(nilg,icc)  !< non-structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nstemmass_s(nilg,icc)   !< structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nrootmass(nilg,icc)     !< root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nrootmass_ns(nilg,icc)  !< non-structural root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nrootmass_s(nilg,icc)   !< structural root nitrogen mass for individual PFTs (\f$g N/m^2\f)
    real, intent(inout) :: nh4_mass(nilg,iccp1)    !< ammonium mass for individual 9 pfts + bare (\f$g N/m^2\f$)
    real, intent(inout) :: no3_mass(nilg,iccp1)    !< nitrate mass for individual 9 pfts + bare (\f$g N/m^2\f$)
    real, intent(inout) :: fcancmx(nilg,icc)    !< max. fractional coverages of ctem's 9 pfts.
    real, intent(inout) :: pfcancmx(nilg,icc)   !< previous max. fractional coverages of ctem's 9 pfts.
    real, intent(inout) :: vgbiomas(nilg)       !< grid averaged vegetation biomass, kg c/m2
    real, intent(inout) :: nvgbiomas(nilg)       !< grid averaged vegetation biomass, kg c/m2
    real, intent(inout) :: soilnmas(nilg,iccp2) !< soil organic nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout) :: nlitrmass(nilg,iccp2) !< litter nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout) :: soilcmas(nilg,iccp2,ignd) !< soil c mass in kg c/m2, for the 9 pfts + bare
    real, intent(inout) :: litrmass(nilg,iccp2,ignd) !< litter mass in kg c/m2, for the 9 pfts + bare

    real, intent(inout) :: gavgltms(nilg)       !< grid averaged litter mass including the LUC product pool, kg c/m2
    real, intent(inout) :: gavgscms(nilg)       !< grid averaged soil c mass including the LUC product pool, kg c/m2
    real, intent(inout) :: nlitrmass_g(nilg)
    real, intent(inout) :: soilnmas_g(nilg)
    real, intent(inout) :: nh4_mass_g(nilg)
    real, intent(inout) :: no3_mass_g(nilg)
    real, intent(inout) :: nfcancmx(nilg,icc)   !< next max. fractional coverages of ctem's 9 pfts.
    real, intent(inout) :: fcanmx(nilg,icp1)    !< fractional coverages of class 4 pfts (these are found based on new fcancmxs)
    real, intent(inout) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    real, intent(out) :: lucemcom(nilg) !< luc related carbon emission losses from combustion u-mol co2/m2.sec
    real, intent(out) :: lucltrin(nilg) !< luc related input to litter pool, u-mol co2/m2.sec
    real, intent(out) :: lucsocin(nilg) !< luc related input to soil carbon pool, u-mol co2/m2.sec
    real, intent(out) :: lucemcomn(nilg) !< luc related carbon emission losses from combustion of N (g N m^-2 day^-1)
    real, intent(out) :: lucltrinn(nilg) !< luc related input to litter pool of N (g N m^-2 day^-1)
    real, intent(out) :: lucsocinn(nilg) !< luc related input to soil nitrogen pool of N (g N m^-2 day^-1)

    ! Local variables
    real :: delfrac(nilg,icc)    !<
    real :: abvgmass(nilg,icc)   !< above-ground biomass
    real :: abvgnmass(nilg,icc)  !< above-ground N biomass
    real :: combustc(nilg,icc)   !< total carbon from deforestation- combustion
    real :: combustn(nilg,icc)   !< total nitrogen from deforestation- combustion
    real :: paperc(nilg,icc)     !< total carbon from deforestation- paper
    real :: papern(nilg,icc)     !< total nitrogen from deforestation- paper
    real :: furnturc(nilg,icc)   !< total carbon from deforestation- furniture
    real :: furnturn(nilg,icc)   !< total nitrogen from deforestation- furniture
    real :: incrlitr             !<
    real :: incrnlitr            !<
    real :: incrsolc             !<
    real :: incrsoln             !<
    real :: incrnh4              !<
    real :: incrno3              !<
    real :: chopedbm(nilg)       !< chopped off biomass
    real :: chopednbm(nilg)      !< chopped off N biomass
    real :: compdelfrac(nilg,icc)!< with competition on, this is the change in pft frac per timestep
    real :: fcancmy(nilg,icc)    !<
    integer :: i, j, k, m, n, q
    real :: redubmas1      !<
    real :: redunbmas1     !<
    real :: term           !<
    real :: barefrac(nilg) !< initialize bare fraction to 1.0
    real :: grsumcom(nilg) !< grid sum of combustion carbon for all pfts that are chopped
    real :: grsumncom(nilg)!< grid sum of combustion nitrogen for all pfts that are chopped
    real :: grsumpap(nilg) !< grid sum of paper carbon for all pfts that are chopped
    real :: grsumnpap(nilg)!< grid sum of paper nitrogen for all pfts that are chopped
    real :: grsumfur(nilg) !< grid sum of furniture carbon for all pfts that are chopped
    real :: grsumnfur(nilg)!< grid sum of furniture nitrogen for all pfts that are chopped
    real :: grsumson(nilg) !< grid sum of soil nitrogen for all pfts that are chopped
    real :: grsumnh4(nilg) !< grid sum of nh4 nitrogen for all pfts that are chopped
    real :: grsumno3(nilg) !< grid sum of no3 nitrogen for all pfts that are chopped
    real :: grdennlit(nilg)!< grid averaged densities for litter nitrogen
    real :: grdenson(nilg) !< grid averaged densities for soil nitrogen
    real :: grsumlit(nilg,ignd) !< grid sum of litter carbon for all pfts that are chopped
    real :: grsumnlit(nilg)!< grid sum of litter nitrogen for all pfts that are chopped
    real :: grsumsoc(nilg,ignd) !< grid sum of soil c carbon for all pfts that are chopped
    real :: grdenlit(nilg,ignd) !< grid averaged densities for litter carbon
    real :: grdensoc(nilg,ignd) !< grid averaged densities for soil c carbon
    real :: pbarefra(nilg) !< initialize previous years's bare fraction to 1.0
    real :: grdencom(nilg) !< grid averaged densities for combustion carbon
    real :: grdenncom(nilg)!< grid averaged densities for combustion nitrogen
    real :: grdenpap(nilg) !< grid averaged densities for paper carbon
    real :: grdennpap(nilg)!< grid averaged densities for paper nitrogen
    real :: grdenfur(nilg) !< grid averaged densities for furniture carbon
    real :: grdennfur(nilg)!< grid averaged densities for furniture nitrogen
    real :: grdennh4(nilg) !< grid averaged densities for nh4 nitrogen
    real :: grdenno3(nilg) !< grid averaged densities for no3 nitrogen

    integer :: fraciord(nilg,icc) !< fractional coverage increase or decrease increase +1, decrease -1
    integer :: treatind(nilg,icc) !< treatment index for combust, paper, & furniture
    integer :: bareiord(nilg)    !< bare fraction increases or decreases
    integer :: lrgstpft(1)       !<
    logical  luctkplc(nilg)   !<
    real :: totcmass(nilg) !< total c mass (live+dead)
    real :: totnmass(nilg) !< total N mass (live+dead)
    real :: totlmass(nilg) !< total c mass (live)
    real :: totlnmass(nilg) !< total N mass (live)
    real :: totdmas1(nilg) !< total c mass (dead) litter
    real :: totdnmas1(nilg)  !< total N mass (dead) litter
    real :: totdnh4mas(nilg) !< total N mass (dead) NH4
    real :: totdno3mas(nilg) !< total N mass (dead) NO3
    real :: ntotcmas(nilg) !< total c mass (live+dead) after luc treatment
    real :: ntotnmas(nilg) !< total c mass (live+dead) after luc treatment
    real :: ntotlmas(nilg) !< total c mass (live) after luc treatment
    real :: ntotlnmas(nilg) !< total c mass (live) after luc treatment
    real :: ntotdms1(nilg) !< total c mass (dead) litter after luc treatment
    real :: ntotdnms1(nilg) !< total N mass (dead) litter after luc treatment
    real :: pvgbioms(nilg) !< Vegetation biomass on entering subroutine
    real :: pvgbiomsn(nilg) !< Vegetation biomass on entering subroutine
    real :: pgavltms(nilg) !< Litter mass on entering subroutine
    real :: pgavscms(nilg) !< Soil C mass on entering subroutine
    real :: pluclitpool(nilg) !< LUC paper pool on entering subroutine
    real :: plucscpool(nilg) !< LUC furniture pool on entering subroutine
    real :: pgavltmsn(nilg)
    real :: pgavscmsn(nilg)
    real :: pgavnh4(nilg)
    real :: pgavno3(nilg)
    real :: pluclitpooln(nilg)
    real :: plucscpooln(nilg)
    real :: redubmas2      !<
    real :: redunbmas2     !<
    real :: totdmas2(nilg) !< total c mass (dead) soil c
    real :: totdnmas2(nilg)   !< total N mass (dead) soil N
    real :: totdnh4mas2(nilg) !< total N mass (dead) nh4
    real :: totdno3mas2(nilg) !< total N mass (dead) no3
    real :: ntotdms2(nilg) !< total c mass (dead) soil c after luc treatment
    real :: ntotdnms2(nilg)   !< total N mass (dead) soil N after luc treatment
    real :: ntotdnh4ms(nilg)  !<
    real :: ntotdno3ms(nilg)  !<
    real :: pftarrays(nilg,icc - numcrops) !<
    integer :: indexpos(nilg,icc - numcrops) !<

    real     :: diff !<
    real     :: sum1 !<
    real     :: sum2 !<
    real     :: rdif !<

    associate( &
    Ncycle_on => c_switch%Ncycle_on                     & !< logical: true if nitrogen cycle is on. 
    )

    ! -------------------------------------

    !> Find/use provided current and previous day's fractional coverage
    !! if competition is on, we will adjust these later.
    if (interpol) then ! perform interpolation
      do j = 1,icc
        do i = il1,il2
          if (PFTCompetition .and. .not. crop(j)) then
            nfcancmx(i,j) = yesfrac(i,j)
            pfcancmx(i,j) = yesfrac(i,j)
          end if
          delfrac(i,j) = nfcancmx(i,j) - pfcancmx(i,j) ! change in fraction

          if (leapnow) then
            delfrac(i,j) = delfrac(i,j)/366.0
          else
            delfrac(i,j) = delfrac(i,j)/365.0
          end if

          fcancmx(i,j) = pfcancmx(i,j) + (real(iday) * delfrac(i,j)) !  current day
          fcancmy(i,j) = pfcancmx(i,j) + (real(iday - 1) * delfrac(i,j)) ! previous day

          if (fcancmx(i,j) < 0.0 .and. abs(fcancmx(i,j)) < 1.0e-05) then
            fcancmx(i,j) = 0.0
          else if (fcancmx(i,j) < 0.0 .and. abs(fcancmx(i,j)) >= 1.0e-05) then
            write(6, * )'fcancmx(',i,',',j,') = ',fcancmx(i,j)
            write(6, * )'fractional coverage cannot be negative'
            call errorHandler('luc', - 1)
          end if

          if (fcancmy(i,j) < 0.0 .and. abs(fcancmy(i,j)) < 1.0e-05) then
            fcancmy(i,j) = 0.0
          else if (fcancmy(i,j) < 0.0 .and. &
            abs(fcancmy(i,j)) >= 1.0e-05) then
            write(6, * )'fcancmy(',i,',',j,') = ',fcancmy(i,j)
            write(6, * )'fractional coverage cannot be negative'
            call errorHandler('luc', - 2)
          end if

        end do ! loop 111
      end do ! loop 110
    else ! use provided values but still check they are not negative
      do j = 1,icc
        do i = il1,il2
          fcancmx(i,j) = todfrac(i,j)
          fcancmy(i,j) = yesfrac(i,j)

          if (fcancmx(i,j) < 0.0 .and. abs(fcancmx(i,j)) < 1.0e-05) then
            fcancmx(i,j) = 0.0
          else if (fcancmx(i,j) < 0.0 .and. &
       abs(fcancmx(i,j)) >= 1.0e-05) then
            write(6, * )'fcancmx(',i,',',j,') = ',fcancmx(i,j)
            write(6, * )'fractional coverage cannot be negative'
            call errorHandler('luc', - 3)
          end if

          if (fcancmy(i,j) < 0.0 .and. abs(fcancmy(i,j)) < 1.0e-05) then
            fcancmy(i,j) = 0.0
          else if (fcancmy(i,j) < 0.0 .and. &
       abs(fcancmy(i,j)) >= 1.0e-05) then
            write(6, * )'fcancmy(',i,',',j,') = ',fcancmy(i,j)
            write(6, * )'fractional coverage cannot be negative'
            call errorHandler('luc', - 4)
          end if

        end do ! loop 116
      end do ! loop 115
    end if
    !>
    !> If competition is on, we need to adjust the other fractions for the increase/decrease
    !> in cropland as only the crops areas is now specified.
    if (PFTCompetition) then

      call adjustFracsComp(il1,il2,nilg,iday,pfcancmx,yesfrac,delfrac,compdelfrac)

      do j = 1,icc
        do i = 1,il1,il2

          if (.not. crop(j)) then
            fcancmx(i,j) = yesfrac(i,j) - compdelfrac(i,j) !  current day
            fcancmy(i,j) = yesfrac(i,j) ! previous day
          end if
        end do

      end do
    end if
    !>
    !! check if this year's fractional coverages have changed or not for any pft
    !!
    luctkplc(:) = .false.  ! did land use change take place for any pft in this grid cell
    do j = 1,icc
      do i = il1,il2
        if ( (abs(fcancmx(i,j) - fcancmy(i,j))) > zero) then
          luctkplc(i) = .true. ! yes,luc did take place in this grid cell
        end if
      end do ! loop 250
    end do ! loop 200

    ! only perform the rest of the subroutine if any luc is actually taking
    ! place,otherwise exit.
    do i = il1,il2
      if (luctkplc(i)) then

        ! initialization
        do j = 1,ican
          fcanmx(i,j) = 0.0 ! fractional coverage of class' pfts
        end do ! loop 260

        do j = 1,icc
          fraciord(i,j) = 0   ! fractional coverage increase or decrease
          !                           ! increase +1, decrease -1
          abvgmass(i,j) = 0.0 ! above-ground biomass
          abvgnmass(i,j) = 0.0  ! above-ground N biomass
          treatind(i,j) = 0   ! treatment index for combust, paper, & furniture
          combustc(i,j) = 0.0 ! total carbon from deforestation- combustion
          combustn(i,j) = 0.0 ! total nitrogen from deforestation- combustion
          paperc(i,j) = 0.0   ! total carbon from deforestation- paper
          papern(i,j)   = 0.0 ! total nitrogen from deforestation- paper
          furnturc(i,j) = 0.0 ! total carbon from deforestation- furniture
          furnturn(i,j) = 0.0 ! total nitrogen from deforestation- furniture
        end do ! loop 270

        pvgbioms(i) = vgbiomas(i)  ! store grid average quantities in
        pvgbiomsn(i) = nvgbiomas(i)  ! store grid average quantities in
        pgavltms(i) = gavgltms(i)  ! temporary arrays
        pgavscms(i) = gavgscms(i)
        pluclitpool(i) = litrmass(i,iccp2,1)  ! LUC product pools are stored in first layer.
        plucscpool(i) = soilcmas(i,iccp2,1) ! LUC product pools are stored in first layer.
        pgavltmsn(i) = nlitrmass_g(i)  ! temporary arrays
        pgavscmsn(i) = soilnmas_g(i)
        pgavnh4(i) = nh4_mass_g(i)
        pgavno3(i) = no3_mass_g(i)
        pluclitpooln(i) = nlitrmass(i,iccp2)  ! LUC product pools are stored in first layer.
        plucscpooln(i) = soilnmas(i,iccp2) ! LUC product pools are stored in first layer.

        vgbiomas(i) = 0.0
        nvgbiomas(i) = 0.0
        gavgltms(i) = 0.0
        gavgscms(i) = 0.0
        nlitrmass_g(i) = 0.0
        soilnmas_g(i) = 0.0
        nh4_mass_g(i) = 0.0
        no3_mass_g(i) = 0.0

        grsumcom(i) = 0.0          ! grid sum of combustion carbon for all
        grsumncom(i) = 0.0         ! grid sum of combustion nitrogen for all
        ! pfts that are chopped
        grsumpap(i) = 0.0          ! similarly for paper,
        grsumnpap(i) = 0.0
        grsumfur(i) = 0.0          ! furniture,
        grsumnfur(i) = 0.0
        grsumnlit(i) = 0.0
        grsumson(i) = 0.0
        grsumnh4(i) = 0.0
        grsumno3(i) = 0.0
        do k = 1,ignd
            grsumlit(i,k)=0.0          ! grid sum of combustion carbon for all chopped PFT's litter, and
            grsumsoc(i,k)=0.0          ! soil c
            grdenlit(i,k)=0.0          ! grid averaged densities for litter, and
            grdensoc(i,k)=0.0          ! soil c
        end do
        grdennlit(i) = 0.0
        grdenson(i)  = 0.0
        grdennh4(i)  = 0.0
        grdenno3(i)  = 0.0

        grdencom(i) = 0.0          ! grid averaged densities for combustion carbon,
        grdenncom(i) = 0.0
        grdenpap(i) = 0.0          ! paper,
        grdennpap(i) = 0.0
        grdenfur(i) = 0.0          ! furniture,
        grdennfur(i) = 0.0

        totcmass(i) = 0.0          ! total c mass (live+dead)
        totnmass(i) = 0.0          ! total N mass (live+dead)
        totlmass(i) = 0.0          ! total c mass (live)
        totlnmass(i) = 0.0         ! total N mass (live)
        totdmas1(i) = 0.0          ! total c mass (dead) litter
        totdnmas1(i) = 0.0         ! total N mass (dead) litter
        totdmas2(i) = 0.0          ! total c mass (dead) soil c
        totdnmas2(i)  = 0.0        ! total N mass (dead) soil N
        totdnh4mas(i) = 0.0
        totdno3mas(i) = 0.0

        barefrac(i) = 1.0
        pbarefra(i) = 1.0


        ntotcmas(i) = 0.0          ! total c mass (live+dead)
        ntotnmas(i) = 0.0          ! total N mass (live+dead)
        ntotlmas(i) = 0.0          ! total c mass (live)
        ntotlnmas(i) = 0.0         ! total N mass (live)
        ntotdms1(i) = 0.0          ! total c mass (dead) litter
        ntotdnms1(i) = 0.0         ! total N mass (dead) litter
        ntotdms2(i) = 0.0          ! total c mass (dead) soil c
        ntotdnms2(i)   = 0.0       ! total N mass (dead) soil N
        ntotdnh4ms(i) = 0.0
        ntotdno3ms(i) = 0.0
        ntotdnh4ms(i) = 0.0
        ntotdno3ms(i) = 0.0

        lucemcom(i) = 0.0
        lucltrin(i) = 0.0
        lucsocin(i) = 0.0
        lucemcomn(i) = 0.0
        lucltrinn(i) = 0.0
        lucsocinn(i) = 0.0

        bareiord(i) = 0
        chopedbm(i) = 0.0
        chopednbm(i) = 0.0

        !> initialization ends

        !> -------------------------------------------------------------------

        !> if land use change has taken place then get fcanmxs for use by class based on the new fcancmxs

        do j = 1,ican
          do m = reindexPFTs(j,1),reindexPFTs(j,2)
            fcanmx(i,j) = fcanmx(i,j) + fcancmx(i,m)
            barefrac(i) = barefrac(i) - fcancmx(i,m)
          end do ! loop 301
        end do ! loop 300
        !>
        !> check if the interpol didn't mess up the barefrac. if so, take the
        !! extra amount from the pft with the largest area.
        !! but you can't take it from crops !
        pftarrays = 0.

        if (barefrac(i) < 0.0) then ! PFTCompetition only needs minbare but it checks later.

          k = 1
          do j = 1,icc
            if (.not. crop(j)) then
              indexpos(i,k) = j
              pftarrays(i,k) = fcancmx(i,j)
              k = k + 1
            end if
          end do

          lrgstpft = maxloc(pftarrays(i,:))
          j = indexpos(i,lrgstpft(1))

          if (PFTCompetition) then
            fcancmx(i,j) = fcancmx(i,j) + barefrac(i) - minbare
            barefrac(i) = minbare
          else
            fcancmx(i,j) = fcancmx(i,j) + barefrac(i)
            barefrac(i) = 0.0
          end if
        end if

        !> find previous day's bare fraction using previous day's fcancmxs
        do j = 1,icc
          pbarefra(i) = pbarefra(i) - fcancmy(i,j)
        end do ! loop 310
        !>
        !> check if the interpol didn't mess up the pbarefra. if so, take the
        !! extra amount from the pft with the largest area.
        !! but you can't take it from crops !
        if (pbarefra(i) < 0.0) then
          k = 1
          do j = 1,icc
            if (.not. crop(j)) then
              indexpos(i,k) = j
              pftarrays(i,k) = fcancmy(i,j)
              k = k + 1
            end if
          end do

          lrgstpft = maxloc(pftarrays(i,:))
          j = indexpos(i,lrgstpft(1))

          if (PFTCompetition) then
            fcancmy(i,j) = fcancmy(i,j) + pbarefra(i) - minbare
            pbarefra(i) = minbare
          else
            fcancmy(i,j) = fcancmy(i,j) + pbarefra(i)
            pbarefra(i) = 0.0
          end if

        end if
        !>
        !> based on sizes of 3 live pools and 2 dead pools we estimate the total amount of c in each grid cell.
        !>
        do j = 1,icc
          totlmass(i) = totlmass(i) &
                        + (fcancmy(i,j) * (gleafmas(i,j) + bleafmas(i,j) &
                        + stemmass(i,j) + rootmass(i,j)) * grclarea(i) &
                        * km2tom2)
          if (Ncycle_on) then
             totlnmass(i) = totlnmass(i) &
                            + (fcancmy(i,j) * (ngleafmas(i,j) + nbleafmas(i,j) &
                            + nstemmass(i,j) + nrootmass(i,j)) * grclarea(i) &
                            * km2tom2)
          end if !for Ncycle_on
        end do ! loop 320

        do j = 1,iccp2

          if (j < iccp1) then
            term = fcancmy(i,j)
          else if (j == iccp1) then
            term = pbarefra(i)
          else if (j == iccp2) then ! LUC product pools cover entire cell.
            term = 1.0
          end if

          do k = 1,ignd
            totdmas1(i)=totdmas1(i)+ &
                       (term*litrmass(i,j,k)*grclarea(i)*km2tom2)
            totdmas2(i)=totdmas2(i)+ &
                       (term*soilcmas(i,j,k)*grclarea(i)*km2tom2)
          end do           

          if (Ncycle_on) then
              totdnmas1(i) = totdnmas1(i) + &
                             (term * nlitrmass(i,j) * grclarea(i) * km2tom2)
              totdnmas2(i) = totdnmas2(i) + &
                             (term * soilnmas(i,j) * grclarea(i) * km2tom2)
              if(j <= iccp1) then
              totdnh4mas(i)= totdnh4mas(i) + &
                             (term * nh4_mass(i,j) * grclarea(i) * km2tom2)
              totdno3mas(i)= totdno3mas(i) + &
                             (term * no3_mass(i,j) * grclarea(i) * km2tom2)
              end if
          end if !for Ncycle_on

        end do ! loop 340

        totcmass(i) = totlmass(i) + totdmas1(i) + totdmas2(i)
        if (Ncycle_on) then
        totnmass(i) = totlnmass(i) + totdnmas1(i) + totdnmas2(i) + &
                      totdnh4mas(i) + totdno3mas(i)
        end if
        ! bare fractions cannot be negative
        if (pbarefra(i) < 0.0 .and. abs(pbarefra(i)) < 1.0e-05) then
          pbarefra(i) = 0.0
        else if (pbarefra(i) < 0.0 .and. abs(pbarefra(i)) >= 1.0e-05) then
          write(6, * )'bare fractions cannot be negative'
          write(6, * )'prev. bare fraction(',i,')  = ',pbarefra(i)
          call errorHandler('luc', - 5)
        end if

        if (barefrac(i) < 0.0 .and. abs(barefrac(i)) < 1.0e-05) then
          write( * , * )'setting bare to zero',barefrac(i)
          barefrac(i) = 0.0
        else if (barefrac(i) < 0.0 .and. abs(barefrac(i)) >= 1.0e-05) then
          write(6, * )'bare fractions cannot be negative'
          write(6, * )'bare fraction(',i,')  = ',barefrac(i)
          call errorHandler('luc', - 6)
        end if

        !> Find above ground biomass and treatment index for combust, paper,
        !> and furniture
        do j = 1,ican
          do m = reindexPFTs(j,1),reindexPFTs(j,2)
            abvgmass(i,m) = gleafmas(i,m) + bleafmas(i,m) + stemmass(i,m)
            if (Ncycle_on) abvgnmass(i,m) = ngleafmas(i,m) + nbleafmas(i,m) + nstemmass(i,m)
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
            case ('BdlSh')  ! FLAG - should shrubs only be bush or shrubs? JM FIXME
              if (abvgmass(i,m) >= bmasthrs(1)) then ! forest
                treatind(i,m) = 1
              else if (abvgmass(i,m) <= bmasthrs(2)) then ! bush
                treatind(i,m) = 3
              else  ! shrubland
                treatind(i,m) = 2
              end if
            case default
              print * ,'Unknown CLASS PFT in LUC ',classpfts(j)
              call errorHandler('luc', - 7)
            end select
          end do ! loop 510
        end do ! loop 500

        !> check if a pft's fractional cover is increasing or decreasing

        do j = 1,icc
          if ( (fcancmx(i,j) > fcancmy(i,j)) .and. &
              (abs(fcancmy(i,j) - fcancmx(i,j)) > zero) ) then
            fraciord(i,j) = 1  ! increasing
          else if ( (fcancmx(i,j) < fcancmy(i,j)) .and. &
                  (abs(fcancmy(i,j) - fcancmx(i,j)) > zero) ) then
            fraciord(i,j) = - 1 ! decreasing
          end if
        end do ! loop 550

        !> check if bare fraction increases of decreases

        if ( (barefrac(i) > pbarefra(i)) .and. &
            (abs(pbarefra(i) - barefrac(i)) > zero) ) then
          bareiord(i) = 1  ! increasing
        else if ( (barefrac(i) < pbarefra(i)) .and. &
                 (abs(pbarefra(i) - barefrac(i)) > zero) ) then
          bareiord(i) = - 1 ! decreasing
        end if

        !>
        !> if the fractional coverage of pfts increases then spread their live & dead biomass uniformly
        !! over the new fraction. this effectively reduces their per m2 c density.
        !!
        do j = 1,icc
          if (fraciord(i,j) == 1) then
            term = fcancmy(i,j)/fcancmx(i,j)
            gleafmas_ns(i,j) = gleafmas_ns(i,j) * term
            gleafmas_s(i,j)  = gleafmas_s(i,j)  * term
            bleafmas(i,j) = bleafmas(i,j) * term
            stemmass_ns(i,j) = stemmass_ns(i,j) *term
            stemmass_s(i,j)  = stemmass_s(i,j)  *term
            rootmass_ns(i,j) = rootmass_ns(i,j) *term
            rootmass_s(i,j)  = rootmass_s(i,j)  *term

            ! refind the total pool sizes:
            gleafmas(i,j) = gleafmas_ns(i,j) + gleafmas_s(i,j)
            stemmass(i,j) = stemmass_ns(i,j) + stemmass_s(i,j)
            rootmass(i,j) = rootmass_ns(i,j) + rootmass_s(i,j)

            do k = 1,ignd
              litrmass(i,j,k)=litrmass(i,j,k)*term
              soilcmas(i,j,k)=soilcmas(i,j,k)*term
            end do 
            if (Ncycle_on) then
               ngleafmas_ns(i,j)= ngleafmas_ns(i,j)* term
               ngleafmas_s(i,j) = ngleafmas_s(i,j) * term
               nbleafmas(i,j)   = nbleafmas(i,j)   * term
               nstemmass_ns(i,j)= nstemmass_ns(i,j)* term
               nstemmass_s(i,j) = nstemmass_s(i,j) * term
               nrootmass_ns(i,j)= nrootmass_ns(i,j)* term
               nrootmass_s(i,j) = nrootmass_s(i,j) * term
               nlitrmass(i,j)   = nlitrmass(i,j)   * term
               soilnmas(i,j)    = soilnmas(i,j)    * term
               nh4_mass(i,j)    = nh4_mass(i,j)    * term
               no3_mass(i,j)    = no3_mass(i,j)    * term

               ! refind the total pool sizes:
               ngleafmas(i,j) = ngleafmas_ns(i,j) + ngleafmas_s(i,j)
               nstemmass(i,j) = nstemmass_ns(i,j) + nstemmass_s(i,j)
               nrootmass(i,j) = nrootmass_ns(i,j) + nrootmass_s(i,j)
               
            end if !for Ncycle_on
            if (useTracer > 0) then ! Now same operation for tracer
              tracerGLeafMass(i,j) = tracerGLeafMass(i,j) * term
              tracerBLeafMass(i,j) = tracerBLeafMass(i,j) * term
              tracerStemMass(i,j) = tracerStemMass(i,j) * term
              tracerRootMass(i,j) = tracerRootMass(i,j) * term
              do k = 1,ignd
                tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) * term
                tracerSoilCMass(i,j,k) = tracerSoilCMass(i,j,k) * term
              end do
            end if

          end if
        end do ! loop 570
        !>
        !> if bare fraction increases then spread its litter and soil c/n
        !> uniformly over the increased fraction
        !>
        if (bareiord(i) == 1) then
          term = pbarefra(i)/barefrac(i)
          do k = 1,ignd
            litrmass(i,iccp1,k)=litrmass(i,iccp1,k)*term
            soilcmas(i,iccp1,k)=soilcmas(i,iccp1,k)*term
            if (useTracer > 0) then ! Now same operation for tracer
              tracerLitrMass(i,iccp1,k) = tracerLitrMass(i,iccp1,k) * term
              tracerSoilCMass(i,iccp1,k) = tracerSoilCMass(i,iccp1,k) * term
            end if
          end do ! loop 582
          ! COMBAK PERLAY
          if (Ncycle_on) then
             nlitrmass(i,iccp1)= nlitrmass(i,iccp1)* term
             soilnmas(i,iccp1) = soilnmas(i,iccp1) * term
             nh4_mass(i,iccp1) = nh4_mass(i,iccp1) * term
             no3_mass(i,iccp1) = no3_mass(i,iccp1) * term
          end if !for Ncycle_on
        end if
        !>
        !> if any of the pfts fractional coverage decreases, then we chop the
        !! aboveground biomass and treat it according to our rules (burn it,
        !! and convert it into paper and furniture). the below ground live
        !! biomass and litter of this pfts gets assimilated into litter of
        !! all pfts (uniformly spread over the whole grid cell), and soil c
        !! from the chopped off fraction of this pft, gets assimilated into
        !! soil c of all existing pfts as well.
        !!
        do j = 1,ican
          do m = reindexPFTs(j,1),reindexPFTs(j,2)
            if (fraciord(i,m) == - 1) then

              !> chop off above ground biomass
              redubmas1 = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i) &
                          * abvgmass(i,m) * km2tom2

              if (redubmas1 < 0.0) then
                write(6, * )'redubmas1 less than zero'
                write(6, * )'fcancmy - fcancmx = ', &
                          fcancmy(i,m) - fcancmx(i,m)
                write(6, * )'grid cell = ',i,' pft = ',m
                call errorHandler('luc', - 8)
              end if
              if (Ncycle_on) then
                 redunbmas1 = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i) &
                              * abvgnmass(i,m) * km2tom2
 
                 if (redunbmas1 < 0.0) then
                   write(6, * )'redunbmas1 less than zero'
                   write(6, * )'fcancmy-fcancmx = ',  &
                               fcancmy(i,m) - fcancmx(i,m)
                   write(6, * )'grid cell = ',i,' pft = ',m
                   !call errorHandler('luc', - 9)
                 end if
              end if !for Ncycle_on

              !> rootmass needs to be chopped as well and all of it goes to the litter/paper pool
              redubmas2 = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i) &
                          * rootmass(i,m) * km2tom2
              if (Ncycle_on) then
                 redunbmas2 = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i)  &
                              * nrootmass(i,m) * km2tom2
              end if !for Ncycle_on

              !> keep adding chopped off biomass for each pft to get the total for a grid cell for diagnostics
              chopedbm(i) = chopedbm(i) + redubmas1 + redubmas2
              if (Ncycle_on) then
                 chopednbm(i) = chopednbm(i) + redunbmas1 + redunbmas2
              end if !for Ncycle_on

              !> find what's burnt,and what's converted to paper & furniture
              combustc(i,m) = combust(treatind(i,m)) * redubmas1
              paperc(i,m) = paper(treatind(i,m)) * redubmas1 + redubmas2
              furnturc(i,m) = furniture(treatind(i,m)) * redubmas1
              if (Ncycle_on) then
                 combustn(i,m) = combust(treatind(i,m))   * redunbmas1
                 papern(i,m)   = paper(treatind(i,m))     * redunbmas1 + redunbmas2
                 furnturn(i,m) = furniture(treatind(i,m)) * redunbmas1
              end if !for Ncycle_on

              !> keep adding all this for a given grid cell
              grsumcom(i) = grsumcom(i) + combustc(i,m)
              grsumpap(i) = grsumpap(i) + paperc(i,m)
              grsumfur(i) = grsumfur(i) + furnturc(i,m)
              if (Ncycle_on) then
                 grsumncom(i) = grsumncom(i) + combustn(i,m)
                 grsumnpap(i) = grsumnpap(i) + papern(i,m)
                 grsumnfur(i) = grsumnfur(i) + furnturn(i,m)
              endif !for Ncycle_on

              !> litter from the chopped off fraction of the chopped
              !> off pft needs to be assimilated, and so does soil c/n from
              !> the chopped off fraction of the chopped pft
              !>
              do k = 1,ignd
                incrlitr=(fcancmy(i,m)-fcancmx(i,m))*grclarea(i) &
                        *litrmass(i,m,k)*km2tom2

                incrsolc=(fcancmy(i,m)-fcancmx(i,m))*grclarea(i) &
                        *soilcmas(i,m,k)*km2tom2

                grsumlit(i,k)=grsumlit(i,k)+incrlitr
                grsumsoc(i,k)=grsumsoc(i,k)+incrsolc
              end do 
              
              if (Ncycle_on) then
                 incrnlitr = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i) &
                             * nlitrmass(i,m) * km2tom2
                 incrsoln  = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i) &
                             * soilnmas(i,m) * km2tom2
                 incrnh4  = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i) &
                             * nh4_mass(i,m) * km2tom2
                 incrno3  = (fcancmy(i,m) - fcancmx(i,m)) * grclarea(i) &
                             * no3_mass(i,m) * km2tom2
                 grsumnlit(i) = grsumnlit(i) + incrnlitr
                 grsumson(i)  = grsumson(i)  + incrsoln
                 grsumnh4(i)  = grsumnh4(i)  + incrnh4
                 grsumno3(i)  = grsumno3(i)  + incrno3
              end if !for Ncycle_on
            end if

          end do ! loop 610
        end do ! loop 600

        !> if bare fraction decreases then chop off the litter and soil c/n
        !! from the decreased fraction and add it to grsumlit & grsumsoc for
        !! putting in the LUC product pools
        if (bareiord(i) == - 1) then
          do k = 1,ignd
            redubmas1=(pbarefra(i)-barefrac(i))*grclarea(i) &
                      *litrmass(i,iccp1,k)*km2tom2

            redubmas2=(pbarefra(i)-barefrac(i))*grclarea(i) &
                      *soilcmas(i,iccp1,k)*km2tom2

            grsumlit(i,k)=grsumlit(i,k)+redubmas1
            grsumsoc(i,k)=grsumsoc(i,k)+redubmas2
          end do
          
          if (Ncycle_on) then
             redunbmas1 = (pbarefra(i) - barefrac(i)) * grclarea(i) &
                          * nlitrmass(i,iccp1) * km2tom2
             redunbmas2 = (pbarefra(i) - barefrac(i)) * grclarea(i) &
                          * soilnmas(i,iccp1) * km2tom2
             grsumnlit(i) = grsumnlit(i) + redunbmas1
             grsumson(i)  = grsumson(i)  + redunbmas2
             redunbmas1 = (pbarefra(i) - barefrac(i)) * grclarea(i) &
                          * nh4_mass(i,iccp1) * km2tom2
             grsumnh4(i) = grsumnh4(i) + redunbmas1
             redunbmas2 = (pbarefra(i) - barefrac(i)) * grclarea(i) &
                          * no3_mass(i,iccp1) * km2tom2
             grsumno3(i) = grsumno3(i) + redunbmas2
          end if !for Ncycle_on
        end if

        !> calculate if the chopped off biomass equals the sum of grsumcom(i), grsumpap(i) & grsumfur(i)

        sum1=grsumcom(i)+grsumpap(i)+grsumfur(i)
        diff=abs(chopedbm(i)-sum1)
        if (chopedbm(i)==0.) then
          if (diff==0.) then
            rdif=0.              ! relative accuracy
          else
            rdif=1.              ! relative accuracy
          end if
        else
          rdif=diff/chopedbm(i)  ! relative accuracy
        end if

        !if (abs(chopedbm(i) - grsumcom(i) - grsumpap(i) - grsumfur(i)) > &
        !    tolrnce1) then
        if (rdif > rel_tolerance_amount) then
          write(6, * )'at grid cell = ',i
          write(6, * )'chopped biomass does not equals sum of total'
          write(6, * )'luc related emissions'
          write(6, * )'chopedbm(i) = ',chopedbm(i)
          write(6, * )'grsumcom(i) = ',grsumcom(i)
          write(6, * )'grsumpap(i) = ',grsumpap(i)
          write(6, * )'grsumfur(i) = ',grsumfur(i)
          write(6, * )'sum of grsumcom,grsumpap,grsumfur(i) = ', &
            grsumcom(i) + grsumpap(i) + grsumfur(i)
          write(6, * )'ABS. DIFF. = ',diff
          write(6, * )'REL. DIFF. = ',rdif,' and should be less than ',rel_tolerance_amount
          call errorHandler('luc', - 10)
        end if

        if ( Ncycle_on .and. abs(chopednbm(i) - grsumncom(i) - &
                                 grsumnpap(i) - grsumnfur(i)) > tolrnce1 ) then
           write(6, * )'at grid cell = ',i
           write(6, * )'chopped N biomass does not equal sum of total'
           write(6, * )'luc related emissions'
           write(6, * )'chopednbm(i) = ',chopednbm(i)
           write(6, * )'grsumncom(i) = ',grsumncom(i)
           write(6, * )'grsumnpap(i) = ',grsumnpap(i)
           write(6, * )'grsumnfur(i) = ',grsumnfur(i)
           write(6, * )'sum of grsumncom, grsumnpap, grsumnfur(i) = ', &
                        grsumncom(i) + grsumnpap(i) + grsumnfur(i)
          call errorHandler('luc', - 11)
        end if
        !>
        !! spread chopped off stuff uniformly over the litter and soil c/n pools of all existing pfts, including the bare fraction.
        !!
        !! convert the available c/n into density
        !!
        grdencom(i) = grsumcom(i)/(grclarea(i) * km2tom2)
        grdenpap(i) = grsumpap(i)/(grclarea(i) * km2tom2)
        grdenfur(i) = grsumfur(i)/(grclarea(i) * km2tom2)
        do k = 1,ignd
          grdenlit(i,k)=grsumlit(i,k)/(grclarea(i)*km2tom2)
          grdensoc(i,k)=grsumsoc(i,k)/(grclarea(i)*km2tom2)
        end do
        
        if (Ncycle_on) then
           grdenncom(i) = grsumncom(i)/(grclarea(i) * km2tom2)
           grdennpap(i) = grsumnpap(i)/(grclarea(i) * km2tom2)
           grdennfur(i) = grsumnfur(i)/(grclarea(i) * km2tom2)
           grdennlit(i) = grsumnlit(i)/(grclarea(i) * km2tom2)
           grdenson(i)  = grsumson(i)/(grclarea(i) * km2tom2)
           grdennh4(i)  = grsumnh4(i)/(grclarea(i) * km2tom2)
           grdenno3(i)  = grsumno3(i)/(grclarea(i) * km2tom2)
        end if !for Ncycle_on

        !>      Now add the C to the gridcell's LUC pools of litter and soil C.
        !!      The fast decaying dead LUC C (paper) and slow (furniture) are kept in
        !!      the first 'soil' layer and iccp2 position. The litter and soil C
        !!      contributions are added to the normal litter and soil C pools below.
        litrmass(i,iccp2,1)=litrmass(i,iccp2,1)+grdenpap(i)
        soilcmas(i,iccp2,1)=soilcmas(i,iccp2,1)+grdenfur(i)
        if (useTracer > 0) then ! Now same operation for tracer
          tracerLitrMass(i,iccp2,1) = tracerLitrMass(i,iccp2,1) + grdenpap(i)
          tracerSoilCMass(i,iccp2,1) = tracerSoilCMass(i,iccp2,1) + grdenfur(i)
        end if
        if (Ncycle_on) then
          nlitrmass(i,iccp2)=nlitrmass(i,iccp2)+grdennpap(i)
          soilnmas(i,iccp2)=soilnmas(i,iccp2)+grdennfur(i)
        end if !for Ncycle_on

        !     Add any adjusted litter and soilc back their respective pools

        do j = 1,icc
          if (fcancmx(i,j) > zero) then
            litrmass(i,j,:)=litrmass(i,j,:)+grdenlit(i,:)
            soilcmas(i,j,:)=soilcmas(i,j,:)+grdensoc(i,:)
            if (Ncycle_on) then
                nlitrmass(i,j) = nlitrmass(i,j) + grdennlit(i)
                soilnmas(i,j)  = soilnmas(i,j)  + grdenson(i)
                nh4_mass(i,j)  = nh4_mass(i,j)  + grdennh4(i)
                no3_mass(i,j)  = no3_mass(i,j)  + grdenno3(i)
            end if
            if (useTracer > 0) then ! Now same operation for tracer
              tracerLitrMass(i,j,:) = tracerLitrMass(i,j,:) + grdenlit(i,:)
              tracerSoilCMass(i,j,:) = tracerSoilCMass(i,j,:) + grdensoc(i,:)
            end if
          else
            gleafmas(i,j) = 0.0
            gleafmas_ns(i,j) = 0.0
            gleafmas_s(i,j)  = 0.0
            bleafmas(i,j) = 0.0
            stemmass(i,j) = 0.0
            stemmass_ns(i,j) = 0.0
            stemmass_s(i,j)  = 0.0
            rootmass(i,j) = 0.0
            rootmass_ns(i,j) = 0.0
            rootmass_s(i,j)  = 0.0
            litrmass(i,j,:)=0.0  ! set all soil layers to 0
            soilcmas(i,j,:)=0.0
            if (Ncycle_on) then
              ngleafmas(i,j)   = 0.0
              ngleafmas_ns(i,j)= 0.0
              ngleafmas_s(i,j) = 0.0
              nbleafmas(i,j)   = 0.0
              nstemmass(i,j)   = 0.0
              nstemmass_ns(i,j)= 0.0
              nstemmass_s(i,j) = 0.0
              nrootmass(i,j)   = 0.0
              nrootmass_ns(i,j)= 0.0
              nrootmass_s(i,j) = 0.0
              nlitrmass(i,j)   = 0.0
              soilnmas(i,j)    = 0.0
              nh4_mass(i,j)    = 0.0
              no3_mass(i,j)    = 0.0
            end if
            if (useTracer > 0) then ! Now same operation for tracer
              tracerGLeafMass(i,j) = 0.
              tracerBLeafMass(i,j) = 0.
              tracerStemMass(i,j) = 0.
              tracerRootMass(i,j) = 0.
              tracerLitrMass(i,j,:) = 0.
              tracerSoilCMass(i,j,:) = 0.
            end if
          end if
        end do ! loop 650

        ! COMBAK PERLAY
        if (barefrac(i) > zero) then
          litrmass(i,iccp1,:)=litrmass(i,iccp1,:)+grdenlit(i,:)
          soilcmas(i,iccp1,:)=soilcmas(i,iccp1,:)+grdensoc(i,:)
          if (Ncycle_on) then
            nlitrmass(i,iccp1) = nlitrmass(i,iccp1) + grdennlit(i)
            soilnmas(i,iccp1)  = soilnmas(i,iccp1)  + grdenson(i)
            nh4_mass(i,iccp1)  = nh4_mass(i,iccp1)  + grdennh4(i)
            no3_mass(i,iccp1)  = no3_mass(i,iccp1)  + grdenno3(i)
          end if
          if (useTracer > 0) then ! Now same operation for tracer
            tracerLitrMass(i,iccp1,:) = tracerLitrMass(i,iccp1,:) + grdenlit(i,:)
            tracerSoilCMass(i,iccp1,:) = tracerSoilCMass(i,iccp1,:) + grdensoc(i,:)
          end if
        else
          litrmass(i,iccp1,:)=0.0 ! set all soil layers to 0
          soilcmas(i,iccp1,:)=0.0
          if (Ncycle_on) then
            nlitrmass(i,iccp1) = 0.0
            soilnmas(i,iccp1)  = 0.0
            nh4_mass(i,iccp1)  = 0.0
            no3_mass(i,iccp1)  = 0.0
          end if
          if (useTracer > 0) then ! Now same operation for tracer
            tracerLitrMass(i,iccp1,:) = 0.0 ! set all soil layers to 0
            tracerSoilCMass(i,iccp1,:) = 0.0
          end if
        end if

        !> the combusted c is used to find the c flux that we can release into the atmosphere.

        lucemcom(i) = grdencom(i) ! this is flux in kg C/m2.day that will be emitted
        lucltrin(i) = grdenpap(i) ! flux in kg C/m2.day
        lucsocin(i) = grdenfur(i) ! flux in kg C/m2.day
        lucemcomn(i) = grdenncom(i)
        lucltrinn(i) = grdennpap(i)
        lucsocinn(i) = grdennfur(i)

        ! If the entire grid cell is bare ground, add lucemcom to the litter C pool:
        if (( 1 - barefrac(i) ) < zero ) then
          do j = 1,icc
            if (fcancmx(i,j) > zero) then
              litrmass(i,j,1)=litrmass(i,j,1)+lucemcom(i)
              if (Ncycle_on) nlitrmass(i,j) = nlitrmass(i,j) + lucemcomn(i)
            end if
          end do
          lucemcom(i) = 0.0
          lucemcomn(i) = 0.0
        end if

        !> convert all land use change fluxes to u-mol CO2-C/m2.sec
        lucemcom(i) = lucemcom(i) * 963.62
        lucltrin(i) = lucltrin(i) * 963.62
        lucsocin(i) = lucsocin(i) * 963.62

        !> and finally we see if the total amount of carbon and nitrogen are conserved

        do j = 1,icc
          ntotlmas(i) = ntotlmas(i) + &
                        (fcancmx(i,j) * (gleafmas(i,j) + bleafmas(i,j) + &
                        stemmass(i,j) + rootmass(i,j)) * grclarea(i) &
                        * km2tom2)
          if (Ncycle_on) then
              ntotlnmas(i) = ntotlnmas(i)+ &
                             (fcancmx(i,j) * (ngleafmas(i,j) + nbleafmas(i,j) + &
                             nstemmass(i,j) + nrootmass(i,j)) * grclarea(i) &
                             * km2tom2)
          end if
        end do ! loop 700

        do j = 1,iccp2
          if (j < iccp1) then
            term = fcancmx(i,j)
          else if (j == iccp1) then
            term = barefrac(i)
          else if (j == iccp2) then ! assumed to be over entire area for LUC product pools
            term = 1.0
          end if
          do k = 1,ignd
            ntotdms1(i)=ntotdms1(i)+ &
                       (term*litrmass(i,j,k)*grclarea(i)*km2tom2)
            ntotdms2(i)=ntotdms2(i)+ &
                       (term*soilcmas(i,j,k)*grclarea(i)*km2tom2)
          end do
          if (Ncycle_on) then
            ntotdnms1(i) = ntotdnms1(i) + &
                           (term * nlitrmass(i,j) * grclarea(i) * km2tom2)
            ntotdnms2(i) = ntotdnms2(i) + &
                           (term * soilnmas(i,j) * grclarea(i) * km2tom2)
            if(j <= iccp1) then
            ntotdnh4ms(i) = ntotdnh4ms(i) + &
                            (term * nh4_mass(i,j) * grclarea(i) * km2tom2)
            ntotdno3ms(i) = ntotdno3ms(i) + &
                            (term * no3_mass(i,j) * grclarea(i) * km2tom2)
            end if
          end if
        end do ! loop 710

        ntotcmas(i) = ntotlmas(i) + ntotdms1(i) + ntotdms2(i)
        if (Ncycle_on) then
           ntotnmas(i) = ntotlnmas(i) + ntotdnms1(i) + ntotdnms2(i) &
                         + ntotdnh4ms(i) + ntotdno3ms(i)
        end if
        !>
        !> total litter mass (before + input from chopped off biomass) and after must be same
        !>
        sum1=totdmas1(i)+grsumpap(i)
        diff=abs(sum1-ntotdms1(i))
        if (ntotdms1(i)==0.) then
          if (diff==0.) then
            rdif=0.             ! relative accuracy
          else
            rdif=1.             ! relative accuracy
          end if
        else
          rdif=diff/ntotdms1(i) ! relative accuracy
        end if

        !if (abs(totdmas1(i) + grsumpap(i) - ntotdms1(i)) > tolrnce1) then
        if (rdif > rel_tolerance_amount) then
          write(6, * )'at grid cell = ',i
          write(6, * )'total litter carbon does not balance after luc'
          write(6, * )'totdmas1(i) = ',totdmas1(i)
          write(6, * )'grsumpap(i) = ',grsumpap(i)
          write(6, * )'totdmas1(i) + grsumpap(i) = ', &
                     totdmas1(i) + grsumpap(i)
          write(6, * )'ntotdms1(i) = ',ntotdms1(i)
          write(6, * )'Absolute diff of ',abs(totdmas1(i) + grsumpap(i) - ntotdms1(i)),' is > than ',tolrnce1
          write(6, * )'ABS. DIFF. = ',diff
          write(6, * )'REL. DIFF. = ',rdif,' and should be less than ',rel_tolerance_amount
          call errorHandler('luc', - 12)
        end if


        if( Ncycle_on .and. abs(totdnmas1(i)+ grsumnpap(i) - ntotdnms1(i)) > tolrnce1 ) then
           write(6, * )'at grid cell = ',i
           write(6, * )'total litter nitrogen does not balance after luc'
           write(6, * )'totdnmas1(i) = ',totdnmas1(i)
           write(6, * )'grsumnpap(i) = ',grsumnpap(i)
           write(6, * )'totdnmas1(i) + grsumnpap(i) = ', &
                      totdnmas1(i) + grsumnpap(i)
           write(6, * )'ntotdnms1(i) = ',ntotdnms1(i)
           call errorHandler('luc', - 13)
        endif


        !>
        !> for conservation totcmass(i) must be equal to ntotcmas(i) plus combustion carbon losses
        !>
        sum1=ntotcmas(i)+grsumcom(i)
        diff=abs(totcmass(i)-sum1)
        if (totcmass(i)==0.) then
          if (diff==0.) then
            rdif=0.             ! relative accuracy
          else
            rdif=1.             ! relative accuracy
          end if
        else
          rdif=diff/totcmass(i) ! relative accuracy
        end if

        !if (abs(totcmass(i) - ntotcmas(i) - grsumcom(i)) > tolrnce1) then
        if (rdif > rel_tolerance_amount) then
          write(6, * )'at grid cell = ',i
          write(6, * )'total carbon does not balance after luc'
          write(6, * )'totcmass(i) = ',totcmass(i)
          write(6, * )'ntotcmas(i) = ',ntotcmas(i)
          write(6, * )'grsumcom(i) = ',grsumcom(i)
          write(6, *)'Absolute diff of ',abs(totcmass(i) - ntotcmas(i) - grsumcom(i)),' is compared to ',tolrnce1
          write(6, * )'ABS. DIFF. = ',diff
          write(6, * )'REL. DIFF. = ',rdif,' and should be less than ',rel_tolerance_amount
          call errorHandler('luc', - 14)
        end if

        if ( Ncycle_on .and. abs(totnmass(i) - ntotnmas(i) - grsumncom(i)) > tolrnce1) then
           write(6, * )'at grid cell = ',i
           write(6, * )'total nitrogen does not balance after luc'
           write(6, * )'totnmass(i) = ',totnmass(i)
           write(6, * )'ntotnmas(i) = ',ntotnmas(i)
           write(6, * )'grsumncom(i)= ',grsumncom(i)
           call errorHandler('luc', - 15)
        end if
        !>
        !> update grid averaged vegetation biomass,and litter and soil c densities
        !>
        do j = 1,icc
          vgbiomas(i) = vgbiomas(i) + fcancmx(i,j) * (gleafmas(i,j) + &
                        bleafmas(i,j) + stemmass(i,j) + rootmass(i,j))
          nvgbiomas(i) = nvgbiomas(i) + fcancmx(i,j) * (ngleafmas(i,j) + &
                        nbleafmas(i,j) + nstemmass(i,j) + nrootmass(i,j))
          do k = 1,ignd
            gavgltms(i)=gavgltms(i)+fcancmx(i,j)*litrmass(i,j,k)
            gavgscms(i)=gavgscms(i)+fcancmx(i,j)*soilcmas(i,j,k)
          end do
          nlitrmass_g(i)=nlitrmass_g(i)+fcancmx(i,j)*nlitrmass(i,j)
          soilnmas_g(i)=soilnmas_g(i)+fcancmx(i,j)*soilnmas(i,j)
          nh4_mass_g(i)=nh4_mass_g(i)+fcancmx(i,j)*nh4_mass(i,j)
          no3_mass_g(i)=no3_mass_g(i)+fcancmx(i,j)*no3_mass(i,j)
        end do ! loop 750
        do k = 1,ignd
          gavgltms(i)=gavgltms(i)+( barefrac(i)*litrmass(i,iccp1,k) )
          gavgscms(i)=gavgscms(i)+( barefrac(i)*soilcmas(i,iccp1,k) )
        end do

        !>
        !> just like total amount of carbon must balance,the grid averaged densities must also balance
        !>

        sum1=pvgbioms(i)+pgavltms(i)+pgavscms(i)+pluclitpool(i)+plucscpool(i)
        sum2=vgbiomas(i)+gavgltms(i)+gavgscms(i)+grdencom(i)+litrmass(i,iccp2,1)+soilcmas(i,iccp2,1)
        diff=abs(sum1-sum2)
        if (sum1==0.) then
          if (diff==0.) then
            rdif=0.              ! relative accuracy
          else
            rdif=1.              ! relative accuracy
          end if
        else
          rdif=diff/sum1         ! relative accuracy
          if (diff < tolrance) rdif=0.
        end if

        if (rdif > rel_tolerance_amount) then
        !if ( abs(pvgbioms(i) + pgavltms(i) + pgavscms(i) &
        !         + pluclitpool(i) + plucscpool(i) &
        !         - vgbiomas(i) - gavgltms(i) - gavgscms(i) &
        !         - litrmass(i,iccp2,1) - soilcmas(i,iccp2,1) &
        !         - grdencom(i)) > tolrance) then
          write(6, * )'at grid cell = ',i
          write(6, * )'total carbon density does not balance after luc'
          write(6, * )'pvgbioms(i) = ',pvgbioms(i)
          write(6, * )'pgavltms(i) = ',pgavltms(i)
          write(6, * )'pgavscms(i) = ',pgavscms(i)
          write(6, * )'pluclitpool(i) = ',pluclitpool(i)
          write(6, * )'plucscpool(i) = ',plucscpool(i)
          write(6, * )'vgbiomas(i) = ',vgbiomas(i)
          write(6, * )'gavgltms(i) = ',gavgltms(i)
          write(6, * )'gavgscms(i) = ',gavgscms(i)
          write(6, * )'grdencom(i) = ',grdencom(i)
          write(6, * )'litrmass(i,iccp2,1) = ',litrmass(i,iccp2,1)
          write(6, * )'soilcmas(i,iccp2,1) = ',soilcmas(i,iccp2,1)
          write(6, * )' '
          write(6, * )'diff = ',abs((pvgbioms(i) + pgavltms(i) + pgavscms(i) + pluclitpool(i) + plucscpool(i)) &
          - (vgbiomas(i) + gavgltms(i) + gavgscms(i) + grdencom(i) + litrmass(i,iccp2,1) + soilcmas(i,iccp2,1) ))
          write(6, * )'tolrance = ',tolrance
          write(6, * )'ABS. DIFF. = ',diff,' and should be less than ',tolrance
          write(6, * )'REL. DIFF. = ',rdif,' and should be less than ',rel_tolerance_amount
          call errorHandler('luc', - 16)
        end if

        do j = 1,icc
          if ((.not. leapnow .and. iday == 365) .or. (leapnow .and. iday == 366)) then
            pfcancmx(i,j) = nfcancmx(i,j)
          end if
        end do ! loop 800

      end if  ! loop to check if any luc took place.
    end do ! loop 255

    end associate
    return

  end subroutine luc
  !! @}
  !=======================================================================
  !> \ingroup landusechange_adjustLucFracs
  !> @{
  !> Adjusts the amount of each pft to ensure that the fraction of
  !> gridcell bare ground is >0.
  !> @author Joe Melton
  subroutine adjustLucFracs (i, onetile_perPFT, nfcancmxrow, &
                             bare_ground_frac, PFTCompetition)

    use classicParams, only : nlat, nmos, icc

    implicit none

    ! arguments:
    integer, intent(in) :: i
    real, dimension(nlat,nmos,icc), intent(inout) :: nfcancmxrow
    real, intent(in) :: bare_ground_frac
    logical, intent(in) :: onetile_perPFT
    logical, intent(in) :: PFTCompetition

    real, dimension(nlat,nmos,icc) :: outnfcrow

    ! local variables:
    integer :: m, j
    real, dimension(icc) :: frac_abv_min_val
    real :: needed_bare, reduce_per_pft, tot
    real :: min_val

    !-------------------------
    tot = 0.0

    if (PFTCompetition) then
      min_val = 0.
    else
      min_val = 0.
    end if

    !> find the amount of needed space in the other pfts
    !> need a minimum bare area of min_val.
    needed_bare = ( - bare_ground_frac) + min_val

    !> get the proportionate amounts of each pft above the min_val lower limit
    do j = 1,icc
      if (onetile_perPFT) then
        m = j
      else
        m = 1
      end if

      frac_abv_min_val(j) = max(0.0,nfcancmxrow(i,m,j) - min_val)
      tot = tot + frac_abv_min_val(j)

    end do

    !> add in the bare ground min fraction of min_val.
    tot = tot + min_val

    !> now reduce the pfts proportional to their fractional area to make
    !> the bare ground fraction be the min_val amount and no other pft be less than min_val
    do j = 1,icc
      if (onetile_perPFT) then
        m = j
      else
        m = 1
      end if
      outnfcrow(i,m,j) = max(min_val,nfcancmxrow(i,m,j) - needed_bare * &
                         frac_abv_min_val(j) / tot)
      nfcancmxrow(i,m,j) = outnfcrow(i,m,j)
    end do

  end subroutine adjustLucFracs
  !! @}
  !=======================================================================
  !> \ingroup landusechange_adjustFracsComp
  !! @{
  !> Used when PFTCompetition = true. It adjusts the amount of each pft
  !> to allow expansion of cropland.
  !> @author Joe Melton

  subroutine adjustFracsComp (il1, il2, nilg, iday, pfcancmx, yesfrac, delfrac, outdelfrac)

    use classicParams, only : icc, crop, zero

    implicit none

    ! arguments:
    integer, intent(in) :: il1
    integer, intent(in) :: il2
    integer, intent(in) :: nilg
    integer, intent(in) :: iday
    real, dimension(nilg,icc), intent(in) :: pfcancmx
    real, dimension(nilg,icc), intent(in) :: yesfrac
    real, dimension(nilg,icc), intent(in) :: delfrac
    real, dimension(nilg,icc), intent(out) :: outdelfrac

    ! local variables:
    integer :: i, j
    real, dimension(nilg) :: chgcrop, cropfrac
    real, dimension(nilg,icc) :: adjus_fracs, fmx, fmy
    real, dimension(nilg) :: barefrac
    real, dimension(nilg,1,icc) :: tmpfcancmx

    real, parameter :: smallnumber = 1.0e-12

    !> -------------------------

    !> Some initializations
    chgcrop = 0.
    cropfrac = 0.
    adjus_fracs = 0.
    outdelfrac = 0.
    barefrac = 1.0

    !> Find how much the crop area changes this timestep. We only care about the total
    !> change,not the per crop PFT change.
    do i = il1,il2
      do j = 1,icc
        if (crop(j)) then
          fmx(i,j) = pfcancmx(i,j) + (real(iday) * delfrac(i,j)) !  current day
          fmy(i,j) = pfcancmx(i,j) + (real(iday - 1) * delfrac(i,j)) ! previous day
          chgcrop(i) = chgcrop(i) + (fmx(i,j) - fmy(i,j))
          cropfrac(i) = cropfrac(i) + fmx(i,j)
        end if
      end do
    end do
    !>
    !> If the crop area changed we have to reduce the other PFT proportional to their
    !> area (if we gained crop area). We don't presently assume anything like grasslands
    !> are converted first. We assume that on the scale of our gridcells, area
    !> is simply converted proportional to the total.
    do i = il1,il2

      if (chgcrop(i) > smallnumber) then

        !> Adjust the non-crop PFT fractions to find their proportional fraction that does not include crops.
        do j = 1,icc
          if (.not. crop(j)) then

            adjus_fracs(i,j) = yesfrac(i,j) / (1. - cropfrac(i))
            outdelfrac(i,j) = chgcrop(i) * adjus_fracs(i,j)

          end if

          if (.not. crop(j)) then
            barefrac(i) = barefrac(i) - (yesfrac(i,j) - outdelfrac(i,j))
            tmpfcancmx(i,1,j) = yesfrac(i,j) - outdelfrac(i,j)
          else  ! crops
            barefrac(i) = barefrac(i) - fmx(i,j)
            tmpfcancmx(i,1,j) = fmx(i,j)
          end if
        end do

        if (barefrac(i) < 0.0) then
          !> tmpfcancmx has a nilg,1,icc dimension as adjustLucFracs expects
          !> a 'row' structure of the array.
          call adjustLucFracs(i,.false.,tmpfcancmx,barefrac(i),.true.)
          do j = 1,icc
            outdelfrac(i,j) = tmpfcancmx(i,1,j) - yesfrac(i,j)
          end do
        end if

      end if

    end do

  end subroutine adjustFracsComp
  !! @}

  !> \namespace landusechange
  !> Central module for all land use change operations
  !!
  !! The land use change (LUC) module of CTEM is based on
  !! (Arora and Boer, 2010) \cite Arora2010-416.
  !! When the area of crop PFTs changes, CTEM generates LUC emissions.
  !! In the simulation where fractional coverage of PFTs is specified, the changes in
  !! fractional coverage of crop PFTs are made consistent with changes in the
  !! fractional coverage of natural non-crop PFTs. That is, an increase or decrease
  !! in the area of crop PFTs is associated with a corresponding decrease or increase
  !! in the area of non-crop PFTs. This approach is taken by Wang et al. (2006) \cite Wang2006-he, which
  !! allows one to reconstruct historical land cover given a spatial data set of
  !! changes in crop area over the historical period and an estimate of potential
  !! natural land cover for the pre-industrial period (as described in Sect.
  !! \ref{methods}). When competition between PFTs for space is allowed, only the
  !! fractional coverage of crop PFTs is specified. Similar to a simulation with
  !! prescribed PFT fractions, when the area of crop PFTs increases, the fractional
  !! coverage of non-crop PFTs is decreased in proportion to their existing coverage
  !! (Wang et al. (2006) \cite Wang2006-he. Alternatively, and in contrast to the
  !! simulation with prescribed
  !! PFT fractions, when the area of crop PFTs decreases then the generated bare
  !! fraction is available for recolonization by non-crop PFTs.
  !!
  !! A decrease in the area of natural non-crop PFTs, associated with an increase in
  !! area of crop PFTs, results in deforested biomass (while the term
  !! \f$\textit{deforested}\f$ implies clearing of forests, the same processes can
  !! occur in grasslands as well and is meant here to imply removal of the biomass).
  !! The deforested biomass is divided into three components: (i) the component that
  !! is combusted or used for fuel wood immediately after natural vegetated is
  !! deforested and which contributes to atmospheric \f$CO_2\f$, (ii) the component
  !! left as slash or used for pulp and paper products and (iii) the component that
  !! is used for long-lasting wood products. The fractions allocated to these three
  !! components depend on whether the PFTs are woody or herbaceous and on their
  !! aboveground vegetation biomass density (see Table 1 of \cite Arora2010-416).
  !! To account for the timescales involved, the fraction allocated to slash or pulp
  !! and paper products is transferred to the model's litter pool and the fraction
  !! allocated to long-lasting wood products is allocated to the model's soil carbon
  !! pool. While they are place in the litter and soil carbon pools, these LUC products
  !! are kept separate to allow accurate accounting. Land use change associated
  !! with a decrease in the area of natural
  !! vegetation thus redistributes carbon from living vegetation to dead litter
  !! and soil carbon pools and emits \f$CO_2\f$ to the atmosphere through direct
  !! burning of the deforested biomass. The net result is positive LUC carbon
  !! emissions from land to the atmosphere.
  !!
  !! When croplands are abandoned, the area of natural PFTs increases. In simulations
  !! with prescribed fractional coverage of PFTs this results in a decreased carbon
  !! density for all model pools as the same amount of carbon is spread over a larger
  !! fraction of the grid cell. This reduced density implies that natural vegetation
  !! is able to take up carbon as it comes into equilibrium with the driving climate
  !! and atmospheric \f$CO_2\f$ concentration. This creates the carbon sink associated
  !! with abandonment of croplands as natural vegetation grows in its place. In
  !! simulations with competition between PFTs, the abandoned land is treated as
  !! bare ground, which is subsequently available for recolonization, as mentioned
  !! above. As natural vegetation expands into bare ground it takes up carbon, again
  !! creating the carbon sink associated with abandonment of croplands. The net
  !! result is negative LUC carbon emissions as carbon is taken from atmosphere to
  !! grow vegetation over the area that was previously a cropland.
  !!

end module landuseChange
