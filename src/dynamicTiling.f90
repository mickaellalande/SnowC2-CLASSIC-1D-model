!> \file
!> Dynamic tiling toolkit
module dynamicTiling

  ! S.R. Curasi April, 2022

  !Later on 
  ! 1. generalize this to make better use of the control vector

  use classicParams,  only : nlat, nmos, ilg, ican, ignd, icp1, icc, iccp2, iccp1,&
                             dynTilingTolrance, convertg2kg, hcpice, hcpsnd, hcpw,&
                             rhoice, rhow, sphice, sphveg, sphw, dynTilingTolranceEnergy,&
                             deltat, percentthreshold, nTilepres

  implicit none

  ! Subroutines contained in this module:
    public  :: moveSplitCopyTiles
    public  :: cleanTiles
    public  :: indexTileAge
    public  :: initialMassEnergyBalance
    public  :: FinalMassEnergyBalance
    public  :: MixTempsRichmann
    public  :: updateHeight
    public  :: majority3d
    public  :: majorityarea3d
    public  :: majorityarea1d
    public  :: majorityarea2dr8
    public  :: weightedAvg2d
    public  :: weightedAvg3d
    public  :: weightedAvg4d
    public  :: weightedAvg3dr8
    public  :: indexSort
    
contains

  ! ------------------------------------------------------------------

  !> \ingroup dynamicTiling_moveSplitCopyTiles
  !! @{
  !> Tile moving, splitting, and copying subroutine is the main feature of this toolkit (see the module-level description for an overview)
  !> @author S.R. Curasi
  subroutine moveSplitCopyTiles(controlVector,inputVector,FAREAROUT,mode,DynTilInitializeFlag,outputIndex, & !these are the main intputs that control the routine
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
                                useTracer, Ncycle_on, DELZ, THPROT, HCPSROT, tileAgerow) !these are switches we need

    use, intrinsic :: iso_fortran_env, only: r8=>real64

    implicit none

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
    real, intent(inout), dimension(nlat,nmos,ignd) :: THPROT !< soil pore volume
    real, intent(inout), dimension(nlat,nmos,ignd) :: HCPSROT !< soil heat capacity
    real, intent(inout), dimension(nlat,nmos) :: tileAgerow 

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

    !declare input variables that control the subroutine
    integer, intent(in), dimension(nlat,nmos) :: controlVector !< the control vector specified in the initializtion file
    integer, intent(in), dimension(nlat,nmos) :: inputVector !< the input vector (1 = operate on tile, 0 = do not operate on tile) 
    character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, copy, or skip) 
    real, intent(in), dimension(nlat) :: FAREAROUT !< a real value greater than zero and less than one which specifies the desired size of the output tile
    logical, intent(inout) :: DynTilInitializeFlag !< a flag that's set if the model needs to be re-initialized

    !declare variables used within the subroutine
    integer :: i, k, j, l
    integer, dimension(nlat) :: inputLength !< the number of tiles to be opperated on (1s in the input vector)
    integer, dimension(nlat) :: headTile !< the index of the head tile
    integer, dimension(nlat), intent(out) :: outputIndex!< the index of the tile where the output will be written
    real, dimension(nlat) :: FAREFINAL !< the fractional area that will be written to the output tile at the end of the operation
    real, dimension(nlat,nmos) :: FAREINPUT !< the fractionl area that will be written to the intput tiles at the end of the operation
    integer, allocatable, dimension(:,:) :: inputIndex !< the index of the input tiles (this is alocatable because the amount of active tiles can vary)

    !declare the variables used to run the mass balance checks 
    real, dimension(nlat) :: vegCMass_initial !< temporary storage variable for the total mass of carbon in vegetation
    real, dimension(nlat) :: soilCMass_initial !< temporary storage variable for the total mass of carbon in soil and litter
    real, dimension(nlat) :: peatCMass_initial !< temporary storage variable for the total mass of carbon in the peatland pools
    real, dimension(nlat) :: lucMass_initial !< temporary storage variable for the total mass of carbon in the LUC pools
    real, dimension(nlat) :: totalCMass_initial !< temporary storage variable for the total mass of carbon in all pools
    real, dimension(nlat) :: soilWaterMass_initial !< temporary storage variable for the total mass of water in the soil
    real, dimension(nlat) :: csWaterMass_initial !< temporary storage variable for the total mass of water in the canopy and on the soil surface  
    real, dimension(nlat) :: totalWaterMass_initial !< temporary storage variable for the total mass of water in all pools
    real, dimension(nlat) :: vegNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, dimension(nlat) :: soilNMass_initial !< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, dimension(nlat) :: lucNmass_initial !< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, dimension(nlat) :: totalNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, dimension(nlat) :: TVSTP_initial !< internal energy of vegetation used for energy balance checks (J)
    real, dimension(nlat) :: TSSTP_initial !< internal energy of snow pack used for energy balance checks (J)
    real, dimension(nlat) :: TSTP_initial !< internal energy of soil used for energy balance checks (J) 
    real, dimension(nlat) :: TSPN_initial !< internal energy ponded water used for energy balance checks (J)
    real, dimension(nlat) :: TTOT_initial !< total internal energy of all pools used for energy balance checks (J)

    !print everything S.R.C used to debug remove later on
    !print *,"-----at startup----"
    !print *,"controlVector = ", controlVector(:,:)
    !print *,"inputVector = ", inputVector(:,:)
    !print *,"FAREINPUT = ", FAREINPUT(:,:)
    !print *,"FAREROT = ", FAREROT(:,:)

    !>First we run a bunch of checks and do some initial calculations
    !!these are just pre-calculations and don't make any modifications

    !>check the the controlVector contains a head tile
    if(any(controlVector(:,1) /= -1) .or. any(controlVector(:,2:nmos) == -1)) then
        write(6, * ) 'moveSplitCopyTiles: missing or incorectly placed head tile'
        call errorHandler('moveSplitCopyTiles', -1)
    else
        headTile(:) = 1
    end if 

    !> check that the values in the inputVector are valid
    if(any(inputVector /= 0 .and. inputVector /= 1)) then
        write(6, * ) 'moveSplitCopyTiles: incorrect inputs in inputVector'
        call errorHandler('moveSplitCopyTiles', -2)
        stop
    end if

    !> check fareout and find the indicies for the input and output tiles
    outputIndex(:) = 0
    do j=1,nlat
        if (mode(j) /= "skip") then
            i = 1
            if (mode(j) == "merge" .and. FAREAROUT(j) == -9999.0) then
                wloop1: do while (outputIndex(j)==0)
                    if(i > nmos) then
                        write(6, * ) 'moveSplitCopyTiles: problem determining output index from inputVector'
                        print *, "mode = ",mode(j)
                        print *, "FAREAROUT = ",FAREAROUT(j)
                        !print *, "inputVector(j,:) = ",inputVector(j,:)
                        !print *, "FAREROT(j,:) = ",FAREROT(j,:)
                        !print *, "controlVector(j,:) = ",controlVector(j,:)
                        call errorHandler('moveSplitCopyTiles', -3)
                        !these exit statement are important or else the code can get stuck here
                        exit wloop1
                    else if(inputVector(j,i) == 1) then
                        outputIndex(j) = i
                    end if
                    i = i + 1
                end do wloop1
            else if (mode(j) == "split" .and. ((FAREAROUT(j) > 0.0) .or. (FAREAROUT(j) <= 1.0))) then
                wloop2: do while (outputIndex(j)==0)
                    if(i > nmos) then
                        write(6, * ) 'moveSplitCopyTiles: problem determining output index from inputVector'
                        print *, "mode = ",mode(j)
                        print *, "FAREAROUT = ",FAREAROUT(j)
                        !print *, "inputVector(j,:) = ",inputVector(j,:)
                        !print *, "FAREROT(j,:) = ",FAREROT(j,:)
                        !print *, "controlVector(j,:) = ",controlVector(j,:)
                        call errorHandler('moveSplitCopyTiles', -3)
                        exit wloop2
                    else if(inputVector(j,i) == 0 .and. FAREROT(j,i) == 0.0 .and. abs(controlVector(j,i)) == 1) then
                        outputIndex(j) = i
                    end if
                    i = i + 1
                end do wloop2
            else if (mode(j) == "move" .and. FAREAROUT(j) == -9999.0) then
                wloop3: do while (outputIndex(j)==0)
                    if(i > nmos) then
                        write(6, * ) 'moveSplitCopyTiles: problem determining output index from inputVector'
                        print *, "mode = ",mode(j)
                        print *, "FAREAROUT = ",FAREAROUT(j)
                        !print *, "inputVector(j,:) = ",inputVector(j,:)
                        !print *, "FAREROT(j,:) = ",FAREROT(j,:)
                        !print *, "controlVector(j,:) = ",controlVector(j,:)
                        call errorHandler('moveSplitCopyTiles', -3)
                        exit wloop3
                    else if(inputVector(j,i) == 0 .and. FAREROT(j,i) == 0.0 .and. abs(controlVector(j,i)) == 1) then
                        outputIndex(j) = i
                    end if
                    i = i + 1
                end do wloop3
            else
                write(6, * ) 'moveSplitCopyTiles: invalid mode or FAREOUT given as input'
                call errorHandler('moveSplitCopyTiles', -3)
            end if
        end if
    end do

    !> find the itterator for the tiles in the input vector and check none of the input tiles are protected
    do i=1,nlat
        if (mode(i) /= "skip") then
            inputLength(i) = sum(inputVector(i,:))
            if (inputLength(i) < 1 .or. inputLength(i) > nmos) then
                write(6, * ) 'moveSplitCopyTiles: invalid number of input tiles given'
                call errorHandler('moveSplitCopyTiles', -4)
                stop
            end if 
            if (mode(i) == "move" .and. inputLength(i) > 1) then
                write(6, * ) 'moveSplitCopyTiles: you cannot move more than one tile'
                call errorHandler('moveSplitCopyTiles', -5)
            end if
        end if
    end do

    !> then this information is used to setup the inputindex vector
    allocate(inputIndex(nlat,maxval(inputLength)))
    do k=1,nlat
        if (mode(k) /= "skip") then
            j = 1
            do i=1,nmos
                if(inputVector(k,i) == 1) then
                    if(controlVector(k,i) == 0) then
                        write(6, * ) 'moveSplitCopyTiles: protected tile given as input'
                        call errorHandler('moveSplitCopyTiles', -6)
                        stop
                    else if (FAREROT(k,i) == 0.0) then
                        write(6, * ) 'moveSplitCopyTiles: inactive tile given as input'
                        call errorHandler('moveSplitCopyTiles', -7)
                        stop
                    end if
                    inputIndex(k,j) = i
                    j = j + 1
                end if  
            end do
        end if
    end do

    !> then we find the fractional area for the output tile and check that it's not larger than the input tile areas
    do i=1,nlat
        if (mode(i) /= "skip") then
            if (mode(i) == "merge" .and. FAREAROUT(i) == -9999.0) then
                FAREFINAL = sum(FAREROT(i,inputIndex(i,:)))
            else if (mode(i) == "move" .and. FAREAROUT(i) == -9999.0) then
                FAREFINAL = sum(FAREROT(i,inputIndex(i,:)))
            else if (mode(i) == "split" .and. (FAREAROUT(i) > 0.0 .or. FAREAROUT(i) < 1.0)) then
                if(FAREAROUT(i) > sum(FAREROT(i,inputIndex(i,:))) .or. FAREAROUT(i) >= 1.0) then
                    write(6, * ) 'moveSplitCopyTiles: either not enough fractional area in input tiles'
                    write(6, * ) 'or the split operation should be a merge because it calls for an output area of 1.0'
                    write(6, * ) 'FAREAROUT(i) = ',FAREAROUT(i)
                    call errorHandler('moveSplitCopyTiles', -8)
                else
                    FAREFINAL(i) = FAREAROUT(i)
                end if       
            else
                write(6, * ) 'moveSplitCopyTiles: incorrect inputs in FAREAROUT'
                call errorHandler('moveSplitCopyTiles', -9)
            end if
        end if
    end do

    !< then we check that fcancmxrow is the same in all the tiles under consideration
    do j=1,nlat
        if (mode(j) /= "skip") then
            do i=1,inputLength(j)
                if(any(fcancmxrow(j,headTile(j),:) /= fcancmxrow(j,inputIndex(j,i),:))) then
                    write(6, * ) 'moveSplitCopyTiles: different fractional covers in input and head tiles'
                    call errorHandler('moveSplitCopyTiles', -10)
                end if    
            end do
        end if
    end do

    !check the fractional area of all tiles
    do i=1,nlat
        if (mode(i) /= "skip") then
            if ((sum(FAREROT(i,:))-1.0 > dynTilingTolrance) .or. (1.0-sum(FAREROT(i,:)) > dynTilingTolrance)) then
                write(6, * ) 'moveSplitCopyTiles: fractional area of the input tiles is not one'
                write(6, * ) 'mode = ', mode
                write(6, * ) 'test threshold = ', dynTilingTolrance
                write(6, * ) 'test value = ', sum(FAREROT(i,:))-1
                !write(6, * ) 'farerot =', FAREROT(i,:)
                call errorHandler('moveSplitCopyTiles', -11)
            end if
        end if
    end do

    FAREINPUT=-9999
    do j=1,nlat
        if (mode(j) /= "skip") then
            !find new fractional areas for the new input tiles
            if (mode(j) == "merge" .and. FAREAROUT(j) == -9999.0) then
                    FAREINPUT(j,inputIndex(j,:)) = 0.0
                    FAREINPUT(j,outputIndex(j))=-9999.0
            else if (mode(j) == "move" .and. FAREAROUT(j) == -9999.0) then
                    FAREINPUT(j,inputIndex(j,:)) = 0.0
                    FAREINPUT(j,outputIndex(j))=-9999.0
            else if (mode(j) == "split" .and. (FAREAROUT(j) > 0.0 .or. FAREAROUT(j) <= 1.0)) then
                do i=1,inputLength(j)
                    FAREINPUT(j,inputIndex(j,i)) = (FAREROT(j,inputIndex(j,i))-((FAREFINAL(j)*FAREROT(j,inputIndex(j,i)))/sum(FAREROT(j,inputIndex(j,:)))))
                end do
            else
                write(6, * ) 'moveSplitCopyTiles: incorrect inputs in FAREAROUT'
                call errorHandler('moveSplitCopyTiles', -12)
            end if
        end if
    end do

    !print everything S.R.C used to debug remove later on
    !print *,"-----after startup----"
    !print *,"headTile = ", headTile(:)
    !print *,"outputIndex = ", outputIndex(:)
    !print *,"FAREAROUT = ", FAREAROUT(:)
    !print *,"FAREFINAL = ", FAREFINAL(:)
    !print *,"inputLength = ", inputLength(:)
    !print *,"inputIndex = ", inputIndex(:,:)

    !> then we perform the intial mass balance checks so we know the model state before any changes were made
    call initialMassEnergyBalance(FAREROT, DLZWROT, fcancmxrow, gleafmasrow, bleafmasrow, stemmassrow, rootmassrow, Cmossmasrow, peatSoilCrow, litrmsmossrow, &
                                soilcmasrow, litrmassrow, nh4_massrow, no3_massrow, nlitrmassrow, soilnmasrow, &
                                ngleafmasrow, nbleafmasrow ,nstemmassrow, nrootmassrow ,RCANROT ,SCANROT, TCANROT, TSNOROT, SNOROT, WSNOROT, TBARROT, &
                                TPNDROT, THICROT, THLQROT, ZPNDROT, vegCMass_initial, soilCMass_initial, peatCMass_initial, lucMass_initial, totalCMass_initial, &
                                soilWaterMass_initial, csWaterMass_initial, totalWaterMass_initial, vegNMass_initial, soilNMass_initial, lucNMass_initial, totalNMass_initial, TVSTP_initial, & 
                                TSSTP_initial, TSTP_initial, TSPN_initial, TTOT_initial, Ncycle_on, ipeatlandrow, DELZ, THPROT, HCPSROT, CMAIROT, mode)

    !> then we move/copy/average information between tiles

    !> temperatures use Richmann’s law (based upon energyWaterBalanceCheck.f90) 
    !> this has to be done first because a lot of other variables are required
    call MixTempsRichmann(FAREROT, DLZWROT, fcancmxrow, DELZ, THPROT, HCPSROT, & 
                            THLQROT, THICROT, RCANROT, SCANROT, SNOROT, WSNOROT, ZPNDROT, TCANROT, TBARROT, &
                            TSNOROT, TPNDROT, CMAIROT, inputLength, outputIndex, inputIndex, mode)

    !> We assign variables which are not read in from the initilization file and assumed to be spatialy invariate
    do i=1,nlat
        if (mode(i) /= "skip") then
            !! the values from the head tile, if they were read in this code might change
            RSMNROT(i,outputIndex(i),:) = RSMNROT(i,headTile(i),:)
            QA50ROT(i,outputIndex(i),:) = QA50ROT(i,headTile(i),:)
            VPDAROT(i,outputIndex(i),:) = VPDAROT(i,headTile(i),:)
            VPDBROT(i,outputIndex(i),:) = VPDBROT(i,headTile(i),:)
            PSGAROT(i,outputIndex(i),:) = PSGAROT(i,headTile(i),:)
            PSGBROT(i,outputIndex(i),:) = PSGBROT(i,headTile(i),:)

            !> Then some variables like DRNROT and SDEPROT that are read in from the init file are assigned values from the head tile
            DRNROT(i,outputIndex(i)) = DRNROT(i,headTile(i))
            SDEPROT(i,outputIndex(i)) = SDEPROT(i,headTile(i))
            SOCIROT(i,outputIndex(i)) = SOCIROT(i,headTile(i))
            SANDROT(i,outputIndex(i),:) = SANDROT(i,headTile(i),:)
            CLAYROT(i,outputIndex(i),:) = CLAYROT(i,headTile(i),:)
            ORGMROT(i,outputIndex(i),:) = ORGMROT(i,headTile(i),:)
            ipeatlandrow(i,outputIndex(i)) = ipeatlandrow(i,headTile(i))
            FCANROT(i,outputIndex(i),:) = FCANROT(i,headTile(i),:) !this gets copied over from the head tile but probably gets overwritten soon after
            slopefracrow(i,outputIndex(i),:) = slopefracrow(i,headTile(i),:)
            DLZWROT(i,outputIndex(i),:) = DLZWROT(i,headTile(i),:) 
            HCPSROT(i,outputIndex(i),:) = HCPSROT(i,headTile(i),:) 
            THPROT(i,outputIndex(i),:) = THPROT(i,headTile(i),:) 

            !> lfstatus could be set to zero so it's value is re-set by phenology.f90
            ! lfstatusrow(i,outputIndex(i),:) = 0
        end if
    end do

    !> other especilly mass based variables like RCANROT and SCANROT used weighted averaged
    RCANROT = weightedAvg2d(RCANROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)   
    SCANROT = weightedAvg2d(SCANROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)  
    SNOROT = weightedAvg2d(SNOROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
    RHOSROT = weightedAvg2d(RHOSROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)  
    !in the modelStateDrivers.f90 consistency checks are done so water in soil should be consistent here as well
    THLQROT = weightedAvg3d(THLQROT,ignd,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    THICROT = weightedAvg3d(THICROT,ignd,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    Cmossmasrow = weightedAvg2d(Cmossmasrow,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    litrmsmossrow = weightedAvg2d(litrmsmossrow,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    peatSoilCrow = weightedAvg2d(peatSoilCrow,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
    litrmassrow = weightedAvg4d(litrmassrow,iccp2,ignd,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    soilcmasrow = weightedAvg4d(soilcmasrow,iccp2,ignd,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    dmossrow = weightedAvg2d(dmossrow,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    lystmmasrow = weightedAvg3d(lystmmasrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    lyrotmasrow = weightedAvg3d(lyrotmasrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    lygleafmasmaxrow = weightedAvg3d(lygleafmasmaxrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    lystemmassmaxrow = weightedAvg3d(lystemmassmaxrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    lyrootmassmaxrow = weightedAvg3d(lyrootmassmaxrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    bleafmasrow = weightedAvg3d(bleafmasrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    rootmasssrow = weightedAvg3d(rootmasssrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    rootmass_NSrow = weightedAvg3d(rootmass_NSrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    stemmasssrow = weightedAvg3d(stemmasssrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    stemmass_NSrow = weightedAvg3d(stemmass_NSrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    gleafmas_NSrow = weightedAvg3d(gleafmas_NSrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    gleafmassrow = weightedAvg3d(gleafmassrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    tileAgerow = weightedAvg2d(tileAgerow,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)

    !these are partial pressures in Pa they should be very similar to the mass related variables
    co2i1cgrow = weightedAvg3d(co2i1cgrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    co2i1csrow = weightedAvg3d(co2i1csrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    co2i2cgrow = weightedAvg3d(co2i2cgrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    co2i2csrow = weightedAvg3d(co2i2csrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)

    !> other temperature related variables
    ZPNDROT = weightedAvg2d(ZPNDROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)  
    maxAnnualActLyrROT = weightedAvg2d(maxAnnualActLyrROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    tymaxlairow = weightedAvg3d(tymaxlairow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)

    !leaf status is set to the majority value
    lfstatusrow = majorityarea3d(lfstatusrow,icc,nlat,nmos,inputLength,outputIndex,inputIndex,FAREROT,mode)

    !these are more complex variables which are again average
    GROROT = weightedAvg2d(GROROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    grwtheffrow = weightedAvg3d(grwtheffrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    flhrloss_nsrow = weightedAvg3d(flhrloss_nsrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    flhrloss_srow = weightedAvg3d(flhrloss_srow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    colddays_harvestrow = INT(weightedAvg2d(REAL(colddays_harvestrow),nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode))
    colddays_leaffallrow = INT(weightedAvg2d(REAL(colddays_leaffallrow),nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode))
    cfluxcgrow = weightedAvg2d(cfluxcgrow,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    cfluxcsrow = weightedAvg2d(cfluxcsrow,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    pandaysrow = INT(weightedAvg3d(REAL(pandaysrow),icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode))
    rothrlosrow = weightedAvg3d(rothrlosrow ,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    stmhrlosrow = weightedAvg3d(stmhrlosrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)

    !these should be set on restary but are coppied here as placeholders
    ALBSROT = weightedAvg2d(ALBSROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    ALICROT = weightedAvg3d(ALICROT,icp1,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    ALVCROT = weightedAvg3d(ALVCROT,icp1,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)

    ROOTROT = weightedAvg3d(ROOTROT,ican,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    PAMNROT = weightedAvg3d(PAMNROT,ican,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    PAMXROT = weightedAvg3d(PAMXROT,ican,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    LNZ0ROT = weightedAvg3d(LNZ0ROT,icp1,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    CMASROT = weightedAvg3d(CMASROT,ican,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)

    !these are not in init file but are setup in modelStatedrivers so they are averaged here as well except for one or two which are reset
    TBASROT = weightedAvg2d(TBASROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    CMAIROT = weightedAvg2d(CMAIROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    WSNOROT = weightedAvg2d(WSNOROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    ZSNLROT = 0.10
    ! ZSNLROT = 0.01
    TSFSROT = weightedAvg3d(TSFSROT,4,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    TACROT = weightedAvg2d(TACROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    QACROT = weightedAvg2d(QACROT,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    ITCTROT = 0

    !then these values can be easily calculated
    gleafmasrow = gleafmas_NSrow + gleafmassrow
    stemmassrow = stemmass_NSrow + stemmasssrow
    rootmassrow = rootmass_NSrow + rootmasssrow
    flhrlossrow = flhrloss_nsrow + flhrloss_srow

    !> if tracers are on we also copy all the related variables 
    if (useTracer > 0) then
        tracerGLeafMassrot = weightedAvg3d(tracerGLeafMassrot,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        tracerBLeafMassrot = weightedAvg3d(tracerBLeafMassrot,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        tracerStemMassrot = weightedAvg3d(tracerStemMassrot,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        tracerRootMassrot = weightedAvg3d(tracerRootMassrot,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        tracerLitrMassrot = weightedAvg4d(tracerLitrMassrot,iccp2,ignd,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        tracerSoilCMassrot = weightedAvg4d(tracerSoilCMassrot,iccp2,ignd,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        tracerMossCMassrot = weightedAvg2d(tracerSoilCMassrot,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        tracerMossLitrMassrot = weightedAvg2d(tracerMossLitrMassrot,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
    end if

    !> if N cycling is on we need to also copy all the related variables
    if (Ncycle_on) then
        !these are all averaged
        ngleafmas_NSrow = weightedAvg3d(ngleafmas_NSrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        ngleafmassrow = weightedAvg3d(ngleafmassrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        nbleafmasrow = weightedAvg3d(nbleafmasrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
        nstemmass_NSrow = weightedAvg3d(nstemmass_NSrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        nstemmasssrow = weightedAvg3d(nstemmasssrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        nrootmass_NSrow = weightedAvg3d(nrootmass_NSrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        nrootmasssrow = weightedAvg3d(nrootmasssrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        bnfnatrow = weightedAvg3d(bnfnatrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        bnfantrow = weightedAvg3d(bnfantrow,icc,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        nh4_massrow = weightedAvg3d(nh4_massrow,iccp2,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        no3_massrow = weightedAvg3d(no3_massrow,iccp2,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        nlitrmassrow = weightedAvg3d(nlitrmassrow,iccp2,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        soilnmasrow = weightedAvg3d(soilnmasrow,iccp2,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode) 
        !ph comes from the head tile
        do i=1,nlat
            if (mode(i) /= "skip") then
                soilpHrow(i,outputIndex(i)) = soilpHrow(i,headTile(i))
            end if
        end do
        ! green leaf mass, stem mass, and root mass can just be calculated now
        ngleafmasrow = ngleafmas_NSrow + ngleafmassrow
        nstemmassrow = nstemmass_NSrow + nstemmasssrow
        nrootmassrow = nrootmass_NSrow + nrootmasssrow
    end if


    !this is some extra code that could be useful later on
    !if the ability to add competition were added these would need to be copied
    !they are assigned values from the first tile in modelStateDrivers.f90 anyway
    !if (PFTCompetition) then
    !   twarmm(:,outputIndex(:),:) = twarmm(:,headTile(:),:)
    !   tcoldm(:,outputIndex(:),:) = tcoldm(:,headTile(:),:)
    !   gdd5(:,outputIndex(:),:) = gdd5(:,headTile(:),:)
    !   aridity(:,outputIndex(:),:) = aridity(:,headTile(:),:)
    !   srplsmon(:,outputIndex(:),:) = srplsmon(:,headTile(:),:)
    !   defctmon(:,outputIndex(:),:) = defctmon(:,headTile(:),:)
    !   anndefct(:,outputIndex(:),:) = anndefct(:,headTile(:),:)
    !   annsrpls(:,outputIndex(:),:) = annsrpls(:,headTile(:),:)
    !   annpcp(:,outputIndex(:),:) = annpcp(:,headTile(:),:)
    !   dry_season_length(:,outputIndex(:),:) = dry_season_length(:,headTile(:),:)
    !end if
    
    ! If the ability to do fire and competition was added in these would also need to be coppied
    !if (dofire .and. PFTCompetition) then
    !   pstemmassrow(:,outputIndex(:),:) = pstemmassrow(:,headTile(:),:)
    !  	pgleafmassrow(:,outputIndex(:),:) = pgleafmassrow(:,headTile(:),:)
    !end if

!change the fractional areas of the input and output tiles
    do j=1,nlat
        if (mode(j) /= "skip") then
            FAREROT(j,inputIndex(j,:)) = FAREINPUT(j,inputIndex(j,:))  
            FAREROT(j,outputIndex(j)) = FAREFINAL(j)
        end if
    end do

    do j=1,nlat
        if (mode(j) /= "skip") then
            !re-check the fractional area of all tiles
            if ((sum(FAREROT(j,:))-1.0 > dynTilingTolrance) .or. (1.0-sum(FAREROT(j,:)) > dynTilingTolrance)) then
                write(6, * ) 'moveSplitCopyTiles: after copying the fractional area of the input tiles is different'
                write(6, * ) 'mode = ', mode
                write(6, * ) 'test threshold = ', dynTilingTolrance
                write(6, * ) 'test value = ', sum(FAREROT(j,:))-1
                !write(6, * ) 'farerot =', FAREROT(j,:)
                call errorHandler('moveSplitCopyTiles', -12)
            end if
        end if
    end do

    !> then we perform the final mass balance checks which ensure that no carbon was lost or gained
    call FinalMassEnergyBalance(FAREROT, DLZWROT, fcancmxrow, gleafmasrow, bleafmasrow, stemmassrow, rootmassrow, Cmossmasrow, peatSoilCrow, litrmsmossrow, &
                                soilcmasrow, litrmassrow, nh4_massrow, no3_massrow, nlitrmassrow, soilnmasrow, &
                                ngleafmasrow, nbleafmasrow ,nstemmassrow, nrootmassrow ,RCANROT ,SCANROT, TCANROT, TSNOROT, SNOROT, WSNOROT, TBARROT, &
                                TPNDROT, THICROT, THLQROT, ZPNDROT, vegCMass_initial, soilCMass_initial, peatCMass_initial, lucMass_initial, totalCMass_initial, &
                                soilWaterMass_initial, csWaterMass_initial, totalWaterMass_initial, vegNMass_initial, soilNMass_initial, lucNMass_initial, totalNMass_initial, TVSTP_initial, & 
                                TSSTP_initial, TSTP_initial, TSPN_initial, TTOT_initial, Ncycle_on, ipeatlandrow, DELZ, THPROT, HCPSROT, CMAIROT, mode)

    !> Then at the end we set a flag so the model gets re-initialized upon re-entry to mainCore.f90
    !> this cuts down on the number of variables that need to be transferred from tile to tile.                        
    DynTilInitializeFlag = .True.

  end subroutine moveSplitCopyTiles
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_cleanTiles
  !! @{
  !> Routine which simplified tiles based on a vegetation height threshold and when space is needed. This is the second main feature in this toolkit (see the module-level description for an overview)
  !> @author S.R. Curasi

  subroutine cleanTiles(controlVector,tneeded,DynTilInitializeFlag,outputIndex, & !these are the main intputs that control the routine
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
    
    use, intrinsic :: iso_fortran_env, only: r8=>real64

    implicit none

    !these are the switches
    logical, intent(in) :: Ncycle_on
    integer, intent(in) :: useTracer

    !setup the intputs for class state vars
    real, intent(inout), dimension(nlat,nmos,ignd) :: DLZWROT !< Permeable thickness of soil layer [m]
    real, intent(in), dimension(ignd) :: DELZ !< Overall thickness of soil layer [m]
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
    real, intent(inout), dimension(nlat,nmos,ignd) :: HCPSROT !< soil heat capacity
    integer, intent(inout), dimension(nlat,nmos,6,50) :: ITCTROT !<
    real, intent(inout), dimension(nlat,nmos,4)  :: TSFSROT !<
    real, intent(inout), dimension(nlat,nmos) :: maxAnnualActLyrROT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, intent(inout), dimension(nlat,nmos,ignd) :: THPROT !< soil pore volume
    real, intent(inout), dimension(nlat,nmos) :: tileAgerow 


    !setup the intputs for ctem state vars
    real, intent(inout), dimension(nlat,nmos) :: litrmsmossrow
    real, intent(inout), dimension(nlat,nmos) :: Cmossmasrow
    real, intent(inout), dimension(nlat,nmos) :: dmossrow
    real, intent(inout), dimension(nlat,nmos) :: peatSoilCrow !< peat soil C mass, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,8)  :: slopefracrow  !< prescribed fraction of wetlands based on slope
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

    !declare the tracer pools variables
    real, intent(inout), dimension(nlat,nmos) :: tracerMossCMassrot !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos) :: tracerMossLitrMassrot   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerGLeafMassrot  !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerBLeafMassrot  !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerStemMassrot   !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(inout), dimension(nlat,nmos,icc) :: tracerRootMassrot   !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(inout), dimension(nlat,nmos,iccp2,ignd) :: tracerLitrMassrot   !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout), dimension(nlat,nmos,iccp2,ignd) :: tracerSoilCMassrot  !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    !declare input variables used specifically by the subroutine
    integer, dimension(nlat) :: outputIndex !< the index of the tile where the output will be written
    real, intent(inout), dimension(nlat,nmos,icc) :: veghghtrow      !< vegetation height(meters)
    real, intent(inout), dimension(nlat) :: tneeded !< the number of tiles needed to continue running the model
    logical, intent(inout) :: DynTilInitializeFlag !< a flag that's set if the model needs to be re-initialized
    integer, intent(in), dimension(nlat,nmos) :: controlVector !< the control vector specified in the initializtion file

    real, dimension(nlat) :: FAREAROUT !< a real value greater than zero and less than one which specifies the desired size of the output tile
    real, dimension(nlat,nmos) :: veghght_avg      !< average vegetation height(meters)
    real, dimension(nlat,nmos,nmos) :: veghght_diff      !< the difference in vegetation height(meters)
    logical, dimension(nlat,nmos,nmos) :: logical_diff      !< masks out tiles that are not valid for a join
    real, dimension(nlat) :: cmin !< The minimum difference in vegetation height
    real, dimension(nlat) :: cjoined !< the number of tiles that have been joined
    character(len=5) :: mode(nlat) !< a character string that specifies the operation (i.e. move, split, copy, or skip) 
    integer, dimension(nlat,nmos) :: inputVector !< the input vector (1 = operate on tile, 0 = do not operate on tile) 
    integer, dimension(nlat,nmos) :: piv1 !< one of the two potential input vectors
    integer, dimension(nlat,nmos) :: piv2 !< one of the two potential input vectors
    integer, dimension(2,nlat) :: cindex !< The index of the tiles with the minimum difference in vegetation height
    integer, dimension(1,nlat) :: hmid !< The index of the tiles with the minimum difference in vegetation height
    integer :: i, k, j, l, k2(nlat), j2(nlat)
    real, dimension(nlat, nmos) :: sortedHeights
    real, dimension(nlat, nTilepres) :: subsetHeights
    integer, dimension(nlat) :: actMax !< the index of the right most active tile (used to shift tiles left at the end) 
    integer, dimension(nlat) :: inactMax !< the index of the right most inactive tiles surrounded by active tiles (used to shift tiles left at the end)
    integer, dimension(nlat) :: mv_cnt !< counter used to track the number of tile that get moved

    !> first we check there are enough tiles
    do i =1,nlat
        if ((tneeded(i) + 1) > nmos) then
            write(6, * ) 'cleanTiles: not enough tiles to run this model scenario'
            call errorHandler('cleanTiles', -1)
        end if
    end do

    !> first we calculate average vegetation height in each tile
    veghght_avg(:,:) = 0.0

    do i =1,nlat
        do j = 1,nmos
            do k = 1,icc
                veghght_avg(i,j) = veghght_avg(i,j) + (veghghtrow(i,j,k)*fcancmxrow(i,j,k))
            end do 
            veghght_avg(i,j) = veghght_avg(i,j)/sum(fcancmxrow(i,j,:))
        end do 
    end do 

    !> then we calculate the difference in average height and mark tiles which can't be joined
    do i =1,nlat
        do j = 1,nmos
            do k = 1,nmos
                ! Note the j and k loops are both over nmos, this creates a matrix to compare all tiles.
                veghght_diff(i,j,k) = ABS(veghght_avg(i,j)-veghght_avg(i,k))
                hmid(:,i) = MAXLOC(veghght_avg(i,:),(.not.(FAREROT(i,:) <= 0.0 .or. controlVector(i,:) == 0 .or. ABS(controlVector(i,:)) /= ABS(controlVector(i,j)))))

                if(veghght_avg(i,hmid(1,i)) <= 0.0) then
                    veghght_diff(i,j,k) = 0.0
                else
                    veghght_diff(i,j,k) = (veghght_diff(i,j,k)/veghght_avg(i,hmid(1,i)))
                end if

                logical_diff(i,j,k) = (.not. (j == k .or. FAREROT(i,j) <= 0.0 .or. FAREROT(i,k) <= 0.0 .or. &
                                        controlVector(i,j) == 0 .or. controlVector(i,k) == 0 .or. &
                                        (ABS(controlVector(i,k)) /= ABS(controlVector(i,j)))))
                if (.not. logical_diff(i,j,k)) veghght_diff(i,j,k) = 0.0
            end do 
        end do 

        !> if nTilespres is set, the smallest n tiles (nTilepres) will be protected
        if(nTilepres > 0) then
            where (.not. (FAREROT(i,:) <= 0.0 .or. controlVector(i,:) == 0))
                sortedHeights(i,:) = veghght_avg(i,:)
            elsewhere
                sortedHeights(i,:) = 9999.0
            end where

            sortedHeights(i,:) = indexSort(sortedHeights(i,:),nmos)
            sortedHeights(i,:) = sortedHeights(i,size(sortedHeights):1:-1)
            subsetHeights(i,:) = sortedHeights(i,1:nTilepres)

            logical_diff(i,:,subsetHeights(i,:)) = .FALSE.
            logical_diff(i,subsetHeights(i,:),:) = .FALSE.

        end if
    end do

    !> now we find the index for the tile with the minimum difference in vegetation height
    do i = 1,nlat
        cindex(:,i) = MINLOC(veghght_diff(i,:,:),logical_diff(i,:,:))
        j2(i) = cindex(1,i)
        k2(i) = cindex(2,i)
        cmin(i) = veghght_diff(i,j2(i),k2(i)) ! set the min diff in veg height
        cjoined(i) = 0 ! set the number of tiles joined to 0
    end do
    
    !> now we enter the main loop. This joins tiles until there are enough for the model to 
    !! continue running and all tiles below the height threshold have been joined.
    !! If there aren’t enough active tiles consumed so far we skip this operation
    wloop5: do while((ANY(cmin < percentthreshold) .or. ANY(cjoined < tneeded)) .and. &
                     ANY(count((FAREROT(:,:) > 0.0 .and. controlVector(:,:) /= 0),dim=2) > (nTilepres + 1)))  

        do i = 1,nlat

            if (ALL(.not. logical_diff(i,:,:)) .and. count((FAREROT(i,:) > 0.0 .and. controlVector(i,:) /= 0)) <= (nTilepres + 1)) then
                !> in the case where we have multiple latitude bands we might need to skip one
                mode(i) = "skip"

            else if (ALL(.not. logical_diff(i,:,:)) .and. count((FAREROT(i,:) > 0.0 .and. controlVector(i,:) /= 0)) > (nTilepres + 1)) then
                write(6, * ) 'cleanTiles: issue selecting input tiles no logical joins'
                write(6, * ) 'vh = ', veghght_diff(:,:,:)
                write(6, * ) 'ld = ', logical_diff(:,:,:)
                write(6, * ) 'FAREA = ', FAREROT(:,:)
                write(6, * ) 'ct = ', count(FAREROT(:,:)>0)
                call errorHandler('cleanTiles', -1)
            
            else
                ! If additional free tiles are needed to continue the run (cjoined(i) < tneeded(i)), but
                ! all the differences in vegetation height exceed the threshold for an automatic join
                ! (veghght_diff(i,j2(i),k2(i)) >= percentthreshold) carry our a forced join of the two
                ! least different tiles so the model can continue running.
                if (veghght_diff(i,j2(i),k2(i)) >= percentthreshold .and. cjoined(i) < tneeded(i)) then
                    inputVector(i,:) = 0
                    FAREAROUT(i) = -9999.
                    mode(i) = "merge"

                    inputVector(i,cindex(:,i)) = 1

                    !write(6, * ) '########################################'
                    !write(6, * ) 'cleanTiles: warning forced join required'
                    !write(6, * ) 'fcancmxrow = ', fcancmxrow(i,:,:)
                    !write(6, * ) 'veghghtrow = ', veghghtrow(i,:,:)
                    !write(6, * ) 'inputVector = ', inputVector(i,:)
                    !write(6, * ) 'cmin = ', cmin(i)
                    !write(6, * ) 'veghght_avg = ', veghght_avg(i,:)
                    !write(6, * ) 'tileage_avg = ', tileAgerow(i,:)
                    !write(6, * ) '########################################'

                ! If a pair of tiles falls below the threshold for an automatic join
                ! (veghght_diff(i,j2(i),k2(i)) >= percentthreshold) join them preemptively.
                else if (veghght_diff(i,j2(i),k2(i)) < percentthreshold) then
                    piv1(i,:) = 0 
                    piv2(i,:) = 0 
                    inputVector(i,:) = 0 
                    FAREAROUT(i) = -9999.0
                    mode(i) = "merge" 

                    do l=1,nmos
                        !< we traverse the difference matrix in the x and y to determine the maximum number of possible joins
                        if((veghght_diff(i,j2(i),l) < percentthreshold) .and. (logical_diff(i,j2(i),l))) then
                            piv1(i,l) = 1
                            piv1(i,j2(i)) = 1
                        end if
                        if((veghght_diff(i,l,k2(i)) < percentthreshold) .and. (logical_diff(i,l,k2(i)))) then
                            piv2(i,l) = 1
                            piv2(i,k2(i)) = 1
                        end if
                    end do	

                    !< then we setup the input vector to join the maximum number of possible tiles
                    if (sum(piv1(i,:)) < sum(piv2(i,:))) then
                        inputVector(i,:) = piv2(i,:)
                    else
                        inputVector(i,:) = piv1(i,:)
                    end if
                
                !If the model is running multiple gridcells (this never happes offline and has not been fully tested)
                !skip the grid cells where no dynamic tiling operations are required
                else if (veghght_diff(i,j2(i),k2(i)) >= percentthreshold .and. cjoined(i) >= tneeded(i) .and. (ANY(cmin < percentthreshold) .or. ANY(cjoined < tneeded))) then
                    mode(i) = "skip"
                !If something odd has happened and no case is selected error out
                else
                    write(6, * ) 'cleanTiles: issue selecting input tiles'
                    call errorHandler('cleanTiles', -2)
                    exit wloop5
                end if
            end if
        end do

        !print *,"m= ", mode, " iv= ",inputVector, "far= ",FAREAROUT, "vhd= ",veghght_diff,"ld= ",logical_diff,"vhm= ",cmin,"cim= ",cindex,"hth= ",heighthresholdIn, "cj= ",cjoined, "tn= ",tneeded
                
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

        !> now we update vegetation height
        call updateHeight(ipeatlandrow, veghghtrow, stemmassrow, gleafmasrow, bleafmasrow)

        !> we re-calculate calculate average vegetation height in each tile
        veghght_avg(:,:) = 0.0

        do i =1,nlat
                do j = 1,nmos
                    do k = 1,icc
                        veghght_avg(i,j) = veghght_avg(i,j) + (veghghtrow(i,j,k)*fcancmxrow(i,j,k))
                    end do 
                    veghght_avg(i,j) = veghght_avg(i,j)/sum(fcancmxrow(i,j,:))
                end do 
        end do 

        !> then we calculate the difference in average height and mark tiles which can't be joined
        do i =1,nlat
                do j = 1,nmos
                    do k = 1,nmos
                        veghght_diff(i,j,k) = ABS(veghght_avg(i,j)-veghght_avg(i,k))
                        hmid(:,i) = MAXLOC(veghght_avg(i,:),(.not.(FAREROT(i,:) <= 0.0 .or. controlVector(i,:) == 0 .or. ABS(controlVector(i,:)) /= ABS(controlVector(i,j)))))
                        veghght_diff(i,j,k) = (veghght_diff(i,j,k)/veghght_avg(i,hmid(1,i)))
                        logical_diff(i,j,k) = (.not. (j == k .or. FAREROT(i,j) <= 0.0 .or. FAREROT(i,k) <= 0.0 .or. controlVector(i,j) == 0 .or. controlVector(i,k) == 0 .or. (ABS(controlVector(i,k)) /= ABS(controlVector(i,j)))))
                        if (.not. logical_diff(i,j,k)) veghght_diff(i,j,k) = 0.0
                    end do 
                end do 
            
            !> if this is set the smallest n tiles will be protected
            if(nTilepres > 0) then
                where (.not. (FAREROT(i,:) <= 0.0 .or. controlVector(i,:) == 0))
                    sortedHeights(i,:) = veghght_avg(i,:)
                elsewhere
                    sortedHeights(i,:) = 9999.0
                end where

                sortedHeights(i,:) = indexSort(sortedHeights(i,:),nmos)
                sortedHeights(i,:) = sortedHeights(i,size(sortedHeights):1:-1)
                subsetHeights(i,:) = sortedHeights(i,1:nTilepres)
                logical_diff(i,:,subsetHeights(i,:)) = .FALSE.
                logical_diff(i,subsetHeights(i,:),:) = .FALSE.
            end if
        end do

        !> now we find the index for the tile with the minimum difference in vegetation height
        do i = 1,nlat
            cindex(:,i) = MINLOC(veghght_diff(i,:,:),logical_diff(i,:,:))
            j2(i) = cindex(1,i)
            k2(i) = cindex(2,i)
            cmin(i) = veghght_diff(i,j2(i),k2(i))
            !sometimes two tiles will be joined in a single operation and this needs to be accounted for
            if(mode(i) == "merge") then
                cjoined(i) = cjoined(i) + (sum(inputVector(i,:)) - 1)
            else if (mode(i) /= "skip" .and. mode(i) /= "merge") then
                write(6, * ) 'cleanTiles: illegal mode argument'
                call errorHandler('cleanTiles', -5)
                exit wloop5
            else
                write(6, * ) 'cleanTiles: illegal mode argument'
                call errorHandler('cleanTiles', -6)
                exit wloop5
            end if
        end do

    end do wloop5

    !< Now we shift all the tiles left to clean things up if needed (this rarely runs)
    do i =1,nlat
        actMax(i) = FINDLOC(FAREROT(i,:) > 0.0,VALUE = .true.,DIM = 1,BACK = .true.)
        if(ANY(FAREROT(i,1:actMax(i)) == 0.0)) then
            inactMax(i) = FINDLOC(FAREROT(i,1:actMax(i)) == 0.0,VALUE = .true.,DIM = 1,BACK = .false.)
        else
            inactMax(i) = nmos
        end if
    end do

    mv_cnt(:) = 0

    wloop6: do while(ANY(inactMax(:) < actMax(:)))
        do i =1,nlat

            if(mv_cnt(i) > nmos) then
                write(6, * ) 'too many attempts to left shift'
                !print *, "move left error: inact = ",inactMax(i)," act = ",actMax(i), " vals = ",FAREROT(i,:), " i = ",i,"mv count = ",mv_cnt(i), "controlVector(i,:) = ", controlVector(i,:)
                call errorHandler('cleanTiles', -7)
                exit wloop6
            end if

            !use to verify these left shifts
            !print *, "move left requred: inact = ",inactMax(i)," act = ",actMax(i), " vals = ",FAREROT(i,:), " i = ",i,"mv count = ",mv_cnt(i)  
            inputVector(i,:) = 0
            FAREAROUT(i) = -9999.0
            if(ANY(inactMax(:) < actMax(:))) then
                mode(i) = "move"
                inputVector(i,actMax(i)) = 1
                mv_cnt(i) = mv_cnt(i) + 1
            else
                mode(i) = "skip"
            end if
       end do
        
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

        do i =1,nlat
            actMax(i) = FINDLOC(FAREROT(i,:) > 0.0,VALUE = .true.,DIM = 1,BACK = .true.)
            if(ANY(FAREROT(i,1:actMax(i)) == 0.0)) then
                inactMax(i) = FINDLOC(FAREROT(i,1:actMax(i)) == 0.0,VALUE = .true.,DIM = 1,BACK = .true.)
            else
                inactMax(i) = nmos
            end if
            !print *, "move left w update: inact = ",inactMax(i)," act = ",actMax(i), " vals = ",FAREROT(i,:), " i = ",i,"mv count = ",mv_cnt(i) 
        end do 
    end do wloop6

  end subroutine cleanTiles
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  ! ------------------------------------------------------------------

  !> \ingroup dynamicTiling_indexTileAge
  !! @{
  !> This subroutine increases the age of the tile at each CTEM timestep
  !> @author S.R. Curasi
  subroutine indexTileAge(FAREGAT,tileAgegat, il1, il2, leapnow)

    integer, intent(in) :: il1   !< il1=1
    integer, intent(in) :: il2   !< il2=ilg
    real, intent(in), dimension(ilg) :: faregat
    real, intent(inout), dimension(ilg) :: tileAgegat
    logical, intent(in) :: leapnow   !< true if this year is a leap year. Only used if the switch 'leap' is true.
    integer :: i

    do i = il1,il2
        if (FAREGAT(i) > 0.0) then
            if (leapnow) then
                tileAgegat(i) = tileAgegat(i) + (deltat/(366./12.0))
            else ! non-leap year
                tileAgegat(i) = tileAgegat(i) + (deltat/(365./12.0))
            end if 
      end if
    end do

  end subroutine indexTileAge
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_initialMassEnergyBalance
  !! @{
  !> This function calculates the initial total mass of carbon, water, etc in a grid cell under consideration
  !> @author S.R. Curasi
  subroutine initialMassEnergyBalance(FAREROT, DLZWROT, fcancmxrow, gleafmasrow, bleafmasrow, stemmassrow, rootmassrow, Cmossmasrow, peatSoilCrow, litrmsmossrow, &
                                soilcmasrow, litrmassrow, nh4_massrow, no3_massrow, nlitrmassrow, soilnmasrow, &
                                ngleafmasrow, nbleafmasrow ,nstemmassrow, nrootmassrow ,RCANROT ,SCANROT, TCANROT, TSNOROT, SNOROT, WSNOROT, TBARROT, &
                                TPNDROT, THICROT, THLQROT, ZPNDROT, vegCMass_initial, soilCMass_initial, peatCMass_initial, lucMass_initial, totalCMass_initial, &
                                soilWaterMass_initial, csWaterMass_initial, totalWaterMass_initial, vegNMass_initial, soilNMass_initial, lucNMass_initial, totalNMass_initial, TVSTP_initial, & 
                                TSSTP_initial, TSTP_initial, TSPN_initial, TTOT_initial, Ncycle_on, ipeatlandrow, DELZ, THPROT, HCPSROT, CMAIROT, mode)

    use, intrinsic :: iso_fortran_env, only: r8=>real64

    implicit none 

    !these are used or calculated within the function
    integer :: i, k, j, l
    character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split or copy)
    real, dimension(nlat,nmos) :: barefrac

    !these are important descriptors of the land surface taken by the function
    real, intent(in), dimension(nlat,nmos) :: FAREROT
    real, intent(in), dimension(nlat,nmos,ignd) :: DLZWROT
    real, intent(in), dimension(nlat,nmos,icc) :: fcancmxrow
    logical, intent(in) :: Ncycle_on
    integer, intent(in), dimension(nlat,nmos) :: ipeatlandrow
    real, intent(in), dimension(ignd) :: DELZ
    real, intent(in), dimension(nlat,nmos,ignd) :: THPROT
    real, intent(in), dimension(nlat,nmos,ignd) :: HCPSROT
    real, intent(in), dimension(nlat,nmos) :: CMAIROT !<

    !vegetation biomass related variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,icc) :: gleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: bleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: stemmassrow
    real, intent(in), dimension(nlat,nmos,icc) :: rootmassrow

    !peatland variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos) :: Cmossmasrow
    real, intent(in), dimension(nlat,nmos) :: peatSoilCrow         !< peat soil C mass, \f$kg C/m^2\f$
    real, intent(in), dimension(nlat,nmos) :: litrmsmossrow

    !soil and litter variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,iccp2,ignd) :: soilcmasrow
    real, intent(in), dimension(nlat,nmos,iccp2,ignd) :: litrmassrow

    !water storage related variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,ignd) :: THLQROT
    real, intent(in), dimension(nlat,nmos,ignd) :: THICROT
    real, intent(in), dimension(nlat,nmos) :: RCANROT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: SCANROT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: SNOROT  !< Mass of snow pack \f$[kg m^{-2}]\f$

    !ncycling related variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,iccp1) :: nh4_massrow
    real, intent(in), dimension(nlat,nmos,iccp1) :: no3_massrow
    real, intent(in), dimension(nlat,nmos,iccp2) :: nlitrmassrow
    real, intent(in), dimension(nlat,nmos,iccp2) :: soilnmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: ngleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: nbleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: nstemmassrow
    real, intent(in), dimension(nlat,nmos,icc) :: nrootmassrow

    !energy budget related variables taken in by mass balance checks
    real, intent(in), dimension(nlat,nmos) :: TCANROT !< Vegetation canopy temperature [K]
    real(r8), intent(in), dimension(nlat,nmos,ignd) :: TBARROT !< Temperature of soil layers [K]
    real, intent(in), dimension(nlat,nmos) :: TSNOROT !< Snowpack temperature [K]
    real, intent(in), dimension(nlat,nmos) :: WSNOROT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: TPNDROT !< Temperature of ponded water [K]
    real, intent(in), dimension(nlat,nmos) :: ZPNDROT !< Depth of ponded water [m]

    !aggregated outputs from mass balance checks
    real, intent(inout), dimension(nlat) :: vegCMass_initial !< temporary storage variable for the total mass of carbon in vegetation
    real, intent(inout), dimension(nlat) :: soilCMass_initial !< temporary storage variable for the total mass of carbon in soil and litter
    real, intent(inout), dimension(nlat) :: peatCMass_initial !< temporary storage variable for the total mass of carbon in the peatland pools
    real, intent(inout), dimension(nlat) :: lucMass_initial !< temporary storage variable for the total mass of carbon in the LUC pools
    real, intent(inout), dimension(nlat) :: totalCMass_initial !< temporary storage variable for the total mass of carbon in all pools

    real, intent(inout), dimension(nlat) :: soilWaterMass_initial !< temporary storage variable for the total mass of water in the soil
    real, intent(inout), dimension(nlat) :: csWaterMass_initial !< temporary storage variable for the total mass of water in the canopy and on the soil surface  
    real, intent(inout), dimension(nlat) :: totalWaterMass_initial !< temporary storage variable for the total mass of water in all pools

    real, intent(inout), dimension(nlat) :: vegNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, intent(inout), dimension(nlat) :: soilNMass_initial !< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, intent(inout), dimension(nlat) :: lucNMass_initial !< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, intent(inout), dimension(nlat) :: totalNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation

    real, intent(inout), dimension(nlat) :: TVSTP_initial !< internal energy of vegetation used for energy balance checks (J)
    real, intent(inout), dimension(nlat) :: TSSTP_initial !< internal energy of snow pack used for energy balance checks (J)
    real, intent(inout), dimension(nlat) :: TSTP_initial !< internal energy of soil used for energy balance checks (J) 
    real, intent(inout), dimension(nlat) :: TSPN_initial !< internal energy ponded water used for energy balance checks (J)
    real, intent(inout), dimension(nlat) :: TTOT_initial !< total internal energy of all pools used for energy balance checks (J)

    vegCMass_initial = 0.0
    soilCMass_initial = 0.0
    peatCMass_initial = 0.0
    lucMass_initial = 0.0
    totalCMass_initial = 0.0

    soilWaterMass_initial = 0.0
    csWaterMass_initial = 0.0
    totalWaterMass_initial = 0.0

    vegNMass_initial = 0.0
    soilNMass_initial = 0.0
    lucNMass_initial = 0.0
    totalNMass_initial = 0.0

    TVSTP_initial = 0.0
    TSSTP_initial = 0.0
    TSTP_initial = 0.0
    TSPN_initial = 0.0
    TTOT_initial = 0.0

    do i=1,nlat
        if (mode(i) /= "skip") then
            do j=1,nmos
                barefrac(i,j) = 1.0 - sum(fcancmxrow(i,j,:))
                csWaterMass_initial(i) = csWaterMass_initial(i) + ((RCANROT(i,j) + SCANROT(i,j) + SNOROT(i,j) + (ZPNDROT(i,j) * RHOW)) * FAREROT(i,j))

                if (ipeatlandrow(i,j) > 0) then
                    peatCMass_initial(i) = peatCMass_initial(i) + ((Cmossmasrow(i,j) + peatSoilCrow(i,j) + litrmsmossrow(i,j)) * FAREROT(i,j))
                end if 


                TVSTP_initial(i) = TVSTP_initial(i) + TCANROT(i,j) * FAREROT(i,j) * (CMAIROT(i,j) * SPHVEG + RCANROT(i,j) * SPHW + SCANROT(i,j) * SPHICE)
                TSSTP_initial(i) = TSSTP_initial(i) + TSNOROT(i,j) * FAREROT(i,j) * (HCPICE * (SNOROT(i,j) / RHOICE) + HCPW * (WSNOROT(i,j) / RHOW))
                TSPN_initial(i) = TSPN_initial(i) + (HCPW * ZPNDROT(i,j) * FAREROT(i,j) * TPNDROT(i,j))

                do k = 1,icc
                    vegCMass_initial(i) =  vegCMass_initial(i) + &
                                            (gleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j)) + &
                                            (bleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j)) + &
                                            (stemmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j)) + &
                                            (rootmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j))
                    if (Ncycle_on) then
                        vegNMass_initial(i) =  vegNMass_initial(i) + &
                                                (ngleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg) + &
                                                (nbleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg) + &
                                                (nstemmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg) + &
                                                (nrootmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg)
                    end if
                end do
                do k = 1,ignd
                    soilWaterMass_initial(i) = soilWaterMass_initial(i) + &
                                            (THLQROT(i,j,k) * FAREROT(i,j) * 1000. * DLZWROT(i,j,k)) + &
                                            (THICROT(i,j,k) * FAREROT(i,j) * 1000. * DLZWROT(i,j,k))
                    
                    TSTP_initial(I) = TSTP_initial(I) + TBARROT(i,j,k) * FAREROT(i,j) * &
                        ((HCPW * THLQROT(i,j,k) + HCPICE * THICROT(i,j,k) + HCPSROT(i,j,k) * (1.0 - THPROT(i,j,k))) * DLZWROT(i,j,k) + &
                        HCPSND * (DELZ(k) - DLZWROT(i,j,k)))
                    do l=1,icc
                        soilCMass_initial(i) = soilCMass_initial(i) + &
                                            (soilcmasrow(i,j,l,k) * fcancmxrow(i,j,l) * FAREROT(i,j)) + &
                                            (litrmassrow(i,j,l,k) * fcancmxrow(i,j,l) * FAREROT(i,j))                   
                    end do
                    soilCMass_initial(i) = soilCMass_initial(i) + &
                                            (soilcmasrow(i,j,iccp1,k) * barefrac(i,j) * FAREROT(i,j)) + &
                                            (litrmassrow(i,j,iccp1,k) * barefrac(i,j) * FAREROT(i,j))  
                end do

                lucMass_initial(i) = lucMass_initial(i) + &
                                    (litrmassrow(i,j,iccp2,1) * FAREROT(i,j)) + &
                                    (soilcmasrow(i,j,iccp2,1) * FAREROT(i,j))  

                if (Ncycle_on) then 
                    do l=1,icc
                        soilNMass_initial(i) = soilNMass_initial(i) + &
                                                (soilnmasrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg) + &
                                                (nlitrmassrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg) + &
                                                (nh4_massrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg) + &
                                                (no3_massrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg)   
                    end do
                        soilNMass_initial(i) = soilNMass_initial(i) + &
                                                (soilnmasrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg) + &
                                                (nlitrmassrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg) + &
                                                (nh4_massrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg) + &
                                                (no3_massrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg)
                        lucNMass_initial(i) = lucMass_initial(i) + &
                                                (soilnmasrow(i,j,iccp2) * FAREROT(i,j) * convertg2kg) + &
                                                (nlitrmassrow(i,j,iccp2) * FAREROT(i,j) * convertg2kg) 
                end if
            end do
        end if
    end do

    totalCMass_initial = vegCMass_initial + soilCMass_initial + peatCMass_initial + lucMass_initial
    totalWaterMass_initial = csWaterMass_initial + soilWaterMass_initial
    totalNMass_initial = soilNMass_initial + vegNMass_initial + lucNmass_initial
    TTOT_initial = TVSTP_initial + TSSTP_initial + TSTP_initial + TSPN_initial

  end subroutine !initialMassEnergyBalance
    !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_finalMassEnergyBalance
  !! @{
  !> This function calculates the final total mass of carbon, water, etc in a grid cell under consideration and compares them to the values from initialMassEnergyBalance
  !> @author S.R. Curasi
  subroutine finalMassEnergyBalance(FAREROT, DLZWROT, fcancmxrow, gleafmasrow, bleafmasrow, stemmassrow, rootmassrow, Cmossmasrow, peatSoilCrow, litrmsmossrow, &
                                soilcmasrow, litrmassrow, nh4_massrow, no3_massrow, nlitrmassrow, soilnmasrow, &
                                ngleafmasrow, nbleafmasrow ,nstemmassrow, nrootmassrow ,RCANROT ,SCANROT, TCANROT, TSNOROT, SNOROT, WSNOROT, TBARROT, &
                                TPNDROT, THICROT, THLQROT, ZPNDROT, vegCMass_initial, soilCMass_initial, peatCMass_initial, lucMass_initial, totalCMass_initial, &
                                soilWaterMass_initial, csWaterMass_initial, totalWaterMass_initial, vegNMass_initial, soilNMass_initial, lucNmass_initial, totalNMass_initial, TVSTP_initial, & 
                                TSSTP_initial, TSTP_initial, TSPN_initial, TTOT_initial, Ncycle_on, ipeatlandrow, DELZ, THPROT, HCPSROT, CMAIROT, mode)

    use, intrinsic :: iso_fortran_env, only: r8=>real64
    
    implicit none 

    !these are used or calcualted within the function
    integer :: i, k, j, l
    character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
    real, dimension(nlat,nmos) :: barefrac
    logical, intent(in) :: Ncycle_on
    integer, intent(in), dimension(nlat,nmos) :: ipeatlandrow
    real, intent(in), dimension(ignd) :: DELZ    !< Overall thickness of soil layer [m]
    real, intent(in), dimension(nlat,nmos,ignd) :: THPROT
    real, intent(in), dimension(nlat,nmos,ignd) :: HCPSROT
    real, intent(in), dimension(nlat,nmos) :: CMAIROT !<

    !these are important descriptors of the land surface taken by the function
    real, intent(in), dimension(nlat,nmos) :: FAREROT
    real, intent(in), dimension(nlat,nmos,ignd) :: DLZWROT
    real, intent(in), dimension(nlat,nmos,icc) :: fcancmxrow
    
    !vegetation biomass related variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,icc) :: gleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: bleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: stemmassrow
    real, intent(in), dimension(nlat,nmos,icc) :: rootmassrow

    !peatland variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos) :: Cmossmasrow
    real, intent(in), dimension(nlat,nmos) :: peatSoilCrow         !< peat soil C mass, \f$kg C/m^2\f$
    real, intent(in), dimension(nlat,nmos) :: litrmsmossrow

    !soil and litter variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,iccp2,ignd) :: soilcmasrow
    real, intent(in), dimension(nlat,nmos,iccp2,ignd) :: litrmassrow

    !water storage related variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,ignd) :: THLQROT
    real, intent(in), dimension(nlat,nmos,ignd) :: THICROT
    real, intent(in), dimension(nlat,nmos) :: RCANROT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: SCANROT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: SNOROT  !< Mass of snow pack \f$[kg m^{-2}]\f$

    !ncycling related variables taken by mass balance checks
    real, intent(in), dimension(nlat,nmos,iccp1) :: nh4_massrow
    real, intent(in), dimension(nlat,nmos,iccp1) :: no3_massrow
    real, intent(in), dimension(nlat,nmos,iccp2) :: nlitrmassrow
    real, intent(in), dimension(nlat,nmos,iccp2) :: soilnmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: ngleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: nbleafmasrow
    real, intent(in), dimension(nlat,nmos,icc) :: nstemmassrow
    real, intent(in), dimension(nlat,nmos,icc) :: nrootmassrow

    !energy budget related variables taken in by mass balance checks
    real, intent(in), dimension(nlat,nmos) :: TCANROT !< Vegetation canopy temperature [K]
    real(r8), intent(in), dimension(nlat,nmos,ignd) :: TBARROT !< Temperature of soil layers [K]
    real, intent(in), dimension(nlat,nmos) :: TSNOROT !< Snowpack temperature [K]
    real, intent(in), dimension(nlat,nmos) :: WSNOROT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: TPNDROT !< Temperature of ponded water [K]
    real, intent(in), dimension(nlat,nmos) :: ZPNDROT !< Depth of ponded water [m]

    !aggregated outputs from initial mass balance checks
    real, intent(in), dimension(nlat) :: vegCMass_initial !< temporary storage variable for the total mass of carbon in vegetation
    real, intent(in), dimension(nlat) :: soilCMass_initial !< temporary storage variable for the total mass of carbon in soil and litter
    real, intent(in), dimension(nlat) :: peatCMass_initial !< temporary storage variable for the total mass of carbon in the peatland pools
    real, intent(in), dimension(nlat) :: lucMass_initial !< temporary storage variable for the total mass of carbon in the LUC pools
    real, intent(in), dimension(nlat) :: totalCMass_initial !< temporary storage variable for the total mass of carbon in all pools

    real, intent(in), dimension(nlat) :: soilWaterMass_initial !< temporary storage variable for the total mass of water in the soil
    real, intent(in), dimension(nlat) :: csWaterMass_initial !< temporary storage variable for the total mass of water in the canopy and on the soil surface  
    real, intent(in), dimension(nlat) :: totalWaterMass_initial !< temporary storage variable for the total mass of water in all pools

    real, intent(in), dimension(nlat) :: vegNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation
    real, intent(in), dimension(nlat) :: soilNMass_initial !< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, intent(in), dimension(nlat) :: lucNMass_initial !< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, intent(in), dimension(nlat) :: totalNMass_initial !< temporary storage variable for the total mass of nitrogen in the vegetation

    real, intent(in), dimension(nlat) :: TVSTP_initial !< internal energy of vegetation used for energy balance checks (J)
    real, intent(in), dimension(nlat) :: TSSTP_initial !< internal energy of snow pack used for energy balance checks (J)
    real, intent(in), dimension(nlat) :: TSTP_initial !< internal energy of soil used for energy balance checks (J) 
    real, intent(in), dimension(nlat) :: TSPN_initial !< internal energy ponded water used for energy balance checks (J)
    real, intent(in), dimension(nlat) :: TTOT_initial !< total internal energy of all pools used for energy balance checks (J)

    !aggregated outputs from initial mass balance checks
    real, dimension(nlat) :: vegCMass_final!< temporary storage variable for the total mass of carbon in vegetation
    real, dimension(nlat) :: soilCMass_final!< temporary storage variable for the total mass of carbon in soil and litter
    real, dimension(nlat) :: peatCMass_final!< temporary storage variable for the total mass of carbon in the peatland pools
    real, dimension(nlat) :: lucMass_final!< temporary storage variable for the total mass of carbon in the LUC pools
    real, dimension(nlat) :: totalCMass_final!< temporary storage variable for the total mass of carbon in all pools

    real, dimension(nlat) :: soilWaterMass_final!< temporary storage variable for the total mass of water in the soil
    real, dimension(nlat) :: csWaterMass_final!< temporary storage variable for the total mass of water in the canopy and on the soil surface  
    real, dimension(nlat) :: totalWaterMass_final!< temporary storage variable for the total mass of water in all pools

    real, dimension(nlat) :: vegNMass_final!< temporary storage variable for the total mass of nitrogen in the vegetation
    real, dimension(nlat) :: soilNMass_final!< temporary storage variable for the total mass of nitrogen in the soil and litter pools
    real, dimension(nlat) :: lucNMass_final!< temporary storage variable for the total mass of nitrogen in the land use change pools
    real, dimension(nlat) :: totalNMass_final!< temporary storage variable for the total mass of nitrogen in the vegetation

    real, dimension(nlat) :: TVSTP_final!< internal energy of vegetation
    real, dimension(nlat) :: TSSTP_final!< internal energy of snow pack
    real, dimension(nlat) :: TSTP_final!< internal energy of soil 
    real, dimension(nlat) :: TSPN_final!< internal energy ponded water
    real, dimension(nlat) :: TTOT_final!< total internal energy of all pools

    vegCMass_final = 0.0
    soilCMass_final = 0.0
    peatCMass_final = 0.0
    lucMass_final = 0.0
    totalCMass_final = 0.0

    soilWaterMass_final = 0.0
    csWaterMass_final = 0.0
    totalWaterMass_final = 0.0

    vegNMass_final = 0.0
    soilNMass_final = 0.0
    lucNMass_final = 0.0
    totalNMass_final = 0.0

    TVSTP_final = 0.0
    TSSTP_final = 0.0
    TSTP_final = 0.0
    TSPN_final = 0.0
    TTOT_final = 0.0

    do i=1,nlat
        if (mode(i) /= "skip") then
            do j=1,nmos
                barefrac(i,j) = 1.0 - sum(fcancmxrow(i,j,:))
                csWaterMass_final(i) = csWaterMass_final(i) + ((RCANROT(i,j) + SCANROT(i,j) + SNOROT(i,j) + (ZPNDROT(i,j) * RHOW)) * FAREROT(i,j))

                if (ipeatlandrow(i,j) > 0) then
                    peatCMass_final(i) = peatCMass_final(i) + ((Cmossmasrow(i,j) + peatSoilCrow(i,j) + litrmsmossrow(i,j)) * FAREROT(i,j))
                end if 

                TVSTP_final(i) = TVSTP_final(i) + TCANROT(i,j) * FAREROT(i,j) * (CMAIROT(i,j) * SPHVEG + RCANROT(i,j) * SPHW + SCANROT(i,j) * SPHICE)
                TSSTP_final(i) = TSSTP_final(i) + TSNOROT(i,j) * FAREROT(i,j) * (HCPICE * (SNOROT(i,j) / RHOICE) + HCPW * (WSNOROT(i,j) / RHOW))
                TSPN_final(i) = TSPN_final(i) + (HCPW * ZPNDROT(i,j) * FAREROT(i,j) * TPNDROT(i,j))

                do k = 1,icc
                    vegCMass_final(i) =  vegCMass_final(i) + &
                                            (gleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j)) + &
                                            (bleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j)) + &
                                            (stemmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j)) + &
                                            (rootmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j))
                    if (Ncycle_on) then
                        vegNMass_final(i) =  vegNMass_final(i) + &
                                                (ngleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg) + &
                                                (nbleafmasrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg) + &
                                                (nstemmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg) + &
                                                (nrootmassrow(i,j,k) * fcancmxrow(i,j,k) * FAREROT(i,j) * convertg2kg)
                    end if
                end do
                do k = 1,ignd
                    soilWaterMass_final(i) = soilWaterMass_final(i) + &
                                            (THLQROT(i,j,k) * FAREROT(i,j) * 1000.0 * DLZWROT(i,j,k)) + &
                                            (THICROT(i,j,k) * FAREROT(i,j) * 1000.0 * DLZWROT(i,j,k))

                    TSTP_final(I) = TSTP_final(I) + TBARROT(i,j,k) * FAREROT(i,j) * &
                        ((HCPW * THLQROT(i,j,k) + HCPICE * THICROT(i,j,k) + HCPSROT(i,j,k) * (1.0 - THPROT(i,j,k))) * DLZWROT(i,j,k) + &
                        HCPSND * (DELZ(k) - DLZWROT(i,j,k)))
                    do l=1,icc
                        soilCMass_final(i) = soilCMass_final(i) + &
                                            (soilcmasrow(i,j,l,k) * fcancmxrow(i,j,l) * FAREROT(i,j)) + &
                                            (litrmassrow(i,j,l,k) * fcancmxrow(i,j,l) * FAREROT(i,j))                    
                    end do
                    soilCMass_final(i) = soilCMass_final(i) + &
                                            (soilcmasrow(i,j,iccp1,k) * barefrac(i,j) * FAREROT(i,j)) + &
                                            (litrmassrow(i,j,iccp1,k) * barefrac(i,j) * FAREROT(i,j))  
                end do

                lucMass_final(i) = lucMass_final(i) + &
                                    (litrmassrow(i,j,iccp2,1) * FAREROT(i,j)) + &
                                    (soilcmasrow(i,j,iccp2,1) * FAREROT(i,j))  
                                    
                if (Ncycle_on) then 
                    do l=1,icc
                        soilNMass_final(i) = soilNMass_final(i) + &
                                                (soilnmasrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg) + &
                                                (nlitrmassrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg) + &
                                                (nh4_massrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg) + &
                                                (no3_massrow(i,j,l) * fcancmxrow(i,j,l) * FAREROT(i,j) * convertg2kg)   
                    end do
                        soilNMass_final(i) = soilNMass_final(i) + &
                                                (soilnmasrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg) + &
                                                (nlitrmassrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg) + &
                                                (nh4_massrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg) + &
                                                (no3_massrow(i,j,iccp1) * barefrac(i,j) * FAREROT(i,j) * convertg2kg)
                        lucNMass_final(i) = lucMass_final(i) + &
                                                (soilnmasrow(i,j,iccp2) * FAREROT(i,j) * convertg2kg) + &
                                                (nlitrmassrow(i,j,iccp2) * FAREROT(i,j) * convertg2kg) 
                end if
            end do
        end if
    end do

    totalCMass_final = vegCMass_final + soilCMass_final + peatCMass_final + lucMass_final
    totalWaterMass_final = csWaterMass_final + soilWaterMass_final
    totalNMass_final = soilNMass_final + vegNMass_final
    TTOT_final = TVSTP_final + TSSTP_final + TSTP_final + TSPN_final

    do i=1,nlat
        if (mode(i) /= "skip") then

            if ((abs(vegCMass_final(i) - vegCMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(soilCMass_final(i) - soilCMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(peatCMass_final(i) - peatCMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(lucMass_final(i) - lucMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(totalCMass_final(i) - totalCMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(soilWaterMass_final(i) - soilWaterMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(csWaterMass_final(i) - csWaterMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(totalWaterMass_final(i) - totalWaterMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(vegNMass_final(i) - vegNMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(soilNMass_final(i) - soilNMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(totalNMass_final(i) - totalNMass_initial(i)) >= dynTilingTolrance) .or. &
                (abs(TVSTP_final(i) - TVSTP_initial(i)) >= dynTilingTolranceEnergy) .or. &
                (abs(TSSTP_final(i) - TSSTP_initial(i)) >= dynTilingTolranceEnergy) .or. &
                (abs(TSTP_final(i) - TSTP_initial(i)) >= dynTilingTolranceEnergy) .or. &
                (abs(TSPN_final(i) - TSPN_initial(i)) >= dynTilingTolranceEnergy) .or. &
                (abs(TTOT_final(i) - TTOT_initial(i)) >= dynTilingTolranceEnergy)) then

                print *, "error in finalMassEnergyBalance checks"
                print *, "---mass balence check results---"
                print *, "nlat = ",i
                print *, "FAREROT = ", FAREROT(i,:)
                print *, "---carbon"
                print *, 'vegCMass_initial = ', vegCMass_initial(i), ' vegCMass_final = ', vegCMass_final(i), 'diff = ', vegCMass_final(i) - vegCMass_initial(i)
                print *, 'soilCMass_initial = ', soilCMass_initial(i), ' soilCMass_final = ', soilCMass_final(i), 'diff = ', soilCMass_final(i) - soilCMass_initial(i)
                print *, 'peatCMass_initial = ', peatCMass_initial(i), ' peatCMass_final = ', peatCMass_final(i), 'diff = ', peatCMass_final(i) - peatCMass_initial(i)
                print *, 'lucMass_initial = ', lucMass_initial(i), ' lucMass_final = ', lucMass_final(i), 'diff = ', lucMass_final(i) - lucMass_initial(i)
                print *, 'totalCMass_initial = ', totalCMass_initial(i), ' totalCMass_final = ', totalCMass_final(i), 'diff = ', totalCMass_final(i) - totalCMass_initial(i)
                print *, "---water"
                print *, 'soilWaterMass_initial = ', soilWaterMass_initial(i), ' soilWaterMass_final = ', soilWaterMass_final(i), 'diff = ', soilWaterMass_final(i) - soilWaterMass_initial(i)
                print *, 'csWaterMass_initial = ', csWaterMass_initial(i), ' csWaterMass_final = ', csWaterMass_final(i), 'diff = ', csWaterMass_final(i) - csWaterMass_initial(i)
                print *, 'totalWaterMass_initial = ', totalWaterMass_initial(i), ' totalWaterMass_final = ', totalWaterMass_final(i), 'diff = ', totalWaterMass_final(i) - totalWaterMass_initial(i)
                if (Ncycle_on) then 
                    print *, "---nitrogen"
                    print *, 'vegNMass_initial = ', vegNMass_initial(i), ' vegNMass_final = ', vegNMass_final(i), 'diff = ', vegNMass_final(i) - vegNMass_initial(i)
                    print *, 'soilNMass_initial = ', soilNMass_initial(i), ' soilNMass_final = ', soilNMass_final(i), 'diff = ', soilNMass_final(i) - soilNMass_initial(i)
                    print *, 'lucNMass_initial = ', lucNMass_initial(i), ' lucNMass_final = ', lucNMass_final(i), 'diff = ', lucNMass_final(i) - lucNMass_initial(i)
                    print *, 'totalNMass_initial = ', totalNMass_initial(i), ' totalNMass_final = ', totalNMass_final(i), 'diff = ', totalNMass_final(i) - totalNMass_initial(i)
                end if
                print *, "---energy"
                print *, 'TVSTP_initial = ', TVSTP_initial(i), ' TVSTP_final = ', TVSTP_final(i), 'diff = ', TVSTP_final(i) - TVSTP_initial(i)
                print *, 'TSSTP_initial = ', TSSTP_initial(i), ' TSSTP_final = ', TSSTP_final(i), 'diff = ', TSSTP_final(i) - TSSTP_initial(i)
                print *, 'TSTP_initial = ', TSTP_initial(i), ' TSTP_final = ', TSTP_final(i), 'diff = ', TSTP_final(i) - TSTP_initial(i)
                print *, 'TSPN_initial = ', TSPN_initial(i), ' TSPN_final = ', TSPN_final(i), 'diff = ', TSPN_final(i) - TSPN_initial(i)
                print *, 'TTOT_initial = ', TTOT_initial(i), ' TTOT_final = ', TTOT_final(i), 'diff = ', TTOT_final(i) - TTOT_initial(i)
                print *, "--------------------------------"
                call errorHandler('finalMassEnergyBalance', - 1)
            end if
        end if
    end do

  end subroutine !finalMassEnergyBalance
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_AvgSoilRichmann
  !! @{
  !> This function calculates the temperatures for new tiles. It applies Richmann’s law as if the material in the target grid cells is being mixed.
  !> @author S.R. Curasi
  subroutine MixTempsRichmann(FAREROT, DLZWROT, fcancmxrow, DELZ, THPROT, HCPSROT, & 
                            THLQROT, THICROT, RCANROT, SCANROT, SNOROT, WSNOROT, ZPNDROT, TCANROT, TBARROT, &
                            TSNOROT, TPNDROT, CMAIROT, inputLength, outputIndex, inputIndex, mode)

    use, intrinsic :: iso_fortran_env, only: r8=>real64

    implicit none 

    !these are vegetation masses used taken by the function
    real, intent(in), dimension(nlat,nmos) :: CMAIROT !<

    !these are important descriptors of the land surface taken by the function
    real, intent(in), dimension(nlat,nmos) :: FAREROT
    real, intent(in), dimension(nlat,nmos,ignd) :: DLZWROT
    real, intent(in), dimension(nlat,nmos,icc) :: fcancmxrow
    real, intent(in), dimension(ignd) :: DELZ
    real, intent(in), dimension(nlat,nmos,ignd) :: THPROT
    real, intent(in), dimension(nlat,nmos,ignd) :: HCPSROT

    !water storage related variables taken by the function
    real, intent(in), dimension(nlat,nmos,ignd) :: THLQROT
    real, intent(in), dimension(nlat,nmos,ignd) :: THICROT
    real, intent(in), dimension(nlat,nmos) :: RCANROT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: SCANROT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: SNOROT  !< Mass of snow pack \f$[kg m^{-2}]\f$

    !energy budget related variables taken in by mass balance checks
    real, intent(in), dimension(nlat,nmos) :: WSNOROT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, intent(in), dimension(nlat,nmos) :: ZPNDROT !< Depth of ponded water [m]

    !temperatures which are eventually overwritted by the function
    real, intent(inout), dimension(nlat,nmos) :: TCANROT !< Vegetation canopy temperature [K]
    real(r8), intent(inout), dimension(nlat,nmos,ignd) :: TBARROT !< Temperature of soil layers [K]
    real, intent(inout), dimension(nlat,nmos) :: TSNOROT !< Snowpack temperature [K]
    real, intent(inout), dimension(nlat,nmos) :: TPNDROT !< Temperature of ponded water [K]

    !< variables from the main subroutine that are required
    integer, intent(in), dimension(nlat) :: inputLength
    integer, intent(in), dimension(nlat) :: outputIndex
    integer, intent(in), dimension(nlat,maxval(inputLength)) :: inputIndex
    character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 

    !< variables used to do calculations internal to the subroutine
    real, dimension(nlat,nmos) :: temp_TSNOROT
    real, dimension(nlat,nmos) :: temp_TCANROT
    real, dimension(nlat,nmos) :: temp_TPNDROT
    real(r8), dimension(nlat,nmos,ignd) :: temp_TBARROT

    real, dimension(nlat,nmos) :: weight_TSNOROT
    real, dimension(nlat,nmos) :: weight_TCANROT
    real, dimension(nlat,nmos) :: weight_TPNDROT
    real(r8), dimension(nlat,nmos,ignd) :: weight_TBARROT

    integer :: i, j, k
    
    temp_TSNOROT = 0.0
    temp_TCANROT = 0.0
    temp_TPNDROT = 0.0
    temp_TBARROT = 0.0

    weight_TSNOROT = 0.0
    weight_TCANROT = 0.0
    weight_TPNDROT = 0.0
    weight_TBARROT = 0.0


    do i=1,nlat
        if (mode(i) /= "skip") then
            do j=1,inputLength(i)

                !snow
                temp_TSNOROT(i,outputIndex(i)) = temp_TSNOROT(i,outputIndex(i)) + TSNOROT(i,inputIndex(i,j)) * FAREROT(i,inputIndex(i,j)) * &
                                                (HCPICE * (SNOROT(i,inputIndex(i,j)) / RHOICE) + HCPW * (WSNOROT(i,inputIndex(i,j)) / RHOW))
                weight_TSNOROT(i,outputIndex(i)) = weight_TSNOROT(i,outputIndex(i)) + FAREROT(i,inputIndex(i,j)) * (HCPICE * &
                                                (SNOROT(i,inputIndex(i,j)) / RHOICE) + HCPW * (WSNOROT(i,inputIndex(i,j)) / RHOW))

                !canopy
                temp_TCANROT(i,outputIndex(i)) = temp_TCANROT(i,outputIndex(i)) + TCANROT(i,inputIndex(i,j)) * FAREROT(i,inputIndex(i,j)) * &
                                                (CMAIROT(i,inputIndex(i,j)) * SPHVEG + RCANROT(i,inputIndex(i,j)) * SPHW + SCANROT(i,inputIndex(i,j)) * SPHICE)
                weight_TCANROT(i,outputIndex(i)) = weight_TCANROT(i,outputIndex(i)) + FAREROT(i,inputIndex(i,j)) * &
                                                (CMAIROT(i,inputIndex(i,j)) * SPHVEG + RCANROT(i,inputIndex(i,j)) * SPHW + SCANROT(i,inputIndex(i,j)) * SPHICE)

                !pond
                temp_TPNDROT(i,outputIndex(i)) = temp_TPNDROT(i,outputIndex(i)) + (HCPW * ZPNDROT(i,inputIndex(i,j)) *&
                                                 FAREROT(i,inputIndex(i,j)) * TPNDROT(i,inputIndex(i,j)))
                weight_TPNDROT(i,outputIndex(i)) = weight_TPNDROT(i,outputIndex(i)) + (HCPW * ZPNDROT(i,inputIndex(i,j)) * FAREROT(i,inputIndex(i,j)))
            
                !soil
                do k=1,ignd
                    temp_TBARROT(i,outputIndex(i),k) = temp_TBARROT(i,outputIndex(i),k) + TBARROT(i,inputIndex(i,j),k) * FAREROT(i,inputIndex(i,j)) * & 
                                                        ((HCPW * THLQROT(i,inputIndex(i,j),k) + HCPICE * THICROT(i,inputIndex(i,j),k) + &
                                                        HCPSROT(i,inputIndex(i,j),k) * (1.0 - THPROT(i,inputIndex(i,j),k))) * DLZWROT(i,inputIndex(i,j),k) + &
                                                        HCPSND * (DELZ(k) - DLZWROT(i,inputIndex(i,j),k)))
                    weight_TBARROT(i,outputIndex(i),k) = weight_TBARROT(i,outputIndex(i),k) + FAREROT(i,inputIndex(i,j)) * ((HCPW * &
                                                        THLQROT(i,inputIndex(i,j),k) + HCPICE * THICROT(i,inputIndex(i,j),k) + HCPSROT(i,inputIndex(i,j),k) *&
                                                        (1.0 - THPROT(i,inputIndex(i,j),k))) * DLZWROT(i,inputIndex(i,j),k) + &
                                                        HCPSND * (DELZ(k) - DLZWROT(i,inputIndex(i,j),k)))
                end do   
            end do

            ! The resulting values get applied. If the weight is zero, for consistency, we take the value representing the largest area.
            if(weight_TSNOROT(i,outputIndex(i)) /= 0.0) then
                TSNOROT(i,outputIndex(i)) = temp_TSNOROT(i,outputIndex(i))/weight_TSNOROT(i,outputIndex(i))
            else
                TSNOROT(i,:) = majorityarea1d(TSNOROT,i,nmos,inputLength,outputIndex,inputIndex,FAREROT,mode)
            end if 

            if(weight_TCANROT(i,outputIndex(i)) /= 0.0) then
                TCANROT(i,outputIndex(i)) = temp_TCANROT(i,outputIndex(i))/weight_TCANROT(i,outputIndex(i))
            else
                TCANROT(i,:) = majorityarea1d(TCANROT,i,nmos,inputLength,outputIndex,inputIndex,FAREROT,mode) 
            end if

            if(weight_TPNDROT(i,outputIndex(i)) /= 0.0) then
                TPNDROT(i,outputIndex(i)) = temp_TPNDROT(i,outputIndex(i))/weight_TPNDROT(i,outputIndex(i))
            else
                TPNDROT(i,:) = majorityarea1d(TPNDROT,i,nmos,inputLength,outputIndex,inputIndex,FAREROT,mode)
            end if 
            ! For TBARROT the weighted value is always used. The weight is always >0 because 
            ! there is always some minearl fraction (1.0 - THPROT(i,inputIndex(i,j),k)))
            TBARROT(i,outputIndex(i),:) = temp_TBARROT(i,outputIndex(i),:)/weight_TBARROT(i,outputIndex(i),:)

        end if
      end do
  end subroutine MixTempsRichmann
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_updateHeight
  !! @{
  !> This function updates the vegetation height of all tiles. This is used by cleanTiles to ensure vegetation height is updated between each comparison/merge
  !> @author S.R. Curasi
  subroutine updateHeight(ipeatlandrow, veghghtrow, stemmassrow, gleafmasrow, bleafmasrow)

    use classicParams, only : fracbofg, reindexPFTs, classpfts

    implicit none 

    integer, intent(inout), dimension(nlat,nmos) :: ipeatlandrow
    real, intent(inout), dimension(nlat,nmos,icc) :: veghghtrow
    real, intent(inout), dimension(nlat,nmos,icc) :: stemmassrow
    real, intent(inout), dimension(nlat,nmos,icc) :: gleafmasrow
    real, intent(inout), dimension(nlat,nmos,icc) :: bleafmasrow

    integer :: i, j, k, l

        do i =1,nlat
        do j = 1,nmos
            do k = 1,ican
                do l = reindexPFTs(k,1),reindexPFTs(k,2) ! loop 260
                    select case (classpfts(k))
                    case ('NdlTr','BdlTr') ! Trees
                        if (ipeatlandrow(i,j) /= 1 .and. ipeatlandrow(i,j) /= 2) then ! Uplands
                        veghghtrow(i,j,k) = min(10.0 * stemmassrow(i,j,k) ** 0.385,45.0)
                        else ! peatland trees have a different relation than normal. Max height 10 m.
                        veghghtrow(i,j,k) = min(3.0 * stemmassrow(i,j,k) ** 0.385,10.0)
                        end if
                    case ('BdlSh')
                        if (ipeatlandrow(i,j) /= 1 .and. ipeatlandrow(i,j) /= 2) then ! Uplands
                        veghghtrow(i,j,k) = min(4.0,0.25 * (stemmassrow(i,j,k) ** 0.2)) !Changed to shrubs like peatland, but they can be taller than 1m, shrub edit**
                        else ! peatland shrubs
                        veghghtrow(i,j,k) = min(1.0,0.25 * (stemmassrow(i,j,k) ** 0.2))
                        end if
                    case ('Crops') ! <Crops
                        veghghtrow(i,j,k) = 1.0 * (stemmassrow(i,j,k) + gleafmasrow(i,j,k)) ** 0.385
                    case ('Grass') ! <Grass
                        if (ipeatlandrow(i,j) /= 1 .and. ipeatlandrow(i,j) /= 2) then ! Uplands
                        veghghtrow(i,j,k) = 3.5 * (gleafmasrow(i,j,k) + fracbofg * bleafmasrow(i,j,k)) ** 0.50
                        else ! peatland grasses and sedges
                        veghghtrow(i,j,k) = min(1.0, (gleafmasrow(i,j,k) + fracbofg * bleafmasrow(i,j,k)) ** 0.3)
                        end if
                    case default
                        print * ,'Unknown CLASS PFT in allometry ',classpfts(j)
                        call errorHandler('dynamic tiling', - 4)
                    end select
                end do
            end do
        end do
        end do

  end subroutine updateHeight
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_majority3d
  !! @{
  !> This function is used by the tile copying subroutine to find the most common value within a vector
  !> @author S.R. Curasi
  function majority3d(val,dim,nlat,nmos,inputLength,outputIndex,inputIndex,mode)

      implicit none

      integer, intent(in) :: nlat,nmos,inputLength(nlat), outputIndex(nlat), inputIndex(nlat,maxval(inputLength)),dim !< variables from the main subroutine that are required
      integer, intent(in), dimension(nlat,nmos,dim) :: val !< the values to be operated on from the tiles in question
      real, dimension(nlat,nmos,dim) :: majority3d !< the output into which calculations are done
      integer :: max_freq = 0, ans = -1, curr_freq
      character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
      integer :: i, k, j, l

      majority3d = val
      
      do i =1,nlat
        if (mode(i) /= "skip") then
            do j = 1,dim
            max_freq = 0
            ans = -1
            if(inputLength(i) == 1) then
                ans = val(i,inputIndex(i,1),j)
            else if (inputLength(i) > 1) then
                    ans = val(i,inputIndex(i,1),j)
                    do k = 1,inputLength(i)
                        curr_freq = 1
                        do l=k+1,inputLength(i)
                            if (val(i,inputIndex(i,k),j) == val(i,inputIndex(i,l),j)) then
                                curr_freq = curr_freq + 1
                            end if
                            if (max_freq < curr_freq) then
                                max_freq = curr_freq
                                ans = val(i,inputIndex(i,k),j)
                            end if
                        end do
                    end do
            else
                print * ,'majority3d function: invalid inputs to the majority3d function'
                call errorHandler('majority3d', - 1)
            end if
            majority3d(i,outputIndex(i),j) = ans
            end do
        end if
      end do
  RETURN
  end function majority3d
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_majorityarea3d
  !! @{
  !> This function is used by the tile copying subroutine to find the value within a vector which represents the largest area
  !> @author S.R. Curasi
  function majorityarea3d(val,dim,nlat,nmos,inputLength,outputIndex,inputIndex,FAREROT,mode)

       implicit none

      integer, intent(in) :: nlat,nmos,inputLength(nlat), outputIndex(nlat), inputIndex(nlat,maxval(inputLength)),dim !< variables from the main subroutine that are required
      integer, intent(in), dimension(nlat,nmos,dim) :: val !< the values to be operated on from the tiles in question
      real, intent(in), dimension(nlat,nmos) :: FAREROT
      real, dimension(nlat,nmos,dim) :: majorityarea3d !< the output into which calculations are done
      integer :: max_freq = 0, ans = -1, curr_freq
      character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
      integer :: i, k, j, l

      majorityarea3d = val
      
      do i =1,nlat
        if (mode(i) /= "skip") then
            do j = 1,dim
            max_freq = 0
            ans = -1
            if(inputLength(i) == 1) then
                ans = val(i,inputIndex(i,1),j)
            else if (inputLength(i) > 1) then
                    ans = val(i,inputIndex(i,1),j)
                    do k = 1,inputLength(i)
                        curr_freq = FAREROT(i,inputIndex(i,k))
                        do l=k+1,inputLength(i)
                            if (val(i,inputIndex(i,k),j) == val(i,inputIndex(i,l),j)) then
                                curr_freq = curr_freq + FAREROT(i,inputIndex(i,k))
                            end if
                            if (max_freq < curr_freq) then
                                max_freq = curr_freq
                                ans = val(i,inputIndex(i,k),j)
                            end if
                        end do
                    end do
            else
                print * ,'majorityarea3d function: invalid inputs to the majorityarea3d function'
                call errorHandler('majorityarea3d', - 1)
            end if
            majorityarea3d(i,outputIndex(i),j) = ans
            end do
        end if
      end do
  RETURN
  end function majorityarea3d
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_majorityarea1d
  !! @{
  !> This function is used by the tile copying subroutine within MixTempsRichmann to find the value within a vector which 
  !> represents the largest area. It is unique in that it operates over latitudes one at a time
  !> @author S.R. Curasi
  function majorityarea1d(val,lat_in,nmos,inputLength,outputIndex,inputIndex,FAREROT,mode)

    implicit none

    integer, intent(in) :: lat_in
    integer, intent(in) :: nmos
    integer, intent(in), dimension(nlat) :: inputLength
    integer, intent(in), dimension(nlat) :: outputIndex
    integer, intent(in), dimension(nlat,maxval(inputLength)) :: inputIndex !< variables from the main subroutine that are required
    real, intent(in), dimension(nlat,nmos) :: val !< the values to be operated on from the tiles in question
    real, intent(in), dimension(nlat,nmos) :: FAREROT
    real, dimension(nmos) :: majorityarea1d !< the output into which calculations are done
    integer :: max_freq = 0
    integer :: curr_freq
    real :: ans = -1 !< initialize this to -1 (this is used as a check and guarantees strange code behavior will error our subsequent dynamic tiling operations)
    character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
    integer :: i, k, j, l
    
    i = lat_in

    majorityarea1d = val(i,:)

    max_freq = 0

    if (inputLength(i) == 1) then
        ans = val(i,inputIndex(i,1))
    else if (inputLength(i) > 1) then
            ans = val(i,inputIndex(i,1))
            do k = 1,inputLength(i)
                curr_freq = FAREROT(i,inputIndex(i,k))
                do l=k+1,inputLength(i)
                    if (val(i,inputIndex(i,k)) == val(i,inputIndex(i,l))) then
                        curr_freq = curr_freq + FAREROT(i,inputIndex(i,k))
                    end if
                    if (max_freq < curr_freq) then
                        max_freq = curr_freq
                        ans = val(i,inputIndex(i,k))
                    end if
                end do
            end do
    else
        print * ,'majorityarea1d function: invalid inputs to the majorityarea1d function'
        call errorHandler('majorityarea1d', - 1)
    end if

    if (ans == -1) then
        print * ,'majorityarea1d function: invalid output from majorityarea1d function'
        call errorHandler('majorityarea1d', - 2)
    end if

    majorityarea1d(outputIndex(i)) = ans
    
    RETURN

  end function majorityarea1d
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_majorityarea2dr8
  !! @{
  !> This function is used by the tile copying subroutine within MixTempsRichmann to find the value within a vector which 
  !> represents the largest area. It is unique in that it operates over latitudes one at a time. It accepts r8 reals
  !> @author S.R. Curasi
  function majorityarea2dr8(val,lat_in,dim,nmos,inputLength,outputIndex,inputIndex,FAREROT,mode)

      use, intrinsic :: iso_fortran_env, only: r8=>real64

      implicit none

      integer, intent(in) :: lat_in,nmos,inputLength(nlat), outputIndex(nlat), inputIndex(nlat,maxval(inputLength)),dim !< variables from the main subroutine that are required
      real(r8), intent(in), dimension(nlat,nmos,dim) :: val !< the values to be operated on from the tiles in question
      real, intent(in), dimension(nlat,nmos) :: FAREROT
      real(r8), dimension(nmos,dim) :: majorityarea2dr8 !< the output into which calculations are done
      integer :: max_freq = 0, curr_freq
      real(r8) :: ans = -1
      character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
      integer :: i, k, j, l

      i = lat_in

      majorityarea2dr8 = val(i,:,:)
      
        do j = 1,dim
            max_freq = 0
            ans = -1
            if(inputLength(i) == 1) then
                ans = val(i,inputIndex(i,1),j)
            else if (inputLength(i) > 1) then
                    ans = val(i,inputIndex(i,1),j)
                    do k = 1,inputLength(i)
                        curr_freq = FAREROT(i,inputIndex(i,k))
                        do l=k+1,inputLength(i)
                            if (val(i,inputIndex(i,k),j) == val(i,inputIndex(i,l),j)) then
                                curr_freq = curr_freq + FAREROT(i,inputIndex(i,k))
                            end if
                            if (max_freq < curr_freq) then
                                max_freq = curr_freq
                                ans = val(i,inputIndex(i,k),j)
                            end if
                        end do
                    end do
            else
                print * ,'majorityarea2dr8 function: invalid inputs to the majorityarea2dr8 function'
                call errorHandler('majorityarea2dr8', - 1)
            end if
            majorityarea2dr8(outputIndex(i),j) = ans
        end do
  RETURN
  end function majorityarea2dr8
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_weightedAvg2d
  !! @{
  !> This function calculates the weighted average of a two-dimensional variable (it takes the variable [val] as an input)
  !> @author S.R. Curasi
  function weightedAvg2d(val,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
      
      implicit none

      integer, intent(in) :: nlat,nmos, inputLength(nlat), headTile(nlat), outputIndex(nlat), inputIndex(nlat,maxval(inputLength)) !< variables from the main subroutine that are required
      real, intent(in), dimension(nlat,nmos) :: val !< the input variable to be averaged and returned
      real, intent(in), dimension(nlat,nmos) :: FAREROT
      real, dimension(nlat,nmos) :: weightedAvg2d !< the variable into which calculations are done
      character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
      integer :: i, j

      weightedAvg2d = val
      do i=1,nlat
        if (mode(i) /= "skip") then
            weightedAvg2d(i,outputIndex(i)) = 0.0
            do j=1,inputLength(i)
                weightedAvg2d(i,outputIndex(i)) =  weightedAvg2d(i,outputIndex(i)) + (val(i,inputIndex(i,j)) * FAREROT(i,inputIndex(i,j)))
            end do
            weightedAvg2d(i,outputIndex(i)) = weightedAvg2d(i,outputIndex(i))/sum(FAREROT(i,inputIndex(i,:)))
        end if
      end do
  RETURN
  end function weightedAvg2d
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_weightedAvg3d
  !! @{
  !> This function calculates the weighted average of a three-dimensional variable (it takes the variable [val] and one additional dimension [dim] as inputs)
  !> @author S.R. Curasi
  function weightedAvg3d(val,dim,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
      
      implicit none

      integer, intent(in) :: nlat,nmos,inputLength(nlat), headTile(nlat), outputIndex(nlat), inputIndex(nlat,maxval(inputLength)) !< variables from the main subroutine that are required
      integer, intent(in) :: dim !< the additional dimension (i.e. not nlat or nmos) for the variable in question
      real, intent(in), dimension(nlat,nmos) :: FAREROT
      real, intent(in), dimension(nlat,nmos,dim) :: val !< the input variable to be averaged and returned
      real, dimension(nlat,nmos,dim) :: weightedAvg3d !< the variable into which calculations are done
      character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
      integer :: i, j

      weightedAvg3d = val
      do i=1,nlat
        if (mode(i) /= "skip") then
            weightedAvg3d(i,outputIndex(i),:) = 0.0
            do j=1,inputLength(i)
                weightedAvg3d(i,outputIndex(i),:) =  weightedAvg3d(i,outputIndex(i),:) + (val(i,inputIndex(i,j),:) * FAREROT(i,inputIndex(i,j)))
            end do
            weightedAvg3d(i,outputIndex(i),:) = weightedAvg3d(i,outputIndex(i),:)/sum(FAREROT(i,inputIndex(i,:)))
        end if
      end do
  RETURN
  end function weightedAvg3d
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup dynamicTiling_weightedAvg4d
  !! @{
  !> This function calculates the weighted average of a four-dimensional variable (it takes the variable [val] and two additional dimensions [dim1,dim2] as inputs)
  !> @author S.R. Curasi
  function weightedAvg4d(val,dim1,dim2,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)
      
      implicit none

      integer, intent(in) :: nlat,nmos,inputLength(nlat), headTile(nlat), outputIndex(nlat), inputIndex(nlat,maxval(inputLength)) !< variables from the main subroutine that are required
      integer, intent(in) :: dim1, dim2 !< the two additional dimensions (i.e. not nlat or nmos) for the variable in question
      real, intent(in), dimension(nlat,nmos) :: FAREROT
      real, intent(in), dimension(nlat,nmos,dim1,dim2) :: val !< the input variable to be averaged and returned
      real, dimension(nlat,nmos,dim1,dim2) :: weightedAvg4d !< the variable into which calculations are done
      character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
      integer :: i, j

      weightedAvg4d = val
      do i=1,nlat
        if (mode(i) /= "skip") then
            weightedAvg4d(i,outputIndex(i),:,:) = 0.0
            do j=1,inputLength(i)
                weightedAvg4d(i,outputIndex(i),:,:) =  weightedAvg4d(i,outputIndex(i),:,:) + (val(i,inputIndex(i,j),:,:) * FAREROT(i,inputIndex(i,j)))
            end do
            weightedAvg4d(i,outputIndex(i),:,:) = weightedAvg4d(i,outputIndex(i),:,:)/sum(FAREROT(i,inputIndex(i,:)))
        end if
      end do
    weightedAvg4d = weightedAvg4d
  RETURN
  end function weightedAvg4d
  !! @}
  ! ---------------------------------------------------------------------------------------------------

!> \ingroup dynamicTiling_weightedAvg3dr8
  !! @{
  !> This function calculates the weighted average of a three-dimensional variable with 64-bit precision (it takes the variable [val] and one additional dimension [dim] as inputs)
  !> @author S.R. Curasi
  function weightedAvg3dr8(val,dim,nlat,nmos,outputIndex,headTile,inputLength,inputIndex,FAREROT,mode)

      use, intrinsic :: iso_fortran_env, only: r8=>real64

      implicit none
      
      integer, intent(in) :: nlat,nmos,inputLength(nlat), headTile(nlat), outputIndex(nlat), inputIndex(nlat,maxval(inputLength)) !< variables from the main subroutine that are required
      integer, intent(in) :: dim !< the additional dimension (i.e. not nlat or nmos) for the variable in question
      real, intent(in), dimension(nlat,nmos):: FAREROT
      real(r8), intent(in), dimension(nlat,nmos,dim) :: val !< the input variable to be averaged and returned
      character(len=5), intent(in), dimension(nlat) :: mode !< a character string that specifies the operation (i.e. move, split, merge, or skip) 
      real(r8), dimension(nlat,nmos,dim) :: weightedAvg3dr8 !< the variable into which calculations are done
      integer :: i, j

      weightedAvg3dr8 = val
      do i=1,nlat
        if (mode(i) /= "skip") then
            weightedAvg3dr8(i,outputIndex(i),:) = 0.0
            do j=1,inputLength(i)
                weightedAvg3dr8(i,outputIndex(i),:) =  weightedAvg3dr8(i,outputIndex(i),:) + (val(i,inputIndex(i,j),:) * FAREROT(i,inputIndex(i,j)))
            end do
            weightedAvg3dr8(i,outputIndex(i),:) = weightedAvg3dr8(i,outputIndex(i),:)/sum(FAREROT(i,inputIndex(i,:)))
        end if
      end do
  RETURN
  end function weightedAvg3dr8
  !! @}

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup tiledDisturbance_indexSort
  !! @{
  !> This function is used to sort tiles by age
  !> @author S.R. Curasi
  function indexSort(vector,dim)

    integer, intent(in) :: dim
    real, intent(in), dimension(dim) :: vector
    real, dimension(dim) :: sorted_vector
    logical, dimension(dim) :: mask_arr
    integer, dimension(dim) :: indexSort
    integer :: i, temp(1)

    mask_arr = .TRUE. ! Set to true so that all elements of the vector are looked at in the first time of the loop.

        do i = 1,dim
            ! Find the max value located in the vector array, excluding all positions that are false in mask_arr
            sorted_vector(i) = MAXVAL(vector,mask_arr)
            ! Then determine the index of the largest value in the vector array, accounting for mask_arr
            temp = MAXLOC(vector,mask_arr)
            ! Set the indexSort to be index found above for this round of the loop.
            indexSort(i) = temp(1)
            ! Set that largest value index to be false in the mask_arr so it is ignored the next round.
            mask_arr(indexSort(i)) = .FALSE.
        end do
        return
  end function indexSort
  !! @}

  ! ---------------------------------------------------------------------------------------------------
  !> \namespace dynamicTiling
  !! Dynamic tiling toolkit 
  !! @author S.R. Curasi
  !!
!! Overview
  !!
  !! These subroutines are designed to allow for dynamic tiling within CLASSIC. That is operations that create or remove tiles from the simulation while preserving
  !! mass balance. These subroutines are currently used to simulate harvest and fire within the classic with stand age and the resulting sub-grid scale heterogeneity
  !! explicitly represented.
  !!
  !! The functions in this toolkit can be divided into a few groups. When dynamic tiling is on cleanTiles runs regularly to ensure the model
  !! runs efficiently and there is generally enough space available for tiles to be manipulated. moveSplitCopyTiles is the major
  !! functions within the tool kit that is used in constructing the simulation. To use this tool an event needs to be broken down into a 
  !! series of operations that move, split or merge some group of tiles.
  !!
  !! For example, if a harvest occurs in a particular time step the dynamic harvest subroutine flags a group of tiles that are going to be
  !! split to create a new tile and the desired area of the resulting tile. Then once the tile had been created and we're certain mass balance is preserved
  !! the dynamic harvest subroutine disturbs the newly created tile, by removing the vegetation and allocating it to the various disturbance pools. 
  !!
  !! In the future, it should be possible to use this toolkit to accomplish more complex operations. For example in the operation above after the new tile is disturbed
  !! it could be converted to farmland by modifying the fractional area of the PFTs and then assigned to a new group in the control vector. Ideally these more complex
  !! changes would occur after the tile has been split into a separate subroutine. As long valid tile control flags are assigned and the assumptions below are followed
  !! this toolkit should interact well with more complex future operations. Leaving tile splitting and changes to the original tile separate helps ensure that mass
  !! water and energy balance is maintained and can be verified.
  !!
  !! General rules and assumptions
  !! 1. Priority should be given to the left-most tiles (i.e. closest to index 1) ideally but not 100% of the time the oldest tiles will be stored as close to index 1 as possible 
  !! this rule is also used to resolve any conflicts and unclear situations involving tiles. It is enshrined in the way the output tile location is determined in moveSplitCopyTiles
  !! and maintined by shifting all active tiles to the left after cleanTiles runs
  !!
  !! 2. To be merged, or split tiles must have the same PFT fractional area as of right now
  !!
  !! 3. Any event like disturbance that moves around mass or alters the fractional area of PFTs should occur separately from tiles being split or merged
  !! to ensure the functions in this module are universally applicable and to prevent mass balance issues
  !!
  !! 4. These subroutines should operate on the smallest set of variables possible (i.e. variables in the initialization file) and then trigger any normal setup, 
  !! to occur afterward
  !!
  !! 5. The control vector uses the following numerical flags
  !!    a. 0 = protected tile (which these subroutines will ignore)
  !!    b. 1 = group 1 modifiable tiles (there can be as many modifiable tiles as memory will allow, these tiles can created destoyed or moved)
  !!    c. -1 = the group 1 head tile (there can only be one head tile, this tile can never be destroyed, and this tile determines the abiotic properties of all subsequent group 1 tiles)
  !!    d. other numerical flags (i.e. 2, 3, 4) and their corresponding head tile flags (i.e. -2, -3, -4) could be used for future tile groups
  !!
  !! Constituent functions and subroutine
  !!
  !! 1. moveSplitCopyTiles: The main function in the dynamic tiling toolkit which allows tiles to be moved, split, or copied. This function can only be called early on in 
  !! main_core_driver because the model needs to be re-initialized. However, it can be called as many times as needed. It accepts the controlVector and a series of flags as
  !! inputs and outputs. These include:    
  !!    - inputVector: a vector of ones and zeros of length nmos which indicates which tiles or tiles should be moved, split, or copied.
  !!    - FAREAROUT: the desired fractional area of the output tile for a merge or copy operation has to be set to -9999 for a split
  !!                it should be a value between 0 and 1 which is not greater than the fractional area of the tiles marked in the inputVector
  !!    - mode: The desired operating mode (i.e. "move", "split", "merge", or "skip"). Move moves a tile to the left-most free location, split 
  !!            splits one or more tiles into a new tile placed into the left-most free location of a user-defined size, merge joins two tiles together,
  !!            and skip does nothing ("skip" was primarily implemented to allow the function to run with nlat > 1)    
  !!   -setup flag: This flag is output to signal the model that it needs to be re-initialized after the subroutine runs
  !!   -outputIndex: This is the index of the tile where moveSplitCopyTiles has placed its output. moveSplitCopyTiles automatically determines the output tile
  !!            from the pool of free tiles available this is often used by subsequent operations
  !! The subroutine carries out these operations by operating on the main state variables present in the initialization file while they are in the ROT structure. This is done
  !! for the sake of simplicity. It uses the other subroutines in this module to carry out the following steps:
  !!   -Calculate the initial mass, water, and energy content of the grid cell for later use using initialMassEnergyBalance
  !!   -Determine how the operation is to be carried out (i.e. what tile to take area from, and where to put the new tile)
  !!   -Calculate the temperatures for the new tile using MixTempsRichmann. This is modeled as if the materials in the tiles are being
  !!    blended using Richmann's law and is done first because of its complexity
  !!   -Copy the rest of the tile variables in most cases (i.e. mass-based variables) using an areally weighted average. 
  !!   -Calculate the final mass balance and check for losses of matter and energy above that expected for a given numerical precision.
  !!
  !! 2. cleanTiles: This will be a high-level subroutine that runs regularly to merge extra tiles for computational efficiency. This is the second major subroutine.
  !! It takes in the number of free tiles needed in a given year [tneeded]. For example, if two tiles are going to be split this year to accommodate disturbance 
  !! two free tiles are needed [tneeded = 2]. The subroutine will forcibly merge tiles as needed to ensure two tiles are free. It does this by calculating the difference
  !! in vegetation height between all tiles and merging the tiles with the smallest difference. If a high threshold greater than zero is set in the parameters file it will
  !! preemptively merge tiles with a height difference less than that threshold
  !!
  !! 3. initialMassEnergyBalance: This subroutine calculates the initial total mass of carbon, water, etc in a grid cell under consideration
  !! it is based heavily upon the per-grid cell output calculations done in prepareOutputs.f90. The calculations are done this way 
  !! because tile copying alters both the mass variables and fractional areas of the tiles and therefore all of this information
  !! needs to be taken into account to avoid any discrepancies. It takes as its inputs info about the model setup [i.e. DLZWROT],
  !! the fractional area of the tiles [FAREROT], the fractional cover of PFTs [fcancmxrow], and all the major mass-based variables [gleafmasrow, THLQROT, THICROT, soilcmas]
  !! it then returns the results of its calculations for each variable to its corresponding "_initial" storage variable [THLQROT_initial, THICROT_initial, gleafmasrow_initial, soilcmasrow_initial]
  !!
  !! 4. finalMassEnergyBalance: This subroutine calculates the final total mass of carbon, water, etc in a grid cell under consideration
  !! and then compares the results to their corresponding "_initial" stored versions output by initialMassEnergyBalance
  !! the inputs, outputs, and much of the calculations done in this function are the same as in initialMassEnergyBalance
  !! the calculations and comparisons use the "_final" variables which are internal to the subroutine
  !! if mass balance is not conserved it indicates an issue and this subroutine will kill the run for the particular cell
  !!
  !! 5. MixTempsRichmann: This subroutine calculated the temperatures for the new tile. It applies Richmann's law to model
  !! the different temperature materials in the tiles being mixed. See energyWaterBalanceCheck.f90 and 
  !! finalMassEnergyBalance for more background.
  !!
  !! 6. majority3d (currently unused): This function is used to determine the most common value within a 3d variable. In cases where 
  !! multiple values occur with the same frequency this function defaults to the leftmost (lowest) index
  !! this function can be used to copy control flags (i.e. lfstatus) its inputs are the variable in question [val] its dimensions [nlat,nmos] 
  !! a dimension which depends upon the variable under consideration [dim], information about which tiles to 
  !! operate on [inputLength,outputIndex,inputIndex]
  !!
  !! 7. majorityarea3d (currently unused): This function is used to determine the value within a 3d variable which represents the largest area
  !! it is very similar to the other majorityarea functions, except it weights the results by the fractional area of the tile in question
  !! this function also defaults to the leftmost (lowest) index this function can be used to copy control flags (i.e. lfstatus)
  !! its inputs are the variable in question [val] its dimensions [nlat,nmos] a dimension that depends upon the variable
  !! under consideration [dim], information about which tiles to operate on [inputLength,outputIndex,inputIndex] and tile area [FAREROT]
  !!
  !! 8. weightedAvg2d: This function is used to determine the area-weighted average of a 2d variable it uses information 
  !! about which tiles should be operated on [outputIndex,headTile,inputLength,inputIndex] the area of the tiles [FAREROT] 
  !! the variable itself [val] and its dimensions [nlat,nmos]
  !!
  !! 9. weightedAvg3d: This function is used to determine the area-weighted average of a 3d variable it uses information 
  !! about which tiles should be operated on [outputIndex,headTile,inputLength,inputIndex] the area of the tiles [FAREROT] 
  !! the variable itself [val], its basic dimensions [nlat,nmos] and one dimension which depends upon the variable under consideration [dim]
  !!
  !! 10. weightedAvg4d: This function is used to determine the area-weighted average of a 4d variable it uses information 
  !! about which tiles should be operated on [outputIndex,headTile,inputLength,inputIndex] the area of the tiles [FAREROT] 
  !! the variable itself [val], its basic dimensions [nlat,nmos], and two dimensions that depend upon the variable under
  !! consideration [dim1, dim2]
  !!
  !! 11. weightedAvg3dr8: This function is the same as weightedAvg3d except for the variable it accepts [val] can have 64bit precision
  !! currently, this function is only used with TBARROT
  !!
  !! 12. indexTileAge: This function is called to increase the age of the tile/grid cell during each CTEM timestep
  !!
  !! 13. majorityarea1d: This function is used by the tile copying subroutine within MixTempsRichmann under certain special cases to find the
  !! value within a vector that represents the largest area. It is different than the other weightedAvg functions in that it operates over
  !! latitudes one at a time. It takes [lat_in] which is the index of the latitude currently being considered as an input.
  !!
  !! 14. majorityarea2dr8 (currently unused): This function is the same as majorityarea1d except it accepts a 64bit precision input
  !! it also excepts a second dimension which depends upon the variable under consideration [dim].
  !!
  !! 15.updateHeight: This function updates the vegetation height of all tiles. This is used by cleanTiles to ensure vegetation height is updated
  !! between each comparison/merge. Because the allometric equations are not defined by parameters and applyAllometry can't be used
  !! here any updates to applyAllometry need to be copied over
  !!
  !! 16. indexSort: The function returns the indices of a vector [vector] of length [dim] sorted from highest to lowest. This is used to sort tiles by age, height, etc.
  !!
  !> \file
end module dynamicTiling
