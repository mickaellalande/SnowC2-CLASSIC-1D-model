!> \file
!> Transfers information between the 'gathered' and 'scattered' form of the CTEM data arrays.
!! @author R. Li, J. Melton, E. Chan
!!
module ctemGatherScatter

  implicit none

  public :: ctemg1
  public :: ctems2
  public :: ctemg2

contains

  !> \ingroup ctemgatherscatter_ctemg1
  !! @{
  !> Performs initial 'gather' operation on CTEM variables for consistency
  !! with physics variables gather operations.
  !!
  !! ctemg1 converts variables from the 'row' format (nlat, nmos, ...)
  !! to the 'gat' format (ilg, ...) which is what the model calculations
  !! are performed on. The ctemg1 subroutine is used to transform the
  !! read in state variables (which come in with the 'row' format from the
  !! various input files).
  !! @author R. Li, J. Melton, A. Asaadi 
  
  subroutine ctemg1 (gleafmasgat,  gleafmasgat_ns, gleafmasgat_s, & ! Out
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
                     bleafmasrow, fcancmxrow, zbtwrot, & ! In
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
                     pgleafmassrow, litrmassrow, soilcmasrow, grwtheffrow, & ! In
                     lygleafmasmaxrow, lystemmassmaxrow, lyrootmassmaxrow, tileAgerow, timharvarearow, prsfirearearow)! In

    use classicParams, only : nlat, nmos, ilg, ignd, ican, icp1, icc, iccp2

    implicit none

    real, intent(out) :: gleafmasgat(ilg,icc)
    real, intent(out) :: gleafmasgat_ns(ilg,icc)
    real, intent(out) :: gleafmasgat_s(ilg,icc)
    real, intent(out) :: bleafmasgat(ilg,icc)
    real, intent(out) :: stemmassgat(ilg,icc)
    real, intent(out) :: stemmassgat_ns(ilg,icc)
    real, intent(out) :: stemmassgat_s(ilg,icc)
    real, intent(out) :: rootmassgat(ilg,icc)
    real, intent(out) :: rootmassgat_ns(ilg,icc)
    real, intent(out) :: rootmassgat_s(ilg,icc)
    real, intent(out) :: fcancmxgat(ilg,icc)
    real, intent(out) :: zbtwgat(ilg,ignd)
    real, intent(out) :: dlzwgat(ilg,ignd)
    real, intent(out) :: sdepgat(ilg)
    real, intent(out) :: grclarea(ilg)
    real, intent(out) :: ailcggat(ilg,icc)
    real, intent(out) :: ailcbgat(ilg,icc)
    real, intent(out) :: ailcgat(ilg,ican)
    real, intent(out) :: zolncgat(ilg,ican)
    real, intent(out) :: rmatcgat(ilg,ican,ignd)
    real, intent(out) :: faregat(ilg)
    real, intent(out) :: rmatctemgat(ilg,icc,ignd)
    real, intent(out) :: slaigat(ilg,icc)
    real, intent(out) :: bmasveggat(ilg,icc)
    real, intent(out) :: cmasvegcgat(ilg,ican)
    real, intent(out) :: veghghtgat(ilg,icc)
    real, intent(out) :: rootdpthgat(ilg,icc)
    real, intent(out) :: alvsctmgat(ilg,ican)
    real, intent(out) :: alirctmgat(ilg,ican)
    real, intent(out) :: paicgat(ilg,ican)
    real, intent(out) :: slaicgat(ilg,ican)
    real, intent(out) :: maxAnnualActLyrGAT(ilg)
    integer, intent(out) :: ipeatlandgat(ilg)
    real, intent(out) :: tracermossCMassgat(ilg)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(out) :: tracermossLitrMassgat(ilg)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(out) :: tracergLeafMassgat(ilg,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerbLeafMassgat(ilg,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerstemMassgat(ilg,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerrootMassgat(ilg,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerlitrMassgat(ilg,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(out) :: tracersoilCMassgat(ilg,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(out) :: twarmmgat(ilg)
    real, intent(out) :: tcoldmgat(ilg)
    real, intent(out) :: gdd5gat(ilg)
    real, intent(out) :: ariditygat(ilg)
    real, intent(out) :: srplsmongat(ilg)
    real, intent(out) :: defctmongat(ilg)
    real, intent(out) :: anndefctgat(ilg)
    real, intent(out) :: annsrplsgat(ilg)
    real, intent(out) :: annpcpgat(ilg)
    real, intent(out) :: dry_season_lengthgat(ilg)
    real, intent(out) :: litrmsmossgat(ilg)
    real, intent(out) :: Cmossmasgat(ilg)
    real, intent(out) :: dmossgat(ilg)
    real, intent(out) :: peatSoilCgat(ilg)
    integer, intent(out) :: pandaysgat(ilg,icc)
    integer, intent(out) :: lfstatusgat(ilg,icc)
    real, intent(out) :: slopefracgat(ilg,8)
    real, intent(out) :: pstemmassgat(ilg,icc)
    real, intent(out) :: pgleafmassgat(ilg,icc)
    real, intent(out)  :: soilcmasgat(ilg,iccp2,ignd)
    real, intent(out)  :: litrmassgat(ilg,iccp2,ignd)
    real, intent(out) :: grwtheffgat(ilg,icc)
    real, intent(out) :: lygleafmasmaxgat(ilg,icc)
    real, intent(out) :: lystemmassmaxgat(ilg,icc)
    real, intent(out) :: lyrootmassmaxgat(ilg,icc)
    real, intent(out) :: tileAgegat(ilg)
    real, intent(out) :: timharvareagat(ilg)
    real, intent(out) :: prsfireareagat(ilg)

    integer, intent(in) :: ilmos (ilg)
    integer, intent(in) :: jlmos  (ilg)
    integer, intent(in) :: nml
    real, intent(in) :: gleafmasrow(nlat,nmos,icc)
    real, intent(in) :: gleafmas_NSrow(nlat,nmos,icc)
    real, intent(in) :: gleafmassrow(nlat,nmos,icc)
    real, intent(in) :: bleafmasrow(nlat,nmos,icc)
    real, intent(in) :: stemmassrow(nlat,nmos,icc)
    real, intent(in) :: stemmass_NSrow(nlat,nmos,icc)
    real, intent(in) :: stemmasssrow(nlat,nmos,icc)
    real, intent(in) :: rootmassrow(nlat,nmos,icc)
    real, intent(in) :: rootmass_NSrow(nlat,nmos,icc)
    real, intent(in) :: rootmasssrow(nlat,nmos,icc)
    real, intent(in) :: fcancmxrow(nlat,nmos,icc)
    real, intent(in) :: zbtwrot(nlat,nmos,ignd)
    real, intent(in) :: dlzwrot(nlat,nmos,ignd)
    real, intent(in) :: sdeprot(nlat,nmos)
    real, intent(in) :: ailcgrow(nlat,nmos,icc)
    real, intent(in) :: ailcbrow(nlat,nmos,icc)
    real, intent(in) :: ailcrow(nlat,nmos,ican)
    real, intent(in) :: zolncrow(nlat,nmos,ican)
    real, intent(in) :: rmatcrow(nlat,nmos,ican,ignd)
    real, intent(in) :: FAREROT(nlat,nmos)
    real, intent(in) :: rmatctemrow(nlat,nmos,icc,ignd)
    real, intent(in) :: slairow(nlat,nmos,icc)
    real, intent(in) :: bmasvegrow(nlat,nmos,icc)
    real, intent(in) :: cmasvegcrow(nlat,nmos,ican)
    real, intent(in) :: veghghtrow(nlat,nmos,icc)
    real, intent(in) :: rootdpthrow(nlat,nmos,icc)
    real, intent(in) :: alvsctmrow(nlat,nmos,ican)
    real, intent(in) :: alirctmrow(nlat,nmos,ican)
    real, intent(in) :: paicrow(nlat,nmos,ican)
    real, intent(in) :: slaicrow(nlat,nmos,ican)
    real, intent(in) :: maxAnnualActLyrROT(nlat,nmos)
    integer, intent(in) ::  ipeatlandrow(nlat,nmos)
    real, intent(in) :: tracermossCMassrot(nlat,nmos)     !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(in) :: tracermossLitrMassrot(nlat,nmos)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(in) :: tracergLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerbLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerstemMassrot(nlat,nmos,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerrootMassrot(nlat,nmos,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerlitrMassrot(nlat,nmos,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in) :: tracersoilCMassrot(nlat,nmos,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in)  :: twarmmrow(nlat)
    real, intent(in)  :: tcoldmrow(nlat)
    real, intent(in)  :: gdd5row(nlat)
    real, intent(in)  :: aridityrow(nlat)
    real, intent(in)  :: srplsmonrow(nlat)
    real, intent(in)  :: defctmonrow(nlat)
    real, intent(in)  :: anndefctrow(nlat)
    real, intent(in)  :: annsrplsrow(nlat)
    real, intent(in)  :: annpcprow(nlat)
    real, intent(in)  :: dry_season_lengthrow(nlat)
    real, intent(in)  :: litrmsmossrow(nlat,nmos)
    real, intent(in)  :: Cmossmasrow(nlat,nmos)
    real, intent(in)  :: dmossrow(nlat,nmos)
    real, intent(in)  :: peatSoilCrow(nlat,nmos)
    integer, intent(in) :: pandaysrow(nlat,nmos,icc)
    integer, intent(in) :: lfstatusrow(nlat,nmos,icc)
    real, intent(in)  :: slopefracrow(nlat,nmos,8)
    real, intent(in)  :: pstemmassrow(nlat,nmos,icc)
    real, intent(in)  :: pgleafmassrow(nlat,nmos,icc)
    real, intent(in)  :: soilcmasrow(nlat,nmos,iccp2,ignd)
    real, intent(in)  :: litrmassrow(nlat,nmos,iccp2,ignd)
    real, intent(in) :: grwtheffrow(nlat,nmos,icc)
    real, intent(in) :: lygleafmasmaxrow(nlat,nmos,icc)
    real, intent(in) :: lystemmassmaxrow(nlat,nmos,icc)
    real, intent(in) :: lyrootmassmaxrow(nlat,nmos,icc)
    
    real, intent(in) :: grclarearow(nlat)
    real, intent(in) :: tileAgerow(nlat,nmos)
    real, intent(in) :: timharvarearow(nlat,nmos)
    real, intent(in) :: prsfirearearow(nlat,nmos)

    ! Local
    integer ::  k, l, m, n

    !----------------------------------------------------------------------
    do k = 1,nml ! loop 100
      grclarea(k) = grclarearow(ilmos(k))
      sdepgat(k) = sdeprot(ilmos(k),jlmos(k))
      ipeatlandgat(k) = ipeatlandrow(ilmos(k),jlmos(k))
      faregat(k) = FAREROT(ilmos(k),jlmos(k))
      maxAnnualActLyrGAT(k) = maxAnnualActLyrROT(ilmos(k),jlmos(k))
      tracerMossCMassgat(k) = tracerMossCMassrot(ilmos(k),jlmos(k))
      tracermossLitrMassgat(k) = tracermossLitrMassrot(ilmos(k),jlmos(k))
      twarmmgat(k)    = twarmmrow(ilmos(k))
      tcoldmgat(k)    = tcoldmrow(ilmos(k))
      gdd5gat(k)      = gdd5row(ilmos(k))
      ariditygat(k)   = aridityrow(ilmos(k))
      srplsmongat(k)  = srplsmonrow(ilmos(k))
      defctmongat(k)  = defctmonrow(ilmos(k))
      anndefctgat(k)  = anndefctrow(ilmos(k))
      annsrplsgat(k)  = annsrplsrow(ilmos(k))
      annpcpgat(k)    = annpcprow(ilmos(k))
      dry_season_lengthgat(k) = dry_season_lengthrow(ilmos(k))
      litrmsmossgat(k) =  litrmsmossrow(ilmos(k),jlmos(k))
      Cmossmasgat(k) = Cmossmasrow(ilmos(k),jlmos(k))
      dmossgat(k) = dmossrow(ilmos(k),jlmos(k))
      peatSoilCgat(k) = peatSoilCrow(ilmos(k),jlmos(k))
      ipeatlandgat(k) = ipeatlandrow(ilmos(k),jlmos(k))
      tileAgegat(k)  =  tileAgerow(ilmos(k),jlmos(k))
      timharvareagat(k)  = timharvarearow(ilmos(k),jlmos(k))
      prsfireareagat(k)  = prsfirearearow(ilmos(k),jlmos(k))
      do n = 1,8
        slopefracgat(k,n) = slopefracrow(ilmos(k),jlmos(k),n)
      end do
    end do ! loop 100

    do l = 1,icc ! loop 101
      do k = 1,nml
        gleafmasgat(k,l) = gleafmasrow(ilmos(k),jlmos(k),l)
        gleafmasgat_ns(k,l) = gleafmas_NSrow(ilmos(k),jlmos(k),l)
        gleafmasgat_s(k,l) = gleafmassrow(ilmos(k),jlmos(k),l)
        bleafmasgat(k,l) = bleafmasrow(ilmos(k),jlmos(k),l)
        stemmassgat(k,l) = stemmassrow(ilmos(k),jlmos(k),l)
        stemmassgat_ns(k,l) = stemmass_NSrow(ilmos(k),jlmos(k),l)
        stemmassgat_s(k,l) = stemmasssrow(ilmos(k),jlmos(k),l)
        rootmassgat(k,l) = rootmassrow(ilmos(k),jlmos(k),l)
        rootmassgat_ns(k,l) = rootmass_NSrow(ilmos(k),jlmos(k),l)
        rootmassgat_s(k,l) = rootmasssrow(ilmos(k),jlmos(k),l)
        fcancmxgat(k,l)  = fcancmxrow(ilmos(k),jlmos(k),l)
        ailcggat(k,l)    = ailcgrow(ilmos(k),jlmos(k),l)
        ailcbgat(k,l)    = ailcbrow(ilmos(k),jlmos(k),l)
        slaigat(k,l)     = slairow(ilmos(k),jlmos(k),l)
        bmasveggat(k,l)  = bmasvegrow(ilmos(k),jlmos(k),l)
        veghghtgat(k,l)  = veghghtrow(ilmos(k),jlmos(k),l)
        rootdpthgat(k,l) = rootdpthrow(ilmos(k),jlmos(k),l)
        tracergLeafMassgat(k,l) = tracergLeafMassrot(ilmos(k),jlmos(k),l)
        tracerbLeafMassgat(k,l) = tracerbLeafMassrot(ilmos(k),jlmos(k),l)
        tracerStemMassgat(k,l) = tracerStemMassrot(ilmos(k),jlmos(k),l)
        tracerRootMassgat(k,l) = tracerRootMassrot(ilmos(k),jlmos(k),l)
        pandaysgat(k,l)  = pandaysrow(ilmos(k),jlmos(k),l)
        lfstatusgat(k,l) = lfstatusrow(ilmos(k),jlmos(k),l)
        pstemmassgat(k,l) = pstemmassrow(ilmos(k),jlmos(k),l)
        pgleafmassgat(k,l) = pgleafmassrow(ilmos(k),jlmos(k),l)
        grwtheffgat(k,l) = grwtheffrow(ilmos(k),jlmos(k),l)
        lygleafmasmaxgat(k,l) = lygleafmasmaxrow(ilmos(k),jlmos(k),l)
        lystemmassmaxgat(k,l) = lystemmassmaxrow(ilmos(k),jlmos(k),l)
        lyrootmassmaxgat(k,l) = lyrootmassmaxrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 101

    do l = 1,ican ! loop 201
      do k = 1,nml
        ailcgat(k,l)     = ailcrow(ilmos(k),jlmos(k),l)
        zolncgat(k,l)    = zolncrow(ilmos(k),jlmos(k),l)
        cmasvegcgat(k,l) = cmasvegcrow(ilmos(k),jlmos(k),l)
        alvsctmgat(k,l)  = alvsctmrow(ilmos(k),jlmos(k),l)
        alirctmgat(k,l)  = alirctmrow(ilmos(k),jlmos(k),l)
        paicgat(k,l)     = paicrow(ilmos(k),jlmos(k),l)
        slaicgat(k,l)    = slaicrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 201

    do l = 1,ignd ! loop 250
      do k = 1,nml
        zbtwgat(k,l) = zbtwrot(ilmos(k),jlmos(k),l)
        dlzwgat(k,l) = dlzwrot(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 250

    do l = 1,icc ! loop 280
      do m = 1,ignd
        do k = 1,nml
          rmatctemgat(k,l,m) = rmatctemrow(ilmos(k),jlmos(k),l,m)
        end do
      end do
    end do ! loop 280

    do l = 1,ican ! loop 290
      do m = 1,ignd
        do k = 1,nml
          rmatcgat(k,l,m) = rmatcrow(ilmos(k),jlmos(k),l,m)
        end do
      end do
    end do ! loop 290

    do l = 1,iccp2
      do k = 1,nml
        do m = 1,ignd
          tracerSoilCMassgat(k,l,m) = tracerSoilCMassrot(ilmos(k),jlmos(k),l,m)
          tracerLitrMassgat(k,l,m) = tracerLitrMassrot(ilmos(k),jlmos(k),l,m)
          litrmassgat(k,l,m) = litrmassrow(ilmos(k),jlmos(k),l,m)
          soilcmasgat(k,l,m) = soilcmasrow(ilmos(k),jlmos(k),l,m)
        end do
      end do
    end do
    return

  end subroutine ctemg1
  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup ctemgatherscatter_ctems2
  !! @{
  !> Performs subsequent scatter operation on biogeochemical variables
  !!
  !! ctems2 converts variables from the 'gat' format to the
  !! 'row' format, which is suitable for writing to output/restart
  !! files. If a variable is not written to either of those files,
  !! there is no need to scatter the variable as it will be in the
  !! correct format for model calclations ('gat').
  !!
  !! @author R. Li, J. Melton, E. Chan

  ! FLAG: Should be able to remove: xdiffus/xdiffusgat 
  
  subroutine ctems2 (fcancmxrow, rmatcrow, zolncrow, paicrow, &
                     ailcrow, ailcgrow, cmasvegcrow, slaicrow, &
                     ailcgsrow, rmatctemrow, &
                     co2concrow, co2i1cgrow, co2i1csrow, co2i2cgrow, &
                     co2i2csrow, xdiffus, slairow, cfluxcgrow, &
                     cfluxcsrow, ancsvegrow, ancgvegrow, rmlcsvegrow, &
                     rmlcgvegrow, canresrow, sdeprow, ch4concrow, &
                     sandrow, clayrow, orgmrow, &
                     anvegrow, rmlvegrow, prbfrhucgrd, &
                     extnprobgrd, pfcancmxrow, nfcancmxrow, &
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
                     neprow, nepCMIProw, hetroresrow, autoresrow, soilresprow, &
                     rmrow, rgrow, nbprow, litresrow, &
                     socresrow, gpprow, dstcemlsrow, litrfallrow, &
                     humiftrsrow, veghghtrow, rootdpthrow, rmlrow, &
                     litrfallvegrow, humiftrsvegrow, &
                     rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                     tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                     afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                     burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                     lucsocinrow, lucemcomnrow, lucltrinnrow, lucsocinnrow, nppvegrow, dstcemls3row, &
                     farerow, gavgscmsrow, &
                     rmlvegaccrow, rmsvegrow, rmrvegrow, rgvegrow, &
                     vgbiomas_vegrow, &
                     gppvegrow, vcmax0row, nepvegrow, &
                     fcanrow, pftexistrow, ccrow, mmrow, &
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
                     anmossrow, rmlmossrow, gppmossrow, armossrow, &
                     nppmossrow, peatdeprow, litrmsmossrow, Cmossmasrow, &
                     dmossrow, peatSoilCrow, pddrow, wetfrac_presrow, & 
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
                     nepgat, nepCMIPgat, hetroresgat, autoresgat, soilrespgat, &
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
                     anmossgat, rmlmossgat, gppmossgat, armossgat, &
                     nppmossgat, peatdepgat, litrmsmossgat, Cmossmasgat, &
                     dmossgat, peatSoilCgat, pddgat, wetfrac_presgat, &
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

    use classicParams, only : nlat, nmos, ilg, ignd, ican, icp1, &
                                icc, iccp2, iccp1

    implicit none

    !
    !
    real, intent(out) :: fcancmxrow(nlat,nmos,icc), rmatcrow(nlat,nmos,ican,ignd), &
                         zolncrow(nlat,nmos,ican), paicrow(nlat,nmos,ican), &
                         ailcrow(nlat,nmos,ican), ailcgrow(nlat,nmos,icc), &
                         cmasvegcrow(nlat,nmos,ican), slaicrow(nlat,nmos,ican), &
                         ailcgsrow(nlat,nmos,icc), rmatctemrow(nlat,nmos,icc,ignd),&
                         co2concrow(nlat,nmos), co2i1cgrow(nlat,nmos,icc), &
                         co2i1csrow(nlat,nmos,icc), co2i2cgrow(nlat,nmos,icc), &
                         co2i2csrow(nlat,nmos,icc), xdiffus(nlat,nmos), &
                         slairow(nlat,nmos,icc), cfluxcgrow(nlat,nmos), &
                         cfluxcsrow(nlat,nmos), ancsvegrow(nlat,nmos,icc), &
                         ancgvegrow(nlat,nmos,icc), rmlcsvegrow(nlat,nmos,icc), &
                         rmlcgvegrow(nlat,nmos,icc), canresrow(nlat,nmos), &
                         sdeprow(nlat,nmos), fcanrow(nlat,nmos,icp1), &
                         ch4concrow(nlat,nmos), wetfrac_presrow(nlat,nmos)
    !
    real, intent(out) :: sandrow(nlat,nmos,ignd), clayrow(nlat,nmos,ignd), &
                         orgmrow(nlat,nmos,ignd)
    !
    real, intent(out) :: anvegrow(nlat,nmos,icc), rmlvegrow(nlat,nmos,icc)
    !
    real, intent(out) :: prbfrhucgrd(nlat), extnprobgrd(nlat), &
                         pfcancmxrow(nlat,nmos,icc), nfcancmxrow(nlat,nmos,icc), &
                         stemmassrow(nlat,nmos,icc), stemmass_NSrow(nlat,nmos,icc), &
                         stemmasssrow(nlat,nmos,icc), rootmassrow(nlat,nmos,icc), &
                         rootmass_NSrow(nlat,nmos,icc),rootmasssrow(nlat,nmos,icc), &
                         pstemmassrow(nlat,nmos,icc), pgleafmassrow(nlat,nmos,icc), &
                         gleafmasrow(nlat,nmos,icc), gleafmas_NSrow(nlat,nmos,icc), &
                         gleafmassrow(nlat,nmos,icc), bleafmasrow(nlat,nmos,icc), &
                         ailcbrow(nlat,nmos,icc), flhrlossrow(nlat,nmos,icc), &
                         flhrloss_nsrow(nlat,nmos,icc), flhrloss_srow(nlat,nmos,icc), &
                         litrmassrow(nlat,nmos,iccp2,ignd), &
                         soilcmasrow(nlat,nmos,iccp2,ignd)

    !
    integer, intent(out) :: pandaysrow(nlat,nmos,icc), lfstatusrow(nlat,nmos,icc), &
                            colddays_leaffallrow(nlat,nmos), colddays_harvestrow(nlat,nmos)

    logical, intent(out) :: pftexistrow(nlat,nmos,icc)
    real, intent(out) :: ccrow(nlat,nmos,icc)
    real, intent(out) :: mmrow(nlat,nmos,icc)
                            !
    real, intent(out) :: grwtheffrow(nlat,nmos,icc), lystmmasrow(nlat,nmos,icc), &
                         lyrotmasrow(nlat,nmos,icc), lmaxtrow(nlat,nmos,icc), &
                         smaxtrow(nlat,nmos,icc), rmaxtrow(nlat,nmos,icc), &
                         lygleafmasmaxrow(nlat,nmos,icc), &
                         lystemmassmaxrow(nlat,nmos,icc), &
                         lyrootmassmaxrow(nlat,nmos,icc), &
                         tymaxlairow(nlat,nmos,icc), &
                         vgbiomasrow(nlat,nmos), &
                         gavgltmsrow(nlat,nmos), &
                         stmhrlosrow(nlat,nmos,icc), bmasvegrow(nlat,nmos,icc), &
                         rothrlosrow(nlat,nmos,icc), &
                         alvsctmrow(nlat,nmos,ican), alirctmrow(nlat,nmos,ican), &
                         gavglairow(nlat,nmos)
    !
    real, intent(out) :: npprow(nlat,nmos), neprow(nlat,nmos), nepCMIProw(nlat,nmos), &
                         hetroresrow(nlat,nmos), autoresrow(nlat,nmos), &
                         soilresprow(nlat,nmos), rmrow(nlat,nmos), rgrow(nlat,nmos), &
                         nbprow(nlat,nmos), litresrow(nlat,nmos), &
                         socresrow(nlat,nmos), gpprow(nlat,nmos), &
                         dstcemlsrow(nlat,nmos), litrfallrow(nlat,nmos), &
                         humiftrsrow(nlat,nmos), veghghtrow(nlat,nmos,icc), &
                         litrfallvegrow(nlat,nmos,icc), &
                         humiftrsvegrow(nlat,nmos,iccp2,ignd), &
                         rootdpthrow(nlat,nmos,icc), rmlrow(nlat,nmos), &
                         rmsrow(nlat,nmos), rmrrow(nlat,nmos), &
                         tltrleafrow(nlat,nmos,icc), tltrstemrow(nlat,nmos,icc), &
                         tltrrootrow(nlat,nmos,icc), leaflitrrow(nlat,nmos,icc), &
                         roottemprow(nlat,nmos,icc), afrleafrow(nlat,nmos,icc), &
                         afrstemrow(nlat,nmos,icc), afrrootrow(nlat,nmos,icc), &
                         wtstatusrow(nlat,nmos,icc), ltstatusrow(nlat,nmos,icc), &
                         burnfracrow(nlat,nmos), smfuncvegrow(nlat,nmos,icc), &
                         lucemcomrow(nlat,nmos), lucltrinrow(nlat,nmos), &
                         lucsocinrow(nlat,nmos), lucemcomnrow(nlat,nmos), lucltrinnrow(nlat,nmos), &
                         lucsocinnrow(nlat,nmos), nppvegrow(nlat,nmos,icc), &
                         dstcemls3row(nlat,nmos)
    !
    !     fire variables
    !
    real, intent(out) :: emit_co2row(nlat,nmos,icc), emit_corow(nlat,nmos,icc), &
                         emit_ch4row(nlat,nmos,icc), emit_nmhcrow(nlat,nmos,icc), &
                         emit_h2row(nlat,nmos,icc), emit_noxrow(nlat,nmos,icc), &
                         emit_n2orow(nlat,nmos,icc), emit_nh3row(nlat,nmos,icc), emit_pm25row(nlat,nmos,icc), &
                         emit_tpmrow(nlat,nmos,icc), emit_tcrow(nlat,nmos,icc), &
                         emit_ocrow(nlat,nmos,icc), emit_bcrow(nlat,nmos,icc), &
                         burnvegfrow(nlat,nmos,icc), btermrow(nlat,nmos,icc), &
                         ltermrow(nlat,nmos), mtermrow(nlat,nmos,icc)

    real, intent(out) :: farerow(nlat,nmos), gavgscmsrow(nlat,nmos)
    !
    real, intent(out) :: rmlvegaccrow(nlat,nmos,icc), rmsvegrow(nlat,nmos,icc), &
                         rmrvegrow(nlat,nmos,icc), rgvegrow(nlat,nmos,icc)
    !
    real, intent(out) :: vgbiomas_vegrow(nlat,nmos,icc)

    real, intent(out) :: tileAgerow(nlat,nmos), timharvarearow(nlat,nmos)

    !
    real, intent(out) :: gppvegrow(nlat,nmos,icc), nepvegrow(nlat,nmos,iccp1), &
                         nbpvegrow(nlat,nmos,iccp1), hetroresvegrow(nlat,nmos,iccp1), &
                         autoresvegrow(nlat,nmos,icc), &
                         leafns2srow(nlat,nmos,icc), stemns2srow(nlat,nmos,icc), &
                         rootns2srow(nlat,nmos,icc), vcmax0row(nlat,nmos,icc), &
                         re_alloc_s2lrow(nlat,nmos,icc), &
                         re_alloc_r2lrow(nlat,nmos,icc), &
                         re_alloc_sr2lrow(nlat,nmos,icc), &
                         litresvegrow(nlat,nmos,iccp2,ignd), &
                         soilcresvegrow(nlat,nmos,iccp2,ignd)

    real, intent(in) :: fcancmxgat(ilg,icc), rmatcgat(ilg,ican,ignd), &
                        zolncgat(ilg,ican), paicgat(ilg,ican), &
                        ailcgat(ilg,ican), ailcggat(ilg,icc), &
                        cmasvegcgat(ilg,ican), slaicgat(ilg,ican), &
                        ailcgsgat(ilg,icc), rmatctemgat(ilg,icc,ignd), &
                        co2concgat(ilg), co2i1cggat(ilg,icc), &
                        co2i1csgat(ilg,icc), co2i2cggat(ilg,icc), &
                        co2i2csgat(ilg,icc), xdiffusgat(ilg), &
                        slaigat(ilg,icc), cfluxcggat(ilg), &
                        cfluxcsgat(ilg), ancsveggat(ilg,icc), &
                        ancgveggat(ilg,icc), rmlcsveggat(ilg,icc), &
                        rmlcgveggat(ilg,icc), canresgat(ilg), &
                        sdepgat(ilg), fcangat(ilg,icp1), &
                        ch4concgat(ilg), wetfrac_presgat(ilg)
    !
    real, intent(in) :: sandgat(ilg,ignd), claygat(ilg,ignd), &
                        orgmgat(ilg,ignd)
    !
    real, intent(in) :: anveggat(ilg,icc), rmlveggat(ilg,icc)
    !
    real, intent(in) :: prbfrhucgat(ilg), extnprobgat(ilg), &
                        pfcancmxgat(ilg,icc), nfcancmxgat(ilg,icc), &
                        stemmassgat(ilg,icc), stemmassgat_ns(ilg,icc), &
                        stemmassgat_s(ilg,icc), rootmassgat(ilg,icc), &
                        rootmassgat_ns(ilg,icc), rootmassgat_s(ilg,icc), &
                        pstemmassgat(ilg,icc), pgleafmassgat(ilg,icc), &
                        gleafmasgat(ilg,icc), gleafmasgat_ns(ilg,icc), &
                        gleafmasgat_s(ilg,icc), bleafmasgat(ilg,icc), &
                        ailcbgat(ilg,icc), flhrlossgat(ilg,icc), &
                        flhrloss_nsgat(ilg,icc), flhrloss_sgat(ilg,icc), &
                        litrmassgat(ilg,iccp2,ignd), &
                        soilcmasgat(ilg,iccp2,ignd)

    integer, intent(in) :: pandaysgat(ilg,icc), lfstatusgat(ilg,icc), &
                           colddays_leaffallgat(ilg), colddays_harvestgat(ilg)

    logical, intent(in) :: pftexistgat(ilg,icc)
    real, intent(in) :: ccgat(ilg,icc)
    real, intent(in) :: mmgat(ilg,icc)
    !
    real, intent(in) :: grwtheffgat(ilg,icc), lystmmasgat(ilg,icc), &
                        lyrotmasgat(ilg,icc), lmaxtgat(ilg,icc), &
                        smaxtgat(ilg,icc), rmaxtgat(ilg,icc), &
                        lygleafmasmaxgat(ilg,icc), lystemmassmaxgat(ilg,icc), &
                        lyrootmassmaxgat(ilg,icc), tymaxlaigat(ilg,icc), &
                        vgbiomasgat(ilg), &
                        gavgltmsgat(ilg), &
                        stmhrlosgat(ilg,icc), bmasveggat(ilg,icc), &
                        rothrlosgat(ilg,icc), alvsctmgat(ilg,ican), &
                        alirctmgat(ilg,ican), gavglaigat(ilg)
    !
    real, intent(in) :: nppgat(ilg), nepgat(ilg), nepCMIPgat(ilg), &
                        hetroresgat(ilg), autoresgat(ilg), &
                        soilrespgat(ilg), rmgat(ilg), rggat(ilg), &
                        nbpgat(ilg), litresgat(ilg), &
                        socresgat(ilg), gppgat(ilg), &
                        dstcemlsgat(ilg), litrfallgat(ilg), &
                        humiftrsgat(ilg), veghghtgat(ilg,icc), &
                        litrfallveggat(ilg,icc), &
                        humiftrsveggat(ilg,iccp2,ignd), &
                        rootdpthgat(ilg,icc), rmlgat(ilg), &
                        rmsgat(ilg), rmrgat(ilg), &
                        tltrleafgat(ilg,icc), tltrstemgat(ilg,icc), &
                        tltrrootgat(ilg,icc), leaflitrgat(ilg,icc), &
                        roottempgat(ilg,icc), afrleafgat(ilg,icc), &
                        afrstemgat(ilg,icc), afrrootgat(ilg,icc), &
                        wtstatusgat(ilg,icc), ltstatusgat(ilg,icc), &
                        burnfracgat(ilg), smfuncveggat(ilg,icc), &
                        lucemcomgat(ilg), lucltringat(ilg), &
                        lucsocingat(ilg), lucemcomngat(ilg), lucltrinngat(ilg), &
                        lucsocinngat(ilg), nppveggat(ilg,icc), &
                        dstcemls3gat(ilg)
    !
    !      fire variables
    real, intent(in) :: emit_co2gat(ilg,icc), emit_cogat(ilg,icc), &
                        emit_ch4gat(ilg,icc), emit_nmhcgat(ilg,icc), &
                        emit_h2gat(ilg,icc), emit_noxgat(ilg,icc), &
                        emit_n2ogat(ilg,icc), emit_nh3gat(ilg,icc), emit_pm25gat(ilg,icc), &
                        emit_tpmgat(ilg,icc), emit_tcgat(ilg,icc), &
                        emit_ocgat(ilg,icc), emit_bcgat(ilg,icc), &
                        burnvegfgat(ilg,icc), btermgat(ilg,icc), &
                        ltermgat(ilg), mtermgat(ilg,icc)

    !
    real, intent(in) :: faregat(ilg)
    real, intent(in) :: gavgscmsgat(ilg)
    !
    real, intent(in) :: rmlvegaccgat(ilg,icc), rmsveggat(ilg,icc), &
                        rmrveggat(ilg,icc), rgveggat(ilg,icc)
    !
    real, intent(in) :: vgbiomas_veggat(ilg,icc)
    !
    real, intent(in) :: gppveggat(ilg,icc), nepveggat(ilg,iccp1), &
                        nbpveggat(ilg,iccp1), hetroresveggat(ilg,iccp1), &
                        autoresveggat(ilg,icc), &
                        leafns2sgat(ilg,icc), stemns2sgat(ilg,icc), &
                        rootns2sgat(ilg,icc), re_alloc_s2lgat(ilg,icc), &
                        re_alloc_r2lgat(ilg,icc), re_alloc_sr2lgat(ilg,icc), &
                        vcmax0gat(ilg,icc), &
                        litresveggat(ilg,iccp2,ignd), &
                        soilcresveggat(ilg,iccp2,ignd)

    !   Methane related variables
    real, intent(out)  :: ch4WetSpecrow(nlat,nmos), wetfdynrow(nlat,nmos), &
                          ch4WetDynrow(nlat,nmos), ch4soillsrow(nlat,nmos)
    real, intent(in)   :: ch4WetSpecgat(ilg), wetfdyngat(ilg), &
                          ch4WetDyngat(ilg), ch4soillsgat(ilg)

    real, intent(out) :: twarmmrow(nlat), tcoldmrow(nlat), &
                         gdd5row(nlat), aridityrow(nlat),  &
                         srplsmonrow(nlat), defctmonrow(nlat),  &
                         anndefctrow(nlat), annsrplsrow(nlat),  &
                         annpcprow(nlat), dry_season_lengthrow(nlat)

    real, intent(in) :: twarmmgat(ilg), tcoldmgat(ilg), gdd5gat(ilg), &
                        ariditygat(ilg), srplsmongat(ilg), defctmongat(ilg), &
                        anndefctgat(ilg), annsrplsgat(ilg), annpcpgat(ilg), &
                        dry_season_lengthgat(ilg)

    !   Peatland variables
    real, intent(out) :: anmossrow(nlat,nmos), rmlmossrow(nlat,nmos),  &
                         gppmossrow(nlat,nmos), armossrow(nlat,nmos),  &
                         nppmossrow(nlat,nmos), peatdeprow(nlat,nmos),  &
                         litrmsmossrow(nlat,nmos), Cmossmasrow(nlat,nmos),  &
                         dmossrow(nlat,nmos), peatSoilCrow(nlat,nmos), pddrow(nlat,nmos)
                         
    real, intent(in) :: anmossgat(ilg), rmlmossgat(ilg), gppmossgat(ilg), &
                        armossgat(ilg), nppmossgat(ilg), peatdepgat(ilg), &
                        litrmsmossgat(ilg), Cmossmasgat(ilg), dmossgat(ilg), &
                        peatSoilCgat(ilg), pddgat(ilg)

    ! allocated with nlat,nmos,...:
    real, intent(out) :: tracermossCMassrot(nlat,nmos)     !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(out) :: tracermossLitrMassrot(nlat,nmos)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(out) :: tracergLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerbLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerstemMassrot(nlat,nmos,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerrootMassrot(nlat,nmos,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerlitrMassrot(nlat,nmos,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(out) :: tracersoilCMassrot(nlat,nmos,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    real, intent(in) :: tracermossCMassgat(ilg)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(in) :: tracermossLitrMassgat(ilg)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(in) :: tracergLeafMassgat(ilg,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerbLeafMassgat(ilg,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerstemMassgat(ilg,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerrootMassgat(ilg,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerlitrMassgat(ilg,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in) :: tracersoilCMassgat(ilg,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    integer, intent(in) :: nml
    integer, intent(in)  :: ilmos(ilg), jlmos(ilg) !     gather-scatter index arrays.
    !
    integer ::  k,l,m

    ! Nitrogen cycle related variables
    real, intent(in) :: bnftotgat(ilg,iccp1), bnffreegat(ilg,iccp1), bnfantgat(ilg,icc), bnfnatgat(ilg,icc), &
                        nstressgat(ilg,icc), & 
                        nh4_massgat(ilg,iccp1), &
                        no3_massgat(ilg,iccp1), &
                        nitrifveggat(ilg,iccp1), &
                        no_nitveggat(ilg,iccp1), &
                        no_denitveggat(ilg,iccp1), &
                        no_nitdenitveggat(ilg,iccp1), &
                        n2o_nitveggat(ilg,iccp1), &
                        n2o_denitveggat(ilg,iccp1), &
                        n2o_nitdenitveggat(ilg,iccp1), &
                        n2_denitveggat(ilg,iccp1), soilpHgat(ilg), &
                        nvolveggat(ilg,iccp1), &
                        nleachveggat(ilg,iccp1), &
                        appl_fertgat(ilg,iccp1), ndep_nh4gat(ilg,iccp1), &
                        ndep_no3gat(ilg,iccp1), ngleafmasgat(ilg,icc), &
                        ngleafmasgat_ns(ilg,icc), ngleafmasgat_s(ilg,icc), &
                        nbleafmasgat(ilg,icc), redcoeff_vcmaxgat(ilg,icc), &
                        nstemmassgat(ilg,icc), nstemmassgat_ns(ilg,icc), &
                        nstemmassgat_s(ilg,icc), nrootmassgat(ilg,icc), &
                        nrootmassgat_ns(ilg,icc), nrootmassgat_s(ilg,icc), &
                        nlitrmassgat(ilg,iccp2), soilnmasgat(ilg,iccp2), &
                        nvgbiomas_veggat(ilg,icc), &
                        ndemandveg_wp_nppgat(ilg,icc), &
                        nuptakeveg_p_nh4gat(ilg,icc), &
                        nuptakeveg_p_no3gat(ilg,icc), &
                        nuptakeveg_a_actl_nh4gat(ilg,icc), &
                        nuptakeveg_a_actl_no3gat(ilg,icc), &
                        nuptakeveggat(ilg,icc), &
                        nallocveg_lgat(ilg,icc), &
                        nallocveg_sgat(ilg,icc), &
                        nallocveg_rgat(ilg,icc), &
                        nresorpedveg_sgat(ilg,icc), &
                        nresorpedveg_rgat(ilg,icc), &
                        nre_allocveg_s2lgat(ilg,icc), &
                        nre_allocveg_r2lgat(ilg,icc), &
                        nleafns2sveggat(ilg,icc), &
                        nstemns2sveggat(ilg,icc), &
                        nrootns2sveggat(ilg,icc), &
                        nlitrveg_lgat(ilg,icc), &
                        nlitrveg_sgat(ilg,icc), &
                        nlitrveg_rgat(ilg,icc), &
                        nlitrveggat(ilg,icc), &
                        gl2bl_grass_nfluxgat(ilg,icc), &
                        c2nveg_lgat(ilg,icc), c2nveg_sgat(ilg,icc), &
                        c2nveg_rgat(ilg,icc), c2nveg_wpgat(ilg,icc), &
                        c2nveg_litrgat(ilg,iccp1), c2nveg_humusgat(ilg,iccp1), &
                        nhumtrsveggat(ilg,iccp1), &
                        nmineralveg_litrgat(ilg,iccp1), &
                        nmineralveg_humusgat(ilg,iccp1), &
                        netnmineralveg_gat(ilg,iccp1), &
                        nimmobilveg_nh4gat(ilg,iccp1), &
                        nimmobilveg_no3gat(ilg,iccp1), &
                        fNnetlandveggat(ilg,iccp1)

    real, intent(out) :: bnftotrow(nlat,nmos,iccp1), &
                         bnffreerow(nlat,nmos,iccp1), bnfantrow(nlat,nmos,icc), bnfnatrow(nlat,nmos,icc), &
                         nstressrow(nlat,nmos,icc), &
                         nh4_massrow(nlat,nmos,iccp1), &
                         no3_massrow(nlat,nmos,iccp1), &
                         nitrifvegrow(nlat,nmos,iccp1), soilpHrow(nlat,nmos), &
                         no_nitvegrow(nlat,nmos,iccp1), &
                         no_denitvegrow(nlat,nmos,iccp1), &
                         no_nitdenitvegrow(nlat,nmos,iccp1), &
                         nvolvegrow(nlat,nmos,iccp1), &
                         n2o_nitvegrow(nlat,nmos,iccp1), &
                         n2o_denitvegrow(nlat,nmos,iccp1), &
                         n2o_nitdenitvegrow(nlat,nmos,iccp1), &
                         n2_denitvegrow(nlat,nmos,iccp1), &
                         nleachvegrow(nlat,nmos,iccp1), &
                         appl_fertrow(nlat,nmos,iccp1), ndep_nh4row(nlat,nmos,iccp1), &
                         ndep_no3row(nlat,nmos,iccp1), ngleafmasrow(nlat,nmos,icc), &
                         ngleafmas_NSrow(nlat,nmos,icc), &
                         ngleafmassrow(nlat,nmos,icc), nbleafmasrow(nlat,nmos,icc), &
                         nstemmassrow(nlat,nmos,icc), nstemmass_NSrow(nlat,nmos,icc), &
                         nstemmasssrow(nlat,nmos,icc), nrootmassrow(nlat,nmos,icc), &
                         nrootmass_NSrow(nlat,nmos,icc), nrootmasssrow(nlat,nmos,icc), &
                         nlitrmassrow(nlat,nmos,iccp2), soilnmasrow(nlat,nmos,iccp2), &
                         nvgbiomas_vegrow(nlat,nmos,icc), &
                         ndemandveg_wp_npprow(nlat,nmos,icc), &
                         nuptakeveg_p_nh4row(nlat,nmos,icc), &
                         nuptakeveg_p_no3row(nlat,nmos,icc), &
                         nuptakeveg_a_actl_nh4row(nlat,nmos,icc), &
                         nuptakeveg_a_actl_no3row(nlat,nmos,icc), &
                         nuptakevegrow(nlat,nmos,icc), &
                         nallocveg_lrow(nlat,nmos,icc), &
                         nallocveg_srow(nlat,nmos,icc), &
                         nallocveg_rrow(nlat,nmos,icc), &
                         nresorpedveg_srow(nlat,nmos,icc), &
                         nresorpedveg_rrow(nlat,nmos,icc), &
                         nre_allocveg_s2lrow(nlat,nmos,icc), &
                         nre_allocveg_r2lrow(nlat,nmos,icc), &
                         nleafns2svegrow(nlat,nmos,icc), &
                         nstemns2svegrow(nlat,nmos,icc), &
                         nrootns2svegrow(nlat,nmos,icc), &
                         nlitrveg_lrow(nlat,nmos,icc), &
                         nlitrveg_srow(nlat,nmos,icc), &
                         nlitrveg_rrow(nlat,nmos,icc), &
                         nlitrvegrow(nlat,nmos,icc), &
                         gl2bl_grass_nfluxrow(nlat,nmos,icc), &
                         c2nveg_lrow(nlat,nmos,icc), c2nveg_srow(nlat,nmos,icc), &
                         c2nveg_rrow(nlat,nmos,icc), c2nveg_wprow(nlat,nmos,icc), &
                         c2nveg_litrrow(nlat,nmos,iccp1), c2nveg_humusrow(nlat,nmos,iccp1), &
                         nhumtrsvegrow(nlat,nmos,iccp1), &
                         nmineralveg_litrrow(nlat,nmos,iccp1), &
                         nmineralveg_humusrow(nlat,nmos,iccp1), &
                         netnmineralveg_row(nlat,nmos,iccp1), &
                         nimmobilveg_nh4row(nlat,nmos,iccp1), &
                         nimmobilveg_no3row(nlat,nmos,iccp1), &
                         fNnetlandvegrow(nlat,nmos,iccp1), &
                         redcoeff_vcmaxrow(nlat,nmos,icc)

    real, intent(in) :: tileAgegat(ilg), timharvareagat(ilg)
    
    !----------------------------------------------------------------------
    do k = 1,nml ! loop 100
      sdeprow(ilmos(k),jlmos(k))        = sdepgat(k)
      co2concrow(ilmos(k),jlmos(k))     = co2concgat(k)
      ch4concrow(ilmos(k),jlmos(k))     = ch4concgat(k)
      cfluxcgrow(ilmos(k),jlmos(k))     = cfluxcggat(k)
      cfluxcsrow(ilmos(k),jlmos(k))     = cfluxcsgat(k)
      canresrow(ilmos(k),jlmos(k))      = canresgat(k)
      xdiffus(ilmos(k),jlmos(k))        = xdiffusgat(k)
      prbfrhucgrd(ilmos(k))             = prbfrhucgat(k)
      extnprobgrd(ilmos(k))             = extnprobgat(k)
      vgbiomasrow(ilmos(k),jlmos(k))    = vgbiomasgat(k)
      gavgltmsrow(ilmos(k),jlmos(k))    = gavgltmsgat(k)
      gavglairow(ilmos(k),jlmos(k))     = gavglaigat(k)
      npprow(ilmos(k),jlmos(k))         = nppgat(k)
      neprow(ilmos(k),jlmos(k))         = nepgat(k)
      nepCMIProw(ilmos(k),jlmos(k))     = nepCMIPgat(k)
      hetroresrow(ilmos(k),jlmos(k))    = hetroresgat(k)
      autoresrow(ilmos(k),jlmos(k))     = autoresgat(k)
      soilresprow(ilmos(k),jlmos(k))    = soilrespgat(k)
      rmrow(ilmos(k),jlmos(k))          = rmgat(k)
      rgrow(ilmos(k),jlmos(k))          = rggat(k)
      nbprow(ilmos(k),jlmos(k))         = nbpgat(k)
      litresrow(ilmos(k),jlmos(k))      = litresgat(k)
      socresrow(ilmos(k),jlmos(k))      = socresgat(k)
      gpprow(ilmos(k),jlmos(k))         = gppgat(k)
      dstcemlsrow(ilmos(k),jlmos(k))    = dstcemlsgat(k)
      litrfallrow(ilmos(k),jlmos(k))    = litrfallgat(k)
      humiftrsrow(ilmos(k),jlmos(k))    = humiftrsgat(k)
      rmlrow(ilmos(k),jlmos(k))         = rmlgat(k)
      rmsrow(ilmos(k),jlmos(k))         = rmsgat(k)
      rmrrow(ilmos(k),jlmos(k))         = rmrgat(k)
      burnfracrow(ilmos(k),jlmos(k))    = burnfracgat(k)
      ltermrow(ilmos(k),jlmos(k))       = ltermgat(k)
      lucemcomrow(ilmos(k),jlmos(k))    = lucemcomgat(k)
      lucltrinrow(ilmos(k),jlmos(k))    = lucltringat(k)
      lucsocinrow(ilmos(k),jlmos(k))    = lucsocingat(k)
      lucemcomnrow(ilmos(k),jlmos(k))    = lucemcomngat(k)
      lucltrinnrow(ilmos(k),jlmos(k))    = lucltrinngat(k)
      lucsocinnrow(ilmos(k),jlmos(k))    = lucsocinngat(k)
      dstcemls3row(ilmos(k),jlmos(k))   = dstcemls3gat(k)
      farerow(ilmos(k),jlmos(k))        = faregat(k)
      gavgscmsrow(ilmos(k),jlmos(k))    = gavgscmsgat(k)
      ch4WetSpecrow(ilmos(k),jlmos(k))     = ch4WetSpecgat(k)
      wetfdynrow(ilmos(k),jlmos(k))     = wetfdyngat(k)
      ch4WetDynrow(ilmos(k),jlmos(k))     = ch4WetDyngat(k)
      ch4soillsrow(ilmos(k),jlmos(k))   = ch4soillsgat(k)
      wetfrac_presrow(ilmos(k),jlmos(k)) = wetfrac_presgat(k)

      twarmmrow(ilmos(k)) = twarmmgat(k)
      tcoldmrow(ilmos(k)) = tcoldmgat(k)
      gdd5row(ilmos(k)) = gdd5gat(k)
      aridityrow(ilmos(k)) = ariditygat(k)
      srplsmonrow(ilmos(k)) = srplsmongat(k)
      defctmonrow(ilmos(k)) = defctmongat(k)
      anndefctrow(ilmos(k)) = anndefctgat(k)
      annsrplsrow(ilmos(k)) = annsrplsgat(k)
      annpcprow(ilmos(k)) = annpcpgat(k)
      dry_season_lengthrow(ilmos(k)) = dry_season_lengthgat(k)

      tracerMossCMassrot(ilmos(k),jlmos(k)) = tracerMossCMassgat(k)
      tracermossLitrMassrot(ilmos(k),jlmos(k)) = tracermossLitrMassgat(k)

      soilpHrow(ilmos(k),jlmos(k))      = soilpHgat(k)
      colddays_leaffallrow(ilmos(k),jlmos(k)) = colddays_leaffallgat(k)
      colddays_harvestrow(ilmos(k),jlmos(k)) = colddays_harvestgat(k)
      tileAgerow(ilmos(k),jlmos(k)) = tileAgegat(k)
      timharvarearow(ilmos(k),jlmos(k)) = timharvareagat(k)
    end do ! loop 100
    !
    do l = 1,icc ! loop 101
      do k = 1,nml
        smfuncvegrow(ilmos(k),jlmos(k),l) = smfuncveggat(k,l)
        mtermrow(ilmos(k),jlmos(k),l)     = mtermgat(k,l)
        btermrow(ilmos(k),jlmos(k),l)     = btermgat(k,l)
        ailcgrow(ilmos(k),jlmos(k),l)     = ailcggat(k,l)
        ailcgsrow(ilmos(k),jlmos(k),l)    = ailcgsgat(k,l)
        co2i1cgrow(ilmos(k),jlmos(k),l)   = co2i1cggat(k,l)
        co2i1csrow(ilmos(k),jlmos(k),l)   = co2i1csgat(k,l)
        co2i2cgrow(ilmos(k),jlmos(k),l)   = co2i2cggat(k,l)
        co2i2csrow(ilmos(k),jlmos(k),l)   = co2i2csgat(k,l)
        slairow(ilmos(k),jlmos(k),l)      = slaigat(k,l)
        anvegrow(ilmos(k),jlmos(k),l)     = anveggat(k,l)
        rmlvegrow(ilmos(k),jlmos(k),l)    = rmlveggat(k,l)
        pfcancmxrow(ilmos(k),jlmos(k),l)  = pfcancmxgat(k,l)
        nfcancmxrow(ilmos(k),jlmos(k),l)  = nfcancmxgat(k,l)
        fcancmxrow(ilmos(k),jlmos(k),l)   = fcancmxgat(k,l)
        stemmassrow(ilmos(k),jlmos(k),l)  = stemmassgat(k,l)
        stemmass_NSrow(ilmos(k),jlmos(k),l)= stemmassgat_ns(k,l)
        stemmasssrow(ilmos(k),jlmos(k),l) = stemmassgat_s(k,l)
        rootmassrow(ilmos(k),jlmos(k),l)  = rootmassgat(k,l)
        rootmass_NSrow(ilmos(k),jlmos(k),l)= rootmassgat_ns(k,l)
        rootmasssrow(ilmos(k),jlmos(k),l) = rootmassgat_s(k,l)
        pstemmassrow(ilmos(k),jlmos(k),l)  = pstemmassgat(k,l)
        pgleafmassrow(ilmos(k),jlmos(k),l)  = pgleafmassgat(k,l)
        gleafmasrow(ilmos(k),jlmos(k),l)  = gleafmasgat(k,l)
        gleafmas_NSrow(ilmos(k),jlmos(k),l)= gleafmasgat_ns(k,l)
        gleafmassrow(ilmos(k),jlmos(k),l) = gleafmasgat_s(k,l)
        bleafmasrow(ilmos(k),jlmos(k),l)  = bleafmasgat(k,l)
        ailcbrow(ilmos(k),jlmos(k),l)     = ailcbgat(k,l)
        flhrlossrow(ilmos(k),jlmos(k),l)  = flhrlossgat(k,l)
        flhrloss_nsrow(ilmos(k),jlmos(k),l)  = flhrloss_nsgat(k,l)
        flhrloss_srow(ilmos(k),jlmos(k),l)  = flhrloss_sgat(k,l)
        pandaysrow(ilmos(k),jlmos(k),l)   = pandaysgat(k,l)
        lfstatusrow(ilmos(k),jlmos(k),l)  = lfstatusgat(k,l)
        grwtheffrow(ilmos(k),jlmos(k),l)  = grwtheffgat(k,l)
        lystmmasrow(ilmos(k),jlmos(k),l)  = lystmmasgat(k,l)
        lyrotmasrow(ilmos(k),jlmos(k),l)  = lyrotmasgat(k,l)
        lmaxtrow(ilmos(k),jlmos(k),l)     = lmaxtgat(k,l)
        smaxtrow(ilmos(k),jlmos(k),l)     = smaxtgat(k,l)
        rmaxtrow(ilmos(k),jlmos(k),l)     = rmaxtgat(k,l)
        lygleafmasmaxrow(ilmos(k),jlmos(k),l) = lygleafmasmaxgat(k,l)
        lystemmassmaxrow(ilmos(k),jlmos(k),l) = lystemmassmaxgat(k,l)
        lyrootmassmaxrow(ilmos(k),jlmos(k),l) = lyrootmassmaxgat(k,l)
        tymaxlairow(ilmos(k),jlmos(k),l)  = tymaxlaigat(k,l)
        stmhrlosrow(ilmos(k),jlmos(k),l)  = stmhrlosgat(k,l)
        bmasvegrow(ilmos(k),jlmos(k),l)   = bmasveggat(k,l)
        rothrlosrow(ilmos(k),jlmos(k),l)  = rothrlosgat(k,l)
        veghghtrow(ilmos(k),jlmos(k),l)   = veghghtgat(k,l)
        rootdpthrow(ilmos(k),jlmos(k),l)  = rootdpthgat(k,l)
        tltrleafrow(ilmos(k),jlmos(k),l)  = tltrleafgat(k,l)
        tltrstemrow(ilmos(k),jlmos(k),l)  = tltrstemgat(k,l)
        tltrrootrow(ilmos(k),jlmos(k),l)  = tltrrootgat(k,l)
        leaflitrrow(ilmos(k),jlmos(k),l)  = leaflitrgat(k,l)
        roottemprow(ilmos(k),jlmos(k),l)  = roottempgat(k,l)
        afrleafrow(ilmos(k),jlmos(k),l)   = afrleafgat(k,l)
        afrstemrow(ilmos(k),jlmos(k),l)   = afrstemgat(k,l)
        afrrootrow(ilmos(k),jlmos(k),l)   = afrrootgat(k,l)
        wtstatusrow(ilmos(k),jlmos(k),l)  = wtstatusgat(k,l)
        ltstatusrow(ilmos(k),jlmos(k),l)  = ltstatusgat(k,l)
        nppvegrow(ilmos(k),jlmos(k),l)    = nppveggat(k,l)
        rmlvegaccrow(ilmos(k),jlmos(k),l) = rmlvegaccgat(k,l)
        rmsvegrow(ilmos(k),jlmos(k),l)    = rmsveggat(k,l)
        rmrvegrow(ilmos(k),jlmos(k),l)    = rmrveggat(k,l)
        rgvegrow(ilmos(k),jlmos(k),l)     = rgveggat(k,l)
        gppvegrow(ilmos(k),jlmos(k),l)    = gppveggat(k,l)
        vcmax0row(ilmos(k),jlmos(k),l)   = vcmax0gat(k,l)
        leafns2srow(ilmos(k),jlmos(k),l) = leafns2sgat(k,l)
        stemns2srow(ilmos(k),jlmos(k),l) = stemns2sgat(k,l)
        rootns2srow(ilmos(k),jlmos(k),l) = rootns2sgat(k,l)
        re_alloc_s2lrow(ilmos(k),jlmos(k),l) = re_alloc_s2lgat(k,l)
        re_alloc_r2lrow(ilmos(k),jlmos(k),l) = re_alloc_r2lgat(k,l)
        re_alloc_sr2lrow(ilmos(k),jlmos(k),l) = re_alloc_sr2lgat(k,l)
        vgbiomas_vegrow(ilmos(k),jlmos(k),l) = vgbiomas_veggat(k,l)
        autoresvegrow(ilmos(k),jlmos(k),l) = autoresveggat(k,l)
        pftexistrow(ilmos(k),jlmos(k),l) = pftexistgat(k,l)
        ccrow(ilmos(k),jlmos(k),l) = ccgat(k,l)
        mmrow(ilmos(k),jlmos(k),l) = mmgat(k,l)
        ancsvegrow(ilmos(k),jlmos(k),l) = ancsveggat(k,l)
        ancgvegrow(ilmos(k),jlmos(k),l) = ancgveggat(k,l)
        rmlcsvegrow(ilmos(k),jlmos(k),l) = rmlcsveggat(k,l)
        rmlcgvegrow(ilmos(k),jlmos(k),l) = rmlcgveggat(k,l)
        litrfallvegrow(ilmos(k),jlmos(k),l) = litrfallveggat(k,l)

        ngleafmasrow(ilmos(k),jlmos(k),l)   = ngleafmasgat(k,l)
        ngleafmas_NSrow(ilmos(k),jlmos(k),l) = ngleafmasgat_ns(k,l)
        ngleafmassrow(ilmos(k),jlmos(k),l)  = ngleafmasgat_s(k,l)
        nbleafmasrow(ilmos(k),jlmos(k),l)   = nbleafmasgat(k,l)
        nstemmassrow(ilmos(k),jlmos(k),l)   = nstemmassgat(k,l)
        nstemmass_NSrow(ilmos(k),jlmos(k),l) = nstemmassgat_ns(k,l)
        nstemmasssrow(ilmos(k),jlmos(k),l)  = nstemmassgat_s(k,l)
        nrootmassrow(ilmos(k),jlmos(k),l)   = nrootmassgat(k,l)
        nrootmass_NSrow(ilmos(k),jlmos(k),l) = nrootmassgat_ns(k,l)
        nrootmasssrow(ilmos(k),jlmos(k),l)  = nrootmassgat_s(k,l)
        nvgbiomas_vegrow(ilmos(k),jlmos(k),l) = nvgbiomas_veggat(k,l)
        ndemandveg_wp_npprow(ilmos(k),jlmos(k),l) = ndemandveg_wp_nppgat(k,l)
        c2nveg_lrow(ilmos(k),jlmos(k),l)    = c2nveg_lgat(k,l)
        c2nveg_srow(ilmos(k),jlmos(k),l)    = c2nveg_sgat(k,l)
        c2nveg_rrow(ilmos(k),jlmos(k),l)    = c2nveg_rgat(k,l)
        c2nveg_wprow(ilmos(k),jlmos(k),l)   = c2nveg_wpgat(k,l)
        nuptakeveg_p_nh4row(ilmos(k),jlmos(k),l) = nuptakeveg_p_nh4gat(k,l)
        nuptakeveg_p_no3row(ilmos(k),jlmos(k),l) = nuptakeveg_p_no3gat(k,l)
        nuptakeveg_a_actl_nh4row(ilmos(k),jlmos(k),l) = nuptakeveg_a_actl_nh4gat(k,l)
        nuptakeveg_a_actl_no3row(ilmos(k),jlmos(k),l) = nuptakeveg_a_actl_no3gat(k,l)
        nuptakevegrow(ilmos(k),jlmos(k),l) = nuptakeveggat(k,l)
        nallocveg_lrow(ilmos(k),jlmos(k),l) = nallocveg_lgat(k,l)
        nallocveg_srow(ilmos(k),jlmos(k),l) = nallocveg_sgat(k,l)
        nallocveg_rrow(ilmos(k),jlmos(k),l) = nallocveg_rgat(k,l)
        nresorpedveg_srow(ilmos(k),jlmos(k),l) = nresorpedveg_sgat(k,l)
        nresorpedveg_rrow(ilmos(k),jlmos(k),l) = nresorpedveg_rgat(k,l)
        nre_allocveg_s2lrow(ilmos(k),jlmos(k),l) = nre_allocveg_s2lgat(k,l)
        nre_allocveg_r2lrow(ilmos(k),jlmos(k),l) = nre_allocveg_r2lgat(k,l)
        nleafns2svegrow(ilmos(k),jlmos(k),l) = nleafns2sveggat(k,l)
        nstemns2svegrow(ilmos(k),jlmos(k),l) = nstemns2sveggat(k,l)
        nrootns2svegrow(ilmos(k),jlmos(k),l) = nrootns2sveggat(k,l)
        nlitrveg_lrow(ilmos(k),jlmos(k),l)   = nlitrveg_lgat(k,l)
        nlitrveg_srow(ilmos(k),jlmos(k),l)   = nlitrveg_sgat(k,l)
        nlitrveg_rrow(ilmos(k),jlmos(k),l)   = nlitrveg_rgat(k,l)
        nlitrvegrow(ilmos(k),jlmos(k),l)     = nlitrveggat(k,l)
        gl2bl_grass_nfluxrow(ilmos(k),jlmos(k),l) = gl2bl_grass_nfluxgat(k,l)
        redcoeff_vcmaxrow(ilmos(k),jlmos(k),l) = redcoeff_vcmaxgat(k,l)
        !         fire variables
        emit_co2row(ilmos(k),jlmos(k),l)    = emit_co2gat(k,l)
        emit_corow(ilmos(k),jlmos(k),l)     = emit_cogat(k,l)
        emit_ch4row(ilmos(k),jlmos(k),l)    = emit_ch4gat(k,l)
        emit_nmhcrow(ilmos(k),jlmos(k),l)   = emit_nmhcgat(k,l)
        emit_h2row(ilmos(k),jlmos(k),l)     = emit_h2gat(k,l)
        emit_noxrow(ilmos(k),jlmos(k),l)    = emit_noxgat(k,l)
        emit_n2orow(ilmos(k),jlmos(k),l)    = emit_n2ogat(k,l)
        emit_nh3row(ilmos(k),jlmos(k),l)    = emit_nh3gat(k,l)
        emit_pm25row(ilmos(k),jlmos(k),l)   = emit_pm25gat(k,l)
        emit_tpmrow(ilmos(k),jlmos(k),l)    = emit_tpmgat(k,l)
        emit_tcrow(ilmos(k),jlmos(k),l)     = emit_tcgat(k,l)
        emit_ocrow(ilmos(k),jlmos(k),l)     = emit_ocgat(k,l)
        emit_bcrow(ilmos(k),jlmos(k),l)     = emit_bcgat(k,l)
        burnvegfrow(ilmos(k),jlmos(k),l)    = burnvegfgat(k,l)

        tracergLeafMassrot(ilmos(k),jlmos(k),l) = tracergLeafMassgat(k,l)
        tracerbLeafMassrot(ilmos(k),jlmos(k),l) = tracerbLeafMassgat(k,l)
        tracerStemMassrot(ilmos(k),jlmos(k),l) = tracerStemMassgat(k,l)
        tracerRootMassrot(ilmos(k),jlmos(k),l) = tracerRootMassgat(k,l)
        bnfantrow(ilmos(k),jlmos(k),l)    = bnfantgat(k,l)
        bnfnatrow(ilmos(k),jlmos(k),l)    = bnfnatgat(k,l)
        nstressrow(ilmos(k),jlmos(k),l)    = nstressgat(k,l)
      end do
    end do ! loop 101
    !
    do l = 1,iccp1 ! loop 102
      do k = 1,nml
        hetroresvegrow(ilmos(k),jlmos(k),l) = hetroresveggat(k,l)
        nepvegrow(ilmos(k),jlmos(k),l)    = nepveggat(k,l)
        nbpvegrow(ilmos(k),jlmos(k),l)    = nbpveggat(k,l)
        nh4_massrow(ilmos(k),jlmos(k),l)  = nh4_massgat(k,l)
        no3_massrow(ilmos(k),jlmos(k),l)  = no3_massgat(k,l)
        bnftotrow(ilmos(k),jlmos(k),l)    = bnftotgat(k,l)
        bnffreerow(ilmos(k),jlmos(k),l)   = bnffreegat(k,l)
        nitrifvegrow(ilmos(k),jlmos(k),l) = nitrifveggat(k,l)
        no_nitvegrow(ilmos(k),jlmos(k),l) = no_nitveggat(k,l)
        no_denitvegrow(ilmos(k),jlmos(k),l) = no_denitveggat(k,l)
        no_nitdenitvegrow(ilmos(k),jlmos(k),l) = no_nitdenitveggat(k,l)
        n2o_nitvegrow(ilmos(k),jlmos(k),l)  = n2o_nitveggat(k,l)
        n2o_denitvegrow(ilmos(k),jlmos(k),l) = n2o_denitveggat(k,l)
        n2o_nitdenitvegrow(ilmos(k),jlmos(k),l) = n2o_nitdenitveggat(k,l)
        n2_denitvegrow(ilmos(k),jlmos(k),l) = n2_denitveggat(k,l)
        nvolvegrow(ilmos(k),jlmos(k),l)   = nvolveggat(k,l)
        nleachvegrow(ilmos(k),jlmos(k),l) = nleachveggat(k,l)
        appl_fertrow(ilmos(k),jlmos(k),l) = appl_fertgat(k,l)
        ndep_nh4row(ilmos(k),jlmos(k),l)  = ndep_nh4gat(k,l)
        ndep_no3row(ilmos(k),jlmos(k),l)  = ndep_no3gat(k,l)
        c2nveg_litrrow(ilmos(k),jlmos(k),l)  = c2nveg_litrgat(k,l)
        c2nveg_humusrow(ilmos(k),jlmos(k),l) = c2nveg_humusgat(k,l)
        nhumtrsvegrow(ilmos(k),jlmos(k),l) = nhumtrsveggat(k,l)
        nmineralveg_litrrow(ilmos(k),jlmos(k),l) = nmineralveg_litrgat(k,l)
        nmineralveg_humusrow(ilmos(k),jlmos(k),l) = nmineralveg_humusgat(k,l)
        netnmineralveg_row(ilmos(k),jlmos(k),l) = netnmineralveg_gat(k,l)
        nimmobilveg_nh4row(ilmos(k),jlmos(k),l) = nimmobilveg_nh4gat(k,l)
        nimmobilveg_no3row(ilmos(k),jlmos(k),l) = nimmobilveg_no3gat(k,l)
        fNnetlandvegrow(ilmos(k),jlmos(k),l) = fNnetlandveggat(k,l)
      end do
    end do ! loop 102

    do l = 1,iccp2 ! loop 103
      do k = 1,nml
        do m = 1,ignd
          litrmassrow(ilmos(k),jlmos(k),l,m) = litrmassgat(k,l,m)
          soilcmasrow(ilmos(k),jlmos(k),l,m) = soilcmasgat(k,l,m)
          litresvegrow(ilmos(k),jlmos(k),l,m) = litresveggat(k,l,m)
          soilcresvegrow(ilmos(k),jlmos(k),l,m)=soilcresveggat(k,l,m)
          humiftrsvegrow(ilmos(k),jlmos(k),l,m) = humiftrsveggat(k,l,m)
          tracerSoilCMassrot(ilmos(k),jlmos(k),l,m) = tracerSoilCMassgat(k,l,m)
          tracerLitrMassrot(ilmos(k),jlmos(k),l,m) = tracerLitrMassgat(k,l,m)
        end do
        nlitrmassrow(ilmos(k),jlmos(k),l) = nlitrmassgat(k,l)
        soilnmasrow(ilmos(k),jlmos(k),l)  = soilnmasgat(k,l)
      end do
    end do ! loop 103
    !
    do l = 1,ican ! loop 201
      do k = 1,nml
        ailcrow(ilmos(k),jlmos(k),l)    = ailcgat(k,l)
        zolncrow(ilmos(k),jlmos(k),l)   = zolncgat(k,l)
        cmasvegcrow(ilmos(k),jlmos(k),l) = cmasvegcgat(k,l)
        paicrow(ilmos(k),jlmos(k),l)    = paicgat(k,l)
        slaicrow(ilmos(k),jlmos(k),l)   = slaicgat(k,l)
        alvsctmrow(ilmos(k),jlmos(k),l) = alvsctmgat(k,l)
        alirctmrow(ilmos(k),jlmos(k),l) = alirctmgat(k,l)
      end do
    end do ! loop 201
    !
    do l = 1,ignd ! loop 250
      do k = 1,nml
        sandrow(ilmos(k),jlmos(k),l)     = sandgat(k,l)
        clayrow(ilmos(k),jlmos(k),l)     = claygat(k,l)
        orgmrow(ilmos(k),jlmos(k),l)     = orgmgat(k,l)
      end do
    end do ! loop 250
    !
    do l = 1,icc ! loop 280
      do m = 1,ignd
        do k = 1,nml
          rmatctemrow(ilmos(k),jlmos(k),l,m) = rmatctemgat(k,l,m)
        end do
      end do
    end do ! loop 280
    !
    do l = 1,ican ! loop 290
      do m = 1,ignd
        do k = 1,nml
          rmatcrow(ilmos(k),jlmos(k),l,m) = rmatcgat(k,l,m)
        end do
      end do
    end do ! loop 290

    !     this class variable is scattered here, but it gathered in classg,
    !     not in ctemg2. jm jan 82013.
    do l = 1,icp1 ! loop 300
      do k = 1,nml
        fcanrow(ilmos(k),jlmos(k),l)     = fcangat(k,l)
      end do
    end do ! loop 300

    !    scatter peatland variables----------------------------------------\
    do k = 1,nml ! loop 400
      anmossrow(ilmos(k),jlmos(k)) = anmossgat(k)
      rmlmossrow(ilmos(k),jlmos(k)) = rmlmossgat(k)
      gppmossrow(ilmos(k),jlmos(k)) = gppmossgat(k)
      armossrow(ilmos(k),jlmos(k)) = armossgat(k)
      nppmossrow(ilmos(k),jlmos(k)) = nppmossgat(k)
      peatdeprow (ilmos(k),jlmos(k))  = peatdepgat(k)
      litrmsmossrow(ilmos(k),jlmos(k)) = litrmsmossgat(k)
      Cmossmasrow(ilmos(k),jlmos(k)) = Cmossmasgat(k)
      dmossrow(ilmos(k),jlmos(k)) = dmossgat(k)
      peatSoilCrow(ilmos(k),jlmos(k)) = peatSoilCgat(k)
      pddrow(ilmos(k),jlmos(k)) = pddgat(k)
    end do ! loop 400

    return

  end subroutine ctems2
  !! @}
  ! ------------------------------------------------------------------------------------
  !> \ingroup ctemgatherscatter_ctemg2
  !! @{
  !> Performs subsequent 'gather' operation on CTEM variables for consistency
  !! with physics variables gather operations.
  !!
  !! ctemg2 takes variables in the 'row' format (nlat, nmos, ...)
  !! and converts them to the 'gat' format (ilg, ...). At present
  !! ctemg2 is bloated with many variables that do not require
  !! gathering. This subroutine should ideally be only used for
  !! state variables that are updated from external files as
  !! the run progresses. Since the model calculations operate
  !! on the 'gat' form, any other variables need not be gathered
  !! as they will already be in the correct format from the previous
  !! model timestep.
  !!
  !! @author R. Li, Y. Wu, E. Chan, J. Melton
  !!
  subroutine ctemg2 (fcancmxgat, ailcgsgat, fcancsgat, fcancgat, &
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
                     gavgltmsgat, stmhrlosgat, colddays_leaffallgat, colddays_harvestgat,  &
                     rothrlosgat,gavglaigat, nppgat, &
                     nepgat, nepCMIPgat, hetroresgat, autoresgat, soilrespgat, &
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
                     rmlcgvegrow, canresrow, sdeprow, ch4concrow, &
                     sandrow, clayrow, orgmrow, &
                     anvegrow, rmlvegrow, prbfrhucrow, &
                     extnprobrow, pfcancmxrow, nfcancmxrow, &
                     stemmassrow, stemmassnsrow, stemmasssrow, rootmassrow, &
                     rootmassnsrow, rootmasssrow, litrmassrow, gleafmasrow, &
                     gleafmasnsrow, gleafmassrow, bleafmasrow, soilcmasrow, flhrlossrow, flhrloss_nsrow, flhrloss_srow, &
                     pandaysrow, lfstatusrow, grwtheffrow, lystmmasrow, &
                     lyrotmasrow, lmaxtrow, smaxtrow, rmaxtrow, &
                     lygleafmasmaxrow, lystemmassmaxrow, lyrootmassmaxrow, &
                     tymaxlairow, vgbiomasrow, &
                     gavgltmsrow, stmhrlosrow, colddays_leaffallrow, colddays_harvestrow, &
                     rothrlosrow, gavglairow, npprow, &
                     neprow, nepCMIProw, hetroresrow, autoresrow, soilresprow, &
                     rmrow, rgrow, nbprow, litresrow, &
                     socresrow, gpprow, dstcemlsrow, litrfallrow, &
                     litrfallvegrow, humiftrsrow, rmlrow, &
                     rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                     tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                     afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                     burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                     lucsocinrow, lucemcomnrow, lucltrinnrow, lucsocinnrow, dstcemls3row, popdinrow, &
                     farerow, gavgscmsrow, rmlvegaccrow, pftexistrow, ccrow, mmrow, &
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
                     ngleafmasrow, ngleafmasnsrow, ngleafmassrow, &
                     nbleafmasrow, nstemmassrow, nstemmassnsrow, nstemmasssrow, &
                     nrootmassrow, nrootmassnsrow, nrootmasssrow, nlitrmassrow, &
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
    !     May 06 2020    Bring in N cycle variables.
    !     A. Asaadi
    !
    !     Dec 232016     Remove thlqaccXXX_m/thicaccXXX_m
    !     Ed Chan
    !     March 192015   Gathering of peatland variables
    !     Yuanqiao Wu
    !
    !     July 122013    Bring in the ctem params use statement
    !     J. Melton
    !     July 282009    Gather operation on CTEM variables.
    !     Rong Li
    !
    use classicParams,      only : nlat, nmos, ilg, ignd, ican, icp1, &
                                   icc, iccp2, iccp1
    implicit none
    !
    !     * integer :: constants.
    !
    integer :: k, l, m, n, j
    integer, intent(in) :: nml
    !
    !     * gather-scatter index arrays.
    !
    integer, intent(in) :: ilmos(ilg), jlmos(ilg)
    !
    !
    real, intent(out) :: fcancmxgat(ilg,icc), &
                         ailcgsgat(ilg,icc), fcancsgat(ilg,icc), &
                         fcancgat(ilg,icc), &
                         co2concgat(ilg), co2i1cggat(ilg,icc), &
                         co2i1csgat(ilg,icc), co2i2cggat(ilg,icc), &
                         co2i2csgat(ilg,icc), xdiffusgat(ilg), &
                         cfluxcggat(ilg), &
                         cfluxcsgat(ilg), ancsveggat(ilg,icc), &
                         ancgveggat(ilg,icc), rmlcsveggat(ilg,icc), &
                         rmlcgveggat(ilg,icc), canresgat(ilg), &
                         sdepgat(ilg), ch4concgat(ilg)
    !
    real, intent(out) :: sandgat(ilg,ignd), claygat(ilg,ignd), &
                         orgmgat(ilg,ignd)
    !
    real, intent(out) :: anveggat(ilg,icc), rmlveggat(ilg,icc)
    !
    real, intent(out) :: prbfrhucgat(ilg), extnprobgat(ilg), &
                         daylgat(ilg), dayl_maxgat(ilg), &
                         pfcancmxgat(ilg,icc), nfcancmxgat(ilg,icc), &
                         stemmassgat(ilg,icc), stemmassgat_ns(ilg,icc), &
                         stemmassgat_s(ilg,icc), rootmassgat(ilg,icc), &
                         rootmassgat_ns(ilg,icc),rootmassgat_s(ilg,icc), &
                         pstemmassgat(ilg,icc), pgleafmassgat(ilg,icc), &
                         gleafmasgat(ilg,icc), gleafmasgat_ns(ilg,icc), &
                         gleafmasgat_s(ilg,icc), bleafmasgat(ilg,icc), &
                         flhrlossgat(ilg,icc), flhrloss_nsgat(ilg,icc), flhrloss_sgat(ilg,icc), &
     			 soilcmasgat(ilg,iccp2,ignd), litrmassgat(ilg,iccp2,ignd)

    integer, intent(out) :: pandaysgat(ilg,icc), lfstatusgat(ilg,icc), &
                            colddays_leaffallgat(ilg), colddays_harvestgat(ilg)
    logical, intent(out) :: pftexistgat(ilg,icc)
    real, intent(out) :: ccgat(ilg,icc)
    real, intent(out) :: mmgat(ilg,icc)
    !
    real, intent(out) :: grwtheffgat(ilg,icc), lystmmasgat(ilg,icc), &
                         lyrotmasgat(ilg,icc), lmaxtgat(ilg,icc), &
                         smaxtgat(ilg,icc), rmaxtgat(ilg,icc), &
                         lygleafmasmaxgat(ilg,icc), lystemmassmaxgat(ilg,icc), &
                         lyrootmassmaxgat(ilg,icc), &
                         tymaxlaigat(ilg,icc), &
                         vgbiomasgat(ilg), &
                         gavgltmsgat(ilg), &
                         stmhrlosgat(ilg,icc), &
                         rothrlosgat(ilg,icc), &
                         gavglaigat(ilg)
    !
    real, intent(out) :: nppgat(ilg), nepgat(ilg), nepCMIPgat(ilg), &
                         hetroresgat(ilg), autoresgat(ilg), &
                         soilrespgat(ilg), rmgat(ilg), rggat(ilg), &
                         nbpgat(ilg), litresgat(ilg), &
                         socresgat(ilg), gppgat(ilg), &
                         dstcemlsgat(ilg), litrfallgat(ilg), &
                         litrfallveggat(ilg,icc), &
                         humiftrsgat(ilg), &
                         rmlgat(ilg), &
                         rmsgat(ilg), rmrgat(ilg), &
                         tltrleafgat(ilg,icc), tltrstemgat(ilg,icc), &
                         tltrrootgat(ilg,icc), leaflitrgat(ilg,icc), &
                         roottempgat(ilg,icc), afrleafgat(ilg,icc), &
                         afrstemgat(ilg,icc), afrrootgat(ilg,icc), &
                         wtstatusgat(ilg,icc), ltstatusgat(ilg,icc), &
                         burnfracgat(ilg), smfuncveggat(ilg,icc), &
                         lucemcomgat(ilg), lucltringat(ilg), &
                         lucsocingat(ilg), lucemcomngat(ilg), lucltrinngat(ilg), lucsocinngat(ilg),&
                         nppveggat(ilg,icc), &
                         dstcemls3gat(ilg), &
                         tracerCO2gat(ilg)
    !
    !      fire emission variables
    real, intent(out) :: emit_co2gat(ilg,icc), emit_cogat(ilg,icc), &
                         emit_ch4gat(ilg,icc), emit_nmhcgat(ilg,icc), &
                         emit_h2gat(ilg,icc), emit_noxgat(ilg,icc), &
                         emit_n2ogat(ilg,icc), emit_nh3gat(ilg,icc), emit_pm25gat(ilg,icc), &
                         emit_tpmgat(ilg,icc), emit_tcgat(ilg,icc), &
                         emit_ocgat(ilg,icc), emit_bcgat(ilg,icc), &
                         burnvegfgat(ilg,icc), btermgat(ilg,icc), &
                         ltermgat(ilg), mtermgat(ilg,icc), &
                         popdingat(ilg)
    !
    real, intent(out) :: faregat(ilg)
    real, intent(out) :: gavgscmsgat(ilg)
    !
    real, intent(out) :: rmlvegaccgat(ilg,icc), rmsveggat(ilg,icc), &
                         rmrveggat(ilg,icc), rgveggat(ilg,icc)
    !
    real, intent(out) :: vgbiomas_veggat(ilg,icc)
    !
    real, intent(out) :: gppveggat(ilg,icc), nepveggat(ilg,iccp1), &
                         nbpveggat(ilg,iccp1), hetroresveggat(ilg,iccp1), &
                         autoresveggat(ilg,icc), &
                         leafns2sgat(ilg,icc), stemns2sgat(ilg,icc), &
                         rootns2sgat(ilg,icc), re_alloc_s2lgat(ilg,icc), &
                         re_alloc_r2lgat(ilg,icc), re_alloc_sr2lgat(ilg,icc), &
                         vcmax0gat(ilg,icc), &
                         litresveggat(ilg,iccp2,ignd), &
                         soilcresveggat(ilg,iccp2,ignd)
    
    real, intent(in) :: fcancmxrow(nlat,nmos,icc), &
                        ailcgsrow(nlat,nmos,icc), fcancsrow(nlat,nmos,icc), &
                        fcancrow(nlat,nmos,icc), &
                        co2concrow(nlat,nmos), co2i1cgrow(nlat,nmos,icc), &
                        co2i1csrow(nlat,nmos,icc), co2i2cgrow(nlat,nmos,icc), &
                        co2i2csrow(nlat,nmos,icc), xdiffus(nlat,nmos), &
                        cfluxcgrow(nlat,nmos), &
                        cfluxcsrow(nlat,nmos), ancsvegrow(nlat,nmos,icc), &
                        ancgvegrow(nlat,nmos,icc), rmlcsvegrow(nlat,nmos,icc), &
                        rmlcgvegrow(nlat,nmos,icc), canresrow(nlat,nmos), &
                        sdeprow(nlat,nmos), ch4concrow(nlat,nmos)
    !
    real, intent(in) :: sandrow(nlat,nmos,ignd), clayrow(nlat,nmos,ignd), &
                        orgmrow(nlat,nmos,ignd)
    !
    real, intent(in) :: anvegrow(nlat,nmos,icc), rmlvegrow(nlat,nmos,icc)
    !
    real, intent(in) :: prbfrhucrow(nlat,nmos), extnprobrow(nlat,nmos), &
                        daylrow(nlat), dayl_maxrow(nlat), &
                        pfcancmxrow(nlat,nmos,icc), nfcancmxrow(nlat,nmos,icc), &
                        stemmassrow(nlat,nmos,icc), stemmassnsrow(nlat,nmos,icc), &
                        stemmasssrow(nlat,nmos,icc), rootmassrow(nlat,nmos,icc), &
                        rootmassnsrow(nlat,nmos,icc), rootmasssrow(nlat,nmos,icc), &
                        pstemmassrow(nlat,nmos,icc), pgleafmassrow(nlat,nmos,icc), &
                        gleafmasrow(nlat,nmos,icc), gleafmasnsrow(nlat,nmos,icc), &
                        gleafmassrow(nlat,nmos,icc), bleafmasrow(nlat,nmos,icc), &
                        flhrlossrow(nlat,nmos,icc), flhrloss_nsrow(nlat,nmos,icc), flhrloss_srow(nlat,nmos,icc), &
                        soilcmasrow(nlat,nmos,iccp2,ignd), &
                        litrmassrow(nlat,nmos,iccp2,ignd)
    !
    integer, intent(in) :: pandaysrow(nlat,nmos,icc), lfstatusrow(nlat,nmos,icc), &
                           colddays_leaffallrow(nlat,nmos), colddays_harvestrow(nlat,nmos)
    logical, intent(in) :: pftexistrow(nlat,nmos,icc)
    real, intent(in) :: ccrow(nlat,nmos,icc)
    real, intent(in) :: mmrow(nlat,nmos,icc)
    !
    real, intent(in) :: grwtheffrow(nlat,nmos,icc), lystmmasrow(nlat,nmos,icc), &
                        lyrotmasrow(nlat,nmos,icc), lmaxtrow(nlat,nmos,icc), &
                        smaxtrow(nlat,nmos,icc), rmaxtrow(nlat,nmos,icc), &
                        lygleafmasmaxrow(nlat,nmos,icc), &
                        lystemmassmaxrow(nlat,nmos,icc), &
                        lyrootmassmaxrow(nlat,nmos,icc), &
                        tymaxlairow(nlat,nmos,icc), &
                        vgbiomasrow(nlat,nmos), &
                        gavgltmsrow(nlat,nmos), &
                        stmhrlosrow(nlat,nmos,icc), &
                        rothrlosrow(nlat,nmos,icc), &
                        gavglairow(nlat,nmos)
    !
    real, intent(in) :: npprow(nlat,nmos), neprow(nlat,nmos), nepCMIProw(nlat,nmos), &
                        hetroresrow(nlat,nmos), autoresrow(nlat,nmos), &
                        soilresprow(nlat,nmos), rmrow(nlat,nmos), &
                        rgrow(nlat,nmos), &
                        nbprow(nlat,nmos), litresrow(nlat,nmos), &
                        socresrow(nlat,nmos), gpprow(nlat,nmos), &
                        dstcemlsrow(nlat,nmos), litrfallrow(nlat,nmos), &
                        litrfallvegrow(nlat,nmos,icc), &
                        humiftrsrow(nlat,nmos), &
                        rmlrow(nlat,nmos), &
                        rmsrow(nlat,nmos), rmrrow(nlat,nmos), &
                        tltrleafrow(nlat,nmos,icc), tltrstemrow(nlat,nmos,icc), &
                        tltrrootrow(nlat,nmos,icc), leaflitrrow(nlat,nmos,icc), &
                        roottemprow(nlat,nmos,icc), afrleafrow(nlat,nmos,icc), &
                        afrstemrow(nlat,nmos,icc), afrrootrow(nlat,nmos,icc), &
                        wtstatusrow(nlat,nmos,icc), ltstatusrow(nlat,nmos,icc), &
                        burnfracrow(nlat,nmos), smfuncvegrow(nlat,nmos,icc), &
                        lucemcomrow(nlat,nmos), lucltrinrow(nlat,nmos), &
                        lucsocinrow(nlat,nmos), lucemcomnrow(nlat,nmos), lucltrinnrow(nlat,nmos), &
                        lucsocinnrow(nlat,nmos), nppvegrow(nlat,nmos,icc), &
                        dstcemls3row(nlat,nmos), tracerCO2rot(nlat,nmos)
    !
    !     fire variables
    real, intent(in) :: emit_co2row(nlat,nmos,icc), emit_corow(nlat,nmos,icc), &
                        emit_ch4row(nlat,nmos,icc), emit_nmhcrow(nlat,nmos,icc), &
                        emit_h2row(nlat,nmos,icc), emit_noxrow(nlat,nmos,icc), &
                        emit_n2orow(nlat,nmos,icc), emit_nh3row(nlat,nmos,icc), emit_pm25row(nlat,nmos,icc), &
                        emit_tpmrow(nlat,nmos,icc), emit_tcrow(nlat,nmos,icc), &
                        emit_ocrow(nlat,nmos,icc), emit_bcrow(nlat,nmos,icc), &
                        burnvegfrow(nlat,nmos,icc), btermrow(nlat,nmos,icc), &
                        ltermrow(nlat,nmos), mtermrow(nlat,nmos,icc), &
                        popdinrow(nlat,nmos)
    !
    real, intent(in) :: farerow(nlat,nmos)
    !
    real, intent(in) :: gavgscmsrow(nlat,nmos)
    !
    real, intent(in) :: rmlvegaccrow(nlat,nmos,icc), rmsvegrow(nlat,nmos,icc), &
                        rmrvegrow(nlat,nmos,icc), rgvegrow(nlat,nmos,icc)
    !
    real, intent(in) :: vgbiomas_vegrow(nlat,nmos,icc)
    !
    real, intent(in) :: gppvegrow(nlat,nmos,icc), nepvegrow(nlat,nmos,iccp1), &
                        nbpvegrow(nlat,nmos,iccp1), hetroresvegrow(nlat,nmos,iccp1), &
                        autoresvegrow(nlat,nmos,icc), &
                        litresvegrow(nlat,nmos,iccp2,ignd), &
                        soilcresvegrow(nlat,nmos,iccp2,ignd), &
                        leafns2srow(nlat,nmos,icc), stemns2srow(nlat,nmos,icc), &
                        rootns2srow(nlat,nmos,icc), re_alloc_s2lrow(nlat,nmos,icc), &
                        re_alloc_r2lrow(nlat,nmos,icc), re_alloc_sr2lrow(nlat,nmos,icc), &
                        vcmax0row(nlat,nmos,icc)
    !   Methane related variables
    real, intent(in)  :: slopefracrow(nlat,nmos,8), ch4WetSpecrow(nlat,nmos), wetfdynrow(nlat,nmos), &
                         ch4WetDynrow(nlat,nmos), ch4soillsrow(nlat,nmos)
    real, intent(out) :: slopefracgat(ilg,8), ch4WetSpecgat(ilg), wetfdyngat(ilg), &
                         ch4WetDyngat(ilg), ch4soillsgat(ilg)
    !    --peatland variables
    real, intent(in)  :: anmossrow(nlat,nmos), rmlmossrow(nlat,nmos), gppmossrow(nlat,nmos), &
                         armossrow(nlat,nmos), nppmossrow(nlat,nmos), peatdeprow(nlat,nmos), &
                         litrmsmossrow(nlat,nmos), Cmossmasrow(nlat,nmos), dmossrow(nlat,nmos), &
                         peatSoilCrow(nlat,nmos), pddrow(nlat,nmos)

    real, intent(out) :: anmossgat(ilg), rmlmossgat(ilg), gppmossgat(ilg), &
                         armossgat(ilg), nppmossgat(ilg), peatdepgat(ilg), &
                         litrmsmossgat(ilg), Cmossmasgat(ilg), dmossgat(ilg), &
                         peatSoilCgat(ilg), pddgat(ilg)

    ! Nitrogen-cycle related variables
    real, intent(out) :: bnftotgat(ilg,iccp1), bnffreegat(ilg,iccp1), bnfantgat(ilg,icc), bnfnatgat(ilg,icc), &
                         nstressgat(ilg,icc), &
                         soilpHgat(ilg), nfertilgat(ilg), &
                         ndepositgat(ilg), nh4_massgat(ilg,iccp1), &
                         no3_massgat(ilg,iccp1), &
                         nitrifveggat(ilg,iccp1), no_nitveggat(ilg,iccp1), &
                         no_denitveggat(ilg,iccp1), &
                         no_nitdenitveggat(ilg,iccp1),  &
                         n2o_nitveggat(ilg,iccp1), &
                         n2o_denitveggat(ilg,iccp1), &
                         n2o_nitdenitveggat(ilg,iccp1), &
                         n2_denitveggat(ilg,iccp1), nvolveggat(ilg,iccp1), &
                         nleachveggat(ilg,iccp1), appl_fertgat(ilg,iccp1), &
                         ndep_nh4gat(ilg,iccp1), ndep_no3gat(ilg,iccp1), &
                         ngleafmasgat(ilg,icc), ngleafmasgat_ns(ilg,icc), &
                         ngleafmasgat_s(ilg,icc), nbleafmasgat(ilg,icc), nstemmassgat(ilg,icc), &
                         nstemmassgat_ns(ilg,icc),nstemmassgat_s(ilg,icc), &
                         nrootmassgat(ilg,icc),   nrootmassgat_ns(ilg,icc), &
                         nrootmassgat_s(ilg,icc), nlitrmassgat(ilg,iccp2), soilnmasgat(ilg,iccp2), &
                         nvgbiomas_veggat(ilg,icc), &
                         ndemandveg_wp_nppgat(ilg,icc), &
                         nuptakeveg_p_nh4gat(ilg,icc), &
                         nuptakeveg_p_no3gat(ilg,icc), &
                         nuptakeveg_a_actl_nh4gat(ilg,icc), &
                         nuptakeveg_a_actl_no3gat(ilg,icc), &
                         nuptakeveggat(ilg,icc), &
                         nallocveg_lgat(ilg,icc), nallocveg_sgat(ilg,icc), &
                         nallocveg_rgat(ilg,icc), &
                         nresorpedveg_sgat(ilg,icc), &
                         nresorpedveg_rgat(ilg,icc), &
                         nre_allocveg_s2lgat(ilg,icc), &
                         nre_allocveg_r2lgat(ilg,icc), &
                         nleafns2sveggat(ilg,icc), nstemns2sveggat(ilg,icc), &
                         nrootns2sveggat(ilg,icc), &
                         nlitrveg_lgat(ilg,icc), nlitrveg_sgat(ilg,icc), &
                         nlitrveg_rgat(ilg,icc), nlitrveggat(ilg,icc), &
                         gl2bl_grass_nfluxgat(ilg,icc), c2nveg_lgat(ilg,icc), &
                         c2nveg_sgat(ilg,icc), c2nveg_rgat(ilg,icc), c2nveg_wpgat(ilg,icc), &
                         c2nveg_litrgat(ilg,iccp1), c2nveg_humusgat(ilg,iccp1), &
                         nhumtrsveggat(ilg,iccp1), &
                         nmineralveg_litrgat(ilg,iccp1), &
                         nmineralveg_humusgat(ilg,iccp1), &
                         netnmineralveg_gat(ilg,iccp1), &
                         nimmobilveg_nh4gat(ilg,iccp1), &
                         nimmobilveg_no3gat(ilg,iccp1), &
                         fNnetlandveggat(ilg,iccp1), &
                         redcoeff_vcmaxgat(ilg,icc)
     real, intent(in) :: bnftotrow(nlat,nmos,iccp1), &
                         bnffreerow(nlat,nmos,iccp1), bnfantrow(nlat,nmos,icc), bnfnatrow(nlat,nmos,icc), &
                         nstressrow(nlat,nmos,icc), &
                         soilpHrow(nlat,nmos), nfertilrow(nlat,nmos), ndepositrow(nlat,nmos), &
                         nh4_massrow(nlat,nmos,iccp1), no3_massrow(nlat,nmos,iccp1), &
                         nitrifvegrow(nlat,nmos,iccp1), &
                         no_nitvegrow(nlat,nmos,iccp1), &
                         no_denitvegrow(nlat,nmos,iccp1), &
                         no_nitdenitvegrow(nlat,nmos,iccp1), &
                         n2o_nitvegrow(nlat,nmos,iccp1), &
                         n2o_denitvegrow(nlat,nmos,iccp1), &
                         n2o_nitdenitvegrow(nlat,nmos,iccp1), & 
                         n2_denitvegrow(nlat,nmos,iccp1), &
                         nvolvegrow(nlat,nmos,iccp1), & 
                         nleachvegrow(nlat,nmos,iccp1), appl_fertrow(nlat,nmos,iccp1), &
                         ndep_nh4row(nlat,nmos,iccp1), ndep_no3row(nlat,nmos,iccp1), &
                         ngleafmasrow(nlat,nmos,icc), ngleafmasnsrow(nlat,nmos,icc), &
                         ngleafmassrow(nlat,nmos,icc), nbleafmasrow(nlat,nmos,icc), &
                         nstemmassrow(nlat,nmos,icc), nstemmassnsrow(nlat,nmos,icc), &
                         nstemmasssrow(nlat,nmos,icc), nrootmassrow(nlat,nmos,icc), &
                         nrootmassnsrow(nlat,nmos,icc), nrootmasssrow(nlat,nmos,icc), &
                         nlitrmassrow(nlat,nmos,iccp2), soilnmasrow(nlat,nmos,iccp2), &
                         nvgbiomas_vegrow(nlat,nmos,icc), & 
                         ndemandveg_wp_npprow(nlat,nmos,icc), &
                         nuptakeveg_p_nh4row(nlat,nmos,icc), &
                         nuptakeveg_p_no3row(nlat,nmos,icc), &
                         nuptakeveg_a_actl_nh4row(nlat,nmos,icc), &
                         nuptakeveg_a_actl_no3row(nlat,nmos,icc), &
                         nuptakevegrow(nlat,nmos,icc), &
                         nallocveg_lrow(nlat,nmos,icc), &
                         nallocveg_srow(nlat,nmos,icc), &
                         nallocveg_rrow(nlat,nmos,icc), &
                         nresorpedveg_srow(nlat,nmos,icc), &
                         nresorpedveg_rrow(nlat,nmos,icc), &
                         nre_allocveg_s2lrow(nlat,nmos,icc), &
                         nre_allocveg_r2lrow(nlat,nmos,icc), &
                         nleafns2svegrow(nlat,nmos,icc), &
                         nstemns2svegrow(nlat,nmos,icc), &
                         nrootns2svegrow(nlat,nmos,icc), &
                         nlitrveg_lrow(nlat,nmos,icc), &
                         nlitrveg_srow(nlat,nmos,icc), &
                         nlitrveg_rrow(nlat,nmos,icc), &
                         nlitrvegrow(nlat,nmos,icc), &
                         gl2bl_grass_nfluxrow(nlat,nmos,icc), c2nveg_lrow(nlat,nmos,icc), &
                         c2nveg_srow(nlat,nmos,icc), c2nveg_rrow(nlat,nmos,icc), &
                         c2nveg_wprow(nlat,nmos,icc), c2nveg_litrrow(nlat,nmos,iccp1), &
                         c2nveg_humusrow(nlat,nmos,iccp1), &
                         nhumtrsvegrow(nlat,nmos,iccp1), &
                         nmineralveg_litrrow(nlat,nmos,iccp1), &
                         nmineralveg_humusrow(nlat,nmos,iccp1), &
                         netnmineralveg_row(nlat,nmos,iccp1), &
                         nimmobilveg_nh4row(nlat,nmos,iccp1), &
                         nimmobilveg_no3row(nlat,nmos,iccp1), & 
                         fNnetlandvegrow(nlat,nmos,iccp1), redcoeff_vcmaxrow(nlat,nmos,icc)
    
    real, intent(in) :: tileAgerow(nlat,nmos)
    real, intent(in) :: timharvarearow(nlat,nmos)
    real, intent(out) :: tileAgegat(ilg)
    real, intent(out) :: timharvareagat(ilg)
    real, intent(in) :: prsfirearearow(nlat,nmos)
    real, intent(out) :: prsfireareagat(ilg)

    !----------------------------------------------------------------------
    do k = 1,nml ! loop 100
      sdepgat(k)      = sdeprow(ilmos(k),jlmos(k))
      co2concgat(k)   = co2concrow(ilmos(k),jlmos(k))
      ch4concgat(k)   = ch4concrow(ilmos(k),jlmos(k))
      cfluxcggat(k)   = cfluxcgrow(ilmos(k),jlmos(k))
      cfluxcsgat(k)   = cfluxcsrow(ilmos(k),jlmos(k))
      canresgat(k)    = canresrow(ilmos(k),jlmos(k))
      xdiffusgat(k)   = xdiffus(ilmos(k),jlmos(k))
      prbfrhucgat(k)  = prbfrhucrow(ilmos(k),jlmos(k))
      extnprobgat(k)  = extnprobrow(ilmos(k),jlmos(k))
      daylgat(k)      = daylrow(ilmos(k))
      dayl_maxgat(k)  = dayl_maxrow(ilmos(k))
      vgbiomasgat(k)  = vgbiomasrow(ilmos(k),jlmos(k))
      gavgltmsgat(k)  = gavgltmsrow(ilmos(k),jlmos(k))
      gavglaigat(k)   = gavglairow(ilmos(k),jlmos(k))
      nppgat(k)       = npprow(ilmos(k),jlmos(k))
      nepgat(k)       = neprow(ilmos(k),jlmos(k))
      nepCMIPgat(k)   = nepCMIProw(ilmos(k),jlmos(k))
      hetroresgat(k)  = hetroresrow(ilmos(k),jlmos(k))
      autoresgat(k)   = autoresrow(ilmos(k),jlmos(k))
      soilrespgat(k)  = soilresprow(ilmos(k),jlmos(k))
      rmgat(k)        = rmrow(ilmos(k),jlmos(k))
      rggat(k)        = rgrow(ilmos(k),jlmos(k))
      nbpgat(k)       = nbprow(ilmos(k),jlmos(k))
      litresgat(k)    = litresrow(ilmos(k),jlmos(k))
      socresgat(k)    = socresrow(ilmos(k),jlmos(k))
      gppgat(k)       = gpprow(ilmos(k),jlmos(k))
      dstcemlsgat(k)  = dstcemlsrow(ilmos(k),jlmos(k))
      litrfallgat(k)  = litrfallrow(ilmos(k),jlmos(k))
      humiftrsgat(k)  = humiftrsrow(ilmos(k),jlmos(k))
      rmlgat(k)       = rmlrow(ilmos(k),jlmos(k))
      rmsgat(k)       = rmsrow(ilmos(k),jlmos(k))
      rmrgat(k)       = rmrrow(ilmos(k),jlmos(k))
      burnfracgat(k)  = burnfracrow(ilmos(k),jlmos(k))
      lucemcomgat(k)  = lucemcomrow(ilmos(k),jlmos(k))
      lucltringat(k)  = lucltrinrow(ilmos(k),jlmos(k))
      lucsocingat(k)  = lucsocinrow(ilmos(k),jlmos(k))
      lucemcomngat(k)  = lucemcomnrow(ilmos(k),jlmos(k))
      lucltrinngat(k)  = lucltrinnrow(ilmos(k),jlmos(k))
      lucsocinngat(k)  = lucsocinnrow(ilmos(k),jlmos(k))
      dstcemls3gat(k) = dstcemls3row(ilmos(k),jlmos(k))
      tracerCO2gat(k) = tracerCO2rot(ilmos(k),jlmos(k))
      faregat(k)      = farerow(ilmos(k),jlmos(k))
      gavgscmsgat(k)  = gavgscmsrow(ilmos(k),jlmos(k))
      soilpHgat(k)    = soilpHrow(ilmos(k),jlmos(k))
      nfertilgat(k)   = nfertilrow(ilmos(k),jlmos(k))
      ndepositgat(k)  = ndepositrow(ilmos(k),jlmos(k))
      do n = 1,8
        slopefracgat(k,n) = slopefracrow(ilmos(k),jlmos(k),n)
      end do
      ch4WetSpecgat(k)   = ch4WetSpecrow(ilmos(k),jlmos(k))
      wetfdyngat(k)   = wetfdynrow(ilmos(k),jlmos(k))
      ch4WetDyngat(k)   = ch4WetDynrow(ilmos(k),jlmos(k))
      ch4soillsgat(k) = ch4soillsrow(ilmos(k),jlmos(k))
      popdingat(k)    = popdinrow(ilmos(k),jlmos(k))
      colddays_leaffallgat(k) = colddays_leaffallrow(ilmos(k),jlmos(k))
      colddays_harvestgat(k) = colddays_harvestrow(ilmos(k),jlmos(k))
      tileAgegat(k)  =  tileAgerow(ilmos(k),jlmos(k))
      timharvareagat(k)  = timharvarearow(ilmos(k),jlmos(k))
      prsfireareagat(k)  = prsfirearearow(ilmos(k),jlmos(k))
    end do ! loop 100

    do l = 1,icc ! loop 101
      do k = 1,nml
        btermgat(k,l)    = btermrow(ilmos(k),jlmos(k),l)
        ltermgat(k)      = ltermrow(ilmos(k),jlmos(k))  ! not per pft but keeping in with the other term vars.
        mtermgat(k,l)    = mtermrow(ilmos(k),jlmos(k),l)
        smfuncveggat(k,l) = smfuncvegrow(ilmos(k),jlmos(k),l)
        fcancmxgat(k,l)  = fcancmxrow(ilmos(k),jlmos(k),l)
        ailcgsgat(k,l)   = ailcgsrow(ilmos(k),jlmos(k),l)
        fcancsgat(k,l)   = fcancsrow(ilmos(k),jlmos(k),l)
        fcancgat(k,l)    = fcancrow(ilmos(k),jlmos(k),l)
        co2i1cggat(k,l)  = co2i1cgrow(ilmos(k),jlmos(k),l)
        co2i1csgat(k,l)  = co2i1csrow(ilmos(k),jlmos(k),l)
        co2i2cggat(k,l)  = co2i2cgrow(ilmos(k),jlmos(k),l)
        co2i2csgat(k,l)  = co2i2csrow(ilmos(k),jlmos(k),l)
        ancsveggat(k,l)  = ancsvegrow(ilmos(k),jlmos(k),l)
        ancgveggat(k,l)  = ancgvegrow(ilmos(k),jlmos(k),l)
        rmlcsveggat(k,l) = rmlcsvegrow(ilmos(k),jlmos(k),l)
        rmlcgveggat(k,l) = rmlcgvegrow(ilmos(k),jlmos(k),l)
        anveggat(k,l)    = anvegrow(ilmos(k),jlmos(k),l)
        rmlveggat(k,l)   = rmlvegrow(ilmos(k),jlmos(k),l)
        pfcancmxgat(k,l) = pfcancmxrow(ilmos(k),jlmos(k),l)
        nfcancmxgat(k,l) = nfcancmxrow(ilmos(k),jlmos(k),l)
        stemmassgat(k,l) = stemmassrow(ilmos(k),jlmos(k),l)
        stemmassgat_ns(k,l) = stemmassnsrow(ilmos(k),jlmos(k),l)
        stemmassgat_s(k,l)  = stemmasssrow(ilmos(k),jlmos(k),l)
        rootmassgat(k,l) = rootmassrow(ilmos(k),jlmos(k),l)
        rootmassgat_ns(k,l) = rootmassnsrow(ilmos(k),jlmos(k),l)
        rootmassgat_s(k,l)  = rootmasssrow(ilmos(k),jlmos(k),l)
        pstemmassgat(k,l) = pstemmassrow(ilmos(k),jlmos(k),l)
        pgleafmassgat(k,l) = pgleafmassrow(ilmos(k),jlmos(k),l)
        gleafmasgat(k,l) = gleafmasrow(ilmos(k),jlmos(k),l)
        gleafmasgat_ns(k,l) = gleafmasnsrow(ilmos(k),jlmos(k),l)
        gleafmasgat_s(k,l)  = gleafmassrow(ilmos(k),jlmos(k),l)
        bleafmasgat(k,l) = bleafmasrow(ilmos(k),jlmos(k),l)
        flhrlossgat(k,l) = flhrlossrow(ilmos(k),jlmos(k),l)
        flhrloss_nsgat(k,l) = flhrloss_nsrow(ilmos(k),jlmos(k),l)
        flhrloss_sgat(k,l) = flhrloss_srow(ilmos(k),jlmos(k),l)
        pandaysgat(k,l)  = pandaysrow(ilmos(k),jlmos(k),l)
        lfstatusgat(k,l) = lfstatusrow(ilmos(k),jlmos(k),l)
        grwtheffgat(k,l) = grwtheffrow(ilmos(k),jlmos(k),l)
        lystmmasgat(k,l) = lystmmasrow(ilmos(k),jlmos(k),l)
        lyrotmasgat(k,l) = lyrotmasrow(ilmos(k),jlmos(k),l)
        lmaxtgat(k,l) = lmaxtrow(ilmos(k),jlmos(k),l)
        smaxtgat(k,l) = smaxtrow(ilmos(k),jlmos(k),l)
        rmaxtgat(k,l) = rmaxtrow(ilmos(k),jlmos(k),l)
        lygleafmasmaxgat(k,l) = lygleafmasmaxrow(ilmos(k),jlmos(k),l)
        lystemmassmaxgat(k,l) = lystemmassmaxrow(ilmos(k),jlmos(k),l)
        lyrootmassmaxgat(k,l) = lyrootmassmaxrow(ilmos(k),jlmos(k),l)
        tymaxlaigat(k,l) = tymaxlairow(ilmos(k),jlmos(k),l)
        stmhrlosgat(k,l) = stmhrlosrow(ilmos(k),jlmos(k),l)
        rothrlosgat(k,l) = rothrlosrow(ilmos(k),jlmos(k),l)
        tltrleafgat(k,l) = tltrleafrow(ilmos(k),jlmos(k),l)
        tltrstemgat(k,l) = tltrstemrow(ilmos(k),jlmos(k),l)
        tltrrootgat(k,l) = tltrrootrow(ilmos(k),jlmos(k),l)
        leaflitrgat(k,l) = leaflitrrow(ilmos(k),jlmos(k),l)
        roottempgat(k,l) = roottemprow(ilmos(k),jlmos(k),l)
        afrleafgat(k,l)  = afrleafrow(ilmos(k),jlmos(k),l)
        afrstemgat(k,l)  = afrstemrow(ilmos(k),jlmos(k),l)
        afrrootgat(k,l)  = afrrootrow(ilmos(k),jlmos(k),l)
        wtstatusgat(k,l) = wtstatusrow(ilmos(k),jlmos(k),l)
        ltstatusgat(k,l) = ltstatusrow(ilmos(k),jlmos(k),l)
        nppveggat(k,l)   = nppvegrow(ilmos(k),jlmos(k),l)
        rmlvegaccgat(k,l) = rmlvegaccrow(ilmos(k),jlmos(k),l)
        rmsveggat(k,l)   = rmsvegrow(ilmos(k),jlmos(k),l)
        rmrveggat(k,l)   = rmrvegrow(ilmos(k),jlmos(k),l)
        rgveggat(k,l)    = rgvegrow(ilmos(k),jlmos(k),l)
        gppveggat(k,l)   = gppvegrow(ilmos(k),jlmos(k),l)
        vcmax0gat(k,l)   = vcmax0row(ilmos(k),jlmos(k),l)
        leafns2sgat(k,l) = leafns2srow(ilmos(k),jlmos(k),l)
        stemns2sgat(k,l) = stemns2srow(ilmos(k),jlmos(k),l)
        rootns2sgat(k,l) = rootns2srow(ilmos(k),jlmos(k),l)
        re_alloc_s2lgat(k,l)  = re_alloc_s2lrow(ilmos(k),jlmos(k),l)
        re_alloc_r2lgat(k,l)  = re_alloc_r2lrow(ilmos(k),jlmos(k),l)
        re_alloc_sr2lgat(k,l) = re_alloc_sr2lrow(ilmos(k),jlmos(k),l)
        autoresveggat(k,l) = autoresvegrow(ilmos(k),jlmos(k),l)
        vgbiomas_veggat(k,l) = vgbiomas_vegrow(ilmos(k),jlmos(k),l)
        pftexistgat(k,l) = pftexistrow(ilmos(k),jlmos(k),l)
        ccgat(k,l) = ccrow(ilmos(k),jlmos(k),l)
        mmgat(k,l) = mmrow(ilmos(k),jlmos(k),l)
        litrfallveggat(k,l)   = litrfallvegrow(ilmos(k),jlmos(k),l)
        ngleafmasgat(k,l)     = ngleafmasrow(ilmos(k),jlmos(k),l)
        ngleafmasgat_ns(k,l)  = ngleafmasnsrow(ilmos(k),jlmos(k),l)
        ngleafmasgat_s(k,l)   = ngleafmassrow(ilmos(k),jlmos(k),l)
        nbleafmasgat(k,l)     = nbleafmasrow(ilmos(k),jlmos(k),l)
        nstemmassgat(k,l)     = nstemmassrow(ilmos(k),jlmos(k),l)
        nstemmassgat_ns(k,l)  = nstemmassnsrow(ilmos(k),jlmos(k),l)
        nstemmassgat_s(k,l)   = nstemmasssrow(ilmos(k),jlmos(k),l)
        nrootmassgat(k,l)     = nrootmassrow(ilmos(k),jlmos(k),l)
        nrootmassgat_ns(k,l)  = nrootmassnsrow(ilmos(k),jlmos(k),l)
        nrootmassgat_s(k,l)   = nrootmasssrow(ilmos(k),jlmos(k),l)
        nvgbiomas_veggat(k,l) = nvgbiomas_vegrow(ilmos(k),jlmos(k),l)
        ndemandveg_wp_nppgat(k,l) = ndemandveg_wp_npprow(ilmos(k),jlmos(k),l)
        c2nveg_lgat(k,l)  = c2nveg_lrow(ilmos(k),jlmos(k),l)
        c2nveg_sgat(k,l)  = c2nveg_srow(ilmos(k),jlmos(k),l)
        c2nveg_rgat(k,l)  = c2nveg_rrow(ilmos(k),jlmos(k),l)
        c2nveg_wpgat(k,l) = c2nveg_wprow(ilmos(k),jlmos(k),l)
        nuptakeveg_p_nh4gat(k,l) = nuptakeveg_p_nh4row(ilmos(k),jlmos(k),l)
        nuptakeveg_p_no3gat(k,l) = nuptakeveg_p_no3row(ilmos(k),jlmos(k),l)
        nuptakeveg_a_actl_nh4gat(k,l) = nuptakeveg_a_actl_nh4row(ilmos(k),jlmos(k),l)
        nuptakeveg_a_actl_no3gat(k,l) = nuptakeveg_a_actl_no3row(ilmos(k),jlmos(k),l)
        nuptakeveggat(k,l) = nuptakevegrow(ilmos(k),jlmos(k),l)
        nallocveg_lgat(k,l) = nallocveg_lrow(ilmos(k),jlmos(k),l)
        nallocveg_sgat(k,l) = nallocveg_srow(ilmos(k),jlmos(k),l)
        nallocveg_rgat(k,l) = nallocveg_rrow(ilmos(k),jlmos(k),l)
        nresorpedveg_sgat(k,l)   = nresorpedveg_srow(ilmos(k),jlmos(k),l)
        nresorpedveg_rgat(k,l)   = nresorpedveg_rrow(ilmos(k),jlmos(k),l)
        nre_allocveg_s2lgat(k,l) = nre_allocveg_s2lrow(ilmos(k),jlmos(k),l)
        nre_allocveg_r2lgat(k,l) = nre_allocveg_r2lrow(ilmos(k),jlmos(k),l)
        nleafns2sveggat(k,l) = nleafns2svegrow(ilmos(k),jlmos(k),l)
        nstemns2sveggat(k,l) = nstemns2svegrow(ilmos(k),jlmos(k),l)
        nrootns2sveggat(k,l) = nrootns2svegrow(ilmos(k),jlmos(k),l)
        nlitrveg_lgat(k,l)   = nlitrveg_lrow(ilmos(k),jlmos(k),l)
        nlitrveg_sgat(k,l)   = nlitrveg_srow(ilmos(k),jlmos(k),l)
        nlitrveg_rgat(k,l)   = nlitrveg_rrow(ilmos(k),jlmos(k),l)
        nlitrveggat(k,l)   = nlitrvegrow(ilmos(k),jlmos(k),l)
        gl2bl_grass_nfluxgat(k,l) = gl2bl_grass_nfluxrow(ilmos(k),jlmos(k),l)
        redcoeff_vcmaxgat(k,l) = redcoeff_vcmaxrow(ilmos(k),jlmos(k),l)
        !         fire emission variables
        emit_co2gat(k,l)  = emit_co2row(ilmos(k),jlmos(k),l)
        emit_cogat(k,l)   = emit_corow(ilmos(k),jlmos(k),l)
        emit_ch4gat(k,l)  = emit_ch4row(ilmos(k),jlmos(k),l)
        emit_nmhcgat(k,l) = emit_nmhcrow(ilmos(k),jlmos(k),l)
        emit_h2gat(k,l)   = emit_h2row(ilmos(k),jlmos(k),l)
        emit_noxgat(k,l)  = emit_noxrow(ilmos(k),jlmos(k),l)
        emit_n2ogat(k,l)  = emit_n2orow(ilmos(k),jlmos(k),l)
        emit_nh3gat(k,l)  = emit_nh3row(ilmos(k),jlmos(k),l)
        emit_pm25gat(k,l) = emit_pm25row(ilmos(k),jlmos(k),l)
        emit_tpmgat(k,l)  = emit_tpmrow(ilmos(k),jlmos(k),l)
        emit_tcgat(k,l)   = emit_tcrow(ilmos(k),jlmos(k),l)
        emit_ocgat(k,l)   = emit_ocrow(ilmos(k),jlmos(k),l)
        emit_bcgat(k,l)   = emit_bcrow(ilmos(k),jlmos(k),l)
        burnvegfgat(k,l)  = burnvegfrow(ilmos(k),jlmos(k),l)
        bnfantgat(k,l)   = bnfantrow(ilmos(k),jlmos(k),l)
        bnfnatgat(k,l)   = bnfnatrow(ilmos(k),jlmos(k),l)
        nstressgat(k,l)  = nstressrow(ilmos(k),jlmos(k),l)
        !
      end do
    end do ! loop 101
    !
    do l = 1,iccp1 ! loop 102
      do k = 1,nml
        hetroresveggat(k,l) = hetroresvegrow(ilmos(k),jlmos(k),l)
        nepveggat(k,l)   = nepvegrow(ilmos(k),jlmos(k),l)
        nbpveggat(k,l)   = nbpvegrow(ilmos(k),jlmos(k),l)
        nh4_massgat(k,l)    = nh4_massrow(ilmos(k),jlmos(k),l)
        no3_massgat(k,l)    = no3_massrow(ilmos(k),jlmos(k),l)
        nitrifveggat(k,l)   = nitrifvegrow(ilmos(k),jlmos(k),l)
        no_nitveggat(k,l)   = no_nitvegrow(ilmos(k),jlmos(k),l)
        no_denitveggat(k,l) = no_denitvegrow(ilmos(k),jlmos(k),l)
        no_nitdenitveggat(k,l)= no_nitdenitvegrow(ilmos(k),jlmos(k),l)
        n2o_nitveggat(k,l)  = n2o_nitvegrow(ilmos(k),jlmos(k),l)
        n2o_denitveggat(k,l) = n2o_denitvegrow(ilmos(k),jlmos(k),l)
        n2o_nitdenitveggat(k,l) = n2o_nitdenitvegrow(ilmos(k),jlmos(k),l)
        n2_denitveggat(k,l) = n2_denitvegrow(ilmos(k),jlmos(k),l)
        nvolveggat(k,l)     = nvolvegrow(ilmos(k),jlmos(k),l)
        nleachveggat(k,l)   = nleachvegrow(ilmos(k),jlmos(k),l)
        appl_fertgat(k,l)   = appl_fertrow(ilmos(k),jlmos(k),l)
        ndep_nh4gat(k,l)    = ndep_nh4row(ilmos(k),jlmos(k),l)
        ndep_no3gat(k,l)    = ndep_no3row(ilmos(k),jlmos(k),l)
        c2nveg_litrgat(k,l) = c2nveg_litrrow(ilmos(k),jlmos(k),l)
        c2nveg_humusgat(k,l) = c2nveg_humusrow(ilmos(k),jlmos(k),l)
        nhumtrsveggat(k,l)  = nhumtrsvegrow(ilmos(k),jlmos(k),l)
        nmineralveg_litrgat(k,l) = nmineralveg_litrrow(ilmos(k),jlmos(k),l)
        nmineralveg_humusgat(k,l) = nmineralveg_humusrow(ilmos(k),jlmos(k),l)
        netnmineralveg_gat(k,l) = netnmineralveg_row(ilmos(k),jlmos(k),l)
        nimmobilveg_nh4gat(k,l) = nimmobilveg_nh4row(ilmos(k),jlmos(k),l)
        nimmobilveg_no3gat(k,l) = nimmobilveg_no3row(ilmos(k),jlmos(k),l)
        fNnetlandveggat(k,l) = fNnetlandvegrow(ilmos(k),jlmos(k),l)
        bnftotgat(k,l)   = bnftotrow(ilmos(k),jlmos(k),l)
        bnffreegat(k,l)  = bnffreerow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 102

    do l = 1,iccp2 ! loop 103
      do k = 1,nml
        do m = 1,ignd
          litrmassgat(k,l,m)=litrmassrow(ilmos(k),jlmos(k),l,m)
          soilcmasgat(k,l,m)=soilcmasrow(ilmos(k),jlmos(k),l,m)
          litresveggat(k,l,m)= litresvegrow(ilmos(k),jlmos(k),l,m)
          soilcresveggat(k,l,m)= soilcresvegrow(ilmos(k),jlmos(k),l,m)
        end do
        nlitrmassgat(k,l)   = nlitrmassrow(ilmos(k),jlmos(k),l)
        soilnmasgat(k,l)    = soilnmasrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 103
    !
    do l = 1,ignd ! loop 250
      do k = 1,nml
        sandgat(k,l) = sandrow(ilmos(k),jlmos(k),l)
        claygat(k,l) = clayrow(ilmos(k),jlmos(k),l)
        orgmgat(k,l) = orgmrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 250
    !
    do k = 1,nml ! loop 300
      ! If peatlands are not present these just gather regardless.
      anmossgat(k) =  anmossrow(ilmos(k),jlmos(k))
      rmlmossgat(k) = rmlmossrow(ilmos(k),jlmos(k))
      gppmossgat(k) = gppmossrow(ilmos(k),jlmos(k))
      armossgat(k) =  armossrow(ilmos(k),jlmos(k))
      nppmossgat(k) = nppmossrow(ilmos(k),jlmos(k))
      peatdepgat(k)    =  peatdeprow(ilmos(k),jlmos(k))
      litrmsmossgat(k) =  litrmsmossrow(ilmos(k),jlmos(k))
      Cmossmasgat(k) = Cmossmasrow(ilmos(k),jlmos(k))
      dmossgat(k) = dmossrow(ilmos(k),jlmos(k))
      peatSoilCgat(k) = peatSoilCrow(ilmos(k),jlmos(k))
      pddgat(k)   = pddrow(ilmos(k),jlmos(k))
    end do ! loop 300
    return
    
  end subroutine ctemg2
  !! @}
  ! ------------------------------------------------------------------------------------

  !> \namespace ctemgatherscatter
  !! Transfers information between the 'gathered' and 'scattered' form of the CTEM data arrays.
  !!

end module ctemGatherScatter
