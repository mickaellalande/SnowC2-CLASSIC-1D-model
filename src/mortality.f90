!> \file
!> Calculates the litter generated from leaves, stem, and root components after
!! vegetation dies due to reduced growth efficiency or due to aging (the intrinsic mortality).
!! Also updates the vegetation pools for changes due to the mortality calculations.
module mortality

  implicit none

  public :: mortalty
  public :: updatePoolsMortality

contains

  !> \ingroup mortality_mortalty
  !! @{
  !> Calculates mortality due to poor plant fitness and aging. Determines the litter generated from these processes.
  !> @author Vivek Arora and Joe Melton
  subroutine mortalty (stemmass, stemmass_ns, stemmass_s, & ! In
                       rootmass, rootmass_ns, rootmass_s, & ! In
                       ailcg, gleafmas, gleafmas_ns, gleafmas_s, & ! In
                       nstemmass, nstemmass_ns, nstemmass_s, & ! In
                       nrootmass, nrootmass_ns, nrootmass_s, & ! In
                       ngleafmas, ngleafmas_ns, ngleafmas_s, & ! In
                       il1, il2, ilg, & ! In
                       leapnow, iday, sort, fcancmx, ipeatland, & ! In
                       useTracer, tracerStemMass, tracerRootMass, tracerGLeafMass, & ! In
                       lystmmas, lyrotmas, tymaxlai, grwtheff, & ! In/Out
                       stemltrm, stemltrm_ns, stemltrm_s, & ! Out
                       rootltrm, rootltrm_ns, rootltrm_s, & ! Out   
                       glealtrm, glealtrm_ns, glealtrm_s, geremort, & ! Out
                       intrmort, tracerStemMort, tracerRootMort, tracerGLeafMort, & ! Out
                       stemltrm_n, stemltrm_ns_n, stemltrm_s_n, & ! Out
                       rootltrm_n, rootltrm_ns_n, rootltrm_s_n, & ! Out
                       glealtrm_n, glealtrm_ns_n, glealtrm_s_n) ! Out
    !
    !     07  Dec 2018  - Pass ilg back in as an argument
    !     V. Arora
    !
    !     17  Jan 2014  - Moved parameters to global file (classicParams.f90)
    !     J. Melton
    !
    !     22  Jul 2013  - Add in module for parameters
    !     J. Melton
    !
    !     24  sep 2012  - add in checks to prevent calculation of non-present
    !     J. Melton       pfts
    !
    !     07  may 2003  - this subroutine calculates the litter generated
    !     V. Arora        from leaves, stem, and root components after
    !                     vegetation dies due to reduced growth efficiency
    !                     or due to aging (the intrinsic mortality)
    !
    use classicParams, only : icc, kk, zero, mxmortge, kmort1, maxage, maxage_peat

    use ctemStatevars, only : c_switch

    implicit none

    integer, intent(in) :: ilg    !< no. of grid cells in latitude circle
    integer, intent(in) :: il1    !< il1=1
    integer, intent(in) :: il2    !< il2=ilg
    integer, intent(in) :: iday   !< day of the year
    integer, intent(in) :: ipeatland(ilg) !< Peatland flag, non-peatlands are 0.
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    logical, intent(in) :: leapnow   !< true if this year is a leap year. Only used if the switch 'leap' is true.
    real, intent(in) :: fcancmx(ilg,icc)  !<
    integer, intent(in) :: sort(icc) !< index for correspondence between ctem 9 pfts and size 12 of parameters vectors
    real, intent(in) :: stemmass(ilg,icc) !< stem mass for each of the 9 ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: stemmass_ns(ilg,icc) !< non-structural stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: stemmass_s(ilg,icc)  !< structural stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: rootmass(ilg,icc) !< root mass for each of the 9 ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: rootmass_ns(ilg,icc) !< non-structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: rootmass_s(ilg,icc)  !< structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: gleafmas(ilg,icc) !< green leaf mass for each of the 9 ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: gleafmas_ns(ilg,icc) !< non-structural green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, intent(in) :: gleafmas_s(ilg,icc)  !< structural green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, intent(in) :: ngleafmas(ilg,icc)      !< N green leaf mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: ngleafmas_ns(ilg,icc)   !< N non-structural green leaf mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: ngleafmas_s(ilg,icc)    !< N structural green leaf mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: nstemmass(ilg,icc)      !< N stem mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: nstemmass_ns(ilg,icc)   !< N non-structural stem mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: nstemmass_s(ilg,icc)    !< N structural stem mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: nrootmass(ilg,icc)      !< N root mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: nrootmass_ns(ilg,icc)   !< N non-structural root mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(in) :: nrootmass_s(ilg,icc)    !< N structural root mass for individual PFTs (\f$g N/m^2\f)
    real, intent(in) :: tracerStemMass(:,:)  !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerRootMass(:,:)  !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerGLeafMass(:,:) !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(in) :: ailcg(ilg,icc)    !< green or live lai
    real, intent(inout) :: grwtheff(ilg,icc) !< growth efficiency. change in biomass per year per unit max. lai (g c/m2)/(m2/m2)
    real, intent(inout) :: lystmmas(ilg,icc) !< stem mass at the end of last year
    real, intent(inout) :: lyrotmas(ilg,icc) !< root mass at the end of last year
    real, intent(inout) :: tymaxlai(ilg,icc) !< this year's maximum lai
    real, intent(out) :: stemltrm(ilg,icc) !< stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: stemltrm_ns(ilg,icc) !< non-structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: stemltrm_s(ilg,icc)  !< structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: rootltrm(ilg,icc) !< root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: rootltrm_ns(ilg,icc) !< non-structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: rootltrm_s(ilg,icc)  !< structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: glealtrm(ilg,icc) !< green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: glealtrm_ns(ilg,icc) !< non-structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: glealtrm_s(ilg,icc)  !< structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: stemltrm_n(ilg,icc) !< N stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: stemltrm_ns_n(ilg,icc) !< N non-structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: stemltrm_s_n(ilg,icc)  !< N structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: rootltrm_n(ilg,icc) !< N root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: rootltrm_ns_n(ilg,icc) !< N non-structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: rootltrm_s_n(ilg,icc)  !< N structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: glealtrm_n(ilg,icc) !< N green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: glealtrm_ns_n(ilg,icc) !< N non-structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: glealtrm_s_n(ilg,icc)  !< N structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(out) :: geremort(ilg,icc) !< growth efficiency related mortality (1/day)
    real, intent(out) :: intrmort(ilg,icc) !< intrinsic mortality (1/day)
    real, intent(out) :: tracerStemMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$
    real, intent(out) :: tracerRootMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$
    real, intent(out) :: tracerGLeafMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$


    integer :: i, j, k, n
    real :: tempVarMort

    associate( &
     Ncycle_on => c_switch%Ncycle_on                     & !< logical: true if nitrogen cycle is on.
    )

    ! Initialize required arrays to zero
    stemltrm = 0.0
    stemltrm_ns = 0.0
    stemltrm_s = 0.0
    rootltrm = 0.0
    rootltrm_ns = 0.0
    rootltrm_s = 0.0
    glealtrm = 0.0
    glealtrm_ns = 0.0
    glealtrm_s = 0.0
    stemltrm_n = 0.0
    stemltrm_ns_n = 0.0
    stemltrm_s_n = 0.0
    rootltrm_n = 0.0
    rootltrm_ns_n = 0.0
    rootltrm_s_n = 0.0
    glealtrm_n = 0.0
    glealtrm_ns_n = 0.0
    glealtrm_s_n = 0.0
    geremort = 0.0
    intrmort = 0.0
    tracerStemMort = 0.0
    tracerRootMort = 0.0
    tracerGLeafMort = 0.0

    ! ------------------------------------------------------------------

    !> At the end of every year, i.e. when iday equals 365/366, we calculate
    !! growth related mortality. rather than using this number to kill
    !! plants at the end of every year, this mortality rate is applied
    !! gradually over the next year.
    do j = 1,icc ! loop 200
      n = sort(j)
      do i = il1,il2
        if (fcancmx(i,j) > 0.0) then

          if (iday == 1) tymaxlai(i,j) = 0.0

          if (ailcg(i,j) > tymaxlai(i,j)) tymaxlai(i,j) = ailcg(i,j)

          if ((.not. leapnow .and. iday == 365) .or. (leapnow .and. iday == 366)) then
            if (tymaxlai(i,j) > zero) then
              grwtheff(i,j) = ((stemmass(i,j) + rootmass(i,j)) &
                              - (lystmmas(i,j) + lyrotmas(i,j))) / tymaxlai(i,j)
            else
              grwtheff(i,j) = 0.0
            end if
            grwtheff(i,j) = max(0.0,grwtheff(i,j)) * 1000.0
            
            lystmmas(i,j) = stemmass(i,j)
            lyrotmas(i,j) = rootmass(i,j)
          end if

          !> Calculate growth related mortality using last year's growth
          !! efficiency or the new growth efficiency if day is 365 and
          !! growth efficiency estimate has been updated above.
          geremort(i,j) = mxmortge(n) / (1.0 + kmort1 * grwtheff(i,j))

          ! convert (1/year) rate into (1/day) rate
          if (leapnow) then
            geremort(i,j) = geremort(i,j) / 366.0
          else
            geremort(i,j) = geremort(i,j) / 365.0
          end if
        end if
      end do ! loop 210
    end do ! loop 200

    !> Calculate intrinsic mortality rate due to aging which implicity includes effects of frost,
    !! hail,wind throw etc. it is assumed that only 1% of the plants exceed maximum age (which is
    !! a pft-dependent parameter). to achieve this some fraction of the plants need to be killed every year.
    do j = 1,icc ! loop 250
      n = sort(j)
      do i = il1,il2
        if (fcancmx(i,j) > 0.0) then
          if (maxage(n) > zero) then
            if (ipeatland(i) /= 1 .and. ipeatland(i) /= 2) then ! Uplands
              intrmort(i,j) = 1.0 - exp( - 4.605 / maxage(n))
            else !peatlands
              intrmort(i,j) = 1.0 - exp( - 4.605 / maxage_peat(n))
            end if 
          else
            intrmort(i,j) = 0.0
          end if

          !> convert (1/year) rate into (1/day) rate
          if (leapnow) then
            intrmort(i,j) = intrmort(i,j) / 366.0
          else
            intrmort(i,j) = intrmort(i,j) / 365.0
          end if
        end if
      end do ! loop 260
    end do ! loop 250

    !> Now that we have both growth related and intrinsic mortality rates,
    !! lets combine these rates for every pft and estimate litter generated
    do j = 1,icc ! loop 300
      do i = il1,il2
        if (fcancmx(i,j) > 0.0) then
          tempVarMort = (1.0 - exp( - 1.0 * (geremort(i,j) + intrmort(i,j))))
          stemltrm(i,j) = stemmass(i,j) * tempVarMort
          stemltrm_ns(i,j) = stemmass_ns(i,j) * tempVarMort
          stemltrm_s(i,j) = stemmass_s(i,j) * tempVarMort
          rootltrm(i,j) = rootmass(i,j) * tempVarMort
          rootltrm_ns(i,j) = rootmass_ns(i,j) * tempVarMort
          rootltrm_s(i,j) = rootmass_s(i,j) * tempVarMort
          glealtrm(i,j) = gleafmas(i,j) * tempVarMort
          glealtrm_ns(i,j) = gleafmas_ns(i,j) * tempVarMort
          glealtrm_s(i,j) = gleafmas_s(i,j) * tempVarMort
          if (useTracer > 0) then
            tracerStemMort(i,j) = tracerStemMass(i,j) * tempVarMort
            tracerRootMort(i,j) = tracerRootMass(i,j) * tempVarMort
            tracerGLeafMort(i,j) = tracerGLeafMass(i,j) * tempVarMort
          end if
          if (Ncycle_on) then
           stemltrm_n(i,j) = nstemmass(i,j) * tempVarMort
           stemltrm_ns_n(i,j) = nstemmass_ns(i,j) * tempVarMort
           stemltrm_s_n(i,j) = nstemmass_s(i,j) * tempVarMort
           rootltrm_n(i,j) = nrootmass(i,j) * tempVarMort
           rootltrm_ns_n(i,j) = nrootmass_ns(i,j) * tempVarMort
           rootltrm_s_n(i,j) = nrootmass_s(i,j) * tempVarMort
           glealtrm_n(i,j) = ngleafmas(i,j) * tempVarMort
           glealtrm_ns_n(i,j) = ngleafmas_ns(i,j) * tempVarMort
           glealtrm_s_n(i,j) = ngleafmas_s(i,j) * tempVarMort
          end if
        end if
      end do ! loop 310
    end do ! loop 300

    end associate

  end subroutine mortalty
  !! @}

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup mortality_updatepoolsmortality
  !! @{
  !> Update leaf, stem, and root biomass pools to take into loss due to mortality, and put the
  !! litter into the litter pool. The mortality for green grasses doesn't generate litter, instead they turn brown.
  !> @author Vivek Arora and Joe Melton
  subroutine updatePoolsMortality (il1, il2, ilg, stemltrm, stemltrm_ns, & ! In
                                   stemltrm_s, rootltrm, rootltrm_ns, rootltrm_s, & ! In
                                   stemltrm_n, stemltrm_ns_n, & ! In
                                   stemltrm_s_n, rootltrm_n, rootltrm_ns_n, rootltrm_s_n, & ! In
                                   useTracer, & ! In
                                   rmatctem, tracerStemMort, tracerRootMort, tracerGLeafMort, & ! In
                                   stemmass, stemmass_ns, stemmass_s, & ! In/Out
                                   rootmass, rootmass_ns, rootmass_s, litrmass, & ! In/Out
                                   glealtrm, glealtrm_ns, glealtrm_s, & ! In/Out
                                   gleafmas, gleafmas_ns, gleafmas_s, bleafmas, & ! In/Out
                                   nstemmass, nstemmass_ns, nstemmass_s, & ! In/Out
                                   nrootmass, nrootmass_ns, nrootmass_s, nlitrmass, & ! In/Out
                                   glealtrm_n, glealtrm_ns_n, glealtrm_s_n, & ! In/Out
                                   ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! In/Out
                                   tracerLitrMass, & ! In/Out
                                   tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass) ! In/Out

    use classicParams, only : ican, nol2pfts, classpfts, ignd, icc, reindexPFTs

    use ctemStatevars, only : c_switch

    implicit none

    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg             !< no. of grid cells/tiles in latitude circle
    real, intent(in)    :: stemltrm(:,:)   !< stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: stemltrm_ns(ilg,icc) !< non-structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: stemltrm_s(ilg,icc)  !< structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: rootltrm(:,:)   !< root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: rootltrm_ns(ilg,icc) !< non-structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: rootltrm_s(ilg,icc)  !< structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: stemltrm_n(:,:)   !< N stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: stemltrm_ns_n(ilg,icc) !< N non-structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: stemltrm_s_n(ilg,icc)  !< N structural stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: rootltrm_n(:,:)   !< N root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: rootltrm_ns_n(ilg,icc) !< N non-structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: rootltrm_s_n(ilg,icc)  !< N structural root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in)    :: rmatctem(:,:,:) !< fraction of roots for each of ctem's 9 pfts in each soil layer
    real, intent(in) :: tracerStemMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$
    real, intent(in) :: tracerRootMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.

    real, intent(inout) :: glealtrm(:,:)   !< green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(inout) :: glealtrm_ns(:,:) !< non-structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(inout) :: glealtrm_s(:,:)  !< structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(inout) :: glealtrm_n(:,:)   !< N green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(inout) :: glealtrm_ns_n(:,:) !< N non-structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(inout) :: glealtrm_s_n(:,:)  !< N structural green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(inout) :: stemmass(:,:)   !< stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: stemmass_ns(:,:) !< non-structural stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: stemmass_s(:,:)  !< structural stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: rootmass(:,:)   !< root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: rootmass_ns(:,:) !< non-structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: rootmass_s(:,:)  !< structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: litrmass(:,:,:) !< litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(inout) :: gleafmas(:,:)   !< green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: gleafmas_ns(:,:) !< non-structural green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, intent(inout) :: gleafmas_s(:,:)  !< structural green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$ 
    real, intent(inout) :: bleafmas(:,:)   !< brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: nstemmass(:,:)   !< N stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: nstemmass_ns(:,:) !< N non-structural stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: nstemmass_s(:,:)  !< N structural stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: nrootmass(:,:)   !< N root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: nrootmass_ns(:,:) !< N non-structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: nrootmass_s(:,:)  !< N structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: nlitrmass(:,:) !< N litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(inout) :: ngleafmas(:,:)   !< N green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: ngleafmas_ns(:,:) !< N non-structural green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, intent(inout) :: ngleafmas_s(:,:)  !< N structural green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$ 
    real, intent(inout) :: nbleafmas(:,:)   !< N brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: tracerStemMass(:,:)  !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(inout) :: tracerRootMass(:,:)  !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(inout) :: tracerGLeafMass(:,:) !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerBLeafMass(:,:) !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout) :: tracerGLeafMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$

    integer :: k1, j, m, k2, i, k

    associate( &
    Ncycle_on => c_switch%Ncycle_on                     & !< logical: true if nitrogen cycle is on.
    )

    !> Update leaf, stem, and root biomass pools to take into loss due to mortality, and put the
    !! litter into the litter pool. the mortality for green grasses doesn't generate litter, instead they turn brown.
    do j = 1,ican
      do m = reindexPFTs(j,1),reindexPFTs(j,2)
        do i = il1,il2

          !stemmass(i,m) = stemmass(i,m) - stemltrm(i,m)
          stemmass_ns(i,m) = stemmass_ns(i,m) - stemltrm_ns(i,m)
          stemmass_s(i,m) = stemmass_s(i,m) - stemltrm_s(i,m)
          stemmass(i,m) = stemmass_ns(i,m) + stemmass_s(i,m)
          
          !rootmass(i,m) = rootmass(i,m) - rootltrm(i,m)
          rootmass_ns(i,m) = rootmass_ns(i,m) - rootltrm_ns(i,m)
          rootmass_s(i,m) = rootmass_s(i,m) - rootltrm_s(i,m)
          rootmass(i,m) = rootmass_ns(i,m) + rootmass_s(i,m)

          if (useTracer > 0) then ! Update the tracer pools if being used.
            tracerStemMass(i,m) = tracerStemMass(i,m) - tracerStemMort(i,m)
            tracerRootMass(i,m) = tracerRootMass(i,m) - tracerRootMort(i,m)
          end if

          if (Ncycle_on) then
            nstemmass(i,m) = nstemmass(i,m) - stemltrm_n(i,m)
            nstemmass_ns(i,m) = nstemmass_ns(i,m) - stemltrm_ns_n(i,m)
            nstemmass_s(i,m) = nstemmass_s(i,m) - stemltrm_s_n(i,m)
            nrootmass(i,m) = nrootmass(i,m) - rootltrm_n(i,m)
            nrootmass_ns(i,m) = nrootmass_ns(i,m) - rootltrm_ns_n(i,m)
            nrootmass_s(i,m) = nrootmass_s(i,m) - rootltrm_s_n(i,m)
          end if

          select case (classpfts(j))

          case ('NdlTr' , 'BdlTr', 'Crops', 'BdlSh')

            ! gleafmas(i,m) = gleafmas(i,m) - glealtrm(i,m)
            gleafmas_ns(i,m) = gleafmas_ns(i,m) - glealtrm_ns(i,m)
            gleafmas_s(i,m) = gleafmas_s(i,m) - glealtrm_s(i,m)

            ! refind the total
            gleafmas(i,m) = gleafmas_ns(i,m) + gleafmas_s(i,m)

            if (useTracer > 0) tracerGLeafMass(i,m) = tracerGLeafMass(i,m) - tracerGLeafMort(i,m)

            if (Ncycle_on) then
             ngleafmas(i,m) = ngleafmas(i,m) - glealtrm_n(i,m)
             ngleafmas_ns(i,m) = ngleafmas_ns(i,m) - glealtrm_ns_n(i,m)
             ngleafmas_s(i,m) = ngleafmas_s(i,m) - glealtrm_s_n(i,m)
            end if

          case ('Grass')    ! grasses

            !gleafmas(i,m) = gleafmas(i,m) - glealtrm(i,m)
            gleafmas_ns(i,m) = gleafmas_ns(i,m) - glealtrm_ns(i,m)
            gleafmas_s(i,m) = gleafmas_s(i,m) - glealtrm_s(i,m)
            ! refind the total
            gleafmas(i,m) = gleafmas_ns(i,m) + gleafmas_s(i,m)
            
            bleafmas(i,m) = bleafmas(i,m) + glealtrm(i,m)
            glealtrm(i,m) = 0.0
            glealtrm_ns(i,m) = 0.0
            glealtrm_s(i,m) = 0.0

            if (useTracer > 0) then ! Update the tracer pools if being used.
              tracerGLeafMass(i,m) = tracerGLeafMass(i,m) - tracerGLeafMort(i,m)
              tracerBLeafMass(i,m) = tracerBLeafMass(i,m) + tracerGLeafMort(i,m)
              tracerGLeafMort(i,m) = 0.0
            end if

            if (Ncycle_on) then
             ngleafmas(i,m) = ngleafmas(i,m) - glealtrm_n(i,m)
             ngleafmas_ns(i,m) = ngleafmas_ns(i,m) - glealtrm_ns_n(i,m)
             ngleafmas_s(i,m) = ngleafmas_s(i,m) - glealtrm_s_n(i,m)
             nbleafmas(i,m) = nbleafmas(i,m) + glealtrm_n(i,m)
             glealtrm_n(i,m) = 0.0
             glealtrm_ns_n(i,m) = 0.0
             glealtrm_s_n(i,m) = 0.0
            end if

          case default

            print * ,'Unknown CLASS PFT in mortality ',classpfts(j)
            call errorHandler('updatePoolsMortality',1)

          end select

          do k = 1,ignd
            if (k == 1) then
              ! The first layer gets the leaf and stem litter. The root litter is given in proportion
              ! to the root distribution
              litrmass(i,m,k) = litrmass(i,m,k) + stemltrm(i,m) &
                                                + rootltrm(i,m) * rmatctem(i,m,k) &
                                                + glealtrm(i,m)
  
              if (useTracer > 0) tracerLitrMass(i,m,k) = tracerLitrMass(i,m,k) + tracerStemMort(i,m) &
                                                + tracerRootMort(i,m) * rmatctem(i,m,k) &
                                                + tracerGLeafMort(i,m)


            else
              litrmass(i,m,k) = litrmass(i,m,k) + rootltrm(i,m) * rmatctem(i,m,k)
  
              if (useTracer > 0) tracerLitrMass(i,m,k) = tracerLitrMass(i,m,k) + tracerRootMort(i,m) * rmatctem(i,m,k)
 
            end if
          end do 

          if (Ncycle_on) then
            nlitrmass(i,m) = nlitrmass(i,m) + stemltrm_n(i,m) &
                                            + rootltrm_n(i,m) &
                                            + glealtrm_n(i,m)
          end if

        end do ! loop 840
      end do ! loop 835
    end do ! loop 830

  end associate

  end subroutine updatePoolsMortality
  !! @}

  !> \namespace mortality
  !!  Calculates the litter generated by leaves, stem and roots due to mortality and updates the C pools.
  !> @author Vivek Arora, Joe Melton
  !! The PFT-dependent mortality rate (\f$day^{-1}\f$),
  !!
  !! \f[ \label{mortality} m_{\alpha} = m_{intr, \alpha} + m_{ge, \alpha} + m_{bioclim, \alpha} + m_{dist, \alpha}, \qquad (Eqn 1) \f]
  !!
  !! reflects the net effect of four different processes: (1) intrinsic- or age-related mortality, \f$m_{intr}\f$, (2)
  !! growth or stress-related mortality, \f$m_{ge}\f$, (3) mortality associated with bioclimatic criteria,
  !! \f$m_{bioclim}\f$ and (4) mortality associated with disturbances, \f$m_{dist}\f$.
  !!
  !! Intrinsic- or age-related mortality uses a PFT-specific maximum age, \f$A_{max}\f$ (see also classicParams.f90),
  !! to calculate an annual mortality rate such that only \f$1\, {\%}\f$ of tree PFTs exceed \f$A_{max}, \alpha\f$.
  !! Intrinsic mortality accounts for processes, whose effect is not explicitly captured in the model including insect damage, hail, wind throw, etc.,
  !!
  !! \f[ \label{intrmort} m_{intr, \alpha} = 1 - \exp(-4.605/A_{max, \alpha}).\qquad (Eqn 2) \f]
  !!
  !! Grasses and crops have \f$m_{intr} = 0\f$. The annual growth-related mortality \f$m_{ge}\f$ is calculated
  !! using growth efficiency of a PFT over the course of the previous year following Prentice et al. (1993)
  !! \cite Prentice1993-xn and Sitch et al. (2003) \cite Sitch2003-847 as
  !!
  !! \f[ \label{mgrow} m_{ge, \alpha} = \frac{m_{{ge}, max, \alpha}}{1 + k_{m} g_{\mathrm{e}, \alpha}}, \qquad (Eqn 3)\f]
  !!
  !! where \f$m_{{ge}, max}\f$ represents the PFT-specific maximum mortality rate when no growth occurs (see also
  !! classicParams.f90). \f$k_{m}\f$ is a parameter set to \f$0.3\, m^{2}\, (g\, C)^{-1}\f$. \f$g_\mathrm{e}\f$
  !! is the growth efficiency of the PFT (\f$g\, C\, m^{-2}\f$) calculated based on the maximum LAI
  !! (\f$L_{\alpha, max}\f$; \f$m^{2}\, m^{-2}\f$) and the increment in stem and root mass over the
  !! course of the previous year (\f$\Delta C_\mathrm{S}\f$ and \f$\Delta C_\mathrm{R}\f$; \f$kg\, C\, m^{-2}\f$, respectively) (Waring, 1983) \cite Waring1983-wc
  !! \f[ g_{\mathrm{e}, \alpha} = 1000\frac{\max(0, (\Delta C_{\mathrm{S}, \alpha}+\Delta C_{\mathrm{R}, \alpha}))}{L_{\alpha, max}}. \qquad (Eqn 4)\f]
  !!
  !! When competition between PFTs is switched on, mortality associated with bioclimatic criteria,
  !! \f$m_{bioclim}\f$ (\f$0.25\, yr^{-1}\f$), is applied when climatic conditions in a grid cell
  !! become unfavourable for a PFT to exist and ensures that PFTs do not exist outside their
  !! bioclimatic envelopes, as explained in competitionScheme :: existence.
  !!
  !! The annual mortality rates for \f$m_{intr}\f$, \f$m_{ge}\f$ and \f$m_{bioclim}\f$ are converted
  !! to daily rates and applied at the daily time step of the model. \f$m_{dist}\f$ is calculated by
  !! the fire module of the model (when switched on) based on daily area burned for each PFT as
  !! summarized in disturbance_scheme :: disturb. In practice, the
  !! \f$\frac{\mathrm{d}f_\alpha}{\mathrm{d}t}=-m_{dist, \alpha}f_\alpha\f$ term of
  !! competitionScheme Eqn 1  is implemented right after area burnt is calculated.
  !!
  !!
end module mortality
