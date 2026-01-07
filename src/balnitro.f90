!>\file
!!Checks nitrogen pools for conservation
!!@author A. Asaadi
!!
!!unless otherwise mentioned all pools are in g.N/m2
!!unless otherwise mentioned all fluxes are in units of g.N/m2.day
!!
subroutine  balnit (il1, il2, iday, leapnow, radj, appl_fert, appl_fert_g, ndeposit, &
                    bnf_free, bnf_free_g, nh4_mass, pnh4_mass,    &
                    nh4_mass_g, pnh4_mass_g, no3_mass, pno3_mass, no3_mass_g, pno3_mass_g, ngleafmas, pngleafmas,         &
                    ngleafmas_g, pngleafmas_g, ngleafmas_ns, pngleafmas_ns, ngleafmas_ns_g, pngleafmas_ns_g, ngleafmas_s, &
                    pngleafmas_s, ngleafmas_s_g, pngleafmas_s_g, nbleafmas, pnbleafmas, nbleafmas_g, pnbleafmas_g,        &
                    nstemmass, pnstemmass, nstemmass_g, pnstemmass_g, nstemmass_ns, pnstemmass_ns, nstemmass_ns_g,        &
                    pnstemmass_ns_g, nstemmass_s, pnstemmass_s, nstemmass_s_g, pnstemmass_s_g, nrootmass, pnrootmass,     &
                    nrootmass_g, pnrootmass_g, nrootmass_ns, pnrootmass_ns, nrootmass_ns_g, pnrootmass_ns_g, nrootmass_s, &
                    pnrootmass_s, nrootmass_s_g, pnrootmass_s_g, nlitrmass, pnlitrmass, nlitrmass_g, pnlitrmass_g,        &
                    soilnmas, psoilnmas, soilnmas_g, psoilnmas_g, nleachveg, nleach, nvolveg, nvol, nitrifveg, nitrif,    &
                    denitveg, denit, no_nitveg, no_nit, n2o_nitveg, n2o_nit, nuptakeveg_nh4, nuptake_nh4, nuptakeveg_no3, &
                    nuptake_no3, ntchlveg, nleafns2sveg, nstemns2sveg, nrootns2sveg, nntchveg_l,                          &
                    nntch_l, nntchveg_s, nntch_s, nntchveg_r, nntch_r,lfstatus, nlitrveg_l, nlitr_l, nlitrveg_s, nlitr_s, &
                    nlitrveg_r, nlitr_r, gl2bl_grass_nflux, gl2bl_grass_nflux_gavg,                                       & 
                    nhumtrsveg, nhumtrs, nmineralveg_litr, nmineral_litr, nmineralveg_humus, nmineral_humus,              &
                    nimmobilveg_nh4, nimmobil_nh4, nimmobilveg_no3, nimmobil_no3)
!--------------------------------------------------------------------------
!
use classicParams, only : ilg, ntolerance, icc, iccp1, iccp2, ctempfts
!
implicit none

! Arguments (in)

integer,intent(in) :: il1  !<il1=1
integer,intent(in) :: il2  !<il2=ilg
integer,intent(in) :: iday !<day of year
logical,intent(in) :: leapnow   !< true if this year is a leap year. Only used if the switch 'leap' is true.
real,dimension(ilg,iccp1),intent(in)    :: nh4_mass         !<pools (after being updated) : ammonium mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nh4_mass_g       !<pools (after being updated) : grid avg. NH4+ mass (\f$g N/m^2\f$)
real,dimension(ilg,iccp1),intent(in)    :: no3_mass         !<pools (after being updated) : nitrate mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: no3_mass_g       !<pools (after being updated) : grid avg. NO3- pool (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: ngleafmas        !<pools (after being updated) : green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: ngleafmas_g      !<pools (after being updated) : grid avg. green leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: ngleafmas_ns     !<pools (after being updated) : non-structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: ngleafmas_ns_g   !<pools (after being updated) : grid avg. non-structural green leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: ngleafmas_s_g    !<pools (after being updated) : grid avg. structural green leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: ngleafmas_s      !<pools (after being updated) : structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: nbleafmas        !<pools (after being updated) : brown leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nbleafmas_g      !<pools (after being updated) : grid avg. brown leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: nstemmass        !<pools (after being updated) : stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nstemmass_g      !<pools (after being updated) : grid avg. stem nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: nstemmass_ns     !<pools (after being updated) : non-structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nstemmass_ns_g   !<pools (after being updated) : grid avg. non-structural stem nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: nstemmass_s      !<pools (after being updated) : structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nstemmass_s_g    !<pools (after being updated) : grid avg. structural stem nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: nrootmass        !<pools (after being updated) : root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nrootmass_g      !<pools (after being updated) : grid avg. root nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: nrootmass_ns     !<pools (after being updated) : non-structural root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nrootmass_ns_g   !<pools (after being updated) : grid avg. non-structural root nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)      :: nrootmass_s      !<pools (after being updated) : structural root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nrootmass_s_g    !<pools (after being updated) : grid avg. structural root nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,iccp2),intent(in)    :: nlitrmass        !<pools (after being updated) : litter nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: nlitrmass_g      !<pools (after being updated) : grid avg. litter nitrogen mass after being updated in n_processes (\f$g N/m^2\f$)
real,dimension(ilg,iccp2),intent(in)    :: soilnmas         !<pools (after being updated) : soil organic nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: soilnmas_g       !<pools (after being updated) : grid avg. soil organic nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,iccp1),intent(in)    :: nitrifveg        !<nitrification flux from nh4_mass to no3_mass pool (\f$g N/m^2 day\f$) for individual PFTs + bareground
real,dimension(ilg),intent(in)          :: nitrif           !<grid avg. nitrification flux from nh4_mass to no3_mass pool (\f$g N/m^2 day\f$
real,dimension(ilg,iccp1),intent(in)    :: denitveg         !<total denitrification outflux from no3_mass pool (\f$g N/m^2 day\f$) for individual PFTs + bareground
real,dimension(ilg),intent(in)          :: denit            !<grid avg. total denitrification outflux from no3_mass pool (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: bnf_free         !< free-living biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
real,dimension(ilg),intent(in)          :: bnf_free_g       !< grid avg. free-living biological nitrogen fixation, \f$g N m^{-2} day^{-1}\f$
real,dimension(ilg),intent(in)          :: ndeposit         !<nitrogen deposition \f$(g N m^{-2} yr^{-1}\f$
real,dimension(ilg),intent(in)          :: radj             !<latitude in radians
real,dimension(ilg,iccp1),intent(in)    :: nvolveg          !<nitrogen volatilization loss from the nh4_mass (i.e., ammonium) pool (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nvol             !<grid avg. nitrogen volatilization loss from the nh4_mass (i.e., ammonium) pool (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: nleachveg        !<nitrogen leaching from nitrate pool (\f$g N/m^2s\f$)
real,dimension(ilg),intent(in)          :: nleach           !<grid avg. nitrogen leaching from nitrate pool (\f$g N/m^2s\f$)
real,dimension(ilg,iccp1),intent(in)    :: no_nitveg        !<nitrification-induced NO loss for individual 9 pfts + bare (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: no_nit           !<grid avg. nitrification-induced NO loss (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: n2o_nitveg       !<nitrification-induced N2O loss for individual 9 pfts + bare (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: n2o_nit          !<grid avg. nitrification-induced N2O loss (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nuptakeveg_nh4   !<nh4+ uptake for individual PFTs (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nuptake_nh4      !<grid avg. nh4+ uptake (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nuptakeveg_no3   !<no3- uptake for individual PFTs (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nuptake_no3      !<grid avg. no3- uptake (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nleafns2sveg     !<nitrogen flux from non-structural to structural leaf pool (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nstemns2sveg     !<nitrogen flux from non-structural to structural stem pool (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nrootns2sveg     !<nitrogen flux from non-structural to structural root pool (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nntchveg_l       !<net change in ngleafmas_ns after allocation, re_allocation & resorption (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nntchveg_s       !<net change in nstemmass_ns after allocation, re_allocation & resorption (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nntchveg_r       !<net change in nrootmass_ns after allocation, re_allocation & resorption (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nntch_l          !<grid avg. net change in ngleafmas_ns after allocation, re_allocation & resorption (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nntch_s          !<grid avg. net change in nstemmass_ns after allocation, re_allocation & resorption (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nntch_r          !<grid avg. net change in nrootmass_ns after allocation, re_allocation & resorption (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nlitrveg_l       !<leaf N litterfall (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nlitrveg_s       !<stem N litterfall (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: nlitrveg_r       !<root N litterfall (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nlitr_l          !<grid avg. leaf N litterfall (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nlitr_s          !<grid avg. stem N litterfall (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nlitr_r          !<grid avg. root N litterfall (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: gl2bl_grass_nflux      !<N flux from ngleafmas to nbleafmas (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: gl2bl_grass_nflux_gavg !<grid avg. N flux from ngleafmas to nbleafmas (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: nhumtrsveg       !<N humification for individual PFTs + bare (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nhumtrs          !<grid avg. N humification from litter mass (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: nmineralveg_litr !<N mineralization from litter pool (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: nmineralveg_humus!<N mineralization from organic soil pool (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: nimmobilveg_nh4  !<N immobilization from NH4+ pool to soilnmas (\f$g N/m^2 day\f$)
real,dimension(ilg,iccp1),intent(in)    :: nimmobilveg_no3  !<N immobilization from NO3- pool to soilnmas (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nmineral_litr    !<grid avg. N mineralization from litter pool (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nmineral_humus   !<grid avg. N mineralization from organic soil pool (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nimmobil_nh4     !<grid avg. N immobilization from NH4+ pool to soilnmas (\f$g N/m^2 day\f$)
real,dimension(ilg),intent(in)          :: nimmobil_no3     !<grid avg. N immobilization from NO3- pool to soilnmas (\f$g N/m^2 day\f$)
real,dimension(ilg,icc),intent(in)      :: ntchlveg         !<net change in gleaf biomass after auto. resp. & allocation
real,dimension(ilg,iccp1),intent(in) :: pnh4_mass        !<pools (before being updated): nitrate mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg,iccp1),intent(in) :: pno3_mass        !<pools (before being updated): nitrate mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pngleafmas       !<pools (before being updated): green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pngleafmas_ns    !<pools (before being updated): non-structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pngleafmas_s     !<pools (before being updated): structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pnbleafmas       !<pools (before being updated): brown leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pnstemmass       !<pools (before being updated): stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pnstemmass_ns    !<pools (before being updated): non-structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pnstemmass_s     !<pools (before being updated): structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pnrootmass       !<pools (before being updated): root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pnrootmass_ns    !<pools (before being updated): non-structural root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,icc),intent(in)   :: pnrootmass_s     !<pools (before being updated): structural root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
real,dimension(ilg,iccp1),intent(in) :: pnlitrmass       !<pools (before being updated): litter nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg,iccp1),intent(in) :: psoilnmas        !<pools (before being updated): soil organic nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnh4_mass_g      !<pools (before being updated): grid avg. NH4 mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pno3_mass_g      !<pools (before being updated): grid avg. NO3- pool (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pngleafmas_g     !<pools (before being updated): grid avg. green leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pngleafmas_ns_g  !<pools (before being updated): grid avg. non-structural green leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pngleafmas_s_g   !<pools (before being updated): grid avg. structural green leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnbleafmas_g     !<pools (before being updated): grid avg. brown leaf nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnstemmass_g     !<pools (before being updated): grid avg. stem nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnstemmass_ns_g  !<pools (before being updated): grid avg. non-structural stem nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnstemmass_s_g   !<pools (before being updated): grid avg. structural stem nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnrootmass_g     !<pools (before being updated): grid avg. root nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnrootmass_ns_g  !<pools (before being updated): grid avg. non-structural root nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnrootmass_s_g   !<pools (before being updated): grid avg. structural root nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: pnlitrmass_g     !<pools (before being updated): grid avg. litter nitrogen mass after being updated in n_processes (\f$g N/m^2\f$)
real,dimension(ilg),intent(in)          :: psoilnmas_g      !<pools (before being updated): grid avg. soil organic nitrogen mass (\f$g N/m^2\f$)
real,dimension(ilg,iccp1),intent(in)    :: appl_fert        !<applied nitrogen fertilizer \f$(g N /m^2 cropland day\f$
real,dimension(ilg),intent(in)          :: appl_fert_g      !<grid avg. applied nitrogen fertilizer \f$(g N /m^2 day\f$
integer,dimension(ilg,icc),intent(in)              :: lfstatus       !<integer indicating leaf status or mode


! Local
character(8) :: pftkind
integer i, j, k
real :: diff1, diff2, diff3
real :: doy
real,dimension(ilg,icc) :: temp_ns        !<temporary variable to store N netchange associated with non-structural pools
real,dimension(ilg,icc) :: temp_s         !<temporary variable to store N netchange associated with structural pools
real,dimension(ilg,icc) :: temp_litr_ns   !<temporary variable to store N litter associated with non-structural pools
real,dimension(ilg,icc) :: temp_litr_s    !<temporary variable to store N litter associated with structural pools
real :: pngleafmas_ns_temp, pngleafmas_s_temp, pngleafmas_temp, pnstemmass_ns_temp, pnstemmass_s_temp, pnstemmass_temp, pnrootmass_ns_temp, pnrootmass_s_temp, pnrootmass_temp

temp_ns(:,:)= 0.0
temp_s(:,:) = 0.0
temp_litr_ns(:,:)= 0.0
temp_litr_s(:,:) = 0.0

if (leapnow) then
  doy = 366.
else
  doy = 365.
end if 

! To check N budget go through each N pool for each PFT:

!!----------NH4+ (i.e., ammonium) pool
do 100 j = 1, iccp1
   do 110 i = il1, il2
        diff1 = 0.0
        diff1 = abs( nh4_mass(i,j) - (  pnh4_mass(i,j)         &
                                      + appl_fert(i,j)         &
                                      + (0.5*ndeposit(i)/doy) &
                                      - nvolveg(i,j)           &
                                      - nitrifveg(i,j)         &
                                      - no_nitveg(i,j)         &
                                      - n2o_nitveg(i,j)        &
                                      + nmineralveg_litr(i,j)  &
                                      + nmineralveg_humus(i,j) &
                                      - nimmobilveg_nh4(i,j)))
        if (j.le.icc) then
           diff1 = abs(diff1 - nuptakeveg_nh4(i,j))
        endif

        if(diff1.gt.ntolerance)then
          write(*,*)''
          write(*,*)'nh4_mass(',i,',',j,')         = ',nh4_mass(i,j)
          write(*,*)'pnh4_mass(',i,',',j,')        = ',pnh4_mass(i,j)
          write(*,*)'ndeposit(',i,',',j,')         = ',0.5*ndeposit(i)/doy
          write(*,*)'appl_fert(',i,',',j,')        = ',appl_fert(i,j)
          write(*,*)'nvolveg(',i,',',j,')          = ',nvolveg(i,j)
          write(*,*)'nitrifveg(',i,',',j,')        = ',nitrifveg(i,j)
          write(*,*)'no_nitveg(',i,',',j,')        = ',no_nitveg(i,j)
          write(*,*)'n2o_nitveg(',i,',',j,')       = ',n2o_nitveg(i,j)
          write(*,*)'nmineralveg_litr(',i,',',j,') = ',nmineralveg_litr(i,j)
          write(*,*)'nmineralveg_humus(',i,',',j,')= ',nmineralveg_humus(i,j)
          write(*,*)'nimmobilveg_nh4(',i,',',j,')  = ',nimmobilveg_nh4(i,j)
          if (j.le.icc) then
             write(*,*)'nuptakeveg_nh4(',i,',',j,')= ',nuptakeveg_nh4(i,j)
          endif
          write(*,*)''
          write(*,*)'imbalance in Ammonium pool at (i)=',i,'pft=',j
          write(*,*)''
          call errorHandler('balnit', - 1)
        endif

110 continue
100 continue


!!---------- NO3- (i.e., nitrate) pool

do 120 j = 1, iccp1
    do 130 i = il1, il2
       diff1 = 0.0
       diff1 = abs(no3_mass(i,j)- (  pno3_mass(i,j)         &
                                   + (0.5*ndeposit(i)/doy) &
                                   + nitrifveg(i,j)         &
                                   - denitveg(i,j)          &
                                   - nleachveg(i,j)         &
                                   - nimmobilveg_no3(i,j)))
       if (j.le.icc) then
          diff1 = abs(diff1 - nuptakeveg_no3(i,j))
       endif

       if(diff1.gt.ntolerance) then
          write(*,*)''
          write(*,*)'no3_mass(',i,',',j,')       = ',no3_mass(i,j)
          write(*,*)'pno3_mass(',i,',',j,')      = ',pno3_mass(i,j)
          write(*,*)'ndeposit(',i,',',j,')       = ',0.5*ndeposit(i)/doy
          write(*,*)'nitrifveg(',i,',',j,')      = ',nitrifveg(i,j)
          write(*,*)'denitveg(',i,',',j,')       = ',denitveg(i,j)
          write(*,*)'nleachveg(',i,',',j,')      = ',nleachveg(i,j)
          write(*,*)'nimmobilveg_no3(',i,',',j,')= ',nimmobilveg_no3(i,j)
          if (j.le.icc) then
           write(*,*)'nuptakeveg_no3(',i,',',j,')= ',nuptakeveg_no3(i,j)
          endif
          write(*,*)''
          write(*,*)'imbalance in Nitrate pool at (i)=',i,'pft=',j  !Note: j=10 is bareground
          write(*,*)''
          call errorHandler('balnit', - 2)
       endif

130 continue
120 continue


!!---------- ngleafmass pool

do 140 j = 1, icc
   pftkind = ctempfts(j)
   do 150 i = il1, il2
        diff3 = 0.0
        if(pftkind == 'GrassC3 ' .or. pftkind == 'GrassC4 ' .or. pftkind == 'Sedge   ') then
          diff3 = ngleafmas(i,j)   - ( pngleafmas(i,j)    + nntchveg_l(i,j) - gl2bl_grass_nflux(i,j))
        else
          diff3 = ngleafmas(i,j)   - ( pngleafmas(i,j)    + nntchveg_l(i,j) - nlitrveg_l(i,j))
        end if
        if (abs(diff3) .gt. ntolerance) then
            write(*,*)''
            write(*,*)'ngleafmas   (',i,',',j,')     = ',ngleafmas(i,j)
            write(*,*)'pngleafmas  (',i,',',j,')     = ',pngleafmas(i,j)
            write(*,*)'nntchveg_l  (',i,',',j,')     = ',nntchveg_l(i,j)
            write(*,*)'nleafns2sveg(',i,',',j,')     = ',nleafns2sveg(i,j)
            write(*,*)'nlitrveg_l  (',i,',',j,')     = ',nlitrveg_l(i,j)
            write(*,*)'gl2bl_grass_nflux(',i,',',j,')= ',gl2bl_grass_nflux(i,j)
            write(*,*)'lfstatus       (',i,',',j,')= ',lfstatus(i,j)
            write(*,*)''
            write(*,*)'imbalance in leaf N pool at (i)=',i,'pft=',j
            write(*,*)''
            call errorHandler('balnit', - 5)
        endif
150 continue
140 continue

!!---------- nbleafmass pool

do 160 j = 1, icc
   pftkind = ctempfts(j)
   do 170 i = il1, il2
        diff1 = 0.0
        if (pftkind == 'GrassC3 ' .or. pftkind == 'GrassC4 ' .or. pftkind == 'Sedge   ') then
              diff1 = abs(nbleafmas(i,j) - (  pnbleafmas(i,j)        &
                                            + gl2bl_grass_nflux(i,j) &
                                            - nlitrveg_l(i,j)))
        endif

        if (diff1 .gt. ntolerance) then
           write(*,*)''
           write(*,*)'nbleafmas(',i,',',j,')        = ',nbleafmas(i,j)
           write(*,*)'pnbleafmas(',i,',',j,')       = ',pnbleafmas(i,j)
           write(*,*)'gl2bl_grass_nflux(',i,',',j,')= ',gl2bl_grass_nflux(i,j)
           write(*,*)'nlitrveg_l(',i,',',j,')       = ',nlitrveg_l(i,j)
           write(*,*)''
           write(*,*)'imbalance in brown leaf N pool at (i)=',i,'pft=',j
           write(*,*)''
           call errorHandler('balnit', - 6)
        endif

170 continue
160 continue


!!---------- nstemmass pool

do 180 j = 1, icc
   do 190 i = il1, il2
          diff3 = 0.0
          diff3 = nstemmass(i,j)   - ( pnstemmass(i,j)    + nntchveg_s(i,j) - nlitrveg_s(i,j))
          if (diff3 .gt. ntolerance) then
             write(*,*)''
             write(*,*)'nstemmass     (',i,',',j,')= ',nstemmass(i,j)
             write(*,*)'pnstemmass    (',i,',',j,')= ',pnstemmass(i,j)
             write(*,*)'nntchveg_s    (',i,',',j,')= ',nntchveg_s(i,j)
             write(*,*)'nstemns2sveg(',i,',',j,')= ',nstemns2sveg(i,j)
             write(*,*)'nlitrveg_s    (',i,',',j,')= ',nlitrveg_s(i,j)
             write(*,*)''
             write(*,*)'imbalance in stem N pool at (i)=',i,'pft=',j
             write(*,*)''
             call errorHandler('balnit', - 9)
          endif
190 continue
180 continue

!!---------- nrootmass

do 200 j = 1, icc
   do 210 i = il1, il2
          diff3 = 0.0
          diff3 = nrootmass(i,j)   - ( pnrootmass(i,j)    + nntchveg_r(i,j) - nlitrveg_r(i,j))
          if (diff3.gt.ntolerance) then
             write(*,*)''
             write(*,*)'nrootmass     (',i,',',j,')= ',nrootmass(i,j)
             write(*,*)'pnrootmass    (',i,',',j,')= ',pnrootmass(i,j)
             write(*,*)'nntchveg_r    (',i,',',j,')= ',nntchveg_r(i,j)
             write(*,*)'nrootns2sveg(',i,',',j,')= ',nrootns2sveg(i,j)
             write(*,*)'nlitrveg_r    (',i,',',j,')= ',nlitrveg_r(i,j)
             write(*,*)''
             write(*,*)'imbalance in root N pool at (i)=',i,'pft=',j
             write(*,*)''
             call errorHandler('balnit', - 12)
          endif

210 continue
200 continue

!!---------- nlitrmass

do 220 j = 1, iccp1
   do 230 i = il1, il2
          diff1 = 0.0
          if (j.le.icc) then
             diff1 = abs(nlitrmass(i,j) - (  pnlitrmass(i,j) &
                                           + nlitrveg_l(i,j) &
                                           + nlitrveg_s(i,j) &
                                           + nlitrveg_r(i,j) &
                                           - nhumtrsveg(i,j) &
                                           - nmineralveg_litr(i,j)))
          else
             diff1 = abs(nlitrmass(i,j) - (  pnlitrmass(i,j) &
                                           - nhumtrsveg(i,j) &
                                           - nmineralveg_litr(i,j)))
          endif
          if(abs(diff1) .gt. ntolerance)then
           write(*,*)''
           write(*,*)'nlitrmass (',i,',',j,')= ',nlitrmass(i,j)
           write(*,*)'pnlitrmass(',i,',',j,')= ',pnlitrmass(i,j)
           if (j.le.icc) then
           write(*,*)'nlitrveg_l(',i,',',j,')= ',nlitrveg_l(i,j)
           write(*,*)'nlitrveg_s(',i,',',j,')= ',nlitrveg_s(i,j)
           write(*,*)'nlitrveg_r(',i,',',j,')= ',nlitrveg_r(i,j)
           endif
           write(*,*)'nhumtrsveg(',i,',',j,')= ',nhumtrsveg(i,j)
           write(*,*)'nmineralveg_litr(',i,',',j,')= ',nmineralveg_litr(i,j)
           write(*,*)''
           write(*,*)'imbalance in litter N pool at (i)=',i,'pft=',j
           write(*,*)''
           call errorHandler('balnit', - 13)
        endif

230 continue
220 continue

!!---------- soilnmas

do 240 j = 1, iccp1
   do 250 i = il1, il2
          diff1 = 0.0
          diff1 = abs(soilnmas(i,j) - (  psoilnmas(i,j)         &
                                       + bnf_free(i,j)          &
                                       + nhumtrsveg(i,j)        &
                                       - nmineralveg_humus(i,j) &
                                       + nimmobilveg_nh4(i,j)   &
                                       + nimmobilveg_no3(i,j)))
          if(diff1 .gt. ntolerance)then
           write(*,*)''
           write(*,*)'soilnmas (',i,',',j,')        = ',soilnmas(i,j)
           write(*,*)'psoilnmas(',i,',',j,')        = ',psoilnmas(i,j)
           write(*,*)'bnf_free(',i,',',j,')        = ',bnf_free(i,j)
           write(*,*)'nhumtrsveg(',i,',',j,')       = ',nhumtrsveg(i,j)
           write(*,*)'nmineralveg_humus(',i,',',j,')= ',nmineralveg_humus(i,j)
           write(*,*)'nimmobilveg_nh4(',i,',',j,')  = ',nimmobilveg_nh4(i,j)
           write(*,*)'nimmobilveg_no3(',i,',',j,')  = ',nimmobilveg_no3(i,j)
           write(*,*)''
           write(*,*)'imbalance in soil N pool at (i)=',i,'pft=',j
           write(*,*)''
           call errorHandler('balnit', - 14)
        endif

250 continue
240 continue

!###############################################################################
!---------- grid averaged pools

do 260 i = il1, il2

        !--- NH4+ (i.e., ammonium) pool
        diff1 = 0.0
        diff1 = abs(nh4_mass_g(i) - (  pnh4_mass_g(i)&
                                     + appl_fert_g(i)         &
                                     + (0.5*ndeposit(i)/doy) &
                                     - nvol(i)                &
                                     - nitrif(i)              &
                                     - no_nit(i)              &
                                     - n2o_nit(i)             &
                                     - nuptake_nh4(i)         &
                                     + nmineral_litr(i)       &
                                     + nmineral_humus(i)      &
                                     - nimmobil_nh4(i)))

        if(diff1 .gt. ntolerance)then
           write(*,*)''
           write(*,*)'nh4_mass_g (',i,')   = ',nh4_mass_g(i)
           write(*,*)'pnh4_mass_g(',i,')   = ',pnh4_mass_g(i)
           write(*,*)'appl_fert_g(',i,')   = ',appl_fert_g(i)
           write(*,*)'ndeposit   (',i,')   = ',0.5*ndeposit(i)/doy
           write(*,*)'nvol       (',i,')   = ',nvol(i)
           write(*,*)'nitrif     (',i,')   = ',nitrif(i)
           write(*,*)'no_nit     (',i,')   = ',no_nit(i)
           write(*,*)'n2o_nit    (',i,')   = ',n2o_nit(i)
           write(*,*)'nuptake    (',i,')   = ',nuptake_nh4(i)
           write(*,*)'nmineral_litr(',i,') = ',nmineralveg_litr(i,j)
           write(*,*)'nmineral_humus(',i,')= ',nmineralveg_humus(i,j)
           write(*,*)'nimmobil_nh4(',i,')  = ',nimmobilveg_nh4(i,j)
           write(*,*)''
           write(*,*)'imbalance in grid averaged ammonium pool at (i)=',i
           write(*,*)''
           call errorHandler('balnit', - 15)
        endif


      	!--- NO3- (i.e., nitrate) pool
        diff1 = 0.0
        diff1 = abs(no3_mass_g(i) - (  pno3_mass_g(i)         &
                                     + (0.5*ndeposit(i)/doy)  &
                                     + nitrif(i)              &
                                     - denit(i)               &
                                     - nleach(i)              &
                                     - nuptake_no3(i)         &
                                     - nimmobil_no3(i)))

        if(diff1 .gt. ntolerance)then
           write(*,*)''
           write(*,*)'no3_mass_g (',i,') = ',no3_mass_g(i)
           write(*,*)'pno3_mass_g(',i,') = ',pno3_mass_g(i)
           write(*,*)'ndeposit   (',i,') = ',0.5*ndeposit(i)/doy
           write(*,*)'nitrif     (',i,') = ',nitrif(i)
           write(*,*)'denit      (',i,') = ',denit(i)
           write(*,*)'nleach     (',i,') = ',nleach(i)
           write(*,*)'nuptake    (',i,') = ',nuptake_no3(i)
           write(*,*)'nimmobil_no3(',i,')= ',nimmobil_no3(i)
           write(*,*)''
           write(*,*)'imbalance in grid averaged nitrate pool at (i)=',i
           write(*,*)''
           call errorHandler('balnit', - 16)
        endif


        !--- ngleafmas & nbleafmas
        !Note that for grasses gl2bl_grass_nflux_gavg(i) is an outfulux from ngleafmas pool and
        !nlitr_l(i) is an outflux for nbleafmas pool. Therefore, in checking grid avg. mass
        !conservation the two pools (ngleafmas & nbleafmas) should be done together.

        diff1 = 0.0
        diff1 = abs(ngleafmas_g(i)+nbleafmas_g(i) - (  pngleafmas_g(i)+pnbleafmas_g(i)  &
                                                     + nntch_l(i)                       &
                                                     - nlitr_l(i)                       &
                                                     - gl2bl_grass_nflux_gavg(i)        &
                                                     + gl2bl_grass_nflux_gavg(i)))

        if(diff1 .gt. ntolerance)then
           write(*,*)''
           write(*,*)'ngleafmas_g (',i,')= ',ngleafmas_g(i)
           write(*,*)'pngleafmas_g(',i,')= ',pngleafmas_g(i)
           write(*,*)'nbleafmas_g(',i,')= ',nbleafmas_g(i)
           write(*,*)'pnbleafmas_g(',i,')= ',pnbleafmas_g(i)
           write(*,*)'nntch_l     (',i,')= ',nntch_l(i)
           write(*,*)'nlitr_l     (',i,')= ',nlitr_l(i)
           write(*,*)'gl2bl_grass_nflux_gavg(',i,')= ',gl2bl_grass_nflux_gavg(i)
           write(*,*)'diff1       (',i,')= ',diff1
           write(*,*)''
           write(*,*)'imbalance in grid averaged leaf N pool at (i)=',i
           write(*,*)''
           call errorHandler('balnit', - 17)
        endif


        !--- nstemmass
        diff1 = 0.0
        diff1 = abs(nstemmass_g(i) - (  pnstemmass_g(i) &
                                      + nntch_s(i)      &
                                      - nlitr_s(i)))

        if(diff1 .gt. ntolerance)then
           write(*,*)''
           write(*,*)'nstemmass_g (',i,')= ',nstemmass_g(i)
           write(*,*)'pnstemmass_g(',i,')= ',pnstemmass_g(i)
           write(*,*)'nntch_s     (',i,')= ',nntch_s(i)
           write(*,*)'nlitr_s     (',i,')= ',nlitr_s(i)
           write(*,*)''
           write(*,*)'imbalance in grid averaged stem N pool at (i)=',i
           write(*,*)''
           call errorHandler('balnit', - 18)
        endif


        !--- nrootmass
        diff1 = 0.0
        diff1 = abs(nrootmass_g(i) - (  pnrootmass_g(i) &
                                      + nntch_r(i)      &
                                      - nlitr_r(i)))

        if(diff1 .gt. ntolerance)then
           write(*,*)''
           write(*,*)'nrootmass_g (',i,')= ',nrootmass_g(i)
           write(*,*)'pnrootmass_g(',i,')= ',pnrootmass_g(i)
           write(*,*)'nntch_r     (',i,')= ',nntch_r(i)
           write(*,*)'nlitr_r     (',i,')= ',nlitr_r(i)
           write(*,*)''
           write(*,*)'imbalance in grid averaged root N pool at (i)=',i
           write(*,*)''
           call errorHandler('balnit', - 19)
        endif


        !--- nlitrmass
        diff1 = 0.0
        diff1 = abs(nlitrmass_g(i) - (  pnlitrmass_g(i) &
                                      + nlitr_l(i)      &
                                      + nlitr_s(i)      &
                                      + nlitr_r(i)      &
                                      - nhumtrs(i)      &
                                      - nmineral_litr(i)))

         if(diff1 .gt. ntolerance)then
           write(*,*)''
           write(*,*)'nlitrmass_g (',i,') = ',nlitrmass_g(i)
           write(*,*)'pnlitrmass_g(',i,') = ',pnlitrmass_g(i)
           write(*,*)'nlitr_l     (',i,') = ',nlitr_l(i)
           write(*,*)'nlitr_s     (',i,') = ',nlitr_s(i)
           write(*,*)'nlitr_r     (',i,') = ',nlitr_r(i)
           write(*,*)'nhumtrs     (',i,') = ',nhumtrs(i)
           write(*,*)'nmineral_litr(',i,')= ',nmineral_litr(i)
           write(*,*)''
           write(*,*)'imbalance in grid averaged litter N pool at (i)=',i
           write(*,*)''
           call errorHandler('balnit', - 20)
         endif


        !--- soilnmas
        diff1 = 0.0
        diff1 = abs(soilnmas_g(i) - (  psoilnmas_g(i)   &
                                     + bnf_free_g(i)    &
                                     + nhumtrs(i)       &
                                     - nmineral_humus(i)&
                                     + nimmobil_nh4(i)  &
                                     + nimmobil_no3(i)))

        if(diff1 .gt. ntolerance)then
          write(*,*)''
          write(*,*)'soilnmas_g (',i,')   = ',soilnmas_g(i)
          write(*,*)'psoilnmas_g(',i,')   = ',psoilnmas_g(i)
          write(*,*)'bnf_free_g (',i,')   = ',bnf_free_g(i)
          write(*,*)'nhumtrs    (',i,')   = ',nhumtrs(i)
          write(*,*)'nmineral_humus(',i,')= ',nmineral_humus(i)
          write(*,*)'nimmobil_nh4(',i,')  = ',nimmobil_nh4(i)
          write(*,*)'nimmobil_no3(',i,')  = ',nimmobil_no3(i)
          write(*,*)''
          write(*,*)'imbalance in grid averaged soil N pool at (i)=',i
          write(*,*)''
          call errorHandler('balnit', - 21)
        endif


260 continue




end
