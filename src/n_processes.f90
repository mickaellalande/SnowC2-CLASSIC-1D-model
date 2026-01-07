!> \file
!! Includes/calls all the processes/subroutines related to the Nitrogen cycle.
!! @author A. Asaadi
!!


subroutine n_processes(il1, il2, spinfast, thliq, thice, thpor, zbotw, tbar, isand, & ! In
                       iday, radj, lfstatus, THFC, THLW, fcancmx, & ! In
                       delzw, sort, ipeatland, ROFB, CFLUX_GA, USTARBS_GA, leapnow, & ! In
                       QFC, rootdpth, soilpH, nfertil, ndeposit, ntchlveg, ntchsveg, & ! In
                       ntchrveg, gleafmas, gleafmas_ns, gleafmas_s, bleafmas, stemmass, & ! In
                       stemmass_ns, stemmass_s, rootmass, rootmass_ns, rootmass_s, & ! In
                       litrmass, soilcmas, leafns2s, stemns2s, rootns2s, & ! In
                       tltrleaf, tltrstem, tltrroot, gl2bl_grass_cflux, bl2ltr_grass_cflux, & ! In
                       re_alloc_sr2l, ailcg, lfthrs, co2conc, & ! In
                       humtrsvg, ltresveg, & ! In
                       scresveg, psisat, bi, & ! In
                       ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! In/Out
                       nstemmass, nstemmass_ns, nstemmass_s, & ! In/Out
                       nrootmass, nrootmass_ns, nrootmass_s, & ! In/Out
                       nlitrmass, soilnmas, & ! In/Out
                       bnf_free, bnf_nat, bnf_ant, bnf_tot, nstress, & !Out 
                       nh4_mass, no3_mass, nvolveg, nleachveg, & ! Out
                       nitrifveg, no_nitveg, no_denitveg, & ! Out
                       no_nitdenitveg, n2o_nitveg, n2o_denitveg, & ! Out
                       n2o_nitdenitveg, n2_denitveg, appl_fert, & ! Out
                       ndep_nh4, ndep_no3, ndemandveg_wp_npp, & ! Out
                       nuptakeveg_p_nh4, nuptakeveg_p_no3, & ! Out
                       nuptakeveg_a_actl_nh4, & ! Out
                       nuptakeveg_a_actl_no3, nuptakeveg, &
                       nallocveg_l, nallocveg_s, & ! Out
                       nallocveg_r, nresorpedveg_s, & ! Out
                       nresorpedveg_r, nre_allocveg_s2l, & ! Out
                       nre_allocveg_r2l, nleafns2sveg, nstemns2sveg, & ! Out
                       nrootns2sveg, nlitrveg_l, nlitrveg_s, & ! Out
                       nlitrveg_r, nlitrveg, gl2bl_grass_nflux, & ! Out
                       c2nveg_l, c2nveg_s, c2nveg_r, & ! Out
                       c2nveg_wp, c2nveg_litr, c2nveg_humus, & ! Out
                       nhumtrsveg, nmineralveg_litr, & ! Out
                       nmineralveg_humus, netnmineralveg, nimmobilveg_nh4, & ! Out
                       nimmobilveg_no3, nvgbiomas_veg, & ! Out
                       fNnetlandveg, redcoeff_vcmax) ! Out

  use classicParams, only : ilg, ignd, ican, icc, iccp1, ctempfts, pi, fertstd, deltat, iccp2, &
                            fertend, zero, c2n_lmax, c2n_smax, c2n_rmax, min_ns2t_l, min_ns2t_s, min_ns2t_r, &
                            gleafmasmax, stemmassmax, rootmassmax, initial_c2n, ntolerance, k_redcoeff_vcmax

  implicit none

  integer,intent(in) :: il1       !< il1=1
  integer,intent(in) :: il2       !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in) :: iday      !< day of year
  integer,intent(in) :: spinfast  !< spinup factor for soil carbon & nitrogen pools
  logical,intent(in) :: leapnow   !< true if this year is a leap year
  integer,intent(in), dimension(icc) :: sort           !< index for correspondence between PFTs and the 12 values in parameters vectors
  integer,intent(in), dimension(ilg) :: ipeatland      !< Peatland flag: 0 = not a peatland, 1= bog, 2 = fen
  integer,intent(in), dimension(ilg,ignd) :: isand     !< flag for bedrock or ice in a soil layer
  integer,intent(in), dimension(ilg,icc) :: lfstatus   !< integer indicating leaf status or mode
  !! 1- max. growth or onset, when all npp is allocated to leaf
  !! 2- normal growth, when npp is allocated to leaf, stem, and root
  !! 3- offset, when zero npp is allocated to leaf 
  !! 4- no leaf
  real,intent(in), dimension(ilg) :: radj               !< latitude in radians
  real,intent(in), dimension(ilg,ignd) :: psisat        !< soil matric potential at saturation, Mpa
  real,intent(in), dimension(ilg,ignd) :: bi            !< parameter B of Clapp and Hornberger 1978
  real,intent(in), dimension(ilg,ignd) :: tbar          !< soil temperature, K 
  real,intent(in), dimension(ilg,ignd) :: thliq         !< volumetric soil moisture, \f$m^3 m^{-3}\f$
  real,intent(in), dimension(ilg,ignd) :: thice         !< frozen volumetric soil moisture, \f$m^3 m^{-3}\f$
  real,intent(in), dimension(ilg,ignd) :: thpor         !< soil total porosity, \f$m^3 m^{-3}\f$
  real,intent(in), dimension(ilg,ignd) :: zbotw         !< bottom of each soil layer, \f$m\f$
  real,intent(in), dimension(ilg,icc) :: fcancmx        !< max. fractional coverage of each PFT
  real,intent(in), dimension(ilg,ignd) :: THFC          !< volumetric soil moisture at field capacity, \f$m^3 m^{-3}\f$
  real,intent(in), dimension(ilg,ignd) :: THLW          !< volumetric soil moisture at wilting point, \f$m^3 m^{-3}\f$
  real,intent(in), dimension(ilg) :: soilpH             !< soil pH
  real,intent(in), dimension(ilg) :: nfertil            !< nitrogen fertilizer, \f$g N m^{-2} yr^{-1}\f$
  real,intent(in), dimension(ilg) :: ndeposit           !< nitrogen deposition, \f$g N m^{-2} yr^{-1}\f$
  real,intent(in), dimension(ilg) :: CFLUX_GA           !< aerodynamic conductance, inverse of aerodynamic resistance, \f$m s^{-1}\f$
  real,intent(in), dimension(ilg) :: USTARBS_GA         !< friction velocity, \f$m s^{-1}\f$
  real,intent(in), dimension(ilg) :: ROFB               !< base flow from bottom of soil column, \f$kg H2O m^{-2} s^{-1}\f$
  real,intent(in), dimension(ilg,ignd) :: QFC           !< water removed from each soil layer by transpiration, \f$kg H2O m^{-2} s^{-1}\f$
  real,intent(in), dimension(ilg,ignd) :: delzw         !< thickness of each soil layer, \f$m\f$
  real,intent(in), dimension(ilg,icc) :: ntchlveg       !< C allocation to leaf, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: ntchsveg       !< C allocation to stem, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: ntchrveg       !< C allocation to root, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: rootdpth       !< rooting depth (contains 99% of rootmass), \f$m\f$
  real,intent(in), dimension(ilg,icc) :: gleafmas       !< green leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: gleafmas_ns    !< non-structural green leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: gleafmas_s     !< structural green leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: bleafmas       !< brown leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: stemmass       !< stem carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: stemmass_ns    !< non-structural stem carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: stemmass_s     !< structural stem carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: rootmass       !< root carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: rootmass_ns    !< non-structural root carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: rootmass_s     !< structural root carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,iccp2,ignd) :: litrmass !< detritus carbon for each PFT + bare soil + LUC product, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,iccp2,ignd) :: soilcmas !< SOM carbon for each PFT + bare soil + LUC product, \f$kg C m^{-2}\f$
  real,intent(in), dimension(ilg,icc) :: tltrleaf       !< total leaf carbon litter fall rate, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: tltrstem       !< total stem carbon litter fall rate, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: tltrroot       !< total root carbon litter fall rate, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: re_alloc_sr2l  !< carbon from stem and root allocated to leaf during leaf onset for each PFT, \f$g C m^{-2} day^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: ailcg          !< green lai
  real,intent(in), dimension(ilg,icc) :: lfthrs         !< threshold lai to determine leaf status
  real,intent(in), dimension(ilg,icc) :: leafns2s       !< carbon flux from non-structural to structural leaf carbon, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: stemns2s       !< carbon flux from non-structural to structural stem carbon, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: rootns2s       !< carbon flux from non-structural to structural root carbon, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,icc) :: gl2bl_grass_cflux  !< carbon flux from green to brown leaf carbon, \f$kg C m^{-2} day^{-1}\f$ 
  real,intent(in), dimension(ilg,icc) :: bl2ltr_grass_cflux !< carbon flux from brown leaf carbon to detritus carbon, \f$kg C m^{-2} day^{-1}\f$
  real,intent(in), dimension(ilg,iccp2,ignd) :: humtrsvg     !< carbon flux from detritus to SOM for individual PFTs + bare soil + LUC product, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,iccp2,ignd) :: ltresveg     !< detritus carbon respiration rate + LUC product, \f$kg C m^{-2} day^{-1}\f$
  real,intent(in), dimension(ilg,iccp2,ignd) :: scresveg     !< SOM carbon respiration rate + LUC product, \f$kg C m^{-2} day^{-1}\f$
  real,intent(in), dimension(ilg) :: co2conc            !< atmospheric CO2 concentration, ppm

  real,intent(out), dimension(ilg,iccp1) :: bnf_tot        !< total biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: bnf_free       !< free-living biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: bnf_nat          !< natural symbiotic biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: bnf_ant          !< anthropogenic symbiotic biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nstress          !< plant nitrogen stress
  real,intent(out), dimension(ilg,iccp1) :: nitrifveg      !< nitrification flux from NH4+ to NO3- for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1)  :: no_nitveg     !< nitrification NO loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: no_denitveg    !< denitrification NO loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: no_nitdenitveg !< total NO loss from denitrification and nitrification for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: n2o_nitveg     !< nitrification N2O loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: n2o_denitveg   !< denitrification N2O loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: n2o_nitdenitveg!< total N2O loss from denitrification and nitrification for individual PFTs + bare, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: n2_denitveg    !< N2 loss from denitrification for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: nvolveg        !< NH3 volatilization loss from NH4+, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: nleachveg      !< nitrogen leaching from NO3-, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: appl_fert      !< daily nitrogen fertilizer for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: ndep_nh4       !< daily nitrogen deposition into NH4+ for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: ndep_no3       !< daily nitrogen deposition into NO3- for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: ndemandveg_wp_npp!< whole plant nitrogen demand based on NPP for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nuptakeveg_p_nh4 !< passive nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nuptakeveg_p_no3 !< passive no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nuptakeveg_a_actl_nh4 !< actual active nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nuptakeveg_a_actl_no3 !< actual active no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nuptakeveg       !< total N uptake (active+passive, NH4+NO3) for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nleafns2sveg     !< nitrogen flux from non-structural to structural leaf nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nstemns2sveg     !< nitrogen flux from non-structural to structural stem nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nrootns2sveg     !< nitrogen flux from non-structural to structural root nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nallocveg_l      !< nitrogen allocation to leaf for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nallocveg_s      !< nitrogen allocation to stem for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nallocveg_r      !< nitrogen allocation to root for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nresorpedveg_s   !< resorbed nitrogen from leaf allocated to stem, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nresorpedveg_r   !< resorbed nitrogen from leaf allocated to root, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nre_allocveg_s2l !< non-structural nitrogen reallocation from stem to leaf during leaf out, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nre_allocveg_r2l !< non-structural nitrogen reallocation from root to leaf during leaf out, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nlitrveg_l       !< total leaf nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nlitrveg_s       !< total stem nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nlitrveg_r       !< total root nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nlitrveg         !< total nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: gl2bl_grass_nflux!< nitrogen flux from green to brown leaf nitrogen for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: c2nveg_l         !< simulated C:N ratio for leaf, \f$g C g N^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: c2nveg_s         !< simulated C:N ratio for stem, \f$g C g N^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: c2nveg_r         !< simulated C:N ratio for root, \f$g C g N^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: c2nveg_wp        !< simulated C:N ratio for the whole plant, \f$g C g N^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: c2nveg_litr    !< simulated C:N ratio for detritus, \f$g C g N^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: c2nveg_humus   !< simulated C:N ratio for SOM, \f$g C g N^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: nhumtrsveg     !< nitrogen flux from detritus to SOM for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: nmineralveg_litr  !< nitrogen mineralization from detritus, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: nmineralveg_humus !< nitrogen mineralization from SOM, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: nimmobilveg_nh4   !< nitrogen immobilization from NH4+ to soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: nimmobilveg_no3   !< nitrogen immobilization from NO3- to soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,iccp1) :: netnmineralveg    !< nitrogen mineralization, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: nvgbiomas_veg    !< vegetation nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(out), dimension(ilg,iccp1) :: fNnetlandveg   !< net terrestrial nitrogen flux, \f$g N m^{-2} day^{-1}\f$
  real,intent(out), dimension(ilg,icc) :: redcoeff_vcmax   !< reduction coefficient of V_c,max passed to the photosynthesis subroutine

  real,intent(inout), dimension(ilg,iccp1) :: nh4_mass     !< NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,iccp1) :: no3_mass     !< NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: ngleafmas      !< green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: ngleafmas_ns   !< non-structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: ngleafmas_s    !< structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: nbleafmas      !< brown leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: nstemmass      !< stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: nstemmass_ns   !< non-structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: nstemmass_s    !< structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: nrootmass      !< root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: nrootmass_ns   !< non-structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,icc) :: nrootmass_s    !< structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,iccp2) :: nlitrmass    !< detritus nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg,iccp2) :: soilnmas     !< SOM nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$

  ! local variables

  integer :: i, j
  character(8) :: pftkind
  real :: temp, temp1, temp2, temp3
  real :: ns2t_nleaf_ratio  !< ratio of non-structural to total leaf nitrogen for each PFT
  real :: ns2t_nstem_ratio  !< ratio of non-structural to total stem nitrogen for each PFT
  real :: ns2t_nroot_ratio  !< ratio of non-structural to total root nitrogen for each PFT
  real,dimension(ilg,iccp2) :: soilcmas_bulk    !< SOM carbon for each PFT + bare soil + LUC product, \f$kg C m^{-2}\f$
  real,dimension(ilg,iccp2) :: litrmass_bulk    !< detritus carbon for each PFT + bare soil + LUC product, \f$kg C m^{-2}\f$
  real,dimension(ilg,iccp2) :: humtrsvg_bulk    !< carbon flux from detritus to SOM for individual PFTs + bare soil + LUC product, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,dimension(ilg,iccp2) :: ltresveg_bulk    !< detritus carbon respiration rate + LUC product, \f$kg C m^{-2} day^{-1}\f$
  real,dimension(ilg,iccp2) :: scresveg_bulk    !< SOM carbon respiration rate + LUC product, \f$kg C m^{-2} day^{-1}\f$
  real,dimension(ilg) :: litrmass_bulk_g        !< grid avg. detritus carbon, \f$kg C m^{-2}\f$
  real,dimension(ilg) :: soilcmas_bulk_g        !< grid avg. SOM carbon, \f$kg C m^{-2}\f$
  real,dimension(ilg,iccp1) :: pno3_mass        !< previous NO3-, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: pnh4_mass        !< previous NH4+, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pngleafmas         !< previous green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pngleafmas_ns      !< previous non-structural green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pngleafmas_s       !< previous structural green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pnbleafmas         !< previous brown leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pnstemmass         !< previous stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pnstemmass_ns      !< previous non-structural stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pnstemmass_s       !< previous structural stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pnrootmass         !< previous root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pnrootmass_ns      !< previous non-structural root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: pnrootmass_s       !< previous structural root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: pnlitrmass       !< previous detritus nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: psoilnmas        !< previous SOM nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnh4_mass_g            !< previous grid avg. NH4+, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pno3_mass_g            !< previous grid avg. NO3-, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pngleafmas_g           !< previous grid avg. green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pngleafmas_ns_g        !< previous grid avg. non-structural green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pngleafmas_s_g         !< previous grid avg. structural green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnbleafmas_g           !< previous grid avg. brown leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnstemmass_g           !< previous grid avg. stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnstemmass_ns_g        !< previous grid avg. non-structural stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnstemmass_s_g         !< previous grid avg. structural stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnrootmass_g           !< previous grid avg. root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnrootmass_ns_g        !< previous grid avg. non-structural root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnrootmass_s_g         !< previous grid avg. structural root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: pnlitrmass_g           !< previous grid avg. detritus nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: psoilnmas_g            !< previous grid avg. SOM nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: barefrac               !< bare fraction of each grid cell
  real,dimension(ilg) :: solctempavg      !< average soil temperature over soil layers in which a process occurs, K
  real,dimension(ilg,iccp1) :: denitveg         !< total denitrification outflux from NO3- for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: denit                  !< grid avg. total denitrification outflux from NO3-, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nuptakeveg_nh4     !< total nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_nh4            !< grid avg. total nh4+ uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nuptakeveg_no3     !< total no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_no3            !< grid avg. total no3- uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: ndemandveg_l_npp   !< leaf nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: ndemandveg_s_npp   !< stem nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: ndemandveg_r_npp   !< root nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nntchveg_l         !< nitrogen net change in ngleafmas_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nntchveg_s         !< nitrogen net change in nstemmass_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nntchveg_r         !< nitrogen net change in nrootmass_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nntch_l                !< grid avg. nitrogen net change in ngleafmas_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nntch_s                !< grid avg. nitrogen net change in nstemmass_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nntch_r                !< grid avg. nitrogen net change in nrootmass_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,iccp1) :: nh4_mass_in      !< previous NH4+ (for N fertilizer & N deposition), \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: no3_mass_in      !< previous NO3- (for N fertilizer & N deposition), \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: ltrestep         !< detritus carbon respiration rate for individual PFTs + bare soil, \f$kg C m^{-2} day^{-1}\f$
  real,dimension(ilg,iccp1) :: screstep         !< SOM carbon respiration rate for individual PFTs + bare soil, \f$kg C m^{-2} day^{-1}\f$
  real,dimension(ilg,icc)   :: exsc2n_l         !< leaf C:N ratio excess above maximum, \f$g C g N^{-1}\f$
  real,dimension(ilg,icc)   :: exsc2n_s         !< stem C:N ratio excess above maximum, \f$g C g N^{-1}\f$
  real,dimension(ilg,icc)   :: exsc2n_r         !< root C:N ratio excess above maximum, \f$g C g N^{-1}\f$
  real,dimension(ilg,icc)   :: nallocf_maxc2n_l !< weight of exsc2n_l to calculate whole plant weighted excess above maximum C:N ratio
  real,dimension(ilg,icc)   :: nallocf_maxc2n_s !< weight of exsc2n_s to calculate whole plant weighted excess above maximum C:N ratio
  real,dimension(ilg,icc)   :: nallocf_maxc2n_r !< weight of exsc2n_r to calculate whole plant weighted excess above maximum C:N ratio
  real,dimension(ilg,icc)   :: wexsc2n_wp       !< whole plant weighted excess above maximum C:N ratio, \f$g C g N^{-1}\f$
  real,dimension(ilg) :: ngleafmas_ns_g         !< grid avg. non-structural green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: ngleafmas_s_g          !< grid avg. structural green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nstemmass_ns_g         !< grid avg. non-structural stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nstemmass_s_g          !< grid avg. structural stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nrootmass_ns_g         !< grid avg. non-structural root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nrootmass_s_g          !< grid avg. structural root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: gleafmas_g             !< grid avg. green leaf carbon, \f$kg C m^{-2}\f$
  real,dimension(ilg) :: bleafmas_g             !< grid avg. brown leaf carbon, \f$kg C m^{-2}\f$
  real,dimension(ilg) :: stemmass_g             !< grid avg. stem carbon, \f$kg C m^{-2}\f$
  real,dimension(ilg) :: rootmass_g             !< grid avg. root carbon, \f$kg C m^{-2}\f$
  real,dimension(ilg) :: bnf_free_g            !< grid avg. free-living biological nitrogen fixation, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nitrif               !< grid avg. nitrification flux, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: no_nit               !< grid avg. nitrification NO loss, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: no_denit             !< grid avg. denitrification NO loss, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: no_nitdenit          !< grid avg. total NO loss from denitrification and nitrification, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: n2o_nit              !< grid avg. nitrification N2O loss, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: n2o_denit            !< grid avg. denitrification N2O loss, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: n2o_nitdenit         !< grid avg. total N2O loss from denitrification and nitrification, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: n2_denit             !< grid avg. N2 loss from denitrification, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nvol                 !< grid avg. nitrogen volatilization loss from NH4+, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nleach               !< grid avg. nitrogen leaching from NO3-, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: ndemand_wp_npp       !< grid avg. whole plant nitrogen demand based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_p_nh4        !< grid avg. passive nh4+ uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_p_no3        !< grid avg. passive no3- uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_a_actl_nh4   !< grid avg. actual active nh4+ uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_a_actl_no3   !< grid avg. actual active no3- uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nleafns2s            !< grid avg. nitrogen flux from non-structural to structural leaf nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nstemns2s            !< grid avg. nitrogen flux from non-structural to structural stem nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nrootns2s            !< grid avg. nitrogen flux from non-structural to structural root nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nalloc_l             !< grid avg. nitrogen allocation to leaf, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nalloc_s             !< grid avg. nitrogen allocation to stem, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nalloc_r             !< grid avg. nitrogen allocation to root, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nresorped_s          !< grid avg. resorbed nitrogen from leaf allocated to stem, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nresorped_r          !< grid avg. resorbed nitrogen from leaf allocated to root, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nre_alloc_s2l        !< grid avg. non-structural nitrogen reallocation from stem to leaf during leaf out, \f$g N m^{-2} day^{-1}\     f$
  real,dimension(ilg) :: nre_alloc_r2l        !< grid avg. non-structural nitrogen reallocation from root to leaf during leaf out, \f$g N m^{-2} day^{-1}\     f$
  real,dimension(ilg) :: nlitr_l              !< grid avg. total leaf nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nlitr_s              !< grid avg. total stem nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nlitr_r              !< grid avg. total root nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: gl2bl_grass_nflux_gavg!< grid avg. nitrogen flux from green to brown leaf nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nhumtrs              !< grid avg. nitrogen flux from detritus to SOM, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nmineral_litr        !< grid avg. nitrogen mineralization from detritus, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nmineral_humus       !< grid avg. nitrogen mineralization from SOM, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nimmobil_nh4         !< grid avg. nitrogen immobilization from NH4+ to soil, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nimmobil_no3         !< grid avg. nitrogen immobilization from NO3- to soil, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: fNnetland            !< grid avg. net terrestrial nitrogen flux, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nvgbiomas          !< grid avg. vegetation nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nh4_mass_g         !< grid avg. NH4+, \f$g N m^{-2}\f$
  real,dimension(ilg) :: no3_mass_g         !< grid avg. NO3+, \f$g N m^{-2}\f$
  real,dimension(ilg) :: ngleafmas_g        !< grid avg. green leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nbleafmas_g        !< grid avg. brown leaf nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nstemmass_g        !< grid avg. stem nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nrootmass_g        !< grid avg. root nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: nlitrmass_g        !< grid avg. detritus nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: soilnmas_g         !< grid avg. SOM nitrogen, \f$g N m^{-2}\f$
  real,dimension(ilg) :: appl_fert_g            !< grid avg. nitrogen fertilizer, \f$g N m^{-2} day^{-1}\f$
 
!***************************************************************************************************************

  ! Convert litrmass and soilcmas from layered to bulk for use in N cycling:
  soilcmas_bulk(:,:)=0.0
  soilcmas_bulk=sum(soilcmas,dim=3)  ! add over all ignd soil layers

  litrmass_bulk(:,:)=0.0 
  litrmass_bulk=sum(litrmass,dim=3)  ! add over all ignd soil layers

  humtrsvg_bulk(:,:)=0.0 
  humtrsvg_bulk=sum(humtrsvg,dim=3)  ! sum over all soil layers

  ltresveg_bulk(:,:)=0.0 
  ltresveg_bulk=sum(ltresveg,dim=3)  ! sum over all soil layers

  scresveg_bulk(:,:)=0.0 
  scresveg_bulk=sum(scresveg,dim=3)  ! sum over all soil layers

  !> **A.** Initialize all variables 
  temp = 0.0
  temp1 = 0.0
  temp2 = 0.0
  temp3 = 0.0
  ns2t_nleaf_ratio = 0.0
  ns2t_nstem_ratio = 0.0
  ns2t_nroot_ratio = 0.0

  do i = il1,il2
     barefrac(i) = 1.0
     nh4_mass_g(i) = 0.0
     no3_mass_g(i) = 0.0
     pnh4_mass_g(i) = 0.0
     pno3_mass_g(i) = 0.0
     ngleafmas_g(i) = 0.0
     ngleafmas_ns_g(i) = 0.0
     ngleafmas_s_g(i) = 0.0
     pngleafmas_g(i) = 0.0
     pngleafmas_ns_g(i) = 0.0
     pngleafmas_s_g(i)  = 0.0
     nbleafmas_g(i) = 0.0
     pnbleafmas_g(i) = 0.0
     nstemmass_g(i) = 0.0
     nstemmass_ns_g(i) = 0.0
     nstemmass_s_g(i) = 0.0
     pnstemmass_g(i) = 0.0
     pnstemmass_ns_g(i) = 0.0
     pnstemmass_s_g(i) = 0.0
     nrootmass_g(i) = 0.0
     nrootmass_ns_g(i) = 0.0
     nrootmass_s_g(i) = 0.0
     pnrootmass_g(i) = 0.0
     pnrootmass_ns_g(i) = 0.0
     pnrootmass_s_g(i) = 0.0
     nlitrmass_g(i) = 0.0
     pnlitrmass_g(i) = 0.0
     soilnmas_g(i) = 0.0
     psoilnmas_g(i) = 0.0
     nvol(i) = 0.0
     nntch_l(i) = 0.0
     nntch_s(i) = 0.0
     nntch_r(i) = 0.0
     nleafns2s(i) = 0.0
     nstemns2s(i) = 0.0
     nrootns2s(i) = 0.0
     nre_alloc_s2l(i) = 0.0
     nre_alloc_r2l(i) = 0.0
     nresorped_s(i) = 0.0
     nresorped_r(i) = 0.0
     ndemand_wp_npp(i) = 0.0
     nuptake_nh4(i) = 0.0
     nuptake_no3(i) = 0.0
     nuptake_p_nh4(i) = 0.0
     nuptake_p_no3(i) = 0.0
     nlitr_l(i) = 0.0
     nlitr_s(i) = 0.0
     nlitr_r(i) = 0.0
     nalloc_l(i) = 0.0
     nalloc_s(i) = 0.0
     nalloc_r(i) = 0.0
     nleach(i) = 0.0
     nitrif(i) = 0.0
     no_nit(i) = 0.0
     n2o_nit(i) = 0.0
     no_denit(i) = 0.0
     n2o_denit(i) = 0.0
     n2_denit(i) = 0.0
     no_nitdenit(i) = 0.0
     n2o_nitdenit(i) = 0.0
     denit(i) = 0.0
     nuptake_a_actl_nh4(i) = 0.0
     nuptake_a_actl_no3(i) = 0.0
     gl2bl_grass_nflux_gavg(i) = 0.0
     nhumtrs(i) = 0.0
     nmineral_litr(i) = 0.0
     nmineral_humus(i) = 0.0
     nimmobil_nh4(i) = 0.0
     nimmobil_no3(i) = 0.0
     nvgbiomas(i) = 0.0
     fNnetland(i) = 0.0
     gleafmas_g(i) = 0.0
     bleafmas_g(i) = 0.0
     stemmass_g(i) = 0.0
     rootmass_g(i) = 0.0
     litrmass_bulk_g(i) = 0.0
     soilcmas_bulk_g(i) = 0.0
     bnf_free_g(i) = 0.0
     appl_fert_g(i) = 0.0
     do j = 1,iccp1
        appl_fert(i,j) = 0.0
        ndep_nh4(i,j) = 0.0
        ndep_no3(i,j) = 0.0
        nh4_mass_in(i,j) = 0.0
        no3_mass_in(i,j) = 0.0
        fNnetlandveg(i,j) = 0.0
        c2nveg_litr(i,j) = initial_c2n
        c2nveg_humus(i,j) = initial_c2n
        no_nitdenitveg(i,j) = 0.0
        n2o_nitdenitveg(i,j) = 0.0
        !! Convert u mol co2/m2.sec -> \f$kg c/m^2\f$ respired over the model time step
        ltrestep(i,j) = ltresveg_bulk(i,j) * (deltat/963.62)
        screstep(i,j) = scresveg_bulk(i,j) * (deltat/963.62)
        bnf_free(i,j) = 0.0
        bnf_tot(i,j) = 0.0
        if (j <= icc) then
          pftkind = ctempfts(j)
          bnf_nat(i,j) = 0.0
          bnf_ant(i,j) = 0.0
          nstress(i,j) = 0.0
          exsc2n_l(i,j) = 0.0
          exsc2n_s(i,j) = 0.0
          exsc2n_r(i,j) = 0.0
          nallocf_maxc2n_l(i,j) = 0.0
          nallocf_maxc2n_s(i,j) = 0.0
          nallocf_maxc2n_r(i,j) = 0.0
          wexsc2n_wp(i,j) = 0.0
          redcoeff_vcmax(i,j) = 0.0
          nvgbiomas_veg(i,j) = 0.0
          nleafns2sveg(i,j) = 0.0
          nstemns2sveg(i,j) = 0.0
          nrootns2sveg(i,j) = 0.0
          nntchveg_l(i,j) = 0.0
          nntchveg_s(i,j) = 0.0
          nntchveg_r(i,j) = 0.0
          c2nveg_l(i,j) = initial_c2n
          c2nveg_s(i,j) = initial_c2n
          c2nveg_r(i,j) = initial_c2n
          c2nveg_wp(i,j) = initial_c2n

          !> Initial C:N ratios are set to an arbitrary large number (initial_c2n) and are updated when N pools
          !> are a reasonable size

          !! C:N for gleaf mass
          if (ngleafmas(i,j) > 1.E-3 .and. lfstatus(i,j) /= 4 .and. &
              1000. * gleafmas(i,j) > 1.E-3) then
              c2nveg_l(i,j) = (gleafmas(i,j) + bleafmas(i,j))* 1000. / (ngleafmas(i,j) + nbleafmas(i,j))
          end if

          !! C:N for stem mass
          select case (pftkind)

          case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
                'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh') !non-grass

            if (nstemmass(i,j) > 1.E-3 .and. 1000. * stemmass(i,j) > 1.E-3) &
                c2nveg_s(i,j) = stemmass(i,j) * 1000. / nstemmass(i,j)

          case ('GrassC3 ','GrassC4 ','Sedge   ')
            ! do nothing
          case default
            print * ,'Unknown CTEM PFT in n_processes ',pftkind
           call errorHandler('n_processes', - 1)
          end select

          !! C:N for root mass
          if (nrootmass(i,j) > 1.E-3 .and. 1000. * rootmass(i,j) > 1.E-3) &
              c2nveg_r(i,j) = rootmass(i,j) * 1000. / nrootmass(i,j)


          !! C:N for the whole plant
          temp1 = ngleafmas(i,j) + nbleafmas(i,j) + nstemmass(i,j) + nrootmass(i,j)
          temp2 = 1000. * (gleafmas(i,j) + bleafmas(i,j) + stemmass(i,j) + rootmass(i,j))
          if (temp1 > 1.E-3 .and. temp2 > 1.E-3) c2nveg_wp(i,j) = temp2 / temp1
          temp1 = 0.0
          temp2 = 0.0

        end if ! (j<=icc)

        !! C:N for litter mass
        if (nlitrmass(i,j) > 1.E-3 .and. 1000. * litrmass_bulk(i,j) > 1.E-3) &
           c2nveg_litr(i,j) = litrmass_bulk(i,j) * 1000. / nlitrmass(i,j)

        !! C:N for humus mass
        if (soilnmas(i,j) > 1.E-3 .and. 1000. * soilcmas_bulk(i,j) > 1.E-3) &
            c2nveg_humus(i,j) = soilcmas_bulk(i,j) * 1000. / soilnmas(i,j)

     end do ! loop 90
  end do ! loop 80

  !> **B.** Save initial pool sizes to check conservation of N
  do j = 1,iccp1
     do i = il1,il2

        if (abs(nh4_mass(i,j))  < zero) nh4_mass(i,j) = 0.0
        if (abs(no3_mass(i,j))  < zero) no3_mass(i,j) = 0.0
        if (j <= icc) then
          if (abs(ngleafmas(i,j)) < zero) ngleafmas(i,j) = 0.0
          if (abs(nstemmass(i,j)) < zero) nstemmass(i,j) = 0.0
          if (abs(nrootmass(i,j)) < zero) nrootmass(i,j) = 0.0
        end if
        if (abs(nlitrmass(i,j)) < zero) nlitrmass(i,j) = 0.0
        if (abs(soilnmas(i,j))  < zero) soilnmas(i,j) = 0.0

        pno3_mass(i,j) = no3_mass(i,j)
        pnh4_mass(i,j) = nh4_mass(i,j)
        pnlitrmass(i,j)= nlitrmass(i,j)
        psoilnmas(i,j) = soilnmas(i,j)
        if (j <= icc) then
           pngleafmas(i,j) = ngleafmas(i,j)
           pngleafmas_ns(i,j) = ngleafmas_ns(i,j)
           pngleafmas_s(i,j) = ngleafmas_s(i,j)
           pnbleafmas(i,j) = nbleafmas(i,j)
           pnstemmass(i,j) = nstemmass(i,j)
           pnstemmass_ns(i,j) = nstemmass_ns(i,j)
           pnstemmass_s(i,j) = nstemmass_s(i,j)
           pnrootmass(i,j) = nrootmass(i,j)
           pnrootmass_ns(i,j) = nrootmass_ns(i,j)
           pnrootmass_s(i,j) = nrootmass_s(i,j)
        end if
     end do ! loop 110
  end do ! loop 100

  !> **C.** Calculate and save grid averaged initial pool sizes to check conservation of N
  do i = il1,il2
     do j = 1,icc
        pnh4_mass_g(i) = pnh4_mass_g(i) + pnh4_mass(i,j) * fcancmx(i,j)
        pno3_mass_g(i) = pno3_mass_g(i) + pno3_mass(i,j) * fcancmx(i,j)
        pngleafmas_g(i) = pngleafmas_g(i) + pngleafmas(i,j) * fcancmx(i,j)
        pngleafmas_ns_g(i) = pngleafmas_ns_g(i) + pngleafmas_ns(i,j) * fcancmx(i,j)
        pngleafmas_s_g(i) = pngleafmas_s_g(i) + pngleafmas_s(i,j) * fcancmx(i,j)
        pnbleafmas_g(i) = pnbleafmas_g(i) + pnbleafmas(i,j) * fcancmx(i,j)
        pnstemmass_g(i) = pnstemmass_g(i) + pnstemmass(i,j) * fcancmx(i,j)
        pnstemmass_ns_g(i) = pnstemmass_ns_g(i) + pnstemmass_ns(i,j) * fcancmx(i,j)
        pnstemmass_s_g(i) = pnstemmass_s_g(i) + pnstemmass_s(i,j) * fcancmx(i,j)
        pnrootmass_g(i) = pnrootmass_g(i) + pnrootmass(i,j) * fcancmx(i,j)
        pnrootmass_ns_g(i) = pnrootmass_ns_g(i) + pnrootmass_ns(i,j) * fcancmx(i,j)
        pnrootmass_s_g(i) = pnrootmass_s_g(i) + pnrootmass_s(i,j) * fcancmx(i,j)
        pnlitrmass_g(i) = pnlitrmass_g(i) + pnlitrmass(i,j) * fcancmx(i,j)
        psoilnmas_g(i)  = psoilnmas_g(i) + psoilnmas(i,j) * fcancmx(i,j)
        barefrac(i) = barefrac(i) - fcancmx(i,j)
     end do ! loop 130
     pnh4_mass_g(i) = pnh4_mass_g(i) + pnh4_mass(i,iccp1) * barefrac(i)
     pno3_mass_g(i) = pno3_mass_g(i) + pno3_mass(i,iccp1) * barefrac(i)
     pnlitrmass_g(i) = pnlitrmass_g(i) + pnlitrmass(i,iccp1)* barefrac(i)
     psoilnmas_g(i) = psoilnmas_g(i) + psoilnmas(i,iccp1) * barefrac(i)
  end do ! loop 120

  !> **D.** Call the biological nitrogen fixation subroutine
  call bnfix(il1, il2, spinfast, thliq,  zbotw, tbar, isand, & ! In
             THFC, THLW, fcancmx, soilcmas_bulk, & ! In
             soilnmas, & ! In/Out
             bnf_free) ! Out

  !> **E.** Calculate N fertilizer and deposition
  !! N fertilizer flux is only to NH4+ and N deposition flux is to both
  !! NH4+ and NO3-. For tropical regions (between 30 N and 30 S), N fertilizer
  !! is applied daily such that the annual rate read from the external file is spread
  !! equally over all days of the year. For non-tropical regions (> 30 N and > 30 S) fertilizer is
  !! applied during spring and autumn equinoxes.
  do j = 1,iccp1
     do i = il1,il2

      ! Calculate N fertilizer per PFT:
      if (j <= icc) then
         pftkind = ctempfts(j)

         select case (pftkind)

         case ('CropC3  ','CropC4  ')

           ! Northern Hemisphere
           if (radj(i) > pi/6. .and. iday >= fertstd .and. iday < fertend) then

              appl_fert(i,j) = nfertil(i) / abs(fertend-fertstd)

           ! Southern Hemisphere
           else if ( radj(i) < -pi/6. .and. (iday >= fertend .or. iday < fertstd)) then

              appl_fert(i,j) = nfertil(i) / abs(fertend-fertstd)

           ! Tropical regions
           else if (abs(radj(i)) <= pi/6.) then
              if (leapnow) then
                 appl_fert(i,j) = nfertil(i) / 366.
              else
                 appl_fert(i,j) = nfertil(i) / 365.
              end if
           end if

         case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
              'GrassC3 ','GrassC4 ','Sedge   ','BdlEvgSh','BdlDCoSh')
           ! other pfts and bareground, do nothing.

         case default
           print * ,'Unknown CTEM PFT in n_processes ',pftkind
           call errorHandler('n_processes', - 2)
         end select

      end if !for j<=icc

      ! Calculate grid-averaged N fertilizer for balnitro:
      if(j<=icc) then
        appl_fert_g(i) = appl_fert_g(i) + appl_fert(i,j) * fcancmx(i,j)
      end if

      ! Calculate N deposition per PFT + bare soil:
      if (leapnow) then
         ndep_nh4(i,j) = 0.5 * ndeposit(i) / 366.
         ndep_no3(i,j) = 0.5 * ndeposit(i) / 366.
      else
         ndep_nh4(i,j) = 0.5 * ndeposit(i) / 365.
         ndep_no3(i,j) = 0.5 * ndeposit(i) / 365.
      end if

      !! Updating NH4+ for N fertilizer & N deposition
      nh4_mass_in(i,j) = nh4_mass(i,j)

         nh4_mass(i,j) = nh4_mass(i,j) + (appl_fert(i,j) &
                         + ndep_nh4(i,j))

      if (nh4_mass(i,j) < 0.0) then
          write(*,*)'nh4_mass lt zero at i = ',i,' for pft=',j
          write(*,*)'nh4_mass(i,j) = ',nh4_mass(i,j)
          write(*,*)'appl_fert(i,j) = ',appl_fert(i,j)
          write(*,*)'ndep_nh4(i,j) = ',ndep_nh4(i,j)
          call errorHandler('n_processes', - 3)
      end if

         if ( (nh4_mass_in(i,j) - nh4_mass(i,j) + &
                    (appl_fert(i,j) + ndep_nh4(i,j))) &
                     > ntolerance) then
             write(*,*)'imbalance in NH4+ at (i)=',i,'pft=',j
             call errorHandler('n_processes', - 4)
         end if

      !! Updating NO3- for N deposition
      no3_mass_in(i,j) = no3_mass(i,j)

         no3_mass(i,j) = no3_mass(i,j) + ndep_no3(i,j) 

      if (no3_mass(i,j) < 0.0) then
          write(*,*)'no3_mass lt zero at i = ',i,' for pft=',j
          write(*,*)'no3_mass(i,j) = ',no3_mass(i,j)
          write(*,*)'ndep_no3(i,j) = ',ndep_no3(i,j)
          call errorHandler('n_processes', - 5)
      end if

         if ((no3_mass_in(i,j) - no3_mass(i,j) + ndep_no3(i,j)) &
                                         > ntolerance) then
             write(*,*)'imbalance in NO3- at (i)=',i,'pft=',j
             call errorHandler('n_processes', - 6)
         end if

     end do ! loop 150
  end do ! loop 140

  !> **F.** Call the NH3 volatilization subroutine
  call volatilnh3(il1, il2, tbar, CFLUX_GA, USTARBS_GA, & ! In
                  soilpH, thliq, delzw, & ! In
                  nh4_mass, & ! In/Out
                  nvolveg) ! Out

  !> **G.** Call the nitrification subroutine 
  call nitrific(il1, il2, zbotw, tbar, thliq, thice, thpor, & ! In
                THFC, sort, isand, ipeatland, & ! In
                psisat, bi,   & ! In
                nh4_mass, no3_mass, & ! In/Out
                solctempavg, nitrifveg, no_nitveg, n2o_nitveg) ! Out

  !> **H.** Call the denitrification subroutine
  call denitrific(il1, il2, zbotw, isand, tbar, thliq, THFC, & ! In
                  THLW, solctempavg, & ! In
                  no3_mass, & ! In/Out
                  no_denitveg, n2o_denitveg, n2_denitveg, denitveg) ! Out

  !> **I.** Call the leaching subroutine 
  call leaching(il1, il2, ROFB, & ! In
                no3_mass, & ! In/Out
                nleachveg) ! Out

  !> **J.** Call the plant N uptake subroutine 
  call nuptake(il1, il2, sort, fcancmx, zbotw, delzw, & ! In
               thliq, rootdpth, QFC, ntchlveg, ntchsveg, ntchrveg, & ! In
               gleafmas, stemmass, rootmass, ngleafmas, nstemmass, & ! In
               lfstatus, iday, & ! In
               nrootmass, nh4_mass, no3_mass, & ! In/Out
               ndemandveg_l_npp, ndemandveg_s_npp, ndemandveg_r_npp, & !Out
               ndemandveg_wp_npp, nuptakeveg_p_nh4, & !Out
               nuptakeveg_p_no3, nuptakeveg_a_actl_nh4, & !Out
               nuptakeveg_a_actl_no3, nuptakeveg_nh4, nuptakeveg_no3) ! Out

  !> **K.** Call the plant N allocation subroutine
  call nallocate(il1, il2, fcancmx, lfstatus, sort, ailcg, & ! In
                 lfthrs, radj, iday, ntchlveg, & ! In
                 ntchsveg, ntchrveg, tltrleaf, ndemandveg_l_npp, & ! In
                 ndemandveg_s_npp, ndemandveg_r_npp, ndemandveg_wp_npp, & ! In
                 nuptakeveg_nh4, & ! In
                 nuptakeveg_no3, gleafmas, rootmass_ns, & ! In
                 re_alloc_sr2l, leafns2s, stemns2s, rootns2s, tbar, isand, & ! In
                 ngleafmas, ngleafmas_ns, ngleafmas_s, & ! In/Out
                 nstemmass, nstemmass_ns, nstemmass_s, nrootmass, nrootmass_ns, & ! In/Out
                 nrootmass_s, nleafns2sveg, nstemns2sveg, nrootns2sveg, & ! In/Out
                 nntchveg_l, nntchveg_s, nntchveg_r, nallocveg_l, nallocveg_s, & ! Out
                 nallocveg_r, nresorpedveg_s, nresorpedveg_r, nre_allocveg_s2l, & ! Out
                 nre_allocveg_r2l, bnf_nat, bnf_ant, nstress) ! Out

  !> **L.** Call the N litter fall subroutine
  call nlitter(il1, il2, sort, fcancmx, lfstatus, iday, tltrleaf, & ! In
               tltrstem, tltrroot, gl2bl_grass_cflux, bl2ltr_grass_cflux, & ! In
               gleafmas, bleafmas, stemmass, rootmass, & ! In
               ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! In/Out
               nstemmass, nstemmass_ns, nstemmass_s, nrootmass, & ! In/Out
               nrootmass_ns, nrootmass_s, nlitrmass, & ! In/Out
               gl2bl_grass_nflux, nlitrveg_l, nlitrveg_s, nlitrveg_r) !Out

  !> **M.** Call the N humification subroutine
  call nhumific(il1, il2, sort, fcancmx, spinfast, humtrsvg_bulk, c2nveg_litr, & ! In
                nlitrmass, soilnmas, & ! In/Out
                nhumtrsveg) ! Out

  !> **N.** Call the N mineralization & immobilization subroutine 
  call nmineralimm(il1, il2, spinfast, soilcmas_bulk, c2nveg_litr, c2nveg_humus, & ! In
                   ltrestep, screstep, litrmass_bulk, co2conc, & ! In
                   nlitrmass, soilnmas, nh4_mass, no3_mass, & ! In/Out
                   nmineralveg_litr, nmineralveg_humus, & ! Out
                   nimmobilveg_nh4, nimmobilveg_no3) ! Out

  !> **O.** Update vegetation N & calculate net terrestrial N flux
  do i = il1,il2
    do j = 1,icc
      nvgbiomas_veg(i,j) = ngleafmas(i,j) + nbleafmas(i,j) &
                           + nstemmass(i,j) + nrootmass(i,j)
      bnf_tot(i,j) = bnf_free(i,j) + bnf_nat(i,j) + bnf_ant(i,j)
      fNnetlandveg(i,j) = bnf_tot(i,j) &
                         + appl_fert(i,j) &
                         + (ndep_nh4(i,j) + ndep_no3(i,j)) - nvolveg(i,j) &
                         - nleachveg(i,j) - no_nitveg(i,j) &
                         - n2o_nitveg(i,j) - no_denitveg(i,j) &
                         - n2o_denitveg(i,j) - n2_denitveg(i,j)
    end do
    bnf_tot(i,iccp1) = bnf_free(i,iccp1)
    fNnetlandveg(i,iccp1) = bnf_tot(i,iccp1) &
                            + appl_fert(i,iccp1) &
                            + (ndep_nh4(i,iccp1) + ndep_no3(i,iccp1)) - nvolveg(i,iccp1) &
                            - nleachveg(i,iccp1) - no_nitveg(i,iccp1) &
                            - n2o_nitveg(i,iccp1) - no_denitveg(i,iccp1) &
                            - n2o_denitveg(i,iccp1) - n2_denitveg(i,iccp1)
  end do


  ! Avoid small N pools: 
  do i = il1,il2
    do j = 1,iccp1
        if (abs(nh4_mass(i,j)) < zero) nh4_mass(i,j) = 0.0
        if (abs(no3_mass(i,j)) < zero) no3_mass(i,j) = 0.0
        if (j <= icc) then
          if (abs(ngleafmas(i,j)) < zero) ngleafmas(i,j) = 0.0
          if (abs(nstemmass(i,j)) < zero) nstemmass(i,j) = 0.0
          if (abs(nrootmass(i,j)) < zero) nrootmass(i,j) = 0.0
        end if
        if (abs(nlitrmass(i,j)) < zero) nlitrmass(i,j) = 0.0
        if (abs(soilnmas(i,j)) < zero) soilnmas(i,j) = 0.0
    end do ! loop 170
  end do ! loop 160

  !> **P.** Calculate grid-averaged pool sizes and fluxes
  do i = il1,il2
     barefrac(i) = 1.0
     do j = 1,icc
        nh4_mass_g(i) = nh4_mass_g(i) + nh4_mass(i,j) * fcancmx(i,j)
        no3_mass_g(i) = no3_mass_g(i) + no3_mass(i,j) * fcancmx(i,j)
        ngleafmas_g(i) = ngleafmas_g(i) + ngleafmas(i,j) * fcancmx(i,j)
        ngleafmas_ns_g(i) = ngleafmas_ns_g(i) + ngleafmas_ns(i,j) * fcancmx(i,j)
        ngleafmas_s_g(i) = ngleafmas_s_g(i) + ngleafmas_s(i,j) * fcancmx(i,j)
        nbleafmas_g(i) = nbleafmas_g(i) + nbleafmas(i,j) * fcancmx(i,j)
        nstemmass_g(i) = nstemmass_g(i) + nstemmass(i,j) * fcancmx(i,j)
        nstemmass_ns_g(i) = nstemmass_ns_g(i) + nstemmass_ns(i,j) * fcancmx(i,j)
        nstemmass_s_g(i) = nstemmass_s_g(i) + nstemmass_s(i,j)  * fcancmx(i,j)
        nrootmass_g(i) = nrootmass_g(i) + nrootmass(i,j) * fcancmx(i,j)
        nrootmass_ns_g(i) = nrootmass_ns_g(i) + nrootmass_ns(i,j) * fcancmx(i,j)
        nrootmass_s_g(i) = nrootmass_s_g(i) + nrootmass_s(i,j) * fcancmx(i,j)
        nlitrmass_g(i) = nlitrmass_g(i) + nlitrmass(i,j) * fcancmx(i,j)
        soilnmas_g(i) = soilnmas_g(i) + soilnmas(i,j) * fcancmx(i,j)
        nvgbiomas(i) = nvgbiomas(i) + ( ngleafmas(i,j) + nbleafmas(i,j) &
                       + nstemmass(i,j) + nrootmass(i,j)) * fcancmx(i,j)
        no_nitdenitveg(i,j) = no_nitveg(i,j) + no_denitveg(i,j)
        n2o_nitdenitveg(i,j) = n2o_nitveg(i,j) + n2o_denitveg(i,j)
        nvol(i) = nvol(i) + nvolveg(i,j) * fcancmx(i,j)
        nitrif(i) = nitrif(i) + nitrifveg(i,j) * fcancmx(i,j)
        no_nit(i) = no_nit(i) + no_nitveg(i,j) * fcancmx(i,j)
        n2o_nit(i) = n2o_nit(i) + n2o_nitveg(i,j) * fcancmx(i,j)
        denit(i) = denit(i) + denitveg(i,j) * fcancmx(i,j)
        no_denit(i) = no_denit(i) + no_denitveg(i,j) * fcancmx(i,j)
        n2o_denit(i) = n2o_denit(i) + n2o_denitveg(i,j) * fcancmx(i,j)
        n2_denit(i) = n2_denit(i) + n2_denitveg(i,j) * fcancmx(i,j)
        nuptake_p_nh4(i) = nuptake_p_nh4(i) + nuptakeveg_p_nh4(i,j) * fcancmx(i,j)
        nuptake_p_no3(i) = nuptake_p_no3(i) + nuptakeveg_p_no3(i,j) * fcancmx(i,j)
        nuptake_a_actl_nh4(i) = nuptake_a_actl_nh4(i) + nuptakeveg_a_actl_nh4(i,j) * fcancmx(i,j)
        nuptake_a_actl_no3(i) = nuptake_a_actl_no3(i) + nuptakeveg_a_actl_no3(i,j) * fcancmx(i,j)
        nuptake_nh4(i) = nuptake_nh4(i) + nuptakeveg_nh4(i,j) * fcancmx(i,j)
        nuptake_no3(i) = nuptake_no3(i) + nuptakeveg_no3(i,j) * fcancmx(i,j)
        nuptakeveg(i,j) = nuptakeveg_p_nh4(i,j) + nuptakeveg_p_no3(i,j) + nuptakeveg_a_actl_nh4(i,j) + nuptakeveg_a_actl_no3(i,j)
        nresorped_s(i) = nresorped_s(i) + nresorpedveg_s(i,j) * fcancmx(i,j)
        nresorped_r(i) = nresorped_r(i) + nresorpedveg_r(i,j) * fcancmx(i,j)
        nre_alloc_s2l(i) = nre_alloc_s2l(i) + nre_allocveg_s2l(i,j) * fcancmx(i,j)
        nre_alloc_r2l(i) = nre_alloc_r2l(i) + nre_allocveg_r2l(i,j) * fcancmx(i,j)
        nntch_l(i) = nntch_l(i) + nntchveg_l(i,j) * fcancmx(i,j)
        nntch_s(i) = nntch_s(i) + nntchveg_s(i,j) * fcancmx(i,j)
        nntch_r(i) = nntch_r(i) + nntchveg_r(i,j) * fcancmx(i,j)
        nleafns2s(i) = nleafns2s(i) + nleafns2sveg(i,j) * fcancmx(i,j)
        nstemns2s(i) = nstemns2s(i) + nstemns2sveg(i,j) * fcancmx(i,j)
        nrootns2s(i) = nrootns2s(i) + nrootns2sveg(i,j) * fcancmx(i,j)
        nlitr_l(i) = nlitr_l(i) + nlitrveg_l(i,j) * fcancmx(i,j)
        nlitr_s(i) = nlitr_s(i) + nlitrveg_s(i,j) * fcancmx(i,j)
        nlitr_r(i) = nlitr_r(i) + nlitrveg_r(i,j) * fcancmx(i,j)
        nlitrveg(i,j) = nlitrveg_l(i,j) + nlitrveg_s(i,j) + nlitrveg_r(i,j)
        gl2bl_grass_nflux_gavg(i) = gl2bl_grass_nflux_gavg(i) &
                                    + gl2bl_grass_nflux(i,j) * fcancmx(i,j)
        nalloc_l(i) = nalloc_l(i) + nallocveg_l(i,j) * fcancmx(i,j)
        nalloc_s(i) = nalloc_s(i) + nallocveg_s(i,j) * fcancmx(i,j)
        nalloc_r(i) = nalloc_r(i) + nallocveg_r(i,j) * fcancmx(i,j)
        nleach(i) = nleach(i) + nleachveg(i,j) * fcancmx(i,j)
        nhumtrs(i) = nhumtrs(i) + nhumtrsveg(i,j) * fcancmx(i,j)
        nmineral_litr(i) = nmineral_litr(i) + nmineralveg_litr(i,j) * fcancmx(i,j)
        nmineral_humus(i) = nmineral_humus(i) + nmineralveg_humus(i,j) * fcancmx(i,j)
        netnmineralveg(i,j) = nmineralveg_litr(i,j) + nmineralveg_humus(i,j) - nimmobilveg_nh4(i,j) - nimmobilveg_no3(i,j)
        nimmobil_nh4(i) = nimmobil_nh4(i) + nimmobilveg_nh4(i,j) * fcancmx(i,j)
        nimmobil_no3(i) = nimmobil_no3(i) + nimmobilveg_no3(i,j) * fcancmx(i,j)
        fNnetland(i) = fNnetland(i) + fNnetlandveg(i,j) * fcancmx(i,j)
        gleafmas_g(i) = gleafmas_g(i) + gleafmas(i,j) * fcancmx(i,j)
        bleafmas_g(i) = bleafmas_g(i) + bleafmas(i,j) * fcancmx(i,j)
        stemmass_g(i) = stemmass_g(i) + stemmass(i,j) * fcancmx(i,j)
        rootmass_g(i) = rootmass_g(i) + rootmass(i,j) * fcancmx(i,j)
        litrmass_bulk_g(i) = litrmass_bulk_g(i) + litrmass_bulk(i,j) * fcancmx(i,j)
        soilcmas_bulk_g(i) = soilcmas_bulk_g(i) + soilcmas_bulk(i,j) * fcancmx(i,j)
        bnf_free_g(i) = bnf_free_g(i) +  bnf_free(i,j) * fcancmx(i,j)
        ndemand_wp_npp(i) = ndemand_wp_npp(i) + ndemandveg_wp_npp(i,j) * fcancmx(i,j)
        barefrac(i) = barefrac(i) - fcancmx(i,j)

     end do ! loop 190

     ! bare soil fraction
     no_nitdenitveg(i,iccp1) = no_nitveg(i,iccp1) + no_denitveg(i,iccp1)
     n2o_nitdenitveg(i,iccp1) = n2o_nitveg(i,iccp1) + n2o_denitveg(i,iccp1)
     nh4_mass_g(i) = nh4_mass_g(i) + nh4_mass(i,iccp1) * barefrac(i)
     no3_mass_g(i) = no3_mass_g(i) + no3_mass(i,iccp1) * barefrac(i)
     nlitrmass_g(i) = nlitrmass_g(i) + nlitrmass(i,iccp1) * barefrac(i)
     soilnmas_g(i) = soilnmas_g(i) + soilnmas(i,iccp1) * barefrac(i)
     nvol(i) = nvol(i) + nvolveg(i,iccp1) * barefrac(i)
     nleach(i) = nleach(i) + nleachveg(i,iccp1) * barefrac(i)
     nhumtrs(i) = nhumtrs(i) + nhumtrsveg(i,iccp1) * barefrac(i)
     nitrif(i) = nitrif(i) + nitrifveg(i,iccp1) * barefrac(i)
     no_nit(i) = no_nit(i) + no_nitveg(i,iccp1) * barefrac(i)
     n2o_nit(i) = n2o_nit(i) + n2o_nitveg(i,iccp1) * barefrac(i)
     denit(i) = denit(i) + denitveg(i,iccp1) * barefrac(i)
     no_denit(i) = no_denit(i) + no_denitveg(i,iccp1) * barefrac(i)
     n2o_denit(i) = n2o_denit(i) + n2o_denitveg(i,iccp1) * barefrac(i)
     n2_denit(i) = n2_denit(i) + n2_denitveg(i,iccp1) * barefrac(i)
     nmineral_litr(i) = nmineral_litr(i) + nmineralveg_litr(i,iccp1) * barefrac(i)
     nmineral_humus(i) = nmineral_humus(i) + nmineralveg_humus(i,iccp1) * barefrac(i)
     netnmineralveg(i,iccp1) = nmineralveg_litr(i,iccp1) + nmineralveg_humus(i,iccp1) - nimmobilveg_nh4(i,iccp1) - nimmobilveg_no3(i,iccp1)
     nimmobil_nh4(i) = nimmobil_nh4(i) + nimmobilveg_nh4(i,iccp1) * barefrac(i)
     nimmobil_no3(i) = nimmobil_no3(i) + nimmobilveg_no3(i,iccp1) * barefrac(i)
     fNnetland(i) = fNnetland(i) + fNnetlandveg(i,iccp1) * barefrac(i)
     no_nitdenit(i) = no_nit(i) + no_denit(i)
     n2o_nitdenit(i) = n2o_nit(i) + n2o_denit(i)
     litrmass_bulk_g(i) = litrmass_bulk_g(i) + litrmass_bulk(i,iccp1) * barefrac(i)
     soilcmas_bulk_g(i) = soilcmas_bulk_g(i) + soilcmas_bulk(i,iccp1) * barefrac(i)
     bnf_free_g(i) = bnf_free_g(i) +  bnf_free(i,iccp1) * barefrac(i)
  end do ! loop 180

  !> **Q.** Calculate updated C:N ratios 
  do i = il1,il2
     barefrac(i) = 1.0
     do j = 1,icc
        pftkind = ctempfts(j)
        barefrac(i) = barefrac(i) - fcancmx(i,j)

      !! leaf C:N ratio
      if (ngleafmas(i,j) > 1.E-3 .and. 1000. * gleafmas(i,j) > 1.E-3 .and. &
                                                    lfstatus(i,j) /= 4) then
         c2nveg_l(i,j) = (gleafmas(i,j) + bleafmas(i,j))* 1000. / (ngleafmas(i,j) + nbleafmas(i,j))
      end if

      !! stem C:N ratio
      select case (pftkind)
      case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr',&
            'BdlDDrTr','CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
         if (nstemmass(i,j) > 1.E-3 .and. 1000.*stemmass(i,j) > 1.E-3) then
            c2nveg_s(i,j) = stemmass(i,j) * 1000. / nstemmass(i,j)
         end if
      case ('GrassC3 ','GrassC4 ','Sedge   ')
           ! do nothing
      case default
           print * ,'Unknown CTEM PFT in n_processes ',pftkind
           call errorHandler('nprocesses', - 7)
      end select

      !! root C:N ratio
      if (nrootmass(i,j) > 1.E-3 .and. 1000. * rootmass(i,j) > 1.E-3) then
         c2nveg_r(i,j) = rootmass(i,j) * 1000. / nrootmass(i,j)
      end if

      !! whole plant C:N ratio
      temp1 = ngleafmas(i,j) + nbleafmas(i,j) + nstemmass(i,j) + nrootmass(i,j)
      temp2 = 1000. * (gleafmas(i,j) + bleafmas(i,j) + stemmass(i,j) + rootmass(i,j))
      if (temp1 > 1.E-3 .and. temp2 > 1.E-3) then
         c2nveg_wp(i,j) = temp2 / temp1
      end if
      temp1 = 0.0
      temp2 = 0.0

      !! litter C:N ratio
      if (nlitrmass(i,j) > 1.E-3 .and. 1000.*litrmass_bulk(i,j) > 1.E-3) then
         c2nveg_litr(i,j) = litrmass_bulk(i,j) * 1000. / nlitrmass(i,j)
      end if

      !! humus C:N ratio
      if (soilnmas(i,j) > 1.E-3 .and. 1000.*soilcmas_bulk(i,j) > 1.E-3) then
         c2nveg_humus(i,j) = soilcmas_bulk(i,j) * 1000. / soilnmas(i,j)
      end if

      !! check for negative ratios
      if (c2nveg_l(i,j) < 0.0) then
         write(*,*)'c2nveg_l lt zero for i= ',i,' and pft= ',j
         call errorHandler('n_processes', - 8)
      else if (c2nveg_s(i,j) < 0.0) then
         write(*,*)'c2nveg_s lt zero for i= ',i,' and pft= ',j
         call errorHandler('n_processes', - 9)
      else if (c2nveg_r(i,j) < 0.0) then
         write(*,*)'c2nveg_r lt zero for i= ',i,' and pft= ',j
         call errorHandler('n_processes', - 10)
      else if (c2nveg_litr(i,j) < 0.0) then
         write(*,*)'c2nveg_litr lt zero for i= ',i,' and pft= ',j
         call errorHandler('n_processes', - 11)
      else if (c2nveg_humus(i,j) < 0.0) then
         write(*,*)'c2nveg_humus lt zero for i= ',i,' and pft= ',j
         call errorHandler('n_processes', - 12)
      end if

     end do ! loop 210

     ! bare soil 
     if (nlitrmass(i,iccp1) > 1.E-3 .and. 1000. * litrmass_bulk(i,iccp1) > 1.E-3) then
        c2nveg_litr(i,iccp1) = litrmass_bulk(i,iccp1) * 1000. / nlitrmass(i,iccp1)
     end if
     if (soilnmas(i,iccp1) > 1.E-3 .and. 1000.*soilcmas_bulk(i,iccp1) > 1.E-3) then
        c2nveg_humus(i,iccp1) = soilcmas_bulk(i,iccp1) * 1000. / soilnmas(i,iccp1)
     end if

  end do ! loop 200

  !> **R.** Calculate reduction coefficient of V_c,max when C:N ratios for the leaves
  !! stem, and root components exceed their maximum specified values.
  !! (see equations 31-33 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ V_{cmax} = \Lambda(\frac{\Gamma_1N_L}{3} + \Gamma_2) \f$
  do j = 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2

        select case (pftkind)

        case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
              'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
           !! Excesses above maximum C:N ratios
           if (c2nveg_l(i,j) < 2000. .and. c2nveg_l(i,j) > c2n_lmax(sort(j)))&
                            exsc2n_l(i,j) = c2nveg_l(i,j) - c2n_lmax(sort(j))
      	   if (c2nveg_s(i,j) < 2000. .and. c2nveg_s(i,j) > c2n_smax(sort(j)))&
                            exsc2n_s(i,j) = c2nveg_s(i,j) - c2n_smax(sort(j))
           if (c2nveg_r(i,j) < 2000. .and. c2nveg_r(i,j) > c2n_rmax(sort(j)))&
                            exsc2n_r(i,j) = c2nveg_r(i,j) - c2n_rmax(sort(j))

           !! Whole plant weighted excess above maximum C:N ratio
           temp3 = (1. / c2n_lmax(sort(j))) + (1. / c2n_smax(sort(j))) + (1. / c2n_rmax(sort(j)))
           if (temp3 > 0.0) then
              nallocf_maxc2n_l(i,j) = (1. / c2n_lmax(sort(j))) / temp3
              nallocf_maxc2n_s(i,j) = (1. / c2n_smax(sort(j))) / temp3
              nallocf_maxc2n_r(i,j) = (1. / c2n_rmax(sort(j))) / temp3
           end if

           !! Whole plant's weighted excess-above max C:N
           wexsc2n_wp(i,j) = (exsc2n_l(i,j) * nallocf_maxc2n_l(i,j)) &
                             + (exsc2n_s(i,j) * nallocf_maxc2n_s(i,j)) &
                             + (exsc2n_r(i,j) * nallocf_maxc2n_r(i,j))

        case ('GrassC3 ','GrassC4 ','Sedge   ')
           !! Excesses above maximum C:N ratios
           if (c2nveg_l(i,j) < 2000. .and. c2nveg_l(i,j) > c2n_lmax(sort(j)))&
                            exsc2n_l(i,j) = c2nveg_l(i,j) - c2n_lmax(sort(j))
           if (c2nveg_r(i,j) < 2000. .and. c2nveg_r(i,j) > c2n_rmax(sort(j)))&
                            exsc2n_r(i,j) = c2nveg_r(i,j) - c2n_rmax(sort(j))

           !! Whole plant weighted excess above maximum C:N ratio
           temp3= (1. / c2n_lmax(sort(j))) + (1. / c2n_rmax(sort(j)))
           if (temp3 > 0.0) then
              nallocf_maxc2n_l(i,j) = (1. / c2n_lmax(sort(j))) / temp3
              nallocf_maxc2n_r(i,j) = (1. / c2n_rmax(sort(j))) / temp3
           end if

           !! Whole plant's weighted excess-above max C:N
           wexsc2n_wp(i,j) = (exsc2n_l(i,j) * nallocf_maxc2n_l(i,j))&
                            + (exsc2n_r(i,j) * nallocf_maxc2n_r(i,j))

        case default
           print * ,'Unknown CTEM PFT in nallocate ',pftkind
           call errorHandler('nprocesses', - 14)
        end select

        redcoeff_vcmax(i,j) = max(0.,min(1.,exp(-k_redcoeff_vcmax * wexsc2n_wp(i,j))))
        temp3 = 0.0
     end do ! loop 220
  end do ! loop 230

  !> **S.** Check if all N pools and fluxes are in balance and we are conserving N to within some specified tolerance
  !> The pools are leaf, stem, root, litter, soil, ammonium, and nitrate 
  if (spinfast == 1) then
     call balnit(il1, il2, iday, leapnow, radj, appl_fert, appl_fert_g, ndeposit,                         &
         bnf_free, bnf_free_g,                                                                            &
         nh4_mass, pnh4_mass, nh4_mass_g, pnh4_mass_g, no3_mass, pno3_mass, no3_mass_g, pno3_mass_g,      &
         ngleafmas, pngleafmas, ngleafmas_g, pngleafmas_g, ngleafmas_ns, pngleafmas_ns, ngleafmas_ns_g,   &
         pngleafmas_ns_g, ngleafmas_s, pngleafmas_s, ngleafmas_s_g, pngleafmas_s_g, nbleafmas,            &
         pnbleafmas, nbleafmas_g, pnbleafmas_g, nstemmass, pnstemmass, nstemmass_g, pnstemmass_g,         &
         nstemmass_ns, pnstemmass_ns, nstemmass_ns_g, pnstemmass_ns_g, nstemmass_s, pnstemmass_s,         &
         nstemmass_s_g, pnstemmass_s_g, nrootmass, pnrootmass, nrootmass_g, pnrootmass_g, nrootmass_ns,   &
         pnrootmass_ns, nrootmass_ns_g, pnrootmass_ns_g, nrootmass_s, pnrootmass_s, nrootmass_s_g,        &
         pnrootmass_s_g, nlitrmass, pnlitrmass, nlitrmass_g, pnlitrmass_g, soilnmas, psoilnmas,           &
         soilnmas_g, psoilnmas_g, nleachveg, nleach, nvolveg, nvol, nitrifveg, nitrif, denitveg, denit,   &
         no_nitveg, no_nit, n2o_nitveg, n2o_nit, nuptakeveg_nh4, nuptake_nh4, nuptakeveg_no3,             &
         nuptake_no3, ntchlveg, nleafns2sveg, nstemns2sveg, nrootns2sveg,                                 &
         nntchveg_l, nntch_l, nntchveg_s, nntch_s, nntchveg_r, nntch_r, lfstatus, nlitrveg_l, nlitr_l,    &
         nlitrveg_s, nlitr_s, nlitrveg_r, nlitr_r, gl2bl_grass_nflux, gl2bl_grass_nflux_gavg,             & 
         nhumtrsveg, nhumtrs, nmineralveg_litr, nmineral_litr, nmineralveg_humus, nmineral_humus,         &
         nimmobilveg_nh4, nimmobil_nh4, nimmobilveg_no3, nimmobil_no3)
  end if

  return
end subroutine n_processes
