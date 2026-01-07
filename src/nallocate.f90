!> \file
!! Performs allocation of nitrogen taken up by the plant into its different tissues.
!> @author A. Asaadi
!!


subroutine nallocate(il1, il2, fcancmx, lfstatus, sort, ailcg, lfthrs, radj, & ! In
                     iday, ntchlveg, ntchsveg, ntchrveg, tltrleaf, & ! In
                     ndemandveg_l_npp, ndemandveg_s_npp, ndemandveg_r_npp, & ! In
                     ndemandveg_wp_npp, nuptakeveg_nh4, & ! In
                     nuptakeveg_no3, gleafmas, rootmass_ns, & ! In
                     re_alloc_sr2l, leafns2s, & ! In
                     stemns2s, rootns2s, tbar, isand, & ! In
                     ngleafmas, ngleafmas_ns, ngleafmas_s, nstemmass, nstemmass_ns, & ! In/Out
                     nstemmass_s, nrootmass, nrootmass_ns, nrootmass_s, nleafns2sveg, & ! In/Out
                     nstemns2sveg, nrootns2sveg, & ! In/Out
                     nntchveg_l, nntchveg_s, nntchveg_r, nallocveg_l, nallocveg_s, & ! Out
                     nallocveg_r, nresorpedveg_s, nresorpedveg_r, nre_allocveg_s2l, & ! Out
                     nre_allocveg_r2l, bnf_nat, bnf_ant, nstress) ! Out


  use classicParams, only : ilg, ignd, icc, iccp1, ntolerance, ctempfts, &
                            deltat, nresorp_coeff, min_ns2t_l, min_ns2t_s, &
                            min_ns2t_r, zero, c2n_lmin, c2n_smin, c2n_rmin, &
                            Fmax_leaf_N, Fmax_stem_N, Fmax_root_N, r_bnf_s, b_bnf, max_frac_bnf, C_cost_bnf_s

  implicit none

  integer,intent(in) :: il1       !< il1=1
  integer,intent(in) :: il2       !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in) :: iday      !< day of year
  integer,intent(in),dimension(ilg,icc) :: lfstatus          !< integer indicating leaf status or mode (max. growth=1, normal growth=2, fall=3, no leaves=4)
  integer,intent(in),dimension(icc) :: sort                  !< index for correspondence between PFTs and the 12 values in parameters vectors
  real,intent(in),dimension(ilg) :: radj                     !< latitude in radians
  real,intent(in),dimension(ilg,icc) :: fcancmx              !< max. fractional coverage of PFTs
  real,intent(in),dimension(ilg,icc) :: ndemandveg_l_npp     !< leaf nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$ 
  real,intent(in),dimension(ilg,icc) :: ndemandveg_s_npp     !< stem nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: ndemandveg_r_npp     !< root nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: ndemandveg_wp_npp    !< whole plant nitrogen demand for each PFT  based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: nuptakeveg_nh4       !< total nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: nuptakeveg_no3       !< total no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: gleafmas             !< green leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc) :: rootmass_ns    !< non-structural root carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc) :: re_alloc_sr2l        !< carbon from stem and root allocated to leaf during leaf onset for each PFT, \f$g C m^{-2} day^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: ailcg                !< green lai
  real,intent(in),dimension(ilg,icc) :: lfthrs               !< threshold lai to determine leaf status 
  real,intent(in),dimension(ilg,icc) :: leafns2s             !< carbon flux from non-structural to structural leaf carbon, \f$umol CO2 m^{-2} sec^{-1}\f$ 
  real,intent(in),dimension(ilg,icc) :: stemns2s             !< carbon flux from non-structural to structural stem carbon, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: rootns2s             !< carbon flux from non-structural to structural root carbon, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: tltrleaf             !< total leaf carbon litter fall rate, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: ntchlveg             !< C allocation to leaf, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: ntchsveg             !< C allocation to stem, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: ntchrveg             !< C allocation to root, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in), dimension(ilg,ignd) :: tbar               !< soil temperature, K
  integer,intent(in), dimension(ilg,ignd) :: isand           !< flag for bedrock or ice in a soil layer

  real,intent(inout),dimension(ilg,icc) :: ngleafmas         !< green leaf nitrogen for each PFT, \f$g N m^{-2}\f$ 
  real,intent(inout),dimension(ilg,icc) :: ngleafmas_ns      !< non-structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: ngleafmas_s       !< structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nstemmass         !< stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nstemmass_ns      !< non-structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nstemmass_s       !< structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nrootmass         !< root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nrootmass_ns      !< non-structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nrootmass_s       !< structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nleafns2sveg      !< nitrogen flux from non-structural to structural leaf nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,intent(inout),dimension(ilg,icc) :: nstemns2sveg      !< nitrogen flux from non-structural to structural stem nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,intent(inout),dimension(ilg,icc) :: nrootns2sveg      !< nitrogen flux from non-structural to structural root nitrogen, \f$g N m^{-2} day^{-1}\f$

  real,intent(out),dimension(ilg,icc) :: nallocveg_l         !< nitrogen allocation to leaf for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nallocveg_s         !< nitrogen allocation to stem for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nallocveg_r         !< nitrogen allocation to root for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nresorpedveg_s      !< nitrogen from leaf allocated to stem during leaf offset for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nresorpedveg_r      !< nitrogen from leaf allocated to root during leaf offset for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nre_allocveg_s2l    !< nitrogen from stem allocated to leaf during leaf onset for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nre_allocveg_r2l    !< nitrogen from root allocated to leaf during leaf onset for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nntchveg_l          !< nitrogen net change in ngleafmas_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nntchveg_s          !< nitrogen net change in nstemmass_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nntchveg_r          !< nitrogen net change in nrootmass_ns after nitrogen allocation and reallocation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: bnf_nat             !< natural symbiotic biological nitrogen fixation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: bnf_ant             !< anthropogenic symbiotic biological nitrogen fixation for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nstress             !< plant nitrogen stress

  ! local variables
  integer :: i, j
  character(8) :: pftkind
  real :: temp_cn             !< temporary denominator of nitrogen allocation (for NPP < 0)
  real :: temp_ndemand        !< temporary denominator of nitrogen allocation (during leaf offset)
  real :: temp_reallocation_onset   !< temporary denominator of nitrogen reallocation during leaf onset
  real :: temp_reallocation_offset  !< temporary denominator of nitrogen reallocation during leaf offset
  real :: temp_reallocation   !< temporary denomination of nitrogen reallocation when leaf C:N ratio less than minimum leaf C:N ratio 
  real :: gleafltr_fr         !< total leaf carbon litter fall as a fraction of green leaf carbon 
  real :: ns2t_nleaf_ratio    !< temporary ratio of non-structural leaf nitrogen to total leaf nitrogen 
  real :: ns2t_nstem_ratio    !< temporary ratio of non-structural stem nitrogen to total stem nitrogen 
  real :: ns2t_nroot_ratio    !< temporary ratio of non-structural root nitrogen to total root nitrogen 
  real :: ngleafmas_resorped  !< nitrogen from leaf allocated to stem and root during leaf offset, \f$g N m^{-2}\f$
  real :: nstemns2t           !< temporary ratio of non-structural stem nitrogen to total stem nitrogen 
  real :: nrootns2t           !< temporary ratio of non-structural root nitrogen to total root nitrogen
  real :: ns_r2sr             !< temporary ratio of non-structural stem nitrogen to non-structural stem and root nitrogen 
  real :: ns_s2sr             !< temporary ratio of non-structural root nitrogen to non-structural stem and root nitrogen 
  real,dimension(ilg,icc) :: nre_allocveg_sr2l !< total nitrogen from stem and root allocated to leaf during leaf onset for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: ngleafmas_in      !< previous green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: ngleafmas_ns_in   !< previous non-structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: ngleafmas_s_in    !< previous structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nstemmass_in      !< previous stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nstemmass_ns_in   !< previous non-structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nstemmass_s_in    !< previous structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nrootmass_in      !< previous root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nrootmass_ns_in   !< previous non-structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nrootmass_s_in    !< previous structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: temp_ns           !< temporary nitrogen flux to non-structural nitrogen 
  real,dimension(ilg,icc) :: temp_s            !< temporary nitrogen flux to structural nitrogen
  temp_cn = 0.0
  temp_ndemand = 0.0
  temp_reallocation_onset = 0.0
  temp_reallocation_offset = 0.0
  temp_reallocation = 0.0
  gleafltr_fr = 0.0
  ngleafmas_resorped = 0.0
  ns2t_nleaf_ratio = 0.0
  ns2t_nstem_ratio = 0.0
  ns2t_nroot_ratio = 0.0
  do i = il1,il2
     do j = 1,icc
        pftkind = ctempfts(j)
        temp_ns(i,j) = 0.0
        temp_s(i,j) = 0.0
        nallocveg_l(i,j) = 0.0
        nallocveg_s(i,j) = 0.0
        nallocveg_r(i,j) = 0.0
        nntchveg_l(i,j) = 0.0
        nntchveg_s(i,j) = 0.0
        nntchveg_r(i,j) = 0.0
        nleafns2sveg(i,j) = 0.0
        nstemns2sveg(i,j) = 0.0
        nrootns2sveg(i,j) = 0.0
        ngleafmas_in(i,j) = 0.0
        ngleafmas_ns_in(i,j) = 0.0
        ngleafmas_s_in(i,j) = 0.0
        nstemmass_in(i,j) = 0.0
        nstemmass_ns_in(i,j) = 0.0
        nstemmass_s_in(i,j) = 0.0
        nrootmass_in(i,j) = 0.0
        nrootmass_ns_in(i,j) = 0.0
        nrootmass_s_in(i,j) = 0.0
        nresorpedveg_s(i,j) = 0.0
        nresorpedveg_r(i,j) = 0.0
        nre_allocveg_sr2l(i,j) = 0.0
        nre_allocveg_s2l(i,j) = 0.0
        nre_allocveg_r2l(i,j) = 0.0
        nstress(i,j) = 0.0
        bnf_nat(i,j) = 0.0
        bnf_ant(i,j) = 0.0
     end do ! loop 110
  end do ! loop 100


 ! Symbiotic BNF:
    do j = 1,icc
       pftkind = ctempfts(j)
       do i = il1,il2

         ! Calculate N stress:
         if(ndemandveg_wp_npp(i,j) > 0.0) then
           nstress(i,j) = max(0.,min((ndemandveg_wp_npp(i,j) - (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j)))/ &
                            ndemandveg_wp_npp(i,j),1.))
         end if

         ! Calculate symbiotic BNF:
         select case (pftkind)

         case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr','GrassC3 ','GrassC4 ','Sedge   ','BdlEvgSh','BdlDCoSh')  
           if (ntchrveg(i,j) > 0.0 .and. isand(i,1) /= -4 .and. isand(i,1) /= -3) then
             if ((tbar(i,1)-273.16)<44.83 .and. (tbar(i,1)-273.16)>1.41) then
               bnf_nat(i,j) = max(0.,(r_bnf_s(sort(j)) * (nstress(i,j) + b_bnf(sort(j))))&
                                     * max(0.,(44.83-(tbar(i,1)-273.16))/(44.83-32.71)*(((tbar(i,1)-273.16)-1.41)/(32.71-1.41))**((32.71-1.41)/(44.83-32.71))))
               if(C_cost_bnf_s > 0.) bnf_nat(i,j) = min(bnf_nat(i,j),max_frac_bnf(sort(j))*1000*rootmass_ns(i,j)/C_cost_bnf_s)
             else
               bnf_nat(i,j) = 0.
             end if
             bnf_ant(i,j) = 0.
           end if

         case ('CropC3  ','CropC4  ')
           if (ntchrveg(i,j) > 0.0 .and. isand(i,1) /= -4 .and. isand(i,1) /= -3) then
             bnf_nat(i,j) = 0.
             if ((tbar(i,1)-273.16)<44.83 .and. (tbar(i,1)-273.16)>1.41) then
               bnf_ant(i,j) = max(0.,(r_bnf_s(sort(j)) * nstress(i,j) + b_bnf(sort(j)))) !&
               if(C_cost_bnf_s > 0.) bnf_ant(i,j) = min(bnf_ant(i,j),max_frac_bnf(sort(j))*1000*rootmass_ns(i,j)/C_cost_bnf_s)
             else
               bnf_ant(i,j) = 0.
             end if
           endif

         case default
                print * ,'Unknown CTEM PFT in nallocate ',pftkind
                call errorHandler('nallocate', - 2)
         end select

      end do
    end do


  ! Nitrogen allocation to tissues:
  do j = 1,icc ! loop 120
     pftkind = ctempfts(j)
     do i = il1,il2 ! loop 130

        !> 1) NITROGEN ALLOCATION TO ALL TISSUES: N uptake is distributed to leaf, stem, and root:
        if (fcancmx(i,j) > 0.0) then

          !! A) in proportion to their N demand (NPP > 0): <br>
          !! (see equations 4 and 12 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147) <br>
          !! \f$ a_i = \frac{\Delta_i}{\Delta_{WP}} \f$ <br>
          !! \f$ A_{R2L} = a_L(U_{NH_4}+U_{NO_3}) \f$ <br>
          !! \f$ A_{R2S} = a_S(U_{NH_4}+U_{NO_3}) \f$ <br>
          if (ndemandveg_wp_npp(i,j) > 0.0) then

             select case (pftkind)             

             case ('NdlDcdTr','BdlDCoTr','BdlDDrTr','BdlDCoSh')

                if (lfstatus(i,j) == 3 .or. lfstatus(i,j) == 4) then

                  nallocveg_l(i,j) = 0.0
                  temp_ndemand = ndemandveg_s_npp(i,j) + ndemandveg_r_npp(i,j)
                  if (temp_ndemand > 0.0) then
                     nallocveg_s(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                        * ndemandveg_s_npp(i,j) / temp_ndemand
                     nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                        * ndemandveg_r_npp(i,j) / temp_ndemand
                  else
                       nallocveg_s(i,j) = 0.0
                       nallocveg_r(i,j) = 0.0
                  end if
                  temp_ndemand = 0.0

                else
                  nallocveg_l(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                     * (ndemandveg_l_npp(i,j) / ndemandveg_wp_npp(i,j))
                  nallocveg_s(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                     * (ndemandveg_s_npp(i,j) / ndemandveg_wp_npp(i,j))
                  nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                     * (ndemandveg_r_npp(i,j) / ndemandveg_wp_npp(i,j))
                end if

             case ('NdlEvgTr','BdlEvgTr','CropC3  ','CropC4  ','BdlEvgSh')
                nallocveg_l(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                   * (ndemandveg_l_npp(i,j) / ndemandveg_wp_npp(i,j))
                nallocveg_s(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                   * (ndemandveg_s_npp(i,j) / ndemandveg_wp_npp(i,j))
                nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                   * (ndemandveg_r_npp(i,j) / ndemandveg_wp_npp(i,j))

              case ('GrassC3 ','GrassC4 ','Sedge   ')
                nallocveg_l(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                   * (ndemandveg_l_npp(i,j) / ndemandveg_wp_npp(i,j))
                nallocveg_s(i,j) = 0.
                nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) &
                                   * (ndemandveg_r_npp(i,j) / ndemandveg_wp_npp(i,j))
              case default
                 print * ,'Unknown CTEM PFT in nallocate ',pftkind
                 call errorHandler('nallocate', - 1)
              end select

          !! B) in proportion to their minimum C:N ratio (NPP < 0): <br>
          !! (see equations 4 and 13 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147) <br>
          !! \f$ a_i = frac{1/C:N_{i,min}}{1/C:N_{L,min}+1/C:N_{S,min}+1/C:N_{R,min}} \f$ <br>
          !! \f$ A_{R2L} = a_L(U_{NH_4}+U_{NO_3}) \f$ <br>
          !! \f$ A_{R2S} = a_S(U_{NH_4}+U_{NO_3}) \f$ <br>
          else

             select case (pftkind)

             case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
                   'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')

                   if (gleafmas(i,j) < zero .or. lfstatus(i,j) == 3 .or. lfstatus(i,j) == 4) then
                      nallocveg_l(i,j) = 0.0
                      temp_cn = (1. / c2n_smin(sort(j))) + (1. / c2n_rmin(sort(j)))
                      nallocveg_s(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) & 
                                         * (1. / c2n_smin(sort(j))) / temp_cn
                      nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) & 
                                         * (1. / c2n_rmin(sort(j))) / temp_cn
                   else
                      temp_cn = (1. / c2n_lmin(sort(j))) + (1. / c2n_smin(sort(j))) &
                        + (1. / c2n_rmin(sort(j)))
                      nallocveg_l(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) & 
                                         * (1. / c2n_lmin(sort(j))) / temp_cn
                      nallocveg_s(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) & 
                                         * (1. / c2n_smin(sort(j))) / temp_cn
                      nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) & 
                                         * (1. / c2n_rmin(sort(j))) / temp_cn
                   end if

             case ('GrassC3 ','GrassC4 ','Sedge   ')

                   temp_cn = (1. / c2n_lmin(sort(j))) + (1. / c2n_rmin(sort(j)))
                   if (gleafmas(i,j) <  zero .or. lfstatus(i,j) == 3 .or. lfstatus(i,j) == 4) then
                      nallocveg_l(i,j) = 0.0
                      nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) 
                   else
                      nallocveg_l(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) & 
                                         * (1. / c2n_lmin(sort(j))) / temp_cn
                      nallocveg_r(i,j) = (nuptakeveg_nh4(i,j) + nuptakeveg_no3(i,j) + bnf_nat(i,j) + bnf_ant(i,j)) & 
                                         * (1. / c2n_rmin(sort(j))) / temp_cn
                   end if

             case default
                print * ,'Unknown CTEM PFT in nallocate ',pftkind
                call errorHandler('nallocate', - 2)
             end select
             temp_cn = 0.0

          end if

        end if !fcancmx


     !> 2) REALLOCATION OF NITROGEN FROM LEAF TO STEM AND ROOT DURING LEAF OFFSET
     !! (see equation 14 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
     !! \f$ R_{L2R} = r_LLF_L\frac{N_{R,NS}}{N_{R,NS}+N_{S,NS}} \f$
     !! \f$ R_{L2S} = r_LLF_L\frac{N_{S,NS}}{N_{R,NS}+N_{S,NS}} \f$
     if (tltrleaf(i,j) > 0.0 .and. & !lfstatus(i,j) == 3 .and. ! leaf status is almost never 3 for PFT4 - this needs to be looked at
         gleafmas(i,j) > 0.0 .and. ngleafmas(i,j) > 0.0) then

        select case (pftkind)

        case ('NdlDcdTr','BdlDCoTr','BdlDDrTr','BdlDCoSh')

            gleafltr_fr = min((deltat/963.62) * tltrleaf(i,j) / ((deltat/963.62) * tltrleaf(i,j) + gleafmas(i,j)), 1.0)

            ngleafmas_resorped = nresorp_coeff(sort(j)) * gleafltr_fr &
                               * ngleafmas(i,j)

            !! Resorped N from leaves to be allocated to stem and roots
            !! in proportion to their existing pool sizes.
            !! For now, we use the following condition to allow resorption
            !! only after summer solstice in each hemisphere (indicated by radj, latitude in radians) 
            !! and before winter solstice.
            !! This avoids N resorption on occasions in the spring season when lfstatus becomes 3.

            temp_reallocation_offset = nstemmass_ns(i,j) + nrootmass_ns(i,j)
            if ( (radj(i) > 0. .and. iday >= 172) .or. &
                              (radj(i) < 0. .and. (iday >= 356 .or. iday <= 172))) then
               if (temp_reallocation_offset > 0.0) then
                  nresorpedveg_s(i,j) = ngleafmas_resorped * nstemmass_ns(i,j) / temp_reallocation_offset
                  nresorpedveg_r(i,j) = ngleafmas_resorped * nrootmass_ns(i,j) / temp_reallocation_offset
                  temp_reallocation_offset = 0.0
               else if (temp_reallocation_offset == 0.0) then
                  nresorpedveg_s(i,j) = ngleafmas_resorped * 0.75
                  nresorpedveg_r(i,j) = ngleafmas_resorped * 0.25
               end if
            end if

        case ('NdlEvgTr','BdlEvgTr','CropC3  ','CropC4  ',&
              'GrassC3 ','GrassC4 ','Sedge   ','BdlEvgSh')
           ! do nothing
        case default
           print * ,'Unknown CTEM PFT in nallocate ',pftkind
           call errorHandler('nallocate', - 4)
        end select

     end if

     end do ! for i, loop 130
  end do ! for j, loop 120

  !> 3) REALLOCATION OF NITROGEN FROM STEM AND ROOT TO LEAVES DURING LEAF ONSET
  !! (see equation 15 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ R_{R2L} = \frac{R_{R2L,C}}{C:N_L} \frac{N_{R,NS}}{N_{R,NS}+N_{S,NS}} \f$
  !! \f$ R_{S2L} = \frac{R_{S2L,C}}{C:N_L} \frac{N_{S,NS}}{N_{R,NS}+N_{S,NS}} \f$
  do j = 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2

        if (&!lfstatus(i,j) == 1 .and. ailcg(i,j) < lfthrs(i,j) .and. & ! removed because re_alloc_sr2l already takes into account leaf status
            re_alloc_sr2l(i,j) > zero .and. gleafmas(i,j) > zero) then

            nre_allocveg_sr2l(i,j) = re_alloc_sr2l(i,j) / c2n_lmin(sort(j))

            select case (pftkind)

            case ('NdlDcdTr','BdlDCoTr','BdlDDrTr','BdlDCoSh') ! for dcd trees

              if (nstemmass(i,j) > 0.0 .and. nrootmass(i,j) > 0.0) then
                  nstemns2t = nstemmass_ns(i,j) / nstemmass(i,j)
                  nrootns2t = nrootmass_ns(i,j) / nrootmass(i,j)

                  !if (nstemns2t < 1.0 .and. nrootns2t < 1.0) then
                     temp_reallocation_onset = nstemmass_ns(i,j) + nrootmass_ns(i,j)
                     if (temp_reallocation_onset > 0.0) then
                        ns_s2sr = nstemmass_ns(i,j) / temp_reallocation_onset
                        ns_r2sr = nrootmass_ns(i,j) / temp_reallocation_onset
                        temp_reallocation_onset = 0.0
                     end if

                      if (nstemns2t > min_ns2t_s(sort(j)) &
                          .and. nrootns2t > min_ns2t_r(sort(j)) &
                          .and. nre_allocveg_sr2l(i,j) * ns_s2sr <= nstemmass_ns(i,j) &
                          .and. nre_allocveg_sr2l(i,j) * ns_r2sr <= nrootmass_ns(i,j)) then

                          nre_allocveg_s2l(i,j) = nre_allocveg_sr2l(i,j) * ns_s2sr
                          nre_allocveg_r2l(i,j) = nre_allocveg_sr2l(i,j) * ns_r2sr

                      else if (nstemns2t > min_ns2t_s(sort(j)) &
                          .and. nrootns2t <= min_ns2t_r(sort(j)) &
                          .and. nre_allocveg_sr2l(i,j) <= nstemmass_ns(i,j)) then

                          nre_allocveg_s2l(i,j) = nre_allocveg_sr2l(i,j)
                          nre_allocveg_r2l(i,j) = 0.0

                      else if (nstemns2t <= min_ns2t_s(sort(j)) &
                          .and. nrootns2t > min_ns2t_r(sort(j)) &
                          .and. nre_allocveg_sr2l(i,j) <= nrootmass_ns(i,j)) then

                          nre_allocveg_s2l(i,j) = 0.0
                          nre_allocveg_r2l(i,j) = nre_allocveg_sr2l(i,j)

                      else if (nstemns2t <= min_ns2t_s(sort(j)) &
                          .and. nrootns2t <= min_ns2t_r(sort(j))) then

                          nre_allocveg_s2l(i,j) = 0.0
                          nre_allocveg_r2l(i,j) = 0.0

                      end if

                      nre_allocveg_sr2l(i,j) = nre_allocveg_s2l(i,j) + nre_allocveg_r2l(i,j)

                  !end if
              end if

            case ('GrassC3 ','GrassC4 ','Sedge   ') ! for grasses

              if (nrootmass(i,j) > 0.0) then
                  nrootns2t = nrootmass_ns(i,j) / nrootmass(i,j)

                  !if (nrootns2t /= 1.0) then

                     if (nrootns2t > min_ns2t_r(sort(j)) .and. &
                         nre_allocveg_sr2l(i,j) < nrootmass_ns(i,j)) then
                         nre_allocveg_s2l(i,j) = 0.0
                         nre_allocveg_r2l(i,j) = nre_allocveg_sr2l(i,j)
                     else
                         nre_allocveg_s2l(i,j) = 0.0
                         nre_allocveg_r2l(i,j) = 0.0
                     end if
                     nre_allocveg_sr2l(i,j) = nre_allocveg_s2l(i,j) + nre_allocveg_r2l(i,j)

                  !end if

              end if

            case ('NdlEvgTr','BdlEvgTr','CropC3  ','CropC4  ','BdlEvgSh')
              ! do nothing

            case default
              print * ,'Unknown CTEM PFT in nallocate ',pftkind
              call errorHandler('nallocate', - 5)

            end select
        end if

     end do ! for i, loop 150
  end do ! for j, loop 140


  !> Calculate nitrogen net change with the above 5 allocation types  
  do j = 1,icc
     do i = il1,il2

        nntchveg_l(i,j) = nallocveg_l(i,j) + nre_allocveg_s2l(i,j) &
                          + nre_allocveg_r2l(i,j) - nresorpedveg_s(i,j) &
                          - nresorpedveg_r(i,j)
        nntchveg_s(i,j) = nallocveg_s(i,j) - nre_allocveg_s2l(i,j) &
                          + nresorpedveg_s(i,j)
        nntchveg_r(i,j) = nallocveg_r(i,j) - nre_allocveg_r2l(i,j) &
                          + nresorpedveg_r(i,j)

    end do ! loop 190
  end do ! loop 180



  !> Update vegetation pools and checking conservation

  !! 1) N GREEN LEAF MASS
  do j = 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2
       select case (pftkind)

       case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
             'CropC3  ','CropC4  ','GrassC3 ','GrassC4 ','Sedge   ','BdlEvgSh','BdlDCoSh')
         
         ngleafmas_in(i,j) = ngleafmas(i,j)
         ngleafmas_ns_in(i,j) = ngleafmas_ns(i,j)
         ngleafmas_s_in(i,j) = ngleafmas_s(i,j)

         !! Nitrogen flux from non-structural to structural leaf nitrogen
         if (nntchveg_l(i,j) > 0.0 .and. ngleafmas(i,j) > 0.0) then
             nleafns2sveg(i,j) = Fmax_leaf_N(sort(j)) * nntchveg_l(i,j) * &
                                 max(0., ngleafmas_ns(i,j) / ngleafmas(i,j) - min_ns2t_l(sort(j)))
         end if

         !! A) Preventing negative ngleafmas_ns(i,j) when nntchveg_l(i,j) is negative
         if ((ngleafmas_ns(i,j) + nntchveg_l(i,j)) < 0.0 .and. ngleafmas(i,j) > 0.0) then

          if ( (ngleafmas_ns(i,j) + nntchveg_l(i,j) * ngleafmas_ns(i,j) / ngleafmas(i,j)) >= 0.0 .and. &
               (ngleafmas_s(i,j) + nntchveg_l(i,j) * ngleafmas_s(i,j) / ngleafmas(i,j)) >= 0.0 ) then

             temp_ns(i,j) = nntchveg_l(i,j) * (ngleafmas_ns(i,j) / ngleafmas(i,j))
             temp_s(i,j) = nntchveg_l(i,j) * (ngleafmas_s(i,j) / ngleafmas(i,j))

          else
              if(nresorpedveg_s(i,j) + nresorpedveg_r(i,j)>0.0) then
                nresorpedveg_s(i,j) = ngleafmas(i,j)*nresorpedveg_s(i,j)/(nresorpedveg_s(i,j) + nresorpedveg_r(i,j))
                nresorpedveg_r(i,j) = ngleafmas(i,j)*nresorpedveg_r(i,j)/(nresorpedveg_s(i,j) + nresorpedveg_r(i,j))
              else
                nresorpedveg_s(i,j) = ngleafmas(i,j)*0.75
                nresorpedveg_r(i,j) = ngleafmas(i,j)*0.25
              end if

              nntchveg_l(i,j) = nallocveg_l(i,j) + nre_allocveg_s2l(i,j) &
                          + nre_allocveg_r2l(i,j) - nresorpedveg_s(i,j) &
                          - nresorpedveg_r(i,j)
              nntchveg_s(i,j) = nallocveg_s(i,j) - nre_allocveg_s2l(i,j) &
                          + nresorpedveg_s(i,j)
              nntchveg_r(i,j) = nallocveg_r(i,j) - nre_allocveg_r2l(i,j) &
                          + nresorpedveg_r(i,j)
              temp_ns(i,j) = nntchveg_l(i,j) * (ngleafmas_ns(i,j) / ngleafmas(i,j))
              temp_s(i,j) = nntchveg_l(i,j) * (ngleafmas_s(i,j) / ngleafmas(i,j))
          end if

         !! B) Preventing negative ngleafmas_ns(i,j) when nntchveg_l(i,j) < nleafns2sveg(i,j)
         else if (ngleafmas_ns(i,j) + nntchveg_l(i,j) < nleafns2sveg(i,j)) then

            nleafns2sveg(i,j) = nntchveg_l(i,j)
            temp_ns(i,j) = nntchveg_l(i,j)
            temp_s(i,j) = 0.0

         !! E) normal situation
         else

            temp_ns(i,j) = nntchveg_l(i,j)
            temp_s(i,j) = 0.0

         end if

       case default
         print * ,'Unknown CTEM PFT in nallocate ',pftkind
         call errorHandler('nallocate', - 7)
       end select

       !! Updating ngleafmas pools
       ngleafmas_ns(i,j) = ngleafmas_ns(i,j) - nleafns2sveg(i,j) + temp_ns(i,j)
       ngleafmas_s(i,j) = ngleafmas_s(i,j) + nleafns2sveg(i,j) + temp_s(i,j)
       ngleafmas(i,j) = ngleafmas_ns(i,j) + ngleafmas_s(i,j)

       if (ngleafmas_ns(i,j) < 0.0 .or. ngleafmas_s(i,j) < 0.0 .or. &
           ngleafmas(i,j) < 0.0) then
           write(*,*)'ngleafmas lt zero in nallocate.f90 at i=',i,' for pft=',j
           !write(*,*)'ngleafmas_ns = ',ngleafmas_ns(i,j)
           !write(*,*)'ngleafmas_s  = ',ngleafmas_s(i,j)
           !write(*,*)'ngleafmas    = ',ngleafmas(i,j)
           call errorHandler('nallocate', - 8)
       end if

       !! Checking conservation for ngleafmas(i,j)
       if ((ngleafmas_in(i,j) - ngleafmas(i,j) + nallocveg_l(i,j) &
            + nre_allocveg_s2l(i,j) + nre_allocveg_r2l(i,j) &
            - nresorpedveg_s(i,j) - nresorpedveg_r(i,j) &
            ) > ntolerance) then
            write(*,*)'imbalance in ngleafmas pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 9)
       end if
       if ((ngleafmas_ns_in(i,j) - ngleafmas_ns(i,j) + temp_ns(i,j) &
            - nleafns2sveg(i,j)) > ntolerance) then
            write(*,*)'ngleafmas_ns_in = ',ngleafmas_ns_in(i,j)
            write(*,*)'ngleafmas_ns = ',ngleafmas_ns(i,j)
            write(*,*)'temp_ns = ',temp_ns(i,j)
            write(*,*)'nleafns2sveg = ',nleafns2sveg(i,j)
            write(*,*)'imbalance in ngleafmas_ns pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 10)
       end if
       if ((ngleafmas_s_in(i,j) - ngleafmas_s(i,j) + temp_s(i,j) &
            + nleafns2sveg(i,j)) > ntolerance) then
            write(*,*)'ngleafmas_s_in = ',ngleafmas_s_in(i,j)
            write(*,*)'ngleafmas_s = ',ngleafmas_s(i,j)
            write(*,*)'temp_s = ',temp_s(i,j)
            write(*,*)'nleafns2sveg = ',nleafns2sveg(i,j)
            write(*,*)'imbalance in ngleafmas_s pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 11)
       end if

       temp_ns(i,j) = 0.0
       temp_s(i,j) = 0.0

     end do ! loop 210
  end do ! loop 200

  !> 2) N STEM MASS

  do j = 1,icc
     do i = il1,il2

       nstemmass_in(i,j) = nstemmass(i,j) 
       nstemmass_ns_in(i,j) = nstemmass_ns(i,j)
       nstemmass_s_in(i,j) = nstemmass_s(i,j)

       !! Nitrogen flux from non-structural to structural stem nitrogen
       if (nntchveg_s(i,j) > 0.0 .and. nstemmass(i,j) > 0.0) then 
          nstemns2sveg(i,j) = Fmax_stem_N(sort(j)) * nntchveg_s(i,j) * &
                              max(0.,nstemmass_ns(i,j) / nstemmass(i,j) - min_ns2t_s(sort(j)))
       end if

       !! A) Preventing negative nstemmass_ns(i,j) resulting from negative nntchveg_s(i,j)
       if ((nstemmass_ns(i,j) + nntchveg_s(i,j)) < 0.0 .and. nstemmass(i,j) > 0.0)  then

          if ( (nstemmass_ns(i,j) + nntchveg_s(i,j) * nstemmass_ns(i,j) / nstemmass(i,j)) >= 0.0 .and. &
               (nstemmass_s(i,j) + nntchveg_s(i,j) * nstemmass_s(i,j) / nstemmass(i,j)) >= 0.0 ) then

              temp_ns(i,j) = nntchveg_s(i,j) * (nstemmass_ns(i,j) / nstemmass(i,j))
              temp_s(i,j) = nntchveg_s(i,j) * (nstemmass_s(i,j) / nstemmass(i,j))

          else
              nre_allocveg_s2l(i,j) = nstemmass(i,j)
              nre_allocveg_sr2l(i,j) = nre_allocveg_s2l(i,j) + nre_allocveg_r2l(i,j)
              nntchveg_l(i,j) = nallocveg_l(i,j) + nre_allocveg_s2l(i,j) &
                          + nre_allocveg_r2l(i,j) - nresorpedveg_s(i,j) &
                          - nresorpedveg_r(i,j)
              nntchveg_s(i,j) = nallocveg_s(i,j) - nre_allocveg_s2l(i,j) &
                          + nresorpedveg_s(i,j)
              nntchveg_r(i,j) = nallocveg_r(i,j) - nre_allocveg_r2l(i,j) &
                          + nresorpedveg_r(i,j)
              nntchveg_s(i,j) = nallocveg_s(i,j) - nre_allocveg_s2l(i,j) + nresorpedveg_s(i,j)
              temp_ns(i,j) = nntchveg_s(i,j) * (nstemmass_ns(i,j) / nstemmass(i,j))
              temp_s(i,j) = nntchveg_s(i,j) * (nstemmass_s(i,j) / nstemmass(i,j))
          endif

       !! B) Preventing negative nstemmass_ns(i,j) when nntchveg_s(i,j) < nstemns2sveg(i,j)
       else if ((nstemmass_ns(i,j) + nntchveg_s(i,j)) < nstemns2sveg(i,j)) then

          nstemns2sveg(i,j) = nntchveg_s(i,j)
          temp_ns(i,j) = nntchveg_s(i,j)
          temp_s(i,j) = 0.0

       !! C) normal situation
       else

          temp_ns(i,j) = nntchveg_s(i,j)
          temp_s(i,j) = 0.0

       end if

       nstemmass_ns(i,j) = nstemmass_ns(i,j) - nstemns2sveg(i,j) + temp_ns(i,j)
       nstemmass_s(i,j) = nstemmass_s(i,j)  + nstemns2sveg(i,j) + temp_s(i,j)
       nstemmass(i,j) = nstemmass_ns(i,j) + nstemmass_s(i,j)

       if (nstemmass_ns(i,j) < 0.0 .and. abs(nstemmass_ns(i,j)) < 1.E-10) then
          nstemmass_ns(i,j) = 0.0
          nstemmass(i,j) = nstemmass_s(i,j)
       end if
       if (nstemmass_s(i,j) < 0.0 .and. abs(nstemmass_s(i,j)) < 1.E-10) then
          nstemmass_s(i,j) = 0.0
          nstemmass(i,j) = nstemmass_ns(i,j)
       end if

       if (nstemmass_ns(i,j) < 0.0 .or. nstemmass_s(i,j) < 0.0 .or. nstemmass(i,j) < 0.0) then
           write(*,*)'nstemmass_ns lt zero in nallocate.f90 at i=',i,' for pft=',j
           ! write(*,*)'nstemmass_ns = ',nstemmass_ns(i,j)
           ! write(*,*)'nstemmass_s  = ',nstemmass_s(i,j)
           ! write(*,*)'nstemmass    = ',nstemmass(i,j)
           call errorHandler('nallocate', - 12)
       end if

       !! Checking conservation for nstemmass(i,j)
       if ((nstemmass_in(i,j) - nstemmass(i,j) + nallocveg_s(i,j) &
            - nre_allocveg_s2l(i,j) + nresorpedveg_s(i,j) &
            ) > ntolerance) then
            write(*,*)'imbalance in nstemmass pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 13)
       end if
       if ((nstemmass_ns_in(i,j) - nstemmass_ns(i,j) + temp_ns(i,j) &
            - nstemns2sveg(i,j)) > ntolerance) then
            write(*,*)'nstemmass_ns_in = ',nstemmass_ns_in(i,j)
            write(*,*)'nstemmass_ns = ',nstemmass_ns(i,j)
            write(*,*)'temp_ns = ',temp_ns(i,j)
            write(*,*)'nstemns2sveg = ',nstemns2sveg(i,j)
            write(*,*)'imbalance in nstemmass_ns pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 10)
       end if
       if ((nstemmass_s_in(i,j) - nstemmass_s(i,j) + temp_s(i,j) &
            + nstemns2sveg(i,j)) > ntolerance) then
            write(*,*)'nstemmass_s_in = ',nstemmass_s_in(i,j)
            write(*,*)'nstemmass_s = ',nstemmass_s(i,j)
            write(*,*)'temp_s = ',temp_s(i,j)
            write(*,*)'nstemns2sveg = ',nstemns2sveg(i,j)
            write(*,*)'imbalance in nstemmass_s pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 11)
       end if

       temp_ns(i,j) = 0.0
       temp_s(i,j) = 0.0

     end do ! loop 230
  end do ! loop 220

  !> 3) N ROOT MASS

  do j = 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2

       nrootmass_in(i,j) = nrootmass(i,j)
       nrootmass_ns_in(i,j) = nrootmass_ns(i,j)
       nrootmass_s_in(i,j) = nrootmass_s(i,j)

       !! Nitrogen flux from non-structural to structural root nitrogen
       if (nntchveg_r(i,j) > 0.0 .and. nrootmass(i,j) > 0.0) then
           nrootns2sveg(i,j) = Fmax_root_N(sort(j)) * nntchveg_r(i,j) * &
                               max(0., nrootmass_ns(i,j) / nrootmass(i,j) - min_ns2t_r(sort(j)))
       end if

       !! A) Preventing negative nrootmass_ns(i,j) resulting from negative nntchveg_r(i,j) 
       if ((nrootmass_ns(i,j) + nntchveg_r(i,j)) < 0.0 .and. nrootmass(i,j) > 0.0)  then

          if ((nrootmass_ns(i,j) + nntchveg_r(i,j) * nrootmass_ns(i,j) / nrootmass(i,j)) >= 0.0 .and. &
              (nrootmass_s(i,j) + nntchveg_r(i,j) * nrootmass_s(i,j) / nrootmass(i,j)) >= 0.0 ) then
              temp_ns(i,j) = nntchveg_r(i,j) * (nrootmass_ns(i,j) / nrootmass(i,j))
              temp_s(i,j) = nntchveg_r(i,j) * (nrootmass_s(i,j) / nrootmass(i,j))
          else
              nre_allocveg_r2l(i,j) = nrootmass(i,j)
              nre_allocveg_sr2l(i,j) = nre_allocveg_s2l(i,j) + nre_allocveg_r2l(i,j)
              nntchveg_l(i,j) = nallocveg_l(i,j) + nre_allocveg_s2l(i,j) &
                          + nre_allocveg_r2l(i,j) - nresorpedveg_s(i,j) &
                          - nresorpedveg_r(i,j)
              nntchveg_s(i,j) = nallocveg_s(i,j) - nre_allocveg_s2l(i,j) &
                          + nresorpedveg_s(i,j)
              nntchveg_r(i,j) = nallocveg_r(i,j) - nre_allocveg_r2l(i,j) &
                          + nresorpedveg_r(i,j)
              temp_ns(i,j) = nntchveg_r(i,j) * (nrootmass_ns(i,j) / nrootmass(i,j))
              temp_s(i,j) = nntchveg_r(i,j) * (nrootmass_s(i,j) / nrootmass(i,j))
          end if

       !! B) Preventing negative nrootmass_ns(i,j) when nntchveg_r(i,j) < nrootns2sveg(i,j)
       else if ((nrootmass_ns(i,j) + nntchveg_r(i,j)) < nrootns2sveg(i,j)) then

          nrootns2sveg(i,j) = nntchveg_r(i,j)
          temp_ns(i,j) = nntchveg_r(i,j)
          temp_s(i,j) = 0.0

       !! C) normal situation
       else

        temp_ns(i,j) = nntchveg_r(i,j)
        temp_s(i,j) = 0.0

       end if

       nrootmass_ns(i,j) = nrootmass_ns(i,j) - nrootns2sveg(i,j) + temp_ns(i,j)
       nrootmass_s(i,j) = nrootmass_s(i,j) + nrootns2sveg(i,j) + temp_s(i,j)
       nrootmass(i,j) = nrootmass_ns(i,j) + nrootmass_s(i,j)

      if (nrootmass_ns(i,j) < 0.0 .and. abs(nrootmass_ns(i,j)) < 1.E-10) then
         nrootmass_ns(i,j) = 0.0
         nrootmass(i,j) = nrootmass_s(i,j)
      end if
      if (nrootmass_s(i,j) < 0.0 .and. abs(nrootmass_s(i,j)) < 1.E-10) then
         nrootmass_s(i,j) = 0.0
         nrootmass(i,j) = nrootmass_ns(i,j)
      end if

      if (nrootmass_ns(i,j) < 0.0 .or. nrootmass_s(i,j) < 0.0 .or. &
          nrootmass(i,j) < 0.0) then
         write(*,*)'nrootmass lt zero in nallocate.f90 at i=',i,' for pft=',j
         ! write(*,*)'nrootmass_ns(i,j)= ',nrootmass_ns(i,j)
         ! write(*,*)'nrootmass_s(i,j) = ',nrootmass_s(i,j)
         ! write(*,*)'nrootmass(i,j)   = ',nrootmass(i,j)
         call errorHandler('nallocate', - 14)
     end if

     !! Checking conservation for nrootmass(i,j)
     if ((nrootmass_in(i,j) - nrootmass(i,j) + nallocveg_r(i,j) &
         - nre_allocveg_r2l(i,j) + nresorpedveg_r(i,j) &
         ) > ntolerance) then
          write(*,*)'imbalance in nrootmass pool in nallocate.f90 subroutine at (i)=',i,'pft=',j
          call errorHandler('nallocate', - 15)
     end if
     if ((nrootmass_ns_in(i,j) - nrootmass_ns(i,j) + temp_ns(i,j) &
            - nrootns2sveg(i,j)) > ntolerance) then
            write(*,*)'nrootmass_ns_in = ',nrootmass_ns_in(i,j)
            write(*,*)'nrootmass_ns = ',nrootmass_ns(i,j)
            write(*,*)'temp_ns = ',temp_ns(i,j)
            write(*,*)'nrootns2sveg = ',nrootns2sveg(i,j)
            write(*,*)'imbalance in nrootmass_ns pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 10)
     end if
     if ((nrootmass_s_in(i,j) - nrootmass_s(i,j) + temp_s(i,j) &
            + nrootns2sveg(i,j)) > ntolerance) then
            write(*,*)'nrootmass_s_in = ',nrootmass_s_in(i,j)
            write(*,*)'nrootmass_s = ',nrootmass_s(i,j)
            write(*,*)'temp_s = ',temp_s(i,j)
            write(*,*)'nrootns2sveg = ',nrootns2sveg(i,j)
            write(*,*)'imbalance in nrootmass_s pool in nallocate subroutine at (i)=',i,'pft=',j
            call errorHandler('nallocate', - 11)
     end if

     temp_ns(i,j) = 0.0
     temp_s(i,j) = 0.0

     end do ! loop 250
  end do ! loop 240


  return
end subroutine nallocate
