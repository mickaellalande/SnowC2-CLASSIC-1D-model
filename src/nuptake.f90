!> \file
!! Calculates plant nitrogen demand and uptake.
!> @author Ali Asaadi
!!


subroutine nuptake(il1, il2, sort, fcancmx, zbotw, delzw, & ! In
                   thliq, rootdpth, QFC, ntchlveg, ntchsveg, ntchrveg, gleafmas, & ! In
                   stemmass, rootmass, ngleafmas, nstemmass, lfstatus, iday, & ! In
                   nrootmass, nh4_mass, no3_mass, & ! In/Out
                   ndemandveg_l_npp, ndemandveg_s_npp, ndemandveg_r_npp, & !Out
                   ndemandveg_wp_npp, nuptakeveg_p_nh4, & !Out
                   nuptakeveg_p_no3, nuptakeveg_a_actl_nh4, nuptakeveg_a_actl_no3, & !Out
                   nuptakeveg_nh4, nuptakeveg_no3) !Out


  use classicParams, only : ilg, ignd, icc, iccp1, ican, ctempfts, ntolerance, &
                            coeff_a, b_nh4, b_no3, kphalf, &
                            c2n_lmin, c2n_smin, c2n_rmin, coeff_p_nh4, &
                            coeff_p_no3, zero, classpfts, reindexPFTs, deltat

  implicit none

  integer,intent(in) :: il1             !< il1=1
  integer,intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in) :: iday            !< day of year
  integer,intent(in),dimension(icc) :: sort             !< index for correspondence between PFTs and the 12 values in parameters vectors
  integer,intent(in),dimension(ilg,icc) :: lfstatus     !< integer indicating leaf status or mode (max. growth=1, normal growth=2, fall=3, no leaves=4)
  real,intent(in),dimension(ilg,icc) :: ntchlveg        !< C allocation to leaf, \f$umol CO2 m^{-2} sec^{-1}\f$ 
  real,intent(in),dimension(ilg,icc) :: ntchsveg        !< C allocation to stem, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: ntchrveg        !< C allocation to root, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc) :: fcancmx         !< max. fractional coverage of PFTs 
  real,intent(in),dimension(ilg,ignd) :: delzw          !< thicknesses of each soil layer, \f$m\f$
  real,intent(in),dimension(ilg,ignd) :: thliq          !< volumetric soil moisture, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg,ignd) :: zbotw          !< bottom of each soil layer, \f$m\f$
  real,intent(in),dimension(ilg,ignd) :: QFC            !< water removed from each soil layer by transpiration, \f$kg H2O m^{-2} s^{-1}\f$ 
  real,intent(in),dimension(ilg,icc) :: rootdpth        !< rooting depth (contains 99% of rootmass), \f$m\f$
  real,intent(in),dimension(ilg,icc) :: gleafmas        !< green leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc) :: stemmass        !< stem carbon for each PFT, \f$kg C m^{-2}\f$ 
  real,intent(in),dimension(ilg,icc) :: rootmass        !< root carbon for each PFT,  \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc) :: ngleafmas       !< green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(in),dimension(ilg,icc) :: nstemmass       !< stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(in),dimension(ilg,icc) :: nrootmass       !< root nitrogen for each PFT, \f$g N m^{-2}\f$
  
  real,intent(inout),dimension(ilg,iccp1) :: nh4_mass   !< NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,iccp1) :: no3_mass   !< NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$

  real,intent(out),dimension(ilg,icc) :: ndemandveg_l_npp       !< leaf nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: ndemandveg_s_npp       !< stem nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: ndemandveg_r_npp       !< root nitrogen demand for each PFT based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: ndemandveg_wp_npp      !< whole plant nitrogen demand for each PFT  based on NPP, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nuptakeveg_nh4         !< total nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nuptakeveg_no3         !< total no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nuptakeveg_p_nh4       !< passive nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nuptakeveg_p_no3       !< passive no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nuptakeveg_a_actl_nh4  !< actual active nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nuptakeveg_a_actl_no3  !< actual active no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$

  ! local variables
  integer :: i, j, k, m
  character(8) :: pftkind
  real :: temp
  real :: temp_mineral                              !< temporary denominator of potential active nitrogen uptake
  real :: total_nuptake_p                           !< passive nitrogen uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real :: total_nuptake_a_ptnl                      !< potential active nitrogen uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: qfc_rz                 !< water removed from rooting zone by transpiration, \f$g H2O m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: rootzw                 !< water content in rooting zone, \f$g H2O m^{-2}\f$
  real,dimension(ilg,icc) :: nconc_rzw_nh4          !< ammonium concentration in rooting zone for each PFT, \f$g N g H2O^{-1}\f$ 
  real,dimension(ilg,icc) :: nconc_rzw_no3          !< nitrate concentration in rooting zone for each PFT, \f$g N g H2O^{-1}\f$
  real,dimension(ilg,iccp1) :: nh4_mass_in          !< previous NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: no3_mass_in          !< previous NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: finerootfr             !< fraction of fine roots in root carbon for each PFT 
  real,dimension(ilg,icc) :: ndemandveg_l           !< leaf nitrogen demand for each PFT based on plant carbon pools, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: ndemandveg_s           !< stem nitrogen demand for each PFT based on plant carbon pools, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: ndemandveg_r           !< root nitrogen demand for each PFT based on plant carbon pools, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nuptakeveg_a_ptnl_nh4  !< potential active nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_a_ptnl_nh4         !< grid avg. potential active nh4+ uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nuptakeveg_a_ptnl_no3  !< potential active no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg) :: nuptake_a_ptnl_no3         !< grid avg. potential active no3- uptake, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: ndemandveg_wp          !< whole plant nitrogen demand for each PFT  based on plant carbon pools, \f$g N m^{-2} day^{-1}\f$
  !
  temp = 0.0
  temp_mineral = 0.0
  total_nuptake_p = 0.0
  total_nuptake_a_ptnl = 0.0

  do i = il1,il2
    nuptake_a_ptnl_nh4(i) = 0.0
    nuptake_a_ptnl_no3(i) = 0.0
    do j = 1,icc
       ndemandveg_l(i,j) = 0.0
       ndemandveg_s(i,j) = 0.0
       ndemandveg_r(i,j) = 0.0
       ndemandveg_l_npp(i,j) = 0.0
       ndemandveg_s_npp(i,j) = 0.0
       ndemandveg_r_npp(i,j) = 0.0
       ndemandveg_wp(i,j) = 0.0
       ndemandveg_wp_npp(i,j) = 0.0
       nconc_rzw_nh4(i,j) = 0.0
       nconc_rzw_no3(i,j) = 0.0
       rootzw(i,j) = 0.0
       qfc_rz(i,j) = 0.0
       nuptakeveg_p_nh4(i,j) = 0.0
       nuptakeveg_p_no3(i,j) = 0.0
       nuptakeveg_a_ptnl_nh4(i,j) = 0.0
       nuptakeveg_a_ptnl_no3(i,j) = 0.0
       nuptakeveg_a_actl_nh4(i,j) = 0.0
       nuptakeveg_a_actl_no3(i,j) = 0.0
       nuptakeveg_nh4(i,j) = 0.0
       nuptakeveg_no3(i,j) = 0.0
       finerootfr(i,j) = 0.0
    end do ! loop 110
    do j = 1,iccp1
      no3_mass_in(i,j) = 0.0
      nh4_mass_in(i,j) = 0.0
    end do ! loop 120
  end do ! loop 100

  !> Calculate nitrogen demand for each PFT based on NPP:
  do j = 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2

        if (fcancmx(i,j) > 0.0) then

           if (ntchlveg(i,j) > 0.0) then
              ndemandveg_l_npp(i,j) = ntchlveg(i,j) * (deltat / 963.62) * 1000. &
                                      / c2n_lmin(sort(j))
           end if

           if (ntchsveg(i,j) > 0.0) then
              select case (pftkind)
              case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr',&
                    'BdlDDrTr','CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
                 ndemandveg_s_npp(i,j) = ntchsveg(i,j) * (deltat / 963.62) * 1000. &
                                       / c2n_smin(sort(j)) ! in (g N/m^2 day)
              case ('GrassC3 ','GrassC4 ','Sedge   ')
                 ! do nothing
              case default
                 print * ,'Unknown CTEM PFT in nuptake ',pftkind
                 call errorHandler('nuptake', - 1)
              end select
           end if

           if (ntchrveg(i,j) > 0.0) then
              ndemandveg_r_npp(i,j) = ntchrveg(i,j) * (deltat / 963.62) * 1000. &
                                      / c2n_rmin(sort(j))
           end if

           ndemandveg_wp_npp(i,j) = ndemandveg_l_npp(i,j) + ndemandveg_s_npp(i,j) &
                                    + ndemandveg_r_npp(i,j)
      end if ! fcancmx
     end do ! loop 140
  end do ! loop 130

  !> Calculate nitrogen demand for each PFT based on plant carbon pools:
  !! This is currently not used - ndemandveg_wp_npp (above) is used.
  do j = 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2

       if (fcancmx(i,j) > 0.0) then

         if (gleafmas(i,j) > 0.0 .and. lfstatus(i,j) /= 3 .and. lfstatus(i,j) /= 4) then
            ndemandveg_l(i,j) = gleafmas(i,j) * 1000. / c2n_lmin(sort(j))
            ndemandveg_l(i,j) = max(0.0, ndemandveg_l(i,j) - ngleafmas(i,j))
         end if

         if (stemmass(i,j) > 0.0 .and. ntchsveg(i,j) > 0.0 ) then
            ndemandveg_s(i,j) = stemmass(i,j) * 1000. / c2n_smin(sort(j))
            ndemandveg_s(i,j) = max(0.0, ndemandveg_s(i,j) - nstemmass(i,j))
         end if

         if (rootmass(i,j) > 0.0 .and. ntchrveg(i,j) > 0.0 ) then
            ndemandveg_r(i,j) = rootmass(i,j) * 1000. / c2n_rmin(sort(j))
            ndemandveg_r(i,j) = max(0.0, ndemandveg_r(i,j) - nrootmass(i,j))
         end if

         ndemandveg_wp(i,j) = ndemandveg_l(i,j) + ndemandveg_s(i,j) &
                              + ndemandveg_r(i,j)

       end if !fcancmx
     end do ! loop 160
  end do ! loop 150

  !> Calculate ammonium and nitrate concentrations in the rooting zone for each PFT
  !! and water removed from rooting zone by transpiration:
  !! (see equation 5 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ NH_4 = \frac{N_{NH_4}}{\sum_{i=1}^{i \leq r_I} 10^6 \theta_i z_i} \f$
  !! \f$ NO_3 = \frac{N_{NO_3}}{\sum_{i=1}^{i \leq r_I} 10^6 \theta_i z_i} \f$
  do j = 1,icc
    do i = il1,il2
      do k = 1, ignd

         if (zbotw(i,k) <= rootdpth(i,j)) then
            rootzw(i,j) = rootzw(i,j) + 1000. * 1000. * thliq(i,k) * delzw(i,k)
            qfc_rz(i,j) = qfc_rz(i,j) + 1000. * 86400. * qfc(i,k)
         end if

      end do ! loop 190

      if (rootzw(i,j) > 0.0) then
        nconc_rzw_nh4(i,j) = nh4_mass(i,j) / b_nh4 / rootzw(i,j)
        nconc_rzw_no3(i,j) = no3_mass(i,j) / b_no3 / rootzw(i,j)
      end if

    end do ! loop 180
  end do ! loop 170

  !> Calculate passive nitrogen uptake
  !! (see equation 6 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ U_{p,NH_4} = 86400 \times 10^3 \beta q_t [NH_4] \f$
  !! \f$ U_{p,NO_3} = 86400 \times 10^3 \beta q_t [NO_3] \f$
  do j = 1,icc
     do i = il1,il2

        if (fcancmx(i,j) > 0.0) then

           nuptakeveg_p_nh4(i,j) = 1.0 * qfc_rz(i,j) * coeff_p_nh4 &
                                   * nconc_rzw_nh4(i,j)
           nuptakeveg_p_no3(i,j) = 1.0 * qfc_rz(i,j) * coeff_p_no3 &
                                   * nconc_rzw_no3(i,j)

          if (nuptakeveg_p_nh4(i,j) < 0.0) then
             write(*,*)'Negative passive NH4+ uptake'
             write(*,*)'nuptakeveg_p_nh4 lt zero at i=',i,' for pft=',j
             write(*,*)'nuptakeveg_p_nh4 = ',nuptakeveg_p_nh4(i,j)
             call errorHandler('nuptake', - 2)
          end if

          if (nuptakeveg_p_no3(i,j) < 0.0) then
             write(*,*)'Negative passive NO3- uptake'
             write(*,*)'nuptakeveg_p_no3 lt zero at i=',i,' for pft=',j
             write(*,*)'nuptakeveg_p_no3 = ',nuptakeveg_p_no3(i,j)
             call errorHandler('nuptake', - 3)
          end if

        end if !fcancmx
     end do ! loop 210
  end do ! loop 200

  !> Calculate fraction of fine roots in root carbon:
  !! (see equation 7 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ f = 1-\frac{C_R}{C_R+0.6}\f$
  do j = 1,ican
     do m = reindexPFTs(j,1),reindexPFTs(j,2)
        do i = il1,il2
           select case (classpfts(j))

           case ('NdlTr','BdlTr','Crops','BdlSh')
               finerootfr(i,m) = 0.16!1.0 - rootmass(i,m) / (0.6 + rootmass(i,m))
           case ('Grass')  ! grasses
               finerootfr(i,m) = 1.0
           case default
              print * ,'Unknown CLASS PFT in nuptake ',classpfts(j)
              call errorHandler('nuptake', - 4)
           end select
        end do ! loop 240
     end do ! loop 230
  end do ! loop 220

  !> Calculate potential active nitrogen uptake (i.e., max.)
  !! (see equation 8 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ U_{a,pot,NH_4} = \frac{\varepsilon f_r C_R N_{NH_4}}{k_{p,1/2}r_d + N_{NH_4} + N_{NO_3}} \f$
  !! \f$ U_{a,pot,NO_3} = \frac{\varepsilon f_r C_R N_{NO_3}}{k_{p,1/2}r_d + N_{NH_4} + N_{NO_3}} \f$
  do j = 1,icc
     do i = il1,il2

       if (fcancmx(i,j) > 0.0) then

          temp_mineral = kphalf + nh4_mass(i,j) / b_nh4 &
                         + no3_mass(i,j) / b_no3

          if (temp_mineral > 0.0) then
             nuptakeveg_a_ptnl_nh4(i,j) = nh4_mass(i,j) / b_nh4 &
                                          * coeff_a(sort(j)) * (finerootfr(i,j) * rootmass(i,j))/ temp_mineral
             nuptakeveg_a_ptnl_no3(i,j) = no3_mass(i,j) / b_no3 &
                                          * coeff_a(sort(j)) * (finerootfr(i,j) * rootmass(i,j))/ temp_mineral
          end if
          temp_mineral = 0.0

          nuptake_a_ptnl_nh4(i) = nuptake_a_ptnl_nh4(i) + nuptakeveg_a_ptnl_nh4(i,j) &
                                                          * fcancmx(i,j)
          nuptake_a_ptnl_no3(i) = nuptake_a_ptnl_no3(i) + nuptakeveg_a_ptnl_no3(i,j) &
                                                          * fcancmx(i,j)

          if (nuptakeveg_a_ptnl_nh4(i,j) < 0.0) then
             write(*,*)'nuptakeveg_a_ptnl_nh4 lt zero at i=',i,' for pft=',j
             write(*,*)'nuptakeveg_a_ptnl_nh4 = ',nuptakeveg_a_ptnl_nh4(i,j)
             call errorHandler('nuptake', - 5)
          end if
          if (nuptakeveg_a_ptnl_no3(i,j) < 0.0) then
             write(*,*)'nuptakeveg_a_ptnl_no3 lt zero at i=',i,' for pft=',j
             write(*,*)'nuptakeveg_a_ptnl_no3 = ',nuptakeveg_a_ptnl_no3(i,j)
             call errorHandler('nuptake', - 6)
          end if

       end if !fcancmx
     end do ! loop 260
  end do ! loop 250

  !> Calculate actual active nitrogen uptake
  !! (see equation 9 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ U_{a,actual,NH_4} = 0, U_{a,actual,NO_3} = 0\f$ if \f$ \Delta_{WP} \leq U_{p,NH_4} + U_{p,NO_3} \f$
  !! \f$ U_{a,actual,NH_4} = (\Delta_{WP} - U_{p,NH_4} - U_{p,NO_3} \frac{U_{a,pot,NH_4}} {U_{a,pot,NH_4} + U_{a,pot,NO_3}}),\f$
  !! \f$ U_{a,actual,NO_3} = (\Delta_{WP} - U_{p,NH_4} - U_{p,NO_3} \frac{U_{a,pot,NO_3}} {U_{a,pot,NH_4} + U_{a,pot,NO_3}})\f$
  !! if \f$ (\Delta_{WP} > U_{p,NH_4} + U_{p,NO_3}) \cap (\Delta_{WP} < U_{p,NH_4} + U_{p,NO_3} + U_{a,pot,NH_4} + U_{a,pot,NO_3}) \f$
  !! \f$ U_{a,actual,NH_4} = U_{a,pot,NH_4}, U_{a,actual,NO_3} = U_{a,pot,NO_3} \f$ if \f$ (\Delta_{WP} \geq U_{p,NH_4} + U_{p,NO_3} + U_{a,pot,NH_4} + U_{a,pot,NO_3})\f$
  do j = 1,icc
     do i = il1,il2

        if (fcancmx(i,j) > 0.0) then

           total_nuptake_p = nuptakeveg_p_nh4(i,j) + nuptakeveg_p_no3(i,j)
           total_nuptake_a_ptnl = nuptakeveg_a_ptnl_nh4(i,j) + nuptakeveg_a_ptnl_no3(i,j)

          if (ndemandveg_wp_npp(i,j) >= (total_nuptake_p+total_nuptake_a_ptnl)) then

              nuptakeveg_a_actl_nh4(i,j) = nuptakeveg_a_ptnl_nh4(i,j)
              nuptakeveg_a_actl_no3(i,j) = nuptakeveg_a_ptnl_no3(i,j)

          else if ((ndemandveg_wp_npp(i,j) >= total_nuptake_p) .and. &
                  (ndemandveg_wp_npp(i,j) < (total_nuptake_p+total_nuptake_a_ptnl))) then

              nuptakeveg_a_actl_nh4(i,j) = (ndemandveg_wp_npp(i,j) - total_nuptake_p) &
                                    * (nuptakeveg_a_ptnl_nh4(i,j) / total_nuptake_a_ptnl)
              nuptakeveg_a_actl_no3(i,j) = (ndemandveg_wp_npp(i,j) - total_nuptake_p) &
                                    * (nuptakeveg_a_ptnl_no3(i,j) / total_nuptake_a_ptnl)

          else if (ndemandveg_wp_npp(i,j) < total_nuptake_p) then

              nuptakeveg_a_actl_nh4(i,j) = 0.0
              nuptakeveg_a_actl_no3(i,j) = 0.0

          end if
          !
          nuptakeveg_nh4(i,j) = nuptakeveg_p_nh4(i,j) + nuptakeveg_a_actl_nh4(i,j)
          nuptakeveg_no3(i,j) = nuptakeveg_p_no3(i,j) + nuptakeveg_a_actl_no3(i,j)

        end if !fcancmx
     end do ! loop 280
  end do ! loop 270

  !> Updating pools and checking the mass conservation
  do j = 1,icc
     do i = il1,il2

        no3_mass_in(i,j) = no3_mass(i,j)
        nh4_mass_in(i,j) = nh4_mass(i,j)

        nh4_mass(i,j) = nh4_mass(i,j) - nuptakeveg_nh4(i,j)
        no3_mass(i,j) = no3_mass(i,j) - nuptakeveg_no3(i,j)

        !! Adjusting to prevent negative nh4_mass(i,j)
        if (nh4_mass(i,j) < 0.0) then
            ! write(*,*)'nh4_mass lt zero in nuptake.f90 at i=',i,' for pft=',j
            ! write(*,*)'nh4_mass_in(i,j) = ',nh4_mass_in(i,j)
            ! write(*,*)'nh4_mass(i,j) = ',nh4_mass(i,j)
            ! write(*,*)'nuptakeveg_nh4(i,j) = ',nuptakeveg_nh4(i,j)
            ! write(*,*)'adjusting to prevent negative nh4_mass(i,j)'
            temp = nuptakeveg_nh4(i,j)
            nuptakeveg_nh4(i,j) = nuptakeveg_nh4(i,j) - nuptakeveg_nh4(i,j) &
                                                  * abs(nh4_mass(i,j)/temp)
            nh4_mass(i,j) = 0.0
            temp = 0.0
        end if

        !! Adjusting to prevent negative no3_mass(i,j)
        if (no3_mass(i,j) < 0.0) then
            ! write(*,*)'no3_mass lt in nuptake.f90 zero at i=',i,' for pft=',j
            ! write(*,*)'no3_mass(i,j) = ',no3_mass(i,j)
            ! write(*,*)'nuptakeveg_no3(i,j) = ',nuptakeveg_no3(i,j)
            ! write(*,*)'adjusting to prevent negative no3_mass(i,j)'
            temp = nuptakeveg_no3(i,j)
            nuptakeveg_no3(i,j) = nuptakeveg_no3(i,j) - nuptakeveg_no3(i,j) &
                                                  * abs(no3_mass(i,j)/temp)
            no3_mass(i,j) = 0.0
            temp = 0.0
        end if

        !! Checking conservation for nh4_mass(i,j)
        if ((nh4_mass_in(i,j) - nh4_mass(i,j) - nuptakeveg_nh4(i,j)) > ntolerance) then
            write(*,*)'imbalance in NH4+ at (i)=',i,'pft=',j
            call errorHandler('nuptake', - 7)
        end if

        !! Checking conservation for no3_mass(i,j)
        if ((no3_mass_in(i,j) - no3_mass(i,j) - nuptakeveg_no3(i,j)) > ntolerance) then
            write(*,*)'imbalance in NO3- at (i)=',i,'pft=',j
            call errorHandler('nuptake', - 8)
        end if

     end do ! loop 300
  end do ! loop 290

  return
end subroutine nuptake
